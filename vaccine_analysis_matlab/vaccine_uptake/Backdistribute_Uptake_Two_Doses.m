clear;
clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Load information for one dose in order to scale relative to the the
% numebr of two doses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
load('At_Least_One_Dose.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Set up arrays for output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
NS=10^6;

Year_Two=[2011:2022];
Age_Two=[16:18];
Dist_At_Least_Two=zeros(NS,length(Year_Two),length(Age_Two));

% Truncate the one dose in order to map towards the second dose
Dist_At_Least_One=Dist_At_Least_One(:,:,ismember(Age_One,[13:15]));
Age_One=Age_One(ismember(Age_One,[13:15]));

Data_At_Least_Two=readtable(['Vaccination_Data.xlsx'],'Sheet','At_least_two_doses');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the logit standard deviation that best fits the confidence intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std_z2=zeros(height(Data_At_Least_Two),1);

for dd=1:height(Data_At_Least_Two)
    p_t=Data_At_Least_Two.Vac_Uptake(dd);
    l_t=Data_At_Least_Two.Lower_Bound(dd);   
    u_t=Data_At_Least_Two.Upper_Bound(dd); 
    if(~isnan(l_t))
        bnds=[log(l_t./(1-l_t)) log(u_t./(1-u_t)) ];
        options=optimoptions('lsqnonlin','FunctionTolerance',10^(-8),'MaxFunctionEvaluations',10^6,'MaxIterations',10^6,'StepTolerance',10^(-9),'OptimalityTolerance',10^(-10));
        std_z2(dd)=lsqnonlin(@(x)(norminv([0.025 0.975],log(p_t./(1-p_t)),x)-bnds),bnds(2)-log(p_t./(1-p_t)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Second Dose among 17 -yr olds from the first year reported to 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Only have age 17 yrs of age from the data
for yy=min(Data_At_Least_Two.Year):2022
    tf_age=Data_At_Least_Two.Year==yy & Data_At_Least_Two.Age==17;
    p_t=Data_At_Least_Two.Vac_Uptake(tf_age);
    r_age=normrnd(log(p_t./(1-p_t)),std_z2(tf_age),NS,1);
    Dist_At_Least_Two(:,Year_Two==yy,Age_Two==17)=1./(1+exp(-r_age));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Infer the earlier years in whch data was not avaialble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Determine prior year update for 2011-2013: Using linear regression as
% applying the ratio would produce a large amount of uncertainty that would
% accumulate over time. 
t_year=[min(Data_At_Least_Two.Year):min(Data_At_Least_Two.Year)+2]; % Use the periods conssitent iwth prior analysis(i.e. 3 yrs) to approximate the linear slope
y=squeeze(Dist_At_Least_Two(:,ismember(Year_Two,t_year),Age_Two==17));
t_year_m=repmat(t_year,NS,1)-min(Data_At_Least_Two.Year);
y=log(y./(1-y));

mdl=fitlm(t_year_m(:),y(:));

p=flip(mdl.Coefficients.Estimate);

t_year_all=[2010:(min(Data_At_Least_Two.Year)-1)];
z_forecast=polyval(p,t_year_all-min(Data_At_Least_Two.Year));
tf= Data_At_Least_Two.Year<=(min(Data_At_Least_Two.Year)+2) & Data_At_Least_Two.Age==17;
std_forecast=mean(std_z2(tf));

for yy=2011:(min(Data_At_Least_Two.Year)-1)
    z_temp=normrnd(z_forecast(t_year_all==yy),std_forecast,NS,1);
    Dist_At_Least_Two(:,Year_Two==yy,Age_Two==17)=1./(1+exp(-z_temp));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Compute the ratio among 14 and 17 and apply for the other age groups
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% projected coverage if it were to start in 2010 among 17 year olds to
z_temp=normrnd(z_forecast(t_year_all==2010),std_forecast,NS,1);
v_17_2010=1./(1+exp(-z_temp));

v_17_2023=Dist_At_Least_Two(:,Year_Two==2022,Age_Two==17).^2./Dist_At_Least_Two(:,Year_Two==2021,Age_Two==17);

for yy=2011:2022
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % 18 yr old coverage
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    if(yy>2011)
        dOne=Dist_At_Least_One(:,Year_One==(yy-3),Age_One==15)-Dist_At_Least_One(:,Year_One==(yy-4),Age_One==14); % Determine the increase in coverage; this will dictate the proportional increase in two-doses for 18 yr olds
        ratio_14_17=Dist_At_Least_Two(:,Year_Two==(yy-1),Age_Two==17)./Dist_At_Least_One(:,Year_One==(yy-4),Age_One==14); % The proportional increase for this specified birth-cohort 
        Dist_At_Least_Two(:,Year_Two==yy,Age_Two==18)=Dist_At_Least_Two(:,Year_Two==(yy-1),Age_Two==17)+dOne.*ratio_14_17; % The increase in coverage among two-doses based on the increase in coveage from 14 to 15
    else
        dOne=Dist_At_Least_One(:,Year_One==(yy-3),Age_One==15)-Dist_At_Least_One(:,Year_One==(yy-4),Age_One==14); % Determine the increase in coverage; this will dictate the proportional increase in two-doses for 18 yr olds
        ratio_14_17=v_17_2010./Dist_At_Least_One(:,Year_One==(yy-4),Age_One==14); % The proportional increase for this specified birth-cohort 
        Dist_At_Least_Two(:,Year_Two==yy,Age_Two==18)=v_17_2010+dOne.*ratio_14_17; % The increase in coverage among two-doses based on the increase in coveage from 14 to 15
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % 16 yr old coverage
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    if(yy<2022)
        dOne=Dist_At_Least_One(:,Year_One==(yy-2),Age_One==14)-Dist_At_Least_One(:,Year_One==(yy-3),Age_One==13); % Determine the increase in coverage; this will dictate the proportional increase in two-doses for 16 yr olds
        ratio_14_17=Dist_At_Least_Two(:,Year_Two==(yy+1),Age_Two==17)./Dist_At_Least_One(:,Year_One==(yy-2),Age_One==14); % The proportional increase for this specified birth-cohort
        Dist_At_Least_Two(:,Year_Two==yy,Age_Two==16)=Dist_At_Least_Two(:,Year_Two==(yy+1),Age_Two==17)-dOne.*ratio_14_17; % The decrease in coverage among two-doses based on the increase in coveage from 13 to 14
    else
        dOne=Dist_At_Least_One(:,Year_One==(yy-2),Age_One==14)-Dist_At_Least_One(:,Year_One==(yy-3),Age_One==13); % Determine the increase in coverage; this will dictate the proportional increase in two-doses for 16 yr olds
        ratio_14_17=v_17_2023./Dist_At_Least_One(:,Year_One==(yy-2),Age_One==14); % The proportional increase for this specified birth-cohort
        Dist_At_Least_Two(:,Year_Two==yy,Age_Two==16)=v_17_2023-dOne.*ratio_14_17; % The decrease in coverage among two-doses based on the increase in coveage from 13 to 14
    end
end

save('At_Least_Two_Doses.mat','Dist_At_Least_Two','Year_Two','Age_Two')

figure('units','normalized','outerposition',[0 0.3 1 0.4]);
for ii=1:3
    
    subplot(1,3,ii);

    patch([Year_Two flip(Year_Two)],100.*[prctile(squeeze(Dist_At_Least_Two(:,:,ii)),2.5,1) flip(prctile(squeeze(Dist_At_Least_Two(:,:,ii)),97.5,1))],'k','FaceAlpha',0.25,'LineStyle','none');
    hold on;
    plot(Year_Two,100.*prctile(squeeze(Dist_At_Least_Two(:,:,ii)),50,1),'k','LineWidth',2);
    box off;
    ylim([0 100])
    set(gca,'linewidth',2,'tickdir','out','FontSize',14)
    xlabel('Year','FontSize',16)
    ylabel({'Vaccine uptake','(at least two doses)'},'FontSize',16)
    ytickformat("percentage")
    title(['Age: ' num2str(15+ii) ' years'])
    xlim([2011 2022])
end