clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Distrbution of dose among ages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% Ages 11-13
Age_1=[11 12 13]; % Added 10.5 to the 11 yr age class
V1=[0.02272727272727287+0.6681818181818182 0.25000000000000017 0.05909090909090911];
V_10=0.02272727272727287;
V1=V1./sum(V1);

Year_One=[2005:2022];
Age_One=[11:18];

% Ages 16 to 18
Age_2=[16 17 18]; % Added the 15.5 to the 16 yr age class; 
V2=[0.03636363636363642+0.3227272727272728 0.34090909090909105 0.3045454545454547]; 
V2=V2./sum(V2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Set up arrays for output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
NS=10^6;
Dist_At_Least_One=zeros(NS,length(Year_One),length(Age_One));

Data_At_Least_One=readtable(['Vaccination_Data.xlsx'],'Sheet','At_least_one_dose');
Data_Val=Data_At_Least_One(strcmp(Data_At_Least_One.Reference,'Digitized from https://www.cdc.gov/vaccines/imz-managers/coverage/teenvaxview/pubs-presentations/NIS-teen-vac-coverage-estimates-2015-2021.html'),:); % Remove the digitized data as it was an estimate and not survey based
Data_At_Least_One=Data_At_Least_One(~strcmp(Data_At_Least_One.Reference,'Digitized from https://www.cdc.gov/vaccines/imz-managers/coverage/teenvaxview/pubs-presentations/NIS-teen-vac-coverage-estimates-2015-2021.html'),:); % Remove the digitized data as it was an estimate and not survey based

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the logit standard deviation that best fits the confidence intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std_z=zeros(height(Data_At_Least_One),1);

for dd=1:height(Data_At_Least_One)
    p_t=Data_At_Least_One.Vac_Uptake(dd);
    l_t=Data_At_Least_One.Lower_Bound(dd);   
    u_t=Data_At_Least_One.Upper_Bound(dd); 
    if(~isnan(l_t))
        bnds=[log(l_t./(1-l_t)) log(u_t./(1-u_t)) ];
        options=optimoptions('lsqnonlin','FunctionTolerance',10^(-8),'MaxFunctionEvaluations',10^6,'MaxIterations',10^6,'StepTolerance',10^(-9),'OptimalityTolerance',10^(-10));
        std_z(dd)=lsqnonlin(@(x)(norminv([0.025 0.975],log(p_t./(1-p_t)),x)-bnds),bnds(2)-log(p_t./(1-p_t)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ages 11-17 for 2006 to 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for yy=2022:-1:2006
    for aa=17:-1:11
        % Search for the data avaialbe (13-17)
        tf_age=Data_At_Least_One.Year==yy & Data_At_Least_One.Age==aa;
        if(sum(tf_age)>0) % If data found, take a random sample from it
             p_t=Data_At_Least_One.Vac_Uptake(tf_age);
             r_age=normrnd(log(p_t./(1-p_t)),std_z(tf_age),NS,1);
             Dist_At_Least_One(:,Year_One==yy,Age_One==aa)=1./(1+exp(-r_age));
        else   
            % If 11 years of age can only go to 2020; as otherwise run out
            % of data to back distribute
            if(yy<=2020 && aa==11) % Can only base this of two-year ahead since vaccination started in 2005
                p_temp=Dist_At_Least_One(:,Year_One==yy+2,Age_One==13);
                Dist_At_Least_One(:,Year_One==yy,Age_One==aa)=p_temp.*V1(Age_1==11); % Scale the data down based on the back distribution
            elseif(yy<=2021 && aa==12) % If 12 years of age, can only go to 2021; as otherwise run out of data to back distribute
                p_temp=Dist_At_Least_One(:,Year_One==yy+1,Age_One==13);
                Dist_At_Least_One(:,Year_One==yy,Age_One==aa)=p_temp.*sum(V1(Age_1<=12)); % Scale the data down based on the back distribution; but be sure to add the prior uptake from the 11 year age class
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Age 11-12 for 2021 and 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Will infer the uptake for 13 yr olds for future years to approximate
% coverage for 11 and 12 yr olds in 2021 and 2022
t_year=[2020:2022]; % Use the pandemic as the point to evalaute the trajectory (this can be chosen based on the analysis done by CDC that there was an observed change)
y=squeeze(Dist_At_Least_One(:,ismember(Year_One,t_year),Age_One==13));
t_year_m=repmat(t_year,NS,1)-2020;
y=log(y./(1-y));

mdl=fitlm(t_year_m(:),y(:));

p=flip(mdl.Coefficients.Estimate);

t_year_all=[2023:2024];
z_forecast=polyval(p,t_year_all-2020);
tf= Data_At_Least_One.Year>=2020 & Data_At_Least_One.Age==13;
std_forecast=mean(std_z(tf));

z_2023=normrnd(z_forecast(t_year_all==2023),std_forecast,NS,1);
Vac_uptake_13_2023=1./(1+exp(-z_2023));

z_2024=normrnd(z_forecast(t_year_all==2024),std_forecast,NS,1);
Vac_uptake_13_2024=1./(1+exp(-z_2024));
% Back distribute for the specified years
for aa=11:12
    for yy=2021:2022
        if(yy==2021 && aa==11) % Age 11 for yr 2021 need to use 2023 for 13 yr old
            p_temp=Vac_uptake_13_2023;
            Dist_At_Least_One(:,Year_One==yy,Age_One==aa)=p_temp.*V1(Age_1==aa);        
        elseif(yy==2022 && aa==11) % Age 11 for yr 2022 need to use 2024 for 13 yr old
            p_temp=Vac_uptake_13_2024;
            Dist_At_Least_One(:,Year_One==yy,Age_One==aa)=p_temp.*V1(Age_1==aa);    
        elseif(yy==2022 && aa==12) % Age 12 for yr 2022 need to use 2023 for 13 yr old           
            p_temp=Vac_uptake_13_2023;   
            Dist_At_Least_One(:,Year_One==yy,Age_One==aa)=p_temp.*sum(V1(Age_1<=12));
        end
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% % 2005 Uptake Data for 11-17
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% Use the prior distribtuin for the back dsitribution for 11 yr olds
Dist_At_Least_One(:,Year_One==2005,Age_One==11)=Dist_At_Least_One(:,Year_One==(2005+1),Age_One==12).*sum(V1(Age_1==11))./sum(sum(V1(Age_1<=12)));
% Use the prior distribtuin for the back dsitribution for 12 yr olds
Dist_At_Least_One(:,Year_One==2005,Age_One==12)=Dist_At_Least_One(:,Year_One==(2005+1),Age_One==13).*sum(V1(Age_1==12))./sum(sum(V1(Age_1>=12)));
% Use the prior distribtuin for the back dsitribution for 16 yr olds
Dist_At_Least_One(:,Year_One==2005,Age_One==16)=Dist_At_Least_One(:,Year_One==(2005+1),Age_One==17).*sum(V1(Age_2==16))./sum(sum(V1(Age_2<=17)));


% Use this wieght to inform the extent of vaccination for 13 year olds
% % https://bmcinfectdis.biomedcentral.com/articles/10.1186/1471-2334-9-175
V_13_16_2005=0.082;

US_Pop_2005=[4149228	4251162	4308851	4353216	4456969	4308145	4239509	4225126];
US_Age=[11:18];
w=[US_Pop_2005(US_Age==13) US_Pop_2005(US_Age==16)];
w=w./(sum(w));
Dist_At_Least_One(:,Year_One==2005,Age_One==13)=(V_13_16_2005-Dist_At_Least_One(:,Year_One==2005,Age_One==16).*w(2))./w(1);

% Bootstrap non-negative samples to replace the negative values
t_neg=find(Dist_At_Least_One(:,Year_One==2005,Age_One==13)<0);
t_pos=find(Dist_At_Least_One(:,Year_One==2005,Age_One==13)>0);
Dist_At_Least_One(t_neg,Year_One==2005,Age_One==13)=Dist_At_Least_One(t_pos(randi(length(t_pos),length(t_neg),1)),Year_One==2005,Age_One==13);


% Approximate what the coverage would possibley be for 14 in 2005
% based on the growth ratio from prior years among the 13 yr old age class
p_2006=Dist_At_Least_One(:,Year_One==2006,Age_One==13);
p_2005=Dist_At_Least_One(:,Year_One==2005,Age_One==13);
r=p_2005./p_2006;

% Construct the distribution for 14 yr olds
Dist_At_Least_One(:,Year_One==2005,Age_One==14)=Dist_At_Least_One(:,Year_One==2006,Age_One==14).*r;

% Use this prior information to construct the coverage for 15 year olds
V_14_15_2005=0.105;
w=[US_Pop_2005(US_Age==14) US_Pop_2005(US_Age==15)];
w=w./(sum(w));
Dist_At_Least_One(:,Year_One==2005,Age_One==15)=(V_14_15_2005-Dist_At_Least_One(:,Year_One==2005,Age_One==14).*w(1))./w(2);

% Bootstrap non-negative samples to replace the negative values
t_neg=find(Dist_At_Least_One(:,Year_One==2005,Age_One==15)<0);
t_pos=find(Dist_At_Least_One(:,Year_One==2005,Age_One==15)>0);
Dist_At_Least_One(t_neg,Year_One==2005,Age_One==15)=Dist_At_Least_One(t_pos(randi(length(t_pos),length(t_neg),1)),Year_One==2005,Age_One==15);

% Approximate what the coverage would possibley be for 17 in 2005
% based on the growth ratio from prior years among the 16 yr old age class
p_2006=Dist_At_Least_One(:,Year_One==2006,Age_One==16);
p_2005=Dist_At_Least_One(:,Year_One==2005,Age_One==16);
r=p_2005./p_2006;
Dist_At_Least_One(:,Year_One==2005,Age_One==17)=Dist_At_Least_One(:,Year_One==2006,Age_One==17).*r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% 18 year old uptake
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for yy=2022:-1:2005
    if(yy>2007 ) % If the year is greater than 2007 we can use information back to they yr they were 15 to propgate coverage
        % The weight in the cited manuscript is based on timing from ~16;
        % so we have to apply the wieght to "NEW VACCINES" starting at 16.
        Dist_At_Least_One(:,Year_One==yy,Age_One==18)=Dist_At_Least_One(:,Year_One==(yy-3),Age_One==15)+(Dist_At_Least_One(:,Year_One==(yy-1),Age_One==17)-Dist_At_Least_One(:,Year_One==(yy-3),Age_One==15))./sum(V2(Age_2<18));
    elseif(yy==2007) % If the year is 2007, can only go to the point in which they are 16
        w_temp=V2(Age_2==17)./sum(V2(Age_2>=17));
        Dist_At_Least_One(:,Year_One==yy,Age_One==18)=Dist_At_Least_One(:,Year_One==(yy-2),Age_One==16)+(Dist_At_Least_One(:,Year_One==(yy-1),Age_One==17)-Dist_At_Least_One(:,Year_One==(yy-2),Age_One==16))./w_temp;
    elseif(yy==2005)
        % Approximate what the coverage would possibley be for 18 in 2005
        % based on the growth ratio from prior years among the 17 yr old age class
        p_2005=Dist_At_Least_One(:,Year_One==yy,Age_One==17);
        p_2006=Dist_At_Least_One(:,Year_One==yy+1,Age_One==17);
        r=p_2005./p_2006;
        Dist_At_Least_One(:,Year_One==yy,Age_One==18)=Dist_At_Least_One(:,Year_One==yy+1,Age_One==18).*r;
    elseif(yy==2006)
        % Use the relative increase in uptake from 2006 17 yr old to 2007 18 yr old
        % to infer increae among 18 yr old for 2006 from 17 yr old 2005
        p_17=Dist_At_Least_One(:,Year_One==(yy),Age_One==17);
        p_18=Dist_At_Least_One(:,Year_One==(yy+1),Age_One==18);
        r=p_18./p_17;
        Dist_At_Least_One(:,Year_One==yy,Age_One==18)=Dist_At_Least_One(:,Year_One==yy-1,Age_One==17).*r;
    end
end
save('At_Least_One_Dose.mat','Dist_At_Least_One','Year_One','Age_One')
close all;

figure('units','normalized','outerposition',[0 0 1 1]);
for ii=1:8
    tf_age=Data_Val.Age==Age_One(ii); % using the 13-yr old age class to determine the possible coverage
    pt=Data_Val.Vac_Uptake(tf_age);
    yr_plot=Data_Val.Year(tf_age);

    subplot(4,2,ii);

    patch([Year_One flip(Year_One)],100.*[prctile(squeeze(Dist_At_Least_One(:,:,ii)),2.5,1) flip(prctile(squeeze(Dist_At_Least_One(:,:,ii)),97.5,1))],'b','FaceAlpha',0.25,'LineStyle','none');
    hold on;
    plot(Year_One,100.*prctile(squeeze(Dist_At_Least_One(:,:,ii)),50,1),'b','LineWidth',2);
    scatter(yr_plot,100.*pt,10,'r','filled');
    box off;
    ylim([0 100])
    set(gca,'linewidth',2,'tickdir','out','FontSize',14)
    xlabel('Year','FontSize',16)
    ylabel({'Vaccine uptake','(at least one dose)'},'FontSize',16)
    ytickformat("percentage")
    title(['Age: ' num2str(10+ii) ' years'])
    xlim([2005 2022])
end


figure('units','normalized','outerposition',[0 0 1 1]);
for ii=1:8
    tf_age=Data_Val.Age==Age_One(ii); % using the 13-yr old age class to determine the possible coverage
    pt=Data_Val.Vac_Uptake(tf_age);
    yr_plot=Data_Val.Year(tf_age);

    subplot(4,2,ii);

    patch([Year_One flip(Year_One)],log([prctile(squeeze(Dist_At_Least_One(:,:,ii)),2.5,1) flip(prctile(squeeze(Dist_At_Least_One(:,:,ii)),97.5,1))]./(1-[prctile(squeeze(Dist_At_Least_One(:,:,ii)),2.5,1) flip(prctile(squeeze(Dist_At_Least_One(:,:,ii)),97.5,1))])),'b','FaceAlpha',0.25,'LineStyle','none');
    hold on;
    plot(Year_One,log(prctile(squeeze(Dist_At_Least_One(:,:,ii)),50,1)./(1-prctile(squeeze(Dist_At_Least_One(:,:,ii)),50,1))),'b','LineWidth',2);
    scatter(yr_plot,log(pt./(1-pt)),10,'r','filled');
    box off;
    ylim([-3 3])
    set(gca,'linewidth',2,'tickdir','out','FontSize',14)
    xlabel('Year','FontSize',16)
    ylabel({'logit Vaccine uptake','(at least one dose)'},'FontSize',16)
%     ytickformat("percentage")
    title(['Age: ' num2str(10+ii) ' years'])
    xlim([2005 2022])
end

