clear;
clc;

% load the data from the excel sheet
ABC=readtable('Vaccination_Incidence_Data.xlsx','Sheet','ABCs');
Nat=readtable('Vaccination_Incidence_Data.xlsx','Sheet','National_Vaccinated');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arrays to store the inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Population.ABC.Unvaccinated.Age_11_15=zeros(height(ABC),1);
Population.ABC.Unvaccinated.Age_16_23=zeros(height(ABC),1);

Population.ABC.Vaccinated_One_Dose.Age_11_15=zeros(height(ABC),1);
Population.ABC.Vaccinated_Two_Dose.Age_11_15=zeros(height(ABC),1);
Population.ABC.Vaccinated_One_Dose.Age_16_23=zeros(height(ABC),1);
Population.ABC.Vaccinated_Two_Dose.Age_16_23=zeros(height(ABC),1);

Population.National.Age_11_15=zeros(height(ABC),1);
Population.National.Age_16_23=zeros(height(ABC),1);

Reported_Cases.ABC.Unvaccinated.Age_11_15=zeros(height(ABC),1);
Reported_Cases.ABC.Unvaccinated.Age_16_23=zeros(height(ABC),1);

Reported_Cases.ABC.Vaccinated_One_Dose.Age_11_15=zeros(height(ABC),1);
Reported_Cases.ABC.Vaccinated_Two_Dose.Age_11_15=zeros(height(ABC),1);
Reported_Cases.ABC.Vaccinated_One_Dose.Age_16_23=zeros(height(ABC),1);
Reported_Cases.ABC.Vaccinated_Two_Dose.Age_16_23=zeros(height(ABC),1);

Reported_Cases.National.Age_11_15=zeros(height(ABC),1);
Reported_Cases.National.Age_16_23=zeros(height(ABC),1);

Vaccine_Coverage.One_Dose.Age_11_15=zeros(height(ABC),1);
Vaccine_Coverage.Two_Dose.Age_11_15=zeros(height(ABC),1);
Vaccine_Coverage.One_Dose.Age_16_23=zeros(height(ABC),1);
Vaccine_Coverage.Two_Dose.Age_16_23=zeros(height(ABC),1);

% 2014 to 2021 Data
Reported_Cases_2014_2021.Unvaccinated.Age_11_15=zeros(height(Nat),1);
Reported_Cases_2014_2021.Vaccinated_One_Dose.Age_11_15=zeros(height(Nat),1);
Reported_Cases_2014_2021.Vaccinated_Two_Dose.Age_11_15=zeros(height(Nat),1);

Reported_Cases_2014_2021.Unvaccinated.Age_16_23=zeros(height(Nat),1);
Reported_Cases_2014_2021.Vaccinated_One_Dose.Age_16_23=zeros(height(Nat),1);
Reported_Cases_2014_2021.Vaccinated_Two_Dose.Age_16_23=zeros(height(Nat),1);

Population_2014_2021.National.Unvaccinated.Age_11_15=zeros(height(Nat),1);
Population_2014_2021.National.Unvaccinated.Age_16_23=zeros(height(Nat),1);

Population_2014_2021.National.Vaccinated_One_Dose.Age_11_15=zeros(height(Nat),1);
Population_2014_2021.National.Vaccinated_Two_Dose.Age_11_15=zeros(height(Nat),1);

Population_2014_2021.National.Vaccinated_One_Dose.Age_16_23=zeros(height(Nat),1);
Population_2014_2021.National.Vaccinated_Two_Dose.Age_16_23=zeros(height(Nat),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read data from tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Year_2005_2013=ABC.Year;
Year_2014_2021=Nat.Year;

Reported_Cases.Additional.Years=[2005:2008];
Reported_Cases.Additional.Vaccinated_2005_to_2008.Age_11_15=1; % Table 1 https://journals.lww.com/pidj/abstract/2011/06000/early_estimate_of_the_effectiveness_of.2.aspx
Reported_Cases.Additional.Vaccinated_2005_to_2008.Age_16_23=13; % Table 1 https://journals.lww.com/pidj/abstract/2011/06000/early_estimate_of_the_effectiveness_of.2.aspx
Population.Additional.Vaccinated_2005_to_2008.Fraction_National=0.54; %https://journals.lww.com/pidj/abstract/2011/06000/early_estimate_of_the_effectiveness_of.2.aspx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABCs site data (2005-2013)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for yy=1:height(ABC)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Age 11 to 15
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Approximate the population size of 11 to 15 for the ABC region
    
    N_ABC_11_15=round(ABC.FractionOfPopulation(yy).*ABC.USAge11To15(yy));
    
    % The vaccination coverage among the age class
    Vaccine_Coverage.One_Dose.Age_11_15(yy)=ABC.CoverageOneDoseAges11To15(yy); 
    Vaccine_Coverage.Two_Dose.Age_11_15(yy)=ABC.CoverageTwoDosesAges11To15(yy); 
        
    % Unvaccincated cases in ABC region
    Reported_Cases.ABC.Unvaccinated.Age_11_15(yy)=ABC.ABCCasesUnvaccinatedAges11To15(yy);

    % Vaccincated cases in ABC region for one dose
    Reported_Cases.ABC.Vaccinated_One_Dose.Age_11_15(yy)=ABC.ABCCasesVaccinatedOneDoseAges11To15(yy)+ABC.ABCCasesVaccinatedUnknownDoseAges11To15(yy);    

    % Vaccincated cases in ABC region for two dose
    Reported_Cases.ABC.Vaccinated_Two_Dose.Age_11_15(yy)=ABC.ABCCasesVaccinatedTwoDoseAges11To15(yy);

    % Total population of age class
    Population.National.Age_11_15(yy)=ABC.USAge11To15(yy); 

    % Total cases of age class
    Reported_Cases.National.Age_11_15(yy)=ABC.NationalCasesAges11To15(yy); 

    % Unvaccinated population in the ABC regions
    Population.ABC.Unvaccinated.Age_11_15(yy)=N_ABC_11_15-round(N_ABC_11_15.*Vaccine_Coverage.One_Dose.Age_11_15(yy))-round(N_ABC_11_15.*Vaccine_Coverage.Two_Dose.Age_11_15(yy));
    % Vaccinated population one dose in the ABC regions
    Population.ABC.Vaccinated_One_Dose.Age_11_15(yy)=round(N_ABC_11_15.*Vaccine_Coverage.One_Dose.Age_11_15(yy));
     % Vaccinated population two dose in the ABC regions
    Population.ABC.Vaccinated_Two_Dose.Age_11_15(yy)=round(N_ABC_11_15.*Vaccine_Coverage.Two_Dose.Age_11_15(yy));
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Age 16 to 23
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Approximate the population size of 16 to 23 for the ABC region
    
    N_ABC_16_23=round(ABC.FractionOfPopulation(yy).*ABC.USAge16To23(yy));
    
    % The vaccination coverage among the age class
    Vaccine_Coverage.One_Dose.Age_16_23(yy)=ABC.CoverageOneDoseAges16To23(yy); 
    Vaccine_Coverage.Two_Dose.Age_16_23(yy)=ABC.CoverageTwoDosesAges16To23(yy); 
        
    % Unvaccincated cases in ABC region
    Reported_Cases.ABC.Unvaccinated.Age_16_23(yy)=ABC.ABCCasesUnvaccinatedAges16To23(yy);

    % Vaccincated cases in ABC region for one dose
    Reported_Cases.ABC.Vaccinated_One_Dose.Age_16_23(yy)=ABC.ABCCasesVaccinatedOneDoseAges16To23(yy)+ABC.ABCCasesVaccinatedUnknownDoseAges16To23(yy);    

    % Vaccincated cases in ABC region for two dose
    Reported_Cases.ABC.Vaccinated_Two_Dose.Age_16_23(yy)=ABC.ABCCasesVaccinatedTwoDoseAges16To23(yy);

    % Total population of age class
    Population.National.Age_16_23(yy)=ABC.USAge16To23(yy); 

    % Total cases of age class
    Reported_Cases.National.Age_16_23(yy)=ABC.NationalCasesAges16To23(yy); 

    % Unvaccinated population in the ABC regions
    Population.ABC.Unvaccinated.Age_16_23(yy)=N_ABC_16_23-round(N_ABC_16_23.*Vaccine_Coverage.One_Dose.Age_16_23(yy))-round(N_ABC_16_23.*Vaccine_Coverage.Two_Dose.Age_16_23(yy));
    % Vaccinated population one dose in the ABC regions
    Population.ABC.Vaccinated_One_Dose.Age_16_23(yy)=round(N_ABC_16_23.*Vaccine_Coverage.One_Dose.Age_16_23(yy));
     % Vaccinated population two dose in the ABC regions
    Population.ABC.Vaccinated_Two_Dose.Age_16_23(yy)=round(N_ABC_16_23.*Vaccine_Coverage.Two_Dose.Age_16_23(yy));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% National data (2014-2021)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Reported_Cases_2014_2021.Unvaccinated.Age_11_15=Nat.NationCasesUnvaccinated11To15;
    Reported_Cases_2014_2021.Vaccinated_One_Dose.Age_11_15=Nat.NationalCasesOneDose11To15+Nat.NationalCasesUnknownDose11To15;
    Reported_Cases_2014_2021.Vaccinated_Two_Dose.Age_11_15=Nat.NationalCasesTwoDoses11To15;
    
    Reported_Cases_2014_2021.Unvaccinated.Age_16_23=Nat.NationCasesUnvaccinated16To23;
    Reported_Cases_2014_2021.Vaccinated_One_Dose.Age_16_23=Nat.NationalCasesOneDose16To23+Nat.NationalCasesUnknownDose16To23;
    Reported_Cases_2014_2021.Vaccinated_Two_Dose.Age_16_23=Nat.NationalCasesTwoDoses16To23;
        
    Population_2014_2021.National.Unvaccinated.Age_11_15=Nat.USAge11To15-round(Nat.USAge11To15.*Nat.CoverageOneDoseAges11To15)-round(Nat.USAge11To15.*Nat.CoverageTwoDosesAges11To15);
    Population_2014_2021.National.Vaccinated_One_Dose.Age_11_15=round(Nat.USAge11To15.*Nat.CoverageOneDoseAges11To15);
    Population_2014_2021.National.Vaccinated_Two_Dose.Age_11_15=round(Nat.USAge11To15.*Nat.CoverageTwoDosesAges11To15);
    
    Population_2014_2021.National.Unvaccinated.Age_16_23=Nat.USAge16To23-round(Nat.USAge16To23.*Nat.CoverageOneDoseAges16To23)-round(Nat.USAge16To23.*Nat.CoverageTwoDosesAges16To23);
    Population_2014_2021.National.Vaccinated_One_Dose.Age_16_23=round(Nat.USAge16To23.*Nat.CoverageOneDoseAges16To23);
    Population_2014_2021.National.Vaccinated_Two_Dose.Age_16_23=round(Nat.USAge16To23.*Nat.CoverageTwoDosesAges16To23);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options=optimoptions('surrogateopt','MaxFunctionEvaluations',500,'useparallel',true);

% Restrictions on the probabilities to ensure that p_vac<=p_inf
A=[-1 0 1 0;
    0 -1 0 1];
b=[0;0];

[x0,fval_0]=surrogateopt(@(x)Objective_Function(x,Year_2005_2013,Population,Reported_Cases,Vaccine_Coverage,Reported_Cases_2014_2021,Population_2014_2021),[0 0 0 0],[1 1 1 1],[],A,b,[],[],options);

options=optimoptions("fmincon",'ConstraintTolerance',10^(-9),'FunctionTolerance',10^(-9),'MaxFunctionEvaluations',10^6,'MaxIterations',10^6,'StepTolerance',10^(-9));
[p_estimate,f_mle]=fmincon(@(x)Objective_Function(x,Year_2005_2013,Population,Reported_Cases,Vaccine_Coverage,Reported_Cases_2014_2021,Population_2014_2021),x0,A,b,[],[],[0 0 0 0],[1 1 1 1],[],options);

if(fval_0<f_mle)
    p_estimate=x0;
end
% Scales the estimates from the optimization 
p_inf_11_15=p_estimate(1)./10^5;
p_inf_16_23=p_estimate(2)./10^5;

p_vac_11_15=p_estimate(3)./10^5;
p_vac_16_23=p_estimate(4)./10^5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the number of cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Searching for a reasonable space to conduct the optimization 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vaccinated_One_Dose_Cases_Age_11_to_15=zeros(length(Year_2005_2013),1);
Vaccinated_Two_Dose_Cases_Age_11_to_15=zeros(length(Year_2005_2013),1);
Vaccinated_One_Dose_Cases_Age_16_to_23=zeros(length(Year_2005_2013),1);
Vaccinated_Two_Dose_Cases_Age_16_to_23=zeros(length(Year_2005_2013),1);

for yy=1:length(Year_2005_2013)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Age 11 to 15
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
     % Compute the log-likelihood for the different possibilities 
    Vac_Case_1=[Reported_Cases.ABC.Vaccinated_One_Dose.Age_11_15(yy):Reported_Cases.National.Age_11_15(yy)];
    Vac_Case_2=[Reported_Cases.ABC.Vaccinated_Two_Dose.Age_11_15(yy):Reported_Cases.National.Age_11_15(yy)];
    [Vac_Case_1,Vac_Case_2]=meshgrid(Vac_Case_1,Vac_Case_2);
    Vac_Case_1=Vac_Case_1(:);
    Vac_Case_2=Vac_Case_2(:);

    tempI=Vac_Case_1+Vac_Case_2;
    Vac_Case_1=Vac_Case_1(tempI<=Reported_Cases.National.Age_11_15(yy));
    Vac_Case_2=Vac_Case_2(tempI<=Reported_Cases.National.Age_11_15(yy));

    L=zeros(length(Vac_Case_2),1);

    for ii=1:length(Vac_Case_2)
        L(ii)=log(binopdf(Reported_Cases.National.Age_11_15(yy)-Vac_Case_1(ii)-Vac_Case_2(ii),Population.National.Age_11_15(yy)-round(Population.National.Age_11_15(yy).*Vaccine_Coverage.One_Dose.Age_11_15(yy))-round(Population.National.Age_11_15(yy).*Vaccine_Coverage.Two_Dose.Age_11_15(yy)),p_inf_11_15));
        L(ii)=L(ii)+log(binopdf(Vac_Case_1(ii),round(Population.National.Age_11_15(yy).*Vaccine_Coverage.One_Dose.Age_11_15(yy)),p_vac_11_15));
        L(ii)=L(ii)+log(binopdf(Vac_Case_2(ii),round(Population.National.Age_11_15(yy).*Vaccine_Coverage.Two_Dose.Age_11_15(yy)),p_vac_11_15));

        if(Vac_Case_1(ii)==0)  
            L(ii)=L(ii)+(round(Population.National.Age_11_15(yy).*Vaccine_Coverage.One_Dose.Age_11_15(yy))).*log(1-p_vac_11_15);
        else
            L(ii)=L(ii)+log(1-poisscdf(Reported_Cases.ABC.Vaccinated_One_Dose.Age_11_15(yy),Vac_Case_1(ii)));
        end

        if(Vac_Case_2(ii)==0)  
            L(ii)=L(ii)+(round(Population.National.Age_11_15(yy).*Vaccine_Coverage.Two_Dose.Age_11_15(yy))).*log(1-p_vac_11_15);
        else
            L(ii)=L(ii)+log(1-poisscdf(Reported_Cases.ABC.Vaccinated_Two_Dose.Age_11_15(yy),Vac_Case_2(ii)));
        end

    end
    
    % Record the cases that maximize the log-likelihood 
    Vaccinated_One_Dose_Cases_Age_11_to_15(yy)=Vac_Case_1(L==max(L));
    Vaccinated_Two_Dose_Cases_Age_11_to_15(yy)=Vac_Case_2(L==max(L));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Age 16 to 23
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
     % Compute the log-likelihood for the different possibilities 
    Vac_Case_1=[Reported_Cases.ABC.Vaccinated_One_Dose.Age_16_23(yy):Reported_Cases.National.Age_16_23(yy)];
    Vac_Case_2=[Reported_Cases.ABC.Vaccinated_Two_Dose.Age_16_23(yy):Reported_Cases.National.Age_16_23(yy)];
    [Vac_Case_1,Vac_Case_2]=meshgrid(Vac_Case_1,Vac_Case_2);
    Vac_Case_1=Vac_Case_1(:);
    Vac_Case_2=Vac_Case_2(:);

    tempI=Vac_Case_1+Vac_Case_2;
    Vac_Case_1=Vac_Case_1(tempI<=Reported_Cases.National.Age_16_23(yy));
    Vac_Case_2=Vac_Case_2(tempI<=Reported_Cases.National.Age_16_23(yy));
    L=zeros(length(Vac_Case_2),1);

    for ii=1:length(Vac_Case_2)
        L(ii)=log(binopdf(Reported_Cases.National.Age_16_23(yy)-Vac_Case_1(ii)-Vac_Case_2(ii),Population.National.Age_16_23(yy)-round(Population.National.Age_16_23(yy).*Vaccine_Coverage.One_Dose.Age_16_23(yy))-round(Population.National.Age_16_23(yy).*Vaccine_Coverage.Two_Dose.Age_16_23(yy)),p_inf_16_23));
        L(ii)=L(ii)+log(binopdf(Vac_Case_1(ii),round(Population.National.Age_16_23(yy).*Vaccine_Coverage.One_Dose.Age_16_23(yy)),p_vac_16_23));
        L(ii)=L(ii)+log(binopdf(Vac_Case_2(ii),round(Population.National.Age_16_23(yy).*Vaccine_Coverage.Two_Dose.Age_16_23(yy)),p_vac_16_23));

        if(Vac_Case_1(ii)==0)  
            L(ii)=L(ii)+(round(Population.National.Age_16_23(yy).*Vaccine_Coverage.One_Dose.Age_16_23(yy))).*log(1-p_vac_16_23);
        else
            L(ii)=L(ii)+log(1-poisscdf(Reported_Cases.ABC.Vaccinated_One_Dose.Age_16_23(yy),Vac_Case_1(ii)));
        end

        if(Vac_Case_2(ii)==0)  
            L(ii)=L(ii)+(round(Population.National.Age_16_23(yy).*Vaccine_Coverage.Two_Dose.Age_16_23(yy))).*log(1-p_vac_16_23);
        else
            L(ii)=L(ii)+log(1-poisscdf(Reported_Cases.ABC.Vaccinated_Two_Dose.Age_16_23(yy),Vac_Case_2(ii)));
        end

    end
    
    % Record the cases that maximize the log-likelihood 
    Vaccinated_One_Dose_Cases_Age_16_to_23(yy)=Vac_Case_1(L==max(L));
    Vaccinated_Two_Dose_Cases_Age_16_to_23(yy)=Vac_Case_2(L==max(L));
end
x_approx=[Vaccinated_One_Dose_Cases_Age_11_to_15(:); Vaccinated_Two_Dose_Cases_Age_11_to_15(:);Vaccinated_One_Dose_Cases_Age_16_to_23(:); Vaccinated_Two_Dose_Cases_Age_16_to_23(:)]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the number of cases through optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp_c=Reported_Cases.ABC.Vaccinated_Two_Dose.Age_16_23;

LBC=[Reported_Cases.ABC.Vaccinated_One_Dose.Age_11_15' Reported_Cases.ABC.Vaccinated_Two_Dose.Age_11_15' Reported_Cases.ABC.Vaccinated_One_Dose.Age_16_23' Reported_Cases.ABC.Vaccinated_Two_Dose.Age_16_23'];

UBC=[Reported_Cases.National.Age_11_15'-Reported_Cases.ABC.Vaccinated_Two_Dose.Age_11_15' Reported_Cases.National.Age_11_15'-Reported_Cases.ABC.Vaccinated_One_Dose.Age_11_15' Reported_Cases.National.Age_16_23'-Reported_Cases.ABC.Vaccinated_Two_Dose.Age_16_23' Reported_Cases.National.Age_16_23'-Reported_Cases.ABC.Vaccinated_One_Dose.Age_16_23'];

A=zeros(2+2.*length(Year_2005_2013),2.*(length(Reported_Cases.ABC.Vaccinated_One_Dose.Age_16_23)+length(Reported_Cases.ABC.Vaccinated_Two_Dose.Age_16_23)));
for yy=1:length(Year_2005_2013)
    A(yy,yy)=1;
    A(yy,yy+length(Reported_Cases.ABC.Vaccinated_One_Dose.Age_11_15))=1;
    A(yy+length(Year_2005_2013),yy+length(Reported_Cases.ABC.Vaccinated_One_Dose.Age_11_15)+length(Reported_Cases.ABC.Vaccinated_Two_Dose.Age_11_15))=1;
    A(yy+length(Year_2005_2013),yy+length(Reported_Cases.ABC.Vaccinated_One_Dose.Age_11_15)+length(Reported_Cases.ABC.Vaccinated_Two_Dose.Age_11_15)+length(Reported_Cases.ABC.Vaccinated_One_Dose.Age_16_23))=1;
end

% Need to specify the lower bound 
temp_11_15=zeros(size(LBC));
temp_16_23=zeros(size(LBC));
temp_11_15(find(ismember(Year_2005_2013,Reported_Cases.Additional.Years)))=-1;
temp_11_15(length(Reported_Cases.ABC.Vaccinated_One_Dose.Age_11_15)+find(ismember(Year_2005_2013,Reported_Cases.Additional.Years)))=-1;
temp_16_23(length(Reported_Cases.ABC.Vaccinated_One_Dose.Age_11_15)+length(Reported_Cases.ABC.Vaccinated_Two_Dose.Age_11_15)+find(ismember(Year_2005_2013,Reported_Cases.Additional.Years)))=-1;
temp_16_23(length(Reported_Cases.ABC.Vaccinated_One_Dose.Age_11_15)+length(Reported_Cases.ABC.Vaccinated_Two_Dose.Age_11_15)+length(Reported_Cases.ABC.Vaccinated_One_Dose.Age_16_23)+find(ismember(Year_2005_2013,Reported_Cases.Additional.Years)))=-1;

A(end-1,:)=temp_11_15;
A(end,:)=temp_16_23;

b=[Reported_Cases.National.Age_11_15(:); Reported_Cases.National.Age_16_23(:) ; -Reported_Cases.Additional.Vaccinated_2005_to_2008.Age_11_15; -Reported_Cases.Additional.Vaccinated_2005_to_2008.Age_16_23];

% Prior convergece;
xt= [0	0	2	3	3	1	1	2	1	0	0	0	0	0	0	0	0	0	0	2	5	5	4	6	5	5	4	0	0	0	1	0	0	0	1	2];
xt2=[0	0	2	3	3	1	1	2	1	0	0	0	0	0	0	0	0	0	0	2	5	5	4	6	5	5	4	0	0	0	1	0	0	0	1	2];
xt3= [0	0	2	3	3	1	1	2	1	0	0	0	0	0	0	0	0	0	0	2	5	5	4	6	5	5	4	0	0	0	1	0	0	0	1	2];
xt4= [1	0	2	2	2	1	1	1	1	0	0	0	0	0	0	0	0	0	0	1	4	7	4	5	6	6	3	0	0	0	1	0	0	2	2	2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% REDUCE THE SIZE OF THE UPPER BOUND DURING THE REFINING PROCESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UBC=min(UBC,10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Approximation
% xt2=[Vaccinated_One_Dose_Cases_Age_11_to_15' Vaccinated_Two_Dose_Cases_Age_11_to_15'  Vaccinated_One_Dose_Cases_Age_16_to_23' Vaccinated_Two_Dose_Cases_Age_16_to_23'];
x0=zeros(1500,length(xt));

% Conduct a Greedy search
f_new=-Inf;
f_start=Objective_Function_Cases(xt,Year_2005_2013,Population,Reported_Cases,Vaccine_Coverage,p_inf_11_15,p_inf_16_23,p_vac_11_15,p_vac_16_23);
x_new=xt;
while(f_new<f_start)
    if(isinf(f_new))
        f_new=f_start;
    else
        f_start=f_new;
    end
    ov=randperm(length(xt));
    for jj=1:length(xt)
        L_temp=zeros(1+UBC(ov(jj))-LBC(ov(jj)),1);
        for ww=LBC(ov(jj)):UBC(ov(jj))
            x_temp=x_new;
            x_temp(ov(jj))=ww;
            L_temp(1+ww-LBC(ov(jj)))=Objective_Function_Cases(x_temp,Year_2005_2013,Population,Reported_Cases,Vaccine_Coverage,p_inf_11_15,p_inf_16_23,p_vac_11_15,p_vac_16_23);
        end
        if(min(L_temp)<f_new)
            f_new=min(L_temp);
            x_temp=x_new;
            varx=LBC(ov(jj)):UBC(ov(jj));
            x_temp(ov(jj))=varx(L_temp==min(L_temp));
            x_new=x_temp;
        end
    end
end

if(max(abs(x_new-xt))==0)
    x_new=xt4;
end

% Pertrub the solutions
for ii=1:500
    x_temp=xt+(randi(5,1,length(xt))-3);
%     x_temp=xt+(randi(3,1,length(xt))-2);
    x_temp(x_temp<LBC)=LBC(x_temp<LBC);
    x_temp(x_temp>UBC)=UBC(x_temp>UBC);
    x0(1+3.*(ii-1),:)=x_temp;

    x_temp=x_new+(randi(7,1,length(x_new))-4);
    x_temp(x_temp<LBC)=LBC(x_temp<LBC);
    x_temp(x_temp>UBC)=UBC(x_temp>UBC);
    x0(2+3.*(ii-1),:)=x_temp;

    x_temp=x_approx+(randi(7,1,length(x_approx))-4);
    x_temp(x_temp<LBC)=LBC(x_temp<LBC);
    x_temp(x_temp>UBC)=UBC(x_temp>UBC);
    x0(3+3.*(ii-1),:)=x_temp;
end

x0=[x_new;xt;x0];


options=optimoptions('surrogateopt','MaxFunctionEvaluations',2.*10^3,'useparallel',true,'InitialPoints',x0,'PlotFcn','surrogateoptplot');

[Vac_Est,f_est]=surrogateopt(@(x)Objective_Function_Cases(x,Year_2005_2013,Population,Reported_Cases,Vaccine_Coverage,p_inf_11_15,p_inf_16_23,p_vac_11_15,p_vac_16_23),LBC,UBC,[1:length(LBC)],A,b,[],[],options);

Vac_Est=Vac_Est(:);

Vaccinated_One_Dose_Cases_Age_11_to_15=Vac_Est(1:length(Year_2005_2013));
Vaccinated_Two_Dose_Cases_Age_11_to_15=Vac_Est(length(Year_2005_2013)+[1:length(Year_2005_2013)]);
Vaccinated_One_Dose_Cases_Age_16_to_23=Vac_Est(2.*length(Year_2005_2013)+[1:length(Year_2005_2013)]);
Vaccinated_Two_Dose_Cases_Age_16_to_23=Vac_Est(3.*length(Year_2005_2013)+[1:length(Year_2005_2013)]);

Unvaccinated_Cases_Age_11_to_15=Reported_Cases.National.Age_11_15-Vaccinated_One_Dose_Cases_Age_11_to_15-Vaccinated_Two_Dose_Cases_Age_11_to_15;
Unvaccinated_Cases_Age_16_to_23=Reported_Cases.National.Age_16_23-Vaccinated_One_Dose_Cases_Age_16_to_23-Vaccinated_Two_Dose_Cases_Age_16_to_23;

% Create Table
Total_Reported_Cases_Age_11_15=Reported_Cases.National.Age_11_15;
Total_Reported_Cases_Age_16_23=Reported_Cases.National.Age_16_23;
T_Cases=table(Year_2005_2013,Total_Reported_Cases_Age_11_15,Unvaccinated_Cases_Age_11_to_15,Vaccinated_One_Dose_Cases_Age_11_to_15,Vaccinated_Two_Dose_Cases_Age_11_to_15,Total_Reported_Cases_Age_16_23,Unvaccinated_Cases_Age_16_to_23,Vaccinated_One_Dose_Cases_Age_16_to_23,Vaccinated_Two_Dose_Cases_Age_16_to_23);

writetable(T_Cases,'Estimated_Cases_Among_Vaccinatied.xlsx');
