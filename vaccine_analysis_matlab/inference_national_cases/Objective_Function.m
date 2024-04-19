function J=Objective_Function(x,Year_ABC,Population,Reported_Cases,Vaccine_Coverage,Reported_Cases_2014_2021,Population_2014_2021)

% Scale the parameters
p_inf_11_15=x(1)./10^5;
p_inf_16_23=x(2)./10^5;

p_vac_11_15=x(3)./10^5;
p_vac_16_23=x(4)./10^5;

% Likelihood arrays
L_11_15_Unvaccinated=zeros(length(Year_ABC),1);
L_16_23_Unvaccinated=zeros(length(Year_ABC),1);

L_11_15_Vaccinated_One_Dose=zeros(length(Year_ABC),1);
L_11_15_Vaccinated_Two_Dose=zeros(length(Year_ABC),1);
L_16_23_Vaccinated_One_Dose=zeros(length(Year_ABC),1);
L_16_23_Vaccinated_Two_Dose=zeros(length(Year_ABC),1);

L_11_15_National=zeros(length(Year_ABC),1);
L_16_23_National=zeros(length(Year_ABC),1);

L_11_15_Vaccination_National=zeros(length(Year_ABC),1);
L_16_23_Vaccination_National=zeros(length(Year_ABC),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% National Data 2014 to 2021 for 11 to 15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
L_11_15_2014to2021=log(binopdf(Reported_Cases_2014_2021.Unvaccinated.Age_11_15,Population_2014_2021.National.Unvaccinated.Age_11_15,p_inf_11_15));
L_11_15_2014to2021=L_11_15_2014to2021+log(binopdf(Reported_Cases_2014_2021.Vaccinated_One_Dose.Age_11_15,Population_2014_2021.National.Vaccinated_One_Dose.Age_11_15,p_vac_11_15));
L_11_15_2014to2021=L_11_15_2014to2021+log(binopdf(Reported_Cases_2014_2021.Vaccinated_Two_Dose.Age_11_15,Population_2014_2021.National.Vaccinated_Two_Dose.Age_11_15,p_vac_11_15));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% National Data 2014 to 2021 for 11 to 15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
L_16_23_2014to2021=log(binopdf(Reported_Cases_2014_2021.Unvaccinated.Age_16_23,Population_2014_2021.National.Unvaccinated.Age_16_23,p_inf_16_23));
L_16_23_2014to2021=L_16_23_2014to2021+log(binopdf(Reported_Cases_2014_2021.Vaccinated_One_Dose.Age_16_23,Population_2014_2021.National.Vaccinated_One_Dose.Age_16_23,p_vac_16_23));
L_16_23_2014to2021=L_16_23_2014to2021+log(binopdf(Reported_Cases_2014_2021.Vaccinated_Two_Dose.Age_16_23,Population_2014_2021.National.Vaccinated_Two_Dose.Age_16_23,p_vac_16_23));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
% Data based on additional report
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555

Count_Additional_Data_11_15=0;
Count_Additional_Data_16_23=0;
Count_National_11_15=0;
Count_National_16_23=0;
for yy=1:length(Reported_Cases.Additional.Years)
    tf=Reported_Cases.Additional.Years(yy)==Year_ABC;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Add the expected number of cases for 11_to 15
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % ONE DOSE
    Count_Additional_Data_11_15=Count_Additional_Data_11_15+binostat(round(Population.Additional.Vaccinated_2005_to_2008.Fraction_National.*Population.National.Age_11_15(tf).*Vaccine_Coverage.One_Dose.Age_11_15(tf)),p_vac_11_15);
    % TWO DOSE
    Count_Additional_Data_11_15=Count_Additional_Data_11_15+binostat(round(Population.Additional.Vaccinated_2005_to_2008.Fraction_National.*Population.National.Age_11_15(tf).*Vaccine_Coverage.Two_Dose.Age_11_15(tf)),p_vac_11_15);

    % ONE DOSE
    Count_National_11_15=Count_National_11_15+binostat(round(Population.National.Age_11_15(tf).*Vaccine_Coverage.One_Dose.Age_11_15(tf)),p_vac_11_15);
    % TWO DOSE
    Count_National_11_15=Count_National_11_15+binostat(round(Population.National.Age_11_15(tf).*Vaccine_Coverage.Two_Dose.Age_11_15(tf)),p_vac_11_15);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Add expected cases 16 to 23
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % ONE DOSE
    Count_Additional_Data_16_23=Count_Additional_Data_16_23+binostat(round(Population.Additional.Vaccinated_2005_to_2008.Fraction_National.*Population.National.Age_16_23(tf).*Vaccine_Coverage.One_Dose.Age_16_23(tf)),p_vac_16_23);
    % TWO DOSE
    Count_Additional_Data_16_23=Count_Additional_Data_16_23+binostat(round(Population.Additional.Vaccinated_2005_to_2008.Fraction_National.*Population.National.Age_16_23(tf).*Vaccine_Coverage.Two_Dose.Age_16_23(tf)),p_vac_16_23);

    % ONE DOSE
    Count_National_16_23=Count_National_16_23+binostat(round(Population.National.Age_16_23(tf).*Vaccine_Coverage.One_Dose.Age_16_23(tf)),p_vac_16_23);
    % TWO DOSE
    Count_National_16_23=Count_National_16_23+binostat(round(Population.National.Age_16_23(tf).*Vaccine_Coverage.Two_Dose.Age_16_23(tf)),p_vac_16_23);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11-15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Likelihood to reflect that estimated number of cases is simialr
L_11_15_Additional=log(poisspdf(Reported_Cases.Additional.Vaccinated_2005_to_2008.Age_11_15,Count_Additional_Data_11_15));
% The report suggest that this number is under reported so opting to use
% the cdf and the reported numebr as a lower bound
L_11_15_Additional=L_11_15_Additional+log(1-poisscdf(Reported_Cases.Additional.Vaccinated_2005_to_2008.Age_11_15,Count_National_11_15));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 16-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Likelihood to reflect that estimated number of cases is simialr
L_16_23_Additional=log(poisspdf(Reported_Cases.Additional.Vaccinated_2005_to_2008.Age_16_23,Count_Additional_Data_16_23));
% The report suggest that this number is under reported so opting to use
% the cdf and the reported numebr as a lower bound
L_16_23_Additional=L_16_23_Additional+log(1-poisscdf(Reported_Cases.Additional.Vaccinated_2005_to_2008.Age_16_23,Count_National_16_23));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
% ABC reported data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555

for yy=1:length(Year_ABC)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Ages 11 to 15
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    L_11_15_Unvaccinated(yy)=log(binopdf(Reported_Cases.ABC.Unvaccinated.Age_11_15(yy),Population.ABC.Unvaccinated.Age_11_15(yy),p_inf_11_15));
    L_11_15_Vaccinated_One_Dose(yy)=log(binopdf(Reported_Cases.ABC.Vaccinated_One_Dose.Age_11_15(yy),Population.ABC.Vaccinated_One_Dose.Age_11_15(yy),p_vac_11_15));
    L_11_15_Vaccinated_Two_Dose(yy)=log(binopdf(Reported_Cases.ABC.Vaccinated_Two_Dose.Age_11_15(yy),Population.ABC.Vaccinated_Two_Dose.Age_11_15(yy),p_vac_11_15));

    p_11_15=(Vaccine_Coverage.One_Dose.Age_11_15(yy)+Vaccine_Coverage.Two_Dose.Age_11_15(yy)).*p_vac_11_15+(1-Vaccine_Coverage.One_Dose.Age_11_15(yy)-Vaccine_Coverage.Two_Dose.Age_11_15(yy)).*p_inf_11_15;
    L_11_15_National(yy)=log(binopdf(Reported_Cases.National.Age_11_15(yy),Population.National.Age_11_15(yy),p_11_15));
    
    % Specify that the number of vaccinated cases reported by ABC is a
    % lower bound for the national number
    temp_v_one_dose_case=binostat(round(Population.National.Age_11_15(yy).*Vaccine_Coverage.One_Dose.Age_11_15(yy)),p_vac_11_15);
    temp_v_two_dose_case=binostat(round(Population.National.Age_11_15(yy).*Vaccine_Coverage.Two_Dose.Age_11_15(yy)),p_vac_11_15);
    L_11_15_Vaccination_National(yy)=log(1-poisscdf(Reported_Cases.ABC.Vaccinated_One_Dose.Age_11_15(yy),temp_v_one_dose_case))+log(1-poisscdf(Reported_Cases.ABC.Vaccinated_Two_Dose.Age_11_15(yy),temp_v_two_dose_case));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % Ages 16 to 23
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    L_16_23_Unvaccinated(yy)=log(binopdf(Reported_Cases.ABC.Unvaccinated.Age_16_23(yy),Population.ABC.Unvaccinated.Age_16_23(yy),p_inf_16_23));
    L_16_23_Vaccinated_One_Dose(yy)=log(binopdf(Reported_Cases.ABC.Vaccinated_One_Dose.Age_16_23(yy),Population.ABC.Vaccinated_One_Dose.Age_16_23(yy),p_vac_16_23));
    L_16_23_Vaccinated_Two_Dose(yy)=log(binopdf(Reported_Cases.ABC.Vaccinated_Two_Dose.Age_16_23(yy),Population.ABC.Vaccinated_Two_Dose.Age_16_23(yy),p_vac_16_23));

    p_16_23=(Vaccine_Coverage.One_Dose.Age_16_23(yy)+Vaccine_Coverage.Two_Dose.Age_16_23(yy)).*p_vac_16_23+(1-Vaccine_Coverage.One_Dose.Age_16_23(yy)-Vaccine_Coverage.Two_Dose.Age_16_23(yy)).*p_inf_16_23;
    L_16_23_National(yy)=log(binopdf(Reported_Cases.National.Age_16_23(yy),Population.National.Age_16_23(yy),p_16_23));

    % Specify that the total number of vaccinated cases reported by ABC is a
    % lower bound for the national number
    temp_v_one_dose_case=binostat(round(Population.National.Age_16_23(yy).*Vaccine_Coverage.One_Dose.Age_16_23(yy)),p_vac_16_23);
    temp_v_two_dose_case=binostat(round(Population.National.Age_16_23(yy).*Vaccine_Coverage.Two_Dose.Age_16_23(yy)),p_vac_16_23);
    L_16_23_Vaccination_National(yy)=log(1-poisscdf(Reported_Cases.ABC.Vaccinated_One_Dose.Age_16_23(yy),temp_v_one_dose_case))+log(1-poisscdf(Reported_Cases.ABC.Vaccinated_Two_Dose.Age_16_23(yy),temp_v_two_dose_case));
end

% Objective function aimed to minimize
J=-sum(L_11_15_2014to2021(:))-sum(L_16_23_2014to2021(:))-L_11_15_Additional-L_16_23_Additional-sum(L_11_15_Unvaccinated)-sum(L_11_15_Vaccinated_One_Dose)-sum(L_11_15_Vaccinated_Two_Dose)-sum(L_11_15_National)-sum(L_11_15_Vaccination_National)-sum(L_16_23_Unvaccinated)-sum(L_16_23_Vaccinated_One_Dose)-sum(L_16_23_Vaccinated_Two_Dose)-sum(L_16_23_National)-sum(L_16_23_Vaccination_National);

end