function J=Objective_Function_Cases(x,Year,Population,Reported_Cases,Vaccine_Coverage,p_inf_11_15,p_inf_16_23,p_vac_11_15,p_vac_16_23)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Estimates to evalaute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
Vac_One_Dose_Case_Age_11_15=x(1:length(Year));
Vac_Two_Dose_Case_Age_11_15=x(length(Year)+[1:length(Year)]);

Vac_One_Dose_Case_Age_16_23=x(2.*length(Year)+[1:length(Year)]);
Vac_Two_Dose_Case_Age_16_23=x(3.*length(Year)+[1:length(Year)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Requirment from the 2005 to 2008 Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cases_2005_2008_11_15=sum(Vac_One_Dose_Case_Age_11_15(ismember(Year,Reported_Cases.Additional.Years))+Vac_Two_Dose_Case_Age_11_15(ismember(Year,Reported_Cases.Additional.Years)));
Cases_2005_2008_16_23=sum(Vac_One_Dose_Case_Age_16_23(ismember(Year,Reported_Cases.Additional.Years))+Vac_Two_Dose_Case_Age_16_23(ismember(Year,Reported_Cases.Additional.Years)));
L_2005_2008=log(1-poisscdf(Reported_Cases.Additional.Vaccinated_2005_to_2008.Age_11_15,Cases_2005_2008_11_15))+log(1-poisscdf(Reported_Cases.Additional.Vaccinated_2005_to_2008.Age_16_23,Cases_2005_2008_16_23));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Age 11 to 15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_11_15=zeros(length(Year),2);

for yy=1:length(Year)
    % Unvaccinated
    L_11_15(yy,1)=log(binopdf(Reported_Cases.National.Age_11_15(yy)-Vac_One_Dose_Case_Age_11_15(yy)-Vac_Two_Dose_Case_Age_11_15(yy),Population.National.Age_11_15(yy)-round(Population.National.Age_11_15(yy).*Vaccine_Coverage.One_Dose.Age_11_15(yy))-round(Population.National.Age_11_15(yy).*Vaccine_Coverage.Two_Dose.Age_11_15(yy)),p_inf_11_15));
    % ONE DOSE
    L_11_15(yy,1)=L_11_15(yy,1)+log(binopdf(Vac_One_Dose_Case_Age_11_15(yy),round(Population.National.Age_11_15(yy).*Vaccine_Coverage.One_Dose.Age_11_15(yy)),p_vac_11_15));
    % TWO DOSE
    L_11_15(yy,1)=L_11_15(yy,1)+log(binopdf(Vac_Two_Dose_Case_Age_11_15(yy),round(Population.National.Age_11_15(yy).*Vaccine_Coverage.Two_Dose.Age_11_15(yy)),p_vac_11_15));
    % ONE DOSE LB
    if(Vac_One_Dose_Case_Age_11_15(yy)>0)
        L_11_15(yy,2)=log(1-poisscdf(Reported_Cases.ABC.Vaccinated_One_Dose.Age_11_15(yy),Vac_One_Dose_Case_Age_11_15(yy)));
    else
        L_11_15(yy,2)=round(Population.National.Age_11_15(yy).*Vaccine_Coverage.One_Dose.Age_11_15(yy)).*log(1-p_vac_11_15);
    end
    % TWO DOSE LB
    if(Vac_Two_Dose_Case_Age_11_15(yy)>0)
        L_11_15(yy,2)=L_11_15(yy,2)+log(1-poisscdf(Reported_Cases.ABC.Vaccinated_Two_Dose.Age_11_15(yy),Vac_Two_Dose_Case_Age_11_15(yy)));
    else
        L_11_15(yy,2)=L_11_15(yy,2)+round(Population.National.Age_11_15(yy).*Vaccine_Coverage.Two_Dose.Age_11_15(yy)).*log(1-p_vac_11_15);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Age 16 to 23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_16_23=zeros(length(Year),2);

for yy=1:length(Year)
    % Unvaccinated
    L_16_23(yy,1)=log(binopdf(Reported_Cases.National.Age_16_23(yy)-Vac_One_Dose_Case_Age_16_23(yy)-Vac_Two_Dose_Case_Age_16_23(yy),Population.National.Age_16_23(yy)-round(Population.National.Age_16_23(yy).*Vaccine_Coverage.One_Dose.Age_16_23(yy))-round(Population.National.Age_16_23(yy).*Vaccine_Coverage.Two_Dose.Age_16_23(yy)),p_inf_16_23));
    % ONE DOSE
    L_16_23(yy,1)=L_16_23(yy,1)+log(binopdf(Vac_One_Dose_Case_Age_16_23(yy),round(Population.National.Age_16_23(yy).*Vaccine_Coverage.One_Dose.Age_16_23(yy)),p_vac_16_23));
    % TWO DOSE
    L_16_23(yy,1)=L_16_23(yy,1)+log(binopdf(Vac_Two_Dose_Case_Age_16_23(yy),round(Population.National.Age_16_23(yy).*Vaccine_Coverage.Two_Dose.Age_16_23(yy)),p_vac_16_23));
    % ONE DOSE LB
    if(Vac_One_Dose_Case_Age_16_23(yy)>0)
        L_16_23(yy,2)=log(1-poisscdf(Reported_Cases.ABC.Vaccinated_One_Dose.Age_16_23(yy),Vac_One_Dose_Case_Age_16_23(yy)));
    else
        L_16_23(yy,2)=round(Population.National.Age_16_23(yy).*Vaccine_Coverage.One_Dose.Age_16_23(yy)).*log(1-p_vac_16_23);
    end
    % TWO DOSE LB
    if(Vac_Two_Dose_Case_Age_16_23(yy)>0)
        L_16_23(yy,2)=L_16_23(yy,2)+log(1-poisscdf(Reported_Cases.ABC.Vaccinated_Two_Dose.Age_16_23(yy),Vac_Two_Dose_Case_Age_16_23(yy)));
    else
        L_16_23(yy,2)=L_16_23(yy,2)+round(Population.National.Age_16_23(yy).*Vaccine_Coverage.Two_Dose.Age_16_23(yy)).*log(1-p_vac_16_23);
    end
end

J=-sum(L_11_15(:))-sum(L_16_23(:))-L_2005_2008;

end