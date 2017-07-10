% CONST_CALL_SWS_95(temperature, Sal, Press):
% Alex gagnon September 2009
% temeprature and salinity dependent constants of carbonate
% system on the seawter pH scale.  Unless otherwise noted the constants 
% are described in Millero, Mar Chem 1995.  I am attemptign to use the same
% set of constants as were evaluated with the GLODAP data by Lee et al.,
% GRL 2000; Lamb et al., Deep-Sea Res II, 2002.
% temperature = temperature in degrees C
% Sal = salinity in paractical salinity units
% Press = pressure in bar (use p=0 for 1atm)
%
%note: when calculating pCO2 at high DIC it may be adventageous
%to correct K_2 as a function of DIC (equation 21 in Mllero review)
%
%eliminated pressure dependence in this set of constants
%
%Alex Gagnon Jan 2009: major update using constants 
%suggested by Millero, Chem. Rev. (2007). an adaptation of my previous
%function const_call.m which was itself inspired  by Zeebe's
%TEMPT.m

function dummy = const_call_SWS_95_calcite_f(temperature, Sal, Press)

global S K_1 K_W K_2 K_B I K_sp K_SO4 K_F B_T F_T SO4_T Ca_T f_coeff_CO2_1atm;
global K_sp_calcite f_h;
global k_plus1 k_minus1 k_plus4 k_minus4 

TK = temperature+273.15;  %kelvin temp
I = 19.924 * Sal /(1000- 1.0049*Sal); %molal ionic strength
R = 83.131; %Gas constant (units? cm^3 bar mol-1 deg-1)

%Conservative seawater components as a function of salinity
B_T = 0.000416*Sal/35; %total boron
F_T = 7e-5*Sal/35; %total HF and F-
%SO4_T = .0293*Sal/35 %total sulfate species from Millero 1995
SO4_T = .02824*Sal/35; %total sulfate species from DOE 1994, and Dickson/Millero's refitting of Mehrbach
Ca_T = .01028*Sal/35;

%1st dissociation constant for seawater as a function of T and S
%Seawater Scale
%Merhrbach et al 1973 as re-fit by Dickson and Millero 1987
tempA = -.0118*Sal + 0.000116*Sal^2;
pK1_star = -62.008 + tempA + 3670.7/TK + 9.7944*log(TK);
K_1_star = 10^(-1*pK1_star);
K_1 = K_1_star;

%2nd dissociation contant for seawater as a function of T and S
%Seawater Scale
%Merhrbach et al 1973 as re-fit by Dickson and Millero 1987
tempA = -.0184*Sal + 0.000118*Sal^2;
pK2_star = 4.777 + tempA + 1394.7/TK;
K_2_star = 10^(-1*pK2_star);
K_2=K_2_star;

%dissociation contant for H2O as reported in Millero 1995
%Seawater Scale
temp1 = 148.9802-13847.26/TK-23.6521*log(TK);
temp2 = (118.67/TK-5.977+1.0495*log(TK))*sqrt(Sal)-0.01615*Sal;
K_W_star = exp(temp1 + temp2);
pKW_star =-1*log10(K_W_star);
K_W=K_W_star;

%Boric Acid
%as reported in Millero 1995 
%from Dickson Deep-Sea Res 1990 (Roy 1993 agrees)
%both are in synthetic seawater without F-
temp1 = (-8966.90-2890.51*sqrt(Sal)-77.942*Sal+1.726*(Sal^(3/2))-0.0993*(Sal^2))/TK;
temp2 = 148.0248+137.1942*sqrt(Sal)+1.62247*Sal;
temp3 = -log(TK)*(24.4344+25.085*sqrt(Sal)+0.2474*Sal)+0.053105*sqrt(Sal)*TK;
temp1+temp2+temp3;
K_B_star = exp(temp1+temp2+temp3);
pKB_star =-1*log10(K_B_star);
K_B = K_B_star;

%Bisulfate dissolution : K_SO4
%K_SO4 = H+ x SO42- / HSO4-
%free scale (of necessity, it doesnt exist in SWS and total scale)
%from Dickson 1990
temp1 = -4276.1/TK + 141.328 - 23.093* log(TK);
temp2 = sqrt(I)*(-13856/TK + 324.57 -47.986*log(TK));
temp3 = I*(35474/TK - 771.54 + 114.723 *log(TK));
temp4 = -2698/TK*(I^(1.5)) + 1776/TK*I^2 + log(1-.001005*Sal);
temp1+temp2+temp3+temp4;
K_SO4 = exp(temp1+temp2+temp3+temp4);

%Hydrogen Fluoride: K_F
%K_F = H+ x F- / HF
%I think this is on the free scale and in mol k solution-1
%I started with the expression in DOE and removed the term which converst
%from the free to total scale.  HF not in the alkalinity expression for
%the seawater scale
temp1 = 1590.2/TK - 12.641 + 1.525 * sqrt(I);
temp2 = log(1-.001005*Sal);
K_F = exp(temp1 + temp2);

%CaCO3 solubility ARAGONITE
%aragonite at 1 atm (Mucci, 1983)
%calculated in natural seawater, using synthetic aragonite,
%carbonate ion was calculated from pH_NBS and total ALK,
%using the originial equilibrium constants of Mehrbach 1973 
temp1 = -171.945-0.077993*TK+2903.293/TK+(71.595*log10(TK));
temp2 = (-0.068393+0.0017276*TK+88.135/TK)*sqrt(Sal);
temp3 = -0.10018*Sal+(0.0059415*(Sal^1.5));
temp1+temp2+temp3;
K_sp = 10^(temp1+temp2+temp3);

%CaCO3 solubility CALCITE
%calcite at 1 atm (Mucci, 1983)
%calculated in natural seawater, using synthetic aragonite,
%carbonate ion was calculated from pH_NBS and total ALK,
%using the originial equilibrium constants of Mehrbach 1973 
%First calculate the thermodynamic solubility constant for calcite
%using the data from Plumer and Busenberg (1982)
temp1 = -171.9065-0.077993*TK+2839.319/TK+(71.595*log10(TK));
temp2 = (-0.77712+0.0028426*TK+178.34/TK)*sqrt(Sal);
temp3 = -0.07711*Sal+(0.0041249*(Sal^1.5));
temp1+temp2+temp3;
K_sp_calcite = 10^(temp1+temp2+temp3);

%solubility of CO2 (Weiss 1974) (pressure effect?)
%does this need to be recalculated with new Millero constants?
temp1 = 9345.17/TK-60.2409+23.3585*log(TK/100);
temp2 = Sal*(0.023517-0.00023656*TK+0.0047036*((TK/100)^2));
S = exp(temp1+temp2);
%log10(S)*-1

%CO2 fugacity coeffcient at 1 atm
BViral = 1e-6*(-1636.75 + 12.0408*TK - 3.27957e-2*(TK^2) + 3.16528e-5*(TK^3)); %first viral coeffceint for CO2
dViral = 1e-6*(57.7-0.118*TK); %mixing cross viral coeffceint for CO2 - air 
press_pascal = 101325; % 1atm
f_coeff_CO2_1atm = exp( (BViral + 2*dViral)*press_pascal/8.314/TK );

%%%%%%%
%proton activity coeffcient at 1 atm
%fit surfcae to data of Mehrbach et al. L&O 1973
%First input experimental data
s=[19.19
25.25
26.75
27.01
31.49
34.76
34.78
34.95
35.1
35.23
42.26
42.67
42.85
42.88];
t=[25
25
35
2
25
13
2
35
25
25
25
35
13
2];
f=[0.677
0.676
0.625
0.772
0.684
0.753
0.805
0.622
0.693
0.692
0.695
0.629
0.788
0.833];
%then fit surface to this data, where the surface is formed from a regualr
%array of points
Salinity_i = [19:43];
Temperature_i = [2:35];
[SI,TI] = meshgrid(Salinity_i,Temperature_i);
FI = griddata(s,t,f,SI,TI);

% %for diagnostic purposes, you can plot the data:
% mesh(SI,TI,FI), hold
% plot3(s,t,f,'o'), hold off

%now interpolate the data form the girdded function to get the activity at
%a particular value

f_h  = interp2(SI,TI,FI,Sal,temperature);

%%%%

dummy=1;

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%I havent critically reviewed constants below this point yet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Hydrogen Fluoride: K_F
%uses the total sulfate and K_SO4
temp1 = 1590.2/TK - 12.641 + 1.525 * sqrt(I);
temp2 = log(1-.001005*Sal) + log(1+SO4_T/K_SO4);
K_F = exp(temp1 + temp2);

%CaCO3 solubility (Mucci)
%aragonite at 1 atm
temp1 = -171.945-0.077993*TK+2903.293/TK+(71.595*log10(TK));
temp2 = (-0.068393+0.0017276*TK+88.135/TK)*sqrt(Sal);
temp3 = -0.10018*Sal+(0.0059415*(Sal^1.5));
K_sp = 10^(temp1+temp2+temp3)
%pressure effect based on summary in Zeebe
tempA0 = 46.00;
tempA1 = .5304;
tempA2 = 0;
Delta_V = -1*tempA0 + tempA1/1e3*temperature + tempA2*(temperature^2); %molal volume change
tempB0 = 11.76;
tempB1 = .3692;
tempB2 = 0;
Delta_k = -1*tempB0/1e3 + tempB1/1e3*temperature + tempB2*(temperature^2); %molal compresability
LN_K_sp_P = log(K_sp)- Delta_V/R/TK*Press + 0.5*Delta_k/R/TK*Press^2; %T in Kelvin
K_sp_P = exp(LN_K_sp_P)
pK_sp_P = -1*log10(K_sp_P)
K_sp=K_sp_P;


%kinetics of hydration/hydroxylation
k_plus1 = exp(1246.98-6.19e4/TK-183.0*log(TK));
%rate enhancement
%k_plus1=k_plus1*1e6;
k_minus1 = k_plus1/K_1;

k_plus4 = 4.7e7*exp(-23.2*1e3/8.314/TK);
%rate enhancement
%k_plus4=k_plus4*1e9;
k_minus4 = k_plus4*K_W/K_1;

dummy=1;
return;

k_plus1 = exp(1246.98-6.19e4/TK-183.0*log(TK));
%rate enhancement
%k_plus1=k_plus1*1e6;

%DOE Recommended from Roy et al.
temp1 = 2.83655-(2307.1266/TK)-(1.5529413*log(TK));
temp3 = -((0.207608410+4.0484/TK)*sqrt(Sal));
temp2 = (0.0846834*Sal)-(0.00654208*(Sal^(3/2)))+(log(1-0.001005*Sal));
K_1 = exp(temp1+temp2+temp3);
k_minus1 = k_plus1/K_1;


temp1 = 148.96502-13847.26/TK-23.6521*log(TK);
temp2 = (118.67/TK-5.977+1.0495*log(TK))*sqrt(Sal)-0.01615*Sal;
K_W = exp(temp1 + temp2);

k_plus4 = 4.7e7*exp(-23.2*1e3/8.314/TK);
%rate enhancement
%k_plus4=k_plus4*1e9;

k_minus4 = k_plus4*K_W/K_1;

temp1 = 9345.17/TK-60.2409+23.3585*log(TK/100);
temp2 = Sal*(0.023517-0.00023656*TK+0.0047036*((TK/100)^2));
S = exp(temp1+temp2);

temp1 = -9.226508-3351.6106/TK-0.2005743*log(TK);
temp2 = -1*(0.106901773+23.9722/TK)*sqrt(Sal);
temp3 = 0.1130822*Sal-0.00846934*(Sal^(3/2))+log(1-0.001005*Sal);
K_2 = exp(temp1+temp2+temp3);

temp1 = (-8966.90-2890.53*sqrt(Sal)-77.942*Sal+1.728*(Sal^(3/2))-0.0996*(Sal^2))/TK;
temp2 = 148.0248+137.1942*sqrt(Sal)+1.62142*Sal;
temp3 = -log(TK)*(24.4344+25.085*sqrt(Sal)+0.2474*Sal)+0.053105*sqrt(Sal)*TK;
K_B = exp(temp1+temp2+temp3);

%aragonite at 1 atm
temp1 = -171.945-0.077993*TK+2903.293/TK+(71.595*log10(TK));
temp2 = (-0.068393+0.0017276*TK+88.135/TK)*sqrt(Sal);
temp3 = -0.10018*Sal+(0.0059415*(Sal^1.5));
K_sp = 10^(temp1+temp2+temp3);

dummy =1;