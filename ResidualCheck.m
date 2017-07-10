%Residual Check Script

global S K_1 K_W K_2 K_B I K_sp K_SO4 K_F B_T F_T SO4_T Ca_T f_coeff_CO2_1atm;
global K_sp_calcite f_h;
global k_plus1 k_minus1 k_plus4 k_minus4 

%%%%%%%%%%%%%%%%%%%%
%USER SET PARAMETERS
%%%%%%%%%%%%%%%%%%%%

%Seawater(culture) conditions
Temperature = 25; %deg C
Pressure = 0; %bar
Salinity =  35;
ALK_SW = 2400e-6; %mol equiv/kg --> when changed to a fucntion these will be an input
DIC_SW = 2000e-6; %mol/kg --> when changed to a fucntion these will be an input
Ca_SW= 10.2e-3; %mol/kg --> when changed to a fucntion these will be an input
F_kzrho = 1e-4;
D_kzrho = 1e-3;

const_call_SWS_95_calcite_f(Temperature, Salinity, Pressure);


%CO2_Cell
pCO2_Cell = 5000e-6; %ppm
CO2_Cell = pCO2_Cell*S;

%Slope of pH rule
slope = 1.6e-5 %this function is defined at end of script and is set by user - AMG changed this value to be consisent with slopes from literature data

TA = ALK_SW;
DIC = DIC_SW
    
%solve for carboante system of surrouding seawater
%(probably good to re-check this, esp the boron part, it has been a while since I wrote it)
fac = conv([1 0],conv([1 K_B],conv([1 K_1 (K_1*K_2)],[(1/K_2) (K_1/K_2) K_1])));  %factor multiplied to general TA equation
temp1= conv([0 -TA],fac);  %TA =
temp2= conv([0 0 (K_1*DIC) 0],deconv(fac, [1 K_1 (K_1*K_2)]));   %+HCO3
temp6= conv([0 0 0 (2*K_1*DIC)],deconv(fac,[(1/K_2) (K_1/K_2) K_1]));  %+2CO3
temp3= conv([0 0 K_W],deconv(fac,[1 0])); % +OH
temp4= conv([-1 0],fac);  % -H
temp5= conv([0 0 (K_B*B_T)],deconv(fac,[1 K_B]));  %+Anionic form of Buffer

poly_to_solve = temp1+temp2+temp3+temp4+temp5+temp6;

%only one positive root
h = max(roots(poly_to_solve)); %this is [H+] in seawater
alpha0 = 1./(1 + (K_1./h) + ((K_1.*K_2)./(h.*h)));
alpha2 = K_1.*K_2 ./ (h.*h + h.*K_1 + K_1.*K_2);
CO2_SW = alpha0.* DIC;    
pCO2_SW = CO2_SW./S;
    
%calculate pH_CF from Seawater Chemistry, using Rule
delta_h = slope*pCO2_SW; % where delta_h is defined as h - h_cf (i.e. [H+]sw - [H+]cf)
h_cf = -delta_h + h; % AMG changed the sign for delta_h to be consistent with the definition in comments above
pH_CF = -log10(h_cf);
    
DIC_CF = -0.01:0.00001:0.01;

hh = 10^(-1*pH_CF);
oh = K_W/hh;
alpha_0 = 1/(1 +K_1/hh +K_1*K_2/hh/hh);
alpha_1 = 1/(1 +hh/K_1 + K_2/hh);
alpha_2 = 1/(1 +hh/K_2 + hh*hh/K_1/K_2);
CO2_CF = DIC_CF.*alpha_0;
HCO3_CF= DIC_CF.*alpha_1;
CO3_CF = DIC_CF.*alpha_2;
BOH4_CF = B_T/(1+ hh/K_B);
ALK_CF = HCO3_CF + 2*CO3_CF + BOH4_CF + oh - hh;
    
%2. Use ALK_CF and dALK/dt to get gamma
% gamma is defined as F/P, the ratio of pumping to ppt
gammaFP = (0.5 - (ALK_CF - ALK_SW)./2./F_kzrho).^-1;
%gammaFP = 2/((ALK_SW-ALK_CF)/F_kzrho + 1);
    
    
%3. Plug guess for DIC_CF into dDIC/dt equation together with calculated CO2_CF and gamma
%this residual will be zero when the DIC_CF guess is correct
residual = D_kzrho.*CO2_Cell - D_kzrho.*CO2_CF + DIC_SW - DIC_CF - F_kzrho./gammaFP;
    
subplot(1,2,1)
plot(DIC_CF, residual)
xlabel('DIC_CF')
ylabel('residual')
    
subplot(1,2,2)
plot(residual, gammaFP)
xlabel('residual')
ylabel('gammaFP')
    
  
   