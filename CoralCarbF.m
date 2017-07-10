%CoralCarbF: Solve for complete state of coral calcifying fluid from the followign parameters:
%F/kzrho and D/kzrho and some empirical realtionship that porvides pH_CF,
%for example pH_CF as a function of seawater chemistry. It is also
%necessary to know ALK_SW, DIC_SW, CO2_cell. If you want to get P, kzrho, or other
%individual fluxes and not just the reduced dynamical paramters then you
%also need to supply a function for precipitation rate as a function of
%omega.
%Alex Gagnon 23 June 2017
%AMG 23 June 2017 changes to lines 42, 50 and 95
%AMG 3 July 2017 edits for addition of C isotopes

function [DIC_CF, pH_CF, ALK_CF, CO3_CF, CO2_CF, P_kzrho, D_kzrhoCO2, pCO2_Cell, pCO2_SW, R1312C_CF] = CoralCarbF(ALK_SW, DIC_SW, F_kzrho, D_kzrho, empirical_slope)

format long %sets the formatting of the output data to be long, i.e. more sig figs are displayed in the command window
tic %start stopwatch timer

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

%Below terms can also be added as inputs to this function
Ca_SW= 10.2e-3; %mol/kg --> when changed to a fucntion these will be an input
pCO2_Cell = 5000e-6; %ppm

%%%%%%%%%%%%%
%CALCULATIONS
%%%%%%%%%%%%%

%Calculation1: load carbonate system constants into global paramters
const_call_SWS_95_calcite_f(Temperature, Salinity, Pressure);
CO2_Cell = pCO2_Cell.*S; %assume same salinity and temperature in "cell" as for seawater to calclaute [CO2]_cell where S here is solubility % AMG changed this equation from 'CO2_sw = pCO2_Cell/S' to 'CO2_sw = pCO2_Cell*S'

%Calculation2: get pH of Calcifying fluis using pHCFRule (function set by user at bottom of script)
[pH_CF pCO2_SW] = pHCFRule(ALK_SW,DIC_SW,empirical_slope);

%Calculation3: solve system of model equations for main model parameters
%using pH_CF. This is done by finding zero of system using itterative guesses for DIC
zerofun = @(x) SolveSSCoral(x, pH_CF, F_kzrho, D_kzrho, ALK_SW, DIC_SW, CO2_Cell);
DIC_CF = fzero(zerofun, 1000e-6);
% below are alternative lines for the fzero calculation including more
% display options for code transparency
%options = optimset('Display','iter'); % show iterations
%DIC_CF = fzero(zerofun, 200e-6, options)

%AMG edited the output of this function to spit out calcifying fluid parameters including
% DIC, ALK, CO3, CO2, as well as other nondimensional parameters gammaFP, etc.
% These lines solve for all of the carbonate system parameters with DIC_CF and H+
hh = 10.^(-1.*pH_CF);  % [H+ cf]
oh = K_W./hh; % [OH - cf]
alpha_0 = 1./(1 +K_1./hh +K_1.*K_2./hh./hh); % alpha0 (and alpha1, alpha2 below
alpha_1 = 1./(1 +hh./K_1 + K_2./hh);
alpha_2 = 1./(1 +hh./K_2 + hh.*hh./K_1./K_2);
CO2_CF = DIC_CF.*alpha_0; %CO2_cf
HCO3_CF= DIC_CF.*alpha_1; %HCO3_cf
CO3_CF = DIC_CF.*alpha_2; %CO3_cf
BOH4_CF = B_T./(1+ hh./K_B); %Borate
ALK_CF = HCO3_CF + 2.*CO3_CF + BOH4_CF + oh - hh; %ALK_CF
gammaFP = (0.5 - (ALK_CF - ALK_SW)./2./F_kzrho).^-1; % Gamma (= F/P)
P_kzrho = F_kzrho./gammaFP; % P/kzrho
D_kzrhoCO2 = D_kzrho.*(CO2_Cell - CO2_CF); %D*CO2/kzrho

%Calculation4. Add in carbon isotopes to this system in order to provide some
%other tracker of processes involved. 
R1312sw = 0.9985; % Seawater ratio of 13/12C
R1312cell = 0.985; % Cell ratio of 13/12C
    
R1312C_CF = (DIC_SW.*R1312sw + D_kzrho.*CO2_Cell.*R1312cell)./(DIC_CF + D_kzrho*CO2_CF + P_kzrho);
%DI12C_CF = -(R12C_TotCsw.*DIC_SW - F_kzrho.*((ALK_SW./2 - BOH4_CF./2 + hh./2 - oh./2)./F_kzrho + 1/2) + R12C_TotCcell*CO2_Cell.*D_kzrho)./(alpha_1./2 + alpha_2 - D_kzrho.*alpha_0 - 1); % calculating 12C-DIC abundance of the calcifying fluid at steady state
%R1312C_CF = DI13C_CF./DI12C_CF; % Calculating the ratio of 13C to 12C in the calcifying fluid at steady state
%ddic = (-1.5)*DIC_SW - DICCF - P/kzrho + D_kzrho(CO2_Cell- DIC_CF.*alpha_0)

%R1312C_CF = -((-1.5).*DIC_SW - F_kzrho.*((ALK_SW./2 - BOH4_CF./2 + hh./2 - oh./2)./F_kzrho + 1/2) + (-15)*CO2_Cell.*D_kzrho)./(alpha_1./2 + alpha_2 - D_kzrho.*alpha_0 - 1); % calculating 13C-DIC abundance of the calcifying fluid at steady state


% 'toc' to count time it takes for these calculations
toc %stop stopwatch timer
end


% Define function 'pHCFRule' to calculate internal calcifying fluid pH and H+. Uses observed
% relationship that delta H+ (sw-cf) is positively correlated with external
% seawater pCO2. 
function [pH_CF pCO2_SW] = pHCFRule(ALK_SW,DIC_SW,slope) 
    global S K_1 K_W  K_2 K_B I K_SO4 K_F B_T F_T SO4_T Ca_T K_sp
    
    % Take input ALK_SW and DIC_SW values and give them the appropriate IDs
    % to fit in with carbonate system calculations below
    TA = ALK_SW;
    DIC = DIC_SW;
    
    % solve for carbonate system of surrouding seawater
    % (probably good to re-check this, esp the boron part, it has been a while since I wrote it)
    fac = conv([1 0],conv([1 K_B],conv([1 K_1 (K_1*K_2)],[(1/K_2) (K_1/K_2) K_1])));  %factor multiplied to general TA equation
    temp1= conv([0 -TA],fac);  %TA =
    temp2= conv([0 0 (K_1*DIC) 0],deconv(fac, [1 K_1 (K_1*K_2)]));   %+HCO3
    temp6= conv([0 0 0 (2*K_1*DIC)],deconv(fac,[(1/K_2) (K_1/K_2) K_1]));  %+2CO3
    temp3= conv([0 0 K_W],deconv(fac,[1 0])); % +OH
    temp4= conv([-1 0],fac);  % -H
    temp5= conv([0 0 (K_B*B_T)],deconv(fac,[1 K_B]));  %+Anionic form of Buffer

    poly_to_solve = temp1+temp2+temp3+temp4+temp5+temp6;

    % only one positive root for carbonate system solution
    h = max(roots(poly_to_solve)); %this is [H+] in seawater
    alpha0 = 1./(1 + (K_1./h) + ((K_1.*K_2)./(h.*h))); %calculate alphas from H+
    alpha2 = K_1.*K_2 ./ (h.*h + h.*K_1 + K_1.*K_2);
    CO2_SW = alpha0.* DIC;    
    pCO2_SW = CO2_SW./S;
    
    %calculate pH_CF from Seawater Chemistry, using empirical slope rule
    %where deltaH+ is proportional to seawater pCO2. 
    delta_h = slope*pCO2_SW; % where delta_h is defined as h - h_cf (i.e. [H+]sw - [H+]cf)
    h_cf = -delta_h + h; % AMG changed the sign for delta_h to be consistent with the definition in comments above
    pH_CF = -log10(h_cf);
end

% Define function 'SolveSSCoral' to solve for DIC of the calcifying fluid
% as a function of ALK_SW, DIC_SW CO2_cell, F_kzrho, and D_kzrho
function [residual] = SolveSSCoral(DIC_CF, pH_CF, F_kzrho, D_kzrho, ALK_SW, DIC_SW, CO2_Cell)
    
    % declare global variables
    global S K_1 K_W  K_2 K_B I K_SO4 K_F B_T F_T SO4_T Ca_T K_sp
    
    %1. Use pH_CF and DIC_CF to calculate ALK_CF and CO2_CF
    hh = 10.^(-1.*pH_CF);
    oh = K_W./hh;
    alpha_0 = 1./(1 +K_1./hh +K_1.*K_2./hh./hh);
    alpha_1 = 1./(1 +hh./K_1 + K_2./hh);
    alpha_2 = 1./(1 +hh./K_2 + hh.*hh./K_1./K_2);
    CO2_CF = DIC_CF.*alpha_0;
    HCO3_CF= DIC_CF.*alpha_1;
    CO3_CF = DIC_CF.*alpha_2;
    BOH4_CF = B_T./(1+ hh./K_B);
    ALK_CF = HCO3_CF + 2.*CO3_CF + BOH4_CF + oh - hh;
    
    %2. Use ALK_CF and dALK/dt to get gamma
    % gamma is defined as F/P, the ratio of pumping to ppt
    gammaFP = (0.5 - (ALK_CF - ALK_SW)./2./F_kzrho).^-1;
    
    %3. Plug guess for DIC_CF into dDIC/dt equation together with calculated CO2_CF and gamma
    %this residual will be zero when the DIC_CF guess is correct. Residual
    %is a linear function of DIC_CF (see 'ResidualCheck.m' code to double
    %check this)
    residual = D_kzrho.*CO2_Cell - D_kzrho.*CO2_CF + DIC_SW - DIC_CF - F_kzrho./gammaFP;
    
end














