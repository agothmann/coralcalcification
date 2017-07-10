%SolveCFchem.m: Explore states of the coral calcifying fluid using
%the functions CoralCarbF.m from Alex which solves for the carbonate chemistry of the
%calcifying fluid and 'const_call_SWS_95_calcite_f.m' from Alex which kicks out
%carbonate system constants. 
%Anne Gothmann 23 June 2017 
%AMG 3 July 2017 edited

global S K_1 K_W K_2 K_B I K_sp K_SO4 K_F B_T F_T SO4_T Ca_T f_coeff_CO2_1atm;
global K_sp_calcite f_h;
global k_plus1 k_minus1 k_plus4 k_minus4 

clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%
% Setting Input Values %
%%%%%%%%%%%%%%%%%%%%%%%%
empirical_slope =  1.3e-5; %this function is defined at end of script and is set by user - AMG changed this value to be consisent with slopes from literature data
f = [1e-4:1e-3:10e-2]; %setting the range of values for F_kzrho - the ratio of alkalinity pumping over the seawater flux
d = [0.1:10:1000]; %setting the range of values for D_kzrho - the ratio of CO2 diffusion over the seawater flux
[F_kzrho,D_kzrho] = meshgrid(f,d); %create a meshgrid of F_kzrho and D_kzrho values for later contour plotting
dim = size(F_kzrho); % determine the dimensions of the meshgrid matrices

%setting DIC and ALK of seawater values...right now, ALK_SW and DIC_SW are
%chosen independently but in the future it might be useful to put these two
%parameters in a single matrix where pairs of DIC_SW and ALK_SW are set.
DIC_SW = [1800e-6]; 
ALK_SW = [2100e-6];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating conditions of the calcifying fluid %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%predefining variables that will be altered by the for loop below
loopcnt = 0;

%for loop where solutions of calcifying fluid parameters are generated for combinations of input ALK_SW, DIC_SW,
%F_kzrho, and D_kzrho values. Loop calls on 'CoralCarbF' function to
%perform calculations. The code within the jj= section of the for loop will
%only work for a single ALK_SW and DIC_SW input combination. Code below
%will have to be edited slightly if multiple DIC_SW and ALK_SW pairs are
%input. 
for AA = 1:length(ALK_SW);
    for DD = 1:length(DIC_SW);
        for ii = 1:dim(1);
            for jj = 1:dim(2);
                loopcnt = loopcnt + 1; % keep track of how many times this loop cycles
                [diccf pHcf alkcf co3cf co2cf pkzrho dkzrhoCO2 pco2cell pco2sw r1312cf] = CoralCarbF(ALK_SW(AA), DIC_SW(AA), F_kzrho(ii,jj), D_kzrho(ii,jj), empirical_slope); % calculate calcifying fluid parameters for a single combination of input variables
                DIC_CF(ii,jj) = diccf; % (this line and below) populate calcifying fluid parameters matrices with values as a function of D_kzrho and F_kzrho.
                pH_CF(ii,jj) = pHcf;
                ALK_CF(ii,jj) = alkcf;
                CO3_CF(ii,jj) = co3cf;
                CO2_CF(ii,jj) = co2cf;
                P_kzrho(ii,jj) = pkzrho;
                D_kzrhoCO2(ii,jj) = dkzrhoCO2;
                pCO2_Cell(ii,jj) = pco2cell;
                pCO2_SW(ii,jj) = pco2sw;
                F_kzrhob(ii,jj) = F_kzrho(ii,jj);
                D_kzrhob(ii,jj) = D_kzrho(ii,jj);
                R1312C_CF(ii,jj) = r1312cf;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove unrealistic cell values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Not all of the values calculated above are reasonable. For example,
% DIC_CF sometimes is negative, F/P is sometimes negative suggesting
% dissolution, etc. The below code assigns a value of '0' to any cells
% where either P_kzrho or D_CF are not greater than zero. 

for l = 1:dim(1)
    for ll = 1:dim(2)
        if P_kzrho(l,ll)>0 && DIC_CF(l,ll)>0
            DIC_CF(l,ll) = DIC_CF(l,ll);
            pH_CF(l,ll) = pH_CF(l,ll);
            ALK_CF(l,ll) = ALK_CF(l,ll);
            CO3_CF(l,ll) = CO3_CF(l,ll);
            CO2_CF(l,ll) =  CO2_CF(l,ll);
            P_kzrho(l,ll) = P_kzrho(l,ll);
            D_kzrhoCO2(l,ll) = D_kzrhoCO2(l,ll);
            pCO2_Cell(l,ll) = pCO2_Cell(l,ll);
            pCO2_SW(l,ll) = pCO2_SW(l,ll);
            F_kzrhob(l,ll) = F_kzrhob(l,ll);
            D_kzrhob(l,ll) = D_kzrhob(l,ll);
            gamma(l,ll) = F_kzrhob(l,ll)./P_kzrho(l,ll);
            psi(l,ll) = D_kzrhoCO2(l,ll)./P_kzrho(l,ll);
            F_SW(l,ll) = F_kzrhob(l,ll)./ALK_SW;
            D_SW(l,ll) = D_kzrhoCO2(l,ll)./DIC_SW;
            Ffactor(l,ll) = P_kzrho(l,ll)./ALK_SW;
            R1312C_CF(l,ll) = R1312C_CF(l,ll);
        else
            DIC_CF(l,ll) = 0;
            pH_CF(l,ll) = 0;
            ALK_CF(l,ll) = 0;
            CO3_CF(l,ll) = 0;
            CO2_CF(l,ll) = 0;
            P_kzrho(l,ll) = 0;
            D_kzrhoCO2(l,ll) = 0;
            pCO2_Cell(l,ll) = 0;
            pCO2_SW(l,ll) = 0;
            F_kzrhob(l,ll) = 0;
            D_kzrhob(l,ll) = 0;
            gamma(l,ll) = 0;
            psi(l,ll) = 0;
            F_SW(l,ll) = 0;
            D_SW(l,ll) = 0;
            Ffactor(l,ll) = 0;
            R1312C_CF(l,ll) = 0;
        end
    end
end

%%

%%%%%%%%%%%%
% Plotting %
%%%%%%%%%%%%
subplot(2,3,1)
contour(F_SW, D_SW, pH_CF,10,'Fill','on')
xlabel('F/SW')
ylabel('D/SW')
ce = colorbar;
ce.Label.String = 'pHcf';
hold on

%subplot(2,3,2)
%contour(D_SW./F_SW, Ffactor,DIC_CF,'Fill','on')
%xlabel('D/F')
%ylabel('P/kzrho[ALK]')
%c = colorbar;
%c.Label.String = 'DIC_CF';
%hold on

subplot(2,3,2)
contour(F_SW,D_SW,R1312C_CF,1000,'Fill','on')
xlabel('F/SW')
ylabel('D/SW')
c = colorbar;
c.Label.String = 'R1312';
hold on


subplot(2,3,3)
plot(DIC_CF, ALK_CF,'ok')
xlabel('DIC_CF')
ylabel('ALK_CF')
hold on

subplot(2,3,4)
contour(F_SW, D_SW, gamma, 1000,'Fill','on')
xlabel('F/SW')
ylabel('D/SW')
c = colorbar;
c.Label.String = 'F/P';
hold on

subplot(2,3,5)
contour(F_kzrhob, D_kzrhob, DIC_CF,10,'Fill','on')
xlabel('F/kzrho')
ylabel('D/kzrho')
c2 = colorbar;
c2.Label.String = 'DIC_CF (umol/kg)';
hold on

subplot(2,3,6)
contour(log10(F_SW), log10(D_SW), Ffactor,10,'Fill','on')
xlabel('F/SW')
ylabel('D/SW')
c3 = colorbar;
c3.Label.String = 'P/kzrho[ALK]';
hold on


%%%%%%%%%%%%%%%%%%%%%%%%%
% Extra code not in use %
%%%%%%%%%%%%%%%%%%%%%%%%%

%figure
%plot(DIC_CF,'o')
%subplot(1,3,2)
%scatter(F_kzrhob, P_kzrho,100, 'o')
%xlabel('F_kzrho')
%ylabel('P_kzrho')
%subplot(1,3,3)
%scatter(MatrixVals(:,3), MatrixVals(:,3)./MatrixVals(:,10),100,MatrixVals(:,8), 'o')
%xlabel('F_kzrho')
%ylabel('gammaFP')

% List of variables that correspond to MatrixVals columns:
% --------------------------------------------------------
% MatrixVals(:,1) = ALK_SW
% MatrixVals(:,2) = DIC_SW
% MatrixVals(:,3) = F_kzrho
% MatrixVals(:,4) = D_kzrho*(CO2_cell - CO2_cf)
% MatrixVals(:,5) = DIC_cf
% MatrixVals(:,6) = pH_cf
% MatrixVals(:,7) = ALK_cf
% MatrixVals(:,8) = CO3_cf
% MatrixVals(:,9) = CO2_cf
% MatrixVals(:,10) = P_kzrho
% MatrixVals(:,11) = D_kzrho
% MatrixVals(:,12) = pCO2_SW

% The ratio of pumping to precipitation, gamma, can be calculated either
% from the Alkalinity steady state equation or by simply dividing F_kzrho
% by P_kzrho. Choose one below: 
%gammaFP = 2./((MatrixVals(:,1)-MatrixVals(:,7))./MatrixVals(:,3) + 1);




%%%%%%%%%%
%Plotting%
%%%%%%%%%%

%pCO2_Cell = 1e6*pCO2_Cell;
%pCO2_SW = 1e6*MatrixVals(:,12);

%subplot(2,3,1)
%plot(MatrixVals(:,4), gammaFP, 'ok')
%xlabel('D*dCO2/kzrho')
%ylabel('gammaFP')
%hold on
%dim = [.2 .5 .3 .3];
%str = str(CO2_cell);
%s = int2str(pCO2_Cell);
%str = strcat('pCO2-cell = ',s,' ppm')
%a = annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none','Fontname','Arial','Fontsize',12,'FontWeight','Bold');



%subplot(2,3,2)
%scatter(MatrixVals(:,2).*1e6, MatrixVals(:,1).*1e6, 100, pCO2_SW, 'filled')
%xlabel('DIC_SW')
%ylabel('ALK_SW')
%hold on

%subplot(2,3,3)
%scatter(MatrixVals(:,4),MatrixVals(:,10), 100, MatrixVals(:,8),'filled')
%xlabel('D*dCO2/kzrho')
%ylabel('P_kzrho')
%hold on

%subplot(2,3,4)
%scatter(MatrixVals(:,5).*1e6, MatrixVals(:,7).*1e6,100, MatrixVals(:,8),'filled')
%xlabel('DIC_cf')
%ylabel('ALK_cf')
%hold on

%subplot(2,3,5)
%scatter(MatrixVals(:,4)./MatrixVals(:,10), MatrixVals(:,3)./MatrixVals(:,10),100, MatrixVals(:,12),'o','filled')
%xlabel('D*dCO2/P')
%ylabel('F/P')
%hold on

%subplot(2,3,6)
%scatter(D_kzrhoCO2./MatrixVals(:,2), MatrixVals(:,3)./MatrixVals(:,1),100, MatrixVals(:,5), 'o','filled')
%xlabel('D*dCO2/kzrho[DIC]sw')
%ylabel('F/kzrho[ALK]sw')
%hold on
