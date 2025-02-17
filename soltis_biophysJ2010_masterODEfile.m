function dydt = soltis_biophysJ2010_masterODEfile(t,y,p)
% This function calls the ode files for EC coupling, CaM reactions, CaMKII
% phosphorylation module, and PKA phosphorylation module

% Dependencies: soltis_biophysJ2010_camODEfile.m,
%               soltis_biophysJ2010_eccODEfile.m
%
% Author: Jeff Saucerman <jsaucerman@virginia.edu>
% Copyright 2008, University of Virginia, All Rights Reserved
%
% References:
% -----------
% JJ Saucerman and BM Bers, Calmodulin mediates differential
% sensitivity of CaMKII and calcineurin to local Ca2+ in cardiac myocytes. 
% Biophys J. 2008 Aug 8. [Epub ahead of print] 
% Please cite the above paper when using this model.
%% Collect params and ICs for each module
% Allocate ICs for each moduel
% y(60) -> y(65) are state transitions for mode 1 junctional LCCs
% y(66) -> y(71) are state transitoins for mode 2 junctional LCCs
% y(72) -> y(77) are state transitions for mode 1 sarcolemmal LCCs
% y(78) -> y(83) are state transitions for mode 2 sarcolemmal LCCs

y_ecc = y(1:83);  % Ca_j is y(36), Ca_sl is y(37), Ca_cytosol is y(38)
y_camDyad = y(83+1:83+15);
y_camSL = y(83+15+1:83+15+15);
y_camCyt = y(83+15+15+1:83+15+15+15);
y_CaMKII = y(83+45+1:83+45+6);      % 6 state vars in CaMKII module
y_BAR = y(83+45+6+1:83+45+6+30);    % 30 state vars in BAR module

% break up parameters from master function into modules
paramsCell=mat2cell(p,ones(size(p,1),1),ones(size(p,2),1));
[cycleLength,recoveryTime,CaMtotDyad,BtotDyad,CaMKIItotDyad,CaNtotDyad,PP1totDyad,...
    CaMtotSL,BtotSL,CaMKIItotSL,CaNtotSL,PP1totSL,...
    CaMtotCyt,BtotCyt,CaMKIItotCyt,CaNtotCyt,PP1totCyt...
    LCCtotDyad,RyRtot,PP1_dyad,PP2A_dyad,OA,PLBtot,LCCtotSL,PP1_SL...
    Ligtot,LCCtotBA,RyRtotBA,PLBtotBA,TnItotBA,IKstotBA,ICFTRtotBA,PP1_PLBtot...
    CKIIOE]=paramsCell{:};

K = 135; % [mM]
Mg = 1;  % [mM]

%% Distribute parameters by module

% CaM module
CaDyad = y(36)*1e3; % from ECC model, *** Converting from [mM] to [uM] ***
compart_dyad = 2;
% ** NOTE: Btotdyad being sent to the dyad camODEfile is set to zero, but is used below for transfer between SL and dyad
pCaMDyad = [K, Mg, CaMtotDyad, 0, CaMKIItotDyad, CaNtotDyad, PP1totDyad, CaDyad, cycleLength, compart_dyad];
CaSL = y(37)*1e3; % from ECC model, *** Converting from [mM] to [uM] ***
compartSL = 1;
pCaMSL = [K, Mg, CaMtotSL, BtotSL, CaMKIItotSL, CaNtotSL, PP1totSL, CaSL, cycleLength,compartSL];
CaCyt = y(38)*1e3; % from ECC model, *** Converting from [mM] to [uM] ***
compartCyt = 0;
pCaMCyt = [K, Mg, CaMtotCyt, BtotCyt, CaMKIItotCyt, CaNtotCyt, PP1totCyt, CaCyt, cycleLength,compartCyt];

% CaMKII phosphorylation module 
CaMKIIact_Dyad = CaMKIItotDyad.*(y(83+8)+y(83+9)+y(83+10)+y(83+11)); % Multiply total by fraction
CaMKIIact_SL = CaMKIItotSL.*(y(83+15+8)+y(83+15+9)+y(83+15+10)+y(83+15+11));
PP1_PLB_avail = y(83+45+6+22)./PP1_PLBtot + .0091;  % Active PP1 near PLB / total PP1 conc + basal value
pCaMKII = [CaMKIIact_Dyad,LCCtotDyad,RyRtot,PP1_dyad,PP2A_dyad,OA,PLBtot,...
           CaMKIIact_SL,LCCtotSL,PP1_SL,...
           PP1_PLB_avail];

LCC_CKdyadp = y(83+45+2)./LCCtotDyad;   % fractional CaMKII-dependent LCC dyad phosphorylation
RyR_CKp = y(83+45+4)./RyRtot;           % fractional CaMKII-dependent RyR phosphorylation
PLB_CKp = y(83+45+5)./PLBtot;           % fractional CaMKII-dependent PLB phosphorylation

% BAR (PKA phosphorylation) module
pBAR = [Ligtot,LCCtotBA,RyRtotBA,PLBtotBA,TnItotBA,IKstotBA,ICFTRtotBA,PP1_PLBtot];
LCCa_PKAp = y(83+45+6+23)./LCCtotBA;
LCCb_PKAp = y(83+45+6+24)./LCCtotBA;
PLB_PKAn = (PLBtotBA - y(83+45+6+19))./PLBtotBA;
RyR_PKAp = y(83+45+6+25)./RyRtotBA;
TnI_PKAp = y(83+45+6+26)./TnItotBA;
IKs_PKAp = y(83+45+6+29)./IKstotBA;
ICFTR_PKAp = y(83+45+6+30)./ICFTRtotBA;

% KO phospho-regulation
% LCC_CKp = 0;
% RyR_CKp = 0;
% PLB_CKp = 0;

% Temporary KO (resting values)
% RyR_CKp = 57.3715/RyRtot;
% PLB_CKp = 2.0751e-322/PLBtot;
% LCC_CKp = 5.2773e-60/LCCtot;

% Temporary Isolation OE (1 Hz) - fixed at normal phos levels - updated 05/10/10
% RyR_CKp = 74.7073/RyRtot;
% LCC_CKp = 15.4709/LCCtot;
% PLB_CKp = .4568/PLBtot;

% Temporary Isolation OE (2 Hz) - fixed at normal phos levels
% USED FOR ANALYSIS OF CAMKII FEEDBACK
% LCC_CKp = 28.0510/LCCtot;
% RyR_CKp = 97.0553/RyRtot;
% PLB_CKp = 1.1833/PLBtot;

% ISO Test (2 Hz) - Fix phos levels at ss normal ISO values
% USED FOR ARRHYTHMIA ANALYSIS
% RyR_CKp = 101.9432/RyRtot;
% LCC_CKp = 28.3810/LCCtot;
% PLB_CKp = 1.3982/PLBtot;

% ISO Test (1 Hz) - using KO value
% RyR_CKp = .15;
% LCC_CKp = 0;

% ECC module
pECC = [cycleLength,LCC_CKdyadp,RyR_CKp,PLB_CKp,...
        LCCa_PKAp,LCCb_PKAp,PLB_PKAn,RyR_PKAp,TnI_PKAp,IKs_PKAp,ICFTR_PKAp,...
        CKIIOE,recoveryTime];

%% Solve dydt in each module
dydt_ecc = soltis_biophysJ2010_eccODEfile(t,y_ecc,pECC);
dydt_camDyad = soltis_biophysJ2010_camODEfile(t,y_camDyad,pCaMDyad);
dydt_camSL = soltis_biophysJ2010_camODEfile(t,y_camSL,pCaMSL);
dydt_camCyt = soltis_biophysJ2010_camODEfile(t,y_camCyt,pCaMCyt);
dydt_CaMKIIDyad = soltis_biophysJ2010_CaMKII_signaling_ODEfile(t,y_CaMKII,pCaMKII);
dydt_BAR = soltis_biophysJ2010_BARsignalling_odefile(t,y_BAR,pBAR);

% incorporate Ca buffering from CaM, convert JCaCyt from uM/msec to mM/msec
global JCaCyt JCaSL JCaDyad;
dydt_ecc(36) = dydt_ecc(36) + 1e-3*JCaDyad;
dydt_ecc(37) = dydt_ecc(37) + 1e-3*JCaSL;
dydt_ecc(38) = dydt_ecc(38) + 1e-3*JCaCyt; 

% incorporate CaM diffusion between compartments
Vmyo = 2.1454e-11;      % [L]
Vdyad = 1.7790e-014;    % [L]
VSL = 6.6013e-013;      % [L]
kDyadSL = 3.6363e-16;	% [L/msec]
kSLmyo = 8.587e-15;     % [L/msec]
k0Boff = 0.0014;        % [s^-1] 
k0Bon = k0Boff/0.2;     % [uM^-1 s^-1] kon = koff/Kd
k2Boff = k0Boff/100;    % [s^-1] 
k2Bon = k0Bon;          % [uM^-1 s^-1]
k4Boff = k2Boff;        % [s^-1]
k4Bon = k0Bon;          % [uM^-1 s^-1]
CaMtotDyad = sum(y_camDyad(1:6))+CaMKIItotDyad*sum(y_camDyad(7:10))+sum(y_camDyad(13:15));
Bdyad = BtotDyad - CaMtotDyad; % [uM dyad]
J_cam_dyadSL = 1e-3*(k0Boff*y_camDyad(1) - k0Bon*Bdyad*y_camSL(1)); % [uM/msec dyad]
J_ca2cam_dyadSL = 1e-3*(k2Boff*y_camDyad(2) - k2Bon*Bdyad*y_camSL(2)); % [uM/msec dyad]
J_ca4cam_dyadSL = 1e-3*(k2Boff*y_camDyad(3) - k4Bon*Bdyad*y_camSL(3)); % [uM/msec dyad]
J_cam_SLmyo = kSLmyo*(y_camSL(1)-y_camCyt(1)); % [umol/msec]
J_ca2cam_SLmyo = kSLmyo*(y_camSL(2)-y_camCyt(2)); % [umol/msec]
J_ca4cam_SLmyo = kSLmyo*(y_camSL(3)-y_camCyt(3)); % [umol/msec]
dydt_camDyad(1) = dydt_camDyad(1) - J_cam_dyadSL;
dydt_camDyad(2) = dydt_camDyad(2) - J_ca2cam_dyadSL;
dydt_camDyad(3) = dydt_camDyad(3) - J_ca4cam_dyadSL;
dydt_camSL(1) = dydt_camSL(1) + J_cam_dyadSL*Vdyad/VSL - J_cam_SLmyo/VSL;
dydt_camSL(2) = dydt_camSL(2) + J_ca2cam_dyadSL*Vdyad/VSL - J_ca2cam_SLmyo/VSL;
dydt_camSL(3) = dydt_camSL(3) + J_ca4cam_dyadSL*Vdyad/VSL - J_ca4cam_SLmyo/VSL;
dydt_camCyt(1) = dydt_camCyt(1) + J_cam_SLmyo/Vmyo;
dydt_camCyt(2) = dydt_camCyt(2) + J_ca2cam_SLmyo/Vmyo;
dydt_camCyt(3) = dydt_camCyt(3) + J_ca4cam_SLmyo/Vmyo;

% Collect all dydt terms
dydt = [dydt_ecc; dydt_camDyad; dydt_camSL; dydt_camCyt; dydt_CaMKIIDyad; dydt_BAR];