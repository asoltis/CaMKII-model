function ydot = soltis_biophysJ2010_BARsignalling_odefile(t,y,pin)
% Re-implemented by Anthony Soltis <ars7h@virginia.edu> for communication
% with Shannon (2004) EC coupling model. Final version 07/21/10

% saucerman_circres2004.m
% coupled beta adrenergic signaling/EC for adult rabbit ventricular
% myocytes, with extensions for yotiao interactions, KCNQ1-G589D mutation
%
% Copyright (2004) The Regents of the University of California
% All Rights Reserved
% Permission to use, copy, and modify, any part of this software for
% academic purposes only, including non-profit  education and research,
% without fee, and without a written agreement is hereby granted, provided
% that the above copyright notice, this paragraph and the following three
% paragraphs appear in all copies.
% The receiving party agrees not to further distribute nor disclose the
% source code to third parties without written permission and not to use
% the software for research under commercial sponsorship or for any other
% any commercial undertaking.  Those desiring to incorporate this software
% into commercial products of to use it for commercial purposes should
% contact the Technology Transfer Office, University of California, San
% Diego, 9500 Gilman Drive, La Jolla, California, 92093-0910, Ph: (619)
% 534 5815, referring to Case SDC98008.
% IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
% FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
% INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF
% THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
% DAMAGE.
% THE SOFTWARE PROVIDED HEREUNDER IS ON AN AS IS BASIS, AND THE
% UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE,
% SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  THE UNIVERSITY OF
% CALIFORNIA MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY
% KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MECHANT ABILITY OR FITNESS FOR A PARTICULAR
% PURPOSE, OR THAT THE USE OF THE MATERIAL WILL NOT INFRINGE ANY PATENT,
% TRADEMARK OR OTHER RIGHTS.
%
% Those using this code should acknowledge Dr. Andrew D. McCulloch and the
% National Biomedical Computation Resource (NBCR), NIH grant #P41 RR08605.
%
% References:
% -----------
% Jeffrey J. Saucerman, Sarah N. Healy, Mary E. Belik, Jose L. Puglisi, and
% Andrew D. McCulloch.  "Proarrhythmic Consequences of a KCNQ1 AKAP-Binding 
% Domain Mutation: Computational Models of Whole Cells and Heterogeneous
% Tissue", Circ. Res., Vol 95: 1216-1224 (2004).
%
% Jeffrey J. Saucerman and Andrew D. McCulloch
% "Mechanistic systems models of cell signaling networks: a case study of
% myocyte adrenergic regulation", Prog. Biophys. Mol. Bio., Vol 84: 261-278 (2004).
%
% Jeffrey J. Saucerman, Laurence L. Brunton, Anushka P. Michailova, and Andrew D. McCulloch
% "Modeling beta-adrenergic control of cardiac myocyte contractility in
% silico", J. Biol. Chem., Vol 278: 47977-48003 (2003).
%
% Last modified: 12/20/2004
% Implemented by: Jeffrey Saucerman <jsaucer@ucsd.edu>
%
% Notes:
% - This code was used for the single cell simulations in the Circ Res ms.
% - In Matlab 7, you can enable cell mode to jump between modules.
% - Please email me for any questions or comments.

%% Assign passed in params
Ltot = pin(1);
LCCtot = pin(2);
RyRtot = pin(3);
PLBtot = pin(4);
TnItot = pin(5);
IKstot = pin(6);
ICFTRtot = pin(7);
PP1_PLBtot = pin(8);

%% Parameters
%% ----- Signaling model parameters -------
% b-AR/Gs module
p(1) = Ltot;    % Ltotmax   [uM] ** apply agonist concentration here **
p(2) = 0.028;   % sumb1AR   [uM]
p(3) = 3.83;    % Gstot     [uM]
p(4) = 0.285;   % Kl        [uM]
p(5) = 0.062;   % Kr        [uM]
p(6) = 33.0;    % Kc        [uM]
p(7) = 1.1e-3;  % k_barkp   [1/sec]
p(8) = 2.2e-3;  % k_barkm   [1/sec]
p(9) = 3.6e-3;  % k_pkap    [1/sec/uM]
p(10) = 2.2e-3; % k_pkam    [1/sec]
p(11) = 16.0;   % k_gact    [1/sec]
p(12) = 0.8;    % k_hyd     [1/sec]
p(13) = 1.21e3; % k_reassoc [1/sec/uM]
% cAMP module
p(14) = 0.047;  % AC_tot    [uM]
p(15) = 5.0e3;  % ATP       [uM]
p(16) = 0.036;  % PDE3tot   [uM]   % Changed from .06 to .036
p(17) = 0.036;  % PDE4tot   [uM]
p(18) = 0.0;    % IBMXtot   [uM]
p(19) = 0.0;    % Fsktot    [uM] (10 uM when used)
p(20) = 0.2;    % k_ac_basal[1/sec]
p(21) = 8.5;    % k_ac_gsa  [1/sec]
p(22) = 7.3;    % k_ac_fsk  [1/sec]
p(23) = 1.03e3; % Km_basal  [uM]
p(24) = 315.0;  % Km_gsa    [uM]
p(25) = 860.0;  % Km_fsk    [uM]
p(26) = 0.4;    % Kgsa      [uM]
p(27) = 44.0;   % Kfsk      [uM]
p(28) = 3.5;    % k_pde3    [1/sec]
p(29) = 0.15;   % Km_pde3   [uM]
p(30) = 5.0;    % k_pde4    [1/sec]
p(31) = 1.3;    % Km_pde4   [uM]
p(32) = 30.0;   % Ki_ibmx   [uM]
% PKA module
p(33) = 0.46;   % PKAItot   [uM]
p(34) = 0.084;  % PKAIItot  [uM]
p(35) = 0.18;   % PKItot    [uM]
p(36) = 9.14;   % Ka        [uM]
p(37) = 1.64;   % Kb        [uM]
p(38) = 4.375;  % Kd        [uM]
p(39) = 0.2e-3; % Ki_pki    [uM]
% PLB module
p(40) = 10;     % epsilon   [none]
p(41) = PLBtot; % PLBtot    [uM]
p(42) = PP1_PLBtot;   % PP1tot    [uM]
p(43) = 0.3;    % Inhib1tot [uM]
p(44) = 54;     % k_pka_plb     [1/sec]
p(45) = 21;     % Km_pka_plb    [uM]
p(46) = 8.5;    % k_pp1_plb     [1/sec]
p(47) = 7.0;    % Km_pp1_plb    [uM]
p(48) = 60;     % k_pka_i1      [1/sec]
p(49) = 1.0;    % Km_pka_i1     [uM]
p(50) = 14.0;   % Vmax_pp2a_i1  [uM/sec]
p(51) = 1.0;    % Km_pp2a_i1    [uM]
p(52) = 1.0e-3; % Ki_inhib1     [uM]
% LCC module
p(53) = LCCtot; % LCCtot        [uM]
p(54) = 0.025;  % PKAIIlcctot   [uM]
p(55) = 0.025;  % PP1lcctot     [uM]
p(56) = 0.025;  % PP2Alcctot    [uM]
p(57) = 54;     % k_pka_lcc     [1/sec]
p(58) = 21;     % Km_pka_lcc    [uM]
p(59) = 8.52;   % k_pp1_lcc     [1/sec]
p(60) = 3;      % Km_pp1_lcc    [uM]
p(61) = 10.1;   % k_pp2a_lcc    [1/sec]
p(62) = 3;      % Km_pp2a_lcc   [uM]
% RyR module
p(63) = RyRtot; % RyRtot        [uM]
p(64) = 0.034;  % PKAIIryrtot   [uM]
p(65) = 0.034;  % PP1ryr        [uM]
p(66) = 0.034;  % PP2Aryr       [uM]
p(67) = 54;     % kcat_pka_ryr  [1/sec]
p(68) = 21;     % Km_pka_ryr    [uM]
p(69) = 8.52;   % kcat_pp1_ryr  [1/sec]
p(70) = 7;      % Km_pp1_ryr    [uM]
p(71) = 10.1;   % kcat_pp2a_ryr [1/sec]
p(72) = 4.1;    % Km_pp2a_ryr   [uM]
% TnI module
p(73) = TnItot; % TnItot        [uM]
p(74) = 0.67;   % PP2Atni       [uM]
p(75) = 54;     % kcat_pka_tni  [1/sec]
p(76) = 21;     % Km_pka_tni    [uM]
p(77) = 10.1;   % kcat_pp2a_tni [1/sec]
p(78) = 4.1;    % Km_pp2a_tni   [uM]
% Iks module
p(79) = IKstot; % Iks_tot       [uM]
p(80) = 0.025;  % Yotiao_tot    [uM]
p(81) = 0.1e-3; % K_yotiao      [uM] ** apply G589D mutation here **
p(82) = 0.025;  % PKAII_ikstot  [uM]
p(83) = 0.025;  % PP1_ikstot    [uM]
p(84) = 54;     % k_pka_iks     [1/sec]
p(85) = 21;     % Km_pka_iks    [uM]
p(86) = 8.52;   % k_pp1_iks     [1/sec]
p(87) = 7;      % Km_pp1_iks    [uM]
% Icftr Module - Added 04/30/10 by Anthony Soltis
p(88) = ICFTRtot;  % CFTR_tot      [uM]
p(89) = 0.025;  % PKAII_CFTRtot [uM]
p(90) = 0.025;  % PP1_CFTRtot   [uM]
p(91) = 54;     % k_pka_CFTR    [1/sec]
p(92) = 8.5;    % Km_pka_CFTR   [uM]
p(93) = 8.52;   % k_pp1_CFTR    [1/sec]
p(94) = 7;      % Km_pp1_CFTR   [uM]
%% -------- SIGNALING MODEL -----------
ydot = zeros(size(y));
%% b-AR module
LR = y(1)*y(2)/p(4);
LRG = LR*y(3)/p(5);
RG = y(2)*y(3)/p(6);
BARKDESENS = p(7)*(LR+LRG);
BARKRESENS = p(8)*y(5);
PKADESENS = p(9)*y(17)*y(4);  
PKARESENS = p(10)*y(6);
GACT = p(11)*(RG+LRG);
HYD = p(12)*y(7);
REASSOC = p(13)*y(8)*y(9);
ydot(1) = p(1)-LR-LRG-y(1);
ydot(2) = y(4)-LR-LRG-RG-y(2);
ydot(3) = p(3)-LRG-RG-y(3);
ydot(4) = (BARKRESENS-BARKDESENS)+(PKARESENS-PKADESENS);
ydot(5) = BARKDESENS-BARKRESENS;
ydot(6) = PKADESENS-PKARESENS;
ydot(7) = GACT-HYD;
ydot(8) = HYD-REASSOC;
ydot(9) = GACT-REASSOC;
% end b-AR module

%% cAMP module
Gsa_gtp_AC = y(10)*y(12)/p(26);
Fsk_AC = y(11)*y(12)/p(27);
AC_ACT_BASAL = p(20)*y(12)*p(15)/(p(23)+p(15));	    
AC_ACT_GSA = p(21)*Gsa_gtp_AC*p(15)/(p(24)+p(15)); 
AC_ACT_FSK = p(22)*Fsk_AC*p(15)/(p(25)+p(15));	   
% PDE3_ACT = p(28)*y(13)*y(16)/(p(29)+y(16));	
% PDE4_ACT = p(30)*y(13)*y(16)/(p(31)+y(16));	
PDE3_ACT = p(28)*p(16)*y(16)/(p(29)*(1+p(18)/p(32))+y(16));	% new PDE3 term w IBMX
PDE4_ACT = p(30)*p(17)*y(16)/(p(31)*(1+p(18)/p(32))+y(16));	% new PDE4 term w IBMX
PDE_IBMX = y(13)*y(14)/p(32);
ydot(10) = y(7)-Gsa_gtp_AC-y(10);
ydot(11) = p(19)-Fsk_AC-y(11);
ydot(12) = p(14)-Gsa_gtp_AC-y(12);  % note: assumes Fsk = 0.  Change Gsa_gtp_AC to Fsk_AC for Forskolin.
% ydot(13) = p(17)-PDE_IBMX-y(13);
% ydot(14) = p(18)-PDE_IBMX-y(14);
ydot(13) = 0;
ydot(14) = 0;
ydot(15) = AC_ACT_BASAL+AC_ACT_GSA+AC_ACT_FSK-PDE3_ACT-PDE4_ACT;
% end cAMP module

%% PKA module
PKI = p(35)*p(39)/(p(39)+y(17)+y(18));
A2RC_I = (y(17)/p(38))*y(17)*(1+PKI/p(39));
A2R_I = y(17)*(1+PKI/p(39));
A2RC_II = (y(18)/p(38))*y(18)*(1+PKI/p(39));
A2R_II = y(18)*(1+PKI/p(39));
ARC_I = (p(36)/y(16))*A2RC_I;
ARC_II = (p(36)/y(16))*A2RC_II;
ydot(16) = y(15)-(ARC_I+2*A2RC_I+2*A2R_I)-(ARC_II+2*A2RC_II+2*A2R_II)-y(16);
PKAtemp = p(36)*p(37)/p(38)+p(36)*y(16)/p(38)+y(16)^2/p(38);
ydot(17) = 2*p(33)*y(16)^2-y(17)*(1+PKI/p(39))*(PKAtemp*y(17)+y(16)^2);
ydot(18) = 2*p(34)*y(16)^2-y(18)*(1+PKI/p(39))*(PKAtemp*y(18)+y(16)^2);
% end PKA module

%% PLB module
PLB = p(41)-y(19);
PLB_PHOSPH = p(44)*y(17)*PLB/(p(45)+PLB);
PLB_DEPHOSPH = p(46)*y(22)*y(19)/(p(47)+y(19));
ydot(19) = PLB_PHOSPH-PLB_DEPHOSPH;
 
Inhib1 = p(43)-y(20);
Inhib1p_PP1 = y(21)*y(22)/p(52);
Inhib1_PHOSPH = p(48)*y(17)*Inhib1/(p(49)+Inhib1); 
Inhib1_DEPHOSPH = p(50)*y(20)/(p(51)+y(20));
ydot(20) = Inhib1_PHOSPH-Inhib1_DEPHOSPH;
ydot(21) = y(20)-Inhib1p_PP1-y(21);
ydot(22) = p(42)-Inhib1p_PP1-y(22);
% end PLB module

%% LCC module
PKAClcc = (p(54)/p(34))*y(18);
LCCa = p(53)-y(23);
LCCa_PHOSPH = p(40)*p(57)*PKAClcc*LCCa/(p(58) + p(40)*LCCa);
LCCa_DEPHOSPH = p(40)*p(61)*p(56)*y(23)/(p(62)+p(40)*y(23));
ydot(23) = LCCa_PHOSPH - LCCa_DEPHOSPH;
 
LCCb = p(53)-y(24);
LCCb_PHOSPH = p(40)*p(57)*PKAClcc*LCCb/(p(58)+p(40)*LCCb);   
LCCb_DEPHOSPH = p(40)*p(59)*p(55)*y(24)/(p(60)+p(40)*y(24));
ydot(24) = LCCb_PHOSPH-LCCb_DEPHOSPH;
% end LCC module

%% RyR module
PKACryr = (p(64)/p(34))*y(18);
RyR = p(63)-y(25);
RyRPHOSPH = p(40)*p(67)*PKACryr*RyR/(p(68)+p(40)*RyR);
RyRDEPHOSPH1 = p(40)*p(69)*p(65)*y(25)/(p(70)+p(40)*y(25));
RyRDEPHOSPH2A = p(40)*p(71)*p(66)*y(25)/(p(72)+p(40)*y(25));
ydot(25) = RyRPHOSPH-RyRDEPHOSPH1-RyRDEPHOSPH2A;
% end RyR module

%% TnI module
TnI = p(73)-y(26);
TnIPHOSPH = p(75)*y(17)*TnI/(p(76)+TnI);
TnIDEPHOSPH = p(77)*p(74)*y(26)/(p(78)+y(26));
ydot(26) = TnIPHOSPH-TnIDEPHOSPH;
% end TnI module

%% Iks module
IksYot = y(27)*y(28)/p(81);           % [uM]
ydot(27) = p(79) - IksYot - y(27);    % [uM]
ydot(28) = p(80) - IksYot - y(28);    % [uM]
PKACiks = (IksYot/p(79))*(p(82)/p(34))*y(18);
PP1iks = (IksYot/p(79))*p(83);
Iks = p(79)-y(29);
IKS_PHOSPH = p(40)*p(84)*PKACiks*Iks/(p(85)+p(40)*Iks);
IKS_DEPHOSPH = p(40)*p(86)*PP1iks*y(29)/(p(87)+p(40)*y(29));
ydot(29) = IKS_PHOSPH-IKS_DEPHOSPH;

%% CFTR module (included 04/30/10)
CFTRn = p(88) - y(30);  % Non-phos = tot - phos
PKAC_CFTR = (p(89)/p(34))*y(18);    % (PKACFTRtot/PKAIItot)*PKAIIact
CFTRphos = p(40)*CFTRn*PKAC_CFTR*p(91)/(p(92)+p(40)*CFTRn);
CFTRdephos = p(90)*p(93)*p(40)*y(30)/(p(94) + p(40)*y(30));
ydot(30) = CFTRphos - CFTRdephos;
%% Gather odes
% Need to convert all ydot terms that are ODEs (not DAEs) to miliseconds
odes = [4,5,6,7,8,9,13,14,15,19,20,23,24,25,26,29,30];
ydot(odes) = ydot(odes).*1e-3;