function yfinal = soltis_biophysJ2010_masterCompute
% This function calls the ode solver and plots results.
% Dependency: soltis_biophysJ2010_masterODEfile.m
%
% Re-implemented by Anthony Soltis <ars7h@virginia.edu> for CaMKII/PKA
% regulation of ECC model
%
% Author: Jeff Saucerman <jsaucerman@virginia.edu>
% Copyright 2008, University of Virginia, All Rights Reserved
%
% Reference: JJ Saucerman and DM Bers, Calmodulin mediates differential
% sensitivity of CaMKII and calcineurin to local Ca2+ in cardiac myocytes. 
% Biophys J. 2008 Aug 8. [Epub ahead of print] 
% Please cite the above paper when using this model.

clear all; 
clc;
%% Parameters for external modules
% ECC and CaM modules
freq = 1.0;                 % [Hz] CHANGE DEPENDING ON FREQUENCY
cycleLength = 1e3/freq;     % [ms]
CaMtotDyad = 418;           % [uM]
BtotDyad = 1.54/8.293e-4;   % [uM]
CaMKIItotDyad = 120;        % [uM] 
CaNtotDyad = 3e-3/8.293e-4; % [uM] 
PP1totDyad = 96.5;          % [uM]
CaMtotSL = 5.65;            % [uM]
BtotSL = 24.2;              % [uM]
CaMKIItotSL = 120*8.293e-4; % [uM]
CaNtotSL = 3e-3;            % [uM]
PP1totSL = 0.57;            % [uM]
CaMtotCyt = 5.65;           % [uM]
BtotCyt = 24.2;             % [uM]
CaMKIItotCyt = 120*8.293e-4;% [uM]
CaNtotCyt = 3e-3;           % [uM] 
PP1totCyt = 0.57;           % [uM]

% ADJUST CAMKII ACTIVITY LEVELS (expression = 'WT', 'OE', or 'KO')
expression = 'WT';
CKIIOE = 0; % Should be zero during 'WT' and 'KO' runs

if strcmp(expression,'OE')
    CKIIOE = 1; % Flag for CKII OE in ECC file (0=WT, 1=OE) - for Ito and INa
    CaMKIItotDyad = 120*6;        % [uM] 
    CaMKIItotSL = 120*8.293e-4*6; % [uM]
    CaMKIItotCyt = 120*8.293e-4*6;% [uM]
elseif strcmp(expression,'KO')
    CaMKIItotDyad = 0;          % [uM] 
    CaMKIItotSL = 0;            % [uM]
    CaMKIItotCyt = 0;           % [uM]
end

% For Recovery from inactivation of LCC
recoveryTime = 10;  % initialize to smallest value

% Parameters for CaMKII module
LCCtotDyad = 31.4*.9;       % [uM] - Total Dyadic [LCC] - (umol/l dyad)
LCCtotSL = 0.0846;          % [uM] - Total Subsarcolemmal [LCC] (umol/l sl)
RyRtot = 382.6;             % [uM] - Total RyR (in Dyad)
PP1_dyad = 95.7;            % [uM] - Total dyadic [PP1]
PP1_SL = 0.57;              % [uM] - Total Subsarcolemmal [PP1]
PP2A_dyad = 95.76;          % [uM] - Total dyadic PP2A
OA = 0;                     % [uM] - PP1/PP2A inhibitor Okadaic Acid
PLBtot = 38;                % [uM] - Total [PLB] in cytosolic units

% Parameters for BAR module
Ligtot = 0;                 % [uM] - SET LIGAND CONCENTRATION HERE
LCCtotBA = 0.025;           % [uM] - [umol/L cytosol]
RyRtotBA = 0.135;           % [uM] - [umol/L cytosol]
PLBtotBA = 38;              % [uM] - [umol/L cytosol]
TnItotBA = 70;              % [uM] - [umol/L cytosol]
IKstotBA = 0.025;           % [uM] - [umol/L cytosol]
ICFTRtotBA = 0.025;         % [uM] - [umol/L cytosol]
PP1_PLBtot = 0.89;          % [uM] - [umol/L cytosol]

%% Collect all parameters and define mass matrix for BAR module
p = [cycleLength,recoveryTime,CaMtotDyad,BtotDyad,CaMKIItotDyad,CaNtotDyad,PP1totDyad,...
    CaMtotSL,BtotSL,CaMKIItotSL,CaNtotSL,PP1totSL,...
    CaMtotCyt,BtotCyt,CaMKIItotCyt,CaNtotCyt,PP1totCyt...
    LCCtotDyad,RyRtot,PP1_dyad,PP2A_dyad,OA,PLBtot,LCCtotSL,PP1_SL...
    Ligtot,LCCtotBA,RyRtotBA,PLBtotBA,TnItotBA,IKstotBA,ICFTRtotBA,PP1_PLBtot,...
    CKIIOE];

% Need to define Mass matrix for BAR portion
m0 = ones(1,83+45+6);   % All state variables in other modules, not BAR
m1 =[0,     0,      0,      1,      1,      1,      1,          1,      1];
m2 =[0,     0,      0,      1,      1,      1];
m3 =[0,     0,      0];
m4 =[1,     1,          0,      0,      1,          1];
m5 =[1,     1,      0,      0,      1,      1];
M = diag([m0,m1,m2,m3,m4,m5]); 

%% Establish and define globals
global tStep tArray I_Ca_store I_to_store I_Na_store I_K1_store ibar_store gates 
global Jserca IKs_store Jleak ICFTR Incx
tStep = 1; tArray = zeros(1,1e6); I_Ca_store=zeros(1,1e6); I_to_store=zeros(3,1e6);
I_Na_store = zeros(1,1e6); I_K1_store = zeros(1,1e6); ibar_store=zeros(1,1e6);
gates = zeros(2,1e6); Jserca = zeros(1,1e6); IKs_store = zeros(1,1e6); Jleak = zeros(1e6,2); ICFTR = zeros(1,1e6);
Incx = zeros(1,1e6);
%% ASSIGN INITIAL CONDITIONS

% WT CaMKII runs, Ltot = 0
% y0n = load('yfinal_0p125Hz_051210.dat');
% y0n = load('yfinal_0p25Hz_051210.dat');
% y0n = load('yfinal_0p5Hz_051210.dat');
y0n = load('yfinal_1Hz_051210.dat');
% y0n = load('yfinal_2Hz_051210.dat');
% y0n = load('yfinal_3Hz_051210.dat');

% OE CaMKII runs, Ltot = 0
% y0n = load('yfinal_0p125Hz_OE_051210.dat');
% y0n = load('yfinal_0p25Hz_OE_051210.dat');
% y0n = load('yfinal_0p5Hz_OE_051210.dat');
% y0n = load('yfinal_1Hz_OE_051210.dat');
% y0n = load('yfinal_2Hz_OE_051210.dat');
% y0n = load('yfinal_3Hz_OE_051210.dat');

% KO CaMKII, Ltot = 0
% y0n=load('0p125Hz_KOCKII_052110.dat');
% y0n=load('0p25Hz_KOCKII_052110.dat');
% y0n=load('yfinal_0p5Hz_KO_CKII_051210.dat');
% y0n=load('yfinal_1Hz_KO_CKII_051210.dat');
% y0n=load('yfinal_2Hz_KO_CKII_051210.dat');
% y0n=load('yfinal_3Hz_KO_CKII_051210.dat');

% VOLTAGE CLAMP (295 K)
% y0n=load('yfinal_min90mVRest_051310.dat');
% y0n=load('yfinal_min140mV_295K_052410.dat');

% VOLTAGE CLAMP (310 K)
% y0n=load('yfinal_min80mV_310K_060110.dat');
% y0n=load('yfinal_min80mV_WTCKII_2Hz_vclamp_SS_310K_060110.dat');
% y0n=load('yfinal_min80mV_OECKII_2Hz_vclamp_SS_310K_060110.dat');
% y0n=load('yfinal_min80mV_KOCKII_2Hz_vclamp_SS_310K_060110.dat');

% LONGER ISO RUNS 
% WT CaMKII, Ltot = 1 uM
% y0n=load('yfinal_0p5Hz_WT_CKII_1uMISO_70s_051210.dat');
% y0n=load('yfinal_1Hz_WT_1uM_ISO_70sRun_051210.dat');
% y0n=load('yfinal_2Hz_WT_1uM_ISO_60sRun_051210.dat');
% OE CaMKII, Ltot = 1 uM
% y0n=load('yfinal_1Hz_OE_CKII_1uMISO_70s_051210.dat');
% y0n=load('yfinal_2Hz_OE_1uM_ISO_60sRun_051210.dat');
% y0n=load('yfinal_3Hz_OE_1uM_ISO_60sRun_051910.dat');
% KO CaMKII, Ltot = 1 uM
% y0n=load('yfinal_0p5Hz_KO_CKII_1uMISO_70s_051210.dat');
% y0n=load('yfinal_1Hz_KO_CKII_1uMISO_60sRun_051210.dat');
% y0n=load('yfinal_1Hz_KO_CKII_1uMISO_70s_051210.dat');
% y0n=load('yfinal_2Hz_KO_CKII_1uMISO_60s_051210.dat');


% INDIVIDUAL TARGET EXPERIMENTS (2 HZ PACING)
% y0n=load('yfinal_2Hz_OE_CKII_LCCEffects_only_051310.dat');
% y0n=load('yfinal_2Hz_OE_CKII_RyREffects_only_051310.dat');
% y0n=load('yfinal_2Hz_OE_CKII_PLBEffects_only_051310.dat');
% y0n=load('yfinal_2Hz_OE_CKII_INa_ItoEffects_only_051310.dat');

% Tests without INa or Ito Changes
% y0n=load('yfinal_1Hz_OE_noItoINaChanges.dat');

% CaMKII activity tests (05/17/10)
% y0n=load('yfinal_2Hz_OE_CKII_noTargetEffects.dat');

% INHIBITOR-1 TESTS
% y0n=load('2Hz_WTCaMKII_no_I1Reg_1uMISO_051810.dat');
% y0n=load('2Hz_OECaMKII_no_I1Reg_1uMISO_051810.dat');

% RYR EFFECT TESTS
% y0n=load('yfinal_1Hz_RyREffect_KO.dat');
% y0n=load('yfinal_2Hz_OECKII_RyREffect_KO.dat');
% y0n=load('yfinal_2Hz_OECKII_RyREffect_KO_1uMISO_60s.dat');
%% Run single simulation
color = 'b';

tspan = [0 3e3]; % [ms]
options = odeset('Mass',M,'RelTol',1e-5,'MaxStep',2); 
[t,y] = ode15s(@soltis_biophysJ2010_masterODEfile,tspan,y0n,options,p);
yfinal = y(end,:)';

tArray = tArray(1:tStep);
I_Ca_store = I_Ca_store(1:tStep);
Ito = I_to_store(1,1:tStep);
Itof = I_to_store(2,1:tStep);
Itos = I_to_store(3,1:tStep);
INa = I_Na_store(1:tStep);
IK1 = I_K1_store(1:tStep);
s1 = gates(1,1:tStep);
k1 = gates(2,1:tStep);
Jserca = Jserca(1:tStep);
Iks = IKs_store(1:tStep);
Jleak = Jleak(1:tStep,:);
ICFTR = ICFTR(1:tStep);
Incx = Incx(1:tStep);
%% Calculate late current integral for INa
% % Run single simulation protocol first, load in ICs
% tStep = 1; tArray = zeros(1,1e6); I_Ca_store=zeros(1,1e6); I_to_store=zeros(3,1e6);
% I_Na_store = zeros(1,1e6); 
% tspan = [0 500];
% [t,y] = ode15s(@soltis_biophysJ2010_masterODEfile,tspan,yfinal,options,p);
% tArray = tArray(1:tStep);
% INa = I_Na_store(1:tStep);
% 
% figure(1)
% plot(tArray,INa,'k'); hold on
% 
% % Normalize by peak current
% INa_norm = INa./min(INa);
% % INa_norm = INa;
% 
% % Determine unique data indices for interpolation
% [tn,indexes,g] = unique(tArray);
% INa_unique = INa_norm(indexes);
% 
% % Roughly find 50 ms index and set new values
% ind50 = find(tn <= 50,1,'last');
% tnew = tn(ind50:end);
% INa_new = INa_unique(ind50:end);
% 
% % Interpolate data over [50 500] ms
% tinterp = [50:.1:500];
% INa_interp = interp1(tnew,INa_new,tinterp);
% 
% % Find area (integrate INa_int over tint)
% INa_integral = trapz(tinterp,INa_interp);
% normIntegral = 100*INa_integral/450;
% disp(['Current integral = ',num2str(INa_integral)])
% disp(['Normalized integral = ',num2str(normIntegral),' %'])
%% ICa Recovery from inactivation
% % Load in ICs and initialize globals
% y0n = load('yfinal_min70mV_295K_040810.dat');
% 
% global tStep tArray I_Ca_store      
% tStep = 1; tArray = zeros(1,1e6); I_Ca_store=zeros(1,1e6);
% 
% % Run preconditioning pulses to load SR
% tspan = [0 12e3];
% options = odeset('RelTol',1e-5,'MaxStep',2); 
% [to,yo] = ode15s(@soltis_biophysJ2010_masterODEfile,tspan,y0n,options,p);
% y0new = yo(end,:)';
% 
% figure(1)
% plot(to./1e3,yo(:,39)); hold on
% xlabel('Time (sec)'); ylabel('Membrane Potential (mV)');
% 
% % Plot total LCC current (ICa)
% figure(2)
% plot(tArray./1e3,I_Ca_store); hold on
% xlabel('Time (sec)'); ylabel('I_C_a (A/F)');
% 
% % Set rest intervals
% restTimes = [10,25,50,75,100,150,250,500,1e3,2e3,3e3];
% % restTimes = [10,25];
% currentMags = zeros(1,length(restTimes));
% for i=1:length(restTimes)
%     tic
%     tStep = 1; tArray = zeros(1,1e6); I_Ca_store=zeros(1,1e6);
%     p(2) = restTimes(i);
%     tend = restTimes(i) + 13e3;
%     tspan2 = [12e3 tend];
% 
%     [t,y] = ode15s(@soltis_biophysJ2010_masterODEfile,tspan2,y0new,options,p);
%     
%     tArray = tArray(1:tStep);
%     I_Ca_store = I_Ca_store(1:tStep);
%     
%     ind = find(tArray<12.5e3+restTimes(i), 1, 'last' );
%     ind2 = find(tArray<tend,1,'last');
%     
%     peakCurrent = min(I_Ca_store(ind:end));
%     residual = max(I_Ca_store(ind2:end));
%     currentMags(i) = peakCurrent-residual;
%     
%     % Plot Membrane potential
% %     figure(1)
% %     plot(t./1e3,y(:,39)); hold on
% %     xlabel('Time (sec)'); ylabel('Membrane Potential (mV)');
% % 
% %     % Plot total LCC current (ICa)
% %     figure(2)
% %     plot(tArray./1e3,I_Ca_store); hold on
% %     xlabel('Time (sec)'); ylabel('I_C_a (A/F)');
% %     
% %     ts = t./1e3;
% %     
% %     % Plot ICa Markov states
% %     Poj_m1 = 1-(y(:,60)+y(:,61)+y(:,62)+y(:,63)+y(:,64)+y(:,65));
% %     Poj_m2 = 1-(y(:,66)+y(:,67)+y(:,68)+y(:,69)+y(:,70)+y(:,71));
% %     figure(3); hold on
% %     subplot(2,1,1)
% %     plot(ts, Poj_m1,ts,y(:,61),ts,y(:,60),ts,y(:,64),ts,y(:,62),ts,y(:,65),ts,y(:,63))
% %     xlabel('Time (ms)'); ylabel('ICaj Mode 1 gates');
% %     legend('Po','C1','C2','I1Ba','I1Ca','I2Ba','I2Ca')
% %     subplot(2,1,2)
% %     plot(ts,Poj_m2,ts,y(:,67),ts,y(:,66),ts,y(:,70),ts,y(:,68),ts,y(:,71),ts,y(:,69))
% %     legend('Po','C1','C2','I1Ba','I1Ca','I2Ba','I2Ca')
% %     xlabel('Time (ms)'); ylabel('ICaj Mode 2 gates');
% 
%     disp(['Completed Simulation # ',num2str(i),' of ',num2str(length(restTimes))])
%     toc
% end
% 
% figure(4)
% % currentMags = currentMags/min(currentMags);   % Divide by min
% currentMags = currentMags/currentMags(end);     % Divide by last pulse
% plot(restTimes,currentMags,'bo-','LineWidth',2);
% xlabel('Rest interval (ms)'); ylabel('% Recovery'); ylim([0 1.1])
% assignin('base','currentMags',currentMags);
% assignin('base','restTimes',restTimes);
%% ICa facilitation
% % Need to use ICs from a long rest
% % Use vclamp protocol
% 
% % VOLTAGE CLAMP (295 K)
% y0n=load('yfinal_min90mVRest_051310.dat');
% 
% tspan = [0 cycleLength];
% options = odeset('RelTol',1e-5,'MaxStep',2); 
% 
% global tStep tArray I_Ca_store      
% 
% numPulses = 21;
% peakCurrent = zeros(1,length(numPulses));
% for i=1:numPulses
%     tic
%     tStep = 1; tArray = zeros(1,1e6); I_Ca_store=zeros(1,1e6);
%     [t,y] = ode15s(@soltis_biophysJ2010_masterODEfile,tspan,y0n,options,p);
%     y0n = y(end,:)'; % Update ICs
%     tspan = [0 cycleLength]; % Update tspan to run another pulse
% 
%     tArray = tArray(1:tStep);
%     I_Ca_store = I_Ca_store(1:tStep);
%     
%     % Plot total LCC current (ICa)
% %     figure(1)
% %     plot(tArray./1e3,I_Ca_store); hold on
% %     xlabel('Time (sec)'); ylabel('I_C_a (A/F)');
%     
%     % Extract peak from 250 ms period
%     ind = find(tArray < t(1) + 250,1,'last');
%     
%     peakCurrent(i) = min(I_Ca_store(1:ind));
%     disp(['Finished simulation # ',num2str(i),' of ',num2str(numPulses)]);
%     toc
% end
% 
% figure(2)
% plot([0:numPulses-1],peakCurrent/peakCurrent(1),'ro-','LineWidth',2); hold on
%% Ito Recovery from inactivation
% % PROTOCOL IN ECC FILE SHOULD BE 'recovery_Ito'
% 
% % Load in ICs and initialize globals
% y0=load('yfinal_min80mV_rest_295K.dat');
% % ADD IN y0 VALUES FOR NEW GATES
% y0n(1:46) = y0(1:46);
% y0n(47) = y0(4);
% y0n(48) = y0(5);
% y0n(49) = y0(6);
% y0n(50:49+15+15+15) = y0(47:end);
% 
% global tStep tArray I_to_store      
% 
% % Run a single simulation to get max peak since it has a long recovery time
% tStep = 1; tArray = zeros(1,1e6); I_to_store=zeros(3,1e6);
% options = odeset('RelTol',1e-5,'MaxStep',2); 
% tspan = [0 499];
% [t,y] = ode15s(@soltis_biophysJ2010_masterODEfile,tspan,y0n,options,p);
% tArray = tArray(1:tStep);
% I_to_store = I_to_store(:,1:tStep);
% 
% maxCurrent = max(I_to_store(1,:));
% 
% restTimes = [10,25,50,75,100,150,200,250,500,750,1e3,1.5e3,2e3,4e3,6e3,10e3,12e3];
% % restTimes = [10,25];
% currentMags = zeros(1,length(restTimes));
% 
% for i=1:length(restTimes)
%     tic
%     tStep = 1; tArray = zeros(1,1e6); I_to_store=zeros(3,1e6);
%     p(2) = restTimes(i);
%     tend = restTimes(i) + 1e3;  % Two 500 ms pulses + restTime
%     tspan2 = [0 tend];
%     
%     options = odeset('RelTol',1e-5,'MaxStep',2); 
%     [t,y] = ode15s(@soltis_biophysJ2010_masterODEfile,tspan2,y0n,options,p);
%     
%     tArray = tArray(1:tStep);
%     I_to_store = I_to_store(:,1:tStep);
%     
%     ind = find(tArray<500+restTimes(i),1,'last');
%     peakCurrent = max(I_to_store(1,ind:end));
%     currentMags(i) = peakCurrent;
%     
%     % Plot Membrane potential
% %     figure(1)
% %     plot(t./1e3,y(:,39)); hold on
% %     xlabel('Time (sec)'); ylabel('Membrane Potential (mV)');
% % 
% %     % Plot total Ito current (Ito)
% %     figure(2)
% %     plot(tArray./1e3,I_to_store(1,:)); hold on
% %     xlabel('Time (sec)'); ylabel('I_t_o (A/F)');
%     toc
% end
% 
% figure(3)
% currentMags = currentMags/maxCurrent;
% plot(restTimes,currentMags,'ro-','LineWidth',2); hold on
% xlabel('Rest interval (ms)'); ylabel('% Recovery'); ylim([0 1.1])
% assignin('base','currentMags',currentMags);
% assignin('base','restTimes',restTimes);
%% INa Recovery from inactivation
% % PROTOCOL IN ECC FILE SHOULD BE 'recovery_INa'
% 
% % Load in ICs and initialize globals
% y0=load('yfinal_min140mV_rest_295K.dat');   % Load -`140 mV ICs @ 22 C (295 K)
% % ADD IN y0 VALUES FOR NEW GATES
% y0n(1:46) = y0(1:46);
% y0n(47) = .9996;      % IC at -140 mV holding potential
% y0n(48:59) = load('yfinal_295K_-140mV_INaMarkovICs.dat'); % -140 mV ICs
% y0n(59+1:59+15+15+15) = y0(47:end);
% 
% global tStep tArray I_Na_store 
% 
% % Run a single simulation to get max peak since it has a long recovery time
% tStep = 1; tArray = zeros(1,1e6); I_Na_store=zeros(1,1e6); 
% options = odeset('RelTol',1e-5,'MaxStep',2); 
% tspan = [0 999];
% [t,y] = ode15s(@soltis_biophysJ2010_masterODEfile,tspan,y0n,options,p);
% tArray = tArray(1:tStep);
% I_Na_store = I_Na_store(1:tStep);
% INa = I_Na_store;
% 
% maxCurrent = min(INa);
% 
% restTimes = [.1,.5,1,2:2:10,20:20:160];
% % restTimes = [.1];
% currentMags = zeros(1,length(restTimes));
% 
% for i=1:length(restTimes)
%     tic
%     tStep = 1; tArray = zeros(1,1e6); I_Na_store=zeros(1,1e6);
%     I_Nal_store  = zeros(1,1e6);
%     p(2) = restTimes(i);
%     tend = restTimes(i) + 1e3 + 10;  % First 1 s pulse + 2nd 10 ms pulse + restTime
%     tspan2 = [0 tend];
%     % Reduce 'MaxStep' option for ode solver when restTime is very small (< 0.7 ms)
%     if restTimes(i) <= 1
%         maxstep = 1e-1;
%     else
%         maxstep = 2;
%     end
%     
%     options = odeset('RelTol',1e-5,'MaxStep',maxstep); 
%     [t,y] = ode15s(@soltis_biophysJ2010_masterODEfile,tspan2,y0n,options,p);
%     
%     tArray = tArray(1:tStep);
%     I_Na_store = I_Na_store(1:tStep);
%     I_Nal_store = I_Nal_store(1:tStep);
%     INa = I_Na_store+I_Nal_store;
%     
%     ind = find(tArray<1e3+restTimes(i),1,'last');
%     peakCurrent = min(INa(ind:end));
%     currentMags(i) = peakCurrent;
%     
%     % Plot Membrane potential
% %     figure(1)
% %     plot(t./1e3,y(:,39)); hold on
% %     xlabel('Time (sec)'); ylabel('Membrane Potential (mV)');
% % 
% %     % Plot total Ito current (Ito)
% %     figure(2)
% %     plot(tArray./1e3,I_Na_store); hold on
% %     xlabel('Time (sec)'); ylabel('I_N_a (A/F)');
% %     
% %     figure(3)
% %     plot(t./1e3,(y(:,1).^3),t./1e3,y(:,2),t./1e3,y(:,3)); hold on
% %     xlabel('Time (sec)'); ylabel('I_N_a gates (A/F)');
% %     legend('m3','h','j');
%     disp(['Completed simulation ', num2str(i), ' of ', num2str(length(restTimes))])
%     toc
% end
% 
% figure(1)
% currentMags = currentMags/maxCurrent;
% plot(restTimes,currentMags,'ro-','LineWidth',2); hold on
% xlabel('Rest interval (ms)'); ylabel('% Recovery'); ylim([0 1.1])
% assignin('base','currentMags',currentMags);
% assignin('base','restTimes',restTimes);
%% Plot Ca4CaM, CaMKII, and CaN timecourses
ts = t/1e3; 
tss = tArray./1e3;

CaMKIIcyt = 100.*(y(:,83+30+8)+y(:,83+30+9)+y(:,83+30+10)+y(:,83+30+11));
CaMKIIdyad = 100.*(y(:,83+8)+y(:,83+9)+y(:,83+10)+y(:,83+11));
CaMKIIsl = 100.*(y(:,83+15+8)+y(:,83+15+9)+y(:,83+15+10)+y(:,83+15+11));

LCCj_CKp = 100*y(:,83+45+2)./LCCtotDyad;
LCCsl_CKp = 100*y(:,83+45+6)./LCCtotSL;
RyR_CKp = 100*y(:,83+45+4)./RyRtot;
PLB_CKp = 100*y(:,83+45+5)./PLBtot;

% CaMKII activity
figure(1)
subplot(1,3,1); plot(ts,CaMKIIcyt,color); hold on; legend('% Cytosolic CaMKII Activity');
subplot(1,3,2); plot(ts,CaMKIIdyad,color); hold on; legend('% Dyadic CaMKII Activity');
subplot(1,3,3); plot(ts,CaMKIIsl,color); hold on; legend('% Subsarcolemmal CaKII Activity');

% Subsarcolemmal LCCp
figure(2)
plot(ts,CaMKIIsl,ts,LCCsl_CKp); legend('Sarcolemmal CaMKII Activity','SL LCCp');
ylabel('% Activity/Phosphorylation');

% CaMKII activity and dyadic phosphorylaiton targets
figure(3)
plot(ts,CaMKIIdyad,ts,LCCj_CKp,ts,RyR_CKp,ts,PLB_CKp)
xlabel('Time (sec)'); ylabel('% Activity/Phosphorylation');
legend('Dyadic CaMKII','LCCp','RyRp','PLBp');

% Membrane potential
figure(4)
plot(ts,y(:,39),color); hold on;
xlabel('Time (sec)'); ylabel('Membrane Potential (mV)');

% LCC current (ICa)
figure(5)
plot(tss,I_Ca_store,color); hold on;
xlabel('Time (sec)'); ylabel('I_C_a (A/F)');

% CaSRT & Caj
figure(6)
subplot(1,3,1)
plot(ts,y(:,30)+y(:,31),color); hold on;
xlabel('Time (sec)'); ylabel('[Ca]_S_R_T (mM)'); 
subplot(1,3,2)
plot(ts,y(:,36).*1e3,color); hold on;
xlabel('Time (sec)'); ylabel('Ca Dyad (\muM)');
subplot(1,3,3)
plot(ts,y(:,37),color); hold on;
xlabel('Time (sec)'); ylabel('Ca sl (mM)');

% Cai 
figure(7)
plot(ts,y(:,38),color); hold on;
xlabel('Time (sec)'); ylabel('[Ca]_i')

% Ito
figure(8)
plot(tss,Ito,color); hold on;
xlabel('Time (sec)'); ylabel('I_t_o'); % legend('I_t_o','I_to_f','I_to_s');
 
% INa 
figure(9)
plot(tss,INa,color); hold on;
xlabel('Time (sec)'); ylabel('I_N_a');

% PKA substrates
figure(10)
subplot(3,2,1)
plot(ts,y(:,83+45+6+23)./LCCtotBA,ts,y(:,83+45+6+24)./LCCtotBA); legend('LCCap','LCCbp')
subplot(3,2,2)
plot(ts,y(:,83+45+6+19)./PLBtotBA); legend('PLBp')
subplot(3,2,3)
plot(ts,y(:,83+45+6+25)./RyRtotBA); legend('RyRp')
subplot(3,2,4)
plot(ts,y(:,83+45+6+26)./TnItotBA); legend('TnIp')
subplot(3,2,5)
plot(ts,y(:,83+45+6+29)./IKstotBA); legend('IKsp')
subplot(3,2,6)
plot(ts,y(:,83+45+6+30)./ICFTRtotBA); legend('ICFTRp')
title('PKA-dependent phosphorylation');

% RyR states
figure(11)
plot(ts,1-sum(y(:,14:16),2),ts,y(:,14),ts,y(:,15),ts,y(:,16));
legend('RI','R','O','I')

% IKs and ICFTR
figure(12)
subplot(2,1,1)
plot(tss,Iks); legend('I_K_s')
subplot(2,1,2);
plot(tss,ICFTR); legend('I_C_F_T_R')

% [Na]
figure(13)
subplot(1,3,1); plot(ts,y(:,32),color); hold on; legend('[Na]_j');
subplot(1,3,2); plot(ts,y(:,33),color); hold on; legend('[Na]_s_l');
subplot(1,3,3); plot(ts,y(:,34),color); hold on; legend('[Na]_i');
ylabel('[Na] (mmol/L relevant compartment');

% I_NCX
figure(14)
plot(tss,Incx,color); hold on; legend('I_N_C_X');

% Integrated JRyR
figure(15)
intJryr = cumtrapz(tss,Jleak(:,1).*1e6);
plot(tss,intJryr,color); ylabel('\intJ_R_y_R (\muM)'); hold on

% RyR fluxes
figure(16)
subplot(3,1,1)
plot(tss,Jleak(:,1),color); hold on; legend('JRyR_t_o_t'); 
subplot(3,1,2)
plot(tss,Jleak(:,2),color); hold on; legend('Passive Leak');
subplot(3,1,3)
plot(tss,Jleak(:,1)-Jleak(:,2),color); hold on; legend('SR Ca release');
