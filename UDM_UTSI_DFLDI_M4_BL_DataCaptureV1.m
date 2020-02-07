%Teledyne Wavesurfer 510 UTSI/UMD M4 2pFLDI V1
%created by Jack Cobourn 20200204 using work from Andrew Ceruzi

%% RUN CONDITIONS

tunnel_section = 1; %section of tunnel FLDI is probing
DailyRunNum = 2;
CampainRunNum = 14;
PositionRunNumber = 2;
x_dist = 16.6;  % streamwise distance from nozzle/test section junction to first beam focus in inches (+-2mm)
y_dist = 2.00;  % vertical distance from tunnel floor to laser points in inch (+-2mm)
z_dist = 0;  % position of tunnel test section centerline relative to FLDI focus (cm)
obj_lens = 'f=-9mm';
Gain = 10;  %gain from amplification (if used)
NDfilter = 'no';
HLfilter = 'none';  %high or low pass Filter type
BPfilter = 'no';  % Bandpass filter
RL = 'none, using photodiodes with built in amplifier'; %[Ohms] terminating resistor
bitRes = 8;
num_diaphrams = 3;  % number of diaphrams in ludweig tube
expected_burst_pressure = '~60Psia'; % tunnel high pressure in Psia
dx1 = '300micron streamwise'; %FLDI beam separation
dx2 = 2.6e-3; %[m] 
model = 'empty tunnel'; %The model inserted in the tunnel during the run
notes = 'chA is ch2';
date = date(); %Comment out for testing on a different day


%% Connect Teledyne Ocilloscope
cd('D:\Jack\Documents\GitHub\Focused_Laser_Dif_interf')
cd('.\oscilloscope')
[VisaInterface, Wavesurfer10] = OscopeVisaConnect(); %get interface and device
VSA = VisaInterface; %simplify code with shorthands
SCP = Wavesurfer10;
ACQ = Wavesurfer10.Acquisition; %simplify code with shorthands
TRG = Wavesurfer10.Trigger;
WVE = Wavesurfer10.Waveform;
cd('..\')
%% Record CHANNELS
vmax(1) = input('Enter the ch A maximum voltage:');
vmin(1) = input('Enter the ch A minimum voltage:');
vmax(2) = input('Enter the ch B maximum voltage:');
vmin(2) = input('Enter the ch B minimum voltage:');
filtenab = 0; %input('High-Pass filter enabled? (Y/N): ');
vavg = (vmax+vmin)./2.0;
%vopts = [10 20 50 100 200 500 1000 2000 5000 10000 20000]; % voltage range options in millivolts
vsep = vmax - vavg; % range of instrument for both channels


ACoupling = query(VisaInterface,'C3:COUPLING?');
BCoupling = query(VisaInterface,'C2:COUPLING?');
AOffset = eval(query(VisaInterface,'C3:OFFSET?'));
BOffset = eval(query(VisaInterface,'C2:OFFSET?'));
AScale = eval(query(VisaInterface,'C3:VOLT_DIV?'))*10;
BScale = eval(query(VisaInterface,'C2:VOLT_DIV?'))*10;

%% GET TIMEBASE
recordtime = get(ACQ,'TimeBase')*10; %seconds
set(ACQ,'Delay',-get(ACQ,'TimeBase').*10./2) %set delay so that nothing from before trigger is saved
numSamples = eval(query(VisaInterface,'MEMORY_SIZE?')); %%% Add to driver properly. also you always get 2 extra points
Fs = numSamples./recordtime; %Hz

%% PAUSE FOR USER INPUT
disp('Paused: Please press any key to collect noise data');
pause;

%% CAPTURE NOISE DATA
 set(ACQ, 'Control', 'single'); %set to single with matlab driver
    invoke(TRG, 'trigger');
    %Note channel A is bigger range.
    [chA_noise, timeMs, chA_noise_info] = invoke(WVE,'readwaveform','C2',true); %grab cha, with the time, and scled to doubles
    [chB_noise, ~, chB_noise_info] = invoke(WVE,'readwaveform','C3',true); %no need to grab time again
    
    
%% PAUSE FOR USER INPUT
disp('Paused: Please press any key to set up for run.');
pause;

%% Prepare to CAPTURE RUN DATA
 set(ACQ, 'Control', 'single'); %set to single with matlab driver
    
 
 %% PAUSE FOR USER INPUT
disp('Paused: Please press any key to collect data after run.');
pause;
%maybe should poll scope for trigger but easier to just click at end of run

%% Collect RUN DATA
    [chA_run, ~, chA_run_info] = invoke(WVE,'readwaveform','C2',true);
    [chB_run, ~, chB_run_info] = invoke(WVE,'readwaveform','C3',true);
    

%% Plot Signal and Noise and Histogram
disp('Plotting Signal')

figure(1)
clf
plot(timeMs,chA_run)
hold on
plot(timeMs,chB_run)
plot(timeMs, chA_noise, 'k')
plot(timeMs, chB_noise, 'k:')
title('Channel A & B run Signal and flow-off Noise');
xlabel('Time (ms)');
ylabel('Voltage (mV)');
% ylim([-vopts(vrange+1) vopts(vrange+1)]);
hold off

% figure(4)
% clf
% plot(timeMs,chB_run);
% hold on
% plot(timeMs, chB_noise, 'k');
% title('Channel B (Trigger) run Signal and Noise');
% xlabel('Time (ms)');
% ylabel('Voltage (mV)');1;
% % ylim([-vopts(vrange+1) vopts(vrange+1)]);
% hold off

% figure(2)
% clf
% hist(chA_noise,50);
% title('Channel A Histogram Noise');
% ylabel('Counts');
% xlabel('Voltage (mV)');

%% Plot Some Spectra

%Note: this should only be run on the steady-state part of the run
%Another note: run this for the noise too so they can be compared

ham=100; %number of divisions for hamming windows
start = 0.1; %[s] trim window start
stop=0.2; %[s] trim window stop
%compute spectra for channel A
chA_run_trim = (chA_run(fix((Fs*start)):fix((Fs*stop)-1)));
% chA_noise_trim = (chA_noise((sps*start)+1:(sps*stop)-1));
[PSDa, f] = pwelch(chA_run_trim-mean(chA_run_trim),ceil(length(chA_run_trim)/ham),[],[],ceil(Fs));
[PSDna, fn] = pwelch((chA_noise-mean(chA_noise)),numSamples/ham,[],[],Fs);
%compute spectra for channel B
chB_run_trim = (chB_run(fix((Fs*start)):fix((Fs*stop))-1));
% chB_noise_trim = (chB_noise((sps*start)+1:(sps*stop)-1));
[PSDb, f] = pwelch(chB_run_trim-mean(chB_run_trim),ceil(length(chB_run_trim)/ham),[],[],ceil(Fs));
[PSDnb, fn] = pwelch((chB_noise-mean(chB_noise)),numSamples/ham,[],[],Fs);

fig3 = figure(3)
clf
loglog(f,PSDa,'linewidth',2)
hold on;
loglog(f,PSDb,'linewidth',2)
loglog(fn,PSDna,'k')
loglog(fn,PSDnb,'k:')

title('Channel A & B spectra')
ylabel('PSD [mV^2/Hz]')
xlabel('Frequency [Hz]')
xlim([0 Fs/2])

legend('Steady Flow: A','Steady Flow: B','Flow Off')
hold off;

%% compute cross correlation
[cross_cor,lag] = xcorr(chA_run_trim-mean(chA_run_trim),-(chB_run_trim-mean(chB_run_trim)), 'coeff');
%find max of xcorr and index
[cc,I] = max(cross_cor);
%check if xcorr peak is high enough?
if cc>0.3
    lagDiff = lag(I);
    %3 pt gaussian interpolation, to get a more exact velocity
    peakOffset = ((log(cross_cor(I-1)))-(log(cross_cor(I+1))))/(2*((log(cross_cor(I-1)))+(log(cross_cor(I+1)))-2*(log(cc))));
    %find delta t
    DT = abs((lagDiff+peakOffset)/Fs);
else
    lagDiff = NaN;
    DT = NaN;
end
%disturbance velocity, [m/s]
Uc = (dx2)/DT

figure(4)
clf
plot(lag,cross_cor)

 %% PAUSE FOR USER INPUT
disp('Paused: Please press any key to collect 2nd noise data at high P.');
disp('Will overwrite run data on scope, make sure its done collecting');
pause;

%% Collect Noise Data2
set(ACQ, 'Control', 'single'); %set to single with matlab driver
    invoke(TRG, 'trigger');
    [chA_noise2, ~, chA_noise_info2] = invoke(WVE,'readwaveform','C2',true);
    [chB_noise2, ~, chB_noise_info2] = invoke(WVE,'readwaveform','C3',true);
    
%% DEVICE DISCONNECTION
disp('Disconnecting')
% Disconnect device object from hardware.
disconnect(Wavesurfer10)

%% SAVE
disp('Saving data...')
% %save the variables
% save F:\FLDI_UMD\200204\UTSI_M4_2ptFLDI_y_2_25_inch_02.mat timeMs
% recordtime chA_run chA_noise chB_run chB_noise chA_run_info chA_noise_info chB_run_info chB_noise_info Fs numSamples bitRes vmax vmin vavg dx1 dx2 tunnel_section x_dist y_dist z_dist Gain HLfilter RL num_diaphrams expected_burst_pressure model NDfilter BPfilter obj_lens notes date DailyRunNum CampainRunNum PositionRunNumber
folderstring='F:\FLDI_UMD\20200206';
savestring='UTSI_M4_2ptFLDI_y_2_00_inch_02';
%mkdir(folderstring);
cd(folderstring)
save([folderstring '\' savestring '.mat'],...
    'chA_run', 'chA_noise','chA_noise2', 'chB_run', 'chB_noise','chB_noise2',...
    'chA_run_info', 'chA_noise_info','chA_noise_info2', 'chB_run_info', 'chB_noise_info','chB_noise_info2',...
    'timeMs','Fs', 'numSamples', 'bitRes', 'vmax', 'vmin', 'vavg', 'dx1', 'dx2', 'tunnel_section', 'x_dist', 'y_dist', 'z_dist',...
    'Gain', 'HLfilter', 'RL', 'num_diaphrams', 'expected_burst_pressure','model', 'NDfilter', 'BPfilter', 'obj_lens',...
    'notes', 'date', 'DailyRunNum', 'CampainRunNum', 'PositionRunNumber','-v7.3','-nocompression');
savefig(fig3,[savestring '.fig']);
saveas(fig3,savestring,'jpeg')
disp('DONE!')