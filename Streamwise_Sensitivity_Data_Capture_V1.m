%% Streamwise_Sensitivity_Data_Capture_V1
% Jack Cobourn cw 20200519
% Using:
%Teledyne Wavesurfer 510, M4, 2pFLDI #1
clear all
close all

%% Data
date = date(); %Comment out for testing on a different day
ISOdate = yyyymmdd(datetime(date));

%position
tunnel_section = 1; %section of tunnel FLDI is probing

z_pos = [6]; %input('enter current position from centerline in 1/64 inches = , positive for catch side');
y_FLDI = [10+28/64];
y_jet = [10+40/64];
x_FLDI = [6+1/8];
X_jet = [4]; %streamwise from FLID

%save Location
folderstring='C:\Users\jack\Documents\FLDI FastWrite';
Runs = length(dir([folderstring filesep '**\*.mat']));

savestring=strrep(sprintf('Sensitivity_2ptFLDI_RunNum=%u_z=%.2f',Runs+1,z_pos),'.',',');

%FLDI Config

dx = [567.9, 605.1 , 909.8, 947.2];
dy = [471.1, 471.2, 473.4, 473.1];
dx1a = 0.00738*(dx(2)-dx(1))/1000; %FLDI beam separation
dy1a = 0.00738*(dy(2)-dy(1))/1000;
dx1b =  0.00738*(dx(4)-dx(3))/1000;
dy1b = 0.00738*(dy(4)-dy(3))/1000;
dx2 =  0.00738*(mean(dx([3,4]))-mean(dx([1,2])))/1000; %[m] 
dy2 = 0.00738*(mean(dy([3,4]))-mean(dy([1,2])))/1000; %[m] 
d2 = sqrt(dx2.^2+dy2.^2);

obj_lens = 'f=-9mm';
Gain = 0;  %gain from amplification (if used)
NDfilter = 'no';
HLfilter = 'none';  %high or low pass Filter type
BPfilter = 'no';  % Bandpass filter
RL = '50'; %[Ohms] terminating resistor
bitRes = 8;

notes = {'chA is ch2, is upstream beam and downstream transducer';'set Both Chanels to 1MOhms'};

%Chanels
vmax(1) = 142; %input('Enter the ch A maximum voltage:');
vmin(1) = 26.2;%input('Enter the ch A minimum voltage:');
vmax(2) = 134;%input('Enter the ch B maximum voltage:');
vmin(2) = 30.9;%input('Enter the ch B minimum voltage:');
filtenab = 0; %input('High-Pass filter enabled? (Y/N): ');
vavg = (vmax+vmin)./2.0;
vsep = vmax - vavg; % range of instrument for both channels


%% Connect Teledyne Ocilloscope
cd('C:\Users\Jack\Documents\GitHub\Focused_Laser_Dif_interf')
cd('.\oscilloscope')
[VisaInterface, Wavesurfer10] = OscopeVisaConnect(); %get interface and device
VSA = VisaInterface; %simplify code with shorthands
SCP = Wavesurfer10;
ACQ = Wavesurfer10.Acquisition; %simplify code with shorthands
TRG = Wavesurfer10.Trigger;
WVE = Wavesurfer10.Waveform;
cd('..\')

% %query(VSA,"VBS? 'return=app.Acquisition.Horizontal.SampleMode'")
% if ~strcmpi(query(VSA,"VBS? 'return=app.Acquisition.Horizontal.SampleMode'"),'RealTime\n')
%     fprintf(VSA, "VBS 'app.Acquisition.Horizontal.SampleMode = RealTime'");
%     error('Incorrect Sample Mode')
% end
%% Record CHANNELS
ACoupling = query(VSA,'C3:COUPLING?');
BCoupling = query(VSA,'C2:COUPLING?');
AOffset = eval(query(VSA,'C3:OFFSET?'));
BOffset = eval(query(VSA,'C2:OFFSET?'));
AScale = eval(query(VSA,'C3:VOLT_DIV?'))*10;
BScale = eval(query(VSA,'C2:VOLT_DIV?'))*10;

%% GET TIMEBASE
recordtime = get(ACQ,'TimeBase')*10; %seconds
set(ACQ,'Delay',-get(ACQ,'TimeBase').*10./2) %set delay so that nothing from before trigger is saved
numSamples = eval(query(VSA,'MEMORY_SIZE?')); %%% Add to driver properly. also you always get 2 extra points
Fs = numSamples./recordtime; %Hz

%% CAPTURE Noise
set(ACQ, 'Control', 'single'); %set to single with matlab driver
    %doesnt seem to work just dont forget to do this by hand.
    %fprintf(VisaInterface,"VBS 'app.Acquisition.Horizontal.SampleMode=RealTime'")
    invoke(TRG, 'trigger');
    
    [chA_noise, ~, chA_noise_info] = invoke(WVE,'readwaveform','C2',true);
    [chB_noise, ~, chB_noise_info] = invoke(WVE,'readwaveform','C3',true);

%% Turn On Air
Pre_Pressure = input('read pressure with air off in millivolts'); %mV
disp('Turn on air');
Pressure = input('read pressure in volts'); %V
%% CAPTURE DATA
set(ACQ, 'Control', 'single'); %set to single with matlab driver
    %doesnt seem to work just dont forget to do this by hand.
    %fprintf(VisaInterface,"VBS 'app.Acquisition.Horizontal.SampleMode=RealTime'")
    invoke(TRG, 'trigger');
    
    [chA_run, timeMs, chA_info] = invoke(WVE,'readwaveform','C2',true);
    [chB_run, ~, chB_info] = invoke(WVE,'readwaveform','C3',true);
    
%% Plot Signal
disp('Plotting Signal')
skipplot = round(linspace(1,length(timeMs),1000));
figure(1)
clf
plot(timeMs(skipplot),chA_run(skipplot),timeMs(skipplot),chB_run(skipplot),timeMs(skipplot), chA_noise(skipplot), 'k',timeMs(skipplot), chB_noise(skipplot), 'k:')
title('Channel A & B run Signal');
xlabel('Time (ms)');
ylabel('Voltage (mV)');
% ylim([-vopts(vrange+1) vopts(vrange+1)]);
hold off

%% Plot Some Spectra

%Note: this should only be run on the steady-state part of the run
%Another note: run this for the noise too so they can be compared

start = 0.1;%input('input start time'); %[s] trim window start
stop= .9;% input('input stop time'); %[s] trim window stop
%compute spectra for channel A
chA_run_trim = (chA_run(fix((Fs*start)):fix((Fs*stop)-1)));
chA_noise_trim = (chA_noise(fix((Fs*start)):fix((Fs*stop)-1)));

[PSDa, f] = pwelch(chA_run_trim-mean(chA_run_trim),...
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(Fs));
[PSDna, fn] =  pwelch(chA_noise_trim-mean(chA_noise_trim),...
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(Fs));

%compute spectra for channel B
chB_run_trim = (chB_run(fix((Fs*start)):fix((Fs*stop))-1));
 chB_noise_trim = (chB_noise(fix((Fs*start)):fix((Fs*stop))-1));
 
[PSDb, f] = pwelch(chB_run_trim-mean(chB_run_trim),...
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(Fs));
[PSDnb, fn] = pwelch(chB_noise_trim-mean(chB_noise_trim),...
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(Fs));
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

%% DEVICE DISCONNECTION
disp('Disconnecting')
% Disconnect device object from hardware.
disconnect(Wavesurfer10)
clear VisaInterface Wavesurfer10 VSA SCP ACQ TRG WVE
%% SAVE
%check if files exist
disp('Saving data...')
if ~isfolder(folderstring)
    mkdir(folderstring);
end
cd(folderstring)
if isfile([savestring '.fig']) || isfile([savestring '.jpg']) || isfile([folderstring '\' savestring '.mat'])
    error('File Exists: cannot overwrite')
end
%save fig 3
savefig(fig3,[savestring '.fig']);
saveas(fig3,[savestring '.jpg'],'jpeg');
clear fig3

%save matfile
save([folderstring '\' savestring '.mat'],'-v7.3');

disp('DONE!')