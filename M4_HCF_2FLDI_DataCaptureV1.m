%Teledyne Wavesurfer 510 M4 2pFLDI V1
%edited for HCF data capture by Jack Cobourn 20200319 based on original code by Andrew Ceruzi

%% RUN CONDITIONS
clear ISOdate date
%date
date = date(); %Comment out for testing on a different day
ISOdate = yyyymmdd(datetime(date));

%run numbers
DailyRunNum = 1;
CampainRunNum = 14;
PositionRunNumber = 1;

%position
tunnel_section = 1; %section of tunnel FLDI is probing
% x_dist = 25.4*(6+6/32);  % streamwise distance model mount to forward laser position (+-2mm)
% rail_pos = 25.4*9.5; % position of rail measured on the stands. should be laser position from back window.
% y_dist = (10+59/64)*25.4;  % vertical distance from tunnel floor to laser points in inch (+-2mm)
 z_dist = 0;  % position of tunnel test section centerline relative to FLDI focus (cm)

xinside = NaN;
xoutside = [6+4.5/64 6+5/64];
yinside = NaN;
youtside = [4+44/64 4+44/64];

%FLDI Config
dx = [453 490 870 908];
dy = [350 350 162 162];
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

%tunnel Conditions
num_diaphrams = 1;  % number of diaphrams in ludweig tube
expected_burst_pressure = '~17Psia'; % tunnel high pressure in Psia

model = 'HCF'; %The model inserted in the tunnel during the run
notes = {'chA is ch2, is upstream beam and downstream transducer';'set Both Chanels to 1MOhms'};


%% Connect Teledyne Ocilloscope
cd('D:\Jack\Documents\GitHub\Focused_Laser_Dif_interf')
cd('.\oscilloscope')
[VisaInterface, Wavesurfer10] = OscopeVisaConnect(); %get interface and device
VSA = VisaInterface; %simplify code with shorthands
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
vmax(1) = input('Enter the ch A maximum voltage:');
vmin(1) = input('Enter the ch A minimum voltage:');
vmax(2) = input('Enter the ch B maximum voltage:');
vmin(2) = input('Enter the ch B minimum voltage:');
filtenab = 0; %input('High-Pass filter enabled? (Y/N): ');
vavg = (vmax+vmin)./2.0;
%vopts = [10 20 50 100 200 500 1000 2000 5000 10000 20000]; % voltage range options in millivolts
vsep = vmax - vavg; % range of instrument for both channels


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

%% PAUSE FOR USER INPUT
disp('Paused: Please press any key to collect noise data');
pause;

%% CAPTURE NOISE DATA
 set(ACQ, 'Control', 'single'); %set to single with matlab driver
    %doesnt seem to work just dont forget to do this by hand.
    %fprintf(VisaInterface,"VBS 'app.Acquisition.Horizontal.SampleMode=RealTime'")
    invoke(TRG, 'trigger');
    
    %Note channel A is bigger range.
    [chA_noise, ~, chA_noise_info] = invoke(WVE,'readwaveform','C2',true); %grab cha, with the time, and scled to doubles
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
    [chA_run, timeMs, chA_run_info] = invoke(WVE,'readwaveform','C2',true);
    [chB_run, ~, chB_run_info] = invoke(WVE,'readwaveform','C3',true);
    

%% Plot Signal
disp('Plotting Signal')

figure(1)
clf
plot(timeMs,chA_run,timeMs,chB_run,timeMs, chA_noise, 'k',timeMs, chB_noise, 'k:')
title('Channel A & B run Signal and flow-off Noise');
xlabel('Time (ms)');
ylabel('Voltage (mV)');
% ylim([-vopts(vrange+1) vopts(vrange+1)]);
hold off

%% PAUSE FOR USER INPUT
disp('Paused: Please press any key to collect 2nd noise data at high P.');
disp('Will overwrite run data on scope, make sure its done collecting');
pause;

%% Collect Noise Data2
set(ACQ, 'Control', 'single'); %set to single with matlab driver
    invoke(TRG, 'trigger');
    [chA_noise2, ~, chA_noise_info2] = invoke(WVE,'readwaveform','C2',true);
    [chB_noise2, ~, chB_noise_info2] = invoke(WVE,'readwaveform','C3',true);
    
%% Plot Some Spectra

%Note: this should only be run on the steady-state part of the run
%Another note: run this for the noise too so they can be compared

start = input('input start time'); %[s] trim window start
stop=input('input stop time'); %[s] trim window stop
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
Uc = (d2)/DT

figure(4)
clf
plot(lag,cross_cor)

 
    
%% DEVICE DISCONNECTION
disp('Disconnecting')
% Disconnect device object from hardware.
disconnect(Wavesurfer10)
clear VisaInterface Wavesurfer10 VSA SCP ACQ TRG WVE
%% SAVE
disp('Saving data...')
folderstring=sprintf('F:\\FLDI_HCF\\Data\\%d',ISOdate);
savestring=sprintf('M4_2ptFLDI_HCF_CRun_%u_DRun_%u',CampainRunNum,DailyRunNum);
if ~isfolder(folderstring)
    mkdir(folderstring);
end
cd(folderstring)
if isfile([savestring '.fig'])
    error('FileExists')
elseif isfile([savestring '.jpg'])
    error('FileExists')
elseif isfile([folderstring '\' savestring '.mat'])
    error('FileExists')
end
savefig(fig3,[savestring '.fig']);
saveas(fig3,savestring,'jpeg');
clear fig3
save([folderstring '\' savestring '.mat'],'-v7.3');
fid = fopen([folderstring '\' savestring '.txt'],'a');
fprintf(fid,'Uc = %.4f',Uc);
fclose(fid);
disp('DONE!')