%% Photo Diode Responce Analysis
% Jack Cobourn cw 20201121
% Using Teledyne Lacroy Wavesurfer 10 and Thorlabs PDA10A2 and PDA36A2

%% Cases and load data
casenum = 3;

switch casenum
    
    case 1 %Unamplified Diode (PDA10A2)
        casestring = 'C2pdresponse100000.dat';
      
    case 2 %Amplified Diode No Gain (PDA36A2)
      casestring = 'C2pdresponse2-nogain00000.dat';
      
    case 3 %Amplified Diode 10 dB Gain (PDA36A2)
      casestring = 'C2pdresponse2-10db00000.dat';
end

Data = readmatrix(casestring);

%% Process Data
fig1 = figure(1);
plot(10^9.*Data(:,1),1000*Data(:,2),'k','LineWidth',1.8)
fig1.Children.FontSize = 20;
ylabel('Voltage [mV]','FontSize',28);
xlabel('Time [ns]','FontSize',28);
bitdepthfilled = numel(unique(Data(:,2)));

fs = 1e10;

%Grab start of falling step respoinse, after 'steady' peak, and end of response.
[x,~] = ginput(2);
xvals(1) = find(10^9.*Data(:,1) > x(1),1);
xvals(2) = find(10^9.*Data(:,1) > x(2),1);
xlim(10^9.*[Data(round(xvals(1)-0.15.*xvals(1))),Data(round(xvals(2)+0.10.*xvals(2)))])
%hold on
%plot([-100,0,1*10^-9,100],[0,0,1,0],'r','LineWidth',1.8)
%plot(Data(xvals(1):xvals(2),1000*Data(xvals(1):xvals(2),2),'k')

%Pwelch
% L = length(Data(xvals(1):xvals(2)));
% two = floor(log2(L));
% Window = hamming(2^(two-3));%ones(1,2^(two-1));
% noverlap = length(Window)./2;
% Nfft = 2^(two);
% 
% [spectrum,f] = pwelch(Data(xvals(1):xvals(2),2),Window,noverlap,Nfft,fs,'power');
% 
% fig2 = figure(2);
% semilogx(1+f,mag2db(spectrum))
% xlim([0,inf])
N = 2^(7+nextpow2(length(Data(xvals(1):xvals(2),2))));
H = 1/N*fftshift(fft(Data(xvals(1):xvals(2),2),N));
df=fs/N; %frequency resolution
sampleIndex = -N/2:N/2-1; %ordered index for FFT plot
f=sampleIndex*df; %x-axis index converted to ordered frequencies

Mag = abs(H);
phase = atan2d(imag(H),real(H));

H2=H;%store the FFT results in another array
%detect noise (very small numbers (eps)) and ignore them
threshold = max(abs(H))/10000; %tolerance threshold
H2(abs(H)<threshold) = 0; %maskout values that are below the threshold
phase=atan2(imag(H2),real(H2))*180/pi; %phase information


fig2 = figure(2);
t = tiledlayout(2,1);

%mag
A = nexttile;
semilogx(f,mag2db(Mag)-max(mag2db(Mag)),'k','LineWidth',1.8);
A.FontSize = 20;
ylabel('Magnitude [dB]','FontSize',28)
xlabel('Frequency [Hz]','FontSize',28)
ylim([min(mag2db(Mag))-3-max(mag2db(Mag)),3])
xlim([-inf,10^8])
ylim([-50,3])



%phase
B = nexttile;
semilogx(f,phase,'k','LineWidth',1.8);
B.FontSize = 20;
ylim([-90,inf])
ylabel('Phase','FontSize',28)
xlabel('Frequency [Hz]','FontSize',28)
xlim([-inf,10^8])
ylim([-90,15])



