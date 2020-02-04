[VisaInterface, Wavesurfer10] = OscopeVisaConnect(); %get interface and device
VSA = VisaInterface; %simplify code with shorthands
SCP = Wavesurfer10;
ACQ = SCP.Acquisition; %simplify code with shorthands
TRG = SCP.Trigger;
WVE = SCP.Waveform;
MEA = SCP.Measurement;

delimiter = ',; ';
form = '%*s %f %*s';


%take 100ms samples
set(ACQ,'Timebase',.1/10); %100ms total
fwrite(VisaInterface,'MSIZ 1e+5') %100 kS => 1MS/s
for i=1:500
    fwrite(VisaInterface,'ARM')
    fwrite(VisaInterface,'FRTR')

    STD = cell2mat(textscan(query(VisaInterface,"C3:PAVA? SDEV"),...
        form,'delimiter',delimiter));
    
    STDVEC1(i) = STD;
    STDSUM1(i) = sum(STDVEC1(1:i))/i;
    
end

%take 1s samples
set(ACQ,'Timebase',1/10); %1 sec total
fwrite(VisaInterface,'MSIZ 1e+6') %1 MS => 1MS/s
for i=1:500
    fwrite(VisaInterface,'ARM')
    fwrite(VisaInterface,'FRTR')

    STD = cell2mat(textscan(query(VisaInterface,"C3:PAVA? SDEV"),...
        form,'delimiter',delimiter));
    
    STDVEC2(i) = STD;
    STDSUM2(i) = sum(STDVEC2(1:i))/i;
    
end

%take 10s samples
set(ACQ,'Timebase',10/10); %10 sec total
fwrite(VisaInterface,'MSIZ 1e+7') %10 MS => 1MS/s
for i=1:500
    fwrite(VisaInterface,'ARM')
    fwrite(VisaInterface,'FRTR')
    pause(1)
    STD = cell2mat(textscan(query(VisaInterface,"C3:PAVA? SDEV"),...
        form,'delimiter',delimiter));
    
    STDVEC3(i) = STD;
    STDSUM3(i) = sum(STDVEC3(1:i))/i;
    
end

%plot
figure(5)
hold on
AVE1 = semilogy(STDSUM1,'k');
points1 = semilogy(STDVEC1,'k*');

AVE2 = semilogy(STDSUM2,'b');
points2 = semilogy(STDVEC2,'b*');

AVE3 = plot(STDSUM3,'c');
points3 = plot(STDVEC3,'c*');

grid on
legend([AVE1 AVE2 AVE3],{'100ms at 1 MHZ','1s at 1 MHz', '10s at 1 MHz'})
xlabel('Measurement number')
ylabel('Standard Deviation')
title({'Wind Off Convergence Study'; 'Cumulative average of Standard Deviation'})
