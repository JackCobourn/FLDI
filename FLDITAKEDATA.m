%% FLDI TAKE DATA
%script to take data at one x/d location
    %everything in mm
%input x
x_loc= 14.6.*10; 
%input d
d=(3/16).*25.4;
JetPos = x_loc./d;
%% Connect scope
[VisaInterface, Wavesurfer510] = OscopeVisaConnect_V2(); %get interface and device
VSA = VisaInterface; %simplify code with shorthands
SCP = Wavesurfer510;
ACQ = Wavesurfer510.Acquisition; %simplify code with shorthands
TRG = Wavesurfer510.Trigger;
WVE = Wavesurfer510.Waveform;
%% Connect and home motor
a=motor.listdevices; %find all attached motors. '27252541' is the y axis
MotorY=motor();        % create a motor object
connect(MotorY,a{1}); %connect. please don't have this step fail
MotorY.home(); %home is at top
MotorY.moveto(50);

%% Position loop

PosStep=0.5; %set step size
HomePos = 0; %function of motor home and cordinate system
EndPos = 50;
CurrentPosition = EndPos;

waves2take = 20;

FLDI_Data = struct('Jetposition',JetPos,'X_loc',x_loc,'d',d);
%take max data
    set(Wavesurfer510.Acquisition(1), 'Control', 'single'); %set to single with matlab driver
    %set(Waveform,'Precision','int16'); %set prcision with matalab driver
    %Already set in new drier
    invoke(fprintf(VisaInterface, "FRTR"); %instead write a low level force trig directly to machine
    [Y,X,YUNIT,XUNIT,HEADER] = invoke(WVE,'readwaveform','channel3');
    FLDI_Data.Max16 = struct('Y',Y,'X',X,'YUNIT',YUNIT,'XUNIT',XUNIT,'Header',HEADER);
    
    set(WVE,'Precision','int8');
    [Y,X,YUNIT,XUNIT,HEADER] = invoke(WVE,'readwaveform','channel3');
    FLDI_Data.Max8 = struct('Y',Y,'X',X,'YUNIT',YUNIT,'XUNIT',XUNIT,'Header',HEADER);


%take min data
    set(Wavesurfer510.Acquisition(1), 'Control', 'single'); %set to single with matlab driver
    set(WVE,'Precision','int16'); %set prcision with matalab driver
    fprintf(VisaInterface, "FRTR"); %instead write a low level force trig directly to machine
    [Y,X,YUNIT,XUNIT,HEADER] = invoke(WVE,'readwaveform','channel3');
    FLDI_Data.Min16 = struct('Y',Y,'X',X,'YUNIT',YUNIT,'XUNIT',XUNIT,'Header',HEADER);
    
    set(WVE,'Precision','int8');
    [Y,X,YUNIT,XUNIT,HEADER] = invoke(WVE,'readwaveform','channel3');
    FLDI_Data.Min8 = struct('Y',Y,'X',X,'YUNIT',YUNIT,'XUNIT',XUNIT,'Header',HEADER);


%set up sweep
FLDI_Data.YSweep = struct();

tic
while CurrentPosition>=HomePos
    %known at start
    CurrentStep = 1+(EndPos-CurrentPosition)./PosStep;
    FLDI_Data.YSweep(CurrentStep).YPosition = CurrentPosition;
    FLDI_Data.YSweep(CurrentStep).TrueYPosition = MotorY.position;
    FLDI_Data.YSweep(CurrentStep).Data = struct();
    
    %DataLoop
    for i = 1:1:waves2take
        set(Wavesurfer510.Acquisition(1), 'Control', 'single'); %set to single with matlab driver
        %invoke(Wavesurfer510.Trigger,'trigger'); still does not work.
        %Double arm command does not trigger as manual says
        fprintf(VisaInterface, "FRTR"); %instead write a low level force trig directly to machine
        
        set(WVE,'Precision','int16'); %set prcision with matalab driver
        [Y,X,YUNIT,XUNIT,HEADER] = invoke(WVE,'readwaveform','channel3');
        FLDI_Data.YSweep(CurrentStep).Data(i).Wave16 = struct('Y',Y,'X',X,'YUNIT',YUNIT,'XUNIT',XUNIT,'Header',HEADER);
        
        
        set(WVE,'Precision','int8'); %set to int8 to do back up data aqq
        [Y,X,YUNIT,XUNIT,HEADER] = invoke(WVE,'readwaveform','channel3');
        FLDI_Data.YSweep(CurrentStep).Data(i).Wave8 = struct('Y',Y,'X',X,'YUNIT',YUNIT,'XUNIT',XUNIT,'Header',HEADER);

    end
    
    %Calculate move
    NextPos = CurrentPosition-PosStep;
    %motor Move
    if NextPos >= HomePos
        MotorY.moveto(NextPos);    
    end
    
    %iterate counter
    CurrentPosition = NextPos;
end
toc
%SaveStr = [sprintf('Jet_FLDI_Data_at_%.0f_Dia_(%ux%ux%u)',JetPos,EndPos./PosStep,waves2take,length(FLDI_Data.YSweep(CurrentStep).Data(i).Wave16(Y))),'.mat'];
SaveStr='20191027_Partial_FLDI_Data';
save(SaveStr,'FLDI_Data','-v7.3','-nocompression')