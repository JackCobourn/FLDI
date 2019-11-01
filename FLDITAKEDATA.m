%% FLDI TAKE DATA
%script to take data at one x/d location
    %everything in mm
%input x
x_loc= 14.6.*10; 
%input d
d=(3/16).*25.4;
JetPos = x_loc./d;
%% Connect scope
[VisaInterface, Wavesurfer510] = OscopeVisaConnect(); %get interface and device
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
Y = int8(1);
PosStep=0.5; %set step size
HomePos = 0; %function of motor home and cordinate system
EndPos = 50;
CurrentPosition = EndPos;
waves2take = 20;

CurrentStep = 0; % intitialize

FLDI_Data = struct('Jetposition',JetPos,'X_loc',x_loc,'d',d);
%take max data
    set(ACQ, 'Control', 'single'); %set to single with matlab driver
    %set(Waveform,'Precision','int16'); %set prcision with matalab driver
    %Already set in new drier
    invoke(TRG, 'trigger'); 
    [Y,~, info] = invoke(WVE,'readwaveform','C3',false);
    FLDI_Data.Max = struct('Y',Y,'INFO',info);
    
%take min data
    set(ACQ, 'Control', 'single'); %set to single with matlab driver
    %set(Waveform,'Precision','int16'); %set prcision with matalab driver
    %Already set in new drier
    invoke(TRG, 'trigger'); 
    [Y,~, info] = invoke(WVE,'readwaveform','C3',false);
    FLDI_Data.Min = struct('Y',Y,'INFO',info);

%take Zero Data
    set(ACQ, 'Control', 'single'); %set to single with matlab driver
    %set(Waveform,'Precision','int16'); %set prcision with matalab driver
    %Already set in new drier
    invoke(TRG, 'trigger'); 
    [Y,~, info] = invoke(WVE,'readwaveform','C3',false);
    FLDI_Data.Zero = struct('Y',Y,'INFO',info);

    
%set up sweep
FLDI_Data.YSweep = struct();
%%%
%initialize waitbar
f = uifigure;
d = uiprogressdlg(f,'Title','Please Wait',...
    'Message','Begining Data ACQ');

tic
d.Title = 'Running data ACQ';
while CurrentPosition>=HomePos
    %known at start
    CurrentStep = 1+(EndPos-CurrentPosition)./PosStep;
    d.Message=sprintf('moving to position %u of 100',CurrentStep);
    FLDI_Data.YSweep(CurrentStep).YPosition = CurrentPosition;
    FLDI_Data.YSweep(CurrentStep).TrueYPosition = MotorY.position;
    FLDI_Data.YSweep(CurrentStep).Data = struct();
    %fprintf('Current Position is %0.3f',MotorY.position)
    %DataLoop
    for i = 1:1:waves2take
        try
            d.Message=sprintf('Running loop %u at position %u of 100',i,CurrentStep);
            d.Value = (20*(CurrentStep-1)+i)./2000;
        set(ACQ, 'Control', 'single'); %set to single with matlab driver
        invoke(TRG, 'trigger'); 
        [Y,~, info] = invoke(WVE,'readwaveform','C3');
        FLDI_Data.YSweep(CurrentStep).Data(waves2take) = struct('Y',Y,'INFO',info);
        
        catch %if it doesnt work try again
        set(ACQ, 'Control', 'single'); %set to single with matlab driver
        pause(1)
        invoke(TRG, 'trigger');
        pause(1)
        [Y,~, info] = invoke(WVE,'readwaveform','C3');
        pause(1)
        FLDI_Data.YSweep(CurrentStep).Data(waves2take) = struct('Y',Y,'INFO',info);
        end
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
close(d)
%SaveStr = [sprintf('Jet_FLDI_Data_at_%.0f_Dia_(%ux%ux%u)',JetPos,EndPos./PosStep,waves2take,length(FLDI_Data.YSweep(CurrentStep).Data(i).Wave16(Y))),'.mat'];
SaveStr='20191027_Partial_FLDI_Data';
save(SaveStr,'FLDI_Data','-v7.3','-nocompression')