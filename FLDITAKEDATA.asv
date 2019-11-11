%% FLDI TAKE DATA
%script to take data at one x/d location
    %everything in mm
%input x
x_loc= 20.6.*10; 
%input d
d=0.203.*25.4;
JetPos = x_loc./d;
%% Connect scope
[VisaInterface, Wavesurfer10] = OscopeVisaConnect(); %get interface and device
VSA = VisaInterface; %simplify code with shorthands
SCP = Wavesurfer10;
ACQ = Wavesurfer10.Acquisition; %simplify code with shorthands
TRG = Wavesurfer10.Trigger;
WVE = Wavesurfer10.Waveform;
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

CurrentStep = 0; % intitialize

FLDI_Data = struct('Jetposition',JetPos,'X_loc',x_loc,'d',d);
%take max data
    set(ACQ, 'Control', 'auto');
    set(ACQ, 'Control', 'single'); %set to single with matlab driver
    invoke(TRG, 'trigger'); 
    [Y,~, info] = invoke(WVE,'readwaveform','C3',false);
    FLDI_Data.Max = struct('Y',Y,'INFO',info);
    
%take min data
    set(ACQ, 'Control', 'auto');
    set(ACQ, 'Control', 'single'); %set to single with matlab driver
    invoke(TRG, 'trigger'); 
    [Y,~, info] = invoke(WVE,'readwaveform','C3',false);
    FLDI_Data.Min = struct('Y',Y,'INFO',info);

%take Spin Data
    set(ACQ, 'Control', 'auto');
    set(ACQ, 'Control', 'stop'); %set to single with matlab driver
    [Y,~, info] = invoke(WVE,'readwaveform','C3',false);
    FLDI_Data.Spin = struct('Y',Y,'INFO',info);
    
%take Zero Data
    set(ACQ, 'Control', 'auto');
    set(ACQ, 'Control', 'single'); %set to single with matlab driver
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
        d.Message=sprintf('Running loop %u at position %u of 101',i,CurrentStep);
        d.Value = (20*(CurrentStep-1)+i)./2020;
        set(ACQ, 'Control', 'single'); %set to single with matlab driver
        invoke(TRG, 'trigger'); 
        [Y,~, info] = invoke(WVE,'readwaveform','C3',false);
        FLDI_Data.YSweep(CurrentStep).Data(i).Y = Y;
        FLDI_Data.YSweep(CurrentStep).Data(i).Info = info;

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
SaveStr='20191110_FLDI_Data_XD30';
save(SaveStr,'FLDI_Data','-v7.3','-nocompression')