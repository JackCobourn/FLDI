function [VisaInterface, Wavesurfer10] = OscopeVisaConnect(~)
instrreset
clear VisaInterface Wavesurfer10 Waveform Wave
VisaInterface = visa('ni','TCPIP0::192.168.137.2::INSTR');
Wavesurfer10 = icdevice('D:\Jack\Documents\GitHub\Focused_Laser_Dif_interf\Oscilloscope\lecroy_WaveSurvfer10_VISA.mdd', VisaInterface);
%Wavesurfer10 = icdevice('lecroy_8600a.mdd', VisaInterface);

connect(Wavesurfer10)
%Waveform=Wavesurfer510.Waveform;

%VisaInterface
%Wavesurfer510
%Waveform

%beep that it worked
invoke(Wavesurfer10,'beep');
%go to defult setup
%invoke(Wavesurfer10,'reset')
%go to saved setting #6
%pause(2);
%invoke(Wavesurfer10.System,'loadstate',6)


%previously saved state
%invoke(Wavesurfer10.System,'savestate',6)

end
