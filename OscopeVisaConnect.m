function [VisaInterface, Wavesurfer510] = OscopeVisaConnect(~)
instrreset
clear VisaInterface Wavesurfer510 Waveform Wave
VisaInterface = visa('ni','TCPIP0::192.168.137.2::INSTR');
Wavesurfer510 = icdevice('D:\Jack\Documents\GitHub\Focused_Laser_Dif_interf\lecroy_WaveSurvfer510_VISA.mdd', VisaInterface);
%Wavesurfer510 = icdevice('lecroy_8600a.mdd', VisaInterface);

connect(Wavesurfer510)
%Waveform=Wavesurfer510.Waveform;

%VisaInterface
%Wavesurfer510
%Waveform
end
