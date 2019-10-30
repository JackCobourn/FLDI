function [VisaInterface, Wavesurfer510, Waveform] = OscopeVisaConnect_V2(~)
instrreset
clear VisaInterface Wavesurfer510 Waveform Wave
VisaInterface = visa('ni','TCPIP0::192.168.137.2::INSTR');
Wavesurfer510 = icdevice('lecroy_WaveSurvfer510_VISA.mdd', VisaInterface);
%Wavesurfer510 = icdevice('lecroy_8600a.mdd', VisaInterface);

connect(Wavesurfer510)
Waveform=Wavesurfer510.Waveform;

VisaInterface
Wavesurfer510
Waveform
end

%Set stuff by hand on Oscope, its easier for now
%set(Waveform,'Precision','int8')
%[Y1,X1,YUNIT1,XUNIT1,HEADER1] = invoke(Waveform,'readwaveform','channel3');
%set(Waveform,'Precision','int16')
%[Y2,X2,YUNIT2,XUNIT2,HEADER2] = invoke(Waveform,'readwaveform','channel3');

%get(Wavesurfer510.Acquisition,'T
%template2 = query(VisaInterface, 'TEMPLATE?');

