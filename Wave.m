function Waveform=Wave(OScopeObj,Channel)
Chanlist= {'channel1','channel2','channel3','channel4'};
Chanstr = Chanlist{Channel};

set(OScopeObj.Acquisition(1), 'Control', 'auto');
pause('on')
pause(0.25)
set(OScopeObj.Acquisition(1), 'Control', 'stop');

WaveObj = get(OScopeObj, 'Waveform');

[Y,X,YUNIT,XUNIT,HEADER] = invoke(WaveObj, 'readwaveform', Chanstr);
Waveform = struct('Y',Y,'X',X,'YUNIT',YUNIT,'XUNIT',XUNIT,'HEADER',HEADER);
end