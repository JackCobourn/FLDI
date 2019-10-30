%%FLDI PROCCESS Data
FLDI_Data.YSweep(CurrentStep).Data(i).RMS = rms(FLDI_Data.YSweep(CurrentStep).Data(i).Wave.Y-mean(FLDI_Data.YSweep(CurrentStep).Data(i).Wave.Y));
FLDI_Data.YSweep(CurrentStep).Data(i).STD = std(FLDI_Data.YSweep(CurrentStep).Data(i).Wave.Y);
FLDI_Data.YSweep(CurrentStep).Data(i).mean = mean(FLDI_Data.YSweep(CurrentStep).Data(i).Wave.Y);
    