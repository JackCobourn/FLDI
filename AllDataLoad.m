%put all mat files in a folder you could even hard link it
cd('F:\FLDI_UMD\Data')
Matfiles = dir('**/*.mat');
for i = 1:length(matfiles)
    