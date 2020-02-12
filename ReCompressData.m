%put all mat files in a folder you could even hard link it
cd('D:\FLDI_UMD\Data')
Matfiles = dir('**/*.mat');
for i = 1:length(Matfiles)
    i
    load([Matfiles(i).folder,'\', Matfiles(i).name])
    C = who('-file',[Matfiles(i).folder,'\', Matfiles(i).name]);
    save([Matfiles(i).folder,'\', Matfiles(i).name],C{:})
end
