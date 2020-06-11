%% Schlieren Height Change Process
% Jack Cobourn cw 20200602

cd('J:\FLDI_HCF_Raw_202004\Schlieren')
rez_files = dir('**/*rez*.tif');
vac_files = dir('**/*vac.mraw');

%% Rez Files
for ii = 1:length(rez_files)
     rez_files_ISO(ii) = regexp(rez_files(ii).name,'(?<Year>\d\d\d\d)(?<Month>\d\d)(?<Day>\d\d)','Names')';
     
     for fn = fieldnames(rez_files_ISO)'
        rez_files(ii).(fn{1}) = rez_files_ISO(ii).(fn{1});
     end
     rez_files(ii).Datetime = datetime(str2num(rez_files(ii).Year),str2num(rez_files(ii).Month),str2num(rez_files(ii).Day));

     rez_files(ii).Tiff = Tiff([rez_files(ii).folder filesep rez_files(ii).name]);
     rez_files(ii).Image = imread([rez_files(ii).folder filesep rez_files(ii).name]);
     fig1 = figure(1);
     Im1 = imshow(rez_files(ii).Image,[]);
     fig1.WindowState = 'maximized';
     disp('Select top left and bottom right of 0-1 square, 2.5 mm')
     [x,y] = ginput(2);
     Im2 = imshow(rez_files(ii).Image(floor(min(y))-25:ceil(max(y))+25,floor(min(x))-25:ceil(max(x))+25),[]);
     fig1.WindowState = 'maximized';
     [x,y] = ginput(2);
     rez_files(ii).HozPixPerMM = (max(x)-min(x))/2.5;
     rez_files(ii).VerpixPerMM = (max(y)-min(y))/2.5;

end
for ii = 1:length(rez_files)
     
end
clear rez_files_ISO x y fig1 Im1 Im2 ii fn
%% Vac Files
for ii = 1:length(vac_files)
     vac_files_ISO(ii) = regexp(vac_files(ii).name,'(?<Year>\d\d\d\d)(?<Month>\d\d)(?<Day>\d\d)','Names')';
     for fn = fieldnames(vac_files_ISO)'
         vac_files(ii).(fn{1}) = vac_files_ISO(ii).(fn{1});
     end
     vac_files(ii).Datetime = datetime(str2num(vac_files(ii).Year),str2num(vac_files(ii).Month),str2num(vac_files(ii).Day));
        
     jj = 1;
     while vac_files(ii).Datetime > rez_files(jj).Datetime || vac_files(ii).Datetime == rez_files(jj).Datetime
         vac_files(ii).rezfilenum = jj;
         jj = jj+1;
         if jj > length(rez_files)
             break
         end
     end    
     I.FirstImage=readmraw(strrep([vac_files(ii).folder filesep vac_files(ii).name],'.mraw',''),1);
     I.LastImage=readmraw(strrep([vac_files(ii).folder filesep vac_files(ii).name],'.mraw',''),I.FirstImage.CameraSetup.Total_Frames);
     
     fig1 = figure(1);
     Im1 = imshow(I.FirstImage.Images.RawImages,[]);
     fig1.WindowState = 'maximized';
     disp('Select edge on far left and far right of iluminated zone')
     [x,y] = ginput(2);
     Im2 = imshow(I.FirstImage.Images.RawImages(floor(min(y))-25:ceil(max(y))+25,floor(min(x)),[]));
     fig1.WindowState = 'maximized';
     [x,y] = ginput(2);
     
     fig2 = figure(2);
     Im3 = imshow(I.LastImage.Images.RawImages,[]);
     fig1.WindowState = 'maximized';
     disp('Select edge on far left and far right of iluminated zone')
     [x,y] = ginput(2);
     Im4 = imshow(I.LastImage.Images.RawImages(floor(min(y))-25:ceil(max(y))+25,floor(min(x)),[]);
     fig1.WindowState = 'maximized';
     [x,y] = ginput(2);
end
%
clear  vac_files_ISO
save('vac_height_change','-v7.3')