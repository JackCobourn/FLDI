%SHCSchlierenProcess.m 
    %    Autor: Jack Cobourn          Creation Date: July 28, 2020
    
    
%load metadata
%% cd to Folder
cd('D:\Jack\Documents\Github\FLDI')
ddrive = 'F:';
folderpath = 'FLDI_SHC_202006';
folderpath_D = 'FLDI_SHC_202006\Data';
folderpath_P = 'FLDI_SHC_202006\Pressure';
folderpath_R = 'FLDI_SHC_202006\Results';
folderpath_S = 'FLDI_SHC_202006\Schlieren';
cd([ddrive filesep folderpath])

%% Load Metadata
SCH = table2struct(readtable('SHCMetadata.xlsb'));

%% Load Pressure Files
cd([ddrive filesep folderpath_P])
load('Pressures.mat');
cd([ddrive filesep folderpath_S])

%% Fix struct and create SCHClass/rez/vac in loop

for ii = 1:length(SCH)
    SCH(ii).Youtside_in_ = [SCH(ii).Youtside_in_,SCH(ii).Var7];
    SCH(ii).Xoutside_in_ = [SCH(ii).Xoutside_in_,SCH(ii).Var10];
end
for ii = 1:length(SCH)
    try
        SCH(ii).Sch = SchClass(SCH(ii).ReleventSch); %schlieren files make sure @SchClass folder is in a folder on path
    catch
        warning('SchClass not loaded, Assigning []');
        SCH(ii).Sch = [];
    end
end

SCH = rmfield(SCH,'Var7'); SCH = rmfield(SCH,'Var10'); 

for ii = 1:length(SCH)
    try
        SCH(ii).Rez = imread(SCH(ii).relScale);
    catch
        warning('Rez not loaded, Assigning []');
        SCH(ii).Rez = [];
    end
    if strcmp(SCH(ii).VacVidOrTiffOrN_A,'Tiff')
            SCH(ii).Vac = true;
        try
            SCH(ii).Vac1 = imread(SCH(ii).tiff1);
        catch
            warning('Vac1 not loaded, Assigning []');
            SCH(ii).Vac1 = [];
        end
        try
            SCH(ii).Vac2 = imread(SCH(ii).tiff2);
        catch
            warning('Vac2 not loaded, Assigning []');
            SCH(ii).Vac2 = [];
        end
    else
        SCH(ii).Vac = false;
        SCH(ii).Vac1 = [];
        SCH(ii).Vac2 = [];
    end
end
%% Measure the rezolution

%% calc movement from the vac

%% Load all data into struct just like to try
for ii= 1:length(SCH) 
    try
        SCH(ii).SCHRun = SCH(ii).Sch.getFrames(floor([PressureData(ii).time_start.*SCH(ii).Sch.FrameRate,PressureData(ii).time_end.*SCH(ii).Sch.FrameRate]));
    catch
        SCH(ii).SCHRun = [];
    end
end
%% Load actual schlieren in loop, don't load more than one at a time (i.e. no struct)
% Could change this to use above but w/e as it only needed to be done once
% test if Sch is usable
for ii= 1:length(SCH) %remember that SCH from run 1 is only ~150 images for some reason
    try
        SCHRun = SCH(ii).Sch.getFrames(floor([PressureData(ii).time_start.*SCH(ii).Sch.FrameRate,PressureData(ii).time_end.*SCH(ii).Sch.FrameRate]));
        MatMovie = implay(SCHRun(:,:,:).*(pow2(16)./pow2(12)),25); %scale the true 12 bits of the camera to 16 bits for implay, play at 25 fps
        good = input('please type 1 for good schlieren, 0 for bad, 0.5 for maybe');
        if good ==1
            SCH(ii).SchGood = 1;
        elseif good == 0.5
            SCH(ii).SchGood = 0.5;
        else 
            SCH(ii).SchGood = 0;
        end
        %clean up loop
        good = [];
        MatMovie.close
        clear MatMovie
    catch
        SCH(ii).SchGood = false;
        good = [];
    end    
end
%% Take ave and Rms
for ii= 1:length(SCH) 
    if SCH(ii).SchGood == 1
        try
            SCH(ii).SCHave = mean(single(SCH(ii).SCHRun),3);
            SCH(ii).SCHrms = rms(single(SCH(ii).SCHRun)-mean(single(SCH(ii).SCHRun),3),3);
        catch
            SCH(ii).SCHave = [];
            SCH(ii).SCHrms = [];
        end
    end
end

%% Stabalize with edgefind (Doesn't work as well)
for ii= 2%1:length(SCH) 
    if SCH(ii).SchGood == 1
        tic
        SCH(ii).SCHedge = ones(size(SCH(ii).SCHRun),'logical');
        for jj= 1:size(SCH(ii).SCHRun,3)
            SCH(ii).SCHedge(:,:,jj) = edge(SCH(ii).SCHRun(:,:,jj),'Prewitt',0.0035);
        end
        toc
        
        %% Try stabilization with phase corrilation
        tic
        for jj= 2:size(SCH(ii).SCHRun,3)
           SCH(ii).tform(jj) = imregcorr(SCH(ii).SCHedge(:,:,jj),SCH(ii).SCHedge(:,:,1),'rigid');
            
        end
        toc
        
        tic
        SCH(ii).SCHemaped(:,:,1) = SCH(ii).SCHedge(:,:,1);
        for jj= 2:size(SCH(ii).SCHRun,3)
            sameAsInput = affineOutputView(size(SCH(ii).SCHRun(:,:,jj)),SCH(ii).tform(jj),'BoundsStyle','SameAsInput');
           [SCH(ii).SCHmaped(:,:,jj),~] = imwarp(SCH(ii).SCHRun(:,:,jj),SCH(ii).tform(jj),'linear','OutputView',sameAsInput);            
        end
        toc
        %% Try stab with intesity corrilation
        
        %% try stabilization with Hugh
         MatMovie = implay(SCH(ii).SCHedge(:,:,:).*(pow2(16)./pow2(12)),25);
         MatMovie = implay(SCH(ii).SCHemaped(:,:,:));
         MatMovie = implay(SCH(ii).SCHmaped(:,:,:).*(pow2(16)./pow2(12)),25);
            SCH(ii).SCHmaprms = rms(single(SCH(ii).SCHmaped)-mean(single(SCH(ii).SCHmaped),3),3);
imshow(SCH(ii).SCHmaprms,[])

    end
end

%% Get Background images
for ii= 1:length(SCH) 
    if SCH(ii).SchGood == 1
        %calc BG as ave of five good wind off frames
        SCH(ii).BG = mean(single(SCH(ii).Sch.getFrames([2:6])),3);
        
        %align to first image of flow on
        transform = imregcorr(SCH(ii).BG,single(SCH(ii).SCHRun(:,:,1)),'rigid');
        sameAsInput = affineOutputView(size(single(SCH(ii).SCHRun(:,:,1))),transform,'BoundsStyle','SameAsInput');
        SCH(ii).BGMapped = imwarp(SCH(ii).BG,transform,'linear','OutputView',sameAsInput);
    end
end
clear transform sameAsInput
%% Stabalize with full image
for ii= 1:length(SCH) 
    if SCH(ii).SchGood == 1
        %% Try stabilization with phase corrilation
        %use full images as background is helpful
        fprintf('Run number %u\n',ii)
        tic
        for jj= 2:size(SCH(ii).SCHRun,3)
           SCH(ii).tform(jj) = imregcorr(SCH(ii).SCHRun(:,:,jj),SCH(ii).SCHRun(:,:,1),'rigid');
            
        end
        toc
        %% Align Images
        tic
        SCH(ii).SCHmaped(:,:,1) = SCH(ii).SCHRun(:,:,1);
        for jj= 2:size(SCH(ii).SCHRun,3)
            sameAsInput = affineOutputView(size(SCH(ii).SCHRun(:,:,jj)),SCH(ii).tform(jj),'BoundsStyle','SameAsInput');
           [SCH(ii).SCHmaped(:,:,jj),~] = imwarp(SCH(ii).SCHRun(:,:,jj),SCH(ii).tform(jj),'linear','OutputView',sameAsInput);            
        end
        toc
     %MatMovie = implay(SCH(ii).SCHRun(:,:,:).*(pow2(16)./pow2(12)),25);
     %MatMovie = implay(SCH(ii).SCHmaped(:,:,:).*(pow2(16)./pow2(12)),25);
    
    end
end
clear sameAsInput
%% Take ave and Rms of mapped
for ii= 1:length(SCH) 
    if SCH(ii).SchGood == 1
            SCH(ii).SCHmappedave = mean(single(SCH(ii).SCHmaped),3);
            SCH(ii).SCHmappedrms = rms(single(SCH(ii).SCHmaped)-mean(single(SCH(ii).SCHmaped),3),3);
    else
        SCH(ii).SCHmappedave = [];
        SCH(ii).SCHmappedrms = [];
    end
end
%% BG sub
thresh = 30;
gain = .2;
for ii= 1:length(SCH) 
    if SCH(ii).SchGood == 1
        %% Remove BG from all images
        for jj= 1:size(SCH(ii).SCHRun,3)
            % Threshold Difference from BG and set all blow thresh to
            % black
            FrameDiff = abs(single(SCH(ii).SCHmaped(:,:,jj)) - single(SCH(ii).BGMapped));
            Foreground_logical = (FrameDiff > thresh); % Tresholding, can be adjusted
            Foreground = (single(SCH(ii).SCHmaped(:,:,jj))).*Foreground_logical;
            
            %rescale image 
            SCH(ii).BGrem(:,:,jj) = Foreground; %(Foreground-min(min(Foreground))).*gain;
        end
     %MatMovie = implay(SCH(ii).SCHRun(:,:,:).*(pow2(16)./pow2(12)),25);
     %MatMovie = .*(pow2(16)./pow2(12)),25);
    
    end
end


%% plot one by one and grab bottom of model and top of bl
[SCH.ydist_Nan] = deal(0);
for ii= 1:length(SCH) 
    if SCH(ii).SchGood == 1
       disp('Select a bunch of points along the model')
       maxplot = imshowpair(SCH(ii).SCHmappedave,SCH(ii).SCHmappedrms,'montage');
       maxplot.Parent.Parent.WindowState = 'maximized';
       [~,y_model]=ginput;
       disp('Select a bunch of points along the BL')
       [~,y_BL]=ginput;
        close(maxplot.Parent.Parent)   
        if numel(y_BL)>0
           SCH(ii).BLH = abs(mean(y_BL)-mean(y_model));
       else
           SCH(ii).ydist_Nan = 1;
        end
    end
end

%% find scale
LPM = 1.12;
Lw = 1/(2*LPM);
Ll = 5*Lw;
for ii= 1:length(SCH) 
    scaleplot = imshow(SCH(ii).Rez,[]);
    scaleplot.Parent.Parent.WindowState = 'maximized';
    scaleplot.Interpolation ='bilinear';
    scaleplot.Parent.InnerPosition = [0 0 1 1 ];
    [Xcorner,Ycorner]=ginput(4);
    Xcorner = sort(Xcorner); Ycorner=sort(Ycorner);
    Xpix = abs(mean(Xcorner(1:2))-mean(Xcorner(3:4)));
    Ypix = abs(mean(Ycorner(1:2))-mean(Ycorner(3:4)));
    SCH(ii).Ypixmm = Ypix./Ll;
    SCH(ii).Xpixmm = Xpix./Ll;                    
        if SCH(ii).SchGood == 1 && SCH(ii).ydist_Nan ~= 1
            SCH(ii).BLHmm = SCH(ii).BLH./SCH(ii).Ypixmm;
            SCH(ii).Festimate = 677*1000/SCH(ii).BLHmm;
        end
    close(scaleplot.Parent.Parent)
end

%% Tiled Print
load([ddrive filesep folderpath_R filesep 'SHCResults.mat'], 'IndexNo');
IndexNo(IndexNo>=5) = IndexNo(IndexNo>=5)+1; %Sch doesn't cut stuff that M does

figtile = figure();
Tiles = tiledlayout(figtile,3,4,'TileSpacing','Compact');
% ylabel(Tiles,'Frequency [Hz]','FontSize',16)
% xlabel(Tiles,'Time [s]','FontSize',16)

for ii = 1:length(IndexNo)
Tb(ii) = nexttile;
imagesc(SCH(IndexNo(ii)).SCHmappedave)
colormap gray

Tb(ii).Title.String = sprintf('Phase Registered Ave, Local Re = %0.3G',SCH(IndexNo(ii)).Re);
Tb(ii).Title.FontSize=16;

grid on
Tb(ii).XMinorGrid='on';Tb(ii).YMinorGrid='on';Tb(ii).GridAlpha=.35;Tb(ii).MinorGridAlpha=.25;

end

figtile = figure();
Tiles = tiledlayout(figtile,3,4,'TileSpacing','Compact');
% ylabel(Tiles,'Frequency [Hz]','FontSize',16)
% xlabel(Tiles,'Time [s]','FontSize',16)

for ii = 1:length(IndexNo)
Tb(ii) = nexttile;
%imagesc(SCH(IndexNo(ii)).SCHmappedrms)
if ~isempty(SCH(IndexNo(ii)).SCHmappedave)
imshowpair(SCH(IndexNo(ii)).SCHmappedave,SCH(IndexNo(ii)).SCHmappedrms,'montage')
end
colormap gray

Tb(ii).Title.String = {'Phase Registered Ave & RMS,';sprintf('Local Re = %0.3G',SCH(IndexNo(ii)).Re)};
Tb(ii).Title.FontSize=16;

grid on
Tb(ii).XMinorGrid='on';Tb(ii).YMinorGrid='on';Tb(ii).GridAlpha=.35;Tb(ii).MinorGridAlpha=.25;

end


%% Save all besides images
vars = whos;
save('SHC_sch_Results',vars(~contains({vars.class},"matlab")).name,'-v7.3'); %Prevents saving of fig handles as they have 'matlab.' in their classname
clear vars
%% extra code
% %% Stabalize schlieren
% for ii= 1:length(SCH) %remember that SCH from run 1 is only ~150 images for some reason
%             SCHRun = SCH(ii).Sch.getFrames(floor([PressureData(ii).time_start.*SCH(ii).Sch.FrameRate,PressureData(ii).time_end.*SCH(ii).Sch.FrameRate]));
%     if SCH(ii).SchGood ==1
%         [centers,radii] = imfindcircles(SCHRun(:,:,1),[100,400],'ObjectPolarity','bright');
%         MatIm = imshow(SCHRun(:,:,1),[]);
%         viscircles(centers, radii,'EdgeColor','b');
%     pause
%     end
% end

%% Get lines
function [x,y] = GetLine(A)
disp('select 3 points along bottom edge, the 3 points along BL')
guinput

end

%% Get Crop