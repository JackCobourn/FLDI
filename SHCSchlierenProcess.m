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

%% Show images and grab ave BL
for ii= 1:length(SCH) 
    if SCH(ii).SchGood == 1
        imshow(SCH(ii).SCHave,[])
        pause
        imshow(SCH(ii).SCHrms,[])
        pause
    end
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
imagesc(SCH(IndexNo(ii)).SCHave)
colormap gray

Tb(ii).Title.String = sprintf('Ave, Local Re = %0.3G',SCH(IndexNo(ii)).Re);
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
imagesc(SCH(IndexNo(ii)).SCHrms)
colormap gray

Tb(ii).Title.String = sprintf('Rms, Local Re = %0.3G',SCH(IndexNo(ii)).Re);
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