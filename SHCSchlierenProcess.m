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

%% Fix struct and create SCHClass in loop

for ii = 1:length(SCH)
    SCH(ii).Youtside_in_ = [SCH(ii).Youtside_in_,SCH(ii).Var7];
    SCH(ii).Xoutside_in_ = [SCH(ii).Xoutside_in_,SCH(ii).Var10];
    
    try
        SCH(ii).Sch = SchClass(SCH(ii).ReleventSch); %schlieren files
    catch
        warning('SchClass not loaded, Assigning []');
        SCH(ii).Sch = [];
    end
end

SCH = rmfield(SCH,'Var7'); SCH = rmfield(SCH,'Var10'); 

for ii = 1:length(SCH)
