%Process large Block of Pressure Data
ddrive = 'F:';
folderpath = 'FLDI_SHC_202006\Pressure';
cd([ddrive filesep folderpath])
PressureFiles = dir('**/*_*')';

for ii = 1:length(PressureFiles)
    PressureData(ii).Fullname = [PressureFiles(ii).folder,'\',PressureFiles(ii).name];
    PressureData(ii).RunNumber = ii;
    
    if PressureFiles(ii).bytes == 277
        PressureData(ii).Useful = 0;
    else
        PressureData(ii).Useful = 1;
    end
end

%handle some weird matlab inport rules
opts = detectImportOptions(PressureData(1).Fullname);
opts.VariableNames = {'DValvePct' 'VValvePct' 'Inlet' 'TotalP' 'StaticP' 'SupplyP' 'VacuumP' 'FatPipeP' 'DriverA' 'DriverB' 'DriverC' 'DriverD' 'CameraTrigger'};
        
fs = 12.5e3;

for ii = 6:length(PressureData)
    %if statement to select runs with useful data
    if PressureData(ii).Useful
        
        %check that the table looks right
        %Temp = preview(PressureData(ii).Fullname,opts)
        %fprintf('Press Enter if table looks ok')
        %pause
        %clc
        %load data, use the camera trigger to normalize t and then plot
        P = readtable(PressureData(ii).Fullname,opts);
        t(:,1) = 0:1/fs:(height(P)-1)/fs;
        Trigger = find(P.CameraTrigger > 0,1);
        t = t - t(Trigger);
        %plot
            fig = figure('Visible','off','Position',[-1000 100 100 100],'WindowStyle','modal','WindowState','maximized');
            plot(t, [P.FatPipeP P.DriverA])
            xlim([-0.10,0.35])
            legend(P.Properties.VariableNames([8 9]));
            grid on
            fig.Visible = 'on';
            fprintf('need to enter 3 points. ')
            fprintf('Preburst, start and stop of steady flow, only x vals matter')
            [x,~] = ginput(3);
            fig.WindowStyle = 'normal';
            close(fig)
            PressureData(ii).BurstP0= P.DriverA(find(t>=x(1),1)); 
            PressureData(ii).RunP0 = mean([P.FatPipeP(find(t>=x(2),1):find(t>=x(3),1))]);
            PressureData(ii).time_start = x(2);
            PressureData(ii).time_end = x(3);
            PressureData(ii).IsDataGood = input('Did that graph look okay? 1 for yes, 0 for no');
    else
        PressureData(ii).BurstP0= 60; 
        PressureData(ii).RunP0 = 0.9*PressureData(ii).BurstP0;
        PressureData(ii).time_start = 0.075;
        PressureData(ii).time_end = 0.15;
        PressureData(ii).IsDataGood = input('Did that graph look okay? 1 for yes, 0 for no');

    end
    clear t Temp Trigger x y P fig
clc
openvar('PressureData')
end
save('Pressures','PressureData','-v7.3')
