%%FLDI PROCCESS Data
load('20191110_FLDI_Data_XD30.mat');

d = FLDI_Data.d;
jet = FLDI_Data.Jetposition;
scopey = []; %introduce variables
info = struct();
RMS = [];
STD = [];
AVE = [];

for i1=1:1:length(FLDI_Data.YSweep) %Sweep through positions
    for j=1:1:length(FLDI_Data.YSweep(1).Data) %sweep through samples
        scopey =  FLDI_Data.YSweep(i1).Data(j).Y; %loop internal scope data
        info =  FLDI_Data.YSweep(i1).Data(j).Info; %loop internal info
        scopey = cast(scopey,'double') .* info.ygain - info.yintercept; %convert to douuble
        FLDI_Results.YSweep(i1).Data(j).RMS = rms(scopey - mean(scopey)); %Calc RMS internal, write to struct
        RMS(j) = FLDI_Results.YSweep(i1).Data(j).RMS; %loop external vector
        FLDI_Results.YSweep(i1).Data(j).STD = std(scopey); %calc std internal and write to struct
        STD(j) = FLDI_Results.YSweep(i1).Data(j).STD; %loop external vector
        FLDI_Results.YSweep(i1).Data(j).AVE = mean(scopey); %calc ave internal and write to struct
        AVE(j) = FLDI_Results.YSweep(i1).Data(j).AVE; %loop external vector
        scopey = []; %clear loop internal scope data
        info = struct();%clear loop internal scope info
    end
    FLDI_Results.YSweep(i1).pos_y = FLDI_Data.YSweep(i1).TrueYPosition; %save position to struct
    
    FLDI_Results.YSweep(i1).RMSofRMS = rms(RMS - mean(RMS)); %calc outer loop rms and write to struct
    FLDI_Results.YSweep(i1).STDofRMS = std(RMS);
    FLDI_Results.YSweep(i1).AVEofRMS = mean(RMS);
    
    FLDI_Results.YSweep(i1).RMSofSTD = rms(STD -mean(STD));
    FLDI_Results.YSweep(i1).STDofSTD = std(STD);
    FLDI_Results.YSweep(i1).AVEofSTD = mean(STD);
    
    FLDI_Results.YSweep(i1).RMSofAVE = rms(AVE -mean(AVE));
    FLDI_Results.YSweep(i1).STDofAVE = std(AVE);
    FLDI_Results.YSweep(i1).AVEofAVE = mean(AVE);
end
for i2=1:1:length(FLDI_Data.YSweep)
    x_plot(i2) = FLDI_Results.YSweep(i2).pos_y; %save values to plot vectors
    y_plot(i2) = FLDI_Results.YSweep(i2).AVEofRMS;
    err_plot(i2) = FLDI_Results.YSweep(i2).STDofRMS;
end

save('20191111_FLDI_Results_XD10_vertical.mat','FLDI_Results','-v7.3','-nocompression');

f1 = figure(1);
P1 = errorbar(x_plot/d-25/d,y_plot,err_plot);
%P1 = plot(x_plot/d-25/d,y_plot);

P1.Color = 'black';
P1.LineWidth=1.5;

axObj=P1.Parent;
axObj.YGrid='on';
axObj.XGrid='on';
xlabel('Jet Location y/d in jet coordinates');
ylabel('Average RMS, (V), N = 20, vertical focii');
axObj.FontSize=26;

legend(sprintf('x/d=%3.3f',jet));
%% XD 40
load('20191110_FLDI_Data_XD40.mat');

scopey = []; %introduce variables
info = struct();
RMS = [];
STD = [];
AVE = [];

for i1=1:1:length(FLDI_Data.YSweep) %Sweep through positions
    for j=1:1:length(FLDI_Data.YSweep(1).Data) %sweep through samples
        scopey =  FLDI_Data.YSweep(i1).Data(j).Y; %loop internal scope data
        info =  FLDI_Data.YSweep(i1).Data(j).Info; %loop internal info
        scopey = cast(scopey,'double') .* info.ygain - info.yintercept; %convert to douuble
        FLDI_Results.YSweep(i1).Data(j).RMS = rms(scopey - mean(scopey)); %Calc RMS internal, write to struct
        RMS(j) = FLDI_Results.YSweep(i1).Data(j).RMS; %loop external vector
        FLDI_Results.YSweep(i1).Data(j).STD = std(scopey); %calc std internal and write to struct
        STD(j) = FLDI_Results.YSweep(i1).Data(j).STD; %loop external vector
        FLDI_Results.YSweep(i1).Data(j).AVE = mean(scopey); %calc ave internal and write to struct
        AVE(j) = FLDI_Results.YSweep(i1).Data(j).AVE; %loop external vector
        scopey = []; %clear loop internal scope data
        info = struct();%clear loop internal scope info
    end
    FLDI_Results.YSweep(i1).pos_y = FLDI_Data.YSweep(i1).TrueYPosition; %save position to struct
    
    FLDI_Results.YSweep(i1).RMSofRMS = rms(RMS - mean(RMS)); %calc outer loop rms and write to struct
    FLDI_Results.YSweep(i1).STDofRMS = std(RMS);
    FLDI_Results.YSweep(i1).AVEofRMS = mean(RMS);
    
    FLDI_Results.YSweep(i1).RMSofSTD = rms(STD -mean(STD));
    FLDI_Results.YSweep(i1).STDofSTD = std(STD);
    FLDI_Results.YSweep(i1).AVEofSTD = mean(STD);
    
    FLDI_Results.YSweep(i1).RMSofAVE = rms(AVE -mean(AVE));
    FLDI_Results.YSweep(i1).STDofAVE = std(AVE);
    FLDI_Results.YSweep(i1).AVEofAVE = mean(AVE);
end
for i2=1:1:length(FLDI_Data.YSweep)
    x_plot2(i2) = FLDI_Results.YSweep(i2).pos_y; %save values to plot vectors
    y_plot2(i2) = FLDI_Results.YSweep(i2).AVEofRMS;
    err_plot2(i2) = FLDI_Results.YSweep(i2).STDofRMS;
end

save('20191110_FLDI_Results_XD40.mat','FLDI_Results','-v7.3','-nocompression');

f2 = figure(2);
EB2 = errorbar(x_plot2,y_plot2,err_plot2);

%% combined plot

f3 = figure(3);
hold on
EB11 = errorbar(x_plot./d,y_plot,err_plot);
EB22 = errorbar(x_plot2./d,y_plot2,err_plot2);

EB11.Color = 'black';%UTKcolors('Smokey');
EB11.LineWidth=1.5;
%EB11.MarkerEdgeColor = UTKcolors('Rock');
EB11.MarkerFaceColor = UTKcolors('Rock');

EB22.Color = 'blue';%UTKcolors('Globe');
EB22.LineWidth=1.5;
EB22.MarkerEdgeColor = UTKcolors('River');
EB22.MarkerFaceColor = UTKcolors('River');

axObj=EB11.Parent;
axObj.YGrid='on';
axObj.XGrid='on';
xlabel('Jet Location y/d in jet coordinates');
ylabel('Average RMS, (V), N = 20');
axObj.FontSize=12;

legend('x/d=30','x/d=40');

% no eb
f4 = figure(4);
hold on
EB111 = plot(x_plot./d,y_plot);
EB222 = plot(x_plot2./d,y_plot2);

EB111.Color = 'black';
EB111.LineWidth=1.5;
EB222.Color = 'blue';
EB222.LineWidth=1.5;

axObj=EB111.Parent;
axObj.YGrid='on';
axObj.XGrid='on';
xlabel('Jet Location y/d in jet coordinates');
ylabel('Average RMS, (V), N = 20');
axObj.FontSize=26;

legend('x/d=30','x/d=40');