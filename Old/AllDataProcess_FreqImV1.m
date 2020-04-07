%% Find Matfiles
close all
cd('D:\FLDI_UMD\Data')
Matfiles = dir('**/*.mat')';
cd('C:\Users\jcobourn\Documents\GitHub\Focused_Laser_Dif_interf')
%cd('D:\Jack\Documents\GitHub\Focused_Laser_Dif_interf')

%% Load Matfiles
%using matfiles() fuction, so they don't load totally into %memory
for ii = 1:length(Matfiles)
Matfiles(ii).Fullname = [Matfiles(ii).folder,'\',Matfiles(ii).name];
M(ii).Matfiles = matfile(Matfiles(ii).Fullname);
end
%add CampainRunNum to initial runs Works out to this pattern asedon date
%and naming
EarlyRunNumbers = [1 2 3 4 5 7 6 8] ;
% for ii = 1:8
%     M(ii).Matfiles.CampainRunNum = EarlyRunNumbers(ii);
% end
emptyIndex = [];
%M.Matfiles = M(~find(M(1:31).Matfiles.num_diaphrams ~= 3).Matfiles
%exclude runs not That have erronious data
for ii = 1:length(Matfiles)
if ii>8  %eliminate by particular campain numbers, below 8 don't have this number;
    if M(ii).Matfiles.CampainRunNum == 27 || M(ii).Matfiles.CampainRunNum == 21 ||...
            M(ii).Matfiles.CampainRunNum == 12
    %M(ii).Matfiles = [];
    emptyIndex(ii) = 1;
    else
        emptyIndex(ii) = 0;
    end
else
    emptyIndex(ii) = 0;
end
end
M = M(~emptyIndex);
clear emptyIndex
for ii = 1:length(M)
    if any(ii==1:8) %Add run number to early runs
        M(ii).RunNumber = EarlyRunNumbers(ii);
    else
        M(ii).RunNumber = M(ii).Matfiles.CampainRunNum;
    end
end
%Put mat files in run order in the script
[~,ordercode] = sort([M.RunNumber]);%reordercode
M = M(ordercode);
clear ordercode EarlyRunNumbers

%% add outside data
start = readmatrix('D:\FLDI_UMD\DFLDI_UMD_UTSI_Campaign_M4_BL.xlsx','Sheet','Data','Range','K3:K32'); %[s] trim window start
stop = readmatrix('D:\FLDI_UMD\DFLDI_UMD_UTSI_Campaign_M4_BL.xlsx','Sheet','Data','Range','L3:L32'); %[s] trim window stop
RunP0 = readmatrix('D:\FLDI_UMD\DFLDI_UMD_UTSI_Campaign_M4_BL.xlsx','Sheet','Data','Range','J3:J32');
for ii = 1:length(M)
    M(ii).start = start(ii);
    M(ii).stop = stop(ii);
    M(ii).RunP0 = RunP0(ii);
end
dx2_opts = num2cell(readmatrix('D:\FLDI_UMD\DFLDI_UMD_UTSI_Campaign_M4_BL.xlsx',...
    'Sheet','Beam Pics','Range','K2:K10'));
BeamPicRel = readmatrix('D:\FLDI_UMD\DFLDI_UMD_UTSI_Campaign_M4_BL.xlsx','Sheet','Data','Range','G3:G32');
[M.dx2] = deal(dx2_opts{BeamPicRel});

clear start stop RunP0 dx2_opts BeamPicRel

%% Load Pressure Files
cd('D:\FLDI_UMD\PressureData')
PressureFiles = dir('**/*_UMD')';
cd('C:\Users\jcobourn\Documents\GitHub\Focused_Laser_Dif_interf')

for ii = 1:length(PressureFiles)
    PressureData(ii).Fullname = [PressureFiles(ii).folder,'\',PressureFiles(ii).name];
    PressureData(ii).RunNumber = ii;
    if any(PressureData(ii).RunNumber == [12 21 27])
        BadRun(ii) = 1;
    else
        BadRun(ii) = 0;
    end
 
    if PressureFiles(ii).bytes == 277
        PressureData(ii).Useful = 0;
    else
        PressureData(ii).Useful = 1;
    end
end %note, need to check that this still works with run in March
PressureData = PressureData(~BadRun);
clear PressureFiles
[M.PressureFile] = deal(PressureData.Fullname); [M.GoodP] = deal(PressureData.Useful);
%handel some wierd matlab inport rules
opts = detectImportOptions(M(1).PressureFile);
opts.VariableNames = {'DValvePct' 'VValvePct' 'Inlet' 'TotalP' 'StaticP'...
    'SupplyP' 'VacuumP' 'FatPipeP' 'DriverA' 'DriverB' 'DriverC' 'DriverD' 'CameraTrigger'};     
fs_pressure = 12.5e3;

%% Load Pressure Curves
for ii = 1:length(M)
    %if statement to select runs with useful data
    if PressureData(ii).Useful
        %load data, use the camera trigger to normalize t and then plot
        M(ii).P = readtable(M(ii).PressureFile,opts);
        t(:,1) = 0:1/fs_pressure:(height(M(ii).P)-1)/fs_pressure;
        Trigger = find(M(ii).P.CameraTrigger > 0,1);
        t = t - t(Trigger);
        M(ii).P.t = t;
        %plot
%             fig = figure();
%             plot(t, [P.FatPipeP P.DriverA])
%             xlim([-0.10,0.35])
%             legend(P.Properties.VariableNames([8 9]));
%             grid on
%            
    end
    clear t Temp Trigger x y P fig
end

%% move runs accosciated with different Re to a different matfile struct
%gets rid of runs at different Re
% for ii = 1:length(M)
% if M(ii).Matfiles.num_diaphrams ~= 3
%     WrongReIndex(ii) = true;
% else
%     WrongReIndex(ii) = false;
% end
% end
% Re = M(WrongReIndex);
% M = M(~WrongReIndex);
% clear WrongReIndex

%% Create Usable Cell Arrays


for ii = 1:length(M)
    Details{ii,1} = whos(M(ii).Matfiles);     
    Y(ii,1) = M(ii).Matfiles.y_dist;
    M(ii).Ydist = Y(ii,1);
    VOLT(ii,:) = [M(ii).Matfiles.vmax(1,1),M(ii).Matfiles.vmin(1,1),...
        M(ii).Matfiles.vmax(1,2),M(ii).Matfiles.vmin(1,2)];
    FS(ii,1) = M(ii).Matfiles.Fs;
    BIT(ii,1) = M(ii).Matfiles.bitRes;
    
    CHA_TRIM{ii,1} = M(ii).Matfiles.chA_run(fix((FS(ii).*M(ii).start)):fix((FS(ii).*M(ii).stop)-1),1);        %CHA{ii}(fix((FS(ii).*start)):fix((FS(ii).*stop)-1));
    CHB_TRIM{ii,1} = M(ii).Matfiles.chB_run(fix((FS(ii).*M(ii).start)):fix((FS(ii).*M(ii).stop)-1),1);
    if BIT(ii) == 8
        CHA_TRIM{ii,1} = CHA_TRIM{ii,1}.*1000; %convert the teleydyne scope results to Volts
        CHB_TRIM{ii,1} = CHB_TRIM{ii,1}.*1000; %note offset doesn't matter due to mean subtracted pwelch
    end
     
end

%% Calculate Re
Mach = 4;
gamma = 1.4;
R = 287.058;
[~, Tratio, Pratio, ~, ~] = flowisentropic(gamma, Mach, 'mach');
P = [M.RunP0].*Pratio.*6894.76; %convert psi to pa
T = 300*Tratio;
rho = P./(R.*T);
a = sqrt(gamma*R*T);
U0 = Mach*a;
mu = 0.00001458*(T)^1.5/(T+110.4); %sotherlans law
Rex =  num2cell(rho.*U0./mu);
[M.Re] = deal(Rex{:});
clear Mach gamma R Tratio Pratio P T rho a mu

%% Spectral Analysis
numPoints = arrayfun(@(x) length(x{:}),CHA_TRIM);
Closestpow2 = pow2(floor(log2(numPoints./100)));
deltaf = FS./Closestpow2;
Fmax = FS/2;
lengthfreqVec = ceil(Fmax./deltaf)+1;
%WelchMatrix = zeroes(cell2mat(CHA_TRIM')';
%PSDa
for ii = 1:length(M)
    ii   
    if BIT(ii) == 15
        [PSDa{ii,1}, f{ii,1}] = pwelch(CHA_TRIM{ii,1}-mean(CHA_TRIM{ii,1}),...
            pow2(13),[], pow2(13)*12.5,fix(FS(ii,1)));
        [MSC{ii,1}, f2{ii,1}] = mscohere(CHA_TRIM{ii,1}-mean(CHA_TRIM{ii,1}),...
            CHB_TRIM{ii,1}-mean(CHB_TRIM{ii,1}),pow2(13),[], pow2(13),fix(FS(ii,1)));
    elseif BIT(ii) == 8
        [PSDa{ii,1}, f{ii,1}] = pwelch(CHA_TRIM{ii,1}-mean(CHA_TRIM{ii,1}),...
            pow2(14),[], pow2(13)*20,fix(FS(ii,1)));
        [MSC{ii,1}, f2{ii,1}] = mscohere(CHA_TRIM{ii,1}-mean(CHA_TRIM{ii,1}),...
            CHB_TRIM{ii,1}-mean(CHB_TRIM{ii,1}),pow2(14),[], pow2(14),fix(FS(ii,1)));

    end
end
for ii = 1:length(M)
    ii
    figure(10)
    semilogx(f2{ii,1},MSC{ii,1})
    drawnow
    pause(2)
end

% cellfun(@(X) mean(diff(X)),f)
% cellfun(@(X) X(end),f)
% cellfun(@(X) length(X),f)
% f15 = f{1,1}; f8 = f{30,1};

for ii = 1:length(M)
    %calc PSD area
    PSDArea(ii,1) = trapz(PSDa{ii,1}).*mean(diff(f{ii,1}));
    PSDNorm{ii,1} = PSDa{ii,1}./PSDArea(ii,1); %normalize to area of 1
        
%         if BIT(ii) == 15 %account for fs difference %Dont do this, just
%         use forier interpolation and norm by area
%             PSDNorm{ii,1} = PSDNorm{ii,1}*(1/12.5);
%             %Inter{ii,1} = griddedInterpolant(f15,PSDNorm{ii,1});
%             %PSDNorm{ii,1} = [PSDNorm{ii,1}'  zeros(size(PSDNorm{30,1},1)-size(PSDNorm{ii,1},1))]'
%         elseif BIT(ii) == 8
%             PSDNorm{ii,1} = PSDNorm{ii,1}*(1/20);
%         end
         NormResults(ii,1) = trapz(PSDNorm{ii,1})*mean(diff(f{ii,1}));
end
for ii = 1:length(M)
if BIT(ii) == 8
    PSDNorm{ii,1} = PSDNorm{ii,1}(1:length(f{1,1}));
    f{ii,1} = f{ii,1}(1:length(f{1,1}));
end
end
%% Create IMageSc waterfall plot
[~,Index] = sortrows(struct2table(struct('Re',[M.Re]','Y',fix(10*[M.Ydist]'))),{'Y','Re'},{'Ascend','descend'});

% for ii = 1:length(M)
% PSDedit{ii,1} = PSDNorm{ii,1}(1:length(PSDa{1,1}))';
% end
%PSDedit(:) = PSDa(:)(1:length(PSDa{1,1}))'; is there a way to do above in
%1 line? could use cell2mat(PSDa) and then remap()
Colors = 1e5*cell2mat(PSDNorm')';
%Colors = cell2mat(PSDa')';
fig3 = figure(3);
imagesc('XData',f{1,1},'CData',log10(Colors(Index,:)))
ax3 = fig3.Children;
grid on

%set colormap
colormap jet
 cb3 = colorbar();
 caxis([-6,-1])
 tickvals = cb3.Ticks;
 for ii = 1:length(tickvals)
     tickname{ii} = sprintf('10^{%.1f}',tickvals(ii));
 end
 cb3.TickLabels = tickname;
cb3.Label.String = 'Power*10^5/Power';
fig3c = fig3.Children; Ax3 = fig3c(2);

%set alphamap

%set y
ylim([0.5,length(M)+0.5])
ylabel('Y [in]')
Ax3.Layer = 'top';
Ax3.YTickMode = 'manual'; Ax3.YTickLabelMode = 'manual'; Ax3.YTick = 1:length(M);
 for ii = 1:length(M)
     yticklab(ii) = {sprintf('y=%.2f, %-2dbit, Re=%8.3G%',Y(Index(ii)),BIT(Index(ii)),M(Index(ii)).Re)};
 end
Ax3.YTickLabel = yticklab;

%set x
xlim([0,2e6])
xlabel('Frequency [Hz]')
Ax3.XAxis.TickDirection = 'out'; Ax3.XAxis.TickLength = [0.005,0.01];
%Did all this with the XDATA command
% tixspace = 5e4;
% Ax3.XTickMode = 'manual'; Ax3.XAxis.TickValues = find(rem(f8,tixspace)==0);
% Ax3.XAxis.Exponent=6; Ax3.XTickLabelMode = 'manual'; 
% Ax3.XAxis.TickLabels = f8(rem(f{30,1},tixspace)==0)/1e6;
% 
% Ax3.XAxis.TickLabelFormat = '%0.3g';

%% Create Waterfall Plot
% %loglog(f,PSDa,'linewidth',2)
% fig1 = figure(1);
% 
% for ii = 1:length(M)
%     h(ii) = waterfall(f{ii,1},Y(ii,1),PSDNorm{ii,1}');
%     hold on
% end
% ax1 = fig1.Children;
% ax1.XScale = 'log';
% %ax.YScale = 'log'
% ax1.ZScale = 'log';
% %ax1.View = [-37.5 30];
% ax1.View = [0 0];

%% Create waterfall scatter plot
% fig3 = figure(3);
% for ii = 1:length(Matfiles)
%     h3(ii) = scatter(f{ii,1},Y(ii,1).*ones(size(f{ii,1})),15,PSDa{ii,1}','filled');
%     hold on
% end
% ax3 = fig3.Children;
% ax3.XScale = 'log';
% xlim([10e3 inf])
% colormap colorcube
% cb3 = colorbar();
% cb3.Ruler.Scale = 'log';
% cb3.Ruler.MinorTick = 'on';

%% Create BL prof
for ii = 1:length(M)
[cross_cor{ii,1},lag{ii,1}] = xcorr(CHA_TRIM{ii,1}-mean(CHA_TRIM{ii,1}),-(CHB_TRIM{ii,1}-mean(CHB_TRIM{ii,1})), 'coeff');
%find max of xcorr and index
[cc(ii,1),I(ii,1)] = max(cross_cor{ii,1});
%check if xcorr peak is high enough?
if cc{ii,1}>0.3
    lagDiff(ii,1) = lag{ii,1}(I{ii,1});
    %3 pt gaussian interpolation, to get a more exact velocity
    peakOffset(ii,1) = ((log(cross_cor{ii,1}(I{ii,1}-1)))-(log(cross_cor{ii,1}(I{ii,1}+1))))/(2*((log(cross_cor{ii,1}(I{ii,1}-1)))+(log(cross_cor{ii,1}(I{ii,1}+1)))-2*(log(cc{ii,1}))));
    %find delta t
    DT(ii,1) = abs((lagDiff(ii,1)+peakOffset(ii,1))/FS(ii,1));
else
    lagDiff = NaN;
    DT = NaN;
end
%disturbance velocity, [m/s]
Uc(ii,1) = dx2(ii)/1000/DT(ii,1);

% figure(4)
% clf
% plot(lag,cross_cor)
end
fig2 = figure(2);
h2 = scatter(Uc(:,1),25.4*Y(:,1),25,[M.Re]','filled');
cmap = colormap(jet);
%cmap = [1  1  1 ; cmap]; %hide low re
colormap(cmap)
cb2 = colorbar();
%caxis([1.4e6,inf])
caxis([-inf,inf])
cb2.Label.String = 'Re/x (m^{-1})';
xlabel('Velocity (m/s)')
ylabel('Height above floor (mm)')

%% Print Voltages and Delta Ts
writematrix(VOLT,'D:\FLDI_UMD\Voltages.xlsx');
writematrix(DT,'D:\FLDI_UMD\Offsets.xlsx');

%% Plot Pressures
fig4 = figure(4);
sgtitle([{'Pressure signals for runs with significant velocity spread'};{'y = 2.50'}]) ;
for ii = 1:length(M)
    subplot(ceil(sqrt(length(M))),ceil(sqrt(length(M))),ii)
    plot(M(ii).P.t, [M(ii).P.FatPipeP M(ii).P.DriverA])
    xlim([-0.10,0.35])
    legend(M(ii).P.Properties.VariableNames([8 9]));
    grid on
    xlabel('t [s]')
    ylabel('Pressure [Psia]')
    title(sprintf('Run Number %g, Ux = %.2f', M(ii).RunNumber,Uc(ii)))
end

%% Plot Signals
fig5 = figure(5);
sgtitle([{'FLDI signals for runs with significant velocity spread'};{'y = 2.50'}]) ;
for ii = 1:length(M)
    subplot(ceil(sqrt(length(M))),ceil(sqrt(length(M))),ii)
    plot((0:(length(CHA_TRIM{ii,1})-1))/FS(ii,1), [CHA_TRIM{ii,1}'; CHB_TRIM{ii,1}'])
    legend('Channel A','Channel B');
    grid on
    xlabel('t [s]')
    ylabel('Voltage [mV]')
    title(sprintf('Run Number %g, Ux = %.2f', M(ii).RunNumber,Uc(ii)))
end

%% Plot Spectra
fig6 = figure(6);
sgtitle([{'FLDI spectra for runs with significant velocity spread'};{'y = 2.50'}]) ;
for ii = 1:length(M)
    subplot(ceil(sqrt(length(M))),ceil(sqrt(length(M))),ii)
    loglog(f{ii,1}, PSDNorm{ii,1})
    legend('Channel A FFT');
    grid on
    xlabel('f [Hz]')
    ylabel('Power')
    title(sprintf('Run Number %g, Ux = %.2f', M(ii).RunNumber,Uc(ii)))
end

%% Plot corrilations


