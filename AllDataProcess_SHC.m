%% Find Matfiles
cd('D:\Jack\Documents\Focused_Laser_Dif_interf')
ddrive = 'F:';
folderpath_D = 'FLDI_SHC_202006\Data';
folderpath_P = 'FLDI_SHC_202006\Pressure';
folderpath_R = 'FLDI_SHC_202006\Results';
cd([ddrive filesep folderpath_D])
Matfiles = dir('**/*.mat')';

%% Load Matfiles
%using matfiles() fuction, so they don't load totally into %memory
for ii = 1:length(Matfiles)
Matfiles(ii).Fullname = [Matfiles(ii).folder,'\',Matfiles(ii).name];
M(ii).Matfiles = matfile(Matfiles(ii).Fullname,'Writable', true);
end

%% Load Pressure Files
cd([ddrive filesep folderpath_P])
load('Pressures.mat');
cd([ddrive filesep folderpath_D])


%% exclude runs that have erronious data
for ii = 1:length(Matfiles)
    if M(ii).Matfiles.CampainRunNum == 5
    %M(ii).Matfiles = [];
    emptyIndex(ii) = 1;
    else
        emptyIndex(ii) = 0;
    end
end
M = M(~emptyIndex);
%PressureData = PressureData(logical([PressureData.IsDataGood]));
PressureData = PressureData([PressureData.RunNumber]~=5);
clear emptyIndex

%% Create Usable Cell Arrays
F2K = @(F) (F-32)./1.8+273.15;
for ii = 1:length(M)
    Details{ii,1} = whos(M(ii).Matfiles);     
%     Yout(ii,:) = M(ii).Matfiles.youtside;
%     Xout(ii,:) = M(ii).Matfiles.xoutside;
%     Yin(ii,1) = M(ii).Matfiles.yinside;
%     Xin(ii,1) = M(ii).Matfiles.xinside;
    VOLT(ii,:) = [M(ii).Matfiles.vmax(1,1),M(ii).Matfiles.vmin(1,1),...
        M(ii).Matfiles.vmax(1,2),M(ii).Matfiles.vmin(1,2)];
    FS(ii,1) = M(ii).Matfiles.Fs;
    BIT(ii,1) = M(ii).Matfiles.bitRes;
        
    CHA_TRIM{ii,1} = M(ii).Matfiles.chA_run(fix((FS(ii).*PressureData(ii).time_start)):fix((FS(ii).*PressureData(ii).time_end)-1),1); 
    CHB_TRIM{ii,1} = M(ii).Matfiles.chB_run(fix((FS(ii).*PressureData(ii).time_start)):fix((FS(ii).*PressureData(ii).time_end)-1),1);
    CHA_NOISE1_TRIM{ii,1} = M(ii).Matfiles.chA_noise(fix((FS(ii).*PressureData(ii).time_start)):fix((FS(ii).*PressureData(ii).time_end)-1),1); 
    CHB_NOISE1_TRIM{ii,1} = M(ii).Matfiles.chB_noise(fix((FS(ii).*PressureData(ii).time_start)):fix((FS(ii).*PressureData(ii).time_end)-1),1);
    CHA_NOISE2_TRIM{ii,1} = M(ii).Matfiles.chA_noise2(fix((FS(ii).*PressureData(ii).time_start)):fix((FS(ii).*PressureData(ii).time_end)-1),1); 
    CHB_NOISE2_TRIM{ii,1} = M(ii).Matfiles.chB_noise2(fix((FS(ii).*PressureData(ii).time_start)):fix((FS(ii).*PressureData(ii).time_end)-1),1);
    if BIT(ii) == 8
        CHA_TRIM{ii,1} = CHA_TRIM{ii,1}.*1000; %convert the teleydyne scope results to Volts
        CHB_TRIM{ii,1} = CHB_TRIM{ii,1}.*1000; %note offset doesn't matter due to mean subtracted pwelch
        CHA_NOISE1_TRIM{ii,1} = CHA_NOISE1_TRIM{ii,1}.*1000;
        CHB_NOISE1_TRIM{ii,1} = CHB_NOISE1_TRIM{ii,1}.*1000;
        CHA_NOISE2_TRIM{ii,1} =  CHA_NOISE2_TRIM{ii,1}.*1000;
        CHB_NOISE2_TRIM{ii,1} = CHB_NOISE2_TRIM{ii,1}.*1000;
    end
    
    if any(strcmpi({Details{ii,1}.name;},'temp_tunnel'))
        T0(ii,1) = F2K(M(ii).Matfiles.temp_tunnel);
    else
        T0(ii,1) = NaN;
    end
end
T0(isnan(T0))=nanmean(T0);

%% Calculate Re
Mach = 4;
gamma = 1.4;
R = 287.058;
[~, Tratio, Pratio, ~, ~] = flowisentropic(gamma, Mach, 'mach');
P = [PressureData.RunP0].*Pratio.*6894.76; %convert psi to pa
T = T0.*Tratio;
rho = P'./(R.*T);
a = sqrt(gamma*R.*T);
U0 = Mach.*a;
mu = ((1.458.*10^-6).*((T).^1.5))./(T+110.4); %sutherlans law
Rex =  num2cell(rho.*U0./mu);
[M.Rex] = deal(Rex{:});
clear Mach gamma R Tratio Pratio P T rho a mu
%% Read Positions From Metadata List (best option with difficulties measureing X every time.
DSLocList = table2array(readtable([ddrive filesep 'FLDI_SHC_202006' filesep 'SHCMetadata.xlsb'],'range','L2:L17'));
DSLocList = num2cell(DSLocList(~isnan(DSLocList)));
[M.DSLocList] = deal(DSLocList{:});
Re = num2cell(arrayfun(@(A,B) A.*B./1000,[M.Rex],[M.DSLocList])');
[M.Re] = deal(Re{:}); 
[~,Index] = sortrows(Re);
IndexNo = Index(Index~=9 & Index~=15 & Index~=1); %remove FS and Out Of BL and Initial Out of BL
%% RMS analysis
for ii = 1:length(M)
    RMSA{ii,1} = std(CHA_TRIM{ii,1}-mean(CHA_TRIM{ii,1}));
    RMSB{ii,1} = std(CHB_TRIM{ii,1}-mean(CHB_TRIM{ii,1}));
    RMSNA{ii,1} = std(CHA_NOISE1_TRIM{ii,1}-mean(CHA_NOISE1_TRIM{ii,1}));
    RMSNB{ii,1} = std(CHB_NOISE1_TRIM{ii,1}-mean(CHB_NOISE1_TRIM{ii,1}));
    RMSN2A{ii,1} = std(CHA_NOISE2_TRIM{ii,1}-mean(CHA_NOISE2_TRIM{ii,1}));
    RMSN2B{ii,1} = std(CHB_NOISE2_TRIM{ii,1}-mean(CHB_NOISE2_TRIM{ii,1}));
end
rmsfig = figure(25);
L1 = plot([Re{IndexNo}],[RMSA{IndexNo}],'r',[Re{IndexNo}], [RMSB{IndexNo}],'b','LineWidth',1.4);
hold on
grid on
L2 = plot([Re{IndexNo}],[RMSNA{IndexNo}],'k',[Re{IndexNo}], [RMSNB{IndexNo}],'k--','LineWidth',1.4);
L3 = plot([Re{IndexNo}],[RMSN2A{IndexNo}],'k-.',[Re{IndexNo}], [RMSN2B{IndexNo}],'k:','LineWidth',1.4);
xlabel('Channel A local Reynolds Number [Re]','FontSize',16);
ylabel('Mean Subtracted RMS [mV]','FontSize',16);
rmsfigax = gca;
rmsfigax.XMinorGrid =  'on'; rmsfigax.XMinorTick='on'; rmsfigax.MinorGridAlpha=.25; rmsfigax.GridAlpha=.35;
rmsfigax.XColor = 'k';
rmsfigax.YColor = 'k';
rmsfigax2 = axes('Position',rmsfigax.Position,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none','YTick',[],'TickLength',[0 0]);
hold on
L4 = plot([DSLocList{IndexNo}],[RMSA{IndexNo}],'k.','LineWidth',.7);
xlim([30,200])
L5 = plot(DSLocList{9},RMSA{9},'r*',DSLocList{9},RMSB{9},'b*','LineWidth',1.4);
L6 = plot(DSLocList{15},RMSA{15},'ro',DSLocList{15},RMSB{15},'bo','LineWidth',1.4);
xlabel('Streamwise Location of Channel A [mm]','FontSize',16);
Leg = legend([L1;L2;L3]','Channel A rms','Channel B rms','Channel A pre-run noise rms',...
    'Channel B pre-run noise rms','Channel A post-run noise rms','Channel B post-run noise rms','Location','northwest'); Leg.FontSize=16; Leg.Title.String='By Re'; Leg.Box='off';
rmsfigax2.Position =rmsfigax.Position;
Leg2 = legend([L4;L5(1);L6(1)]','A Channel Streamwise Locations','Freestream A/B','Outside BL A/B','Location','northeast'); Leg2.FontSize=16; Leg2.Title.String='By Position'; Leg2.Box='off';
ylim(rmsfigax.YLim);
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
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
        [PSDb{ii,1}, ~] = pwelch(CHB_TRIM{ii,1}-mean(CHB_TRIM{ii,1}),...
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
        [MSC{ii,1}, f2{ii,1}] = mscohere(CHA_TRIM{ii,1},...
            CHB_TRIM{ii,1},hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
    elseif BIT(ii) == 8
        [PSDa{ii,1}, f{ii,1}] = pwelch(CHA_TRIM{ii,1}-mean(CHA_TRIM{ii,1}),...
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
        [PSDb{ii,1}, ~] = pwelch(CHB_TRIM{ii,1}-mean(CHB_TRIM{ii,1}),...
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
        [MSC{ii,1}, f2{ii,1}] = mscohere(CHA_TRIM{ii,1},...
            CHB_TRIM{ii,1},hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
        
        [PSDNa{ii,1}, ~] = pwelch(CHA_NOISE1_TRIM{ii,1}-mean(CHA_NOISE1_TRIM{ii,1}),...
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
        [PSDNb{ii,1}, ~] = pwelch(CHB_NOISE1_TRIM{ii,1}-mean(CHB_NOISE1_TRIM{ii,1}),...
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
        [MSCN{ii,1}, f2{ii,1}] = mscohere(CHA_NOISE1_TRIM{ii,1},...
            CHB_NOISE1_TRIM{ii,1},hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
        
        [PSDN2a{ii,1}, ~] = pwelch(CHA_NOISE2_TRIM{ii,1}-mean(CHA_NOISE2_TRIM{ii,1}),...
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
        [PSDN2b{ii,1}, ~] = pwelch(CHB_NOISE2_TRIM{ii,1}-mean(CHB_NOISE2_TRIM{ii,1}),...
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
        [MSCN2{ii,1}, f2{ii,1}] = mscohere(CHA_NOISE2_TRIM{ii,1},...
            CHB_NOISE2_TRIM{ii,1},hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));

    end
end


%% Create IMageSc waterfall plot
%[~,Index] = sortrows(struct2table(struct('Re',[M.Re]','Y',fix(10*[M.Ydist]'))),{'Y','Re'},{'Ascend','descend'});

% for ii = 1:length(M)
% PSDedit{ii,1} = PSDNorm{ii,1}(1:length(PSDa{1,1}))';
% end
%PSDedit(:) = PSDa(:)(1:length(PSDa{1,1}))'; is there a way to do above in
%1 line? could use cell2mat(PSDa) and then remap()
Colors = cell2mat(PSDa(IndexNo)')';
%Colors = cell2mat(PSDa')';
fig3 = figure(3);
imagesc('XData',f{1,1}.*10^-6,'CData',log10(Colors))
ax3 = fig3.Children;
grid on

%set colormap
colormap jet
 cb3 = colorbar();
 caxis([-8,-3.5])
 tickvals = cb3.Ticks;
 for ii = 1:length(tickvals)
     tickname{ii} = sprintf('10^{%.1f}',tickvals(ii));
 end
 cb3.TickLabels = tickname;
cb3.Label.String = 'PSD [mV^2/Hz]';
fig3c = fig3.Children; Ax3 = fig3c(2);
Ax3.FontSize = 15;
cb3.Label.FontSize = 20;
Ax3.YLabel.FontSize = 24;
Ax3.XLabel.FontSize = 24;


%set alphamap

%set y
ylim([0.5,length(IndexNo)+0.5])
ylabel('Channel A Local R_e')
Ax3.Layer = 'top';
Ax3.YTickMode = 'manual'; Ax3.YTickLabelMode = 'manual'; Ax3.YTick = 1:length(M);
Ax3.YAxis.TickDirection = 'out';
 for ii = 1:length(IndexNo)
     yticklab(ii) = {sprintf('%0.3G',M(IndexNo(ii)).Re)};
 end
Ax3.YTickLabel = yticklab;

%set x
xlim([0,2])
xlabel('Frequency [MHz]')
Ax3.XAxis.TickDirection = 'out'; Ax3.XAxis.TickLength = [0.005,0.01];
%Did all this with the XDATA command
% tixspace = 5e4;
% Ax3.XTickMode = 'manual'; Ax3.XAxis.TickValues = find(rem(f8,tixspace)==0);
% Ax3.XAxis.Exponent=6; Ax3.XTickLabelMode = 'manual'; 
% Ax3.XAxis.TickLabels = f8(rem(f{30,1},tixspace)==0)/1e6;
% 
% Ax3.XAxis.TickLabelFormat = '%0.3g';

%% Create Waterfall Plot
%loglog(f,PSDa,'linewidth',2)
fig1 = figure(1);

for ii = 1:length(IndexNo)
    h(ii) = waterfall(f{IndexNo(ii),1},M(IndexNo(ii)).Re,PSDa{IndexNo(ii),1}');
    hold on
end
ax1 = fig1.Children;
ax1.XScale = 'log';
%ax.YScale = 'log'
ax1.ZScale = 'log';
%ax1.View = [37.5 30];
ax1.View = [0 0];

%% Tiled Layout
figtile = figure();
Tiles = tiledlayout(figtile,3,4,'TileSpacing','Compact');
ylabel(Tiles,'PSD [mV^2]','FontSize',28)
xlabel(Tiles,'Frequency [Hz]','FontSize',28)
IndexNo2 = IndexNo;
IndexNo2(IndexNo>=5) = IndexNo(IndexNo>=5)+1; %Sch doesn't cut stuff that M does

for ii = 1:length(IndexNo)
Ta(ii) = nexttile;
loglog(f{IndexNo(ii),1},PSDa{IndexNo(ii),1}','r',f{IndexNo(ii),1},PSDb{IndexNo(ii),1},'b');
hold on
loglog(f{IndexNo(ii),1},PSDNa{IndexNo(ii),1}','k',f{IndexNo(ii),1},PSDNb{IndexNo(ii),1},'k--');
loglog(f{IndexNo(ii),1},PSDN2a{IndexNo(ii),1}','k-.',f{IndexNo(ii),1},PSDN2b{IndexNo(ii),1},'k:');
Ta(ii).Title.String = sprintf('Local Re = %0.3G',M(IndexNo(ii)).Re);
Ta(ii).FontSize=18;
Ta(ii).Title.FontSize=20;
Ta(ii).XAxis.Label.FontSize=28;
Ta(ii).YAxis.Label.FontSize=28;

% if ii==4
%     l1 = legend('Channel A','Channel B','Ch A pre-run noise',...
%     'Ch B pre-run noise','Ch A post-run noise','Ch B post-run noise');
% 	l1.Title.String=sprintf('Power Spectral\nDensity [mV^2]');
%     l1.Title.FontSize=16;
%     l1.FontSize=14;
%     l1.Location='northeastoutside';
% end
if any(IndexNo2(ii) == [2,3,9,11,12])
    xline(SCH(IndexNo2(ii)).Festimate,'Color','green','LineWidth',1.8)
end

grid on
end
l1 = legend(Ta(1),{'Channel A','Channel B','Ch A pre-run noise',...
    'Ch B pre-run noise','Ch A post-run noise','Ch B post-run noise','U_{edge}/y_{BL}'});
l1.Title.String=sprintf('Power Spectral\nDensity [mV^2]');
    l1.Title.FontSize=16;
    l1.FontSize=14;
    l1.Location='northeastoutside';
linkaxes([Ta(:)])
% xlim([3e3,5e6])
% ylim([10^(-8.5),10^(-3.5)])
%print([ddrive filesep folderpath_r filesep 'Spectra_V1')


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
    %savefig(fig3,[savestring '.fig']);
    %saveas(fig3,savestring,'jpeg');
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

%% Plot Spectragram
for ii = 1:length(M)
    [~, fspc{ii,1},tspc{ii,1},SPCa{ii,1}] = spectrogram(CHA_TRIM{ii,1}-mean(CHA_TRIM{ii,1}),...
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
    [~, fspc{ii,1},tspc{ii,1},SPCb{ii,1}] = spectrogram(CHB_TRIM{ii,1}-mean(CHB_TRIM{ii,1}),...
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
    [~, fspc{ii,1},tspc{ii,1},SPCNa{ii,1}] = spectrogram(CHA_NOISE1_TRIM{ii,1}-mean(CHA_NOISE1_TRIM{ii,1}),...
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
    [~, fspc{ii,1},tspc{ii,1},SPCNb{ii,1}] = spectrogram(CHB_NOISE1_TRIM{ii,1}-mean(CHB_NOISE1_TRIM{ii,1}),...
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
    [~, fspc{ii,1},tspc{ii,1},SPCN2a{ii,1}] = spectrogram(CHA_NOISE2_TRIM{ii,1}-mean(CHA_NOISE2_TRIM{ii,1}),...
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
    [~, fspc{ii,1},tspc{ii,1},SPCN2b{ii,1}] = spectrogram(CHB_NOISE2_TRIM{ii,1}-mean(CHB_NOISE2_TRIM{ii,1}),...
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
end

for ii = 1:length(M)
    ii
    [~, fspc_Full{ii,1},tspc_Full{ii,1},SPCa_Full{ii,1}] = spectrogram(M(ii).Matfiles.chA_run-mean(M(ii).Matfiles.chA_run),...
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
end

figtile = figure();
Tiles = tiledlayout(figtile,4,3,'TileSpacing','Compact');
ylabel(Tiles,'Frequency [MHz]','FontSize',16)
xlabel(Tiles,'Time [s]','FontSize',16)
for ii = 1:length(IndexNo)
Tb(ii) = nexttile;
imagesc('XData',tspc{IndexNo(ii)},'YData',fspc{IndexNo(ii)},'CData',mag2db(SPCa{IndexNo(ii)}))
colormap jet
Tb(ii).Title.String = sprintf('Local Re = %0.3G',M(IndexNo(ii)).Re);
Tb(ii).Title.FontSize=16;
ylim([-inf,1e6])
     cb(ii) = colorbar();
     cb(ii).Label.String = 'PSD [dB/Hz]';
     caxis([-200,0])
%     l1 = legend('Channel A','Channel B','Ch A pre-run noise',...
%     'Ch B pre-run noise','Ch A post-run noise','Ch B post-run noise');
% 	l1.Title.String=sprintf('Power Spectral\nDensity [mV^2]');
%     l1.Title.FontSize=16;
%     l1.FontSize=14;
%     l1.Location='northeastoutside';
 
grid on
Tb(ii).XMinorGrid='on';Tb(ii).YMinorGrid='on';Tb(ii).GridAlpha=.35;Tb(ii).MinorGridAlpha=.25;

end
linkaxes([Tb(:)])
print([ddrive filesep folderpath_r filesep 'SpcGram_V1'],'-dsvg');

figtile = figure();
Tiles = tiledlayout(figtile,4,3,'TileSpacing','Compact','Padding','normal');
ylabel(Tiles,'Frequency [MHz]','FontSize',28)
xlabel(Tiles,'Time [s]','FontSize',28)
for ii = 1:length(IndexNo)
Tb(ii) = nexttile;
imagesc('XData',tspc_Full{IndexNo(ii)},'YData',fspc_Full{IndexNo(ii)}./10^6,'CData',log10(SPCa_Full{IndexNo(ii)}))
hold on
colormap jet
xline(PressureData(IndexNo(ii)).time_start,'Color','black','LineWidth',2.5)
xline(PressureData(IndexNo(ii)).time_end,'Color','black','LineWidth',2.5)
Tb(ii).Title.String = sprintf('Local Re = %0.3G',M(IndexNo(ii)).Re);
Tb(ii).FontSize=18;
Tb(ii).Title.FontSize=20;
ylim([-inf,2.5])
xlim([0.059526209677419,0.165372983870968])
caxis([-15,-10])
if any(ii == [3,6,9,12])
     cb(ii) = colorbar();
     tickvals = cb(ii).Ticks;
  for jj = 1:length(tickvals)
      tickname{jj} = sprintf('10^{%.1f}',tickvals(jj));
  end
  cb(ii).TickLabels = tickname;
 cb(ii).Label.String = 'PSD [mV^2/Hz]';
 cb(ii).Label.FontSize = 18;
end
grid on
Tb(ii).XMinorGrid='on';Tb(ii).YMinorGrid='on';Tb(ii).GridAlpha=.35;Tb(ii).MinorGridAlpha=.25;
for ll=0.5:0.5:2.5
    yline(ll)
end
for ll=0.06:0.02:0.16
    xline(ll)
end
end
save('SpcGram_V3.mat','figtile','-v7.3')
linkaxes([Tb(:)])
print([ddrive filesep folderpath_r filesep 'SpcGram_V2'],'-dsvg')

%% Save everything but images.
vars = whos;
save('SHCResults',vars(~contains({vars.class},"matlab")).name,'-v7.3'); %Prevents saving of fig handles as they have 'matlab.' in their classname
clear vars
