%% Find Matfiles
close all
cd('C:\Users\jack\Documents\FLDI FastWrite')
%cd('L:\FLDI_Sensitivity_study\Data');
Matfiles = dir('**/*.mat')';
%cd('C:\Users\jcobourn\Documents\GitHub\Focused_Laser_Dif_interf') %work
%cd('C:\Users\Jack\Documents\GitHub\Focused_Laser_Dif_interf') %laptop
%cd('D:\Jack\Documents\Focused_Laser_Dif_interf') %homepc
PCString = 'C:\Users\Jack\Documents\GitHub\Focused_Laser_Dif_interf';
cd(PCString);
%% Load Matfiles
%using matfiles() fuction, so they don't load totally into %memory
for ii = 1:length(Matfiles)
Matfiles(ii).Fullname = [Matfiles(ii).folder,'\',Matfiles(ii).name];
M(ii).Matfiles = matfile(Matfiles(ii).Fullname,'Writable', true);
end

%% Create Usable Cell Arrays

for ii = 1:length(M)
    Details{ii,1} = whos(M(ii).Matfiles);     
    Z(ii,1) = M(ii).Matfiles.z_pos;
    %M(ii).Ydist = Y(ii,1);
    VOLT(ii,:) = [M(ii).Matfiles.vmax(1,1),M(ii).Matfiles.vmin(1,1),...
        M(ii).Matfiles.vmax(1,2),M(ii).Matfiles.vmin(1,2)];
    FS(ii,1) = M(ii).Matfiles.Fs;
    BIT(ii,1) = M(ii).Matfiles.bitRes;
        
    CHA_TRI M{ii,1} = M(ii).Matfiles.chA_run(1:end,1); 
    %CHA{ii}(fix((FS(ii).*start)):fix((FS(ii).*stop)-1));
    CHB_TRIM{ii,1} = M(ii).Matfiles.chB_run(1:end,1);
    if BIT(ii) == 8
        CHA_TRIM{ii,1} = CHA_TRIM{ii,1}.*1000; %convert the teleydyne scope results to Volts
        CHB_TRIM{ii,1} = CHB_TRIM{ii,1}.*1000; %note offset doesn't matter due to mean subtracted pwelch
    end
     
end

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
            hann(pow2(18)),0.75*pow2(18), pow2(18),fix(FS(ii,1)));
        [PSDb{ii,1}, f{ii,1}] = pwelch(CHB_TRIM{ii,1}-mean(CHB_TRIM{ii,1}),...
            hann(pow2(18)),0.75*pow2(18), pow2(18),fix(FS(ii,1)));
        [MSC{ii,1}, f2{ii,1}] = mscohere(CHA_TRIM{ii,1},...
            CHB_TRIM{ii,1},hann(pow2(18)),0.75*pow2(18), pow2(18),fix(FS(ii,1)));
    elseif BIT(ii) == 8
        [PSDa{ii,1}, f{ii,1}] = pwelch(CHA_TRIM{ii,1}-mean(CHA_TRIM{ii,1}),...
            hann(pow2(18)),0.75*pow2(18), pow2(18),fix(FS(ii,1)));
        [PSDb{ii,1}, f{ii,1}] = pwelch(CHB_TRIM{ii,1}-mean(CHB_TRIM{ii,1}),...
            hann(pow2(18)),0.75*pow2(18), pow2(18),fix(FS(ii,1)));
        [MSC{ii,1}, f2{ii,1}] = mscohere(CHA_TRIM{ii,1},...
            CHB_TRIM{ii,1},hann(pow2(18)),0.75*pow2(18), pow2(18),fix(FS(ii,1)));

    end
end

%% Create IMageSc waterfall plot
[~,Index] = sortrows(struct2table(struct('Z',[Z])),{'Z'},{'descend'});

for ii = 1:length(M)
PSDedit{ii,1} = mean([PSDa{ii,1}(1:length(PSDa{1,1}))';PSDb{ii,1}(1:length(PSDb{1,1}))'],1);
end
%PSDedit(:) = PSDa(:)(1:length(PSDa{1,1}))'; is there a way to do above in
%1 line? could use cell2mat(PSDa) and then remap()
Colors = cell2mat(PSDedit);
%Colors = cell2mat(PSDa')';
fig3 = figure(3);
imagesc('XData',f{1,1},'CData',log10(Colors(Index,:)));
ax3 = fig3.Children;
grid on

%set colormap
colormap jet
 cb3 = colorbar();
 %caxis([-8,-4.5])
 tickvals = cb3.Ticks;
 for ii = 1:length(tickvals)
     tickname{ii} = sprintf('10^{%.1f}',tickvals(ii));
 end
 cb3.TickLabels = tickname;
cb3.Label.String = 'Power/Power';
fig3c = fig3.Children; Ax3 = fig3c(2);

%set alphamap

%set y
ylim([0.5,length(M)+0.5])
ylabel('Run Number')
Ax3.Layer = 'top';
Ax3.YTickMode = 'manual'; Ax3.YTickLabelMode = 'manual'; Ax3.YTick = 1:length(M);
Ax3.YAxis.TickDirection = 'out';
 for ii = 1:length(M)
     yticklab(ii) = {sprintf('Run %u, z=%.2f',M(ii).Matfiles.Runs+1,M(ii).Matfiles.z_pos)};
 end
Ax3.YTickLabel = yticklab(Index);

%set x
xlim([0,2e5])
xlabel('Frequency [Hz]')
Ax3.XAxis.TickDirection = 'out'; Ax3.XAxis.TickLength = [0.005,0.01];

%% Create Waterfall Plot
%loglog(f,PSDa,'linewidth',2)
fig1 = figure(1);

for ii = 1:length(M)
    h(ii) = waterfall(f{ii,1},M(ii).Matfiles.Runs+1,PSDa{ii,1}');
    hold on
end
ax1 = fig1.Children;
ax1.XScale = 'log';
%ax.YScale = 'log'
ax1.ZScale = 'log';
%ax1.View = [-37.5 30];
ax1.View = [0 0];


%% Plot Signals
fig5 = figure(5);
sgtitle({'FLDI signals from Sensitivity study') ;
for ii = 1:length(M)
    subplot(ceil(sqrt(length(M))),ceil(sqrt(length(M))),ii)
    plot((0:(length(CHA_TRIM{ii,1})-1))/FS(ii,1), [CHA_TRIM{ii,1}'; CHB_TRIM{ii,1}'])
    legend('Channel A','Channel B');
    grid on
    xlabel('t [s]')
    ylabel('Voltage [mV]')
    title(sprintf('Run Number %g, Z=%.2f', M(Index(ii)).Matfiles.Runs+1,Z(Index(ii))))
end

%% Plot Spectra
fig6 = figure(6);
sgtitle('FLDI spectra from Sensitivity study') ;
for ii = 1:length(M)
    subplot(ceil(sqrt(length(M))),ceil(sqrt(length(M))),ii)
    loglog(f{Index(ii),1}, PSDa{Index(ii),1})
    legend('Channel A FFT');
    grid on
    xlabel('f [Hz]')
    ylabel('Power')
    title(sprintf('Run Number %g, Z=%.2f', M(Index(ii)).Matfiles.Runs+1,Z(Index(ii))))
end

%% Plot corrilations


