%% Find Matfiles
close all
cd('I:\FLDI_HCF_Workable\Data')
Matfiles = dir('**/*.mat')';
%cd('C:\Users\jcobourn\Documents\GitHub\Focused_Laser_Dif_interf') %work
%cd('D:\Jack\Documents\GitHub\Focused_Laser_Dif_interf') %laptop
%cd('D:\Jack\Documents\Focused_Laser_Dif_interf') %homepc
PCString = 'D:\Jack\Documents\Focused_Laser_Dif_interf';
cd(PCString);
%% Load Matfiles
%using matfiles() fuction, so they don't load totally into %memory
for ii = 1:length(Matfiles)
Matfiles(ii).Fullname = [Matfiles(ii).folder,'\',Matfiles(ii).name];
M(ii).Matfiles = matfile(Matfiles(ii).Fullname,'Writable', true);
end

%% Load Pressure Files
cd('I:\FLDI_HCF_Workable\Pressure')
load('Pressures.mat');
cd(PCString)


%% exclude runs not That have erronious data
for ii = 1:length(Matfiles)
    if M(ii).Matfiles.CampainRunNum == 4 || M(ii).Matfiles.CampainRunNum == 5
    %M(ii).Matfiles = [];
    emptyIndex(ii) = 1;
    else
        emptyIndex(ii) = 0;
    end
end
M = M(~emptyIndex);
PressureData = PressureData(~emptyIndex);
clear emptyIndex

%% Create Usable Cell Arrays

for ii = 1:length(M)
    Details{ii,1} = whos(M(ii).Matfiles);     
    %Y(ii,1) = M(ii).Matfiles.y_dist;
    %M(ii).Ydist = Y(ii,1);
    VOLT(ii,:) = [M(ii).Matfiles.vmax(1,1),M(ii).Matfiles.vmin(1,1),...
        M(ii).Matfiles.vmax(1,2),M(ii).Matfiles.vmin(1,2)];
    FS(ii,1) = M(ii).Matfiles.Fs;
    BIT(ii,1) = M(ii).Matfiles.bitRes;
        
    CHA_TRIM{ii,1} = M(ii).Matfiles.chA_run(fix((FS(ii).*PressureData(ii).time_start)):fix((FS(ii).*PressureData(ii).time_end)-1),1); 
    %CHA{ii}(fix((FS(ii).*start)):fix((FS(ii).*stop)-1));
    CHB_TRIM{ii,1} = M(ii).Matfiles.chB_run(fix((FS(ii).*PressureData(ii).time_start)):fix((FS(ii).*PressureData(ii).time_end)-1),1);
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
P = [PressureData.RunP0].*Pratio.*6894.76; %convert psi to pa
T = 300*Tratio;
rho = P./(R.*T);
a = sqrt(gamma*R*T);
U0 = Mach*a;
mu = 0.00001458*(T)^1.5/(T+110.4); %sutherlans law
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
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
        [PSDb{ii,1}, f{ii,1}] = pwelch(CHB_TRIM{ii,1}-mean(CHB_TRIM{ii,1}),...
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
        [MSC{ii,1}, f2{ii,1}] = mscohere(CHA_TRIM{ii,1},...
            CHB_TRIM{ii,1},hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
    elseif BIT(ii) == 8
        [PSDa{ii,1}, f{ii,1}] = pwelch(CHA_TRIM{ii,1}-mean(CHA_TRIM{ii,1}),...
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
        [PSDb{ii,1}, f{ii,1}] = pwelch(CHB_TRIM{ii,1}-mean(CHB_TRIM{ii,1}),...
            hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));
        [MSC{ii,1}, f2{ii,1}] = mscohere(CHA_TRIM{ii,1},...
            CHB_TRIM{ii,1},hann(pow2(13)),0.75*pow2(13), pow2(13),fix(FS(ii,1)));

    end
end
% fig10 = figure(10)
% fig10.Position = [-1000,250,500,500]; fig10.WindowStyle = 'normal'; fig10.WindowState = 'Maximized';
% for ii = 1:length(M)
%     ii
%     clf(fig10)
%      
%     ax1 = axes('Units','normalized','Position',[0.05 0.05 0.45 .90]);
%     loglog(f{ii,1},PSDa{ii,1},f{ii,1},PSDb{ii,1})
%     x=logspace(4,6,10000); y=2e5*x.^(-5/3); hold on; loglog(x,y,'k--'); hold off;
%     %ax1.YAxis.Scale='log'; ax1.XAxis.Scale='log';
%     
%     ax2 = axes('Units','normalized','Position',[0.525 0.05 0.45 .90]);
%     semilogx(f2{ii,1},MSC{ii,1})
%     yline(.7,'k--');
%     freqCutOff = f2{ii,1}(find((f2{ii,1}<1e6 & MSC{ii,1}>0.7),1,'last'));
%     xline(freqCutOff,'k');
%     
%     axes(ax1);
%     xline(freqCutOff,'k');
%     
%     fig10.WindowStyle = 'normal'; fig10.WindowState = 'Maximized';
%     drawnow
%     pause(3)
%    
% end
% 
% % cellfun(@(X) mean(diff(X)),f)
% % cellfun(@(X) X(end),f)
% % cellfun(@(X) length(X),f)
% % f15 = f{1,1}; f8 = f{30,1};
% 
% for ii = 1:length(M)
%     %calc PSD area
%     PSDArea(ii,1) = trapz(PSDa{ii,1}).*mean(diff(f{ii,1}));
%     PSDNorm{ii,1} = PSDa{ii,1}./PSDArea(ii,1); %normalize to area of 1
%         
% %         if BIT(ii) == 15 %account for fs difference %Dont do this, just
% %         use forier interpolation and norm by area
% %             PSDNorm{ii,1} = PSDNorm{ii,1}*(1/12.5);
% %             %Inter{ii,1} = griddedInterpolant(f15,PSDNorm{ii,1});
% %             %PSDNorm{ii,1} = [PSDNorm{ii,1}'  zeros(size(PSDNorm{30,1},1)-size(PSDNorm{ii,1},1))]'
% %         elseif BIT(ii) == 8
% %             PSDNorm{ii,1} = PSDNorm{ii,1}*(1/20);
% %         end
%          NormResults(ii,1) = trapz(PSDNorm{ii,1})*mean(diff(f{ii,1}));
% end
% for ii = 1:length(M)
% if BIT(ii) == 8
%     PSDNorm{ii,1} = PSDNorm{ii,1}(1:length(f{1,1}));
%     f{ii,1} = f{ii,1}(1:length(f{1,1}));
% end
% end
%% Create IMageSc waterfall plot
%[~,Index] = sortrows(struct2table(struct('Re',[M.Re]','Y',fix(10*[M.Ydist]'))),{'Y','Re'},{'Ascend','descend'});

% for ii = 1:length(M)
% PSDedit{ii,1} = PSDNorm{ii,1}(1:length(PSDa{1,1}))';
% end
%PSDedit(:) = PSDa(:)(1:length(PSDa{1,1}))'; is there a way to do above in
%1 line? could use cell2mat(PSDa) and then remap()
Colors = 1e5*cell2mat(PSDa')';
%Colors = cell2mat(PSDa')';
fig3 = figure(3);
imagesc('XData',f{1,1},'CData',log10(Colors))
ax3 = fig3.Children;
grid on

%set colormap
colormap jet
 cb3 = colorbar();
 caxis([-3.5,-1])
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
ylabel('Run Number')
Ax3.Layer = 'top';
Ax3.YTickMode = 'manual'; Ax3.YTickLabelMode = 'manual'; Ax3.YTick = 1:length(M);
Ax3.YAxis.TickDirection = 'out';
 for ii = 1:length(M)
     yticklab(ii) = {sprintf('%u',M(ii).Matfiles.CampainRunNum)};
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
%loglog(f,PSDa,'linewidth',2)
fig1 = figure(1);

for ii = 1:length(M)
    h(ii) = waterfall(f{ii,1},M(ii).Matfiles.CampainRunNum,PSDa{ii,1}');
    hold on
end
ax1 = fig1.Children;
ax1.XScale = 'log';
%ax.YScale = 'log'
ax1.ZScale = 'log';
%ax1.View = [-37.5 30];
ax1.View = [0 0];

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
if cc(ii,1)>0.3
    lagDiff(ii,1) = lag{ii,1}(I(ii,1));
    %3 pt gaussian interpolation, to get a more exact velocity
    peakOffset(ii,1) = ((log(cross_cor{ii,1}(I(ii,1)-1)))-(log(cross_cor{ii,1}(I(ii,1)+1))))/(2*((log(cross_cor{ii,1}(I(ii,1)-1)))+(log(cross_cor{ii,1}(I(ii,1)+1)))-2*(log(cc(ii,1)))));
    %find delta t
    DT(ii,1) = abs((lagDiff(ii,1)+peakOffset(ii,1))/FS(ii,1));
else
    lagDiff = NaN;
    DT = NaN;
end
%disturbance velocity, [m/s]
Uc(ii,1) = M(ii).dx2/1000/DT(ii,1);
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

%%% Overplot with BL Calcs
%Non-dim
Uc_star = Uc/mean(Uc(arrayfun(@(X) X == 3.24,[M.Ydist])));
Y_star = Y/3.24;
cfit1 = fit(Uc_star,Y_star,'power1');
ft2 = fittype(@(a,x) a*x.^(7));
cfit2 = fit(Uc_star,Y_star,ft2,'StartPoint',1);
ft3 = fittype(@(a,x) a*x.^(9));
cfit3 = fit(Uc_star,Y_star,ft3,'StartPoint',1);

fig7 = figure(7);
xlabel('Velocity [m/s]','FontSize',30)
ylabel('Height above floor [mm]','FontSize',30)
hold on
u = linspace(0.01,1,1001);
Uvec = u*mean(Uc(arrayfun(@(X) X == 3.24,[M.Ydist])));
h2_cf1 = plot(Uvec,3.24*25.4*cfit1(u),'k--','LineWidth',2);
h2_cf2 = plot(Uvec,3.24*25.4*cfit2(u),'k:','LineWidth',2);
h2_cf3 = plot(Uvec,3.24*25.4*cfit3(u),'k-.','LineWidth',2);
ylim([0,100]);
yline(3.24*25.4*cfit1(0.9900),'--')
text(300,5+3.24*25.4*cfit1(0.9900),['\delta_{99}=' sprintf('%.2fmm',3.24*25.4*cfit1(0.9900))],'FontSize',30)
grid on
h2 = scatter(Uc,25.4*Y,25,[M.Re]','filled');
cmap = colormap(jet);
%cmap = [1  1  1 ; cmap]; %hide low re
colormap(cmap)
cb2 = colorbar();
%caxis([1.4e6,inf])
caxis([-inf,inf])
cb2.Label.String = 'Re/x (m^{-1})';
C1 = coeffvalues(cfit1); C2 = coeffvalues(cfit2); C3 = coeffvalues(cfit3);
cfstring1 = ['$$y^*=' sprintf('%.2f',C1(1)) 'u^{* ' sprintf('%.2f',C1(2)) '}$$'];
cfstring2 = ['$$y^*=' sprintf('%.2f',C2(1)) 'u^{*7}$$'];
cfstring3 = ['$$y^*=' sprintf('%.2f',C3(1)) 'u^{*9}$$'];
l3 = legend([h2,h2_cf2,h2_cf3,h2_cf1],'Boundary Layer Profile',cfstring2,cfstring3,cfstring1,'Location','Northwest','Interpreter','latex','FontSize',30);

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


