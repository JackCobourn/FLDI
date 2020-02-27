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
    if M(ii).Matfiles.CampainRunNum == 27 || M(ii).Matfiles.CampainRunNum == 21 || M(ii).Matfiles.CampainRunNum == 12
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
start = readmatrix('D:\FLDI_UMD\DFLDI_UMD_UTSI_Campaign_M4_BL.xlsx','Sheet','Data','Range','K3:K32'); %[s] trim window start %will eventuall calc dynammically from p data
stop = readmatrix('D:\FLDI_UMD\DFLDI_UMD_UTSI_Campaign_M4_BL.xlsx','Sheet','Data','Range','L3:L32'); %[s] trim window stop
RunP0 = readmatrix('D:\FLDI_UMD\DFLDI_UMD_UTSI_Campaign_M4_BL.xlsx','Sheet','Data','Range','J3:J32');

for ii = 1:length(M)
    Details{ii,1} = whos(M(ii).Matfiles);     
    Y(ii,1) = M(ii).Matfiles.y_dist;
    VOLT(ii,:) = [M(ii).Matfiles.vmax(1,1),M(ii).Matfiles.vmin(1,1),M(ii).Matfiles.vmax(1,2),M(ii).Matfiles.vmin(1,2)];
    FS(ii,1) = M(ii).Matfiles.Fs;
    BIT(ii,1) = M(ii).Matfiles.bitRes;
    
    CHA_TRIM{ii,1} = M(ii).Matfiles.chA_run(fix((FS(ii).*start(ii))):fix((FS(ii).*stop(ii))-1),1);        %CHA{ii}(fix((FS(ii).*start)):fix((FS(ii).*stop)-1));
    CHB_TRIM{ii,1} = M(ii).Matfiles.chB_run(fix((FS(ii).*start(ii))):fix((FS(ii).*stop(ii))-1),1);
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
P = RunP0.*Pratio.*6894.76; %convert psi to pa
T = 300*Tratio;
rho = P./(R.*T);
a = sqrt(gamma*R*T);
U0 = Mach*a;
mu = 0.00001458*(T)^1.5/(T+110.4); %sotherlans law
Rex = rho.*U0./mu;
clear Mach gamma R Tratio Pratio P T rho a mu

%% Spectral Analysis
% 
% for ii = 1:length(M)
%     [PSDa{ii,1}, f{ii,1}] = pwelch(CHA_TRIM{ii,1}-mean(CHA_TRIM{ii,1}),ceil(length(CHA_TRIM{ii,1})/ham),[],[],fix(FS(ii,1)));
%     %calc PSD area
%     PSDArea(ii,1) = trapz(PSDa{ii,1}).*mean(diff(f{ii,1}));
%     PSDNorm{ii,1} = PSDa{ii,1}./PSDArea(ii,1); %normalize to area of 1
%         
% %         if BIT(ii) == 15 %account for fs difference
% %             PSDNorm{ii,1} = PSDNorm{ii,1}*(20/12.5);
% %         elseif BIT(ii) == 8
% %             PSDNorm{ii,1} = PSDNorm{ii,1}*(12.5/20);
% %         end
%         NormResults(ii,1) = trapz(PSDNorm{ii,1}(1:find(f{ii,1}>=f{1,1}(end),1)));
% end
% 
%% Create IMageSc waterfall plot
[~,Index] = sort(Y(:,1));
for ii = 1:length(M)
PSDedit{ii,1} = PSDNorm{ii,1}(1:length(PSDa{1,1}))';
end
%PSDedit(:) = PSDa(:)(1:length(PSDa{1,1}))'; is there a way to do above in
%1 line? could use cell2mat(PSDa) and then remap()
Colors = cell2mat(PSDedit);

fig3 = figure(3);
imagesc('CData',log10(Colors(Index,:)))
ax3 = fig3.Children;
%xlim([1000,6000] )
xlim([0,inf])
ylim([0.5,length(M)+0.5])
colormap jet
 cb3 = colorbar();
 tickvals = cb3.TickLabels;
 for ii = 1:length(tickvals)
     tickname{ii} = sprintf('10^{%d}',str2double(tickvals{ii}));
 end
 cb3.TickLabels = tickname;
 for ii = 1:length(Matfiles)
     text(100,ii,sprintf('y = %.2f, %d bit',Y(Index(ii)),BIT(Index(ii))))
 end
 
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
dx2_opts = readmatrix('D:\FLDI_UMD\DFLDI_UMD_UTSI_Campaign_M4_BL.xlsx','Sheet','Beam Pics','Range','K2:K10'); %need to calc dynamically 
BeamPicRel = readmatrix('D:\FLDI_UMD\DFLDI_UMD_UTSI_Campaign_M4_BL.xlsx','Sheet','Data','Range','G3:G32');
dx2 = dx2_opts(BeamPicRel);

for ii = 1:length(M)
[cross_cor{ii,1},lag{ii,1}] = xcorr(CHA_TRIM{ii,1}-mean(CHA_TRIM{ii,1}),-(CHB_TRIM{ii,1}-mean(CHB_TRIM{ii,1})), 'coeff');
%find max of xcorr and index
[cc{ii,1},I{ii,1}] = max(cross_cor{ii,1});
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
h2 = scatter(Uc(:,1),25.4*Y(:,1),25,Rex,'filled');
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