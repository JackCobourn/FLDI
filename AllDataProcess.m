%% Find Matfiles
cd('D:\FLDI_UMD\Data')
Matfiles = dir('**/*.mat')';
cd('C:\Users\jcobourn\Documents\GitHub\Focused_Laser_Dif_interf')

%% Load Matfiles
%using matfiles() fuction, so they don't load totally into %memory

for i = 1:length(Matfiles)
Matfiles(i).Fullname = [Matfiles(i).folder,'\',Matfiles(i).name];
M(i).Matfiles = matfile(Matfiles(i).Fullname);
end
EarlyRunNumbers = [1 2 3 4 5 7 6 8] ;

%% Create Waterfall Plot
ham=100; %number of divisions for hamming windows
start = 0.1; %[s] trim window start %will eventuall calc dynammically from p data
stop=0.2; %[s] trim window stop

for i = 1:length(Matfiles)
    if any(i==1:8) %Add run number to early runs
        RunNumber(i,1) = EarlyRunNumbers(i);
    else
        RunNumber(i,1) = M(i).Matfiles.CampainRunNum;
    end
    Numpoints =     
    Y(i,1) = M(i).Matfiles.y_dist;
    VOLT(i,:) = [M(i).Matfiles.vmax,M(i).Matfiles.vmin];
    FS(i,1) = M(i).Matfiles.Fs;
    %CHA{i,1} = M(i).Matfiles.chA_run; lets avoid loading the whole thing
    CHA_TRIM{i,1} = M(i).Matfiles.chA_run(fix((FS(i).*start)):fix((FS(i).*stop)-1),1);        %CHA{i}(fix((FS(i).*start)):fix((FS(i).*stop)-1));
    CHB_TRIM{i,1} = M(i).Matfiles.chB_run(fix((FS(i).*start)):fix((FS(i).*stop)-1),1);
%     if any(i==[8:length(Matfiles)])
%         CHA_TRIM{i,1} = CHA_TRIM{i,1}.*1000-500;
%     end
end
[~,ordercode] = sort(RunNumber);%reordercode
% figure out how to deal with noise not recorded for 8A(i) = is(who(M(i).Matfiles,'chA_noise'))
% calc spectra

for i = 1:length(Matfiles)
    [PSDa{i,1}, f{i,1}] = pwelch(CHA_TRIM{i,1}-mean(CHA_TRIM{i,1}),ceil(length(CHA_TRIM{i,1})/ham),[],[],fix(FS(i,1)));
end

%loglog(f,PSDa,'linewidth',2)
fig1 = figure(1)

for i = 1:length(Matfiles)
    h(i) = waterfall(f{i,1},Y(i,1),PSDa{i,1}')
    hold on
end
ax1 = fig1.Children;
ax1.XScale = 'log';
%ax.YScale = 'log'
ax1.ZScale = 'log';
%ax1.View = [-37.5 30];
ax1.View = [0 0];

%% Create waterfall scatter plot
fig3 = figure(3);
for i = 1:length(Matfiles)
    h3(i) = scatter(f{i,1},Y(i,1).*ones(size(f{i,1})),15,PSDa{i,1}','filled');
    hold on
end
ax3 = fig3.Children;
ax3.XScale = 'log';
xlim([10e3 inf])
colormap colorcube
cb3 = colorbar();
cb3.Ruler.Scale = 'log';
cb3.Ruler.MinorTick = 'on';


%% Create BL prof
dx2 = 2.6e-3; %need to calc dynamically 
for i = 1:length(Matfiles)
[cross_cor{i,1},lag{i,1}] = xcorr(CHA_TRIM{i,1}-mean(CHA_TRIM{i,1}),-(CHB_TRIM{i,1}-mean(CHB_TRIM{i,1})), 'coeff');
%find max of xcorr and index
[cc{i,1},I{i,1}] = max(cross_cor{i,1});
%check if xcorr peak is high enough?
if cc{i,1}>0.3
    lagDiff(i,1) = lag{i,1}(I{i,1});
    %3 pt gaussian interpolation, to get a more exact velocity
    peakOffset(i,1) = ((log(cross_cor{i,1}(I{i,1}-1)))-(log(cross_cor{i,1}(I{i,1}+1))))/(2*((log(cross_cor{i,1}(I{i,1}-1)))+(log(cross_cor{i,1}(I{i,1}+1)))-2*(log(cc{i,1}))));
    %find delta t
    DT(i,1) = abs((lagDiff(i,1)+peakOffset(i,1))/FS(i,1));
else
    lagDiff = NaN;
    DT = NaN;
end
%disturbance velocity, [m/s]
Uc(i,1) = (dx2)/DT(i,1);

% figure(4)
% clf
% plot(lag,cross_cor)
end
fig2 = figure(2);
h2 = scatter(Uc(:,1),25.4*Y(:,1),'*');
hold on
xlabel('Velocity (m/s)')
ylabel('Height above floor (mm)')