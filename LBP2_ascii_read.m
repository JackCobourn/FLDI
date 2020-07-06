%% LBP2 Jpeg read
% Jack Cobourn cw 20200601
% Using: Newport LBP2-VIS2

[BP_File,BP_path] = uigetfile({'*.csv','Comma Separated Values)'},'Select a File');
BP = csvread([BP_path filesep BP_File]);
imshow(BP,[])

[centers,radii] = imfindcircles(BP,[2,10],'Sensitivity',0.9,'EdgeThreshold',.3);
viscircles(centers, radii,'Color','b');

[~,I] = sort(centers,1);
sortedC=centers(I(:,1),:)';

dx = sortedC(1,:);
dy = sortedC(2,:);
dx1a = 0.00738*(dx(2)-dx(1))/1000; %FLDI beam separation
dy1a = 0.00738*(dy(2)-dy(1))/1000;
dx1b =  0.00738*(dx(4)-dx(3))/1000;
dy1b = 0.00738*(dy(4)-dy(3))/1000;
dx2 =  0.00738*(mean(dx([3,4]))-mean(dx([1,2])))/1000; %[m] 
dy2 = 0.00738*(mean(dy([3,4]))-mean(dy([1,2])))/1000; %[m] 
d2 = sqrt(dx2.^2+dy2.^2);