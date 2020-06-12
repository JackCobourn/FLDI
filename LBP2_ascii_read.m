%% LBP2 Jpeg read
% Jack Cobourn cw 20200601
% Using: Newport LBP2-VIS2

BP_File = 'F:\FLDI_SHC_202006\BP\20200612\20200612_PreRun1_0001.ascii.csv';
BP = csvread(BP_File);
imshow(BP,[])

[centers,radii] = imfindcircles(BP,[2,10],'Sensitivity',0.9,'EdgeThreshold',.3);
viscircles(centers, radii,'Color','b');

[~,I] = sort(centers,1);
sortedC=centers(I(:,1),:)';