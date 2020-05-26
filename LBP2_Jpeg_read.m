%% LBP2 Jpeg read
% Jack Cobourn cw 20200519
% Using: Newport LBP2-VIS2

BP_File = 'C:\Users\Jack\Desktop\20200526_BeamPic_estCenter_noLines_0001.jpeg';
BP = imread(BP_File);
BP2 = reshape(BP,[],3);
BP2Vals = unique(BP2,'rows');
BP = imcrop=[161 0 width 655]
BP(327:328,:,:) = repmat(0,[2,1207,3]);

[centers,radii] = imfindcircles(BP,[2,10],'Sensitivity',0.9,'EdgeThreshold',.3);
imshow(BP);
viscircles(centers, radii,'Color','b');

[~,I] = sort(centers,1);
sortedC=centers(I(:,1),:);