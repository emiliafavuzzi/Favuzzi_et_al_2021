function [mythreshold]=recountwithZ(thatpath)
%thatpath='control1_brain2/slice1_side1';
savepath=fullfile('/broad/hptmp/lbinan/microglia/',thatpath);
mypath=fullfile(savepath,'/merfish_mosaics');
dapipath=fullfile(savepath,'/smFISH_mosaics');
% mask_Tmem119=imread(fullfile(mypath,'mask_Tmem119.tif'));
mask_DAPI=imread(fullfile(dapipath,'mosaic_DAPI_3.tif'));
% nuclei=bwareaopen(imbinarize(mask_DAPI,'adaptive','Sensitivity',0.001),200);
load(fullfile(savepath,'matlab.mat'));
csvfiles = dir(strcat(savepath,'/*fastCounts.csv'));
counttable=table2array(readtable(fullfile(csvfiles(1).folder,csvfiles(1).name)));
% csvfiles = dir(strcat(savepath,'/*fastIntensitiesallZ.csv'));
% intensitytable=table2array(readtable(fullfile(csvfiles(1).folder,csvfiles(1).name)));
newcounts=counttable;
newcounts(:,2)=0;
newcounts(:,3)=0;
Gabbr1_0=newcounts(:,2);
Gabbr1_1=newcounts(:,2);
Gabbr1_2=newcounts(:,2);
Gabbr1_3=newcounts(:,2);
Gabbr1_4=newcounts(:,2);
Gabbr1_5=newcounts(:,2);
Gabbr1_6=newcounts(:,2);
Gabbr2_0=newcounts(:,3);
Gabbr2_1=newcounts(:,3);
Gabbr2_2=newcounts(:,3);
Gabbr2_3=newcounts(:,3);
Gabbr2_4=newcounts(:,3);
Gabbr2_5=newcounts(:,3);
Gabbr2_6=newcounts(:,3);
numberofZ=newcounts(:,3);

% newintensities=[];
mythreshold=7;


%   thisImage=zeros(size(cellmaskTmem119));
  mask_Tmem119_0=imread(fullfile(mypath,'Tmem119_0.tif'));
  mask_Tmem119_1=imread(fullfile(mypath,'Tmem119_1.tif'));
  mask_Tmem119_2=imread(fullfile(mypath,'Tmem119_2.tif'));
  mask_Tmem119_3=imread(fullfile(mypath,'Tmem119_3.tif'));
  mask_Tmem119_4=imread(fullfile(mypath,'Tmem119_4.tif'));
  mask_Tmem119_5=imread(fullfile(mypath,'Tmem119_5.tif'));
  mask_Tmem119_6=imread(fullfile(mypath,'Tmem119_6.tif'));
  mask_Tmem119=(mask_Tmem119_0+mask_Tmem119_1+mask_Tmem119_2+mask_Tmem119_3+mask_Tmem119_4+mask_Tmem119_5+mask_Tmem119_6);
  SE=strel('disk',6);
thisImage=bwareaopen(imdilate(mask_Tmem119,SE),750);
SE=strel('disk',3);
thisImage=imdilate(thisImage,SE);
thisImage=bwareaopen(thisImage,800);
thisImage=thisImage.*mask;
thisImage=bwareaopen(thisImage,750);
% figure, imshow(thisImage);
cellmaskTmem119=im2bw(thisImage);
% mycompt=0;%%
% disp('here')
stats=regionprops(cellmaskTmem119,'area','PixelIdxList');
  mask_Gabbr1_0=imread(fullfile(mypath,'Gabbr1_0.tif'));
  mask_Gabbr1_1=imread(fullfile(mypath,'Gabbr1_1.tif'));
  mask_Gabbr1_2=imread(fullfile(mypath,'Gabbr1_2.tif'));
  mask_Gabbr1_3=imread(fullfile(mypath,'Gabbr1_3.tif'));
  mask_Gabbr1_4=imread(fullfile(mypath,'Gabbr1_4.tif'));
  mask_Gabbr1_5=imread(fullfile(mypath,'Gabbr1_5.tif'));
  mask_Gabbr1_6=imread(fullfile(mypath,'Gabbr1_6.tif'));
  mask_Gabbr2_0=imread(fullfile(mypath,'Gabbr2_0.tif'));
  mask_Gabbr2_1=imread(fullfile(mypath,'Gabbr2_1.tif'));
  mask_Gabbr2_2=imread(fullfile(mypath,'Gabbr2_2.tif'));
  mask_Gabbr2_3=imread(fullfile(mypath,'Gabbr2_3.tif'));
  mask_Gabbr2_4=imread(fullfile(mypath,'Gabbr2_4.tif'));
  mask_Gabbr2_5=imread(fullfile(mypath,'Gabbr2_5.tif'));
  mask_Gabbr2_6=imread(fullfile(mypath,'Gabbr2_6.tif'));
  mask_Fcrls_0=imread(fullfile(mypath,'Fcrls_0.tif'));
  mask_Fcrls_1=imread(fullfile(mypath,'Fcrls_1.tif'));
  mask_Fcrls_2=imread(fullfile(mypath,'Fcrls_2.tif'));
  mask_Fcrls_3=imread(fullfile(mypath,'Fcrls_3.tif'));
  mask_Fcrls_4=imread(fullfile(mypath,'Fcrls_4.tif'));
  mask_Fcrls_5=imread(fullfile(mypath,'Fcrls_5.tif'));
  mask_Fcrls_6=imread(fullfile(mypath,'Fcrls_6.tif'));
   mask_P2ry12_0=imread(fullfile(mypath,'P2ry12_0.tif'));
  mask_P2ry12_1=imread(fullfile(mypath,'P2ry12_1.tif'));
  mask_P2ry12_2=imread(fullfile(mypath,'P2ry12_2.tif'));
  mask_P2ry12_3=imread(fullfile(mypath,'P2ry12_3.tif'));
  mask_P2ry12_4=imread(fullfile(mypath,'P2ry12_4.tif'));
  mask_P2ry12_5=imread(fullfile(mypath,'P2ry12_5.tif'));
  mask_P2ry12_6=imread(fullfile(mypath,'P2ry12_6.tif'));
%      
%     for i=1:size(stats,1)
%     for i=1:150
parfor (i=1:size(stats,1),12)
    thiscell=zeros(size(cellmaskTmem119));
    thiscell(stats(i).PixelIdxList)=1;
%     TmemthisZ=bsxfun(@times, mask_Tmem119, cast(thiscell, 'like', mask_Tmem119));
%     mynucleus=sum(sum(thisnucleus));
Z0microglia=(sum(sum(bsxfun(@times, mask_Tmem119_0, cast(thiscell, 'like', mask_Tmem119_0))))+sum(sum(bsxfun(@times, mask_P2ry12_0, cast(thiscell, 'like', mask_P2ry12_0))))+sum(sum(bsxfun(@times, mask_Fcrls_0, cast(thiscell, 'like', mask_Fcrls_0)))))/255;
Z1microglia=(sum(sum(bsxfun(@times, mask_Tmem119_1, cast(thiscell, 'like', mask_Tmem119_1))))+sum(sum(bsxfun(@times, mask_P2ry12_1, cast(thiscell, 'like', mask_P2ry12_1))))+sum(sum(bsxfun(@times, mask_Fcrls_1, cast(thiscell, 'like', mask_Fcrls_1)))))/255;
Z2microglia=(sum(sum(bsxfun(@times, mask_Tmem119_2, cast(thiscell, 'like', mask_Tmem119_2))))+sum(sum(bsxfun(@times, mask_P2ry12_2, cast(thiscell, 'like', mask_P2ry12_2))))+sum(sum(bsxfun(@times, mask_Fcrls_2, cast(thiscell, 'like', mask_Fcrls_2)))))/255;
Z3microglia=(sum(sum(bsxfun(@times, mask_Tmem119_3, cast(thiscell, 'like', mask_Tmem119_3))))+sum(sum(bsxfun(@times, mask_P2ry12_3, cast(thiscell, 'like', mask_P2ry12_3))))+sum(sum(bsxfun(@times, mask_Fcrls_3, cast(thiscell, 'like', mask_Fcrls_3)))))/255;
Z4microglia=(sum(sum(bsxfun(@times, mask_Tmem119_4, cast(thiscell, 'like', mask_Tmem119_4))))+sum(sum(bsxfun(@times, mask_P2ry12_4, cast(thiscell, 'like', mask_P2ry12_4))))+sum(sum(bsxfun(@times, mask_Fcrls_4, cast(thiscell, 'like', mask_Fcrls_4)))))/255;
Z5microglia=(sum(sum(bsxfun(@times, mask_Tmem119_5, cast(thiscell, 'like', mask_Tmem119_5))))+sum(sum(bsxfun(@times, mask_P2ry12_5, cast(thiscell, 'like', mask_P2ry12_5))))+sum(sum(bsxfun(@times, mask_Fcrls_5, cast(thiscell, 'like', mask_Fcrls_5)))))/255;
Z6microglia=(sum(sum(bsxfun(@times, mask_Tmem119_6, cast(thiscell, 'like', mask_Tmem119_6))))+sum(sum(bsxfun(@times, mask_P2ry12_6, cast(thiscell, 'like', mask_P2ry12_6))))+sum(sum(bsxfun(@times, mask_Fcrls_6, cast(thiscell, 'like', mask_Fcrls_6)))))/255;
% disp('Z0microglia');disp(Z0microglia)
% disp('Z1microglia');disp(Z1microglia)
% disp('Z2microglia');disp(Z2microglia)
% disp('Z3microglia');disp(Z3microglia)
% disp('Z4microglia');disp(Z4microglia)
% disp('Z5microglia');disp(Z5microglia)
% disp('Z6microglia');disp(Z6microglia)
Z0Gabbr1=(sum(sum(bsxfun(@times, mask_Gabbr1_0, cast(thiscell, 'like', mask_Gabbr1_0)))))/255;
Z1Gabbr1=(sum(sum(bsxfun(@times, mask_Gabbr1_1, cast(thiscell, 'like', mask_Gabbr1_1)))))/255;
Z2Gabbr1=(sum(sum(bsxfun(@times, mask_Gabbr1_2, cast(thiscell, 'like', mask_Gabbr1_2)))))/255;
Z3Gabbr1=(sum(sum(bsxfun(@times, mask_Gabbr1_3, cast(thiscell, 'like', mask_Gabbr1_3)))))/255;
Z4Gabbr1=(sum(sum(bsxfun(@times, mask_Gabbr1_4, cast(thiscell, 'like', mask_Gabbr1_4)))))/255;
Z5Gabbr1=(sum(sum(bsxfun(@times, mask_Gabbr1_5, cast(thiscell, 'like', mask_Gabbr1_5)))))/255;
Z6Gabbr1=(sum(sum(bsxfun(@times, mask_Gabbr1_6, cast(thiscell, 'like', mask_Gabbr1_6)))))/255;
Z0Gabbr2=(sum(sum(bsxfun(@times, mask_Gabbr2_0, cast(thiscell, 'like', mask_Gabbr2_0)))))/255;
Z1Gabbr2=(sum(sum(bsxfun(@times, mask_Gabbr2_1, cast(thiscell, 'like', mask_Gabbr2_1)))))/255;
Z2Gabbr2=(sum(sum(bsxfun(@times, mask_Gabbr2_2, cast(thiscell, 'like', mask_Gabbr2_2)))))/255;
Z3Gabbr2=(sum(sum(bsxfun(@times, mask_Gabbr2_3, cast(thiscell, 'like', mask_Gabbr2_3)))))/255;
Z4Gabbr2=(sum(sum(bsxfun(@times, mask_Gabbr2_4, cast(thiscell, 'like', mask_Gabbr2_4)))))/255;
Z5Gabbr2=(sum(sum(bsxfun(@times, mask_Gabbr2_5, cast(thiscell, 'like', mask_Gabbr2_5)))))/255;
Z6Gabbr2=(sum(sum(bsxfun(@times, mask_Gabbr2_6, cast(thiscell, 'like', mask_Gabbr2_6)))))/255;
% disp('Z0Gabbr1');disp(Z0Gabbr1)
% disp('Z1Gabbr1');disp(Z1Gabbr1)
% disp('Z2Gabbr1');disp(Z2Gabbr1)
% disp('Z3Gabbr1');disp(Z3Gabbr1)
% disp('Z4Gabbr1');disp(Z4Gabbr1)
% disp('Z5Gabbr1');disp(Z5Gabbr1)
% disp('Z6Gabbr1');disp(Z6Gabbr1)
% newGabbr2=0;
% newGabbr1=0;
count=0;
    if Z0microglia>=mythreshold
        Gabbr1_0(i+1)=Z0Gabbr1;
        Gabbr2_0(i+1)=Z0Gabbr2;
        count=count+1;
    end
    if Z1microglia>=mythreshold
        Gabbr1_1(i+1)=Z1Gabbr1;
        Gabbr2_1(i+1)=Z1Gabbr2;
        count=count+1;
    end    
    if Z2microglia>=mythreshold
        Gabbr1_2(i+1)=Z2Gabbr1;
        Gabbr2_2(i+1)=Z2Gabbr2;
        count=count+1;
    end    
    if Z3microglia>=mythreshold
        Gabbr1_3(i+1)=Z3Gabbr1;
        Gabbr2_3(i+1)=Z3Gabbr2;
        count=count+1;
    end    
    if Z4microglia>=mythreshold
        Gabbr1_4(i+1)=Z4Gabbr1;
        Gabbr2_4(i+1)=Z4Gabbr2;
        count=count+1;
    end    
    if Z5microglia>=mythreshold
        Gabbr1_5(i+1)=Z5Gabbr1;
        Gabbr2_5(i+1)=Z5Gabbr2;
        count=count+1;
    end    
    if Z6microglia>=mythreshold
        Gabbr1_6(i+1)=+Z6Gabbr1;
        Gabbr2_6(i+1)=Z6Gabbr2;
        count=count+1;
    end
    numberofZ(i+1)=count;
end
    newcounts(:,2)=Gabbr1_0+Gabbr1_1+Gabbr1_2+Gabbr1_3+Gabbr1_4+Gabbr1_5+Gabbr1_6;
     newcounts(:,3)=Gabbr2_0+Gabbr2_1+Gabbr2_2+Gabbr2_3+Gabbr2_4+Gabbr2_5+Gabbr2_6;
      
writematrix([newcounts,numberofZ],fullfile(savepath,strcat(thatpath(17:end),'fastCountsFixedZCounts.csv')));
% writematrix(newintensities,fullfile(savepath,strcat(thatpath(17:end),'fastIntensitiesFixedZCounts.csv')));

