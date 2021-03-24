function [cellmaskTmem119]=functionmaskeMasksOnClusterMoreTmemSIZE(thatpath)
savepath=fullfile('/broad/hptmp/lbinan/microglia/',thatpath);
mypath=fullfile(savepath,'/merfish_mosaics');
dapipath=fullfile(savepath,'/smFISH_mosaics');
mask_Tmem119=imread(fullfile(mypath,'mask_Tmem119.tif'));
mask_DAPI=imread(fullfile(dapipath,'mosaic_DAPI_3.tif'));
nuclei=bwareaopen(imbinarize(mask_DAPI,'adaptive','Sensitivity',0.001),200);
load(fullfile(savepath,'matlab.mat'));
%mask_Fcrls=imread(fullfile(mypath,'mask_Fcrls.tif'));
% mask_Gabbr1=imread(fullfile(mypath,'mask_Gabbr1.tif'));
% mask_Gabbr2=imread(fullfile(mypath,'mask_Gabbr2.tif'));

% SE=strel('disk',6);
% thisimage=imdilate(mask_Fcrls,SE);
% thisimage=imerode(thisimage,strel('disk',2));
% thisimage=bwareaopen(thisimage,420);
% thisimage=imdilate(thisimage,SE);
% %imshow(thisimage)
% cellmaskFcrls=thisimage;


%imwrite(imresize(cellmaskFcrls,10),fullfile(savepath,'analysis','cellmaskFcrlsMoreTmem.png'))

SE=strel('disk',6);
thisImage=bwareaopen(imdilate(mask_Tmem119,SE),750);
SE=strel('disk',3);
thisImage=imdilate(thisImage,SE);
thisImage=bwareaopen(thisImage,800);
thisImage=thisImage.*mask;
thisImage=bwareaopen(thisImage,750);
% figure, imshow(thisImage);
cellmaskTmem119=im2bw(thisImage);
mycompt=0;%%

    
    
%    removed=200*double((imdilate(removed,strel('disk',9))-removed));
%    removed=reduceImage(removed);
%    nuclei=reduceImage(nuclei);
%    savethisrgb=cat(3,removed,removed,nuclei);
 % cellmaskTmem119=thisImage;
%%
%imwrite(savethisrgb,fullfile(savepath,'analysis','removedonDAPI.png'))
%figure, imshow(thisImage);
%imshowpair(cellmaskTmem119,cellmaskFcrls);

% SE=strel('disk',4);
% thisImage=imdilate(mask_Gabbr1,SE).*imdilate(mask_Gabbr2,SE);
% thisImage=bwareaopen(thisImage,300);
% thisImage=imdilate(thisImage,strel('disk',3));
% thisImage=bwareaopen(thisImage,400);
% cellmaskGabbr1and2=thisImage;
% imshow(cellmaskGabbr1and2)
% stats=regionprops(cellmaskTmem119,'PixelIdxList');
% cellmaskallTmem_Gabbr1_2=zeros(size(cellmaskTmem119));
% for i=1:size(stats,1)
%     thiscell=zeros(size(cellmaskTmem119));
%     thiscell(stats(i).PixelIdxList)=1;
%     if sum(sum(thiscell.*double(imbinarize(mask_Gabbr1))))>10&& sum(sum(thiscell.*double(imbinarize(mask_Gabbr2))))>10
%         cellmaskallTmem_Gabbr1_2(stats(i).PixelIdxList)=1;
%     end
% end
% if contains(thatpath,'control1_brain2')
% % negacircle=200*(imdilate(cellmaskTmem119,strel('disk',13))-cellmaskTmem119);
% % posicircle=200*(imdilate(cellmaskallTmem_Gabbr1_2,strel('disk',9))-cellmaskallTmem_Gabbr1_2);
% % size(regionprops(cellmaskallTmem_Gabbr1_2))
% % size(regionprops(cellmaskTmem119))
% % 
% %  rgb=cat(3,negacircle,posicircle,mask_Tmem119);
%   imwrite(cellmaskTmem119,fullfile(savepath,'analysis','newsegmentation.tif'));
% end
%cellmaskallTmem_Gabbr1_2=imdilate(cellmaskTmem119.*cellmaskGabbr1and2,strel('disk',14));

%imwrite(imresize(cellmaskallTmem_Gabbr1_2,10),fullfile(savepath,'analysis','cellmaskallTmem_Gabbr1_2.png'))

%imshow(cellmaskallTmem_Gabbr1_2)
