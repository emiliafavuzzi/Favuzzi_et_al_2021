function [cellmaskTmem119]=functionmaskeMasksOnClusterMoreTmem(thatpath)
savepath=fullfile('/broad/hptmp/lbinan/microglia/',thatpath);
mypath=fullfile(savepath,'/merfish_mosaics');
mask_Tmem119=imread(fullfile(mypath,'mask_Tmem119.tif'));
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
% figure, imshow(thisImage);
cellmaskTmem119=im2bw(thisImage);
%%
disp('here')
stats=regionprops(cellmaskTmem119,'area','PixelIdxList');
  thisImage=zeros(size(cellmaskTmem119));
    for i=1:size(stats,1)
    thiscell=zeros(size(cellmaskTmem119));
    thiscell(stats(i).PixelIdxList)=1;
    thisarea=stats(i).Area;
%     while thisarea<1700
%         thiscell=imdilate(thiscell,strel('disk',4));
%         thisarea=sum(sum(im2bw(thiscell)));
%     end
    while thisarea>3000
        thiscell=imerode(thiscell,strel('disk',4));
        newstats=regionprops(im2bw(bwareaopen(thiscell,400)),'area');
        thisarea=newstats(1).Area
        disp(thisarea)
    end 
    disp('cell done')
    thisImage=thisImage+thiscell;
    end
 cellmaskTmem119=thisImage;
