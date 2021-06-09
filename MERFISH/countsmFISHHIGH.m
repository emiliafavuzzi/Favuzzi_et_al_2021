function [mythreshold]=countsmFISHHIGH(thatpath)
% Loic Binan
%lbinan@broadinstitute.org
%3/26/2021
%this functions edits the counts generate by the main script to replace intensity values with counts for C1qc, Tms4bx,Rsp29,Ftl1
% thatpath='KO2_brain2/slice2_side1';
savepath=fullfile('/broad/hptmp/lbinan/microglia/',thatpath);
mypath=fullfile(savepath,'/merfish_mosaics');
dapipath=fullfile(savepath,'/smFISH_mosaics');
% mask_Tmem119=imread(fullfile(mypath,'mask_Tmem119.tif'));
mask_DAPI=imread(fullfile(dapipath,'mosaic_DAPI_3.tif'));
nuclei=bwareaopen(imbinarize(mask_DAPI,'adaptive','Sensitivity',0.001),200);

% nuclei=bwareaopen(imbinarize(mask_DAPI,'adaptive','Sensitivity',0.001),200);
% load(fullfile(savepath,'matlab.mat'));
csvfiles = dir(strcat(savepath,'/*astIntensitiesallZ.csv'));
counttable=table2array(readtable(fullfile(csvfiles(1).folder,csvfiles(1).name)));
% csvfiles = dir(strcat(savepath,'/*fastIntensitiesallZ.csv'));
% intensitytable=table2array(readtable(fullfile(csvfiles(1).folder,csvfiles(1).name)));
newcounts=counttable;
newcounts(:,5)=0;
newcounts(:,6)=0;
newcounts(:,7)=0;
newcounts(:,8)=0;
C1qc=newcounts(:,5);
%   thisImage=zeros(size(cellmaskTmem119));
  mask_Tmem119_0=imread(fullfile(mypath,'Tmem119_0.tif'));
  mask_Tmem119_1=imread(fullfile(mypath,'Tmem119_1.tif'));
  mask_Tmem119_2=imread(fullfile(mypath,'Tmem119_2.tif'));
  mask_Tmem119_3=imread(fullfile(mypath,'Tmem119_3.tif'));
  mask_Tmem119_4=imread(fullfile(mypath,'Tmem119_4.tif'));
  mask_Tmem119_5=imread(fullfile(mypath,'Tmem119_5.tif'));
  mask_Tmem119_6=imread(fullfile(mypath,'Tmem119_6.tif'));
  mask_Tmem119=im2bw(mask_Tmem119_0+mask_Tmem119_1+mask_Tmem119_2+mask_Tmem119_3+mask_Tmem119_4+mask_Tmem119_5+mask_Tmem119_6);
  SE=strel('disk',6);
thisImage=bwareaopen(imdilate(mask_Tmem119,SE),750);
SE=strel('disk',3);
thisImage=imdilate(thisImage,SE);
thisImage=bwareaopen(thisImage,800);
% thisImage=thisImage.*mask;
thisImage=bwareaopen(thisImage,750);
% figure, imshow(thisImage);
cellmaskTmem119=im2bw(thisImage);
cellMask=cellmaskTmem119;
mycompt=0;%%
% disp('here')
stats=regionprops(cellmaskTmem119,'area','PixelIdxList');
disp('here')
% while(size(C1qc,1)<size(stats,1))
%     C1qc=[C1qc;0];
% end
  Gabbr2_1=imread(fullfile(dapipath,'mosaic_bit13_compiled.tif'));
  Tspan7_1=imread(fullfile(dapipath,'mosaic_bit15_compiled.tif'));

  mask_C1qc=imread(fullfile(dapipath,'mosaic_bit17_compiled.tif'));
  mask_Tms4bx=imread(fullfile(dapipath,'mosaic_bit18_compiled.tif'));
  mask_Rps29=imread(fullfile(dapipath,'mosaic_bit19_compiled.tif'));
  mask_Ftl1=imread(fullfile(dapipath,'mosaic_bit20_compiled.tif'));

mylastGood=1;
thiscount=1;
i=1;
mypixels=find(nuclei==1);
C1qcThresh=mean(prctile(mask_C1qc(mypixels),95));
Tms4bxThresh=mean(prctile(mask_Tms4bx(mypixels),95));
Rps29Thresh=mean(prctile(mask_Rps29(mypixels),95));
Ftl1Thresh=mean(prctile(mask_Ftl1(mypixels),95));

disp('there')
while i<=size(stats,1)
         disp(i)
    thiscell=zeros(size(cellMask));
    thiscell(stats(i).PixelIdxList)=1;
    thisImage=bsxfun(@times, Gabbr2_1, cast(thiscell, 'like', Gabbr2_1));
    thisintensity=sum(sum(thisImage));
    if thisintensity==counttable(thiscount,1)
    thisImage=bsxfun(@times, Tspan7_1, cast(thiscell, 'like', Tspan7_1));
    thisintensity=sum(sum(thisImage));  
    disp('almost')
       if thisintensity==counttable(thiscount,3) 
           disp('found one')
           thisImage=bsxfun(@times, mask_C1qc, cast(thiscell, 'like', mask_C1qc));
        C1qc=size(pkfnd(thisImage,C1qcThresh,4),1);
        thisImage=bsxfun(@times, mask_Tms4bx, cast(thiscell, 'like', mask_Tms4bx));
        Tms4bx=size(pkfnd(thisImage,Tms4bxThresh,4),1);
        thisImage=bsxfun(@times, mask_Rps29, cast(thiscell, 'like', mask_Rps29));
        Rps29=size(pkfnd(thisImage,Rps29Thresh,4),1);
        thisImage=bsxfun(@times, mask_Ftl1, cast(thiscell, 'like', mask_Ftl1));
        Ftl1=size(pkfnd(thisImage,Ftl1Thresh,4),1);
        newcounts(thiscount,5)=C1qc;
        newcounts(thiscount,6)=Tms4bx;
        newcounts(thiscount,7)=Rps29;
        newcounts(thiscount,8)=Ftl1;
        thiscount=thiscount+1;
         mylastGood=i;
    end
    end

    if thiscount<size(counttable,1)&& i==size(stats,1)
        i=mylastGood+1;
        thiscount=thiscount+1;
    else
        i=i+1;
    end
    if thiscount>size(counttable,1)
        break
    end
end
writematrix(newcounts,fullfile(savepath,strcat(thatpath(17:end),'counted4genesNUCLEITHRESHHIGH.csv')))

% writematrix(newintensities,fullfile(savepath,strcat(thatpath(17:end),'fastIntensitiesFixedZCounts.csv')));

