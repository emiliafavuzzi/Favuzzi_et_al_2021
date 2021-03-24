function fastintensityCompiled(thatpath,cellmaskTmem119)
SE=strel('disk',8);
savepath=fullfile('/broad/hptmp/lbinan/microglia/',thatpath);
mypath=fullfile(savepath,'/smFISH_mosaics');
%%
%cellmaskallTmem_Gabbr1_2=imread(fullfile(savepath,'analysis','cellmaskallTmem_Gabbr1_2.png'));
%cellmaskTmem119=imread(fullfile(savepath,'analysis','cellmaskTmem119.png'));
%cellMask=im2bw(imdilate(cellmaskallTmem_Gabbr1_2,SE));
cellMask=im2bw(cellmaskTmem119);
%   SE=strel('disk',2);
%  white=ones(size(cellMask));
%   negativeCells=bwareaopen(cellmaskTmem119.*(white-imdilate(cellMask,SE)),650);
% %imshow(negativeCells);
% %%
% cellMask=im2bw(negativeCells);
%imshow(cellMask);
Gabbr2_1_3=(imread(fullfile(mypath,'mosaic_bit13_compiled.tif')));
%imwrite(imresize(Gabbr2_1_3,10),fullfile(mypath,'filteredsmFISHGabbr2_1.png'));
Gabbr2_2_3=(imread(fullfile(mypath,'mosaic_bit14_compiled.tif')));
%imwrite(imresize(Gabbr2_2_3,10),fullfile(mypath,'filteredsmFISHGabbr2_2.png'));
Tspan7_1_3=(imread(fullfile(mypath,'mosaic_bit15_compiled.tif')));
%imwrite(imresize(Tspan7_1_3,10),fullfile(mypath,'filteredsmFISHTspan7_1.png'));
Tspan7_2_3=(imread(fullfile(mypath,'mosaic_bit16_compiled.tif')));
%imwrite(imresize(Tspan7_2_3,10),fullfile(mypath,'filteredsmFISHTspan7_2.png'));
C1qc_3=(imread(fullfile(mypath,'mosaic_bit17_compiled.tif')));
%imwrite(imresize(C1qc_3,10),fullfile(mypath,'filteredsmFISHC1qc.png'));
Tms4bx_3=(imread(fullfile(mypath,'mosaic_bit18_compiled.tif')));
%imwrite(imresize(Tms4bx_3,10),fullfile(mypath,'filteredsmFISHTms4bx.png'));
Rps29_3=(imread(fullfile(mypath,'mosaic_bit19_compiled.tif')));
%imwrite(imresize(Rps29_3,10),fullfile(mypath,'filteredsmFISHRps29.png'));
Ftl1_3=(imread(fullfile(mypath,'mosaic_bit20_compiled.tif')));
%imwrite(imresize(Ftl1_3,10),fullfile(mypath,'filteredsmFISHFtl1.png'));
mt_Nd2_3=(imread(fullfile(mypath,'mosaic_bit21_compiled.tif')));
%imwrite(imresize(mt_Nd2_3,10),fullfile(mypath,'filteredsmFISHmt_Nd2.png'));
mt_Co3_3=(imread(fullfile(mypath,'mosaic_bit22_compiled.tif')));
%imwrite(imresize(mt_Co3_3,10),fullfile(mypath,'filteredsmFISHmt_Co3.png'));

stats=regionprops(cellMask,'image','area','PixelIdxList');
thiscell=zeros(size(cellMask));
intensity4Gabbr2_1=zeros([1+size(stats,1),1]);
intensity4Gabbr2_1(1)="Gabbr2_1";
intensity4Gabbr2_2=zeros([1+size(stats,1),1]);
intensity4Gabbr2_2(1)="Gabbr2_2";
intensity4Tspan7_1=zeros([1+size(stats,1),1]);
intensity4Tspan7_1(1)="Tspan7_1";
intensity4Tspan7_2=zeros([1+size(stats,1),1]);
intensity4Tspan7_2(1)="Tspan7_2";
intensity4C1qc=zeros([1+size(stats,1),1]);
intensity4C1qc(1)="C1qc";
intensity4Tms4bx=zeros([1+size(stats,1),1]);
intensity4Tms4bx(1)="Tms4bx";
intensity4Rps29=zeros([1+size(stats,1),1]);
intensity4Rps29(1)="Rps29";
intensity4Ftl1=zeros([1+size(stats,1),1]);
intensity4Ftl1(1)="Ftl1";
intensity4mt_Nd2=zeros([1+size(stats,1),1]);
intensity4mt_Nd2(1)="mt_Nd2";
intensity4mt_Co3=zeros([1+size(stats,1),1]);
intensity4mt_Co3(1)="mt_Co3";
normalizedintensity4Gabbr2_1=zeros([1+size(stats,1),1]);
normalizedintensity4Gabbr2_1(1)="Gabbr2_1";
normalizedintensity4Gabbr2_2=zeros([1+size(stats,1),1]);
normalizedintensity4Gabbr2_2(1)="Gabbr2_2";
normalizedintensity4Tspan7_1=zeros([1+size(stats,1),1]);
normalizedintensity4Tspan7_1(1)="Tspan7_1";
normalizedintensity4Tspan7_2=zeros([1+size(stats,1),1]);
normalizedintensity4Tspan7_2(1)="Tspan7_2";
normalizedintensity4C1qc=zeros([1+size(stats,1),1]);
normalizedintensity4C1qc(1)="C1qc";
normalizedintensity4Tms4bx=zeros([1+size(stats,1),1]);
normalizedintensity4Tms4bx(1)="Tms4bx";
normalizedintensity4Rps29=zeros([1+size(stats,1),1]);
normalizedintensity4Rps29(1)="Rps29";
normalizedintensity4Ftl1=zeros([1+size(stats,1),1]);
normalizedintensity4Ftl1(1)="Ftl1";
normalizedintensity4mt_Nd2=zeros([1+size(stats,1),1]);
normalizedintensity4mt_Nd2(1)="mt_Nd2";
normalizedintensity4mt_Co3=zeros([1+size(stats,1),1]);
normalizedintensity4mt_Co3(1)="mt_Co3";


parfor (i=1:size(stats,1),16)
    thiscell=zeros(size(cellMask));
    thiscell(stats(i).PixelIdxList)=1;
%     cellsize=stats(i).Area;
%     if stats(i).Area<2500
%         thiscell=imdilate(thiscell,strel('disk',4));
%     end
    thisImage=bsxfun(@times, Gabbr2_1_3, cast(thiscell, 'like', Gabbr2_1_3));
    thisintensity=sum(sum(thisImage));
    intensity4Gabbr2_1(i+1)=(thisintensity);normalizedintensity4Gabbr2_1(i+1)=(thisintensity/sum(sum(thiscell)));
    thisImage=bsxfun(@times, Gabbr2_2_3, cast(thiscell, 'like', Gabbr2_2_3));
    thisintensity=sum(sum(thisImage));
    intensity4Gabbr2_2(i+1)=(thisintensity);normalizedintensity4Gabbr2_2(i+1)=(thisintensity/sum(sum(thiscell)));
    thisImage=bsxfun(@times, C1qc_3, cast(thiscell, 'like', C1qc_3));
    thisintensity=sum(sum(thisImage));
    intensity4C1qc(i+1)=(thisintensity);    normalizedintensity4C1qc(i+1)=(thisintensity/sum(sum(thiscell)));
     thisImage=bsxfun(@times, Ftl1_3, cast(thiscell, 'like', Ftl1_3));
    thisintensity=sum(sum(thisImage));
    intensity4Ftl1(i+1)=(thisintensity);normalizedintensity4Ftl1(i+1)=(thisintensity/sum(sum(thiscell)));
     thisImage=bsxfun(@times, mt_Nd2_3, cast(thiscell, 'like', mt_Nd2_3));
    thisintensity=sum(sum(thisImage));
    intensity4mt_Nd2(i+1)=(thisintensity);normalizedintensity4mt_Nd2(i+1)=(thisintensity/sum(sum(thiscell)));
         thisImage=bsxfun(@times, mt_Co3_3, cast(thiscell, 'like', mt_Co3_3));
    thisintensity=sum(sum(thisImage));
    intensity4mt_Co3(i+1)=(thisintensity);normalizedintensity4mt_Co3(i+1)=(thisintensity/sum(sum(thiscell)));
    thisImage=bsxfun(@times, Rps29_3, cast(thiscell, 'like', Rps29_3));
    thisintensity=sum(sum(thisImage));
    intensity4Rps29(i+1)=(thisintensity);normalizedintensity4Rps29(i+1)=(thisintensity/sum(sum(thiscell)));
    thisImage=bsxfun(@times, Tms4bx_3, cast(thiscell, 'like', Tms4bx_3));
    thisintensity=sum(sum(thisImage));
    intensity4Tms4bx(i+1)=(thisintensity);  normalizedintensity4Tms4bx(i+1)=(thisintensity/sum(sum(thiscell)));
    thisImage=bsxfun(@times, Tspan7_1_3, cast(thiscell, 'like', Tspan7_1_3));
    thisintensity=sum(sum(thisImage));
    intensity4Tspan7_1(i+1)=(thisintensity);normalizedintensity4Tspan7_1(i+1)=(thisintensity/sum(sum(thiscell)));
    thisImage=bsxfun(@times, Tspan7_2_3, cast(thiscell, 'like', Tspan7_2_3));
    thisintensity=sum(sum(thisImage));
    intensity4Tspan7_2(i+1)=(thisintensity);normalizedintensity4Tspan7_2(i+1)=(thisintensity/sum(sum(thiscell)));

end

intensitytable=[intensity4Gabbr2_1 intensity4Gabbr2_2 intensity4Tspan7_1 intensity4Tspan7_2 intensity4C1qc intensity4Tms4bx intensity4Rps29 intensity4Ftl1 intensity4mt_Nd2 intensity4mt_Co3];
 writematrix(intensitytable,fullfile(savepath,strcat(thatpath(17:end),'fastIntensitiesallZ.csv')))
