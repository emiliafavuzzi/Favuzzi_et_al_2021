function fastCounts(thatpath,cellmaskTmem119)
% Loic Binan
%lbinan@broadinstitute.org
%3/26/2021
savepath=fullfile('/broad/hptmp/lbinan/microglia/',thatpath);
mypath=fullfile(savepath,'/merfish_mosaics');
SE=strel('disk',6);
%cellmaskallTmem_Gabbr1_2=imread(fullfile(savepath,'analysis','cellmaskallTmem_Gabbr1_2.png'));
%cellmaskTmem119=imread(fullfile(savepath,'analysis','cellmaskTmem119.png'));

%cellMask=im2bw(imdilate(cellmaskallTmem_Gabbr1_2,SE));
cellMask=im2bw(cellmaskTmem119);

%imwrite(cellMask,fullfile(savepath,'analysis','positiveCellsMAskMoreTmem.png'))

%%
%  SE=strel('disk',2);
% white=ones(size(cellMask));
%  negativeCells=bwareaopen(cellmaskTmem119.*(white-imdilate(cellMask,SE)),650);
% imwrite(negativeCells,fullfile(savepath,'analysis','negativeCellsMAskMoreTmem.png'))

 %imshow(negativeCells);
%%
% cellMask=im2bw(negativeCells);
%imshow(cellMask);
stats=regionprops(cellMask,'image','area','PixelIdxList');
thiscell=zeros(size(cellMask));
countsforTmem119=zeros([1+size(stats,1),1]);
countsforTmem119(1)="Tmem119";
countsforGabbr1=zeros([1+size(stats,1),1]);
countsforGabbr1(1)="Gabbr1";
countsforGabbr2=zeros([1+size(stats,1),1]);
countsforGabbr2(1)="Gabbr2";
countsforCd164=zeros([1+size(stats,1),1]);
countsforCd164(1)="Cd164";
countsforClec4a3=zeros([1+size(stats,1),1]);
countsforClec4a3(1)="Clec4a3";

countsforEcscr=zeros([1+size(stats,1),1]);
countsforEcscr(1)="Ecscr";

countsforFcrls=zeros([1+size(stats,1),1]);
countsforFcrls(1)="Fcrls";

countsforGpr34=zeros([1+size(stats,1),1]);
countsforGpr34(1)="Gpr34";
countsforLaptm4a=zeros([1+size(stats,1),1]);
countsforLaptm4a(1)="Laptm4a";
countsforLaptm5=zeros([1+size(stats,1),1]);
countsforLaptm5(1)="Laptm5";
countsforP2ry12=zeros([1+size(stats,1),1]);
countsforP2ry12(1)="P2ry12";
countsforP2ry13=zeros([1+size(stats,1),1]);
countsforP2ry13(1)="P2ry13";
countsforPpib=zeros([1+size(stats,1),1]);
countsforPpib(1)="Ppib";
countsforSelenok=zeros([1+size(stats,1),1]);
countsforSelenok(1)="Selenok";
countsforSelenop=zeros([1+size(stats,1),1]);
countsforSelenop(1)="Selenop";
countsforTmem14c=zeros([1+size(stats,1),1]);
countsforTmem14c(1)="Tmem14c";
countsforTrem2=zeros([1+size(stats,1),1]);
countsforTrem2(1)="Trem2";
countsforTspan3=zeros([1+size(stats,1),1]);
countsforTspan3(1)="Tspan3";
countsforTspan4=zeros([1+size(stats,1),1]);
countsforTspan4(1)="Tspan4";
countsforTspan7=zeros([1+size(stats,1),1]);
countsforTspan7(1)="Tspan7";
thissize=zeros([1+size(stats,1),1]);
thissize(1)="cell size";
maskCd164=imread(fullfile(mypath,'mask_Cd164.tif'));
maskTmem119=imread(fullfile(mypath,'mask_Tmem119.tif'));
maskGabbr1=imread(fullfile(mypath,'mask_Gabbr1.tif'));
maskGabbr2=imread(fullfile(mypath,'mask_Gabbr2.tif'));

maskFcrls=imread(fullfile(mypath,'mask_Fcrls.tif'));
maskEcscr=imread(fullfile(mypath,'mask_Ecscr.tif'));
maskClec4a3=imread(fullfile(mypath,'mask_Clec4a3.tif'));
maskFcrls=imread(fullfile(mypath,'mask_Fcrls.tif'));
maskGpr34=imread(fullfile(mypath,'mask_Gpr34.tif'));
maskLaptm4a=imread(fullfile(mypath,'mask_Laptm4a.tif'));
maskLaptm5=imread(fullfile(mypath,'mask_Laptm5.tif'));
maskP2ry12=imread(fullfile(mypath,'mask_P2ry12.tif'));
maskP2ry13=imread(fullfile(mypath,'mask_P2ry13.tif'));
maskPpib=imread(fullfile(mypath,'mask_Ppib.tif'));
        maskSelenok=imread(fullfile(mypath,'mask_Selenok.tif'));
    maskSelenop=imread(fullfile(mypath,'mask_Selenop.tif'));    
    maskTmem14c=imread(fullfile(mypath,'mask_Tmem14c.tif'));    
    maskTrem2=imread(fullfile(mypath,'mask_Trem2.tif'));
    maskTspan3=imread(fullfile(mypath,'mask_Tspan3.tif'));    
    maskTspan4=imread(fullfile(mypath,'mask_Tspan4.tif'));    
    maskTspan7=imread(fullfile(mypath,'mask_Tspan7.tif'));
parfor (i=1:size(stats,1),12)
    thiscell=zeros(size(cellMask));
    thissize(i+1)=stats(i).Area;
    thiscell(stats(i).PixelIdxList)=1;
%     if stats(i).Area<2500
%         thiscell=imdilate(thiscell,strel('disk',4));
%     end

    thisImage=thiscell.*double(imbinarize(maskTmem119));
    thesecounts=sum(sum(thisImage));
    countsforTmem119(i+1)=(thesecounts);
  
        thisImage=thiscell.*double(imbinarize(maskGabbr1));
    thesecounts=sum(sum(thisImage));
    countsforGabbr1(i+1)=(thesecounts);
            thisImage=thiscell.*double(imbinarize(maskGabbr2));
    thesecounts=sum(sum(thisImage));
    countsforGabbr2(i+1)=(thesecounts);
thisImage=thiscell.*double(imbinarize(maskCd164));
    thesecounts=sum(sum(thisImage));
    countsforCd164(i+1)=(thesecounts);
    thisImage=thiscell.*double(imbinarize(maskFcrls));
    thesecounts=sum(sum(thisImage));
    countsforFcrls(i+1)=(thesecounts);
    thisImage=thiscell.*double(imbinarize(maskEcscr));
    thesecounts=sum(sum(thisImage));
    countsforEcscr(i+1)=(thesecounts);
    thisImage=thiscell.*double(imbinarize(maskClec4a3));
    thesecounts=sum(sum(thisImage));
    countsforClec4a3(i+1)=(thesecounts);
    thisImage=thiscell.*double(imbinarize(maskFcrls));
    thesecounts=sum(sum(thisImage));
    countsforFcrls(i+1)=(thesecounts);
    thisImage=thiscell.*double(imbinarize(maskGpr34));
    thesecounts=sum(sum(thisImage));
    countsforGpr34(i+1)=(thesecounts); 
    thisImage=thiscell.*double(imbinarize(maskLaptm4a));
    thesecounts=sum(sum(thisImage));
    countsforLaptm4a(i+1)=(thesecounts); 
    thisImage=thiscell.*double(imbinarize(maskLaptm5));
    thesecounts=sum(sum(thisImage));
    countsforLaptm5(i+1)=(thesecounts);     
    thisImage=thiscell.*double(imbinarize(maskP2ry12));
    thesecounts=sum(sum(thisImage));
    countsforP2ry12(i+1)=(thesecounts);   
    thisImage=thiscell.*double(imbinarize(maskP2ry13));
    thesecounts=sum(sum(thisImage));
    countsforP2ry13(i+1)=(thesecounts);  
     thisImage=thiscell.*double(imbinarize(maskPpib));
    thesecounts=sum(sum(thisImage));
    countsforPpib(i+1)=(thesecounts);   
    thisImage=thiscell.*double(imbinarize(maskSelenok));
    thesecounts=sum(sum(thisImage));
    countsforSelenok(i+1)=(thesecounts); 
    thisImage=thiscell.*double(imbinarize(maskSelenop));
    thesecounts=sum(sum(thisImage));
    countsforSelenop(i+1)=(thesecounts);
    thisImage=thiscell.*double(imbinarize(maskTmem14c));
    thesecounts=sum(sum(thisImage));
    countsforTmem14c(i+1)=(thesecounts);
        thisImage=thiscell.*double(imbinarize(maskTrem2));
    thesecounts=sum(sum(thisImage));
    countsforTrem2(i+1)=(thesecounts);
            thisImage=thiscell.*double(imbinarize(maskTspan3));
    thesecounts=sum(sum(thisImage));
    countsforTspan3(i+1)=(thesecounts);
    thisImage=thiscell.*double(imbinarize(maskTspan4));
    thesecounts=sum(sum(thisImage));
    countsforTspan4(i+1)=(thesecounts);
        thisImage=thiscell.*double(imbinarize(maskTspan7));
    thesecounts=sum(sum(thisImage));
    countsforTspan7(i+1)=(thesecounts);
end
    
counttable=[countsforTmem119 countsforGabbr1 countsforGabbr2 countsforCd164 countsforClec4a3 countsforEcscr countsforFcrls countsforGpr34 countsforLaptm4a countsforLaptm5 countsforP2ry12 countsforP2ry13 countsforPpib countsforSelenok countsforSelenop countsforTmem14c countsforTrem2 countsforTspan3 countsforTspan4 countsforTspan7 thissize];
% counttable=num2str(counttable);

writematrix(counttable,fullfile(savepath,strcat(thatpath(17:end),'fastCountsNoSmallCell.csv')));
% 


clear maskCd164;
    clear maskFcrls;
    clear maskEcscr;
    clear maskClec4a3; 
    clear maskFcrls;
    clear maskGpr34; 
    clear maskLaptm4a;
    clear maskLaptm5;
    clear maskP2ry12;
    clear maskP2ry13;
     clear maskPpib;
    clear maskSelenok;
    clear maskSelenop;
    clear maskTmem14c;
    clear Trem2;
    clear maskTspan3;
    clear maskTspan4;
        clear maskTspan7;


