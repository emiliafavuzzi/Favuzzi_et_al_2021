function [thisImage]=groupSMfish(inputpath, outputpath)

for i=13:22
mypath=fullfile(strcat(inputpath,'/mosaic_bit',num2str(i),'_0','.tif'));
savepath=fullfile(strcat(outputpath,'/smFISH_mosaics/mosaic_bit',num2str(i),'_compiled','.tif'));

thisImage=myHighPass(imread(mypath));

% if i==13
% %       myzeros=zeros(size(thisImage));
% %     for count=1:10
% %         myzeros=zeros(size(thisImage));
% %     mosaic(:,:,count)=myzeros;
% %     disp(count)
% %     end
% end
for a=1:6
    mypath=fullfile(strcat(inputpath,'/mosaic_bit',num2str(i),'_',num2str(a),'.tif'));
    addthis=myHighPass(imread(mypath));
    thisImage=thisImage+addthis;
end
% mosaic(:,:,i-12)=thisImage;
imwrite(thisImage,fullfile(savepath));
% save(fullfile(outputpath,'smfishmosaics.mat'));
end
