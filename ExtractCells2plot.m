% Loic Binan
%lbinan@broadinstitute.org
%3/26/2021

function [myList]=ExtractCells2plot(myclustergreen,myslice)
load('listsOfcells.mat')
MyTable(:,1)=table2array(ctrkoinfoseparateCTRKO(:,1));
MyTable(:,2)=table2array(ctrkoinfoseparateCTRKO(:,3));
MyTable(:,3)=table2array(ctrkoinfoseparateCTRKO(:,4));
MyTable(:,4)=table2array(ctrkoinfoseparateCTRKO(:,22));
MyTable(:,5)=table2array(ctrkoinfoseparateCTRKO(:,34));
MyTable(:,6)=table2array(ctrkoinfoseparateCTRKO(:,35));
MyTable(:,7)=table2array(ctrkoinfoseparateCTRKO(:,36));


myList=[];
for i=1:size(MyTable,1)
    if MyTable(i,5)==myslice%slice1side2
        if MyTable(i,6)==myclustergreen
            myList=[myList;MyTable(i,1:7)];
        end
    end
end
myMinimum=0;
for ii=1:size(MyTable,1)
    if MyTable(ii,5)==myslice && myMinimum==0
        myMinimum=ii;
    end
end
myList(:,1)=myList(:,1)-myMinimum+1;
