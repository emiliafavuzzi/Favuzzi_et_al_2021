%script used to generate a mask of the region of interest, avoiding counts from cell outside of it. save the variable space as matlab.mat after running it
% Loic Binan
%lbinan@broadinstitute.org
%3/26/2021
figure; imshow(5*mosaic_DAPI_3, []);
roi = images.roi.AssistedFreehand;
draw(roi);mask = createMask(roi);imshow(mask)
