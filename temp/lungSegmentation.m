clc;
clear;
close all;
%% MATLAB EXAMPLE
% https://www.mathworks.com/help/images/segment-lungs-from-3-d-chest-mri-data.html

load chestVolume
whos
V1 = im2single(V);
whos
% volumeViewer(V1);
XY = V1(:,:,160);
XZ = squeeze(V1(256,:,:));
figure, imshow(XY, [],'Border','tight');
figure, imshow(XZ, [],'Border','tight');
% imageSegmenter(XY);
% imageSegmenter(XZ);

%% Image Segmenter XY and XZ
[BW,maskedImageXY] = lungSegmentationXY_graphcut(XY);
[BW2,maskedImageXZ] = lungSegmentationXZ_graphcut(XZ);

figure, imshow(maskedImageXY,'Border','tight');
figure, imshow(maskedImageXZ,'Border','tight');

% imageSegmenter(maskedImageXY);
% imageSegmenter(maskedImageXZ);

% figure, imshow(maskedImageXYnodules, [],'Border','tight');
% figure, imshow(maskedImageXZnodules, [],'Border','tight');

%% Create a 3D mask

% Create a logical 3-D volume the same size as the input volume and insert mask_XY and mask_XZ at the appropriate spatial locations.
mask = false(size(V1));
mask(:,:, 160) = maskedImageXY;
mask(256, :, :) = mask(256, :, :)|reshape(maskedImageXZ, [1, 512, 318]);

% % Volumetric 3d image
% array3D = squeeze(mask);
% sliceNo = 4;  % Show desired no of Slices
% img_out = array3D(:,:,sliceNo);


%% Active Contour using Chan-Vese

% Using this 3-D seed mask, segment the lungs in the 3-D volume using the active contour method. This operation can take a few minutes. To get a quality segmentation, use histeq to spread voxel values over the available range.
V1 = (histeq(V1));

BW  = activecontour(V1,mask, 100,'Chan-Vese');

segmentedImage = V1.*single(BW);

volumeViewer(segmentedImage);

%% Compute the Volume of the Segmented Lungs
