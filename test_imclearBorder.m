
clc;
clear;
close all;

%% imclearBorder



dInfo = dicominfo('000003.dcm');
dReference = imread('abnormal1.jpg');
% dImage = uint8(dicomread(dInfo));
dImage = dicomread(dInfo);
I = imhistmatch(dImage, dReference);
binaryImage = imbinarize(I);
figure, imshow(binaryImage, []), title('Original Image');

binaryImage = zeros(size(I)+1,'uint8');    %// initialize a padded matrix
% binaryImage = binaryImage(:,:,1:3);                  %// (we didn't need to pad the third-dimension)
binaryImage(2:end,2:end,:) = I;            %// assign I to the lower-right of J



% Draw lines around the border of your binary image:
binaryImage(:,1) = true;
binaryImage(:,end) = true;
binaryImage(1,:) = true;
binaryImage(end,:) = true;

binaryImage = imfill(binaryImage, 'holes');      
figure, imshow(binaryImage),title('fill holes');

% Now do an imclearborder to get rid of all that and leave just the lung edges
binaryImage = imclearborder(binaryImage);
figure, imshow(binaryImage);

%Now do a fill again to get the lungs filled.
lungsMask = imfill(binaryImage, 'holes');
figure, imshow(lungsMask)

test_imclearBorder2();


%%%%%%%%%%%%%%%%%%%%%%%%%
% cropped = I(90:400,40:460);
% figure, imshow(cropped, []), title('Cropped Image');
% 
% %%Threshold to isolate lung tissue
% thresholded = cropped < 86;
% %%Remove artifacts attached to border
% clearThresh = imclearborder(thresholded);
% %%Remove objects less than 40 pixels in size
% imgBigTissue = bwareaopen(clearThresh,40);
% figure, imshow(imgBigTissue), title('Lung tissue');

%%%%%%%%%%%%%%%%%%%%%%%%%
% J = zeros(size(I)+1,'uint8');    %// initialize a padded matrix
% J = J(:,:,1:3);                  %// (we didn't need to pad the third-dimension)
% J(2:end,2:end,:) = I;            %// assign I to the lower-right of J
% 
% C = imclearborder(J,4);          %// now run imclearborder, using connectedness option 4
% C = C(2:end,2:end,:);            %// remove the padding