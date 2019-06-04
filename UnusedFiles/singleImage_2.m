clc;
clear;
close all;

%% Read single DICOM Image
dInfo = dicominfo('000082.dcm');
dReference = imread('abnormal1.jpg');
% dImage = uint8(dicomread(dInfo));
dImage = dicomread(dInfo);
dImage = imhistmatch(dImage, dReference);
figure, imshow(dImage, []), title('Original Image');

%extract size for planeXY, XZ, YZ from meta data
pixel_spacing = dInfo.PixelSpacing;
per_pixel_area = pixel_spacing(1)*pixel_spacing(2);
% num_white_pixels = nnz(binarized_img);
% total_white_area = num_white_pixels * per_pixel_area;

%% Smoothing - Apply median filter 
I_t = medfilt2(dImage);

%% Smoothing - Gaussian filter
I_t = imgaussfilt(I_t,2);
imageSegmenter(dImage);

%% Sharpen Image
% I_t = imsharpen(I_t, 'Radius', 1.5, 'Amount', 10);
% figure, imshow(I_t), title('Sharpened Image');

%% Adaptive histogram - Not using*
I_t = adapthisteq(I_t);
figure, imshow(I_t), title('Adaptive Histogram Equalization');

%% Image Enhancement - Gabor filter
lambda  = 2;
theta   = 0;
bw      = 3;
psi     = [0 0];
gamma   = 2; 
N       = 8;
img_out = zeros(size(dImage,1), size(dImage,2), N);
for n=1:N
    gb = gabor_fn(bw,gamma,psi(1),lambda,theta)...
        + gabor_fn(bw,gamma,psi(2),lambda,theta);
    img_out(:,:,n) = imfilter(dImage, gb, 'symmetric');
    figure, imshow(img_out(:,:)), title('Gabor');
    theta = theta + pi/N;
end
img_out_disp = sum(abs(uint8(img_out)), 3);
img_out_disp = img_out_disp./max(img_out_disp(:));
figure, imshow(img_out_disp, []), title('Sum of Gabor at all orientations (theta)');

%% Image Segmentation - Erosion and Dilation to get the binarized image's threshold value (Using Otsu's thresholding)
se = strel('disk', 2);
% Use img_out_disp to in  imeroed and imreconstruct to visualize Gabor's
% output --- more branches (lighter shades will be included)
Ie = imerode(I_t, se);
Iobr = imreconstruct(Ie, I_t);
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
BW = imbinarize(Iobrcbr, graythresh(Iobr));
figure, imshow(BW, []), title('Segmented Image using Otsu Thresholding');

%% Marker Controlled Watershed - Not using*
% fgm = imregionalmax(Iobrcbr);
% I2 = labeloverlay(dImage,fgm);
% se2 = strel(ones(5,5));
% fgm2 = imclose(fgm,se2);
% fgm3 = imerode(fgm2,se2);
% fgm4 = bwareaopen(fgm3,20);
% I3 = labeloverlay(dImage,fgm4);
% bw = imbinarize(Iobrcbr);
% D = bwdist(bw);
% DL = watershed(D);
% bgm = DL == 0;
% gmag2 = imimposemin(gmag, bgm | fgm4);
% L = watershed(gmag2);
% labels = imdilate(L==0,ones(3,3)) + 2*bgm + 3*fgm4;
% I4 = labeloverlay(dImage,labels);
% imshow(I4)

%%

% Get rid of stuff touching the border
holes = imclearborder(BW);
holesAccurate = bwareafilt(holes, [50 1000]);% tumour are usually larger than 100
figure, imshow(holesAccurate), title('Clear Border');

% Show tumour boundaries
boundary = bwboundaries(holesAccurate);
figure, imshow(img_out_disp)
hold on
visboundaries(holesAccurate, 'Color', 'b');

labeledImage = bwlabel(holesAccurate, 8); % Label each blob so we can make measurements of it
blobMeasurements = regionprops(labeledImage,dImage,'all'); % Get all the blob properties
stage_one = {};
stage_two = {};
stage_three = {};
stage_four = {};
tumour = [];
for i = 1:length(blobMeasurements)
    thisBlob = blobMeasurements(i);
    blobBoundary = (thisBlob.BoundingBox);
    blobArea = thisBlob.Area * per_pixel_area
    blobPerimeter = thisBlob.PerimeterOld * pixel_spacing
    blobMajorAxis = thisBlob.MajorAxisLength * pixel_spacing
    blobEccentricity = thisBlob.Eccentricity;
    
    % Less Circular when result deviates far from 1.
    blobCircularity = (blobPerimeter.^2) ./ (4*pi*blobArea)
   
    if (blobArea>206 && blobArea<=341) && (blobPerimeter > 54 && blobPerimeter <=77)
        tumour = [tumour blobMeasurements{i}];
        stage_one = {stage_one blobMeasurements{i}};
    end
    
    if (blobArea>341 && blobArea<=491) && (blobPerimeter > 77 && blobPerimeter <=94)
        tumour = [tumour blobMeasurements{i}];
        stage_two = {stage_two blobMeasurements{i}};
    end
    
    if (blobArea>491 && blobArea<=608) && (blobPerimeter > 94 && blobPerimeter <=109)
        tumour = [tumour blobMeasurements{i}];
        stage_three = {stage_three blobMeasurements{i}};
    end
    
    if (blobArea>608) && (blobPerimeter > 109)
        tumour = [tumour blobMeasurements{i}];
        stage_four = {stage_four blobMeasurements{i}};
    end
        
end
 
 

%% Feature Extration - Find possible tumours

% % Get rid of stuff touching the border
% lungsOnly = imclearborder(BW);
% figure, imshow(lungsOnly), title('Clear Border');
% detect = bwareafilt(imclearborder(BW), [50 1000]);% tumour are usually larger than 100
% figure, imshow(detect); title('Possible cancer nodule location');
% [B,~] = bwboundaries(detect,'noholes');
% stats = regionprops(BW,'Area','centroid');
% tumour = 0;
% for j = 1:length(B)
%   metric = 4*pi*stats(j).Area/sum(sqrt(sum(diff(B{j}).^2,2)))^2;
%   %disp(metric);
%     
%   if metric > 3
%       tumour = B{j};
%   end
% end

%% Binarization for black pixels - Not using*
% tmp =ones(512,512);
% black=0;
%
% for i=1:512;
%     for j=1:512;
%         tmp(i,j)=img_out_disp(i,j);
%     end
% end
% 
% for i=1:512;
%     for j=1:512;
%         if tmp(i,j)<=0.12
%                 black=black+1;
%             end
%     end
% end
% if black>=17000
%     ('not a lung cancer patient')
% else
%     ('lung cancer patient')
% end

%% Feature Extration - Plot tumour location
if tumour ~= 0
    figure, imshow(dImage), title('Possible cancer nodule on original image - with Gabor filter')
    hold on;
    x = tumour(:,2);
    y = tumour(:,1);
    plot(x, y, 'color','r','linestyle','-','linewidth', 3)
    disp(metric);
    disp('Possible cancer nodule found - Category: Cancer Patient');
else
    disp('No lung cancer nodule found - Category: Normal Patient');
end