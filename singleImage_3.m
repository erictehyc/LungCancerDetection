clc;
clear;
close all;
%% Eric

fileFolder = fullfile(pwd, 'LIDC-IDRI-0001','01-01-2000-30178','3000566-03192');
files = dir(fullfile(fileFolder, '*.dcm'));%specify data file diectory
fileNames = {files.name};


%% Read single DICOM Image
dInfo = dicominfo('000096.dcm');
%dReference = imread('abnormal1.jpg');
% dImage = uint8(dicomread(dInfo));
dImage = dicomread(dInfo);
%img_in = imhistmatch(dImage, dReference);
img_in = dImage;
figure, imshow(img_in, []), title('Original Image');

%extract size for planeXY, XZ, YZ from meta data
voxel_size = [dInfo.PixelSpacing; dInfo.SliceThickness];

%% Smoothing - Apply median filter 
img_in = medfilt2(dImage);

%% Smoothing - Gaussian filter
img_in = imgaussfilt(img_in,2);
figure, imshow(img_in, []), title('Gaussian filtering');

%% Smoothing - Anisotropic diffusion filtering of images - preserves the sharpness of edges better than Gaussian blurring.
img_in_2 = imdiffusefilt(dImage);
figure, imshow(img_in_2, []), title('Anisotropic diffusion filtering');


%% Adaptive histogram - Not using*
% img_in = adapthisteq(img_in);

%% Image Enhancement - Gabor filter
lambda  = 25;
theta   = 0;
bw      = 3;
psi     = [0 0];
gamma   = 2;
N       = 4;
img_out = zeros(size(img_in,1), size(img_in,2), N);
for n=1:N
    gb = gabor_fn(bw,gamma,psi(1),lambda,theta)...
        + gabor_fn(bw,gamma,psi(2),lambda,theta);
    img_out(:,:,n) = imfilter(img_in, gb, 'symmetric');
    theta = theta + pi/4;
end

img_out_disp = sum(abs(img_out).^2, 3).^0.5;
img_out_disp = img_out_disp./max(img_out_disp(:));
figure, imshow(img_out_disp), title('gabor output, normalized');

%% Eric - get lung volume

% % Built-in Otsu Global Thresholding to find threshold value
% T = graythresh(uint16(img_in));         % Threshold value
% BW_OG = imbinarize(uint16(img_in),T);      % Threshold mask
% % figure,imshow(BW_OG, []), title('Mask using Built-in Otsu Threshold');
% 
% % Invert mask
% BW = imcomplement(BW_OG);
% 
% % Fill holes
% BW = imfill(BW, 'holes');
% % figure,imshow(BW, []), title('Mask using Built-in Otsu Threshold');
% 
% 
% % Clear borders
% BW = imclearborder(BW);
% % figure,imshow(BW, []), title('Mask using Built-in Otsu Threshold after refining');
% 
% maskedImage = img_in;
% maskedImage(~BW) = 0;
% figure,imshow(maskedImage, []), title('Lung Volume using Built-in Otsu Threshold');

%% Eric - get nodules
% holes = imclearborder(BW_OG);
% holesAccurate = bwareafilt(holes, [0, 1000]);
% labeledImage = bwlabel(holesAccurate, 8);
% blobMeasurements = regionprops(labeledImage, img_in,'all');
% figure, imshow(labeledImage), title('Possible Tumors Mask');
% 
% maskedImage = img_in;
% maskedImage(~labeledImage) = 0;
% nnz(maskedImage)
% figure, imshow(maskedImage, []), title('Masked Image showing Possible Tumors');

%% This is marker controlled watershed using masking
Ir = img_out_disp;
se = strel('disk', 20);
Ie = imerode(Ir, se);
Iobr = imreconstruct(Ie, Ir);
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
bw = imbinarize(Iobrcbr, graythresh(Iobr));
figure, imshow(bw), title('Binarized')
Ia=bw;
m = zeros(size(Ia,1),size(Ia,2));
m(200:320,95:180) = 1;
m(186:321,348:410) = 1;
seg = region_seg(Ia, m, 800); % Run segmentation with 800 iteration
figure, imshow(seg); title('Global Region-Based Segmentation - MarkerControlledWatershed Mask')

%% Binarization for image classification (masked Image)
tmp=ones(512,512);
black=0;
for i=1:512
    for j=1:512
        if seg(i,j)==1
            tmp(i,j)=img_out_disp(i,j);
            if tmp(i,j)<=0.12
                black=black+1;
            end
        else
            tmp(i,j)=1;
        end
    end
end
tmp = imclearborder(tmp);
figure, imshow(tmp), title('Lung volume applied with MarkerControlledWatershed Mask')


%%
tmpBW = imbinarize(tmp); % Built-in Otsu Thresholding to get nodules - blobs will clump together
figure, imshow(tmpBW), title('BW Segmented lung nodules mask')

%% This is marker controlled watershed using masking

I_eq = adapthisteq(tmp);

tmpBW = imbinarize(tmp, graythresh(I_eq)); % Built-in Otsu Thresholding to get nodules - blobs will clump together
figure, imshow(tmpBW), title('Otsu BW Segmented lung nodules mask')
bw2 = imfill(tmpBW,'holes');
bw3 = imopen(bw2, ones(5,5));
bw4 = bwareaopen(bw3, 40);
bw4_perim = bwperim(bw4);
overlay1 = imoverlay(I_eq, bw4_perim, [.3 1 .3]);
figure, imshow(overlay1), title('Overlay');

mask_em = imextendedmax(I_eq, 30);
figure, imshow(mask_em)

mask_em = imclose(mask_em, ones(5,5));
mask_em = imfill(mask_em, 'holes');
mask_em = bwareaopen(mask_em, 40);
overlay2 = imoverlay(I_eq, bw4_perim | mask_em, [.3 1 .3]);
figure, imshow(overlay2), title('Overlay2');

I_eq_c = imcomplement(I_eq);
I_mod = imimposemin(I_eq_c, ~bw4 | mask_em);

L = watershed(I_mod);
figure, imshow(label2rgb(L)); title('Global Region-Based Segmentation to get nodules - MarkerControlledWatershed Mask')

%% Eric to extract nodules

% Show tumour boundaries
boundary = bwboundaries(tmpBW);
figure, imshow(dImage, []), title('Possible Nodules location')
hold on
visboundaries(boundary, 'Color', 'b');

%% Tumor parameters to reduce false positives

holesAccurate = bwareafilt(tmpBW, [50 1000]);% malignant tumour are usually larger than 100
figure, imshow(holesAccurate), title('Possible Nodules location (area considered)');
hold on
visboundaries(boundary, 'Color', 'r');

labeledImage = bwlabel(holesAccurate); % Label each blob so we can make measurements of it
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

%%
if black>=17179
    ('normal lung')
else
    ('lung cancer')
end