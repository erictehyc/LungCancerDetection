clc;
clear;
close all;
%% Eric open folder for 3d volume - not used

% fileFolder = fullfile(pwd, 'LIDC-IDRI-0001','01-01-2000-30178','3000566-03192');
% files = dir(fullfile(fileFolder, '*.dcm'));%specify data file diectory
% fileNames = {files.name};


%% Read single DICOM Image
dInfo = dicominfo(fullfile('Patient1','000034.dcm'));
%dReference = imread('abnormal1.jpg');
% dImage = uint8(dicomread(dInfo));
dImage = dicomread(dInfo);
%img_in = imhistmatch(dImage, dReference);

%%
img_in = dImage;
figure, imshow(img_in, []), title('Original Image');

%extract size for planeXY, XZ, YZ from meta data
voxel_size = [dInfo.PixelSpacing; dInfo.SliceThickness];

%% Smoothing - Apply median filter 
% img_in = medfilt2(dImage);

%% Smoothing - Gaussian filter
% img_in = imgaussfilt(img_in,2);
% figure, imshow(img_in, []), title('Gaussian filtering');

%% Smoothing - Anisotropic diffusion filtering of images - preserves the sharpness of edges better than Gaussian blurring.
% img_in_2 = imdiffusefilt(dImage);
% figure, imshow(img_in_2, []), title('Anisotropic diffusion filtering');


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

%% Active Contour using masking to get Lung Volume
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

% Closing and filling up holes to obtain nodules close to pleural wall
seg = imclose(seg, strel('disk', 20));
seg = imfill(seg, 'holes');

figure, imshow(seg); title('Global Region-Based Segmentation - RegionActiveContour Mask')

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
figure, imshow(tmp), title('Lung volume applied with RegionActiveContour Mask')

%% Eric to extract nodules using Watershed Segmentation
% https://www.mathworks.com/company/newsletters/articles/the-watershed-transform-strategies-for-image-segmentation.html
% This is marker controlled watershed using masking to get nodules

I_eq = adapthisteq(tmp);
% figure, imhist(tmp), title('Histogram of lung volume image');
% figure, imhist(I_eq), title('Histogram of lung volume image after equalization');

tmpBW = imbinarize(tmp, graythresh(I_eq)); % Built-in Otsu Thresholding to get nodules - blobs will clump together
figure, imshow(tmpBW), title('Otsu BW Segmented lung nodules mask')

watershed_nodule = watershedTransform(tmpBW); % Watershed to separate out blobs (obtained from MathWorks)
watershed_nodule = imbinarize(watershed_nodule);
figure, imshow(watershed_nodule), title('Watershed BW Segmented lung nodules mask');

maskedNodule = img_in;
maskedNodule(~watershed_nodule) = 0;
figure,imshow(maskedNodule, []), title('Nodules using WatershedSegmentation');


%% Visualize the nodules on original image
% Show tumour boundaries
boundary = bwboundaries(watershed_nodule);
figure, imshow(dImage, []), title('Possible Nodules location')
hold on
visboundaries(boundary, 'Color', 'b');

%% Tumor parameters to reduce false positives - hard coded (does not stage tumor based on TNM)

% holesAccurate = bwareafilt(watershed_nodule, [50 1000]); % malignant tumour are usually larger than 50
% boundary = bwboundaries(holesAccurate);
% figure, imshow(dImage, []), title('Possible Nodules location (area considered)');
% hold on
% visboundaries(boundary, 'Color', 'r');

%% Tumor parameters to reduce false positives

holesAccurate = bwareafilt(watershed_nodule, [50 1000]); % malignant tumour are usually larger than 50
labeledImage = bwlabel(holesAccurate); % Label each blob so we can make measurements of it
blobMeasurements = regionprops(labeledImage,dImage,'all'); % Get all the blob properties

%% Tumor parameters to reduce false positives - hard coded - more parameters involved

% Binarized image of nodules using watershed segmentation
BW = watershed_nodule;

% Find connected components (blobs)
cc = bwconncomp(BW); 

% Find statistics/properties of each blob
stats = regionprops(cc, 'all'); 

%% Create new column (dervied) in regionprops (stats) structure array 

% DICOM header - extract pixel size for planeXY, XZ, YZ from DICOM meta data
pixel_spacing = dInfo.PixelSpacing;
per_pixel_area = pixel_spacing(1)*pixel_spacing(2);

% Area based on per pixel area
actual_area = num2cell([stats.Area]*[per_pixel_area]);
[stats.ActualArea] = actual_area{:};

% Actual perimeter based on pixel spacing
actual_perimeter = num2cell([stats.Perimeter]*[pixel_spacing(1)]);
[stats.ActualPerimeter] = actual_perimeter{:};

% Actual major axis based on pixel spacing
actual_diameter = num2cell([stats.MajorAxisLength]*[pixel_spacing(1)]);
[stats.ActualMajorAxisLength] = actual_diameter{:};

% Mean Gray Level Intensity of nodule
for k = 1 : length(stats)           % Loop through all blobs.
        thisBlobsPixels = stats(k).PixelIdxList;        % Get list of pixels in current blob.
        MeanIntensity = mean(dImage(thisBlobsPixels));  % Find mean intensity (in original image!)
        stats(k).MeanIntensity = MeanIntensity;
end
    
%% Find desired Parameters
% Get index in regionprops stucture array that satisfy property according
% to TNM 8th edition

% Stage t1a
if (length(stats) ~= 0)
idx = find([stats.ActualMajorAxisLength] > 3 & [stats.ActualMajorAxisLength] <= 10 & [stats.Eccentricity] < 0.8 & [stats.MeanIntensity] < 800);

t1a = ismember(labelmatrix(cc), idx);  
% figure, imshow(t1a), title('Stage t1a mask - Discard due to high false positive');
t1a_stats = stats(idx);

% Stage t1b
idx = find([stats.ActualMajorAxisLength] > 10 & [stats.ActualMajorAxisLength] <= 20 & [stats.Eccentricity] < 0.8 & [stats.MeanIntensity] < 800);

t1b = ismember(labelmatrix(cc), idx);  
% figure, imshow(t1b), title('Stage t1b mask');
t1b_stats = stats(idx);

% Stage t1c
idx = find([stats.ActualMajorAxisLength] > 20 & [stats.ActualMajorAxisLength] <= 30 & [stats.Eccentricity] < 0.8 & [stats.MeanIntensity] < 800);

t1c = ismember(labelmatrix(cc), idx);  
% figure, imshow(t1c), title('Stage t1c mask');
t1c_stats = stats(idx);

% Stage 2
idx = find([stats.ActualMajorAxisLength] > 30 & [stats.ActualMajorAxisLength] <= 50 & [stats.Eccentricity] < 0.8 & [stats.MeanIntensity] < 800);

t2 = ismember(labelmatrix(cc), idx);  
% figure, imshow(t2), title('Stage t2 mask');
t2_stats = stats(idx);


% Stage 3
idx = find([stats.ActualMajorAxisLength] > 50 & [stats.Eccentricity] < 0.8 & [stats.MeanIntensity] < 800);

t3 = ismember(labelmatrix(cc), idx);  
% figure, imshow(t3), title('Stage t3 mask');
t3_stats = stats(idx);
end

%% Display Result on Command Window

% Discard stage t1a due to high false positive rate
% If there are no tumors > stage t1a, we consider as normal
if length(stats) == 0
    ('normal lung')

else
    if (length(t1b_stats) + length(t1c_stats) + length(t2_stats) == 0)
    ('normal lung')
    
    
    % Visualize output
    fprintf('\nEarly Detection saves lives! Possible nodules:\n\n');
    if length(t1a_stats) ~= 0
        t1a_holes = bwlabel(t1a); 
        boundary = bwboundaries(t1a_holes);
        figure, imshow(dImage, []), title('Probable T1a Tumors');
        hold on
        visboundaries(boundary, 'Color', 'r');

        maskedImageT1a = dImage;
        maskedImageT1a(~t1a_holes) = 0;
        
        disp('T1a\n');
        figure, imshow(maskedImageT1a, []), title(' Probable T1a Tumors');
        showNoduleStats(t1a_stats);
    end
        
    else
        ('tumor detected')

        % Visualize output of t1a nodules
        if length(t1a_stats) ~= 0
            t1a_holes = bwlabel(t1a); 
            boundary = bwboundaries(t1a_holes);
            figure, imshow(dImage, []), title('Probable T1a Tumors');
            hold on
            visboundaries(boundary, 'Color', 'r');
    
            maskedImageT1a = dImage;
            maskedImageT1a(~t1a_holes) = 0;
            
            disp('T1a');
            figure, imshow(maskedImageT1a, []), title(' Probable T1a Tumors');
            showNoduleStats(t1a_stats);
        end

        % Visualize output of t1b nodules
        if length(t1b_stats) ~= 0
            t1b_holes = bwlabel(t1b); 
            boundary = bwboundaries(t1b_holes);
            figure, imshow(dImage, []), title('T1b Tumors');
            hold on
            visboundaries(boundary, 'Color', 'r');

            maskedImageT1b = dImage;
            maskedImageT1b(~t1b_holes) = 0;
            
            disp('T1b');
            figure, imshow(maskedImageT1b, []), title('T1b Tumors');
            showNoduleStats(t1b_stats);
        end

        % Visualize output
        if length(t1c_stats) ~= 0
            t1c_holes = bwlabel(t1c); 
            boundary = bwboundaries(t1c_holes);
            imshow(dImage, []), title('T1c Tumors');
            hold on
            visboundaries(boundary, 'Color', 'g');

            maskedImageT1c = dImage;
            maskedImageT1c(~t1c_holes) = 0;
            figure, imshow(maskedImageT1c, []), title('T1c Tumors');

            disp('T1c');
            showNoduleStats(t1c_stats);
        end

        % Visualize output
        if length(t2_stats) ~= 0
            t2_holes = bwlabel(t2); 
            boundary = bwboundaries(t2_holes);
            figure, imshow(dImage, []), title('T2 Tumors');
            hold on
            visboundaries(boundary, 'Color', 'c');

            maskedImageT2 = dImage;
            maskedImageT2(~t2_holes) = 0;
            figure, imshow(maskedImageT2, []), title('T2 Tumors');
        
            disp('T2');
            showNoduleStats(t2_stats)
        end
        
          % Visualize output
        if length(t3_stats) ~= 0
            t3_holes = bwlabel(t3); 
            boundary = bwboundaries(t3_holes);
            figure, imshow(dImage, []), title('T3 Tumors');
            hold on
            visboundaries(boundary, 'Color', 'y');

            maskedImageT3 = dImage;
            maskedImageT3(~t3_holes) = 0;
            figure, imshow(maskedImageT3, []), title('T3 Tumors');

            disp('T3');
            showNoduleStats(t3_stats)
        end
    end
end
