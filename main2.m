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
img_in = uint16(dImage);
figure, imshow(img_in, []), title('Original Image');

%extract size for planeXY, XZ, YZ from meta data
voxel_size = [dInfo.PixelSpacing; dInfo.SliceThickness];

%% Smoothing - Apply median filter 
img_in = medfilt2(img_in);

%% Smoothing - Gaussian filter
img_in = imgaussfilt(img_in,2);

%% Adaptive histogram - Not using*
%img_in = adapthisteq(img_in);

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

% Built-in Otsu Global Thresholding to find threshold value
T = graythresh(img_in);         % Threshold value
BW = imbinarize(img_in,T);      % Threshold mask
figure,imshow(BW, []), title('Mask using Built-in Otsu Threshold');

% Invert mask
BW = imcomplement(BW);

% Fill holes
BW = imfill(BW, 'holes');

% Clear borders
BW = imclearborder(BW);

figure,imshow(BW, []), title('Mask using Built-in Otsu Threshold After Refining');

maskedImage = img_in;
maskedImage(~BW) = 0;
figure,imshow(maskedImage, []), title('Lung Volume using Built-in Otsu Threshold');

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
figure, imshow(seg); title('Global Region-Based Segmentation')

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
figure, imshow(tmp), title('Segmented lung')

%% Eric
[B,L] = bwboundaries(tmp,'noholes');
figure, imshow(label2rgb(L, @jet, [.5 .5 .5]))
hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end

%%
if black>=17179
    ('normal lung')
else
    ('lung cancer')
end
