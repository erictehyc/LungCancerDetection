clc;
clear;

%% Read image
I = imread('abnormal1.jpg');
figure, imshow(I), title('Original Image');

%% Smoothing - Apply median filter 
I_t = medfilt2(I);

%% Smoothing - Gaussian filter
I_t = imgaussfilt(I_t,2);

%% Adaptive histogram - Not using*
% I_t = adapthisteq(I_t);

%% Image Enhancement - Gabor filter
lambda  = 2;
theta   = 0;
bw      = 3;
psi     = [0 0];
gamma   = 2;
N       = 4;
img_out = zeros(size(I,1), size(I,2), N);
for n=1:N
    gb = gabor_fn(bw,gamma,psi(1),lambda,theta)...
        + gabor_fn(bw,gamma,psi(2),lambda,theta);
    img_out(:,:,n) = imfilter(I, gb, 'symmetric');
    theta = theta + pi/4;
end
img_out_disp = sum(abs(img_out).^2, 3).^0.5;
img_out_disp = img_out_disp./max(img_out_disp(:));

%% Image Segmentation - Erosion and Dilation to get the binarized image (Using Otsu's thresholding)
se = strel('disk', 2);
Ie = imerode(I_t, se);
Iobr = imreconstruct(Ie, I_t);
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
BW = imbinarize(Iobrcbr, graythresh(Iobr));

%% Marker Controlled Watershed - Not using*
% fgm = imregionalmax(Iobrcbr);
% I2 = labeloverlay(I,fgm);
% se2 = strel(ones(5,5));
% fgm2 = imclose(fgm,se2);
% fgm3 = imerode(fgm2,se2);
% fgm4 = bwareaopen(fgm3,20);
% I3 = labeloverlay(I,fgm4);
% bw = imbinarize(Iobrcbr);
% D = bwdist(bw);
% DL = watershed(D);
% bgm = DL == 0;
% gmag2 = imimposemin(gmag, bgm | fgm4);
% L = watershed(gmag2);
% labels = imdilate(L==0,ones(3,3)) + 2*bgm + 3*fgm4;
% I4 = labeloverlay(I,labels);
% imshow(I4)

%% Feature Extration - Find possible tumours
detect = bwareafilt(imclearborder(BW), [100 1000]);% tumour are usually larger than 100
figure, imshow(detect); title('Possible cancer nodule location');
[B,~] = bwboundaries(detect,'noholes');
stats = regionprops(BW,'Area','centroid');
tumour = 0;
for j = 1:length(B)
  metric = 4*pi*stats(j).Area/sum(sqrt(sum(diff(B{j}).^2,2)))^2;
  %disp(metric);
    
  if metric > 3
      tumour = B{j};
  end
end

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
    figure, imshow(I), title('Possible cancer nodule on original image - with Gabor filter')
    hold on;
    x = tumour(:,2);
    y = tumour(:,1);
    plot(x, y, 'color','r','linestyle','-','linewidth', 3)
    disp(metric);
    disp('Possible cancer nodule found - Category: Cancer Patient');
else
    disp('No lung cancer nodule found - Category: Normal Patient');
end