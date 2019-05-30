%%
clc;
clear;
close all;
%%
% function L =  watershedSample(img)

I = uint16(imread('pears.png'));
I = tmp;
I = adapthisteq(I);
figure, imshow(I)
text(732,501,'Image courtesy of Corel(R)',...
     'FontSize',7,'HorizontalAlignment','right')
 
 
%% Step 2: Use the Gradient Magnitude as the Segmentation Function
 gmag = imgradient(I);
 figure, imshow(gmag,[]), title('Gradient Magnitude');

L = watershed(gmag);
Lrgb = label2rgb(L);
figure, imshow(Lrgb), title('Watershed Transform of Gradient Magnitude');

%% Step 3: Mark the Foreground Objects

se = strel('disk',3);
Io = imopen(I,se);
figure, imshow(Io), title('Opening');

%Next compute the opening-by-reconstruction using imerode and imreconstruct.
Ie = imerode(I,se);
Iobr = imreconstruct(Ie,I);
figure, imshow(Iobr), title('Opening-by-Reconstruction');

% Following the opening with a closing can remove the dark spots and stem marks. Compare a regular morphological closing with a closing-by-reconstruction. First try imclose:
Ioc = imclose(Io,se);
figure, imshow(Ioc), title('Opening-Closing');

%Now use imdilate followed by imreconstruct. Notice you must complement the image inputs and output of imreconstruct.
Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
figure, imshow(Iobrcbr), title('Opening-Closing by Reconstruction');

% Calculate the regional maxima of Iobrcbr to obtain good foreground markers.
fgm = imregionalmax(Iobrcbr);
figure, imshow(fgm), title('Regional Maxima of Opening-Closing by Reconstruction');

% To help interpret the result, superimpose the foreground marker image on the original image.
I2 = labeloverlay(I,fgm);
figure, imshow(I2), title('Regional Maxima Superimposed on Original Image');

%some of the mostly-occluded and shadowed objects are not marked.
% clean the edges of the marker blobs and then shrink them a bit. You can do this by a closing followed by an erosion.
se2 = strel(ones(5,5));
fgm2 = imclose(fgm,se2);
fgm3 = imerode(fgm2,se2);

% stray isolated pixels that must be removed. You can do this using bwareaopen, which removes all blobs that have fewer than a certain number of pixel
fgm4 = bwareaopen(fgm3,20);
I3 = labeloverlay(I,fgm4);
figure, imshow(I3), title('Modified Regional Maxima Superimposed on Original Image');

%% Step 4: Compute Background Markers
bw = imbinarize(Iobrcbr);
figure, imshow(bw), title('Thresholded Opening-Closing by Reconstruction');

% We don't want the background markers to be too close to the edges of the objects we are trying to segment. We'll "thin" the background by computing the "skeleton by influence zones", or SKIZ, of the foreground of bw
% omputing the watershed transform of the distance transform of bw
D = bwdist(bw);
DL = watershed(D);
bgm = DL == 0;
figure, imshow(bgm), title('Watershed Ridge Lines)');

%% Step 5: Compute the Watershed Transform of the Segmentation Function.

% imimposemin can be used to modify an image so that it has regional minima only in certain desired locations
gmag2 = imimposemin(gmag, bgm | fgm4);
figure, imshow(gmag2), title('gmag2');

% y we are ready to compute the watershed-based segmentation
L = watershed(gmag2);
figure, imshow(L), title('L');


%% Step 6: Visualize the Result

labels = L==0 + 2*bgm + 3*fgm4;
I4 = labeloverlay(I,labels);
figure, imshow(I4), title('Markers and Object Boundaries Superimposed on Original Image');

Lrgb = label2rgb(L,'jet','w','shuffle');
figure, imshow(Lrgb), title('Colored Watershed Label Matrix');

% use transparency to superimpose this pseudo-color label matrix on top of the original intensity image.
figure, imshow(I)
hold on
himage = imshow(Lrgb);
himage.AlphaData = 0.3;
title('Colored Labels Superimposed Transparently on Original Image')
% end