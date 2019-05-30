clc;
clear;
[V,s,d] = dicomreadVolume(fullfile('LIDC-IDRI-0001','01-01-2000-30178', '3000566-03192'));

fileFolder = fullfile('LIDC-IDRI-0001','01-01-2000-30178', '3000566-03192');
% imds = imageDatastore('LIDC-IDRI-0001\01-01-2000-30178\3000566-03192','FileExtensions',{'.dcm'});
files = dir(fullfile(fileFolder, '*.dcm'));%specify data file diectory
fileNames = {files.name};

%examine file header meta datafrom dicom stack
info = dicominfo(fullfile(fileFolder, fileNames{1}));

%extract size for planeXY, XZ, YZ from meta data
voxel_size = [info.PixelSpacing; info.SliceThickness];

% hWaitBar=  waitbar(0,'Reading DICOM files');
% for i=1:length(fileNames)
for i=1:15
    fname = fullfile(fileFolder, fileNames{i});
    dInfo = dicominfo(fname);
    I = dicomread(dInfo);
    figure, imshow(I, []);    
%     waitbar((length(fileNames):-1+1)/length(fileNames))
end

