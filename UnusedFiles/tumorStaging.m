%% Run singleImage3.m first
BW = watershed_nodule;

cc = bwconncomp(BW); 

stats = regionprops(cc, 'all'); 
%% Create new column (dervied) in regionprops structure array 

% Area based on per pixel area
actual_area = num2cell([stats.Area]*[per_pixel_area]);
[stats.ActualArea] = actual_area{:};

% Actual perimeter based on pixel spacing
actual_perimeter = num2cell([stats.Perimeter]*[pixel_spacing(1)]);
[stats.ActualPerimeter] = actual_perimeter{:};

% Actual major axis based on pixel spacing
actual_diameter = num2cell([stats.MajorAxisLength]*[pixel_spacing(1)]);
[stats.ActualMajorAxisLength] = actual_diameter{:};

% Actual major axis based on pixel spacing
actual_diameter = num2cell([stats.MajorAxisLength]*[pixel_spacing(1)]);
[stats.ActualMajorAxisLength] = actual_diameter{:};

%% Find desired Parameters
idx = find([stats.ActualMajorAxisLength] > 5 & [stats.ActualMajorAxisLength] < 10);
figure, imshow(idx), title('idx');

BW2 = ismember(labelmatrix(cc), idx);  
figure, imshow(BW2), title('Stage t1a mask');

% Get all desired tumor info from index retrieved from find() on stats
t1a = stats(idx);
