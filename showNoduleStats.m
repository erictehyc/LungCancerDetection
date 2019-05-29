%% Show nodule characteristics on the command Window and label each nodule with numbers on image
function showNoduleStats(noduleStage_stats)
    for k = 1 : length(noduleStage_stats)           % Loop through all blobs.

        blobArea = noduleStage_stats(k).ActualArea;		% Get area.
        blobPerimeter = noduleStage_stats(k).ActualPerimeter;		% Get perimeter.
        blobCentroid = noduleStage_stats(k).Centroid;		% Get centroid one at a time
        blobDiameter = noduleStage_stats(k).ActualMajorAxisLength;		% Get Diameter (ActualMajorAxisLength)
        meanGL = noduleStage_stats(k).MeanIntensity; % Mean Pixel Intensity of Blob
        
    %     x = ['Mean Pixel Internsity', 'Actual Area', 'Actual Perimeter', 'Centroid Coordinates'];
    %     y = [meanGL, blobArea, blobPerimeter, blobCentroid];

        fprintf('Nodule #%2d\n', k);
        fprintf('Diameter of nodule/tumor = %5.1f mm\n', blobDiameter);
        fprintf('Mean Intensity of Pixels = %8.1f\n', meanGL);
        fprintf('Blob Area (mm^2) = %15.1f\n', blobArea);
        fprintf('Blob Perimeter (mm) = %12.1f\n', blobPerimeter);
        fprintf('Centroid Coordinates (x,y) = (%.1f, %.1f)\n\n', blobCentroid(1), blobCentroid(2));
        
        textFontSize = 20;	% Used to control size of "blob number" labels put atop the image.
        labelShiftX = 2;	% Used to align the labels in the centers of the nodules.

        % Put the "blob number" labels on the "boundaries" grayscale image.
        text(blobCentroid(1) + labelShiftX, blobCentroid(2), num2str(k), 'FontSize', textFontSize, 'Color', 'r');
    end
    fprintf('\n');
end