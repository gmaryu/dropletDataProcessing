function labelDifference = compareAgainstToolbox(brightFieldImage, legacyParameters, toolboxParameters, savePath)
    %   compareAgainstToolbox Compares segmentation results between the legacy code and the toolbox
    %
    %   compareAgainstToolbox(birghtFieldImage, legacyParameters, toolboxParameters)
    %
    %   Inputs:
    %       brightFieldImage (MxN) numeric - Bright field image to segment
    %       legacyParameters legacyCode.Parameters - Segmentation parameters to use for the legacy code
    %       toolboxParameters segmentation.Parameters - Segmentation parameters to use for the toolbox
    %      savePath (1x1) string - Path to save the comparison figure to. If empty, no figure is saved. Please include full path, including file name and extension
    %
    %   Outputs:
    %       labelDifference (MxN) numeric - Label difference between the two segmentations. 1 = Toolbox only, 2 = Legacy only
    arguments
        brightFieldImage (:,:) {mustBeNumeric}
        legacyParameters (1,1) legacyCode.Parameters
        toolboxParameters (1,1) segmentation.Parameters
        savePath (1,1) string
    end
    segmentationResult = segmentation.segmentBrightFieldImage(brightFieldImage, toolboxParameters);
    legacyResult = legacyCode.segmentBrightFieldImage(brightFieldImage, legacyParameters);
    % Comparison
    toolboxOverlay = labeloverlay(histeq(brightFieldImage), segmentationResult.labeledImage);
    legacyOverlay = labeloverlay(histeq(brightFieldImage), legacyResult.labeledImage);
    labelDifference = double(logical(segmentationResult.labeledImage)) - double(logical(legacyResult.labeledImage));
    labelDifference(~bwareaopen(labelDifference, 100)) = 0; % Remove small objects
    labelDifference(labelDifference == 1) = 1; % Blue, on toolbox not on legacy
    labelDifference(labelDifference == -1) = 2; % Red, on legacy not on toolbox

    if ~isempty(savePath)
        figure('Visible', 'off');
        subplot(1, 3, 1);
        imshow(toolboxOverlay);
        title('Toolbox');
        subplot(1, 3, 2);
        imshow(legacyOverlay);
        title('Legacy Code');
        subplot(1, 3, 3);
        overlay = labeloverlay(histeq(brightFieldImage), labelDifference, "Colormap", "flag"); % Red, Blue
        imshow(overlay);
        title({'Label Difference', '\color{blue}Toolbox-only, \color{red}Legacy-only'})
        exportgraphics(gcf, savePath, 'Resolution', 500);
        close;
    end
end