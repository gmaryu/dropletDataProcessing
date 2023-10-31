function debugSegmentation(brightFieldImage, saveFilePath, segmentationParameters)
    %   debugSegmentation Create a 6 plot figure with every segmentation crucial step
    %
    %   debugSegmentation(brightFieldImage, saveFilePath, segmentationParameters)
    %
    %   Inputs:
    %       brightFieldImage (MxN) numeric - Bright field image to debug
    %       savePath (1,1) string - Full path including filename for saving the resulting figure. Please include the .png extension. If left empty, the figure will not be saved but only displayed.
    %       segmentationParameters (1,1) segmentation.Parameters - Parameters to use for segmentation
    arguments
        brightFieldImage (:,:) {mustBeNumeric}
        saveFilePath (1,1) string
        segmentationParameters (1,1) segmentation.Parameters
    end

    preProcessed = segmentation.preProcessBrightFieldImage(brightFieldImage, segmentationParameters);
    enhancedBordersImage = segmentation.enhanceDropletBorders(preProcessed, segmentationParameters);
    watershedMask = segmentation.getMaskForWatershed(preProcessed, enhancedBordersImage, segmentationParameters);
    labeledImage = segmentation.labelObjects(watershedMask, segmentationParameters);
    [filteredLabeledImage, ~] = segmentation.filterLabels(labeledImage, segmentationParameters);

    if saveFilePath ~= ""
        figure('Visible', 'off');
    end
    subplot(2, 3, 1);
    imshow(histeq(brightFieldImage));
    title('Raw Image');
    subplot(2, 3, 2);
    imshow(preProcessed);
    title('Pre-processed');
    subplot(2, 3, 3);
    imshow(enhancedBordersImage);
    title('Border enhancement')
    subplot(2, 3, 4);
    imshow(watershedMask);
    title('Watershed mask')
    subplot(2, 3, 5);
    allLabelOverlay = labeloverlay(preProcessed, labeledImage);
    imshow(allLabelOverlay);
    title('All labels')
    subplot(2, 3, 6);
    filteredLabelsOverlay = labeloverlay(preProcessed, filteredLabeledImage);
    imshow(filteredLabelsOverlay);
    title('Filtered labels')
    if saveFilePath ~= ""
        exportgraphics(gcf, saveFilePath, 'Resolution', 500);
        close;
    else
        % Display current figure
        shg;
    end
end