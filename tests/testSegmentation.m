classdef testSegmentation < matlab.unittest.TestCase
    %   Parametrized test to compare legacy code segmentation results with current toolbox
    %
    %   Notes:
    %       - It is assumed that images to test are named 'test_image_X.tif' with X an integer
    %       - Ground truth results containing a labeled image and region properties are used for comparison
    %         and are expected to be named 'ground_truth_result_image_X.mat' with X an integer
    properties
        testImagePath = 'images/segmentation_test_images/'; % Location of bright field images to segment
        groundTruthPath = 'tests/segmentation_ground_truth/'; 
        segmentationParameters = segmentation.Parameters(1, 21.6915, 68.5945); % To match legacy thresholds
    end

    properties (TestParameter)
        imageFileNames = num2cell(strcat('test_image_', string(1:14)));
    end

    methods (Test)
        function testAgainstLegacyCode(testCase, imageFileNames)
            % Result with toolbox
            filePath = pwd() + "/" + testCase.testImagePath + imageFileNames + ".tif";
            brightFieldImage = imread(filePath);
            segmentationResult = segmentation.segmentBrightFieldImage(brightFieldImage, ...
                                                                      testCase.segmentationParameters);
            % Legacy result
            fileNameComponents = split(imageFileNames, "_");
            groundTruthFileName = "ground_truth_result_image_" + fileNameComponents(end) + ".mat";
            legacyResult = load(testCase.groundTruthPath + groundTruthFileName).obj;
            % Comparison
            toolboxOverlay = labeloverlay(histeq(brightFieldImage), segmentationResult.labeledImage);
            legacyOverlay = labeloverlay(histeq(brightFieldImage), legacyResult.labeledImage);
            labelDifference = logical(segmentationResult.labeledImage) - logical(legacyResult.labeledImage);
            labelDifference = bwareaopen(labelDifference, 100); % Remove small objects
            differenceRegions = bwlabel(labelDifference);
            numRegions = length(unique(differenceRegions));
            percentDiff = 100 * numRegions / size(legacyResult.regionProperties, 1);

            figure('visible', 'off');
            subplot(1, 3, 1);
            imshow(toolboxOverlay);
            title('Toolbox');
            subplot(1, 3, 2);
            imshow(legacyOverlay);
            title('Legacy Code');
            subplot(1, 3, 3);
            imshow(labelDifference);
            title('Label Difference')
            figSavePath = pwd() + "/debug/test_segmentation/";
            fileName = figSavePath + "overlay_" + fileNameComponents(end) + ".png";
            exportgraphics(gcf, fileName, 'Resolution', 300);
            close;

            disp('Fininshed with ' + imageFileNames + '.tif. Difference percentage: ' + string(percentDiff))
        end
    end
end