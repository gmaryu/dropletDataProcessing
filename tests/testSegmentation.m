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
            labelDifference = double(logical(segmentationResult.labeledImage)) - double(logical(legacyResult.labeledImage));
            labelDifference(~bwareaopen(labelDifference, 100)) = 0; % Remove small objects
            labelDifference(labelDifference == 1) = 1; % Blue, on toolbox not on legacy
            labelDifference(labelDifference == -1) = 2; % Red, on legacy not on toolbox
            differenceRegions = bwlabel(labelDifference);
            numRegions = length(unique(differenceRegions));
            percentDiff = 100 * numRegions / size(legacyResult.regionProperties, 1);

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
            figSavePath = pwd() + "/debug/test_segmentation/";
            fileName = figSavePath + "comparison_" + fileNameComponents(end) + ".png";
            exportgraphics(gcf, fileName, 'Resolution', 500);
            close;

            disp('Fininshed with ' + imageFileNames + '.tif. Difference percentage: ' + string(percentDiff))
        end
    end
end