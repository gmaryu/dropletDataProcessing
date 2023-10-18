function tests = testSegmentation
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    brightFieldImage = imread("images/brightFieldExample.tif");
    segmentationParameters = segmentation.Parameters();
    testCase.TestData.brightFieldImage = brightFieldImage;
    testCase.TestData.segmentationParameters = segmentationParameters;
end

function testPipeline(testCase)
    % TODO: Implement actual tests. This code is just a placeholder
    % See https://www.mathworks.com/help/matlab/matlab_prog/write-function-based-unit-tests.html
    % See https://www.mathworks.com/help/matlab/matlab_prog/write-test-using-setup-and-teardown-functions.html
    brightFieldImage = testCase.TestData.brightFieldImage;
    segmentationParameters = testCase.TestData.segmentationParameters;
    segmentationResult = segmentation.segmentBrightFieldImage(brightFieldImage, segmentationParameters);
    labeledImage = segmentationResult.labeledImage;
    verifyEqual(testCase, size(labeledImage), size(brightFieldImage));
end