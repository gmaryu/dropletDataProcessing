function tests = testTracking
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    segmentationParameters = segmentation.Parameters();
    trackingParameters = tracking.Parameters();

    brightFieldImage0 = imread("images/tracking_example_data/BF_0.tif");
    brightFieldImage1 = imread("images/tracking_example_data/BF_1.tif");

    segmentationResult0 = segmentation.segmentBrightFieldImage(brightFieldImage0, segmentationParameters);
    segmentationResult1 = segmentation.segmentBrightFieldImage(brightFieldImage1, segmentationParameters);

    testCase.TestData.trackingParameters = trackingParameters;
    testCase.TestData.segmentationResult0 = segmentationResult0;
    testCase.TestData.segmentationResult1 = segmentationResult1;
end

function testLinkSegmentedObjects(testCase)
    segmentationResult0 = testCase.TestData.segmentationResult0;
    segmentationResult1 = testCase.TestData.segmentationResult1;
    trackingParameters = testCase.TestData.trackingParameters;
    linkResult = tracking.linkSegmentedObjects(segmentationResult0, segmentationResult1, ...
                                               trackingParameters);
    assertEqual(testCase, size(linkResult.assigned), [454, 2])
    assertEqual(testCase, size(linkResult.unassignedFromFirstInput), [12, 1])
    assertEqual(testCase, size(linkResult.unassignedFromSecondInput), [13, 1])
end