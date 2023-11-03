function calculateGroundTruthForTests(testImagesPath, saveResultsPath, totalImagesToTest)
    %   calculateGroundTruthForTests Creates labeled images and regions to be used for testing the segmentation algorithm in this toolbox
    %
    %   calculateGroundTruthForTests(testImagesPath, saveResultsPath, totalImagesToTest)
    %
    %   Inputs:
    %       testImagesPath string - Full path to test images. Located by default at /images/segmentation_test_images
    %       saveResultsPath string - Path where results will be saved. Located by default at /tests/segmentation_ground_truth
    %       totalImagesToTest integer - Number of images to include for testing. Maximum is 14
    arguments
        testImagesPath (1,1) string
        saveResultsPath (1,1) string
        totalImagesToTest (1,1) {mustBeInteger}
    end
    if ~endsWith(testImagesPath, '/')
        testImagesPath = testImagesPath + '/';
    end
    if ~endsWith(saveResultsPath, '/')
        saveResultsPath = saveResultsPath + '/';
    end
    segmentationParameters = legacyCode.Parameters();
    for idx=1:totalImagesToTest
        filePath = testImagesPath + 'test_image_' + string(idx) + '.tif';
        brightFieldImage = imread(filePath);
        segmentationResult = legacyCode.segmentBrightFieldImage(brightFieldImage, segmentationParameters);
        saveName = 'ground_truth_result_image_' + string(idx);
        segmentationResult.save(saveResultsPath, saveName);
        disp("Finished with test image " + string(idx) + " of " + string (totalImagesToTest) + ".");
    end
end