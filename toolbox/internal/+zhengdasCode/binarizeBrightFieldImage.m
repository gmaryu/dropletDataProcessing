function binaryMask = binarizeBrightFieldImage(brightFieldImage, segmentationParameters)
    % TODO: Add description
    arguments
        brightFieldImage (:,:) {mustBeNumeric}
        segmentationParameters zhengdasCode.Parameters
    end

    largeFilteredImage = imgaussfilt(brightFieldImage, segmentationParameters.gaussianFilterSizeLarge);
    smallFilteredImage = imgaussfilt(brightFieldImage, segmentationParameters.gaussianFilterSizeSmall);
    gimg2 = largeFilteredImage - smallFilteredImage;
    gimg3 = gimg2;
    gimg3(brightFieldImage < 0) = 0;

    [gDxx,gDxy,gDyy] = zhengdasCode.hessian2D(gimg2, segmentationParameters.hessianSigma);
    % Correct for scale
    gDxx = (segmentationParameters.hessianSigma^2)*gDxx;
    gDxy = (segmentationParameters.hessianSigma^2)*gDxy;
    gDyy = (segmentationParameters.hessianSigma^2)*gDyy;

    % Calculate (abs sorted) eigenvalues and vectors
    [gLambda2, gLambda1, ~, ~] = zhengdasCode.eig2image(gDxx,gDxy,gDyy);
    gS2 = gLambda1.^2;% + Lambda2.^2;
    gI2 = (ones(size(gimg2)) - exp(-gS2/segmentationParameters.c)) .* (gLambda1 < 0) .* gimg3;

    tmpthresh = graythresh(gI2) * segmentationParameters.thresholdFactor;
    binaryMask = gI2 > tmpthresh;

    ngimgl = bwlabel(binaryMask);
    gS = regionprops(ngimgl,'Area','Eccentricity');
    % pick-up of junk labels
    condition1 = [gS.Area] < (50 * segmentationParameters.resAmp ^ 2);
    condition2 = [gS.Area] < 1000 * segmentationParameters.resAmp ^ 2;
    condition3 = [gS.Eccentricity] < segmentationParameters.maxEccentricity;
    delidx = condition1|((condition2)& condition3); 
    % replace labels to zero
    for iij=1:length(delidx)
        if delidx(iij)
            ngimgl(ngimgl==iij)=0;
        end
    end
    ngimg = logical(ngimgl); 
    binaryMask = gimg3 & ngimg;
end