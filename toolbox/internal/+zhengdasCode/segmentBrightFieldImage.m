function segmentationResult = segmentBrightFieldImage(brightFieldImage, segmentationParameters)
    % TODO: Add description
    arguments (Input)
        brightFieldImage (:,:) {mustBeNumeric}
        segmentationParameters zhengdasCode.Parameters
    end

    arguments (Output)
        segmentationResult zhengdasCode.Result
    end
    % Pre-process image to make it suitable for thresholding and segmentation
    preProcessedImage = zhengdasCode.preProcessBrightFieldImage(brightFieldImage, segmentationParameters);
    % Convert image into binary mask
    binaryMask = zhengdasCode.binarizeBrightFieldImage(preProcessedImage, segmentationParameters);
    % Apply watershed algorithm (see https://www.mathworks.com/help/images/ref/watershed.html)
    seeds = imextendedmax(bwdist(binaryMask), 2);
    gimg = double(brightFieldImage);
    gimg4 = max(gimg(:)) - gimg;
    gss1 = watershed(imimposemin(gimg4, seeds));
    gss2 = watershed(imimposemin(gimg4, gss1==0|seeds));
    gss2(gss2==1) = 0;
    gss3 = logical(gss2);
    gss3 = imopen(gss3,strel('disk',round(2 * segmentationParameters.resAmp)));
    gss2 = double(gss2) .* double(gss3);

    S1 = regionprops(gss2, gimg, 'Centroid', 'Area', 'Perimeter', ...
                     'PixelIdxList', 'PixelValues', 'Eccentricity', 'EulerNumber');

    % collection and overwrite of junk label
    areaT=[1e3*segmentationParameters.resAmp^2, 1e4*segmentationParameters.resAmp^2]; % min and max threshold of area size
    delid = [S1.Area]>areaT(2)|[S1.Area]<areaT(1); % too big or too small
    delid = delid|[S1.Eccentricity]>segmentationParameters.maxEccentricity; % roundness
    delid = delid|[S1.Perimeter].^2./[S1.Area]>(1+segmentationParameters.prthresh)*4*pi; % mismutch of perimeter and area ?
    delid = delid|[S1.EulerNumber]>1|[S1.EulerNumber]<1;
    S1(delid==1) = [];
    
    % apply delid to watershed result (cleaned segmentation result)
    ss5=gss2;
    for ij=1:max(ss5(:))
        if delid(ij)
            ss5(ss5==ij)=0;
        end
    end
    
    ss5l = unique(ss5);
    uqss5 = num2cell(ss5l(2:end));
    [S1.Label] = uqss5{:};

    % Create segmentation result
    filteredLabeledImage = ss5;
    filteredRegionProperties = S1;
    segmentationResult = zhengdasCode.Result(filteredLabeledImage, filteredRegionProperties);
end