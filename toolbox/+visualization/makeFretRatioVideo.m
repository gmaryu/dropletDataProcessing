function makeFretRatioVideo(dataDirectory, fretChannelName, outputDirectory, outputName, videoSize, autoDetectFretLimits, fretRatioLimits, videoFPS, resizeMethod, verbose)
    %   Makes a fret ratio video from a directory containing CFP and FRET images.
    %
    %   makeFretRatioVideo(dataDirectory, fretChannelName, outputDirectory, outputName);
    % 
    %   Required Inputs:
    %       dataDirectory - Directory containing CFP and FRET images
    %       fretChannelName - Name of the fret channel (e.g. 'FRET' or 'Custom')
    %       outputDirectory - Directory to save the output video
    %       outputName - Name of the output video without file extension
    %   Optional Inputs:
    %       videoSize - Size of the output video in pixels (default [576 576])
    %       autoDetectFretLimits - Automatically detect the fret ratio limits (default true)
    %       fretRatioLimits - Fret ratio limits to use when autoDetectFretLimits is false (default [NaN NaN])
    %       videoFPS - Frames per second of the output video (default 20)
    %       resizeMethod - Method used to resize the images (default 'bilinear')
    %       verbose - Display progress messages (default true)
    arguments
        % Required
        dataDirectory (1,1) string
        fretChannelName (1,1) string
        outputDirectory (1,1) string
        outputName (1,1) string
        % Optional
        videoSize (1,2) {mustBeInteger} = [576 576]
        autoDetectFretLimits (1,1) logical = true
        fretRatioLimits (1,2) {mustBeNumeric} = [-Inf Inf]
        videoFPS (1,1) {mustBeInteger} = 20
        resizeMethod (1,1) string = "bilinear"
        verbose (1,1) logical = true
    end
    if ~autoDetectFretLimits && any(isnan(fretRatioLimits))
        error("fretRatioLimits must be specified when autoDetectFretLimits is false")
    end
    % Locate CFP and FRET images
    CFP_files = dir(dataDirectory + "/*CFP*.tif");
    FRET_files = dir(dataDirectory + "/*" + fretChannelName + "*.tif");
    % Initialize global limits and video info arrays
    if autoDetectFretLimits
        globalLowRatio = 10;
        globalHighRatio = 0;
    else
        globalLowRatio = fretRatioLimits(1);
        globalHighRatio = fretRatioLimits(2);
    videoImages = cell(1, length(CFP_files));
    backgroundMasks = cell(1, length(CFP_files));
    end

    for i = 1:length(CFP_files)
        if (mod(i, 100) == 0) && verbose
            disp("Processing image " + i + " of " + length(CFP_files))
        end
        % Read CFP and FRET images
        CFP = imread(fullfile(CFP_files(i).folder, CFP_files(i).name));
        FRET = imread(fullfile(FRET_files(i).folder, FRET_files(i).name));
        % Resize images
        resizedCFP = imresize(CFP, videoSize, resizeMethod);
        resizedFRET = imresize(FRET, videoSize, resizeMethod);
        % Create droplet mask
        dropletMaskCFP = imquantize(resizedCFP, multithresh(resizedCFP, 1), [0 1]);
        dropletMaskFRET = imquantize(resizedFRET, multithresh(resizedFRET, 1), [0 1]);
        dropletMask = dropletMaskCFP | dropletMaskFRET;
        % Calculate fret ratio
        fretRatio = double(resizedFRET) ./ double(resizedCFP);
        maskedRatio = fretRatio .* dropletMask;
        % Store images
        videoImages{i} = maskedRatio;
        backgroundMasks{i} = ~ dropletMask;
        if autoDetectFretLimits
            % Update global limits
            globalLowRatio = min(globalLowRatio, quantile(nonzeros(maskedRatio), 0.01));
            globalHighRatio = max(globalHighRatio, quantile(nonzeros(maskedRatio), 0.99));
        end
    end

    % Create video
    cmap = jet(256);
    video = VideoWriter(fullfile(outputDirectory, outputName + ".avi"));
    video.FrameRate = videoFPS;
    open(video);
    for i = 1:length(videoImages)
        % Normalize image
        normalizedImage = (videoImages{i} - globalLowRatio) ./ (globalHighRatio - globalLowRatio);
        normalizedImage(normalizedImage < 0) = 0;
        normalizedImage(normalizedImage > 1) = 1;
        % Create video frame
        indFrame = gray2ind(normalizedImage, 256);
        rgbFrame = ind2rgb(indFrame, cmap);
        % Add background mask
        background = repmat(backgroundMasks{i}, [1, 1, 3]);
        rgbFrame(background) = 0;
        % Write frame
        writeVideo(video, rgbFrame);
    end
    close(video);
end