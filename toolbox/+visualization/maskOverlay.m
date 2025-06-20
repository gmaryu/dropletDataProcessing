function maskOverlay(database, position, dropletID, numFrames, FPChannel, saveDir, alphaValue)

for k = 1:length(database)
    posID = database{k}.posId;
    if position == posID
        croppedDir = database{k}.croppedImages;
        matDir = database{k}.maskMatFiles;
        mkdir(fullfile(saveDir, sprintf('Pos%d', position)));
        break;
    end
end

% Loop over each item folder
% Construct the folder name (e.g. 'droplet_001')
folderName = fullfile(croppedDir, sprintf('droplet_%03d', dropletID));
% Construct GIF output file name
gifName    = fullfile(saveDir, sprintf('Pos%d', position), sprintf('overlay_%03d.gif', dropletID));

% --- 1) Attempt to load binaryMask from the .mat file ---
maskFilePath = fullfile(matDir,sprintf('nuclear_%03d.mat', dropletID));
hasMask      = false;  % default to no mask
if exist(maskFilePath, 'file') == 2
    loadedData = load(maskFilePath);
    % Check if the variable nuclearMask is present
    if isfield(loadedData, 'nuclearMask')
        binaryMask = loadedData.nuclearMask;  % [W x H x T]
        hasMask    = true;
    end
end

% If there's no valid mask, skip creating a GIF
if ~hasMask
    warning('No mask found for item %03d. Skipping GIF creation.', dropletID);
    return;
end

% --- 2) Gather TIFF frames and sort by numeric index ---
fn = sprintf('Pos%d_%s_*.tif', position, FPChannel);
imgFiles = dir(fullfile(folderName, fn));
if isempty(imgFiles)
    warning('No images found in %s. Skipping...', folderName);
    return;
end

% Parse numeric part of each filename, e.g. 'Pos7_CFP_15.tif' -> 15
fileNumbers = zeros(size(imgFiles));
for f = 1:numel(imgFiles)
    thisName  = imgFiles(f).name;
    % Extract all digits
    numStr    = regexp(thisName, '\d+', 'match');
    % Convert the first match to a number
    fileNumbers(f) = str2double(numStr{1});
end

% Sort by these numeric values
[~, sortIdx] = sort(fileNumbers);
imgFiles = imgFiles(sortIdx);

% Ensure we only process up to numFrames frames
if length(imgFiles) < numFrames
    numFrames = length(imgFiles);
end
imgFiles = imgFiles(1:numFrames);

% --- 3) Build frames for the GIF ---
for tIdx = 1:numFrames

    % Read the image
    imgName = fullfile(folderName, imgFiles(tIdx).name);
    grayRaw = imread(imgName);

    % Convert to double, then rescale intensities to [0..255]
    grayDouble = double(grayRaw);
    maxGrayVal = max(grayDouble(:));
    if maxGrayVal == 0
        % Avoid divide-by-zero if the image is entirely 0
        grayDouble = zeros(size(grayDouble));
    else
        grayDouble = (grayDouble ./ maxGrayVal) * 255;
    end

    % If it's multi-channel, reduce to grayscale
    if ndims(grayDouble) == 3
        grayDouble = rgb2gray(grayDouble);
    end

    % Create an RGB "canvas" in double
    overlayRGB = cat(3, grayDouble, grayDouble, grayDouble);

    % We must ensure we don’t exceed the mask’s size in the T dimension
    if tIdx <= size(binaryMask, 3)
        maskSlice = binaryMask(:,:,tIdx);

        % Blend color [0, 0, 255] (blue) in masked regions
        % newValue = (1-alpha)*oldValue + alpha*colorValue
        R = overlayRGB(:,:,1);
        G = overlayRGB(:,:,2);
        B = overlayRGB(:,:,3);

        idxMask = (maskSlice == 1);
        R(idxMask) = (1 - alphaValue)*R(idxMask) + alphaValue*0;
        G(idxMask) = (1 - alphaValue)*G(idxMask) + alphaValue*0;
        B(idxMask) = (1 - alphaValue)*B(idxMask) + alphaValue*255;

        overlayRGB = cat(3, R, G, B);
    end

    % Convert blended overlay to uint8 for GIF
    overlayRGB = uint8(overlayRGB);

    % Convert the overlaid image to an indexed image for GIF
    [indImg, cm] = rgb2ind(overlayRGB, 256);

    % Write or append to the GIF file
    if tIdx == 1
        imwrite(indImg, cm, gifName, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(indImg, cm, gifName, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end

fprintf('Finished GIF for folder %s -> %s\n', folderName, gifName);
end