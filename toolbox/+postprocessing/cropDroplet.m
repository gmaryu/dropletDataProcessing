function cropDroplet(trackmateOutput, pathDatabase, ignoreFile, labelPath, pathOutput, imageNameFormats, outputPrefix, outputChannelNames, dx)
% cropDroplet2  Crop individual droplet images based on tracking data.
%
%   cropDroplet2(trackmateOutput, pathDatabase, ignoreFile, pathOutput, ...
%                imageNameFormats, outputPrefix, outputChannelNames, dx)
%
% This function reads droplet tracking data from a TrackMate CSV file, optionally
% filters out droplets specified in an ignore file, computes the maximum droplet radius
% (plus a safety margin), and crops individual droplets from multi-channel microscopy images.
%
% Inputs:
%   trackmateOutput    - (1x1 string) Path to the TrackMate output CSV file containing droplet
%                        tracking data.
%   pathDatabase       - (1x1 string) Directory containing the full frame images. The last folder 
%                        in this path must denote the position (e.g., "Pos0").
%   ignoreFile         - (1x1 string) Path to the "force_ignore.csv" file used to filter out droplets.
%   labelPath          - (1x1 string) Path to the "label" directory used to filter out droplets. 
%   pathOutput         - (1x1 string) Directory where the cropped droplet images will be saved.
%   imageNameFormats   - (1xN string array) Format strings for the image filenames for each channel.
%                        Example: ["img_%09d_4-BF_000.tif", "img_%09d_1-DAPI_000.tif", ...].
%   outputPrefix       - (1x1 string) A prefix for the output file names (e.g., "Pos0").
%   outputChannelNames - (1xN string array) Short names for each channel (e.g., ["BF", "DAPI", ...]).
%   dx                 - (1x1 double) Physical length per pixel (e.g., in micrometers).
%
% Behavior:
%   For each frame in the tracking data, the function:
%     1. Loads the corresponding multi-channel image.
%     2. Computes a crop region for each droplet based on its maximum radius (with a margin).
%     3. Saves each cropped droplet image into a subfolder named after the droplet ID.
%
% Notes:
%   - Tracking data must include the fields: TRACK_ID, FRAME, POSITION_X, POSITION_Y, and RADIUS.
%   - Droplet IDs are assumed to start at 0.
%
% Example:
%   cropDroplet("exports/20250328_Nocodazole/Pos0_segmented_spots.csv", ...
%                "raw/20250328_Nocodazole/Pos0", ...
%                "exports/20250328_Nocodazole/force_ignore.csv", ...
%                "exports/20250328_Nocodazole/cropped_pos0", ...
%                ["img_%09d_4-BF_000.tif", "img_%09d_1-DAPI_000.tif", "img_%09d_5-CFP_000.tif", ...
%                 "img_%09d_6-YFP_000.tif", "img_%09d_8-Custom_000.tif"], ...
%                "Pos0", ["BF", "DAPI", "CFP", "YFP", "FRET"], 2.649);

    arguments
        trackmateOutput (1,1) string
        pathDatabase    (1,1) string
        ignoreFile      (1,1) string
        labelPath       (1,1) string
        pathOutput      (1,1) string
        imageNameFormats (1,:) string
        outputPrefix    (1,1) string
        outputChannelNames (1,:) string
        dx              (1,1) double {mustBePositive}
    end

    %% Parameters
    radiusMargin = 0.01;  % Additional margin added to droplet radius (fraction)

    %% Create output folder if it does not exist
    if ~exist(pathOutput, "dir")
        mkdir(pathOutput);
    end

    %% Load tracking data
    opts = detectImportOptions(trackmateOutput);
    opts.DataLines = 5;  % Adjust if necessary based on file format
    opts.VariableNamesLine = 1;
    trackData = readtable(trackmateOutput, opts);
    trackData = trackData(~isnan(trackData.TRACK_ID), :);

    %% Filter out droplets using force_ignore.csv if available
    if exist(ignoreFile, 'file')
        ignoreTable = readtable(ignoreFile);
        % Extract current position from the last folder of pathDatabase (e.g., "Pos0")
        [~, currentPosID, ~] = fileparts(pathDatabase);
        posNum = sscanf(currentPosID, "Pos%d");
        ignoreRows = ignoreTable.PosID == posNum;
        ignoreTablePos = ignoreTable(ignoreRows, :);
        % Assuming the ignore table has a column named 'DropID' listing droplet IDs.
        ignoreDroplets = ignoreTablePos{:, 'DropID'};
        trackData = trackData(~ismember(trackData.TRACK_ID, ignoreDroplets), :);
    end

    %% Create lookup table for maximum droplet radius (over all frames)
    uniqueDropletIds = unique(trackData.TRACK_ID);
    % Preallocate lookup table; using id+1 for zero-indexed droplet IDs.
    maxRadiusLookupTable = zeros(max(uniqueDropletIds) + 1, 1);
    for dropletIdx = 1:length(uniqueDropletIds)
        id = uniqueDropletIds(dropletIdx);
        % Calculate the maximum droplet radius over frames, add margin, and convert to pixels.
        maxRadius = max(trackData(trackData.TRACK_ID == id, :).("RADIUS"));
        maxRadiusLookupTable(id + 1) = round(maxRadius * (1 + radiusMargin) / dx);
    end

    %% Process each frame
    %maxFrame = max(trackData.FRAME);
    maxFrame = 2;
    tic;
    for frameIdx = 1:(maxFrame + 1)
        frame = frameIdx - 1;
        elapsedTime = toc;
        % Print progress information.
        lineLength = fprintf("Processing frame %d / %d, elapsed time: %.1f sec (estimated total: %.1f sec)", ...
            frame, maxFrame, elapsedTime, elapsedTime * (maxFrame + 1) / frameIdx);

        % Get tracking data for the current frame and sort by TRACK_ID.
        frameData = sortrows(trackData(trackData.FRAME == frame, :), "TRACK_ID");

        % Extract droplet IDs and positions (converted to pixel coordinates).
        ids = frameData.TRACK_ID;
        x = round(frameData.POSITION_X / dx);
        y = round(frameData.POSITION_Y / dx);

        % Loop over each channel.
        for channelIdx = 1:length(imageNameFormats)
            % Read the full frame image for the given channel.
            imagePath = fullfile(pathDatabase, sprintf(imageNameFormats(channelIdx), frame));
            fullImage = imread(imagePath)';
            [maxX, maxY] = size(fullImage);

            % Process each droplet in the current frame.
            for dropletIdx = 1:length(ids)
                id = ids(dropletIdx);
                % Determine crop radius for the droplet.
                r = maxRadiusLookupTable(id + 1);
                % Define crop boundaries ensuring they are within the image dimensions.
                xStart = max(x(dropletIdx) - r, 1);
                xEnd   = min(x(dropletIdx) + r, maxX);
                yStart = max(y(dropletIdx) - r, 1);
                yEnd   = min(y(dropletIdx) + r, maxY);
                croppedImage = fullImage(xStart:xEnd, yStart:yEnd);

                % Create subfolder for the droplet if it doesn't exist.
                dropletDir = fullfile(pathOutput, sprintf("droplet_%03d", id));
                if ~exist(dropletDir, "dir")
                    mkdir(dropletDir);
                end

                % Save the cropped image with an appropriate file name.
                outputFileName = sprintf("%s_%s_%03d.tif", outputPrefix, outputChannelNames(channelIdx), frame);
                imwrite(croppedImage, fullfile(dropletDir, outputFileName));
            end

            
        end

        %% label image crop
        % Read the full label image for the given channel.
        labelNameFormat = 'raw_label_img_BF_%03d.tif';
        labelimagePath = fullfile(labelPath, sprintf(labelNameFormat, frame));
        if ~exist(labelimagePath,'file')
            fprintf('Label Not Found. Time or labelNameFromat might be different.')
            return
        else
            fullLabel = imread(labelimagePath)';
            [maxX, maxY] = size(fullLabel);
        end
        
        % Process each droplet in the current frame.
        for dropletIdx = 1:length(ids)
            id = ids(dropletIdx);
            labelid = fullLabel(x(dropletIdx),y(dropletIdx));
            % Determine crop radius for the droplet.
            r = maxRadiusLookupTable(id + 1);
            % Define crop boundaries ensuring they are within the image dimensions.
            xStart = max(x(dropletIdx) - r, 1);
            xEnd   = min(x(dropletIdx) + r, maxX);
            yStart = max(y(dropletIdx) - r, 1);
            yEnd   = min(y(dropletIdx) + r, maxY);
            croppedLabel = fullLabel(xStart:xEnd, yStart:yEnd);
            croppedLabel(croppedLabel ~= labelid) = 0;

            % Save the cropped image with an appropriate file name.
            dropletDir = fullfile(pathOutput, sprintf("droplet_%03d", id));
            outputFileName = sprintf("%s_label_%03d.tif", outputPrefix, frame);
            imwrite(croppedLabel, fullfile(dropletDir, outputFileName));
        end        

        % Erase the previous progress message.
        fprintf(repmat('\b', 1, lineLength));
    end
    fprintf("Done.\n");
end
