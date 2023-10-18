function pathToFile = getPathToFile(pathToFolder, channel, frame, microscope)
    % TODO: Add description
    arguments
        pathToFolder string % Full path to folder containing data
        channel string % Channel name
        frame {mustBeInteger} % Frame number
        microscope string % Which epifluorescence microscope was used. Must be "new" or "old"
    end
    if strcmp(microscope, 'new')
        channelID = MicroscopeInfo.newEpiChannelID(channel);
    elseif strcmp(microscope, 'old')
        channelID = MicroscopeInfo.oldEpiChannelID(channel);
    else
        error('Invalid microscope type. Must be "new" or "old"')
    end
    formattedFrame = sprintf("%09d", frame);
    pathToFile = [pathToFolder, "/img_", ...
                    formattedFrame, "_", channelID, ...
                    "-", channel, "_000.tif"];
    pathToFile = strjoin(pathToFile, '');
end