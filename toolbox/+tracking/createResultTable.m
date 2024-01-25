function resultTable = createResultTable(segmentationResultArray, trackingResult, channelsOfInterest, timeDelta, calculateRatio)
    %   createResultTable Combines the segmentation and tracking results into a single table for downstream analysis
    %
    %   resultTable = tracking.createResultTable(segmentationResultArray, trackingResult, channelsOfInterest, timeDelta, calculateRatio)
    %
    %   Inputs:
    %       segmentationResultArray (1,:) segmentation.Result - Array of segmentation results
    %       trackingResult (1,1) tracking.Result - Tracking result
    %       channelsOfInterest (1,:) string - Array of channel names to include in the result table. Use the channel name as it appears on the image files
    %       timeDelta (1,1) double - Time between frames in minutes. Used to calculate the time column of the result table
    %       calculateRatio (1,1) logical - If true, the ratio between FRET and CFP channels will be calculated and added to the result table
    %
    %   Output:
    %       resultTable (1,:) table - Table containing the combined segmentation and tracking results. Each entry of the table corresponds to a single droplet
    arguments
        segmentationResultArray (1,:) segmentation.Result
        trackingResult (1,1) tracking.Result
        channelsOfInterest (1,:) string
        timeDelta (1,1) double
        calculateRatio (1,1) logical
    end
    tableColumns = ["Frame", "Label", "x", "y", "Area", "Time"];
    tableColumns = [tableColumns, channelsOfInterest];
    if calculateRatio
        tableColumns = [tableColumns, "Ratio"];
    end

    variableUnits = ["", "", "pixels", "pixels", "pixels^2", "minutes"];
    variableUnits = [variableUnits, repmat("a.u.", 1, length(channelsOfInterest))];
    if calculateRatio
        variableUnits = [variableUnits, "dimensionless"];
    end

    totalDroplets = length(trackingResult.dropletArray);
    totalFrames = length(segmentationResultArray);
    resultTable = cell(1, totalDroplets);

    for dropletIndex = 1:totalDroplets
        droplet = trackingResult.dropletArray(dropletIndex);
        frames = 1:totalFrames;
        time = (frames - 1) * timeDelta;
        frames = frames(~isnan(droplet.segmentationLabels));
        time = time(~isnan(droplet.segmentationLabels));
        labels = droplet.segmentationLabels;
        labels = labels(~isnan(labels));
        x = zeros(size(frames));
        y = zeros(size(frames));
        area = zeros(size(frames));
        channelIntensities = zeros(length(frames), length(channelsOfInterest));
        if calculateRatio
            ratio = zeros(size(frames));
        end

        for frameIndex = 1:length(frames)
            frame = frames(frameIndex);
            label = labels(frameIndex);
            segmentationResult = segmentationResultArray(frame);
            regionProperties = segmentationResult.regionProperties;
            dropletArea = regionProperties(label).Area;
            dropletX = regionProperties(label).Centroid(1); 
            dropletY = regionProperties(label).Centroid(2);
            x(frameIndex) = dropletX;
            y(frameIndex) = dropletY;
            area(frameIndex) = dropletArea;
            for channelIndex = 1:length(channelsOfInterest)
                channelName = channelsOfInterest(channelIndex);
                channelIntensities(frameIndex, channelIndex) = regionProperties(label).(channelName);
            end
            if calculateRatio
                fretIndex = find(channelsOfInterest == "FRET" | channelsOfInterest == "Custom");
                cfpIndex = find(channelsOfInterest == "CFP");
                if ~fretIndex || ~cfpIndex
                    error("Could not find 'FRET' or 'CFP' channel in the list of channels of interest");
                end
                ratio(frameIndex) = channelIntensities(frameIndex, fretIndex) / channelIntensities(frameIndex, cfpIndex);
            end
        end
        data = [frames', labels', x', y', area', time', channelIntensities, ratio'];
        resultTable{dropletIndex} = array2table(data, 'VariableNames', tableColumns);
        resultTable{dropletIndex}.Properties.VariableUnits = variableUnits;
    end

end