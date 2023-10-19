function dropletArray = createDropletArray(linkResultArray, trackingParameters)
    % TODO: Add description
    % TODO: Threshold by the minimum length of track desired
    % TODO: Rename variables to make the code more understandable
    arguments
        linkResultArray (1,:) tracking.LinkResult
        trackingParameters (1,1) tracking.Parameters
    end

    referenceLinkResult = linkResultArray(1);
    % The total number of tracks corresponds to all those segments that are assigned in the first
    % frame plus all those segments that are unassigned in the subsequent frames. Note that the
    % unassigned segments on the second frame are found in the unassignedFromSecondInput property
    % of the first linkResult.
    totalTrackNumber = size(referenceLinkResult.assigned, 1);
    for i = 1:length(linkResultArray)-1
        linkResult = linkResultArray(i);
        totalTrackNumber = totalTrackNumber + length(linkResult.unassignedFromSecondInput);
    end
    trackList = NaN(totalTrackNumber, length(linkResultArray)+1);
    % Tracks that start on the first frame
    for selectedRow = 1:size(referenceLinkResult.assigned, 1)
        label = referenceLinkResult.assigned(selectedRow,1);
        associatedLabel = referenceLinkResult.assigned(selectedRow,2);
        trackList(selectedRow, 1) = label;
        for i = 2:length(linkResultArray)
            linkResult = linkResultArray(i);
            assigned = linkResult.assigned;
            if any(assigned(:,1) == associatedLabel)
                trackList(selectedRow, i) = associatedLabel;
                row = assigned(:,1) == associatedLabel;
                associatedLabel = assigned(row,2);
                % If we are at the last linkResult, then add the associatedLabel
                if i == length(linkResultArray)
                    trackList(selectedRow, i+1) = associatedLabel;
                end
            else
                trackList(selectedRow, i) = associatedLabel;
                break
            end
        end
    end
    % Tracks that don't start on the first frame
    offset = size(referenceLinkResult.assigned, 1);
    for startingLinkResult = 1:length(linkResultArray)-1
        if startingLinkResult > 1
            previousLinkResult = linkResultArray(startingLinkResult-1);
            previousRows = length(previousLinkResult.unassignedFromSecondInput);
            offset = offset + previousRows;
        end
        referenceLinkResult = linkResultArray(startingLinkResult);
        for selectedRow = 1:length(referenceLinkResult.unassignedFromSecondInput)
            associatedLabel = referenceLinkResult.unassignedFromSecondInput(selectedRow);
            for i = startingLinkResult+1:length(linkResultArray)
                linkResult = linkResultArray(i);
                assigned = linkResult.assigned;
                if any(assigned(:,1) == associatedLabel)
                    trackList(selectedRow + offset, i) = associatedLabel;
                    row = assigned(:,1) == associatedLabel;
                    associatedLabel = assigned(row,2);
                    % If we are at the last linkResult, then add the associatedLabel
                    if i == length(linkResultArray)
                        trackList(selectedRow + offset, i+1) = associatedLabel;
                    end
                else
                    trackList(selectedRow + offset, i) = associatedLabel;
                    break
                end
            end
        end
    end
    % Remove tracks that are shorter than the minimum length
    trackLengths = sum(~isnan(trackList), 2);
    trackList(trackLengths < trackingParameters.minTrackLength, :) = [];
    % If the tracking must start on the first frame, then remove tracks that don't
    if trackingParameters.mustStartOnFirstFrame
        trackList(isnan(trackList(:,1)), :) = [];
    end
    % Create the dropletArray
    dropletArray = tracking.Droplet.empty(length(trackList), 0);
    for i = 1:length(trackList)
        id = i;
        labels = trackList(i, :);
        dropletArray(i) = tracking.Droplet(id, labels);
    end
end