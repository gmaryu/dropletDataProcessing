function dropletArray = createDropletArray(linkResultArray, trackingParameters)
    %   createDropletArray Connects linkResults together forming tracks for individual droplets
    %
    %   dropletArray = tracking.createDropletArray(linkResultArray, trackingParameters);
    %
    %   Notes:
    %       Tracks of length 1 frame that start on the first frame are not incuded in the result
    %
    %   Inputs:
    %       linkResultArray (1,:) tracking.LinkResult - Array with connections among labels for each frame
    %       trackingParameters (1,1) tracking.Parameters - Parameters for the tracking algorithm
    %
    %   Output:
    %       dropletArray (1,:) - Array containing tracking.Droplet objects
    arguments
        % Array with connections among labels for each frame
        linkResultArray (1,:) tracking.LinkResult
        % Parameters for the tracking algorithm
        trackingParameters (1,1) tracking.Parameters
    end
    referenceLinkResult = linkResultArray(1);
    % Total number of tracks = 
    % Assigned tracks on the first frame + 
    % Unassigned tracks from second inputs (these create new tracks in the next link)
    totalTrackNumber = size(referenceLinkResult.assigned, 1);
    for i = 1:length(linkResultArray)-1
        linkResult = linkResultArray(i);
        totalTrackNumber = totalTrackNumber + length(linkResult.unassignedFromSecondInput);
    end
    trackList = NaN(totalTrackNumber, length(linkResultArray)+1);
    % Calculation for tracks assigned on the first frame
    for selectedRow = 1:size(referenceLinkResult.assigned, 1)
        label = referenceLinkResult.assigned(selectedRow, 1);
        associatedLabel = referenceLinkResult.assigned(selectedRow, 2);
        trackList(selectedRow, 1) = label;
        for i = 2:length(linkResultArray)
            linkResult = linkResultArray(i);
            assigned = linkResult.assigned;
            if any(assigned(:,1) == associatedLabel)
                trackList(selectedRow, i) = associatedLabel;
                row = assigned(:,1) == associatedLabel;
                associatedLabel = assigned(row, 2);
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
    % Calculation for unassigned tracks from second inputs
    if ~trackingParameters.mustStartOnFirstFrame
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
    end
    % Remove tracks that are shorter than the minimum length
    trackLengths = sum(~isnan(trackList), 2);
    trackList(trackLengths < trackingParameters.minTrackLength, :) = [];
    % Create the dropletArray
    dropletArray = tracking.Droplet.empty(length(trackList), 0);
    for i = 1:length(trackList)
        id = i;
        labels = trackList(i, :);
        dropletArray(i) = tracking.Droplet(id, labels);
    end
end