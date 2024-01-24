function trackingResult = trackDroplets(segmentationResultArray, trackingParameters, verbose)
    %   trackDroplets Link segmented droplets into tracks.
    %
    %   trackingResult = trackDroplets(segmentationResultArray, trackingParameters)
    %
    %   Required Inputs:
    %       segmentationResultArray (1,:) segmentation.Result - Segmentation result for each frame
    %       trackingParameters (1,1) tracking.Parameters - Tracking parameters
    %
    %   Optional Inputs:
    %       verbose (1,1) logical = true - Enable verbose output
    %
    %   Output:
    %       trackingResult (1,1) tracking.Result - Tracking result containing linkResultArray and dropletArray
    arguments
        % Segmentation result for each frame
        segmentationResultArray (1,:) segmentation.Result
        % Tracking parameters
        trackingParameters (1,1) tracking.Parameters
        % Whether to enable verbose output
        verbose (1,1) logical = true
    end
        totalFrames = length(segmentationResultArray);
        linkResultArray = tracking.LinkResult.empty(totalFrames-1, 0);
        if verbose
            fig = uifigure;
            set(fig, 'Visible', 'on');
            d = uiprogressdlg(fig,'Title','Please Wait',...
                                'Message','Opening the application');
        end
        for frame = 1:totalFrames-1
            if verbose
                d.Value = frame/totalFrames;
                d.Message = sprintf('Tracking frame %d/%d...', frame, totalFrames);
            end
            segmentationResult1 = segmentationResultArray(frame);
            segmentationResult2 = segmentationResultArray(frame + 1);
            linkResult = tracking.linkSegmentedObjects(segmentationResult1, ...
                                                       segmentationResult2, ...
                                                       trackingParameters);
            linkResultArray(frame) = linkResult;
        end
        if verbose
            close(fig);
        end
        dropletArray = tracking.createDropletArray(linkResultArray, ...
                                                   trackingParameters);
        trackingResult = tracking.Result(linkResultArray, dropletArray);
end