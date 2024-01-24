function linkResult = linkSegmentedObjects(segmentationResult1, segmentationResult2, ...
                                           trackingParameters)
    %   linkSegmentedObjects Links two segmentation results by matching labels from the first result to the second result. The matching is done using the Hungarian assignment algorithm.
    %
    %   linkResult = tracking.linkSegmentedObjects(segmentationResult1, segmentationResult2, trackingParameters);
    %
    %   Inputs:
    %       segmentationResult1 (1,1) segmentation.Result - A segmentation.Result object containing the segmentation result for a particular frame.
    %       segmentationResult2 (1,1) segmentation.Result - A segmentation.Result object containing the segmentation result for the consecutive frame.
    %       trackingParameters (1,1) tracking.Parameters - A tracking.Parameters object containing the parameters for the tracking algorithm.
    %
    %   Outputs:
    %       linkResult: A tracking.LinkResult object containing the label assignment as well as those labels that were not assigned.
    arguments
        segmentationResult1 (1,1) segmentation.Result
        segmentationResult2 (1,1) segmentation.Result
        trackingParameters (1,1) tracking.Parameters
    end
    predictions = vertcat(segmentationResult1.regionProperties.Centroid);
    detections = vertcat(segmentationResult2.regionProperties.Centroid);
    % Calculate the cost matrix (euclidean distance)
    cost = pdist2(predictions, detections);
    % Threshold the cost matrix to avoid linking detections that are too far away
    cost(cost > trackingParameters.maxCost) = Inf;
    % Link results using the Hungarian assignment algorithm (see https://www.mathworks.com/help/vision/ref/assigndetectionstotracks.html#d126e215950)
    [assigned, unassignedTracks, unassignedDetections] = assignDetectionsToTracks(cost, trackingParameters.costOfNonAssignment);
    linkResult = tracking.LinkResult(assigned, unassignedTracks, unassignedDetections);
end