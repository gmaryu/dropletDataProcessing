function [assignment, unassigned] = linkSegmentationResults(segmentationResult1, ...
                                                            segmentationResult2, ...
                                                            trackingParameters)
    % TODO: Add description
    arguments
        segmentationResult1 segmentation.Result
        segmentationResult2 segmentation.Result
        trackingParameters tracking.Parameters
    end
    regionProperties1 = segmentationResult1.regionProperties;
    regionProperties2 = segmentationResult2.regionProperties;
    predictions = vertcat(regionProperties1.Centroid);
    detections = vertcat(regionProperties2.Centroid);
    % Calculate the cost matrix
    cost = zeros(size(predictions, 1), size(detections, 1));
    for i = 1:size(predictions, 1)
        diff = detections - repmat(predictions(i,:), [size(detections, 1), 1]);
        cost(i,:) = sqrt(sum(diff .^ 2, 2));
    end
    % Threshold the cost matrix to avoid linking detections that are too far away
    cost(cost > trackingParameters.maxCost) = Inf;
    % Link results using the Hungarian assignment algorithm (see https://www.mathworks.com/help/vision/ref/assigndetectionstotracks.html#d126e215950)
    [assignment, unassignedTracks, unassignedDetections] = assignDetectionsToTracks(cost, ...
                                                                                    trackingParameters.costOfNonAssignment);
    unassigned = struct('firstInput', unassignedTracks, 'secondInput', unassignedDetections);
end