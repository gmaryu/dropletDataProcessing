classdef MicroscopePosition < handle
    % TODO: Add description
    properties
        id (1,1) {mustBeInteger} % Unique identifier for this microscope position
        channels (1,:) string % List of channels
        resolution {mustBeNumeric} % 1.0 for 1x1, 0.5 for 2x2, etc.
        totalFrames {mustBeNumeric} % Total number of frames
        timeDelta {mustBeNumeric} % Time between frames in minutes
        microscope (1,1) string % 'new' or 'old' referring to which epifluorescence microscope was used
        pathToData (1,1) string % Full path to the folder containing the data
        positionMetadata (1,1) string % Metadata for this microscope position
        segmentationResultArray (1,:) segmentation.Result % Array of segmentation results for each frame
        segmentationParameters (1,1) segmentation.Parameters % Parameters used for segmentation
        trackingParameters (1,1) tracking.Parameters % Parameters used for tracking
        trackingResult (1,1) tracking.Result % Tracking result for this microscope position
    end

    methods
        function obj = MicroscopePosition(id, channels, resolution, totalFrames, ...
                                          timeDelta, microscope, pathToData, positionMetadata)
            % TODO: Add description
            arguments
                id (1,1) {mustBeInteger}
                channels (1,:) string
                resolution (1,1) double
                totalFrames (1,1) {mustBeInteger}
                timeDelta (1,1) double
                microscope string
                pathToData string
                positionMetadata (1,1) string
            end
            obj.id = id;
            obj.channels = channels;
            obj.resolution = resolution;
            obj.totalFrames = totalFrames;
            obj.timeDelta = timeDelta;
            obj.microscope = microscope;
            obj.pathToData = pathToData;
            obj.positionMetadata = positionMetadata;
        end

        function pathToFile = getPathToFile(obj, channel, frame)
            % TODO: Add description
            arguments
                obj MicroscopePosition
                channel string
                frame (1,1) {mustBeInteger}
            end
            pathToFile = fileIO.getPathToFile(obj.pathToData, channel, frame, obj.microscope);
        end

        function obj = segment(obj, segmentationParameters)
            % TODO: Add description
            arguments
                obj MicroscopePosition
                segmentationParameters segmentation.Parameters
            end
            obj.segmentationParameters = segmentationParameters;
            allSegmentationResults = segmentation.Result.empty(obj.totalFrames, 0);
            for frame = 1:obj.totalFrames
                pathToFile = obj.getPathToFile('BF', frame);
                image = imread(pathToFile);
                allSegmentationResults(frame) = segmentation.segmentBrightFieldImage(image, segmentationParameters);
            end
            obj.segmentationResultArray = allSegmentationResults;
        end

        function plotSegmentation(obj, channel, frame)
            % TODO: Add description
            arguments
                obj MicroscopePosition
                channel string
                frame (1,1) {mustBeInteger}
            end
            pathToFile = obj.getPathToFile(channel, frame);
            image = imread(pathToFile);
            segmentationResult = obj.segmentationResultArray(frame);
            segmentationResult.plot(image)
        end

        function track(obj, trackingParameters)
            % TODO: Add description
            arguments
                obj MicroscopePosition
                trackingParameters tracking.Parameters
            end
            % Check that segmentation has been performed
            if isempty(obj.segmentationResultArray)
                error('Segmentation has not been performed for this microscope position. Please run segment() first.')
            end
            obj.trackingParameters = trackingParameters;

            linkResultArray = tracking.LinkResult.empty(obj.totalFrames-1, 0);
            % Link segmentation results
            for frame = 1:obj.totalFrames-1
                segmentationResult1 = obj.segmentationResultArray(frame);
                segmentationResult2 = obj.segmentationResultArray(frame+1);
                linkResult = tracking.linkSegmentedObjects(segmentationResult1, ...
                                                           segmentationResult2, ...
                                                           trackingParameters);
                linkResultArray(frame) = linkResult;
            end
            dropletArray = tracking.createDropletArray(linkResultArray, ...
                                                       trackingParameters);
            obj.trackingResult = tracking.Result(linkResultArray, dropletArray);
        end
    end
end