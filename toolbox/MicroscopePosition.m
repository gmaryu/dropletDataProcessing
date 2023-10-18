classdef MicroscopePosition < handle
    % TODO: Add description
    properties
        channels (1,:) string % List of channels
        resolution {mustBeNumeric} % 1 for 1x1, 0.5 for 2x2, etc.
        totalFrames {mustBeNumeric} % Total number of frames
        timeDelta {mustBeNumeric} % Time between frames in minutes
        microscope string % 'new' or 'old' referring to which epifluorescence microscope was used
        pathToData string % Full path to the folder containing the data
        segmentationResultArray (1,:) segmentation.Result % Array of segmentation results for each frame
        segmentationParameters segmentation.Parameters % Parameters used for segmentation
    end

    methods
        function obj = MicroscopePosition(channels, resolution, totalFrames, ...
                                          timeDelta, microscope, pathToData)
            % TODO: Add description
            arguments
                channels (1,:) string
                resolution (1,1) double
                totalFrames (1,1) {mustBeInteger}
                timeDelta (1,1) double
                microscope string
                pathToData string
            end
            obj.channels = channels;
            obj.resolution = resolution;
            obj.totalFrames = totalFrames;
            obj.timeDelta = timeDelta;
            obj.microscope = microscope;
            obj.pathToData = pathToData;
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
    end
end