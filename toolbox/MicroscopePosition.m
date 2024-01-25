classdef MicroscopePosition < handle
    %   MicroscopePosition Class to handle all the results for the collection of images in a single microscope position
    %
    %   Constructor:
    %       position = MicroscopePosition(id, channels, totalFrames, timeDelta, microscope, ...
    %                                     dataFolder, saveFolder, positionMetadata)
    %
    %   Properties:
    %       id (1,1) {mustBeInteger} % Unique identifier for this microscope position
    %       channels (1,:) string % List of channels of interest
    %       totalFrames {mustBeNumeric} % Total number of frames
    %       timeDelta {mustBeNumeric} % Time between frames in minutes
    %       microscope (1,1) string % 'new' or 'old' referring to which epifluorescence microscope was used
    %       dataFolder (1,1) string % Full path to the folder containing the data
    %       saveFolder (1,1) string % Full path to the folder where results will be saved
    %       positionMetadata (1,1) string % Metadata for this microscope position
    %       segmentationParameters (1,1) segmentation.Parameters % Parameters used for segmentation
    %       segmentationResultArray (1,:) segmentation.Result % Segmentation result for each frame
    %       trackingParameters (1,1) tracking.Parameters % Parameters used for tracking
    %       trackingResult (1,1) tracking.Result % Tracking result for this microscope position
    %
    %   Methods:
    %       pathToFile = getPathToFile(obj, channel, frame) % Returns the full path to the file for the given channel and frame
    %       segment(obj, segmentationParameters) % Performs segmentation on all frames in the microscope position
    %       plotSegmentation(obj, channel, frame) % Plots the segmentation result for the given channel and frame as a label overlay
    properties
        % Unique identifier for this microscope position
        id (1,1) {mustBeInteger} 
        % List of channels of interest
        channels (1,:) string 
        % Microscope resolution relative to 1x1 binning with 4x objective. Use 1 for 1x1 binning, 0.5 for 2x2 binning, etc.
        resolution (1,1) double
        % Total number of frames
        totalFrames {mustBeNumeric} 
        % Time between frames in minutes
        timeDelta {mustBeNumeric} 
        % 'new' or 'old' referring to which epifluorescence microscope was used
        microscope (1,1) string 
        % Full path to the folder containing the data
        dataFolder (1,1) string 
        % Full path to the folder where results will be saved
        saveFolder (1,1) string 
        % Metadata for this microscope position
        positionMetadata (1,1) string 
        % Parameters used for segmentation
        segmentationParameters (1,1) segmentation.Parameters 
        % Segmentation result for each frame
        segmentationResultArray (1,:) segmentation.Result 
        % Parameters used for tracking
        % trackingParameters (1,1) tracking.Parameters
    end

    methods
        function obj = MicroscopePosition(id, channels, resolution, totalFrames, ...
                                          timeDelta, microscope, dataFolder, saveFolder, positionMetadata)
            %   Creates a new MicroscopePosition object
            %
            %   Inputs:
            %    id (1,1) {mustBeInteger} % Unique identifier for this microscope position
            %    channels (1,:) string % List of channels of interest
            %    resolution (1,1) double % Microscope resolution relative to 1x1 binning with 4x objective. Use 1 for 1x1 binning, 0.5 for 2x2 binning, etc.
            %    totalFrames {mustBeNumeric} % Total number of frames
            %    timeDelta {mustBeNumeric} % Time between frames in minutes
            %    microscope (1,1) string % 'new' or 'old' referring to which epifluorescence microscope was used
            %    dataFolder (1,1) string % Full path to the folder containing the data
            %    saveFolder (1,1) string % Full path to the folder where results will be saved
            %    positionMetadata (1,1) string % Metadata for this microscope position
            arguments
                id (1,1) {mustBeInteger}
                channels (1,:) string
                resolution (1,1) double
                totalFrames (1,1) {mustBeInteger}
                timeDelta (1,1) double
                microscope string
                dataFolder string
                saveFolder string
                positionMetadata (1,1) string
            end
            obj.id = id;
            obj.channels = channels;
            obj.resolution = resolution;
            obj.totalFrames = totalFrames;
            obj.timeDelta = timeDelta;
            obj.microscope = microscope;
            obj.dataFolder = dataFolder;
            obj.saveFolder = saveFolder;
            obj.positionMetadata = positionMetadata;
        end

        function pathToFile = getPathToFile(obj, channel, frame)
            %   getPathToFile Returns the full path to the file for the given channel and frame
            %
            %   Inputs:
            %    channel string % Channel of interest
            %    frame (1,1) {mustBeInteger} % Frame of interest
            arguments
                obj MicroscopePosition
                channel string
                frame (1,1) {mustBeInteger}
            end
            pathToFile = fileIO.getPathToFile(obj.dataFolder, channel, frame, obj.microscope);
        end

        function segment(obj, segmentationParameters)
            %   segment Performs segmentation on all frames in the microscope position
            %
            %   Inputs:
            %    segmentationParameters segmentation.Parameters % Parameters used for segmentation
            %
            %   Notes:
            %    - This function will display a progress bar
            %    - The resolution of the microscope position must be set before calling this function
            arguments
                obj MicroscopePosition
                segmentationParameters segmentation.Parameters
            end
            segmentationParameters.resolution = obj.resolution;
            obj.segmentationParameters = segmentationParameters;
            allSegmentationResults = segmentation.Result.empty(obj.totalFrames, 0);
            fig = uifigure;
            set(fig, 'Visible', 'on');
            d = uiprogressdlg(fig,'Title','Please Wait',...
                              'Message','Opening the application');
            for frame = 1:obj.totalFrames
                d.Value = frame/obj.totalFrames;
                d.Message = sprintf('Segmenting frame %d/%d...', frame, obj.totalFrames);
                pathToFile = obj.getPathToFile('BF', frame);
                image = imread(pathToFile);
                segmentationResult = segmentation.segmentBrightFieldImage(image, segmentationParameters);
                % Add data for channels of interest
                for channel = obj.channels
                    pathToFile = obj.getPathToFile(channel, frame);
                    image = imread(pathToFile);
                    segmentationResult.addAverageChannelIntensity(image, channel);
                end
                allSegmentationResults(frame) = segmentationResult;
            end
            obj.segmentationResultArray = allSegmentationResults;
            close(fig);
        end

        function plotSegmentation(obj, channel, frame)
            %   plotSegmentation Plots the segmentation result for the given channel and frame as a label overlay
            %
            %   Inputs:
            %    channel string % Channel of interest
            %    frame (1,1) {mustBeInteger} % Frame of interest
            arguments
                obj MicroscopePosition
                channel string
                frame (1,1) {mustBeInteger}
            end
            figure;
            set(gcf, 'Visible', 'on');
            pathToFile = obj.getPathToFile(channel, frame);
            image = imread(pathToFile);
            segmentationResult = obj.segmentationResultArray(frame);
            segmentationResult.plot(image)
        end
        
        function track(obj, trackingParameters)
            % TODO: Implement
            %{
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
            fig = uifigure;
            set(fig, 'Visible', 'on');
            d = uiprogressdlg(fig,'Title','Please Wait',...
                              'Message','Opening the application');
            for frame = 1:obj.totalFrames-1
                d.Value = frame/obj.totalFrames;
                d.Message = sprintf('Tracking frame %d/%d...', frame, obj.totalFrames);
                segmentationResult1 = obj.segmentationResultArray(frame);
                segmentationResult2 = obj.segmentationResultArray(frame+1);
                linkResult = tracking.linkSegmentedObjects(segmentationResult1, ...
                                                           segmentationResult2, ...
                                                           trackingParameters);
                linkResultArray(frame) = linkResult;
            end
            close(fig);
            dropletArray = tracking.createDropletArray(linkResultArray, ...
                                                       trackingParameters);

            obj.trackingResult = tracking.Result(linkResultArray, dropletArray);
            %}
            obj.trackingParameters = trackingParameters;
        end
    end
end