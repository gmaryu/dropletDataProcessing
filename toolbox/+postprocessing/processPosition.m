function data_output = processPosition(db, frameToMin, pixelToUm, initialPeakTimeBound, forceIgnore, ...
                              spermCondition, hoechstCondition, nucChannel, dnaChannel, overwriteNucMask, overwriteDNAInfo, ...
                              automaticNucleiCount, hoechstoffset, FRETNumerator, FRETDenominator)
% processPosition Process signal quantification and segmentation data for a single experimental position.
%
%   db = processPosition(db, frameToMin, pixelToUm, initialPeakTimeBound, forceIgnore, ...
%                          spermCondition, hoechstCondition, nucChannel, dnaChannel, overwriteNucMask, overwriteDNAInfo, ...
%                          automaticNucleiCount, hoechstoffset, FRETNumerator, FRETDenominator)
%
% Description:
%   This function performs a comprehensive processing workflow for droplet tracking data
%   associated with a specific experimental position contained in the structure "db." The workflow includes:
%     1. Analyzing TrackMate output to obtain droplet tracking data, droplet oscillation peaks,
%        and droplets with no detected peaks.
%     2. Ensuring that a manual sperm count CSV file exists and reading its data.
%     3. Performing nuclear segmentation and DNA intensity quantification if spermCondition is true.
%        This is achieved by calling the nuclearSegmentation routine in the postprocessing package,
%        which in turn utilizes routines such as cropBrightChunk and cropDNAMask.
%     4. If oscillation peaks (trackPeaks) are present, processing the droplet-level data (oscillation dynamics)
%        via the processDroplets function to extract time series data, cycle (oscillation) characteristics, and
%        droplet-specific information.
%     5. Updating the input structure "db" with the new processed information:
%            - db.info: A table with droplet-level summary (including position, droplet ID, sperm count, cycle number and median diameter).
%            - db.timeSeries: A table containing the droplet time series (tracking dynamics).
%            - db.cycle: A table containing details of the oscillation cycles.
%            - db.noOcillation: Droplets for which no oscillation peaks were detected.
%
% Inputs:
%   db                 - (1×1 struct) Structure containing all data for one experimental position.
%                        The structure must include fields such as posId, trackMate data path, sperm count CSV path, and croppedImages.
%   frameToMin         - (double scalar) Conversion factor to convert frame numbers to minutes.
%   pixelToUm          - (double scalar) Conversion factor to convert pixels to micrometers.
%   initialPeakTimeBound - (double scalar) Time (in minutes) beyond which droplet oscillation peaks are ignored.
%   forceIgnore        - (table) Table listing droplets to ignore in the analysis (e.g., due to short tracking).
%   spermCondition     - (logical) Flag indicating whether to perform nuclear segmentation and DNA intensity analysis.
%   hoechstCondition   - (logical) Flag indicating whether Hoechst intensity measurements are used in the analysis.
%   nucChannel         - (string scalar) Name of the channel used for nuclear staining (e.g., "CFP").
%   dnaChannel         - (string scalar) Name of the channel used for DNA (Hoechst) staining (e.g., "DAPI").
%   overwriteNucMask   - (logical) If true, any existing nuclear segmentation masks will be overwritten.
%   overwriteDNAInfo   - (logical) If true, any existing DNA intensity data will be recalculated.
%   automaticNucleiCount - (logical) Flag indicating whether automatic nuclei counting is performed.
%   hoechstoffset      - (logical) Flag to determine whether to apply an offset correction to Hoechst intensity.
%   FRETNumerator      - (string scalar) Name of the field used as the numerator for FRET calculations.
%   FRETDenominator    - (string scalar) Name of the field used as the denominator for FRET calculations.
%
% Outputs:
%   db                 - (1×1 struct) The input structure "db" updated with processed data. New fields added:
%         .info          - Table containing droplet-level summary information (e.g., position ID, droplet ID, sperm count, median diameter).
%         .timeSeries    - Table containing the droplet tracking (time series) data.
%         .cycle         - Table containing oscillation (cycle) data.
%         .noOcillation  - Table (or array) for droplets that did not show any oscillation peaks.
%
% Processing Steps:
%   1. Calls postprocessing.analyzeTrackMate to extract droplet tracking information and identify oscillation peaks.
%   2. Ensures that the sperm count file exists; if not, creates an empty table and writes it.
%   3. When spermCondition is true, calls postprocessing.nuclearSegmentation for nuclear mask generation and
%      DNA quantification.
%   4. If oscillation peaks are detected, calls postprocessing.processDroplets to process droplet-level dynamics,
%      which returns time series, cycle, and droplet summary data.
%   5. Updates the input structure "db" with all newly computed data.
%
% Example:
%   % Given an existing structure "db" for a position, and parameters defined as:
%   frameToMin         = 6;
%   pixelToUm          = 2.649;
%   initialPeakTimeBound = 100;
%   forceIgnore        = <your forceIgnore table>;
%   spermCondition     = true;
%   hoechstCondition   = true;
%   nucChannel         = "CFP";
%   dnaChannel         = "DAPI";
%   overwriteNucMask   = true;
%   overwriteDNAInfo   = true;
%   automaticNucleiCount = true;
%   hoechstoffset      = true;
%   FRETNumerator      = "MEAN_INTENSITY_CH5";
%   FRETDenominator    = "MEAN_INTENSITY_CH3";
%
%   % Process the position:
%   db = processPosition(db, frameToMin, pixelToUm, initialPeakTimeBound, forceIgnore, ...
%                        spermCondition, hoechstCondition, nucChannel, dnaChannel, overwriteNucMask, overwriteDNAInfo, ...
%                        automaticNucleiCount, hoechstoffset, FRETNumerator, FRETDenominator);
%
% See also: postprocessing.analyzeTrackMate, postprocessing.nuclearSegmentation, postprocessing.processDroplets


arguments
    db (1,1) struct
    frameToMin double
    pixelToUm double
    initialPeakTimeBound double
    forceIgnore table
    spermCondition logical
    hoechstCondition logical
    nucChannel string
    dnaChannel string
    overwriteNucMask logical
    overwriteDNAInfo logical
    automaticNucleiCount logical
    hoechstoffset logical
    FRETNumerator string
    FRETDenominator string
end

    %% Analyze the TrackMate data for this position.
    [trackMate, trackPeaks, trackNoPeaks] = postprocessing.analyzeTrackMate(db, FRETNumerator, FRETDenominator, frameToMin, forceIgnore);

    % Ensure sperm count file exists.
    if ~isfile(db.spermCountCsv)
        spc_table = table([], [], 'VariableNames', {'DropID','Count'});
        writetable(spc_table, db.spermCountCsv);
    else
        spermRef = readtable(db.spermCountCsv);
    end
    spermRef = readtable(db.spermCountCsv);

    %% Generate nuclear masks and DNA mask mat files.
    if spermCondition
        % Here we call nuclearQuantification on just this position to segment nuclei and DNA.
        % cropBrightChunk
        % cropDNAMask
        segmentation.nuclearSegmentation(db, nucChannel, dnaChannel, overwriteNucMask, overwriteDNAInfo);
    end
    

    %% Analyze oscillation dynamics
    if ~isempty(trackPeaks)
        % Process droplet-level data.
        [timeSeriesData, cycleData, dropletInfo] = postprocessing.processDroplets(db, trackMate, trackPeaks, spermRef, db.posId, frameToMin, ...
            pixelToUm, initialPeakTimeBound, forceIgnore, spermCondition, hoechstCondition, automaticNucleiCount, hoechstoffset);

        % Save results into the database.
        data_output.posId = db.posId;
        data_output.info = [array2table(db.posId * ones(height(dropletInfo),1), 'VariableNames', {'POS_ID'}), dropletInfo];
        data_output.timeSeries = timeSeriesData;
        data_output.cycle = cycleData;
        data_output.noOcillation = trackNoPeaks;
    else
        data_output.posId = db.posId;
        data_output.info = [];
        data_output.timeSeries = [];
        data_output.cycle = [];
        data_output.noOcillation = trackNoPeaks;
    end
    
    
end