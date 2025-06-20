function [timeSeriesData, cycleData, dropletInfo] = processDroplets(db, trackMate, trackPeaks, spermRef, posId, frameToMin, pixelToUm, initialPeakTimeBound, forceIgnore, ...
                                                                 spermCondition, hoechstCondition, automaticNucleiCount, hoechstoffset)
% processDroplets  Perform per‐droplet quantification and cycle analysis.
%
%   [timeSeriesData, cycleData, dropletInfo] = processDroplets(db, trackMate, trackPeaks, spermRef, ...
%       posId, frameToMin, pixelToUm, initialPeakTimeBound, forceIgnore, ...
%       spermCondition, hoechstCondition, automaticNucleiCount, hoechstoffset)
%
% Description:
%   For a single experimental position (`posId`), this function iterates over each droplet
%   identified in `trackPeaks`, applies optional nuclear (NLS) and DNA (Hoechst) quantification,
%   detects oscillation cycles, and assembles:
%     • a time series table (`timeSeriesData`),
%     • a cycle summary table (`cycleData`), and
%     • a droplet‐level summary (`dropletInfo`).
%
% Inputs:
%   db                    – Struct for one position, containing at least:
%                            .croppedImages     path to cropped image subfolders
%                            .spermCountCsv     (if used by getNuclearData/getDNAData)
%   trackMate             – Table of per‐frame spot data (fields include TRACK_ID, FRAME, RADIUS, etc.)
%   trackPeaks            – Table of oscillation peaks (fields include TRACK_ID, START_FRAME, etc.)
%   spermRef              – Table with manual sperm or nuclei counts (fields: DropID, Count)
%   posId                 – Scalar position identifier (must match db.posId)
%   frameToMin            – Conversion factor: frames → minutes
%   pixelToUm             – Conversion factor: pixels → μm
%   initialPeakTimeBound  – Ignore peaks whose first START_FRAME × frameToMin exceeds this bound
%   forceIgnore           – Table of droplets to skip (fields: PosID, DropID)
%   spermCondition        – Logical; if true, perform nuclear quantification via getNuclearData
%   hoechstCondition      – Logical; if true, perform DNA quantification via getDNAData
%   automaticNucleiCount  – Logical; if true, enable automatic nuclei‐count logic in getNuclearData
%   hoechstoffset         – Logical; if false, apply offset correction in getDNAData
%
% Outputs:
%   timeSeriesData  – Table concatenating all per‐frame tracking rows (possibly augmented with nuclei/DNA fields)
%   cycleData       – Table concatenating oscillation cycle rows for all droplets
%   dropletInfo     – Table with one row per droplet, columns:
%                        TRACK_ID       droplet identifier
%                        NUCLEI_COUNT   count from getNuclearData or NaN
%                        SPERM_COUNT    count from getDNAData or NaN
%                        CYCLE_NUMBER   number of detected cycles (rows in cycleData per droplet)
%                        MEDIAN_DIAMETER median droplet diameter (2×RADIUS)
%
% Workflow:
%   1. Extract the unique droplet IDs from `trackPeaks`.
%   2. Build a list of `ignoredDroplets` from `forceIgnore` for this `posId`.
%   3. For each dropletID:
%       a. Skip if in `ignoredDroplets` or if first peak starts after `initialPeakTimeBound`.
%       b. Extract its per‐frame data (`tm`) and peak data (`tp`).
%       c. Compute `medDiam` = median(tm.RADIUS * 2).
%       d. If `spermCondition`, call getNuclearData to load/compute nuclear mask results.
%       e. If `hoechstCondition`, call getDNAData to load/compute DNA intensity results.
%       f. Call processCycleData to detect cycles and append results.
%       g. Accumulate into `timeSeriesData`, `cycleData`, and a row in `dropletInfo`.
%
% Example:
%   [ts, cy, info] = processDroplets(db, trackMateTable, trackPeaksTable, spermRefTable, ...
%                                     db.posId, 6, 2.649, 100, forceIgnoreTable, ...
%                                     true, true, true, true);
%
% See also:
%   postprocessing.getSpermCount, postprocessing.getNuclearData, postprocessing.getDNAData, postprocessing.processCycleData

arguments
    db                      (1,1) struct
    trackMate               table
    trackPeaks              table
    spermRef                table
    posId                   double
    frameToMin              double
    pixelToUm               double
    initialPeakTimeBound    double
    forceIgnore             table
    spermCondition          logical
    hoechstCondition        logical
    automaticNucleiCount    logical
    hoechstoffset           logical
end

    % Initialize output containers.
    timeSeriesData = table();
    cycleData = table();
    dropletInfoRows = [];  % later converted to a table
    
    % Get unique droplet IDs (from peaks, for example).
    uniqueDropletIDs = unique(trackPeaks.TRACK_ID);
    
    % Get droplets to ignore from the forceIgnore table.
    ignoredDroplets = forceIgnore.DropID(forceIgnore.PosID == posId);
   
    % Main Part: loop with droplet ID
    for j = 1:length(uniqueDropletIDs)
        dropletID = uniqueDropletIDs(j);
        
        % Print information.
        fprintf(" - Droplet %d of Pos %d", dropletID, posId);

        % Data loading for the target droplet
        % Extract data for this droplet.
        tm = trackMate(trackMate.TRACK_ID == dropletID, :);
        tp = trackPeaks(trackPeaks.TRACK_ID == dropletID, :);
        
        % Get sperm count from spermRef. (To import manual counting information)
        spermCount = postprocessing.getSpermCount(spermRef, dropletID);
        nucleiCount = postprocessing.getSpermCount(spermRef, dropletID);
        
        % Skip analysis
        % flagged in the force ignore list.
        if ismember(dropletID, ignoredDroplets)
            fprintf(" - Already Ignored\n");
            continue;
        end
        % Check if the droplet oscillation starts too late.
        if tp.START_FRAME(1) * frameToMin > initialPeakTimeBound
            fprintf(" - Ignored. Very late oscillation.\n");
            continue;
        end

        
        % Calculate median diameter.
        medDiam = median(tm.RADIUS * 2);
        
        % Quantification of nuclear (NLS) data
        % If spermCondition true, perform nuclear quantification.
        if spermCondition
            try
                % (Assume nuclearQuantification already processes the necessary .mat files.)
                % [tm_updated, tp_updated, nucleiCount] = postprocessing.getNuclearData(db.croppedImages, dropletID, tm, tp, nucleiCount, automaticNucleiCount);
                [tm, tp, nucleiCount] = postprocessing.getNuclearData(db, dropletID, tm, tp, nucleiCount, automaticNucleiCount);
               
            catch
                %fprintf(" - Failed. getNuclearData\n");
                %disp(ME.identifier);
                %rethrow(ME);
                %fprintf(" - Ignored. .mat file not found.\n");

                continue;
            end
        else
            nucleiCount = NaN;
            fprintf(" - cytoplasm only");
        end

        % If hoechstCondition true, perform DNA quantification.
        if hoechstCondition
            try
                % (Assume nuclearQuantification already processes the necessary .mat files.)
                % run detectMultiNuclei
                % [tm_updated, tp_updated, spermCount] = postprocessing.getDNAData(db.croppedImages, dropletID, db.posId, tm_updated, tp_updated, spermCount, hoechstoffset);
                [tm, tp, spermCount, nucleiCount] = postprocessing.getDNAData(db, dropletID, db.posId, tm, tp, spermCount, hoechstoffset);
                
            catch ME
                fprintf(" - Failed. getDNAData\n");
                %fprintf(" - Ignored. .mat file not found.\n");
                disp(ME.identifier);
                %rethrow(ME);
                
                
                continue;
            end

        else
            spermCount = NaN;
            fprintf(" - No DNA staining -");
        end

        
        % Process cycle data for current droplet.
       
        %[tp_updated, cycleMetrics] = postprocessing.processCycleData(tp_updated, tm_updated, frameToMin, pixelToUm, spermCondition);
        [tp, cycleMetrics] = postprocessing.processCycleData(tp, tm, frameToMin, pixelToUm, spermCondition);
        
        % Gather processed data.
        timeSeriesData = [timeSeriesData; tm];
        % timeSeriesData = [timeSeriesData; tm_updated];
        cycleData = [cycleData; tp];
        %cycleData = [cycleData; tp_updated];
        dropletInfoRows = [dropletInfoRows; [dropletID, nucleiCount, spermCount, size(tp,1)-1, medDiam]]; %#ok<AGROW>
        %dropletInfoRows = [dropletInfoRows; [dropletID, nucleiCount, spermCount, size(tp_updated,1)-1, medDiam]]; %#ok<AGROW>
        
        fprintf(" - \n");
    end

    if ~isempty(dropletInfoRows)
        % Convert droplet info into a table with appropriate variable names.
        dropletInfo = array2table(dropletInfoRows, 'VariableNames', {'TRACK_ID','NUCLEI_COUNT','SPERM_COUNT','CYCLE_NUMBER','MEDIAN_DIAMETER'});
    else
        dropletInfo = array2table(zeros(0,5), 'VariableNames', {'TRACK_ID','NUCLEI_COUNT','SPERM_COUNT','CYCLE_NUMBER','MEDIAN_DIAMETER'});
    end
end
