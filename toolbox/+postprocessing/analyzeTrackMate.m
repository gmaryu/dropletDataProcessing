function [trackMate, trackPeaks, trackNoPeaks] = analyzeTrackMate(db, ratNum, ratDen, frameToMin, forceIgnore)
% analyzeTrackMate Processes TrackMate spots data and detects periodic peaks.
%
%   [trackMate, trackPeaks, trackNoPeaks] = analyzeTrackMate(db, ratNum, ratDen, frameToMin, forceIgnore)
%
% This function reads a TrackMate spots CSV file (skipping the first row) from the
% database structure, computes a main signal as the ratio of the fields given by ratNum
% and ratDen, and then detects periodic peaks in the signal. Tracks with abnormal area
% variations are filtered out. Droplets marked for ignore (via forceIgnore) are skipped.
%
% Inputs:
%   db         - Structure with at least the fields:
%                   .trackMateSpotsCsv: Path to the CSV file.
%                   .posId: Position identifier.
%   ratNum     - (1,1) string for the numerator channel field name.
%   ratDen     - (1,1) string for the denominator channel field name.
%   frameToMin - Positive scalar conversion factor from frame to minutes.
%   forceIgnore- Table indicating droplets to ignore (fields: PosID and DropID).
%
% Outputs:
%   trackMate   - Table of processed tracking data.
%   trackPeaks  - Table with detected periodic peak data.
%   trackNoPeaks- Table with tracks that show no valid peaks.
%
% Example:
%   [tm, tp, tn] = postprocessing.analyzeTrackMate2(db, "MEAN_INTENSITY_CH5", "MEAN_INTENSITY_CH3", 6, forceIgnore);

    arguments
        db struct
        ratNum (1,1) string
        ratDen (1,1) string
        frameToMin (1,1) double {mustBePositive}
        forceIgnore table
    end


    trackMateFile = db.trackMateSpotsCsv;
    target_position = db.posId;
    
    % Read and sort tracking data (skipping first row).
    trackMate = readtable(trackMateFile);
    trackMate = sortrows(sortrows(trackMate(2:end, :), "FRAME"), "TRACK_ID");

    % Calculate FRET ratio
    flagFRET = 0;
    try
        trackMate.MAIN_SIGNAL = trackMate.(ratNum) ./ trackMate.(ratDen);
    catch
        % Return sorted TrackMate reuslt without ratio calculation
        fprintf('No FRET denominator and neumerator');
        flagFRET = 1;
         
    end

    % Peak detection
    ids = unique(trackMate.TRACK_ID);
    cnt = 0;
    trackPeaks = table();
    trackNoPeaks = table();
    
    for i = 1:length(ids)

        track = sortrows(trackMate(trackMate.TRACK_ID == ids(i), :), "FRAME");
                
        t = track.FRAME;
        preIgnored = forceIgnore.DropID(forceIgnore.PosID == target_position);
        if ~ismember(ids(i), preIgnored)
            try
                pidx = postprocessing.findPeriodicPeaks(track.MAIN_SIGNAL, frameToMin);
            catch
                if flagFRET
                    fprintf("ID:%d - No ratio column\n", ids(i));
                else
                    fprintf("ID:%d - Invalid peaks\n", ids(i));
                end
                trackNoPeaks = [trackNoPeaks; track];
                continue;
            
            end
            
            if isempty(pidx) || all(isnan(pidx(:)))
                fprintf("ID:%d  - No peaks\n", ids(i));
                trackNoPeaks = [trackNoPeaks; track];
                continue;
            end
            
            tpoints = array2table([pidx(:,1), pidx(:,2), pidx(:,3), t(pidx(:,1)), t(pidx(:,2)), t(pidx(:,3))], ...
                                   'VariableNames', {'START_INDEX','END_INDEX','TROUGH_INDEX','START_FRAME','END_FRAME','TROUGH_FRAME'});
            labels = array2table([repmat(ids(i), height(tpoints), 1), (1:height(tpoints))'], ...
                                 'VariableNames', {'TRACK_ID','CYCLE_ID'});
            trackPeaks = [trackPeaks; [labels, tpoints]];
            cnt = cnt + 1;
            %fprintf("\n");
        else
            %fprintf("ID:%d - found in force_ignored \n", ids(i));
        end
    end
    fprintf("%d / %d droplets with valid signals\n", cnt, length(ids));
end
