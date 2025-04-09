function [trackMate, trackPeaks, trackNoPeaks] = analyzeTrackMate(db, ratNum, ratDen, frameToMin, forceIgnore)
% analyzeTrackMate2 Processes TrackMate spots data and detects periodic peaks.
%
%   [trackMate, trackPeaks, trackNoPeaks] = analyzeTrackMate2(db, ratNum, ratDen, frameToMin, forceIgnore)
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

    debug = false;

    trackMateFile = db.trackMateSpotsCsv;
    target_position = db.posId;
    
    % Read and sort tracking data (skipping first row).
    trackMate = readtable(trackMateFile);
    trackMate = sortrows(sortrows(trackMate(2:end, :), "FRAME"), "TRACK_ID");
    trackMate.MAIN_SIGNAL = trackMate.(ratNum) ./ trackMate.(ratDen);
    
    ids = unique(trackMate.TRACK_ID);
    cnt = 0;
    trackPeaks = table();
    trackNoPeaks = table();
    
    for i = 1:length(ids)
        if debug
            fprintf("%d", ids(i));
        end
        track = sortrows(trackMate(trackMate.TRACK_ID == ids(i), :), "FRAME");
        
        % Filter tracks with abnormal area variations.
        areaAnomaly = abs((track.AREA - median(track.AREA)) / median(track.AREA));
        if sum(areaAnomaly > 0.05) > 0.01 * numel(areaAnomaly)
            validIdx = areaAnomaly < 0.05;
            origLen = height(track);
            track = track(validIdx, :);
            if height(track) < 0.8 * origLen
                fprintf("ID:%d - Area anomaly (%.3f, %.3f, %d)\n", ids(i), median(areaAnomaly), std(areaAnomaly), sum(areaAnomaly > 0.1));
                continue;
            end
        end
        
        t = track.FRAME;
        preIgnored = forceIgnore.DropID(forceIgnore.PosID == target_position);
        if ~ismember(ids(i), preIgnored)
            try
                pidx = findPeriodicPeaks(track.MAIN_SIGNAL, frameToMin);
            catch
                fprintf("ID:%d - Invalid peaks\n", ids(i));
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
            fprintf("ID:%d - found in force_ignored \n", ids(i));
        end
    end
    fprintf("%d / %d droplets with valid signals\n", cnt, length(ids));
end

function retv = findPeriodicPeaks(signal, frameToMin)
% findPeriodicPeaks Detects periodic peaks in a signal.
%
%   retv = findPeriodicPeaks(signal, frameToMin)
%
% The function uses findpeaks on both the signal and its negative to determine peak and trough
% positions. It returns a matrix with each row as [start_index, end_index, trough_index] if the peaks
% and troughs match the expected pattern; otherwise, it returns NaN.

    p = 0.1;
    maxw = 60 / frameToMin;  % Maximum expected peak width in frames
    
    [~, ip] = findpeaks(signal, "MinPeakProminence", p, "MaxPeakWidth", maxw);
    [~, it] = findpeaks(-signal, "MinPeakProminence", p);
    
    if numel(ip) == numel(it)
        if all(it - ip > 0)
            retv = [ip(1:end-1), ip(2:end), it(1:end-1)];
        elseif all(it - ip < 0)
            retv = [ip(1:end-1), ip(2:end), it(2:end)];
        else
            retv = nan;
        end
    elseif numel(ip) == numel(it) + 1
        if all(it - ip(1:end-1) > 0)
            retv = [ip(1:end-1), ip(2:end), it(1:end)];
        else
            retv = nan;
        end
    else
        retv = nan;
    end
end
