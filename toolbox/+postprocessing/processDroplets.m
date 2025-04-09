function [timeSeriesData, cycleData, dropletInfo] = processDroplets(db, trackMate, trackPeaks, spermRef, posId, frameToMin, pixelToUm, initialPeakTimeBound, forceIgnore, ...
                                                                 spermCondition, nucChannel, dnaChannel, automaticSpermCount, hoechstoffset)
    % Initialize output containers.
    timeSeriesData = table();
    cycleData = table();
    dropletInfoRows = [];  % later converted to a table
    
    % Get unique droplet IDs (from peaks, for example).
    uniqueDropletIDs = unique(trackPeaks.TRACK_ID);
    
    % Get droplets to ignore from the forceIgnore table.
    ignoredDroplets = forceIgnore.DropID(forceIgnore.PosID == posId);
    
    for j = 1:length(uniqueDropletIDs)
        dropletID = uniqueDropletIDs(j);
        
        % Extract data for this droplet.
        tm = trackMate(trackMate.TRACK_ID == dropletID, :);
        tp = trackPeaks(trackPeaks.TRACK_ID == dropletID, :);
        
        % Get sperm count from spermRef.
        spermCount = postprocessing.getSpermCount(spermRef, dropletID);
        
        % Calculate median diameter.
        medDiam = median(tm.RADIUS * 2);
        
        % Print information.
        fprintf(" - Droplet %d of Pos %d", dropletID, posId);
        
        % Skip droplet if flagged in the force ignore list.
        if ismember(dropletID, ignoredDroplets)
            fprintf(" - short tracking frames\n");
            continue;
        end

        % Check if the droplet oscillation starts too late.
        if tp.START_FRAME(1) * frameToMin > initialPeakTimeBound
            fprintf(" - Ignored. Very late oscillation.\n");
            continue;
        end
        
        % If spermCondition true, perform nuclear quantification.
        if spermCondition
            try
                % (Assume nuclearQuantification already processes the necessary .mat files.)
                nuclearData = postprocessing.getNuclearData(db.croppedImages, posId, dropletID, nucChannel, dnaChannel, automaticSpermCount, hoechstoffset);
                % Append the obtained data to the tracking table.
                tm.NPIXEL_NUC = nuclearData.nuclearArea';
                tm.NPIXEL_DNA = nuclearData.hoechstNPixels';
                tm.SUMINTENSITY_DNA = nuclearData.hoechstSum';
                fprintf(" - Nuclear mask obtained");
            catch
                fprintf(" - Ignored. .mat file not found.\n");
                continue;
            end
        else
            fprintf(" - cytoplasm only -");
        end

        
        % Process cycle data for current droplet.
        [tp_updated, cycleMetrics] = postprocessing.processCycleData(tp, tm, frameToMin, pixelToUm, spermCondition);
        
        % Gather processed data.
        timeSeriesData = [timeSeriesData; tm];
        cycleData = [cycleData; tp_updated];
        dropletInfoRows = [dropletInfoRows; [dropletID, spermCount, medDiam]]; %#ok<AGROW>
        
        fprintf(" - \n");
    end
    
    % Convert droplet info into a table with appropriate variable names.
    dropletInfo = array2table(dropletInfoRows, 'VariableNames', {'TRACK_ID','SPERM_COUNT','MEDIAN_DIAMETER'});
end
