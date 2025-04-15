function [timeSeriesData, cycleData, dropletInfo] = processDroplets(db, trackMate, trackPeaks, spermRef, posId, frameToMin, pixelToUm, initialPeakTimeBound, forceIgnore, ...
                                                                 spermCondition, hoechstCondition, automaticNucleiCount, hoechstoffset)

arguments
        db (1,1) struct
        trackMate table
        trackPeaks table
        spermRef table
        posId double
        frameToMin double
        pixelToUm double
        initialPeakTimeBound double
        forceIgnore table
        spermCondition logical
        hoechstCondition logical
        automaticNucleiCount logical
        hoechstoffset logical
    end

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
        
        % Get sperm count from spermRef. (To import manual counting information)
        spermCount = postprocessing.getSpermCount(spermRef, dropletID);
        nucleiCount = postprocessing.getSpermCount(spermRef, dropletID);
        
        % Calculate median diameter.
        medDiam = median(tm.RADIUS * 2);
        
        % Print information.
        fprintf(" - Droplet %d of Pos %d", dropletID, posId);
        
        % Skip droplet if flagged in the force ignore list.
        if ismember(dropletID, ignoredDroplets)
            fprintf(" - Already Ignored\n");
            continue;
        end

        % Check if the droplet oscillation starts too late.
        if tp.START_FRAME(1) * frameToMin > initialPeakTimeBound
            fprintf(" - Ignored. Very late oscillation.\n");
            continue;
        end
        
        % Process cycle data for current droplet.
        [tp_updated, cycleMetrics] = postprocessing.processCycleData(tp, tm, frameToMin, pixelToUm, spermCondition);

        % If spermCondition true, perform nuclear quantification.
        if spermCondition
            try
                % (Assume nuclearQuantification already processes the necessary .mat files.)
                [tm, tp_updated, nucleiCount] = postprocessing.getNuclearData(db.croppedImages, dropletID, tm, tp_updated, nucleiCount, automaticNucleiCount);
                fprintf(" - Nuclear mask obtained");
            catch
                %fprintf(" - Ignored. .mat file not found.\n");
                fprintf(" - Failed. getNuclearData\n");
                continue;
            end
        else
            nucleiCount = NaN;
            fprintf(" - cytoplasm only -");
        end

        % If spermCondition true, perform nuclear quantification.
        if hoechstCondition
            try
                % (Assume nuclearQuantification already processes the necessary .mat files.)
                % run detectMultiNuclei
                [tm, tp_updated, spermCount] = postprocessing.getDNAData(db.croppedImages, dropletID, tm, tp_updated, spermCount, hoechstoffset);
                fprintf(" - DNA mask obtained");
            catch
                %fprintf(" - Ignored. .mat file not found.\n");
                fprintf(" - Failed. getDNAData\n");
                
                continue;
            end
        else
            spermCount = NaN;
            fprintf(" - cytoplasm only -");
        end
        
        
        % Gather processed data.
        timeSeriesData = [timeSeriesData; tm];
        cycleData = [cycleData; tp_updated];
        dropletInfoRows = [dropletInfoRows; [dropletID, nucleiCount, medDiam]]; %#ok<AGROW>
        
        fprintf(" - \n");
    end
    
    % Convert droplet info into a table with appropriate variable names.
    dropletInfo = array2table(dropletInfoRows, 'VariableNames', {'TRACK_ID','NUCLEI_COUNT','MEDIAN_DIAMETER'});
end
