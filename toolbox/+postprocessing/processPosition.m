function db = processPosition(db, frameToMin, pixelToUm, initialPeakTimeBound, forceIgnore, ...
                              spermCondition, nucChannel, dnaChannel, overwriteNucMask, overwriteDNAInfo, ...
                              automaticSpermCount, hoechstoffset, FRETNumerator, FRETDenominator)
    % Analyze the TrackMate data for this position.
    [trackMate, trackPeaks, trackNoPeaks] = postprocessing.analyzeTrackMate(db, FRETNumerator, FRETDenominator, frameToMin, forceIgnore);
    
    % Ensure sperm count file exists.
    if ~isfile(db.spermCountCsv)
        spc_table = table([], [], 'VariableNames', {'DropID','Count'});
        writetable(spc_table, db.spermCountCsv);
    else
        spermRef = readtable(db.spermCountCsv);
    end
    spermRef = readtable(db.spermCountCsv);

    % Generate nuclear masks and DNA content intensity data mat files.
    if spermCondition
        % Here we call nuclearQuantification on just this position.
        postprocessing.nuclearQuantification(db, nucChannel, dnaChannel, overwriteNucMask, overwriteDNAInfo);
    end
    
    % Analyze oscillation dynamics
    if ~isempty(trackPeaks)
        % Process droplet-level data.
        [timeSeriesData, cycleData, dropletInfo] = postprocessing.processDroplets(db, trackMate, trackPeaks, spermRef, db.posId, frameToMin, pixelToUm, ...
            initialPeakTimeBound, forceIgnore, spermCondition, nucChannel, dnaChannel, automaticSpermCount, hoechstoffset);

        % Save results into the database.
        db.info = [array2table(db.posId * ones(height(dropletInfo),1), 'VariableNames', {'POS_ID'}), dropletInfo];
        db.timeSeries = timeSeriesData;
        db.cycle = cycleData;
        db.noOcillation = trackNoPeaks;
    else
        db.info = [];
        db.timeSeries = [];
        db.cycle = [];
        db.noOcillation = trackNoPeaks;
    end
    
end