function db = processPosition(db, frameToMin, pixelToUm, initialPeakTimeBound, forceIgnore, ...
                              spermCondition, nucChannel, dnaChannel, overwriteNucMask, overwriteDNAInfo, ...
                              automaticSpermCount, hoechstoffset, FRETNumerator, FRETDenominator)
    % Analyze the TrackMate data for this position.
    [trackMate, trackPeaks, trackNoPeaks] = postprocessing.analyzeTrackMate(db, FRETNumerator, FRETDenominator, frameToMin, forceIgnore);
    
    % Ensure sperm count file exists.
    if spermCondition && ~isfile(db.spermCountCsv)
        spc_table = table([], [], 'VariableNames', {'DropID','Count'});
        writetable(spc_table, db.spermCountCsv);
    end
    spermRef = readtable(db.spermCountCsv);
    
    % Process droplet-level data.
    [timeSeriesData, cycleData, dropletInfo] = postprocessing.processDroplets(trackMate, trackPeaks, spermRef, db.posId, frameToMin, pixelToUm, initialPeakTimeBound, forceIgnore, spermCondition, nucChannel, dnaChannel, overwriteNucMask, overwriteDNAInfo, automaticSpermCount, hoechstoffset);
    
    % Save results into the database.
    db.info = [array2table(db.posId * ones(height(dropletInfo),1), 'VariableNames', {'POS_ID'}), dropletInfo];
    db.timeSeries = timeSeriesData;
    db.cycle = cycleData;
end