function data = mergeDatabase(database, totalPositions, frameToMin, pixelToUm, initialPeakTimeBound)
    % Initialize empty tables for final merged data.
    mergedInfo = table();
    mergedTimeSeries = table();
    mergedCycle = table();
    
    for i = 1:length(database)
        db = database{i};
        if ismember(db.posId, totalPositions)
            posInfo = array2table(db.posId * ones(height(db.info),1), 'VariableNames', {'POS_ID'});
            mergedInfo = [mergedInfo; posInfo];
            mergedTimeSeries = [mergedTimeSeries; [array2table(db.posId * ones(height(db.timeSeries),1), 'VariableNames', {'POS_ID'}), db.timeSeries]];
            mergedCycle = [mergedCycle; [array2table(db.posId * ones(height(db.cycle),1), 'VariableNames', {'POS_ID'}), db.cycle]];
        end
    end
    
    data.info = mergedInfo;
    data.timeSeries = mergedTimeSeries;
    data.cycle = mergedCycle;
    data.FrameToMin = frameToMin;
    data.PixelToUm = pixelToUm;
    data.InitialPeakTimeBound = initialPeakTimeBound;
end
