function data = mergeDatabase(database, totalPositions, frameToMin, pixelToUm, initialPeakTimeBound)
    % Initialize empty tables for final merged data.
    mergedInfo = table();
    mergedTimeSeries = table();
    mergedCycle = table();
    mergednoOsci = table();
    
    for i = 1:length(database)
        db = database{i};
        if ismember(db.posId, totalPositions)
            %posInfo = array2table(db.posId * ones(height(db.info),1), 'VariableNames', {'POS_ID'});
            mergedInfo = [mergedInfo; db.info];
            mergedTimeSeries = [mergedTimeSeries; [array2table(db.posId * ones(height(db.timeSeries),1), 'VariableNames', {'POS_ID'}), db.timeSeries]];
            if ~isempty(db.cycle)
                mergedCycle = [mergedCycle; [array2table(db.posId * ones(height(db.cycle),1), 'VariableNames', {'POS_ID'}), db.cycle]];
            end
            if ~isempty(db.noOcillation)
                mergednoOsci = [mergednoOsci; [array2table(db.posId * ones(height(db.noOcillation),1), 'VariableNames', {'POS_ID'}), db.noOcillation]];
            end
        end
    end
    
    data.info = mergedInfo;
    data.timeSeries = mergedTimeSeries;
    data.cycle = mergedCycle;
    data.FrameToMin = frameToMin;
    data.PixelToUm = pixelToUm;
    data.InitialPeakTimeBound = initialPeakTimeBound;
    data.noOscillation = mergednoOsci;
end
