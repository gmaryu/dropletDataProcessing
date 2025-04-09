function data = pipelineProcess2(database, totalPositions, frameToMin, pixelToUm, initialPeakTimeBound, forceIgnore, ...
                                spermCondition, nucChannel, dnaChannel, overwriteNucMask, overwriteDNAInfo, ...
                                automaticSpermCount, hoechstoffset, FRETNumerator, FRETDenominator)
    % Loop through each database entry and process only the selected positions.
    for i = 1:length(database)
        db = database{i};
        if ismember(db.posId, totalPositions)
            db = postprocessing.processPosition(db, frameToMin, pixelToUm, initialPeakTimeBound, forceIgnore, ...
                                 spermCondition, nucChannel, dnaChannel, overwriteNucMask, overwriteDNAInfo, ...
                                 automaticSpermCount, hoechstoffset, FRETNumerator, FRETDenominator);
            database{i} = db; % update the database structure with results
        end
    end

    % Merge results across positions into final data structure.
    data = postprocessing.mergeDatabase(database, totalPositions, frameToMin, pixelToUm, initialPeakTimeBound);
end