function dataSet = dropPosition(dataSet, targetPos)
% dropPosition - Mark all entries with a given POS_ID as ignored
%
% Inputs:
%   dataSet   : struct with fields .info, .cycle, .timeSeries (all tables)
%   targetPos : numeric POS_ID to ignore
%
% Output:
%   dataSet   : same struct, with IGNORED flag updated

% Ensure an IGNORED column exists in each table (initialize with zeros if missing)
if ~ismember("IGNORED", dataSet.info.Properties.VariableNames)
    dataSet.info.IGNORED = zeros(size(dataSet.info,1),1);
end
if ~ismember("IGNORED", dataSet.cycle.Properties.VariableNames)
    dataSet.cycle.IGNORED = zeros(size(dataSet.cycle,1),1);
end
if ~ismember("IGNORED", dataSet.timeSeries.Properties.VariableNames)
    dataSet.timeSeries.IGNORED = zeros(size(dataSet.timeSeries,1),1);
end

% Find rows in info table that match the target POS_ID
indexInfo = dataSet.info.POS_ID == targetPos;
% Mark those rows as ignored (logical OR to preserve previous flags)
dataSet.info.IGNORED = dataSet.info.IGNORED | indexInfo;

% Repeat the same process for cycle table
indexCycle = dataSet.cycle.POS_ID == targetPos;
dataSet.cycle.IGNORED = dataSet.cycle.IGNORED | indexCycle;

% Repeat the same process for timeSeries table
indexTimeSeries = dataSet.timeSeries.POS_ID == targetPos;
dataSet.timeSeries.IGNORED = dataSet.timeSeries.IGNORED | indexTimeSeries;

end