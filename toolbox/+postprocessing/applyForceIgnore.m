function dataSet = applyForceIgnore(dataSet, fi_table)

tmpInfo = dataSet.info;
tmpCycle = dataSet.cycle;
tmpTS = dataSet.timeSeries;

for f = 1:size(fi_table,1)
    targetPos = fi_table.PosID(f);
    targetTrack = fi_table.DropID(f);

    rmIndInfo = tmpInfo.POS_ID == targetPos & tmpInfo.TRACK_ID == targetTrack;
    rmIndCycle = tmpCycle.POS_ID == targetPos & tmpCycle.TRACK_ID == targetTrack;
    rmIndTS = tmpTS.POS_ID == targetPos & tmpTS.TRACK_ID == targetTrack;

    if ~ismember("IGNORED", tmpInfo.Properties.VariableNames)
        disp("IGNORED is temporally added to dataSet.Info");
        tmpInfo.IGNORED = rmIndInfo;
    end
    if ~ismember("IGNORED", tmpCycle.Properties.VariableNames)
        disp("IGNORED is temporally added to dataSet.cycle");
        tmpCycle.IGNORED = rmIndCycle;
    end
    if ~ismember("IGNORED", tmpTS.Properties.VariableNames)
        disp("IGNORED is temporally added to dataSet.timeSeries");
        tmpTS.IGNORED = rmIndTS;
    end

    tmpInfo.IGNORED = tmpInfo.IGNORED | rmIndInfo;
    tmpCycle.IGNORED = tmpCycle.IGNORED | rmIndCycle;
    tmpTS.IGNORED = tmpTS.IGNORED | rmIndTS;
end

dataSet.info = tmpInfo;
dataSet.cycle = tmpCycle;
dataSet.timeSeries = tmpTS;

end