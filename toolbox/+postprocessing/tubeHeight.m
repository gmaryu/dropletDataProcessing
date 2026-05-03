function dataSet = tubeHeight(dataSet, pos_100, pos_200)
%{
tubeHeight: Inject glass tube height information to a quantified results
INPUTS:
    dataSet: quantified data structure includnig "info", "cycle", "timeSeries", etc...
    pos_100: vector for positions of 100 um height glass tubes
    pos_200: vector for positions of 200 um height glass tubes
OUTPUT:
    dataSet: added "TUBE_HEIGHT" column dataSet

EXAMPLE:
cntrl_noSperm = 0:2;
qc = postprocessing.tubeHeight(qc, cntrl_noSperm, []);
%}
dataSet.info.TUBE_HEIGHT = NaN(length(dataSet.info.POS_ID),1);
dataSet.timeSeries.TUBE_HEIGHT = NaN(length(dataSet.timeSeries.POS_ID),1);

positions = unique(dataSet.info.POS_ID);

for i = 1:length(positions)
    index_info = dataSet.info.POS_ID == positions(i);
    index_ts = dataSet.timeSeries.POS_ID == positions(i);
    if ismember(positions(i), pos_100)
       dataSet.info.TUBE_HEIGHT(index_info) = 100;
       dataSet.timeSeries.TUBE_HEIGHT(index_ts) = 100;
    elseif ismember(positions(i), pos_200)
       dataSet.info.TUBE_HEIGHT(index_info) = 200;
       dataSet.timeSeries.TUBE_HEIGHT(index_ts) = 200;
    end
end
end