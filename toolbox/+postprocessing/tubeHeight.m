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
for i = 1:length(dataSet.info.POS_ID)
    tmp_p = dataSet.info.POS_ID(i);
    if ismember(tmp_p, pos_100)
       dataSet.info.TUBE_HEIGHT(i) = 100;
    elseif ismember(tmp_p, pos_200)
       dataSet.info.TUBE_HEIGHT(i) = 200;
    end
end

end