function dataSet = predReplicationLimit(dataSet)
% Extract dilution rate
d_rate = 0.7;

%dataSet = control_nuc;
info = dataSet.info;
cycle = dataSet.cycle;

info.rmax = nan*ones(size(info,1),1);
info.rmax_DNACR = nan*ones(size(info,1),1);
%%
for i = 1:size(info,1)
    if info.IGNORED(i) == 0
        tmpP = info.POS_ID(i);
        tmpD = info.TRACK_ID(i);
        tmpDia = info.ORIGINAL_MED_DIAMETER(i);

        tmpCycles = cycle(cycle.POS_ID == tmpP & cycle.TRACK_ID == tmpD, :);
        DNACR1 = tmpCycles.DNACR(1);

        D = 1000;
        N = 12;
        m = info.SPERM_COUNT(i);
        k = log2(1+(2/m)*(tmpDia/D).^3*((2.^N)-1)*d_rate);

        info.rmax(i) = k;
        info.rmax_DNACR(i) = (2^k)*DNACR1;
    else
        info.rmax(i) = nan;
        info.rmax_DNACR(i) = nan;
    end
end
dataSet.info = info;
end