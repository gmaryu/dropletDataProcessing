function dataSet = predReplicationLimit(dataSet)

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
        k = log2(1+2*(tmpDia/D).^3*((2.^N)-1));

        info.rmax(i) = k;
        info.rmax_DNACR(i) = k*DNACR1;
    else
        info.rmax(i) = nan;
        info.rmax_DNACR(i) = nan;
    end
end
dataSet.info = info;
end