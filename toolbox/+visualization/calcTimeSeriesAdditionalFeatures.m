function dataSet = calcTimeSeriesAdditionalFeatures(dataSet)
dnarenormfactor = 1e7; % 1e7 factor to make numbers not too big

tm = dataSet.timeSeries;


%%
tm.START_MINUTE = tm.FRAME * dataSet.FrameToMin; % minute
tm.VOLUMEUM3 = visualization.convertAreaPixelsToVolume(tm.AREA, dataSet.PixelToUm);
tm.NCVR = power(tm.NPIXEL_NUC./ tm.AREA, 3/2); % n.d.
tm.DNACR = tm.SUM_NUCLEUS_HOECHST_INT ./ visualization.convertAreaPixelsToVolume(tm.AREA, dataSet.PixelToUm) / dnarenormfactor; % a.u (sum px) / um^3
%%
dataSet.timeSeries = tm;
end