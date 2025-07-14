function retv = calcCycleAdditionalFeatures(cycle, FrameToMin, PixelToUm)
dnarenormfactor = 1e7;
% 1e7 factor to make numbers not to big

retv = cycle;

retv.START_MINUTE = retv.START_FRAME * FrameToMin; % minute
retv.DURATION = (retv.END_FRAME - retv.START_FRAME) * FrameToMin; % minute
retv.DNA_INC_RATE_COEFF = retv.DNA_INC_RATE_COEFF / log10(dnarenormfactor); % a.u. (sum px) / minute
retv.NCVR = power(retv.NUC_NPIXELS_Q90 ./ retv.AREA_NPIXELS_MEDIAN, 3 / 2); % n.d.
retv.DNACR = retv.DNA_SUM_INT_Q90 ./ visualization.convertAreaPixelsToVolume(retv.AREA_NPIXELS_MEDIAN, PixelToUm) / dnarenormfactor; % a.u (sum px) / um^3

retv.VOLUMEUM3 = visualization.convertAreaPixelsToVolume(retv.AREA_NPIXELS_MEDIAN, PixelToUm);
retv.MARKERSIZE = (log10(visualization.convertAreaPixelsToVolume(retv.AREA_NPIXELS_MEDIAN, PixelToUm)) - 5.25) * 3; % log um3 volume roughly within 5-7
retv.MARKERSIZE(retv.MARKERSIZE < 0.5, :) = 0.5;
end