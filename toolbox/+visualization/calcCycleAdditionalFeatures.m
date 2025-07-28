function retv = calcCycleAdditionalFeatures(cycle, timeseries, FrameToMin, PixelToUm)
dnarenormfactor = 1e7;
% 1e7 factor to make numbers not to big

retv = cycle;
tp = cycle;
tm = timeseries;

if ~isfield(retv, 'NUC_INC_RATE_COEFF')
    disp('Added a new field "NUC_INC_RATE_COEFF"');
    % Collect sum of nuclear area
    nCycles = height(tp);
    nucAreaIncRateCoeff  = nan(nCycles, 1);

    for k = 1:nCycles
        posid = tp.POS_ID(k);
        trackid = tp.TRACK_ID(k);
                
        tmp_tm = tm(bitand(tm.POS_ID==posid, tm.TRACK_ID==trackid),:);
        startidx = find(tmp_tm.FRAME == tp.INTERPHASE_START_FRAME(k));
        endidx = find(tmp_tm.FRAME == tp.INTERPHASE_END_FRAME(k));
        tmp_time = startidx:endidx;

        if k == 1651
            disp(k);
        end
        
        if ~isnan(tmp_time)
            cycleData = tmp_tm(tmp_time, :);

            % Compute time vector in minutes.
            %t = FrameToMin * cycleData.FRAME;
            t = tmp_time;

            % Collect sum of nuclear area
            nucAreas = cycleData.NPIXEL_NUC;
            %nucAreas = cycleData.NPIXEL_NUC / PixelToUm^2;
            
            %valid = ~isnan(nucAreas);
            %t = t(valid);
            %nucAreas = nucAreas(valid);

            % Compute Nuclear Area increase rate via linear fitting if enough
            % data points exist.
            if numel(nucAreas) > 4
                p = polyfit(cycleData.FRAME, nucAreas, 1);
                if p(1) > 0
                    nucAreaIncRateCoeff(k) = p(1);
                else
                    nucAreaIncRateCoeff(k) = nan;
                end
            else
                nucAreaIncRateCoeff(k) = nan;
            end
        else
           nucAreaIncRateCoeff(k) = nan;
        end
    end
    retv.NUC_INC_RATE_COEFF = nucAreaIncRateCoeff;
end


retv.START_MINUTE = retv.START_FRAME * FrameToMin; % minute
retv.DURATION = (retv.END_FRAME - retv.START_FRAME) * FrameToMin; % minute
retv.DNA_INC_RATE_COEFF = retv.DNA_INC_RATE_COEFF / log10(dnarenormfactor); % a.u. (sum px) / minute
retv.NCVR = power(retv.NUC_NPIXELS_Q90 ./ retv.AREA_NPIXELS_MEDIAN, 3 / 2); % n.d.
retv.DNACR = retv.DNA_SUM_INT_Q90 ./ visualization.convertAreaPixelsToVolume(retv.AREA_NPIXELS_MEDIAN, PixelToUm) / dnarenormfactor; % a.u (sum px) / um^3

retv.VOLUMEUM3 = visualization.convertAreaPixelsToVolume(retv.AREA_NPIXELS_MEDIAN, PixelToUm);
retv.MARKERSIZE = (log10(visualization.convertAreaPixelsToVolume(retv.AREA_NPIXELS_MEDIAN, PixelToUm)) - 5.25) * 3; % log um3 volume roughly within 5-7
retv.MARKERSIZE(retv.MARKERSIZE < 0.5, :) = 0.5;
end