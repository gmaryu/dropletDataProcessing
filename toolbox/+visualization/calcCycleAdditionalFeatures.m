function dataSet = calcCycleAdditionalFeatures(dataSet)
dnarenormfactor = 1; % 1e7 factor to make numbers not too big

% --- make IGNORE column if necessary ---
if ~ismember("IGNORED", dataSet.info.Properties.VariableNames)
    disp("IGNORED is added to dataSet.Info");
    dataSet.info.IGNORED = zeros(size(dataSet.info,1),1);
end

if ~ismember("IGNORED", dataSet.cycle.Properties.VariableNames)
    disp("IGNORED is added to dataSet.Info");
    dataSet.cycle.IGNORED = zeros(size(dataSet.cycle,1),1);
end

if ~ismember("IGNORED", dataSet.cycle.Properties.VariableNames)
    disp("IGNORED is added to dataSet.Info");
    dataSet.timeSeries.IGNORED = zeros(size(dataSet.timeSeries,1),1);
end

% --- extract variables ---
tp = dataSet.cycle;
tm = dataSet.timeSeries;
FrameToMin = dataSet.FrameToMin;
PixelToUm = dataSet.PixelToUm;

%% Calculate new features related to nucleus if nuclear mask data exists 
if any(~isnan(tp.("NUC_NPIXELS_Q90")))

    % --- check if a new column has to be calculated ----
    flag_rateCoff = 0;
    if ~ismember('NUC_INC_RATE_COEFF',tp.Properties.VariableNames)
        flag_rateCoff = 1;
        disp('Added a new field "NUC_INC_RATE_COEFF"');
    end

    flag_NUC_NPIXELS_MOD_Q90 = 0;
    if ~ismember('NUC_NPIXELS_MOD_Q90',tp.Properties.VariableNames)
        flag_NUC_NPIXELS_MOD_Q90 = 1;
        disp('Added a new field "NUC_NPIXELS_MOD_Q90"');
    end

    flag_DNA_SUM_INT_MOD_Q90 = 0;
    if ~ismember('DNA_SUM_INT_MOD_Q90',tp.Properties.VariableNames)
        flag_DNA_SUM_INT_MOD_Q90 = 1;
        disp('Added a new field "DNA_SUM_INT_MOD_Q90"');
    end

    flag_NSURF_Q90 = 0;
    if ~ismember('NSURF_Q90',tp.Properties.VariableNames)
        flag_NSURF_Q90 = 1;
        disp('Added a new field "NSURF_Q90"');
    end

    flag_NUCVOL_Q90 = 0;
    if ~ismember('NUC_VOLUMEUM3_Q90',tp.Properties.VariableNames)
        flag_NUCVOL_Q90 = 1;
        disp('Added a new field "NUC_VOLUMEUM3_Q90"');
    end


    % --- number of loop and data allocation ---
    nCycles              = height(tp); % number of peaks
    
    nucAreaIncRateCoeff  = nan(nCycles, 1);
    nucNpixelMODQ90      = nan(nCycles, 1);
    nucDNAintMODQ90      = nan(nCycles, 1);
    nucSurfQ90           = nan(nCycles, 1);
    nucVolQ90            = nan(nCycles, 1);

    for k = 1:nCycles
        if ~tp.IGNORED(k)
            % specify droplet
            posid = tp.POS_ID(k);
            trackid = tp.TRACK_ID(k);

            % time course data for the target droplet
            tmp_tm = tm(tm.POS_ID==posid & tm.TRACK_ID==trackid,:);

            % define START and END time. Use whole cycle if there is no nuclear mask.
            if ~isnan(tp.INTERPHASE_START_FRAME(k))
                startidx = find(tmp_tm.FRAME == tp.INTERPHASE_START_FRAME(k));
            else
                startidx = find(tmp_tm.FRAME == tp.START_FRAME(k));
            end
            if ~isnan(tp.INTERPHASE_END_FRAME(k))
                endidx = find(tmp_tm.FRAME == tp.INTERPHASE_END_FRAME(k));
            else
                endidx = find(tmp_tm.FRAME == tp.END_FRAME(k));
            end

            % extract data for analysis
            cycle_time = startidx:endidx;
            cycleData = tmp_tm(cycle_time, :);

            % --- Calc stats ---
            if flag_rateCoff
                % Compute Nuclear Area increase rate via linear fitting if enough data points exist.
                p = polyfit(cycleData.FRAME, cycleData.NPIXEL_NUC_MOD, 1);
                if p(1) > 0
                    nucAreaIncRateCoeff(k) = p(1);
                else
                    nucAreaIncRateCoeff(k) = nan;
                end
            else
                nucAreaIncRateCoeff(k) = nan;
            end

            if flag_NUC_NPIXELS_MOD_Q90
                nucNpixelMODQ90(k) = quantile(cycleData.NPIXEL_NUC_MOD, 0.9);
            else
                nucNpixelMODQ90(k) = nan;
            end

            if flag_DNA_SUM_INT_MOD_Q90
                nucDNAintMODQ90(k) = quantile(cycleData.SUM_NUCLEUS_HORCHST_INT_MOD, 0.9);
            else
                nucDNAintMODQ90(k) = nan;
            end

            if flag_NSURF_Q90
                nucSurfQ90(k) = quantile(cycleData.NUC_SURF_AREA, 0.9);
            else
                nucSurfQ90(k) = nan;
            end

            if flag_NUCVOL_Q90
                nucVolQ90(k) = quantile(cycleData.NUC_VOLUMEUM3, 0.9);
            else
                nucVolQ90(k) = nan;
            end
        end

    end
    
    % --- add column to the cycle information table ---
    if flag_rateCoff
        tp.NUC_INC_RATE_COEFF = nucAreaIncRateCoeff;
    end
    if flag_NUC_NPIXELS_MOD_Q90
        tp.NUC_NPIXELS_MOD_Q90 = nucNpixelMODQ90;
    end
    if flag_DNA_SUM_INT_MOD_Q90
        tp.DNA_SUM_INT_MOD_Q90 = nucDNAintMODQ90;
    end
    if flag_NSURF_Q90
        tp.NSURF_Q90 = nucSurfQ90;
    end
    if flag_NUCVOL_Q90
        tp.NUC_VOLUMEUM3_Q90 = nucVolQ90;
    end

end

tp.START_MINUTE = tp.START_FRAME * FrameToMin; % minute
tp.DURATION = (tp.END_FRAME - tp.START_FRAME) * FrameToMin; % minute
tp.DNA_INC_RATE_COEFF = tp.DNA_INC_RATE_COEFF / log10(dnarenormfactor); % a.u. (sum px) / minute
if  ismember('NUC_NPIXELS_Q90',tp.Properties.VariableNames) && any(~isnan(tp.NUC_NPIXELS_Q90))
    tp.NCVR_ORI = power(tp.NUC_NPIXELS_Q90 ./ tp.AREA_NPIXELS_MEDIAN, 3 / 2); % n.d.
    tp.NCVR = power(tp.NUC_NPIXELS_MOD_Q90 ./ tp.AREA_NPIXELS_MEDIAN, 3 / 2); % n.d.
    tp.DNACR_ORI = tp.DNA_SUM_INT_Q90 ./ visualization.convertAreaPixelsToVolume(tp.AREA_NPIXELS_MEDIAN, PixelToUm) / dnarenormfactor; % a.u (sum px) / um^3
    tp.DNACR = tp.DNA_SUM_INT_MOD_Q90 ./ visualization.convertAreaPixelsToVolume(tp.AREA_NPIXELS_MEDIAN, PixelToUm) / dnarenormfactor; % a.u (sum px) / um^3
    tp.FC_DNA = visualization.foldChangeDNA(tp);
end
tp.VOLUMEUM3 = visualization.convertAreaPixelsToVolume(tp.AREA_NPIXELS_MEDIAN, PixelToUm);
tp.MARKERSIZE = (log10(visualization.convertAreaPixelsToVolume(tp.AREA_NPIXELS_MEDIAN, PixelToUm)) - 5.25) * 3; % log um3 volume roughly within 5-7
tp.MARKERSIZE(tp.MARKERSIZE < 0.5, :) = 0.5;
dataSet.cycle = tp;
end