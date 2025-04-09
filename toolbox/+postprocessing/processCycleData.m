function [tp_updated, cycleMetrics] = processCycleData(tp, tm, frameToMin, pixelToUm, spermCondition)
% processCycleData Process cycle metrics for a single droplet.
%
%   [tp_updated, cycleMetrics] = processCycleData(tp, tm, frameToMin, pixelToUm, spermCondition)
%
% This function processes cycle data defined in the table tp (each row representing a cycle)
% by extracting a segment of the droplet tracking data (tm) for each cycle. It computes:
%   - INTERPHASE_START_FRAME and INTERPHASE_END_FRAME,
%   - AREA_NPIXELS_MEDIAN: median droplet area (converted to µm²),
%   - NUC_NPIXELS_Q90: 90th percentile of nuclear pixel counts,
%   - DNA_SUM_INT_Q90: 90th percentile of Hoechst intensity,
%   - DNA_NPIXELS_Q90: 90th percentile of DNA pixel counts,
%   - DNA_INC_RATE_COEFF: the coefficient (log10) from a linear fit of Hoechst intensity.
%
% Inputs:
%   tp             - Table with cycle definitions (must contain fields START_INDEX, END_INDEX, FRAME, etc.).
%   tm             - Table with droplet tracking data (must contain FRAME, AREA, and if spermCondition is true:
%                    NPIXEL_NUC, SUMINTENSITY_DNA, NPIXEL_DNA).
%   frameToMin     - Scalar conversion factor to convert frame numbers to minutes.
%   pixelToUm      - Scalar conversion factor to convert pixels to micrometers.
%   spermCondition - Logical flag indicating whether nuclear/DNA data are present.
%
% Outputs:
%   tp_updated   - The updated cycle table with new computed columns added.
%   cycleMetrics - A structure containing the computed vectors for each metric.
%
% Example:
%   [tp_updated, metrics] = postprocessing.processCycleData(tp, tm, 6, 2.649, true);

    arguments
        tp table
        tm table
        frameToMin (1,1) double {mustBePositive}
        pixelToUm (1,1) double {mustBePositive}
        spermCondition logical
    end
    
    nCycles = height(tp);
    
    interphaseStartFrame = nan(nCycles, 1);
    interphaseEndFrame   = nan(nCycles, 1);
    areaPixelsMedian     = nan(nCycles, 1);
    nucPixelsQ90         = nan(nCycles, 1);
    dnaSumIntQ90         = nan(nCycles, 1);
    dnaPixelsQ90         = nan(nCycles, 1);
    dnaIncRateCoeff      = nan(nCycles, 1);
    
    for k = 1:nCycles
        startidx = tp.START_INDEX(k);
        endidx = tp.END_INDEX(k);
        cycleData = tm(startidx:endidx-1, :);
        
        if spermCondition && ismember('NPIXEL_NUC', cycleData.Properties.VariableNames)
            intstart = find(cycleData.NPIXEL_NUC > 0, 1);
            revIdx = find(flipud(cycleData.NPIXEL_NUC > 0), 1);
            if ~isempty(revIdx)
                intend = height(cycleData) - revIdx + 1;
            else
                intend = [];
            end
        else
            intstart = [];
            intend = [];
        end
        
        % Compute median area (convert to µm²).
        areaPixelsMedian(k) = median(cycleData.AREA / pixelToUm^2);
        
        % Compute time vector in minutes.
        t = frameToMin * cycleData.FRAME;
        
        if spermCondition && ismember('SUMINTENSITY_DNA', cycleData.Properties.VariableNames)
            hoechst = cycleData.SUMINTENSITY_DNA;
            valid = ~isnan(hoechst);
            t = t(valid);
            hoechst = hoechst(valid);
        else
            hoechst = zeros(size(cycleData.FRAME));
        end
        
        % Compute DNA increase rate via linear fitting if enough data points exist.
        if numel(t) > 5
            p = polyfit(t, hoechst, 1);
            if p(1) > 0
                dnaIncRateCoeff(k) = log10(p(1));
            else
                dnaIncRateCoeff(k) = nan;
            end
        else
            dnaIncRateCoeff(k) = nan;
        end
        
        if isempty(intstart)
            interphaseStartFrame(k) = nan;
            interphaseEndFrame(k)   = nan;
            nucPixelsQ90(k)         = nan;
            dnaSumIntQ90(k)         = nan;
            dnaPixelsQ90(k)         = nan;
        else
            interphaseStartFrame(k) = cycleData.FRAME(intstart);
            intData = cycleData(cycleData.FRAME >= interphaseStartFrame(k), :);
            nucPixelsQ90(k) = quantile(intData.NPIXEL_NUC, 0.9);
            dnaSumIntQ90(k) = quantile(intData.SUMINTENSITY_DNA, 0.9);
            dnaPixelsQ90(k) = quantile(intData.NPIXEL_DNA, 0.9);
            if ~isempty(intend) && intend > intstart
                interphaseEndFrame(k) = cycleData.FRAME(intend);
            else
                interphaseEndFrame(k) = nan;
            end
        end
    end
    
    % Update the cycle table with computed metrics.
    tp.INTERPHASE_START_FRAME = interphaseStartFrame;
    tp.INTERPHASE_END_FRAME = interphaseEndFrame;
    tp.AREA_NPIXELS_MEDIAN = areaPixelsMedian;
    tp.NUC_NPIXELS_Q90 = nucPixelsQ90;
    tp.DNA_SUM_INT_Q90 = dnaSumIntQ90;
    tp.DNA_NPIXELS_Q90 = dnaPixelsQ90;
    tp.DNA_INC_RATE_COEFF = dnaIncRateCoeff;
    
    tp_updated = tp;
    
    % Also return cycleMetrics in a structure.
    cycleMetrics.interphaseStartFrame = interphaseStartFrame;
    cycleMetrics.interphaseEndFrame = interphaseEndFrame;
    cycleMetrics.areaPixelsMedian = areaPixelsMedian;
    cycleMetrics.nucPixelsQ90 = nucPixelsQ90;
    cycleMetrics.dnaSumIntQ90 = dnaSumIntQ90;
    cycleMetrics.dnaPixelsQ90 = dnaPixelsQ90;
    cycleMetrics.dnaIncRateCoeff = dnaIncRateCoeff;
end
