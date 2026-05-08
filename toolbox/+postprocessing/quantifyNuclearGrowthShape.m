function result = quantifyNuclearGrowthShape(nucVol, nucleusDetected, cycleStartIdx, cycleEndIdx, frameToMin)
% quantifyNuclearGrowthShape
%
% Quantifies nuclear growth shape within each Cdk1-defined cell cycle.
%
% INPUT
%   NCratio         : vector of N/C ratio
%   nucleusDetected : logical vector, true when nuclei are detected in NLS image
%   cycleStartIdx   : vector of cycle start indices
%   cycleEndIdx     : vector of cycle end indices
%   frameToMin      : minutes per frame
%
% OUTPUT
%   result : table with slopes, plateau index, and growth type

nucVol = nucVol(:);
nucleusDetected = nucleusDetected(:);
t = (1:length(nucVol))';

% Parameters
smoothWindow = 3;
minFitPoints = 5;
earlyFraction = 0.15;
lateFraction = 0.70;
minSegmentLength = 5;
minPiecewiseImprovement = 0.25;
minRobustPoints = 5;


% Slope threshold for plateau classification
% This should be tuned using control data.
minSlopePerMin = 1e-4;

NucSmooth = smoothdata(nucVol, 'movmedian', smoothWindow);

rows = [];

for c = 1:length(cycleStartIdx)

    idxCycle = cycleStartIdx(c):cycleEndIdx(c);

    % Use only frames where nuclei were detected
    idxNuc = idxCycle(nucleusDetected(idxCycle));

    if length(idxNuc) < minFitPoints
        continue
    end

    % Optional: keep the longest continuous detected segment
    idxNuc = keepLongestContinuousSegment(idxNuc);

    if length(idxNuc) < minFitPoints
        continue
    end

    tNuc = t(idxNuc);
    yRaw = nucVol(idxNuc);
    ySmooth = NucSmooth(idxNuc);

    % Normalize time within this nuclear period
    tau = (tNuc - tNuc(1)) * frameToMin;
    durationMin = tau(end) - tau(1);

    if durationMin <= 0
        continue
    end

    % -------------------------
    % Overall robust slope
    % -------------------------
    overallSlope = safeSlopeFit(tau, yRaw, minRobustPoints);

    % -------------------------
    % Early slope
    % -------------------------
    nPts = length(idxNuc);

    if nPts < 4
        growthType = "tooShort";
        earlySlope = NaN;
        lateSlope = NaN;
        overallSlope = safeSlopeFit(tau, yRaw, 5);
    end

    nEarly = round(earlyFraction * nPts);
    nLate  = round(lateFraction  * nPts);

    nEarly = max(2, nEarly);
    nLate  = max(2, nLate);

    nEarly = min(nEarly, floor(nPts / 2));
    nLate  = min(nLate,  floor(nPts / 2));

    earlyIdxLocal = 1:nEarly;
    lateIdxLocal  = (nPts - nLate + 1):nPts;

    if nEarly >= nPts
        nEarly = floor(nPts/2);
    end
    if nLate >= nPts
        nLate = floor(nPts/2);
    end

    earlyIdxLocal = 1:nEarly;
    lateIdxLocal = (nPts - nLate + 1):nPts;


    earlySlope = safeSlopeFit(tau(earlyIdxLocal), yRaw(earlyIdxLocal), minRobustPoints);
    lateSlope  = safeSlopeFit(tau(lateIdxLocal),  yRaw(lateIdxLocal),  minRobustPoints);

    % -------------------------
    % Plateau index
    % -------------------------
    earlySlopePos = max(earlySlope, eps);
    lateSlopePos = max(lateSlope, 0);

    plateauIndex = 1 - lateSlopePos / earlySlopePos;

    % Keep within a reasonable range
    plateauIndex = max(0, min(1, plateauIndex));

    % -------------------------
    % Growth fraction
    % -------------------------
    dy = diff(ySmooth) / frameToMin;

    % A relative threshold based on this cycle's amplitude
    amp = max(ySmooth) - min(ySmooth);
    derivativeThreshold = 0.05 * amp / durationMin;

    growthFraction = mean(dy > derivativeThreshold);

    % -------------------------
    % Piecewise linear fit
    % -------------------------
    if nPts < 10
        piecewiseImprovement = NaN;
        breakIdx = NaN;
        slope1 = NaN;
        slope2 = NaN;
    else
        [piecewiseImprovement, breakIdx, slope1, slope2] = ...
            bestTwoPhaseFit(tau, yRaw, idxNuc, minSegmentLength);
    end

    % -------------------------
    % Shape metrics for NucVolume
    % -------------------------

    yStart = ySmooth(1);
    yEnd   = ySmooth(end);
    yMax   = max(ySmooth);
    yMin   = min(ySmooth);

    totalGrowth = yMax - yStart;
    dynamicRange = yMax - yMin;

    % Avoid division by zero
    if totalGrowth <= 0 || dynamicRange <= 0
        growthType = "noGrowth";
    else
        nPts = length(ySmooth);

        earlyEndLocal = max(minFitPoints, round(0.30 * nPts));
        midLocal      = max(minFitPoints, round(0.50 * nPts));

        earlyEndValue = mean(ySmooth(max(1, earlyEndLocal-1):min(nPts, earlyEndLocal+1)));
        midValue      = mean(ySmooth(max(1, midLocal-1):min(nPts, midLocal+1)));
        endValue      = mean(ySmooth(max(1, nPts-2):nPts));

        earlyGain = earlyEndValue - yStart;
        lateGain  = endValue - midValue;

        earlyGainFraction = earlyGain / totalGrowth;
        lateGainFraction  = lateGain  / totalGrowth;

        [~, maxLocalIdx] = max(ySmooth);
        timeToMaxFraction = maxLocalIdx / nPts;

        overallGrowthRate = totalGrowth / durationMin;

        if overallGrowthRate <= 0
            lateSlopeNorm = NaN;
            earlySlopeNorm = NaN;
        else
            earlySlopeNorm = earlySlope / overallGrowthRate;
            lateSlopeNorm  = lateSlope  / overallGrowthRate;
        end

        % -------------------------
        % Classification
        % -------------------------

        minClassifyPoints = 8;

        if nPts < minClassifyPoints
            growthType = "tooShort";
        else
            isPlateau = ...
                lateSlopeNorm < -0.15 && ...
                lateGainFraction < 0.20 && ...
                timeToMaxFraction < 0.65;

            isTwoPhase = ...
                piecewiseImprovement > minPiecewiseImprovement && ...
                earlySlopeNorm > 1.5 && ...
                lateSlopeNorm > -0.15 && ...
                lateGainFraction >= 0.20;

            isContinuous = ...
                lateGainFraction >= 0.20 || ...
                timeToMaxFraction >= 0.65;

            if isPlateau
                growthType = "plateau";
            elseif isTwoPhase
                growthType = "twoPhase";
            elseif isContinuous
                growthType = "continuous";
            else
                growthType = "weakGrowth";
            end
        end


        
    end

    % -------------------------
    % Store result
    % -------------------------
    row.cycle = c;
    row.startIdx = idxNuc(1);
    row.endIdx = idxNuc(end);
    row.durationMin = durationMin;

    row.overallSlope = overallSlope;
    row.earlySlope = earlySlope;
    row.lateSlope = lateSlope;
    row.plateauIndex = plateauIndex;
    row.growthFraction = growthFraction;

    row.piecewiseImprovement = piecewiseImprovement;
    row.piecewiseBreakIdx = breakIdx;
    row.piecewiseEarlySlope = slope1;
    row.piecewiseLateSlope = slope2;

    row.maxNC = max(yRaw);
    row.minNC = min(yRaw);
    row.amplitudeNC = max(yRaw) - min(yRaw);

    row.growthType = growthType;

    rows = [rows; row];

end

result = struct2table(rows);
end

function idxOut = keepLongestContinuousSegment(idxIn)
% Keeps the longest consecutive segment from a vector of indices.

if isempty(idxIn)
    idxOut = idxIn;
    return
end

breaks = find(diff(idxIn) > 1);

segStarts = [1; breaks(:) + 1];
segEnds = [breaks(:); length(idxIn)];

segLengths = segEnds - segStarts + 1;

[~, bestSeg] = max(segLengths);

idxOut = idxIn(segStarts(bestSeg):segEnds(bestSeg));
end

function [bestImprovement, bestBreakIdx, bestSlope1, bestSlope2] = bestTwoPhaseFit(t, y, originalIdx, minSegmentLength)
% Finds the best two-phase robust linear fit.

n = length(y);

bestImprovement = NaN;
bestBreakIdx = NaN;
bestSlope1 = NaN;
bestSlope2 = NaN;

if n < 2 * minSegmentLength
    return
end

% Single robust fit
bSingle = robustfit(t, y);
yPredSingle = bSingle(1) + bSingle(2) * t;
RSSsingle = sum((y - yPredSingle).^2);

if RSSsingle <= 0
    return
end

bestRSS = Inf;

for b = minSegmentLength:(n - minSegmentLength)

    t1 = t(1:b);
    y1 = y(1:b);

    t2 = t((b+1):end);
    y2 = y((b+1):end);

    b1 = robustfit(t1, y1);
    b2 = robustfit(t2, y2);

    yPred1 = b1(1) + b1(2) * t1;
    yPred2 = b2(1) + b2(2) * t2;

    RSS = sum((y1 - yPred1).^2) + sum((y2 - yPred2).^2);

    if RSS < bestRSS
        bestRSS = RSS;
        bestBreakIdx = originalIdx(b);
        bestSlope1 = b1(2);
        bestSlope2 = b2(2);
    end
end

bestImprovement = (RSSsingle - bestRSS) / RSSsingle;
end

function slope = safeSlopeFit(t, y, minRobustPoints)
% safeSlopeFit
% Returns slope using robustfit when possible.
% Falls back to polyfit for short segments.
% Returns NaN if fewer than 2 valid points are available.

    if nargin < 3
        minRobustPoints = 5;
    end

    t = t(:);
    y = y(:);

    valid = isfinite(t) & isfinite(y);
    t = t(valid);
    y = y(valid);

    if numel(y) < 2
        slope = NaN;
        return
    end

    % If all time values are identical, slope cannot be estimated
    if numel(unique(t)) < 2
        slope = NaN;
        return
    end

    if numel(y) >= minRobustPoints
        try
            b = robustfit(t, y);
            slope = b(2);
        catch
            p = polyfit(t, y, 1);
            slope = p(1);
        end
    else
        p = polyfit(t, y, 1);
        slope = p(1);
    end
end