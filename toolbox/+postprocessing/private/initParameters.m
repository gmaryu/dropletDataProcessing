% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function P = initParameters()

    %% Smoothing
    P.detectSmoothWindow = 1;   % retained for compatibility; no independent peak detection is done
    P.featureSmoothWindow = 1;  % 1 = no smoothing for feature extraction

    %% Cycle QC
    P.minCycleLength = 5;
    P.maxCycleLength = inf;
    P.minAmplitude = 0.05;
    P.minTroughPhase = 0.05;
    P.maxTroughPhase = 0.95;
    P.maxNaNFractionTotal = 0.20;

    %% Peak-to-peak trough-based windows
    P.interphaseOffsetAfterTrough = 35;
    P.lateOffsetAfterTrough = 35; % retained for compatibility/documentation
    P.minWindowPoints = 3;

    %% Interphase end detection from NLS-defined nuclear volume
    P.useNucVolForInterphaseEnd = true;
    P.nucVolLossThreshold = 0;
    P.useRelativeNucVolLossThreshold = false;
    P.relativeNucVolLossThreshold = 0.05;
    P.nucVolLossConsecutiveN = 2;
    P.minFramesAfterTroughForLossSearch = 1;
    P.excludeLastNPhasePointsFromInterphase = 3;

    %% State 1: cytoplasmic-like waveform
    P.useCorrForState1 = true;
    P.useRMSEForState1 = true;
    P.state1UseOrLogic = true;     % true: corr OR RMSE, false: corr AND RMSE
    P.state1CorrMultiplier = 1.5;
    P.state1RMSEMultiplier = 1.5;
    P.useMedianCytoTemplateForDistance = false;

    %% State 2: high or gradually increasing interphase nuclear Cdk1
    P.nucRatioHighThreshold = 1.10;
    P.useDataDrivenNucRatioThreshold = false;
    P.useP90ForState2 = true;
    P.useLateInterphaseMeanForState2 = true;
    P.useIncreaseForState2 = true;
    P.useNucMinusCytoForState2 = false;
    P.nucRatioIncreaseThreshold = 0.05;
    P.nucMinusCytoHighThreshold = 0.00;
    P.interphasePercentileForState2 = 90; % historical name kept as P90 output column for compatibility

    %% State 3: late nuclear activation + context gate
    % State 3 is defined as:
    %   late nuclear activation
    %   AND (high N/C volume ratio OR within-droplet long cycle)
    % A global/cytoplasm-reference long-cycle gate is intentionally not used.
    P.lateActivationThreshold = 0.25;
    P.lateTimingOffsetFromTrough = 0.20;  % phase units
    P.useState3ContextGate = true;
    P.useDataDrivenState3ContextThresholds = true;
    P.highNCvolRatioThreshold = 0.05;      % if NaN and data-driven fails, N/C criterion is ignored
    P.state3ContextMADMultiplier = 3;
    P.useInterphaseNCvolRatioForState3 = true;

    %% Within-droplet long-cycle context for State 3
    P.withinDropletFoldChangeThreshold = 1.50;
    P.withinDropletDeltaThreshold = 20;
    P.withinDropletZThreshold = 5.0;

    %% Optional state-sequence post-processing
    P.applyConsecutiveSmoothing = false;
    P.consecutiveN = 2;

    %% Output
end

