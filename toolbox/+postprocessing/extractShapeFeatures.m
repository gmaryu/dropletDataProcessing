function [nucCycleInfo, waveforms] = extractShapeFeatures(nucData, dataSet, traceInfo, refPath)
%EXTRACTSHAPEFEATURES  Preprocessing + segmentation + Cdk1 waveform feature
% extraction only (no state classification). Preprocessing/segmentation/
% feature-extraction local functions live canonically in
% toolbox/+postprocessing/private/ (originally lifted verbatim from
% toolbox/internal/+archived/classifyCdk1dynamics.m, now archived/superseded,
% so features are identical to the original circular classifier's).
%
% Returns nucCycleInfo with shape-only features (CytoDistance*, Interphase*,
% Late*, PeakPhase_nuc, etc.) plus context columns (NCvolRatio, CycleLength)
% that Stage A must NOT use for classification.
%
% Optional 2nd output `waveforms` (struct) returns the phase-resampled cycle
% waveforms aligned row-for-row with nucCycleInfo, for supplementary plotting:
%   .normNuc/.normCyto/.normTotal  = phase-resampled, amplitude PRESERVED
%   .shapeNuc/.shapeCyto/.shapeTotal = phase-resampled AND amplitude-normalized
%   .phase (phase grid), .nPhasePoints

P = initParameters();
[cytoRef, phase, nPhasePoints] = loadCytoplasmicReference(refPath);
raw = getNuclearRawMatrices(nucData);
validateRawMatrices(raw);
[nDroplets, nTime] = size(raw.total);
[~, cycleBoundaryTable, traceInfo] = prepareCycleBoundaryTable(dataSet, traceInfo, nDroplets); %#ok<ASGLU>
featuresRaw = smoothAndDeriveSignals(raw, P);
[nucCycleInfo, cycles] = segmentPeakToPeakCycles(cycleBoundaryTable, featuresRaw, nDroplets, nTime);
nucCycleInfo = applyCycleQC(nucCycleInfo, P);
[normCycles, shapeCycles] = normalizeCycleSignals(cycles, phase, nPhasePoints);
[nucCycleInfo, ~] = extractCycleFeatures(nucCycleInfo, normCycles, shapeCycles, featuresRaw, cytoRef, phase, nPhasePoints, P);
nucCycleInfo = addWithinDropletCycleLengthContext(nucCycleInfo, P);
if nargout > 1
    waveforms = struct( ...
        'normNuc', normCycles.nuc, 'normCyto', normCycles.cyto, 'normTotal', normCycles.total, ...
        'shapeNuc', shapeCycles.nuc, 'shapeCyto', shapeCycles.cyto, 'shapeTotal', shapeCycles.total, ...
        'phase', phase, 'nPhasePoints', nPhasePoints);
end
fprintf('extractShapeFeatures: %d cycles, %d valid\n', height(nucCycleInfo), sum(nucCycleInfo.IsValidCycle));
end
