% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function [cytoRef, phase, nPhasePoints] = loadCytoplasmicReference(refPathOrDir)

    % Accept either a full path to the reference .mat, or a folder containing
    % the default-named file. This lets each condition/batch keep its own
    % cytoplasmic reference without overwriting.
    if endsWith(string(refPathOrDir), ".mat", "IgnoreCase", true) && isfile(refPathOrDir)
        refFile = char(refPathOrDir);
    else
        refFile = fullfile(refPathOrDir, "cytoplasm_only_reference_peak2peak.mat");
    end

    load(refFile, ...
        "cytoTemplateByCycle_peak2peak", ...
        "cytoTemplateMedianByCycle_peak2peak", ...
        "cytoThresholdCorrByCycle_peak2peak", ...
        "cytoThresholdRMSEByCycle_peak2peak", ...
        "cytoReferenceSummary_peak2peak", ...
        "phase", ...
        "nPhasePoints");

    cytoRef.templateMean = cytoTemplateByCycle_peak2peak;
    cytoRef.templateMedian = cytoTemplateMedianByCycle_peak2peak;
    cytoRef.thresholdCorr = cytoThresholdCorrByCycle_peak2peak;
    cytoRef.thresholdRMSE = cytoThresholdRMSEByCycle_peak2peak;
    cytoRef.summary = cytoReferenceSummary_peak2peak;

    fprintf("Loaded cytoplasmic reference: %s\n", refFile);
end

