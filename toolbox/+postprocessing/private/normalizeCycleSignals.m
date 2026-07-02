% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function [normCycles, shapeCycles] = normalizeCycleSignals(cycles, phase, nPhasePoints)

    nCyclesTotal = numel(cycles.total);
    names = fieldnames(cycles);

    for i = 1:numel(names)
        name = names{i};
        normCycles.(name) = nan(nCyclesTotal, nPhasePoints);
        for c = 1:nCyclesTotal
            normCycles.(name)(c, :) = phaseResampleAllowNaN(cycles.(name){c}, phase);
        end
    end

    shapeCycles.total = nan(nCyclesTotal, nPhasePoints);
    shapeCycles.nuc = nan(nCyclesTotal, nPhasePoints);
    shapeCycles.cyto = nan(nCyclesTotal, nPhasePoints);

    for c = 1:nCyclesTotal
        shapeCycles.total(c, :) = amplitudeNormalize(normCycles.total(c, :));
        shapeCycles.nuc(c, :) = amplitudeNormalize(normCycles.nuc(c, :));
        shapeCycles.cyto(c, :) = amplitudeNormalize(normCycles.cyto(c, :));
    end
end

%% =========================
%  Local functions: feature extraction
%  =========================

