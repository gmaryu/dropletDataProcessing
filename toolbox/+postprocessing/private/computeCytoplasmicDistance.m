% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function F = computeCytoplasmicDistance(F, c, k, y, cytoRef, P)

    if isnan(k) || k < 1 || k > size(cytoRef.templateMean, 1)
        return;
    end

    if P.useMedianCytoTemplateForDistance
        template = cytoRef.templateMedian(k, :);
    else
        template = cytoRef.templateMean(k, :);
    end

    if all(isnan(template))
        return;
    end

    r = corr(y(:), template(:), "Rows", "complete");
    F.CytoDistanceCorr(c) = 1 - r; % close to 0, if the dynamics looks like cytoplasmic activation
    F.CytoDistanceRMSE(c) = sqrt(mean((y - template).^2, "omitnan"));
    F.CytoThresholdCorr(c) = cytoRef.thresholdCorr(k);
    F.CytoThresholdRMSE(c) = cytoRef.thresholdRMSE(k);
end

