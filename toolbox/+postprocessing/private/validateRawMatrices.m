% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function validateRawMatrices(raw)

    baseSize = size(raw.total);
    fields = fieldnames(raw);

    for i = 1:numel(fields)
        f = fields{i};
        assert(isequal(baseSize, size(raw.(f))), ...
            "%s must have the same size as raw.total.", f);
    end
end

