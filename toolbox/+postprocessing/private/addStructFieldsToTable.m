% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function T = addStructFieldsToTable(T, S)
    names = fieldnames(S);
    for i = 1:numel(names)
        T.(names{i}) = S.(names{i});
    end
end

