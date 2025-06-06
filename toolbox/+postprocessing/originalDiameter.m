function D = originalDiameter(d, h)
    % originalDiameter  Estimate initial sphere diameter
    %   Skips elements where d < h by returning NaN for those indices.
    %
    %   D = originalDiameter(d, h)
    %
    %   d, h : equal-sized vectors (or scalars)
    %   D    : same size, NaN where d < h
    
    mask = (d < h);                 % keep only plausible cases
    D = nthroot(d.^2 .* h, 3);
    D(mask) = d(mask);
end