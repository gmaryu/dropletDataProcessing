function retv = normalizedCycleAUC(t, s)
    % input
    % s: time series fragment for a neighboring peak-to-peak window (1xN or Nx1)
    % t: corresponding time (same dimensions)

    tn = t(~isnan(s));
    sn = s(~isnan(s));

    % min max normalization
    mm = max(sn);
    mn = min(sn);

    sn = (sn - mn) / (mm - mn);

    % time rescaling (into 0-1 window)
    tn = tn - tn(1);
    tn = tn / tn(end);

    % factor 0.5 is taken from int_0^1 \cos^2(\pi x) = int_0^1 (1 + \cos(2\pi x)) / 2 = 0.5
    % number close to 1 - sinusoid
    % lower AUC - potentially relaxation oscillator-like waveforms
    retv = trapz(tn, sn) / 0.5;
end