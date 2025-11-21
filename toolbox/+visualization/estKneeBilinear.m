function estKneeBilinear(x, y)
x = x(:);
y = y(:);

% ==========
% Parameters
% ==========
smoothing = 'sgolay';   % 'none' | 'sgolay' | 'malowess' | 'rlowess'
sgolay_order = 3;       % sgolay poly order
sgolay_frame = 7;       % sgolay frame length (odd)
loess_span   = 0.07;    % fraction (0-1), adjust if needed
do_bootstrap = true;    % set false to skip CI (faster)
B = 500;                % bootstrap iterations

% =========================
% Safety + sorting + log-x
% =========================
assert(all(x>0),'x must be positive for log-axis analysis.');
[x,ord] = sort(x(:));
y = y(:); y = y(ord);
lx = log10(x);

% ==================
% Optional smoothing
% ==================
ys = y;
switch lower(smoothing)
    case 'none'
        % leave as is
    case 'sgolay'
        ys = sgolayfilt(y, sgolay_order, sgolay_frame);
    case 'malowess'
        ys = smooth(lx, y, loess_span, 'loess');   % distance-weighted local poly
    case 'rlowess'
        ys = smooth(lx, y, loess_span, 'rlowess'); % robust loess
    otherwise
        error('Unknown smoothing option.');
end

% =========================
% Method 1: Bilinear (cont.)
% =========================
out1 = knee_bilinear_continuous(lx, ys);   % lx = log10(x)

% Optional bootstrap CI
if do_bootstrap
    ci1 = knee_bootstrap_CI(lx, ys, B);
else
    ci1 = [NaN NaN];
end

% ==============================
% Method 2: Slope change-point
% ==============================
out2 = knee_by_slope(lx, ys);

% =========
% Reporting
% =========
fprintf('--- Smoothing: %s ---\n', smoothing);
fprintf('[Bilinear]   knee x  = %.6g (idx=%d)\n', 10.^out1.lx_knee, out1.idx_knee);
fprintf('             95%% CI  = [%.6g, %.6g]\n', 10.^ci1(1), 10.^ci1(2));
fprintf('[Slope-CP]   knee x  = %.6g (idx=%d)\n', 10.^out2.lx_knee, out2.idx_knee);

% =====
% Plot
% =====
figure('Color','w'); 
semilogx(x, y, '.', 'MarkerSize', 10); hold on;
% Bilinear fit line
semilogx(x, out1.yhat, '-', 'LineWidth', 1.5);
% Vertical lines for knees
xline(10.^out1.lx_knee, '--', 'Bilinear knee');
%xline(10.^out2.lx_knee, ':', 'Slope knee');

xlabel('NCVR (log scale)'); ylabel('Period (min)');
title(sprintf('Knee detection (%s smoothing)', smoothing));
grid on; legend({'data','bilinear fit','bilinear knee','slope knee'}, 'Location','best');
end
%% ----------------- subfunctions -----------------
function out = knee_bilinear_continuous(lx, y)
% Knee by continuous 2-segment linear fit on lx=log10(x)
% Returns: idx_knee, lx_knee, params, yhat

    n = numel(lx);
    kmin = 2; kmax = n-2;

    best.SSE = inf;
    for k = kmin:kmax
        % Left line: y = a1*lx + b1 (ordinary least squares)
        p1 = polyfit(lx(1:k), y(1:k), 1);
        yk = polyval(p1, lx(k));

        % Right line with continuity at (lx(k), yk):
        % y = a2 * (lx - lx(k)) + yk  => estimate a2 by LS
        dxR = lx(k+1:n) - lx(k);
        dyR = y(k+1:n) - yk;
        denom = sum(dxR.^2);
        if denom <= eps, a2 = 0; else, a2 = sum(dxR .* dyR) / denom; end

        yhat = [ polyval(p1, lx(1:k));
                 a2*(lx(k+1:n) - lx(k)) + yk ];

        SSE = sum((y - yhat).^2);
        if SSE < best.SSE
            best.SSE = SSE;
            best.k   = k;
            best.p1  = p1;
            best.a2  = a2;
            best.yhat= yhat;
        end
    end

    out.idx_knee = best.k;
    out.lx_knee  = lx(best.k);
    out.yhat     = best.yhat;
    out.left     = best.p1;   % [a1 b1]
    out.right_a  = best.a2;   % slope on the right in shifted coords
end

function ci = knee_bootstrap_CI(lx, y, B)
% Simple bootstrap CI for the bilinear knee on lx
    n  = numel(lx);
    ks = zeros(B,1);
    for b = 1:B
        idx = randsample(n, n, true);         % resample pairs
        lxB = lx(idx);
        yB  = y(idx);
        % keep the joint distribution roughly by sorting lx within bootstrap sample
        [lxB, ordB] = sort(lxB); 
        yB = yB(ordB);

        outB = knee_bilinear_continuous(lxB, yB);
        ks(b) = outB.lx_knee;
    end
    ci = quantile(ks, [0.025 0.975]);
end

function out = knee_by_slope(lx, y)
% Knee via change point in slope dy/d(lx)
    dx  = diff(lx);
    dy  = diff(y);
    slope = dy ./ dx;

    % Robustify slope a bit
    win = max(3, round(numel(slope)/50));
    slope = movmedian(slope, win);

    % Use findchangepts if available, otherwise MAD threshold
    try
        k = findchangepts(slope, 'MaxNumChanges', 1, 'Statistic', 'mean');
        idx = k + 1;                   % back to x-index
    catch
        % MAD-based first up-jump
        b = median(slope);
        s = 1.4826 * median(abs(slope - b));
        alpha = 3; % sensitivity
        k = find(slope > b + alpha*s, 1, 'first');
        if isempty(k), k = ceil(numel(slope)/2); end
        idx = k + 1;
    end

    out.idx_knee = idx;
    out.lx_knee  = lx(idx);
    out.slope    = slope;
end





