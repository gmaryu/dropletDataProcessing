function result = quantifySingleCdk1Inactivation(signal, varargin)
% quantifySingleCdk1Inactivation Quantify one Cdk1 inactivation trace.
%
% RESULT = quantifySingleCdk1Inactivation(SIGNAL) removes NaNs, treats the
% first valid point as the peak, treats the last valid point as the local
% baseline/trough, normalizes the decline to peak = 1 and trough = 0, then
% fits a decreasing Hill/sigmoid curve:
%
%   y = 1 / (1 + (t / hill_t50)^hill_coefficient)
%
% The Hill coefficient reports the steepness of the inactivation transition.
% The function also reports threshold crossing times, maximum decline rate,
% AUC, and an exponential-fit k/half-life for optional comparison.
%
% RESULT = quantifySingleCdk1Inactivation(SIGNAL, "TimeStep", 5) uses a
% sampling interval of 5, for example 5 minutes between frames.
%
% Example:
%   signal = [1 0.988 0.972 0.940 0.875 0.821 0.803 0.799 0.796 NaN];
%   result = quantifySingleCdk1Inactivation(signal, "TimeStep", 1, ...
%       "MakePlot", true);
%   disp(result.metrics)

opts = parseOptions(varargin{:});

signal = signal(:);
signal = signal(isfinite(signal));

if numel(signal) < 4
    error("At least four finite signal values are required.");
end

time = (0:numel(signal)-1)' .* opts.TimeStep;
peakValue = signal(1);
troughValue = signal(end);
amplitude = peakValue - troughValue;

if amplitude <= 0
    error("The first valid point must be larger than the last valid point.");
end

normalizedSignal = (signal - troughValue) ./ amplitude;
normalizedSignal = min(max(normalizedSignal, 0), 1);

expFit = fitExponentialDecay(time, normalizedSignal, opts);
hillFit = fitHillDecay(time, normalizedSignal, opts);

metrics = table( ...
    hillFit.hill_coefficient, ...
    hillFit.hill_t50, ...
    hillFit.r_squared, ...
    crossingTime(time, normalizedSignal, 0.50), ...
    crossingTime(time, normalizedSignal, 0.25), ...
    crossingTime(time, normalizedSignal, 0.10), ...
    maxDeclineRate(time, normalizedSignal), ...
    trapz(time, normalizedSignal), ...
    expFit.k, ...
    expFit.half_life, ...
    expFit.r_squared, ...
    peakValue, ...
    troughValue, ...
    amplitude, ...
    numel(signal), ...
    'VariableNames', {'hill_coefficient','hill_t50','hill_fit_r_squared', ...
    't50','t25','t10','max_decline_rate','auc','exp_k', ...
    'exp_half_life','exp_fit_r_squared','peak_value', ...
    'trough_value','amplitude','n_points'});

trace = table(time, signal, normalizedSignal, ...
    hillPredict(time, hillFit.hill_t50, hillFit.hill_coefficient), ...
    exp(-expFit.k .* time), ...
    'VariableNames', {'time','raw_signal','normalized_signal', ...
    'hill_fit_signal','exp_fit_signal'});
trace.hill_fit_signal(~isfinite(trace.hill_fit_signal)) = NaN;
trace.exp_fit_signal(~isfinite(trace.exp_fit_signal)) = NaN;

result = struct();
result.metrics = metrics;
result.trace = trace;
result.options = opts;

if opts.MakePlot
    makePlot(result);
end
end

function opts = parseOptions(varargin)
parser = inputParser;
parser.addParameter("TimeStep", 1);
parser.addParameter("FitLowerBound", 0.05);
parser.addParameter("FitUpperBound", 0.95);
parser.addParameter("MinFitPoints", 3);
parser.addParameter("MaxIterations", 10000);
parser.addParameter("MakePlot", false);
parser.parse(varargin{:});
opts = parser.Results;
end

function fitStats = fitExponentialDecay(time, y, opts)
fitStats = struct("k", NaN, "half_life", NaN, "r_squared", NaN);

mask = isfinite(time) & isfinite(y) & ...
    y > opts.FitLowerBound & y < opts.FitUpperBound;

if nnz(mask) < opts.MinFitPoints
    mask = isfinite(time) & isfinite(y) & y > 0 & y < 1;
end

if nnz(mask) < opts.MinFitPoints
    return;
end

x = time(mask);
logY = log(y(mask));
p = polyfit(x, logY, 1);
k = -p(1);

if k <= 0 || ~isfinite(k)
    return;
end

pred = polyval(p, x);
ssRes = sum((logY - pred).^2);
ssTot = sum((logY - mean(logY)).^2);

if ssTot == 0
    rSquared = NaN;
else
    rSquared = 1 - ssRes / ssTot;
end

fitStats.k = k;
fitStats.half_life = log(2) / k;
fitStats.r_squared = rSquared;
end

function fitStats = fitHillDecay(time, y, opts)
fitStats = struct("hill_coefficient", NaN, "hill_t50", NaN, ...
    "r_squared", NaN);

mask = isfinite(time) & isfinite(y) & ...
    y > opts.FitLowerBound & y < opts.FitUpperBound & time > 0;

if nnz(mask) < opts.MinFitPoints
    mask = isfinite(time) & isfinite(y) & y > 0 & y < 1 & time > 0;
end

if nnz(mask) < opts.MinFitPoints
    return;
end

x = time(mask);
yy = y(mask);
initialT50 = crossingTime(time, y, 0.50);
if ~isfinite(initialT50) || initialT50 <= 0
    initialT50 = median(x);
end
initialHill = 2;

initialParams = log([initialT50, initialHill]);
objective = @(logParams) hillObjective(logParams, x, yy);
searchOptions = optimset("Display", "off", ...
    "MaxIter", opts.MaxIterations, "MaxFunEvals", opts.MaxIterations);
bestParams = fminsearch(objective, initialParams, searchOptions);

params = exp(bestParams);
hillT50 = params(1);
hillCoefficient = params(2);
pred = hillPredict(x, hillT50, hillCoefficient);

ssRes = sum((yy - pred).^2);
ssTot = sum((yy - mean(yy)).^2);
if ssTot == 0
    rSquared = NaN;
else
    rSquared = 1 - ssRes / ssTot;
end

fitStats.hill_coefficient = hillCoefficient;
fitStats.hill_t50 = hillT50;
fitStats.r_squared = rSquared;
end

function sse = hillObjective(logParams, time, y)
params = exp(logParams);
pred = hillPredict(time, params(1), params(2));
sse = sum((y - pred).^2);
end

function yHat = hillPredict(time, hillT50, hillCoefficient)
yHat = 1 ./ (1 + (time ./ hillT50) .^ hillCoefficient);
end

function tCross = crossingTime(time, y, level)
tCross = NaN;
idx = find(y <= level, 1, "first");

if isempty(idx)
    return;
end

if idx == 1
    tCross = time(1);
    return;
end

t1 = time(idx - 1);
t2 = time(idx);
y1 = y(idx - 1);
y2 = y(idx);

if y1 == y2
    tCross = t2;
else
    tCross = t1 + (level - y1) * (t2 - t1) / (y2 - y1);
end
end

function rate = maxDeclineRate(time, y)
dy = diff(y);
dt = diff(time);
rate = min(dy ./ dt);
end

function makePlot(result)
figure;
plot(result.trace.time, result.trace.normalized_signal, "ko-", ...
    "LineWidth", 1.2, "MarkerFaceColor", "k");
hold on;
plot(result.trace.time, result.trace.hill_fit_signal, "r-", "LineWidth", 1.5);
plot(result.trace.time, result.trace.exp_fit_signal, "b--", "LineWidth", 1.0);
xlabel("Time");
ylabel("Normalized Cdk1 activity");
ylim([-0.05 1.05]);
legend({"data", "Hill fit", "exp fit"}, "Location", "best");
title("Cdk1 inactivation");
end
