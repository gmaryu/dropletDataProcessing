function out = classifyPhasesShapeOnly(nucCycleInfo, opts)
%% DEPRECATED / ARCHIVED -- do not call.
% Earlier standalone Stage-A version (GMM + rule + agreement), functionally
% absorbed by postprocessing.runPhaseClassification. Kept here for provenance
% only.
%CLASSIFYPHASESSHAPEONLY  Stage A: phase classification from Cdk1 waveform SHAPE ONLY.
%
% Excludes N/C volume ratio, absolute volumes, cycle length / period, and any
% "long cycle" gate. Those are reserved for Stage B post-hoc characterization.
%
% Runs two independent classifiers and reports their agreement:
%   (A1) Unsupervised GMM (fitgmdist, full covariance, shared-nothing) with
%        automatic k selection by BIC over k = 1..kMax. Discovery: do 3 groups emerge?
%   (A2) N/C-free, period-free rule-based classifier. Confirmation.
%
% opts fields (all optional):
%   .featNames   string array of shape features to use (default below)
%   .kMax        max clusters for BIC scan (default 6)
%   .nReplicates fitgmdist replicates (default 20)
%   .regVal      covariance regularization (default 1e-4)
%   .applyModel  a previously returned out.gmm to APPLY instead of refit
%                (fixes the "k=3 breaks when a phase is absent" failure mode)
%
% out fields:
%   .featNames .X .validIdx .gmm .labelGMM .labelRule .agreement .bic .table

if nargin < 2, opts = struct; end
if ~isfield(opts,'featNames')
    opts.featNames = [ ...
        "CytoDistanceCorr","CytoDistanceRMSE", ...
        "InterphaseNucCytoRatioP90","InterphaseNucCytoRatioSlope", ...
        "InterphaseNucCytoRatioIncrease", ...
        "LateActivationScore_nucShape","TimeToHalfActivation_nuc", ...
        "PeakPhase_nuc","MaxSlopeLate_nucShape"];
end
if ~isfield(opts,'kMax'),        opts.kMax = 6;       end
if ~isfield(opts,'nReplicates'), opts.nReplicates = 20; end
if ~isfield(opts,'regVal'),      opts.regVal = 1e-4;  end

featNames = opts.featNames;
vn = string(nucCycleInfo.Properties.VariableNames);
miss = setdiff(featNames, vn);
assert(isempty(miss), "Missing shape features: %s", strjoin(cellstr(miss), ", "));

% ---- assemble feature matrix over valid, complete-feature cycles ----
nAll = height(nucCycleInfo);
Xraw = nan(nAll, numel(featNames));
for j = 1:numel(featNames)
    Xraw(:,j) = nucCycleInfo.(char(featNames(j)));
end
rowok = nucCycleInfo.IsValidCycle & all(~isnan(Xraw),2);
Xv = Xraw(rowok,:);

% ---- standardize (store transform for apply-mode) ----
if isfield(opts,'applyModel') && ~isempty(opts.applyModel)
    mu_ = opts.applyModel.mean_; sd_ = opts.applyModel.std_;
else
    mu_ = mean(Xv,1); sd_ = std(Xv,0,1); sd_(sd_==0) = 1;
end
X = (Xv - mu_) ./ sd_;

% ---- GMM: fit-or-apply ----
if isfield(opts,'applyModel') && ~isempty(opts.applyModel)
    gmm = opts.applyModel;
    rawlbl = cluster(gmm.model, X);
    labelGMM = gmm.phaseMap(rawlbl);
    bic = [];
else
    [gmm.model, bic, allmodels] = fitGMMbyBIC(X, opts.kMax, opts.nReplicates, opts.regVal); %#ok<ASGLU>
    gmm.kBest = gmm.model.NumComponents;
    rawlbl = cluster(gmm.model, X);
    gmm.phaseMap = orderClustersByShape(X, rawlbl, gmm.model.NumComponents, featNames);
    labelGMM = gmm.phaseMap(rawlbl);
    gmm.mean_ = mu_; gmm.std_ = sd_; gmm.featNames = featNames;
end
labelGMM = labelGMM(:);

% ---- rule-based (shape only, no N/C, no period) ----
labelRule = ruleClassifyShapeOnly(nucCycleInfo(rowok,:));

% ---- agreement ----
agreement = agreementStats(labelGMM, labelRule);

% ---- write back ----
T = nucCycleInfo;
T.PhaseGMM  = nan(nAll,1);
T.PhaseRule = nan(nAll,1);
T.PhaseGMM(rowok)  = labelGMM;
T.PhaseRule(rowok) = labelRule;

out = struct('featNames',{featNames},'X',X,'validIdx',rowok, ...
    'gmm',gmm,'labelGMM',labelGMM,'labelRule',labelRule, ...
    'agreement',agreement,'bic',bic,'table',T);
end

% =====================================================================
function [best, bic, models] = fitGMMbyBIC(X, kMax, nRep, regVal)
bic = inf(kMax,1); models = cell(kMax,1);
opt = statset('MaxIter',500);
for k = 1:kMax
    try
        gm = fitgmdist(X,k,'CovarianceType','full','RegularizationValue',regVal, ...
            'Replicates',nRep,'Options',opt);
        bic(k) = gm.BIC; models{k} = gm;
    catch ME
        warning("k=%d failed: %s", k, ME.message);
    end
end
[~,kBest] = min(bic);
best = models{kBest};
end

function phaseMap = orderClustersByShape(X, rawlbl, k, featNames)
% Order clusters into phases 1..3 by a shape severity score (higher = more
% phase-3 / G2-M-like). Uses standardized feature means; NO N/C, NO period.
fn = cellstr(featNames);
gi = @(nm) find(strcmp(fn,nm));
score = zeros(k,1);
for c = 1:k
    m = mean(X(rawlbl==c,:),1);
    s = 0;
    s = s + getf(m,gi('TimeToHalfActivation_nuc'));
    s = s + getf(m,gi('LateActivationScore_nucShape'));
    s = s + getf(m,gi('MaxSlopeLate_nucShape'));
    s = s + getf(m,gi('PeakPhase_nuc'));
    s = s + getf(m,gi('CytoDistanceCorr'));            % phase1 ~ cyto (low corr distance)
    s = s - getf(m,gi('InterphaseNucCytoRatioP90'));   % phase2 high interphase nuc
    score(c) = s;
end
[~,order] = sort(score,'ascend');
phaseMap = zeros(k,1);
if k <= 3
    for r=1:k, phaseMap(order(r)) = r; end
else
    edges = round(linspace(0,k,4));
    for r=1:k
        for t=1:3
            if r>edges(t) && r<=edges(t+1), phaseMap(order(r))=t; break; end
        end
    end
end
end

function v = getf(vec, idx)
if isempty(idx), v = 0; else, v = vec(idx); end
end

function lbl = ruleClassifyShapeOnly(T)
n = height(T); lbl = nan(n,1);
for c = 1:n
    cytoLike = T.CytoDistanceCorr(c) <= 0.08 ...
        & T.LateActivationScore_nucShape(c) < 0.20 ...
        & T.TimeToHalfActivation_nuc(c) < 0.30;
    g2m = T.LateActivationScore_nucShape(c) >= 0.40 ...
        & T.TimeToHalfActivation_nuc(c) >= 0.60 ...
        & T.InterphaseNucCytoRatioP90(c) < 1.12;
    nuc = T.InterphaseNucCytoRatioP90(c) >= 1.12 ...
        | T.InterphaseNucCytoRatioIncrease(c) >= 0.05;
    if g2m
        lbl(c) = 3;
    elseif nuc
        lbl(c) = 2;
    elseif cytoLike
        lbl(c) = 1;
    else
        lbl(c) = 2;
    end
end
end

function A = agreementStats(a, b)
C = zeros(3);
for i = 1:numel(a)
    if isnan(a(i))||isnan(b(i)), continue; end
    C(a(i),b(i)) = C(a(i),b(i)) + 1;
end
N = sum(C(:));
po = sum(diag(C))/N;
pe = sum(sum(C,2).*sum(C,1)')/N^2;
A.confusion = C;
A.pctAgree = po;
A.kappa = (po-pe)/(1-pe+eps);
A.adjRand = adjustedRand(a,b);
end

function ari = adjustedRand(a,b)
va = unique(a(~isnan(a))); vb = unique(b(~isnan(b)));
C = zeros(numel(va),numel(vb));
for i=1:numel(va)
    for j=1:numel(vb)
        C(i,j) = sum(a==va(i) & b==vb(j));
    end
end
n = sum(C(:));
sc = sum(sum(C.*(C-1)/2));
ai = sum(C,2); bj = sum(C,1)';
sa = sum(ai.*(ai-1)/2); sb = sum(bj.*(bj-1)/2);
expected = sa*sb/(n*(n-1)/2);
maxIdx = (sa+sb)/2;
ari = (sc-expected)/(maxIdx-expected+eps);
end
