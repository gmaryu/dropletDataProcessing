function out = runPhaseClassification(nucData, dataSet, traceInfo, resultDataPath, batchName, conditionLabel, saveDir, applyModel)
%RUNPHASECLASSIFICATION  End-to-end Stage A phase classification for one batch.
%
% Shape-only (N/C-, period-, condition-free) 3-state classification of Cdk1
% cycles, with an unsupervised GMM (primary) and a data-driven rule (independent
% confirmation). Post-hoc quantities (N/C, period, cycle, sperm, diameter) are
% attached for Stage B but NEVER used for classification. Results are saved so
% multiple batches / conditions can be pooled later.
%
% Label convention: cycleTable carries both PhaseRule and PhaseGMM.
% PhaseRule is the primary/reported label; PhaseGMM is the supporting
% unsupervised label used to derive the rule thresholds and to compute the
% agreement (kappa) statistic. Downstream phase-dependent analyses should use
% PhaseRule (see README_phase_classification.md sec 3).
%
% INPUTS
%   nucData, dataSet, traceInfo, resultDataPath : as for extractShapeFeatures
%   batchName        : e.g. "20240821_Control"
%   conditionLabel   : e.g. "control" | "CHIR124_1uM"
%   saveDir          : folder to write <batchName>_<conditionLabel>_phases.mat
%   applyModel (opt) : a previously returned out.gmm to APPLY (fixed-model
%                      transfer) instead of fitting a fresh GMM. Features are
%                      z-scored WITHIN this batch first, then the fixed cluster
%                      geometry is applied (structure transfer, not absolute
%                      thresholds).
%
% OUTPUT struct with fields: batchName, conditionLabel, featNames, cycleTable
%   (one row per valid cycle: DropletID, CycleInDroplet, 3 shape features,
%    PhaseGMM, PhaseRule, NCratio, Period, SpermCount, Diameter), gmm, rule,
%    agreement, bic, containment.

if nargin < 8, applyModel = []; end
featNames = ["CytoDistanceCorr","InterphaseNucCytoRatioP90","TimeToHalfActivation_nuc"];

% ---- features ----
T = postprocessing.extractShapeFeatures(nucData, dataSet, traceInfo, resultDataPath);

valid = T.IsValidCycle;
X = nan(height(T), numel(featNames));
for j = 1:numel(featNames), X(:,j) = T.(char(featNames(j))); end
rowok = valid & all(~isnan(X),2);
Xv = X(rowok,:);
viok = find(rowok);

% ---- within-batch standardization (store transform) ----
if ~isempty(applyModel)
    % structure transfer: z-score within THIS batch, then apply fixed geometry
    mu_ = mean(Xv,1); sd_ = std(Xv,0,1); sd_(sd_==0)=1;
    Xz = (Xv - mu_)./sd_;
    gm = applyModel.model;
    rawlbl = cluster(gm, Xz);
    phaseGMM = applyModel.phaseMap(rawlbl);
    bic = []; gmmOut = applyModel; gmmOut.transferredTo = batchName;
else
    mu_ = mean(Xv,1); sd_ = std(Xv,0,1); sd_(sd_==0)=1;
    Xz = (Xv - mu_)./sd_;
    % BIC scan 1..6 (for the record)
    bic = inf(6,1);
    for k=1:6
        g = fitgmdist(Xz,k,'CovarianceType','full','RegularizationValue',1e-4, ...
            'Replicates',30,'Options',statset('MaxIter',1000));
        bic(k)=g.BIC;
    end
    % 3-state model (justified separately by containment analysis)
    gm = fitgmdist(Xz,3,'CovarianceType','full','RegularizationValue',1e-4, ...
        'Replicates',50,'Options',statset('MaxIter',1000));
    rawlbl = cluster(gm, Xz);
    % Order clusters into phases 1..3 by cytoplasmic-template DISTANCE alone.
    % Rationale: CytoDistanceCorr is (a) computed against a within-batch cyto-only
    % reference (batch-normalized), and (b) the only shape axis that is monotone
    % across phases in every batch tested. Ordering by it is circularity-free
    % (uses no N/C, period, or condition) and defines Phase 1 = most cytoplasm-like,
    % consistent with the biological definition. Other axes (interphase nuclear
    % activation p90, activation timing tth) are non-monotone across batches and
    % are therefore NOT used to number phases.
    cm = arrayfun(@(c) median(Xv(rawlbl==c,1)), 1:3);  % raw cyto distance median
    [~,ord] = sort(cm,'ascend');
    phaseMap = zeros(3,1); for n=1:3, phaseMap(ord(n))=n; end
    phaseGMM = phaseMap(rawlbl);
    gmmOut = struct('model',gm,'phaseMap',phaseMap,'mean_',mu_,'std_',sd_, ...
        'featNames',featNames,'cytoOrderMedians',cm);
end

% ---- data-driven rule (reproduces GMM geometry; no hardcoded literals) ----
rawFeat = Xv;  % unstandardized [cyto, p90, tth]
med = @(lab,col) median(rawFeat(phaseGMM==lab,col));
thr.cyto12 = mean([med(1,1), med(2,1)]);
thr.cyto23 = mean([med(2,1), med(3,1)]);
thr.tth23  = mean([med(2,3), med(3,3)]);
thr.p9023  = mean([med(2,2), med(3,2)]);
N = size(rawFeat,1); phaseRule = zeros(N,1);
for i=1:N
    isP3 = rawFeat(i,1) >= thr.cyto23 && rawFeat(i,3) >= thr.tth23 && rawFeat(i,2) <= thr.p9023;
    if isP3, phaseRule(i)=3;
    elseif rawFeat(i,1) < thr.cyto12, phaseRule(i)=1;
    else, phaseRule(i)=2; end
end

% ---- agreement ----
C=zeros(3); for a=1:3,for b=1:3,C(a,b)=sum(phaseGMM==a & phaseRule==b);end;end
po=sum(diag(C))/N; pe=sum(sum(C,2).*sum(C,1)')/N^2;
agreement=struct('confusion',C,'pctAgree',po,'kappa',(po-pe)/(1-pe));

% ---- post-hoc quantities (Stage B only) ----
NCratio = T.InterphaseNCvolRatioMean(viok);
Period  = T.CycleLength(viok);
Cycle   = T.CycleInDroplet(viok);
DropletID = string(T.DropletID(viok));
% sperm & diameter from dataSet.info via POS/TRACK
Sperm = nan(N,1); Diam = nan(N,1);
posC = T.Position(viok); trkC = T.Droplet(viok);
for i=1:N
    r = dataSet.info.POS_ID==posC(i) & dataSet.info.TRACK_ID==trkC(i);
    if any(r)
        Sperm(i)=dataSet.info.SPERM_COUNT(find(r,1));
        Diam(i)=dataSet.info.ORIGINAL_MED_DIAMETER(find(r,1));
    end
end

cycleTable = table(DropletID, Cycle, rawFeat(:,1), rawFeat(:,2), rawFeat(:,3), ...
    phaseGMM, phaseRule, NCratio, Period, Sperm, Diam, ...
    'VariableNames', {'DropletID','CycleInDroplet','CytoDistanceCorr', ...
    'InterphaseNucCytoRatioP90','TimeToHalfActivation_nuc','PhaseGMM','PhaseRule', ...
    'NCratio','Period','SpermCount','Diameter'});
cycleTable.BatchName = repmat(string(batchName), N, 1);
cycleTable.Condition = repmat(string(conditionLabel), N, 1);

% ---- containment k=4,5 -> k=3 (only when fitting) ----
containment = struct();
if isempty(applyModel)
    for kk=[4 5]
        gk=fitgmdist(Xz,kk,'CovarianceType','full','RegularizationValue',1e-4, ...
            'Replicates',30,'Options',statset('MaxIter',1000));
        lk=cluster(gk,Xz);
        % which k3 phase each k-cluster maps to
        M=zeros(kk,3); for a=1:kk,for b=1:3,M(a,b)=sum(lk==a & phaseGMM==b);end;end
        containment.(sprintf('k%d',kk))=M;
    end
end

out = struct('batchName',string(batchName),'conditionLabel',string(conditionLabel), ...
    'featNames',{featNames},'cycleTable',cycleTable,'gmm',gmmOut,'rule',thr, ...
    'agreement',agreement,'bic',bic,'containment',containment);

% ---- save ----
if ~isempty(saveDir)
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    fn = fullfile(saveDir, sprintf('%s_%s_phases.mat', batchName, conditionLabel));
    save(fn, 'out');
    writetable(cycleTable, fullfile(saveDir, sprintf('%s_%s_cycles.csv', batchName, conditionLabel)));
    fprintf('saved: %s\n', fn);
end

fprintf('%s / %s : n=%d cycles | GMM phases %d/%d/%d | rule kappa=%.3f\n', ...
    batchName, conditionLabel, N, sum(phaseGMM==1),sum(phaseGMM==2),sum(phaseGMM==3), agreement.kappa);
end
