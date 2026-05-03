function output = foldChangeDNA(tm, tp, varargin)


% ---------- parse options ----------
p = inputParser;
addParameter(p,'StartIndexMin',5,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'RespectIgnored',true,@islogical);
parse(p,varargin{:});
thresh = p.Results.StartIndexMin;
respectIgnored = p.Results.RespectIgnored;

T = tp;
%%
% ---------- safety: ensure required columns exist ----------
req = ["TRACK_ID","START_INDEX","CYCLE_ID","DNA_SUM_INT_MOD_Q90"];
missing = req(~ismember(req, string(T.Properties.VariableNames)));
assert(isempty(missing), 'Missing required column(s): %s', strjoin(cellstr(missing), ', '));

% ensure IGNORED column exists and is logical
if ~ismember("IGNORED", string(T.Properties.VariableNames))
    T.IGNORED = false(height(T),1);
elseif ~islogical(T.IGNORED)
    T.IGNORED = logical(T.IGNORED);
end


% ---------- stable re-numbering per TRACK_ID ----------
% Keep original order index, then sort by TRACK_ID, START_INDEX to define sequence

uniquePairs = unique(tm(:, {'POS_ID', 'TRACK_ID'}), 'rows');

for i = 1:size(uniquePairs,1)
    % time series data
    dna_int = tm.SUM_SPERM_HOECHST_INT(tm.POS_ID == uniquePairs.POS_ID(i) & tm.TRACK_ID == uniquePairs.TRACK_ID(i));

    % peak data
    tmpDNAQ90 = T.DNA_SUM_INT_MOD_Q90(T.POS_ID == uniquePairs.POS_ID(i) & T.TRACK_ID == uniquePairs.TRACK_ID(i)); 

    firstPeakIdx = min(T.START_INDEX(T.POS_ID == uniquePairs.POS_ID(i) & T.TRACK_ID == uniquePairs.TRACK_ID(i)));
    if length(tmpDNAQ90) > 1 & firstPeakIdx < 20 % short interphase before the first peak
        fc = [1; tmpDNAQ90(2:end)./tmpDNAQ90(1:end-1)];
        T.FC_DNA(T.POS_ID == uniquePairs.POS_ID(i) & T.TRACK_ID == uniquePairs.TRACK_ID(i)) = fc;
    elseif length(tmpDNAQ90) > 1 & firstPeakIdx >= 20 % enough frames before the first peak
        dna_before_1st = median(dna_int(1:firstPeakIdx));
        fc = [tmpDNAQ90(1)/dna_before_1st; tmpDNAQ90(2:end)./tmpDNAQ90(1:end-1)];
        T.FC_DNA(T.POS_ID == uniquePairs.POS_ID(i) & T.TRACK_ID == uniquePairs.TRACK_ID(i)) = fc;
    else
        T.FC_DNA(T.POS_ID == uniquePairs.POS_ID(i) & T.TRACK_ID == uniquePairs.TRACK_ID(i)) = NaN;
    end
end
T.CYCLE_ID(T.IGNORED) = NaN;
output = T.FC_DNA;
end
