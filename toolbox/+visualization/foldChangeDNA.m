function output = foldChangeDNA(dataCycle, varargin)
% remapCycleIdByStartIndex - Ignore early peaks (START_INDEX<=thresh) and
%                            re-label CYCLE_ID per TRACK_ID in START_INDEX order.
%
% Effect:
%   - Rows with START_INDEX <= thresh (default 5) are flagged IGNORED (union with existing flags).
%   - Within each TRACK_ID, non-ignored rows are assigned CYCLE_ID = 1,2,3,... by
%     ascending START_INDEX.
%   - IGNORED rows receive CYCLE_ID = NaN (so downstream code can skip them).
%
% Inputs:
%   dataSet.cycle : table with columns TRACK_ID, START_INDEX, CYCLE_ID (numeric), IGNORED (logical/0-1)
%
% Name-Value:
%   'StartIndexMin' (5)      : threshold; START_INDEX <= this is ignored
%   'RespectIgnored' (true)  : if true, existing IGNORED==true stays ignored
%
% Output:
%   dataSet.cycle : updated table with CYCLE_ID overwritten

% ---------- parse options ----------
p = inputParser;
addParameter(p,'StartIndexMin',5,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'RespectIgnored',true,@islogical);
parse(p,varargin{:});
thresh = p.Results.StartIndexMin;
respectIgnored = p.Results.RespectIgnored;

T = dataCycle;
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

% Identify first row of each track after sorting (group boundary)
isFirstOfTrack = [true; diff(T.TRACK_ID)~=0];
grp = cumsum(isFirstOfTrack);            % group id per row

ids = unique(grp);
for i = 1:numel(ids)
    % disp(i);
    tmpDNAQ90 = T.DNA_SUM_INT_MOD_Q90(grp==ids(i));
    if length(tmpDNAQ90) > 1
        fc = [1; tmpDNAQ90(2:end)./tmpDNAQ90(1:end-1)];
        T.FC_DNA(grp==ids(i)) = fc;
    else
        T.FC_DNA(grp==ids(i)) = 1;
    end
end
T.CYCLE_ID(T.IGNORED) = NaN;
output = T.FC_DNA;
end
