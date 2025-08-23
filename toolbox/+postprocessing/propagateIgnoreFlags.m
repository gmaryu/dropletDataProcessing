function dataOut = propagateIgnoreFlags(dataSet, mask, varargin)
% propagateIgnoreFlags - Propagate IGNORED flags from info to cycle/timeSeries via (POS_ID, TRACK_ID).
%
% Usage:
%   dataOut = propagateIgnoreFlags(dataIn, mask)
%       % mask: logical vector of length height(dataIn.info)
%
%   dataOut = propagateIgnoreFlags(dataIn, mask, 'KeyFields', {'POS_ID','TRACK_ID'})
%       % customize key field names if necessary
%
% Inputs:
%   dataIn.info, dataIn.cycle, dataIn.timeSeries : tables (must contain key fields)
%   mask : logical vector, same length as dataIn.info, specifying rows to ignore
%
% Outputs:
%   dataOut : copy of dataIn, with IGNORED column updated in all tables

% --- parse options ---
p = inputParser;
addParameter(p, 'KeyFields', {'POS_ID','TRACK_ID'}, @(c)iscellstr(c) || all(cellfun(@isstring,c)));
parse(p, varargin{:});
opt = p.Results;
keys  = cellstr(string(opt.KeyFields));

% --- ensure IGNORED column exists and is logical ---
infoT = ensureIgnoredVar(dataSet.info);
cycT  = ensureIgnoredVar(dataSet.cycle);
tsT   = ensureIgnoredVar(dataSet.timeSeries);

% --- apply mask to info ---
assert(islogical(mask) && numel(mask)==height(infoT), ...
    'Mask must be logical with length equal to height(info).');
infoT.IGNORED = mask;

% --- prepare key matrices ---
keyInfo = infoT{:, keys};
keyCyc  = cycT{:,  keys};
keyTS   = tsT{:,   keys};

% --- set of keys to ignore ---
keyIgnore = keyInfo(mask, :);

% --- vectorized matching ---
ixC = ismember(keyCyc, keyIgnore, 'rows');
ixT = ismember(keyTS, keyIgnore, 'rows');

cycT.IGNORED(ixC) = true;
tsT.IGNORED(ixT)  = true;

% --- pack output ---
dataOut = dataSet;
dataOut.info       = infoT;
dataOut.cycle      = cycT;
dataOut.timeSeries = tsT;
end

% --- helper: add IGNORED column if missing ---
function T = ensureIgnoredVar(T)
    if ~ismember("IGNORED", T.Properties.VariableNames)
        T.IGNORED = false(height(T),1);
    elseif ~islogical(T.IGNORED)
        T.IGNORED = logical(T.IGNORED);
    end
end
