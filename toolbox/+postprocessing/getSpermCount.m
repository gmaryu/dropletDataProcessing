function spermCount = getSpermCount(spermRef, dropletID)
% getSpermCount Retrieve the sperm count for a given droplet.
%
%   spermCount = getSpermCount(spermRef, dropletID)
%
% Inputs:
%   spermRef  - A table containing at least the columns 'DropID' and 'Count'.
%   dropletID - A scalar numeric identifier for the droplet.
%
% Output:
%   spermCount - The sperm count value (if found) or NaN if not available.
%
% Example:
%   sc = getSpermCount(spermRef, 5);

    arguments
        spermRef table
        dropletID (1,1) double
    end

    row = spermRef(spermRef.DropID == dropletID, :);
    if isempty(row)
        spermCount = [];
    else
        spermCount = row.Count;
    end
end
