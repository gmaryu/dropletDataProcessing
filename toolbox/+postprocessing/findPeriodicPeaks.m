function retv = findPeriodicPeaks(signal, frameToMin)
% findPeriodicPeaks Detects periodic peaks in a signal.
%
%   retv = findPeriodicPeaks(signal, frameToMin)
%
% The function uses findpeaks on both the signal and its negative to determine peak and trough
% positions. It returns a matrix with each row as [start_index, end_index, trough_index] if the peaks
% and troughs match the expected pattern; otherwise, it returns NaN.

    p = 0.08;
    maxw = 60 / frameToMin;  % Maximum expected peak width in frames
    
    [~, ip] = findpeaks(signal, "MinPeakProminence", p, "MaxPeakWidth", maxw);
    [~, it] = findpeaks(-signal, "MinPeakProminence", p);
    
    if numel(ip) == numel(it)
        if all(it - ip > 0)
            retv = [ip(1:end-1), ip(2:end), it(1:end-1)];
        elseif all(it - ip < 0)
            retv = [ip(1:end-1), ip(2:end), it(2:end)];
        else
            retv = nan;
        end
    elseif numel(ip) == numel(it) + 1
        if all(it - ip(1:end-1) > 0)
            retv = [ip(1:end-1), ip(2:end), it(1:end)];
        else
            retv = nan;
        end
    elseif numel(ip) +1 == numel(it)
        % disp('first peak is too close to 0');
        % discard first trough
        if all(it(2:end) - ip(1:end) > 0)
            retv = [ip(1:end-1), ip(2:end), it(2:end-1)];
        else
            retv = nan;
        end
    else
        retv = nan;
    end
end