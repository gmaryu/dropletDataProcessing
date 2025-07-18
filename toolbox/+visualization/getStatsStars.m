function retv = getStatsStars(y1, y2)
    p = ranksum(y1, y2);

    if p > 0.05
        retv = "n.s.";
    elseif p > 0.01
        retv = "*";
    elseif p > 0.001
        retv = "**";
    else
        retv = "***";
    end
end