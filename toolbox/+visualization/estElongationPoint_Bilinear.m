function dataSet = estElongationPoint_Bilinear(dataSet, varargin)

% ---------- parse options ----------
p = inputParser;
addParameter(p,'plotFig',false,@islogical);
addParameter(p,'savePath',"",@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(p,'Xvar',"NCVR",@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(p,'Yvar',"DURATION",@(x)ischar(x)||(isstring(x)&&isscalar(x)));
parse(p,varargin{:});
opt = p.Results;
Xvar = opt.Xvar;
Yvar = opt.Yvar;

% --- figure setting ---
if opt.plotFig & isempty(opt.savePath)
    disp("savePath is required for saving plots")
    return
end
savePath = fullfile(opt.savePath, 'KneeFitting');

%%
info = dataSet.info;
cycle = dataSet.cycle;


% --- data allocation ---
knee_vals = nan*zeros(size(info,1),1);
columnNameK = sprintf('KneeVlas_%s',Xvar);
%%
for i = 1:size(info,1)
    if ~info.IGNORED(i)
        %disp(i)
        tmp_POS = info.POS_ID(i);
        tmp_DID = info.TRACK_ID(i);
        tmp_cycles = cycle(cycle.POS_ID==tmp_POS & cycle.TRACK_ID == tmp_DID,:);

        % sample data
        ydata = tmp_cycles.(Yvar)(~isnan(tmp_cycles.(Xvar)));
        xdata = tmp_cycles.(Xvar)(~isnan(tmp_cycles.(Xvar)));
        %{
        if strcmp(Xvar,"NCVR") | strcmp(Xvar,"NCVR_ORI") | strcmp(Xvar,"DNACR")
            xdata = log10(tmp_cycles.(Xvar)(~isnan(tmp_cycles.(Xvar))));
        else
            xdata = tmp_cycles.(Xvar)(~isnan(tmp_cycles.(Xvar)));
        end
        %}

        % fitting
        if numel(xdata) > 4
            xi = xdata;
            yi = ydata;
            if ~all(xi > 0)
                fprintf('Skipping index %d: negative or zero values detected.\n', i);
                continue; 
            end

            [xi,ord] = sort(xi(:));
            yi = yi(:);
            yi = yi(ord);
            if strcmp(Xvar,"NCVR") | strcmp(Xvar,"NCVR_ORI") | strcmp(Xvar,"DNACR")
                lxi = log10(xi);
            else
                break
            end
            
            % ---- Knee detection (bilinear fit) ----
            out = knee_bilinear_continuous(lxi, yi);
            knee_vals(i) = 10.^out.lx_knee;

            if opt.plotFig

                f = figure('Visible','off');
                f.Units = 'centimeters';
                f.Position = [0,0,5,5];

                semilogx(xi, yi, '.', 'MarkerSize', 8); hold on;
                semilogx(xi, out.yhat, '-', 'LineWidth', 1.5);
                xlabel(Xvar,'Interpreter','none');
                ylabel(Yvar,'Interpreter','none');
                grid on;


                if ~exist(savePath,'dir')
                    mkdir(savePath)
                end
                exportgraphics(f, fullfile(savePath,sprintf('Pos%d_TRACK_ID_%03d.png',tmp_POS, tmp_DID)));
            end
        else
            fprintf('Pos%d TRACK_ID_%03d does not have enough data point\n',tmp_POS, tmp_DID);
            knee_vals(i) = nan;
            
        end
    end

end

dataSet.info.(columnNameK) = knee_vals;


%% --------- subfunction ---------
function out = knee_bilinear_continuous(lx, y)
% Knee by continuous 2-segment linear fit on lx=log10(x)
    n = numel(lx);
    kmin = 2; kmax = n-2;
    best.SSE = inf;
    for k = kmin:kmax
        % Left side
        p1 = polyfit(lx(1:k), y(1:k), 1);
        yk = polyval(p1, lx(k));
        % Right side
        dxR = lx(k+1:n) - lx(k);
        dyR = y(k+1:n) - yk;
        denom = sum(dxR.^2);
        if denom <= eps, a2 = 0; else, a2 = sum(dxR.*dyR)/denom; end
        yhat = [polyval(p1,lx(1:k)); a2*(lx(k+1:n)-lx(k))+yk];
        SSE = sum((y - yhat).^2);
        if SSE < best.SSE
            best.SSE = SSE; best.k=k; best.p1=p1; best.a2=a2; best.yhat=yhat;
        end
    end
    if isempty(k)
        out.idx_knee = NaN;
        out.lx_knee  = NaN;
        out.yhat     = NaN;
        disp('error')
    else
        out.idx_knee = best.k;
        out.lx_knee  = lx(best.k);
        out.yhat     = best.yhat;
    end
end

end