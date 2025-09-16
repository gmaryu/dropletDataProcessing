function dataSet = estElongationPoint(dataSet, varargin)

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
savePath = fullfile(opt.savePath, 'fitting');

%%
info = dataSet.info;
cycle = dataSet.cycle;

%% fitting-equation
modelfun = @(b,x) b(1)+(b(2)./(1+exp(-b(3).*(x-b(4)))));

% --- data allocation ---
fitting_K = nan*zeros(size(info,1),1);
fitting_rmid = nan*zeros(size(info,1),1);
columnNameK = sprintf('fitting_K_%s',Xvar);
columnNameRmid = sprintf('fitting_rmid_%s',Xvar);
%%
for i = 1:size(info,1)
    if ~info.IGNORED(i)
        disp(i)
        tmp_POS = info.POS_ID(i);
        tmp_DID = info.TRACK_ID(i);
        tmp_cycles = cycle(cycle.POS_ID==tmp_POS & cycle.TRACK_ID == tmp_DID,:);

        % sample data
        ydata = tmp_cycles.(Yvar)(~isnan(tmp_cycles.(Xvar)));
        if strcmp(Xvar,"NCVR") | strcmp(Xvar,"NCVR_ORI") | strcmp(Xvar,"DNACR")
            xdata = log10(tmp_cycles.(Xvar)(~isnan(tmp_cycles.(Xvar))));
        else
            xdata = tmp_cycles.(Xvar)(~isnan(tmp_cycles.(Xvar)));
        end

        % fitting
        if numel(xdata) > 4
            p0 = [min(ydata), 2*(max(ydata)-min(ydata)), 1.0, 0.5];
            logistic_transition = fittype("a + b ./ (1 + exp(-c * (x - d)))", "independent", "x", "dependent", ...
                "y", "coefficients", {'a', 'b', 'c', 'd'});
            fitresult = fit(xdata, ydata, logistic_transition,"Robust","Bisquare", "StartPoint", p0);
            fitting_K(i) = fitresult.c;
            fitting_rmid(i) = fitresult.d;


            % %        beta0 = [min(ydata), 2*(max(ydata)-min(ydata)), 1.0, (min(xdata)+max(xdata))/2];
            %    beta0 = [min(ydata), 2*(max(ydata)-min(ydata)), 1.0, 0.5];
            %    opts=optimoptions("lsqcurvefit", 'MaxFunctionEvaluations',1e4, 'MaxIterations',1e4);
            %    lb = [0,0,0,0]; ub = [200,200,10,1];
            %    beta = lsqcurvefit(modelfun,beta0,xdata,ydata,lb,ub,opts)
            %    %beta = lsqcurvefit(modelfun,beta0,xdata,ydata);
            %    fitting_K(i) = beta(3);
            %    fitting_rmid(i) = beta(4);

            if opt.plotFig

                f = figure('Visible','off');
                f.Units = 'centimeters';
                f.Position = [0,0,5,5];

                if strcmp(Xvar,"NCVR") | strcmp(Xvar,"NCVR_ORI")| strcmp(Xvar,"DNACR")
                    fitX = log10(min(tmp_cycles.(Xvar)(~isnan(tmp_cycles.(Xvar))))...
                        :0.001:max(tmp_cycles.(Xvar)(~isnan(tmp_cycles.(Xvar)))));
                else
                    fitX = linspace(min(tmp_cycles.(Xvar)(~isnan(tmp_cycles.(Xvar)))),...
                        max(tmp_cycles.(Xvar)(~isnan(tmp_cycles.(Xvar)))),100);
                end

                yhat = fitresult(fitX);

                hold on
                plot(tmp_cycles.(Xvar), tmp_cycles.(Yvar),'o');
                if strcmp(Xvar,"NCVR") | strcmp(Xvar,"NCVR_ORI")| strcmp(Xvar,"DNACR")
                    plot(power(10,fitX), yhat, '-');
                else
                    plot(fitX, yhat, '-');
                end

                xlabel(Xvar,'Interpreter','none');
                ylabel(Yvar,'Interpreter','none');
                title(sprintf('K:%.2g, Rmid:%.2g',fitresult.c,power(10,fitresult.d)));
                if strcmp(Xvar,"NCVR") | strcmp(Xvar,"NCVR_ORI")| strcmp(Xvar,"DNACR")
                    xscale('log')
                end
                %xlim([0,0.15])
                hold off

                if ~exist(savePath,'dir')
                    mkdir(savePath)
                end
                exportgraphics(f, fullfile(savePath,sprintf('Pos%d_TRACK_ID_%03d.png',tmp_POS, tmp_DID)));
            end
        else

            fitting_K(i) = nan;
            fitting_rmid(i) = nan;
        end
    end

end

dataSet.info.(columnNameK) = fitting_K;
dataSet.info.(columnNameRmid) = power(10,fitting_rmid);


%%
% figure();
% mask = control_nuc.cycle.POS_ID == 10 & control_nuc.cycle.TRACK_ID==31;
% adroplet_cycle_period = control_nuc.cycle.DURATION(mask);
% adroplet_cycle_DNA = control_nuc.cycle.DNA_SUM_INT_Q90(mask);
% adroplet_cycle_NCVR = control_nuc.cycle.NCVR(mask);
%
% %xlim([0,0.1])
%
%
% xdata = log10(adroplet_cycle_NCVR(~isnan(adroplet_cycle_NCVR)))
% ydata = adroplet_cycle_period(~isnan(adroplet_cycle_NCVR))
% modelfun = @(b,x) b(1)+(b(2)./(1+exp(-b(3).*(x-b(4)))));
% beta0 = [min(y), 2*(max(y)-min(y)), 1.0, 0.5];
% bata = lsqcurvefit(modelfun,beta0,xdata,ydata)
%
%
%
% figure();
% hold on
% plot(adroplet_cycle_NCVR, adroplet_cycle_period,'o');
% plot(power(10,xdata), yhat,'*')
% hold off
%
% power(10,-1.18)
end

