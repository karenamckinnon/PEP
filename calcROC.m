function [roc_score, TPR, FPR, cutoffs] = ...
    calcROC(GHCND_clustered, C, cone_ratio, lag_range, predictand, ...
    summer_doy, makeplots, fig_folder, domain, cutoffPerc)
     
% cutoffPerc = 0:10:100;
Ncutoff = length(cutoffPerc) - 1;
startVec = double(ismember(GHCND_clustered.time, predictand));
origSummer = startVec(ismember(GHCND_clustered.doy, summer_doy));
fullSummerLength = length(origSummer);
summerDates = GHCND_clustered.time(ismember(GHCND_clustered.doy, summer_doy));

TPR = NaN(length(lag_range), Ncutoff);
FPR = TPR;
roc_score = NaN(size(lag_range));

for jj = 1:length(lag_range)
	cutoffs = prctile(C(jj,:), cutoffPerc);
	lag = lag_range(jj);
	% window_size = abs(floor(cone_ratio * lag));
	if lag_range > 30 | year(GHCND_clustered.time(1)) == 1950 % latter for HadSST3 analysis
		window_size = 3;
	else
		window_size = 0;
	end
	ts = runmean(startVec, window_size);
	tsSummer = ts(ismember(GHCND_clustered.doy, summer_doy));

	for ct = 1:length(cutoffs)

		predUse = C(jj, :) > cutoffs(ct);
		TP = sum(predUse == 1 & tsSummer > 0);
		FP = sum(predUse == 1 & tsSummer == 0);
		totalPos = sum(tsSummer > 0);
		totalNeg = sum(tsSummer == 0);
		TPR(jj,ct) = TP/totalPos;
		FPR(jj,ct) = FP/totalNeg;

		
	end


	[x,ind] = sort(FPR(jj, :),'ascend');
	y = TPR(jj, :);
	y = y(ind);
	dx = diff(x);
	roc_score(jj) = sum(dx.*y(1:end-1) + 1/2*dx.*(y(2:end) - y(1:end-1)));

end



% Plot ROC scores for different lead times
if makeplots
    cutoffUse = 0.5;
    plotLags = [20 30 40 50];
    for ct = 1:length(plotLags)
        idx = find(lag_range == plotLags(ct));
        clf
        plot([0 1],[0 1],'--','color',0.8*[1 1 1],'linewidth',2)
        hold on
        plot(FPR(idx, :), TPR(idx, :) ,'-sk','linewidth',2)
        grid on
        set(gca,'layer','top')
        set(gca,'box','on')
        title(['ROC score: ' num2str(round(100*roc_score(idx))/100) ''])
        orient landscape
        print('-dpsc',[fig_folder '/ROCplots_' domain '_window_' num2str(cone_ratio) '_lag_' num2str(plotLags(ct)) '.ps'])
    end

    clf
    plot([0 80],[0.5 0.5],'--','color',0.8*[1 1 1],'linewidth',2)
    hold on

    plot(lag_range, roc_score0, '-','color',0.5*[1 1 1],'linewidth',2)
    plot(lag_range, roc_score, '-k', 'linewidth',2)
    xlim([0 80])
    % ylim([0.4 0.7])
    xlabel('Lead time')
    ylabel('ROC score')
    set(gca,'xdir','reverse')
    orient landscape
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0.25 0.25 10 5]);
    print('-dpsc',[fig_folder '/ROCscore_' domain '_window_' num2str(cone_ratio) '.ps'])

end

