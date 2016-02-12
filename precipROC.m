function [roc_score, TPRprecip, FPRprecip, cutoffs, CPrecip] = precipROC(lag_range, lagsToPlot, varsToPlot, names, precip_comp, summer_doy, predictand, ...
	predictandName, GHCND_clustered, fig_folder);
  
cutoffPerc = 0:10:100;
% also use precip as predictor
summerPrecip = find(ismember(precip_comp.doy, summer_doy));
clear TPRprecip FPRprecip roc_score
for jj = 1:length(varsToPlot)
	fld = varsToPlot{jj};
	CPrecip = NaN(length(lag_range), length(summerPrecip));
	for kk = 1:length(lag_range)
		try
			CPrecip(kk,:) = fld(summerPrecip - lag_range(kk));
		end
	end
	CPrecip = -bsxfun(@times, CPrecip, 1./nanstd(CPrecip,[],2));
	makeplots = 0;

	[roc_score(:,jj), TPRprecip(:,:,jj), FPRprecip(:,:,jj), cutoffs] = ...
		calcROC(GHCND_clustered, CPrecip, [], lag_range, predictand, summer_doy, makeplots, ...
		fig_folder, [names{jj}], cutoffPerc);
end

OR = TPRprecip.*(1-FPRprecip)./(FPRprecip.*(1-TPRprecip));

if ~isempty(lagsToPlot)

	cmap = flipud(lbmap(length(lagsToPlot), 'brownblue'));
	for ct = 1:length(varsToPlot)
		figname = ['ROCplots_' names{ct} '_' predictandName '.ps'];
		clf
		hold on
		plot([0 1],[0 1],'--','color',0.8*[1 1 1],'linewidth',2)
		clear h legendText
		for kk = 1:length(lagsToPlot)
			loc = find(lag_range == lagsToPlot(kk));
			h(kk) = plot(FPRprecip(loc,:,ct), TPRprecip(loc,:,ct),'-s','linewidth',2,'color',cmap(kk,:));
			legendText{kk} = [num2str(lagsToPlot(kk)) ' days  (' num2str(round(roc_score(loc,ct)*100)/100) ')'];
		end
		LEG = legend(h,strtrim(cellstr(legendText)),'location','southeast');
		set(LEG,'fontsize',24)
		set(gca,'fontsize',24)
		set(gca,'layer','top')
		set(gca,'box','on')
		title(names{ct})
	
		axis square
		orient landscape
		set(gcf,'Renderer','zbuffer')
		%set(gcf, 'PaperPositionMode', 'manual');
		%set(gcf, 'PaperUnits', 'inches');
		%set(gcf, 'PaperPosition', [0.25 0.25 10 10]);
		print('-dpsc',[fig_folder '/' figname])
	
	
		% also plot odds ratio
	
		figname = ['ORplots_' names{ct} '_' predictandName '.ps'];
		clf
		hold on
		plot([50 90],[1 1],'--','color',0.8*[1 1 1],'linewidth',2)
		clear h legendText
		for kk = 1:length(lagsToPlot)
			loc = find(lag_range == lagsToPlot(kk));
			h(kk) = plot(cutoffPerc, OR(loc,:, ct),'-s','linewidth',2,'color',cmap(kk,:));
			if lagsToPlot(kk) == 1
				legendText{kk} = [num2str(lagsToPlot(kk)) ' month'];
			elseif lagsToPlot(kk) == 2
				legendText{kk} = [num2str(lagsToPlot(kk)) ' months'];
			else
				legendText{kk} = [num2str(lagsToPlot(kk)) ' days'];
			end
		end
		LEG = legend(h,strtrim(cellstr(legendText)),'location','eastoutside');
		set(LEG,'fontsize',24)
		xlim([50 90])
		ylim([0 6])
	
		set(gca,'fontsize',24)
		set(gca,'layer','top')
		set(gca,'box','on')
		title(names{ct})
	
		orient landscape
		set(gcf, 'PaperPositionMode', 'manual');
		set(gcf, 'PaperUnits', 'inches');
		set(gcf, 'PaperPosition', [0.25 0.25 10 5]);
		print('-dpsc',[fig_folder '/' figname])
	
	end
end