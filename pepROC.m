function [roc_score,TPR, FPR, cutoffs] = pepROC(lag_range, lagsToPlot, C, names, ...
    summer_doy, predictand, predictandName, ...
    GHCND_clustered, fig_folder)
 
 
makeplots = 0;
cutoffPerc = 0:10:100;
if strcmp(names,'NINO') % calculate for with and without 1988
     [roc_score, TPR, FPR, cutoffs] = ... % with 1988
        calcROC(GHCND_clustered, squeeze(C(1, :, :)), [], lag_range, predictand, ...
        summer_doy, makeplots, fig_folder, names, cutoffPerc);
        
    [roc_score2 TPR2, FPR2, cutoffs2] = ... % without 1988
        calcROC(GHCND_clustered, squeeze(C(2, :, :)), [], lag_range, predictand, ...
        summer_doy, makeplots, fig_folder, names, cutoffPerc);
        
    if ~isempty(lagsToPlot)
        cmap = flipud(lbmap(length(lagsToPlot), 'brownblue'));
        
        figname = ['ROCplots_' names '_' predictandName '.ps'];
        clf
        hold on
        plot([0 1],[0 1],'--','color',0.8*[1 1 1],'linewidth',2)
        clear h legendText
        for kk = 1:length(lagsToPlot)
            loc = find(lag_range == lagsToPlot(kk));
            h(kk) = plot(FPR(loc,:), TPR(loc,:),'-s','linewidth',2,'color',cmap(kk,:));
            if lagsToPlot(kk) == 1
                legendText{kk} = [num2str(lagsToPlot(kk)) ' month (' num2str(round(roc_score(loc)*100)/100) ')'];
            elseif lagsToPlot(kk) == 2
                legendText{kk} = [num2str(lagsToPlot(kk)) ' months (' num2str(round(roc_score(loc)*100)/100) ')'];
            else
                legendText{kk} = [num2str(lagsToPlot(kk)) ' days (' num2str(round(roc_score(loc)*100)/100) ')'];
            end
        end
        
        for kk = 1:length(lagsToPlot)
            loc = find(lag_range == lagsToPlot(kk));
            h(kk + length(lagsToPlot)) = plot(FPR2(loc,:), TPR2(loc,:),'--s','linewidth',2,'color',cmap(kk,:));
            if lagsToPlot(kk) == 1
                legendText{kk + length(lagsToPlot)} = [num2str(lagsToPlot(kk)) ' month (' num2str(round(roc_score(loc)*100)/100) ')'];
            elseif lagsToPlot(kk) == 2
                legendText{kk + length(lagsToPlot)} = [num2str(lagsToPlot(kk)) ' months (' num2str(round(roc_score(loc)*100)/100) ')'];
            else
                legendText{kk + length(lagsToPlot)} = [num2str(lagsToPlot(kk)) ' days (' num2str(round(roc_score(loc)*100)/100) ')'];
            end
        end
        
        
        LEG = legend(h,strtrim(cellstr(legendText)),'location','northwest');
        set(LEG,'fontsize',20)
        set(gca,'fontsize',24)
        set(gca,'layer','top')
        set(gca,'box','on')
        title(names)
        
        axis square
        orient landscape
        %set(gcf, 'PaperPositionMode', 'manual');
        %set(gcf, 'PaperUnits', 'inches');
        %set(gcf, 'PaperPosition', [0.25 0.25 10 10]);
        print('-dpsc',[fig_folder '/' figname])
    end
    
else
    
    [roc_score, TPR, FPR, cutoffs] = ...
        calcROC(GHCND_clustered, C, [], lag_range, predictand, ...
        summer_doy, makeplots, fig_folder, names, cutoffPerc);
    
     
    if ~isempty(lagsToPlot)
        cmap = flipud(lbmap(length(lagsToPlot), 'brownblue'));
        
        figname = ['ROCplots_' names '_' predictandName '.ps'];
        clf
        hold on
        plot([0 1],[0 1],'--','color',0.8*[1 1 1],'linewidth',2)
        clear h legendText
        for kk = 1:length(lagsToPlot)
            loc = find(lag_range == lagsToPlot(kk));
            h(kk) = plot(FPR(loc,:), TPR(loc,:),'-s','linewidth',2,'color',cmap(kk,:));
            if lagsToPlot(kk) == 1
                legendText{kk} = [num2str(lagsToPlot(kk)) ' month (' num2str(round(roc_score(loc)*100)/100) ')'];
            elseif lagsToPlot(kk) == 2
                legendText{kk} = [num2str(lagsToPlot(kk)) ' months (' num2str(round(roc_score(loc)*100)/100) ')'];
            else
                legendText{kk} = [num2str(lagsToPlot(kk)) ' days (' num2str(round(roc_score(loc)*100)/100) ')'];
            end
        end
        if strcmp(names, 'ATL')
        	LEG = legend(h,strtrim(cellstr(legendText)),'location','northwest');
        	set(LEG,'fontsize',20)
        else
        	LEG = legend(h,strtrim(cellstr(legendText)),'location','southeast');
        	set(LEG,'fontsize',24)
        end
        
        set(gca,'fontsize',24)
        set(gca,'layer','top')
        set(gca,'box','on')
        title(names)
        
        axis square
        orient landscape
        %set(gcf, 'PaperPositionMode', 'manual');
        %set(gcf, 'PaperUnits', 'inches');
        %set(gcf, 'PaperPosition', [0.25 0.25 10 10]);
        print('-dpsc',[fig_folder '/' figname])
        
        % also plot odds ratio
        OR = TPR.*(1-FPR)./(FPR.*(1-TPR));
        figname = ['ORplots_' names '_' predictandName '.ps'];
        clf
        hold on
        plot([50 90],[1 1],'--','color',0.8*[1 1 1],'linewidth',2)
       	clear h legendText
        for kk = 1:length(lagsToPlot)
        	loc = find(lag_range == lagsToPlot(kk));
            h(kk) = plot(cutoffPerc, OR(loc,:),'-s','linewidth',2,'color',cmap(kk,:));
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
        title(names)
        
        orient landscape
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperPosition', [0.25 0.25 10 5]);
        print('-dpsc',[fig_folder '/' figname])
            
        
    end
    
    
    
end
