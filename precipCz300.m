function precipCz300(C, lag_range, GHCND_clustered, precip_comp, hgt300_comp, ...
	sst_comp, summerDays, precip_cutoff, fig_folder)
  
% Look at relationship between precip and long-lead PEP
CEarly = mean(C(lag_range >= 15 & lag_range <= 45, :), 1);
highCDays = CEarly > 1;
highDates = GHCND_clustered.time(summerDays);
highDates =  highDates(highCDays);

% Want to look at averages over the 45-15 days before these dates
highDatesMat = sort(unique(bsxfun(@minus, highDates, [15:45]')));

% Get SPI30 at a lead time of 15 days (synchronous)
spiIdx = ismember(precip_comp.time, highDates - 15);
spiComp = mean(precip_comp.SPI30(:,:,spiIdx), 3);
spiCompEUS = mean(precip_comp.spi30(spiIdx));

% Get average z300 over same time
z300Comp = NaN(size(hgt300_comp.anom,1), size(hgt300_comp.anom,2), length(highDates));
sstComp = NaN(size(sst_comp.anom, 1), size(sst_comp.anom, 2), length(highDates));
for ct = 1:length(highDates)
	z300Idx = find(hgt300_comp.time == highDates(ct));
	z300Comp(:,:,ct) = mean(hgt300_comp.anom(:,:,(z300Idx - 45):(z300Idx - 15)), 3);
	sstIdx = find(sst_comp.time == highDates(ct));
	sstComp(:,:,ct) = mean(sst_comp.anom(:,:,(sstIdx - 45):(sstIdx - 15)), 3);
end 

sstComp = mean(sstComp, 3);
z300Comp = mean(z300Comp, 3);
 

% precip_cutoff = 0.5;
lowPrecip = spiComp';
lowPrecip(lowPrecip > -precip_cutoff) = NaN;
highPrecip = spiComp';
highPrecip(highPrecip < precip_cutoff) = NaN;

dry = rgb('sienna');
wet = rgb('darkgreen');

clf
m_proj('miller','lat',[20 65],'lon',[120 330]);
hold on

m_pcolorKM(sst_comp.lon, sst_comp.lat, sstComp');
shading interp
h = colorbar;
set(h,'fontsize',16)
colormap(flipud(lbmap(21,'redblue')));
caxis([-1 1])
hold on
%m_coast('line','color','k');
%m_grid('box','fancy');
m_coast('patch',[.7 .7 .7],'edgecolor','none');
m_grid('box','fancy');
[LON LAT] = meshgrid(precip_comp.lon, precip_comp.lat);
m_plot(LON(~isnan(lowPrecip)), LAT(~isnan(lowPrecip)), 'o', 'color', dry, 'linewidth', 1, 'markersize', 3)
m_plot(LON(~isnan(highPrecip)), LAT(~isnan(highPrecip)), '+','color', wet,'linewidth', 1, 'markersize', 3)

clines = linspace(-80,80,17);
m_contour(hgt300_comp.lon,hgt300_comp.lat,z300Comp',clines(clines>0),'linewidth',2,'col','k','linestyle','-');
m_contour(hgt300_comp.lon,hgt300_comp.lat,z300Comp',clines(clines<0),'linewidth',2,'col','k','linestyle','--');

set(gca,'fontsize',16)


set(gcf,'Renderer','zbuffer')
orient landscape
print('-dpsc',[fig_folder '/SST_z300_SPI.ps'])
