%% Retrieve fCO2 global maps from SOCAT, NEMO-PlankTOM12, and GOBMs
clearvars

ywanted = 1990:2022;
YWANTED = sort(repmat(ywanted,1,12));

% read SOCAT
file = 'data/SOCATv2023_tracks_gridded_monthly.nc';

SOCAT = ncread(file,'fco2_ave_unwtd');
time  = ncread(file,'tmnth');
lat   = ncread(file,'ylat');
lon   = ncread(file,'xlon');
time  = round(datenum(1970,1,1) + time);

YYYY  = str2num(datestr(time,'yyyy'));
MM    = str2num(datestr(time,'mm'));

% rearrange longitude (same as PlankTOM) and select time
lon = lon([181:end 1:180]);
lon(lon < 0) = lon(lon < 0) + 360; 
SOCAT = SOCAT([181:end 1:180],:,ismember(YYYY, ywanted));
[LAT, LON] = meshgrid(lat,lon);
LAT = repmat(LAT,1,1,size(SOCAT,3));

clearvars -except ywanted YWANTED LON LAT COAST AREA SOCAT

% --- PlankTOM
% read PlankTOM regridded
pfolder = 'data/PlankTOM12/TOM12_RW_RNB0/';
pfiles = dir(pfolder);
pfile = pfiles(contains({pfiles.name},'sfco2')).name;
SFCO2 = ncread([pfolder,pfile],'sfco2');
yyyy = double(ncread([pfolder,pfile],'TIME'));
yyyy = str2num(datestr(round(datenum(1950,1,1) + (yyyy/86400)),'yyyy'));
SFCO2_PlankTOM = SFCO2(:,:,ismember(yyyy, ywanted));

clearvars -except ywanted YWANTED LON LAT COAST AREA SOCAT SFCO2_PlankTOM

% --- Ocean Model from GCB2023
files = dir('data/zenodo_model/*.nc');
files(~cellfun('isempty',regexpi({files.name},'Breakdown'))) = [];
files(~cellfun('isempty',regexpi({files.name},'PlankTOM12'))) = [];
ywanted = 1990:2022; % Years wanted for times series

oceanmodel   = nan(360,180,length(ywanted)*12,length(files));
for m = 1:length(files)
    filename = [files(m).folder,'/',files(m).name];
    time  = datenum(1959,1,1) + double(ncread(filename,'time'));
    sfco2 = double(ncread(filename,'sfco2'));
    yyyy = str2num(datestr(time,'yyyy'));
    
    sfco2 = sfco2(:,:,:,1);
    
    oceanmodel(:,:,:,m) = sfco2(:,:,ismember(yyyy,ywanted));
end

SFCO2_GOBM = oceanmodel;
SFCO2_GOBM(SFCO2_GOBM > 1000) = NaN;
SFCO2_PlankTOM(isnan(repmat(mean(mean(cat(4,SFCO2_GOBM,SFCO2_PlankTOM),3),4),1,1,length(ywanted)*12))) = NaN; % Keep only pixels where every model has a value
SFCO2_GOBM(isnan(repmat(mean(mean(cat(4,SFCO2_GOBM,SFCO2_PlankTOM),3),4),1,1,length(ywanted)*12,size(SFCO2_GOBM,4)))) = NaN; % Keep only pixels where every model has a value


clearvars -except ywanted YWANTED LON LAT COAST AREA SOCAT SFCO2_PlankTOM...
    SFCO2_GOBM


%%

txt_title = {'(a) SOCAT','(b) NEMO-PlankTOM12.1','(c) GOBMs'};
for cpt = 1:5
    lat = LAT(:,:,1)';
    lon = LON(:,:,1)';

    if cpt == 1
        DIFF = nanmean(SOCAT,3)';
        ax = subplot(3,4,[2 3]);
    elseif cpt == 2
        DIFF = nanmean(SFCO2_PlankTOM,3)';
        ax = subplot(3,4,[5 6]);
    elseif cpt == 4
        DIFF = nanmean(SFCO2_PlankTOM - SOCAT,3)';
        ax = subplot(3,4,[7 8]);
    elseif cpt == 3
        DIFF = nanmean(nanmean(SFCO2_GOBM,3),4)';
        ax = subplot(3,4,[9 10]);
    elseif cpt == 5
        DIFF = nanmean(nanmean( SFCO2_GOBM - repmat(SOCAT,1,1,1,size(SFCO2_GOBM,4)) ,3),4)';
        ax = subplot(3,4,[11 12]);
    end
    
    lonw = 40;

    % Rearrange data to lie in the longitude limits I give for the projection
    ind=[lonw+1:360 1:lonw]; % Move left side to right
    DIFF=DIFF(:,ind);
    lat=lat(:,ind);
    lon=lon(:,ind);lon(lon>lonw)=lon(lon>lonw)-360; %...and subtract 360 to some longitudes

    m_proj('robinson','lon',[-360+lonw lonw]);
    m_pcolor(lon,lat,DIFF);
    m_coast('patch',[.5 .5 .5],'edgecolor','k');
    m_grid('fontsize',8,'xtick','');

    if cpt <= 3
        set(gca,'Clim',[280 420])
        col = cbrewer('seq','YlOrRd',7);
        colormap(ax, col)
        hcb = colorbar('southoutside');
        ylabel(hcb,'fCO_2 (μatm)')
        title(txt_title(cpt),'fontweight','normal')
    else
        set(gca,'Clim',[-50 50])
        col = flipud(cbrewer('div','RdBu',10));
        colormap(ax, col)
        hcb = colorbar('southoutside');
        ylabel(hcb,'Bias (μatm)')
        set(hcb,'ytick',-50:10:50)
    end
    
end

% set(gcf,'PaperPosition',[0 0 22 30])
% print('bias_maps.jpeg','-djpeg','-r300')