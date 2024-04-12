%% Trend with dataprod
clearvars
cd('data/zenodo_dataprod/')
files = dir('*.nc');
files(~cellfun('isempty',regexpi({files.name},'Watson'))) = [];

% open mask
load('../RECCAP2_modified.mat')
ocmask(ocmask > 3 | ocmask == 1) = NaN;


YYYY = 2010:2019;
FMAP = [];
alpha = 0.05;

for m = 1:length(files)
    time  = datenum(1959,1,1) + double(ncread(files(m).name,'time'));
    fgco2 = double(ncread(files(m).name,'fgco2'));
    ywanted = str2num(datestr(time(1),'yyyy')):str2num(datestr(time(end),'yyyy'));

    FGOBM = reshape(fgco2,[],size(fgco2,3));
    FGOBM = FGOBM';
    FGOBM = reshape(FGOBM,12,size(fgco2,3)/12,[]);
    FGOBM = mean(FGOBM);
    FGOBM = permute(FGOBM,[3,2,1]);
    
    F   = FGOBM *(86400*365)*12.0107;
    if m == 3
        F = -F;
    end

    idx = find(~isnan(mean(F(:,ismember(ywanted,YYYY)),2)));
    idy = find(ismember(ywanted,YYYY));
    
    DecTrend = nan(size(F,1),3);
    for px = idx'
        datain = [YYYY',F(px,idy)'];
        
        x = datain(:,1);
        y = datain(:,2);
        [N, c]=size(datain(:,1));
        Comb = combnk(1:N,2);
        deltay=diff(y(Comb),1,2);
        deltax=diff(x(Comb),1,2);
        theil=diff(y(Comb),1,2)./diff(x(Comb),1,2);
        senD=median(theil);
        DecTrend(px,1) = senD*10;

    end

    Fmap = DecTrend(:,1);
    Fmap = reshape(Fmap,360,180);
    Fmap(isnan(ocmask)) = NaN;
    FMAP = cat(3,FMAP,Fmap);
end

clearvars -except FMAP ocmask


%% Number of data in SOCAT
cd('../../')
YYYY = 2010:2019;
lon  = ncread('data/SOCATv2023_tracks_gridded_monthly.nc','xlon');
lat  = ncread('data/SOCATv2023_tracks_gridded_monthly.nc','ylat');
time = ncread('data/SOCATv2023_tracks_gridded_monthly.nc','tmnth');
fco2_count_nobs = ncread('data/SOCATv2023_tracks_gridded_monthly.nc','fco2_count_nobs');

time  = round(time+datenum(1970,1,1))-1;
ytime = str2num(datestr(time,'yyyy'));

idt   = ismember(ytime,YYYY);

SOCAT = reshape(fco2_count_nobs,[],size(fco2_count_nobs,3));
SOCAT(SOCAT ~= 0) = 1;
SOCAT = reshape(SOCAT',12,[],size(SOCAT,1));
SOCAT = permute(sum(SOCAT),[3,2,1]);
SOCAT = SOCAT(:,ismember(unique(ytime),YYYY));
SOCAT = reshape(SOCAT,360,180,size(SOCAT,2));
SOCAT = SOCAT([181:end 1:180],:,:);


[LG,LT]    = meshgrid(lon([181:end 1:180]),lat);
LG(LG < 0) = LG(LG < 0) + 360;

for f = 2
    if f == 1
        B = mean(SOCAT,3)';
    elseif f == 2
        B = median(SOCAT,3)';
    end
    B(isnan(ocmask')) = NaN;
end

clearvars -except B LG LT FMAP ocmask

%%
ax = subplot(1,2,1);
Fmap = mean(FMAP(:,:,[2 5 6 7]),3)-mean(FMAP(:,:,1:7),3);

m_proj('stereographic','lat',90,'long',0,'radius',60);
m_pcolor(LG-.5,LT-.5,Fmap');
m_coast('patch',[.7 .7 .7],'edgecolor','none');
m_grid('tickdir','out','ytick',40:20:80,'yticklabels',[],'xtick',[-120 -60 0 60 120]);

col = flipud(cbrewer('div','RdBu',10));
colormap(ax,col)
set(gca,'Clim',[-10 10])
cmb = colorbar();
ylabel(cmb,'CO_2 flux trend (gC m^{-2} yr^{-1} decade^{-1})')
title("(a) Decadal trend in the 2010s",'fontweight','normal')


%% Maps with Trend and number SOCAT
ax = subplot(1,2,2);
Fmap = B;

m_proj('stereographic','lat',90,'long',0,'radius',60);
m_pcolor(LG-.5,LT-.5,Fmap);
m_coast('patch',[.7 .7 .7],'edgecolor','none');
m_grid('tickdir','out','ytick',40:20:80,'yticklabels',[],'xtick',[-120 -60 0 60 120]);

col = cbrewer('seq','Greens',13);
col(1,:) = [1 1 1];
col(2:4,:) = repmat(col(4,:),3,1);
col(5:7,:) = repmat(col(7,:),3,1);
col(8:10,:) = repmat(col(10,:),3,1);
col(11:13,:) = repmat(col(13,:),3,1);

colormap(ax,col)
set(gca,'Clim',[-0.5 12.5])
cmb = colorbar();
set(cmb,'Ytick',0:12)
ylabel(cmb,'Nb. of months with observations each year')
title("(b) SOCAT observations in the 2010s",'fontweight','normal')
