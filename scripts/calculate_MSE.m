%% Read fgco2

ywanted = 1970:2022;
fgco2 = nan(length(ywanted),4);

% read global and regional PlankTOM fgco2
pfolder = 'data/PlankTOM12/TOM12_RW_RVY4/';
pfiles = dir(pfolder);
pfile = pfiles(contains({pfiles.name},'integrated_timelines')).name;
fgco2(:,1) = mean(reshape(ncread([pfolder,pfile],'fgco2_glob'),12,[]))';

fgco2_reg = ncread([pfolder,pfile],'fgco2_reg')';
fgco2(:,2:4) = permute(mean(reshape(fgco2_reg(:,[3 2 1]),12,[],3)),[2,3,1]);

%% Calculate MSE globally and for a particular region, with the coast
clearvars

ywanted = 1980:2022;

% read PlankTOM regridded
pfolder = 'data/PlankTOM12/TOM12_RW_RNB0/';
pfiles = dir(pfolder);
pfile = pfiles(contains({pfiles.name},'sfco2')).name;
SFCO2 = ncread([pfolder,pfile],'sfco2');
yyyy = double(ncread([pfolder,pfile],'TIME'));
yyyy = str2num(datestr(round(datenum(1950,1,1) + (yyyy/86400)),'yyyy'));
SFCO2 = SFCO2(:,:,ismember(yyyy, ywanted));

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

% Calculate MSE
[LT, ~]    = meshgrid(double(lat),double(lon));
LT = repmat(LT,1,1,size(SOCAT,3));
DIFF = SOCAT - SFCO2;
YYYY      = sort(repmat(ywanted,1,12));

MSE = nan(length(ywanted),4);
N   = nan(length(ywanted),4);
for r = 1:4
    % ---- Select region
    diff = DIFF;  
    if r == 2
        diff(LT <= 30) = NaN;  
    elseif r == 3
        diff(LT > 30 | LT < -30) = NaN;  
    elseif r == 4
        diff(LT >= -30) = NaN;  
    end
    
    % ---- Select time
    for Y = ywanted
        value = diff(:,:,ismember(YYYY,Y));
        value(isnan(value)) = [];
        MSE(Y == ywanted, r) = mean(value.^2);
        N(Y == ywanted, r)   = length(value);
    end
end

%% Calculate MSE globally and for a particular region, without the coast
clearvars

ywanted = 1970:2022;

% read PlankTOM regridded
pfolder = 'data/PlankTOM12/TOM12_RW_RNB4/';
pfiles = dir(pfolder);
pfile = pfiles(contains({pfiles.name},'sfco2')).name;
SFCO2 = ncread([pfolder,pfile],'sfco2');
yyyy = double(ncread([pfolder,pfile],'TIME'));
yyyy = str2num(datestr(round(datenum(1950,1,1) + (yyyy/86400)),'yyyy'));
SFCO2 = SFCO2(:,:,ismember(yyyy, ywanted));

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

% open mask
coast  = double(ncread('~/Documents/UEA-postdoc/Mirror/projects/4C/data/ice_extent_fronts/RECCAP2_region_masks_all_v20210412.nc','coast'));
% load('~/Documents/UEA-postdoc/Mirror/projects/4C/data/ice_extent_fronts/coast_fco2_products.mat')
% coast = coast2;
COAST  = repmat(coast,1,1,size(SOCAT,3));

% Calculate MSE
[LT, ~]    = meshgrid(double(lat),double(lon));
LT = repmat(LT,1,1,size(SOCAT,3));
DIFF = SOCAT - SFCO2;
YYYY      = sort(repmat(ywanted,1,12));

MSE = nan(length(ywanted),4);
N   = nan(length(ywanted),4);
for r = 1:4
    % ---- Select region
    diff = DIFF;  
    if r == 2
        diff(LT <= 30) = NaN;  
    elseif r == 3
        diff(LT > 30 | LT < -30) = NaN;  
    elseif r == 4
        diff(LT >= -30) = NaN;  
    end
    
    diff(COAST == 1) = NaN;
    
    % ---- Select time
    for Y = ywanted
        value = diff(:,:,ismember(YYYY,Y));
        value(isnan(value)) = [];
        MSE(Y == ywanted, r) = mean(value.^2);
        N(Y == ywanted, r)   = length(value);
    end
end


%% Calculate MSE by 5° latitude, with coast
clearvars
file = 'data/SOCATv2023_tracks_gridded_monthly.nc';

SOCAT = ncread(file,'fco2_ave_unwtd');
time  = ncread(file,'tmnth');
lat   = ncread(file,'ylat');
lon   = ncread(file,'xlon');
time  = round(datenum(1970,1,1) + time);

YYYY  = str2num(datestr(time,'yyyy'));
MM    = str2num(datestr(time,'mm'));

% rearrange longitude (same as PlankTOM) and select time
[~, ~, grp] = histcounts(lat,-90:5:90);
lwanted = unique(grp);
ywanted = 2000:2022;
lon = lon([181:end 1:180]);
lon(lon < 0) = lon(lon < 0) + 360; 
SOCAT = SOCAT([181:end 1:180],:,ismember(YYYY, ywanted));

area = double(ncread('/Users/nicolasmayot/Documents/UEA-postdoc/Mirror/projects/4C/data/GOBM/plankTOM/AncillaryData.nc','AREA'));
area = repmat(area,1,1,size(SOCAT,3));

modelruns = {'TOM12_RW_RNB1','TOM12_RW_RNA0','TOM12_RW_RNB0','TOM12_RW_RNB2','TOM12_RW_RNB4'}; 
MSE  = nan(length(lwanted),length(modelruns));
CFLX = nan(length(lwanted),length(modelruns));
N    = nan(length(lwanted),length(modelruns));

for M = 1:length(modelruns)
    % read PlankTOM regridded
    pfolder = ['data/PlankTOM12/',cell2mat(modelruns(M)),'/'];
    pfiles = dir(pfolder);

    pfile = pfiles(contains({pfiles.name},'sfco2')).name;
    SFCO2 = ncread([pfolder,pfile],'sfco2');
    
    pfile = pfiles(contains({pfiles.name},'fgco2')).name;
    FGCO2 = ncread([pfolder,pfile],'fgco2');

    yyyy = double(ncread([pfolder,pfile],'TIME'));
    yyyy = str2num(datestr(datenum(1950,1,1) + (yyyy/86400),'yyyy'));
    SFCO2 = SFCO2(:,:,ismember(yyyy, ywanted));
    FGCO2 = FGCO2(:,:,ismember(yyyy, ywanted));

    DIFF = SOCAT - SFCO2;
    
    cflx = FGCO2;
    cflx = cflx .* area .* (86400*365*12.0107) /10^15;
    cflx = mean(cflx,3);

    for L = lwanted'
        value = DIFF(:,ismember(grp,L),:);
        value(isnan(value)) = [];
        MSE(L,M) = mean(value.^2);
        N(L,M)   = length(value);
        CFLX(L,M) = sum(cflx(:,ismember(grp,L)),[1 2],'omitnan')/5;
    end
end

%% Calculate MSE by 5° latitude, without coast
clearvars
file = '/Users/nicolasmayot/Documents/UEA-postdoc/Mirror/projects/GCB2022/GCBocean/data/cache/SOCATv2022_tracks_gridded_monthly.nc.zip.unzip/SOCATv2022_tracks_gridded_monthly.nc';

SOCAT = ncread(file,'fco2_ave_unwtd');
time  = ncread(file,'tmnth');
lat   = ncread(file,'ylat');
lon   = ncread(file,'xlon');
time  = round(datenum(1970,1,1) + time);

YYYY  = str2num(datestr(time,'yyyy'));
MM    = str2num(datestr(time,'mm'));

% rearrange longitude (same as PlankTOM) and select time
[~, ~, grp] = histcounts(lat,-90:5:90);
lwanted = unique(grp);
ywanted = 1970:2019;
lon = lon([181:end 1:180]);
lon(lon < 0) = lon(lon < 0) + 360; 
SOCAT = SOCAT([181:end 1:180],:,ismember(YYYY, ywanted));

area = double(ncread('/Users/nicolasmayot/Documents/UEA-postdoc/Mirror/projects/4C/data/GOBM/plankTOM/AncillaryData.nc','AREA'));
area = repmat(area,1,1,size(SOCAT,3));
coast  = double(ncread('~/Documents/UEA-postdoc/Mirror/projects/4C/data/ice_extent_fronts/RECCAP2_region_masks_all_v20210412.nc','coast'));
COAST  = repmat(coast,1,1,size(SOCAT,3));

modelruns = {'TOM12_ET_PIHQ','TOM12_ET_PIHM','TOM12_ET_PIHR','TOM12_ET_PIHT'}; 
MSE  = nan(length(lwanted),length(modelruns));
CFLX = nan(length(lwanted),length(modelruns));
N    = nan(length(lwanted),length(modelruns));

for M = 1:length(modelruns)
    % read PlankTOM regridded
    load(['/Users/nicolasmayot/Documents/UEA-postdoc/Mirror/projects/4C/data/GOBM/plankTOM/',cell2mat(modelruns(M)),'/',cell2mat(modelruns(M)),'_SFCO2.mat'])
    load(['/Users/nicolasmayot/Documents/UEA-postdoc/Mirror/projects/4C/data/GOBM/plankTOM/',cell2mat(modelruns(M)),'/',cell2mat(modelruns(M)),'_Cflx.mat'])

    DIFF = SOCAT - SFCO2;
    DIFF(COAST == 1) = NaN;
    
    cflx = FGCO2;
    cflx = cflx .* area .* (86400*365*12.0107);
    cflx = mean(cflx,3);

    for L = lwanted'
        value = DIFF(:,ismember(grp,L),:);
        value(isnan(value)) = [];
        MSE(L,M) = mean(value.^2);
        N(L,M)   = length(value);
        CFLX(L,M) = sum(cflx(:,ismember(grp,L)),[1 2],'omitnan');
    end
end