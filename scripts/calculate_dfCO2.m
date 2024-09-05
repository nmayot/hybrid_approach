%% calculate the fugacity coefficient to convert fCO2 (from SOCAT, dataprod, and model) into pCO2
clearvars

cd('/Users/nicolasmayot/Documents/UEA-postdoc/Mirror/redaction/20220505 - Erik/hybrid_approach/data')

% Time and location wanted
ywanted = 1990:2019;
lat_wanted = -89.5:89.5;
lon_wanted = 0.5:359.5;
[lon_wanted,lat_wanted] = ndgrid(lon_wanted,lat_wanted);

% open SST data
sst     = double(ncread('sst-OISST1-2.nc','sst'));
nanval  = ncreadatt('sst-OISST1-2.nc','sst','missing_value');
sst(sst == nanval) = NaN;

time = ncread('sst-OISST1-2.nc','time'); 
time = time + datenum(1800,1,1);
lon_sst = ncread('sst-OISST1-2.nc','lon');
lat_sst = ncread('sst-OISST1-2.nc','lat');

% regrid and subset SST
[x,y]  = ndgrid(lon_sst,lat_sst);
gdata  = griddedInterpolant(x,y,sst);
sst_r  = gdata(lon_wanted,lat_wanted);
sst_r  = sst_r(:,:,ismember(str2num(datestr(time,'yyyy')),ywanted));

% regrid and subset SSS
load('EN4_SSS');
lon_sss = lon;
lat_sss = lat;

[x,y]  = ndgrid(lon_sss,lat_sss);
gdata  = griddedInterpolant(x,y,sss_monthly);
sss_r  = gdata(lon_wanted,lat_wanted);

% open SP data
pres  = ncread('surface_air_pressure-ERA5.nc','sp',[1 1 1 1],[inf inf 1 inf]); % read Atmospheric Pressure
pres = permute(pres,[1,2,4,3]);

time = double(ncread('surface_air_pressure-ERA5.nc','time')); 
time = (time/24) + datenum(1900,1,1);
lon_pres = ncread('surface_air_pressure-ERA5.nc','longitude');
lat_pres = ncread('surface_air_pressure-ERA5.nc','latitude');

% regrid and subset
[x,y]   = ndgrid(lon_pres,flipud(lat_pres));
gdata   = griddedInterpolant(x,y,fliplr(pres));
pres_r  = gdata(lon_wanted,lat_wanted);
pres_r  = pres_r(:,:,ismember(str2num(datestr(time,'yyyy')),ywanted));

% read atmospheric CO2 from GCB22 time series
T = readtable('global_co2_merged_2022.txt');
time = table2array(T(:,1));
xco2 = table2array(T(:,2));

xco2(isnan(time)) = [];
time(isnan(time)) = [];

xco2 = xco2(ismember(floor(time), ywanted));
xco2 = permute(repmat(xco2,1,360,180),[2,3,1]);

% Compute vapor pressure of seawater (atm)
TK = sst_r + 273.15;
% pH2O = exp(24.4543 - 67.4509*(100/TK) - 4.8489*log(TK/100) - 0.000544*sss_r);

cv1 =  -7.85951783; 
cv2 =   1.84408259;
cv3 =  -11.7866497;
cv4 =   22.6807411;
cv5 =  -15.9618719;
cv6 =   1.80122502; 
co0 =   0.90799;   
co1 =  -0.08992;  
co2 =   0.18458; 
co3 =  -0.07395;
co4 =  -0.00221;

ztc = 647.096;
zpc = 22.064 / 101325.0e-6; 

zrt = 1 - TK/ztc;
zrt05 = sqrt(zrt);
zrt15 = zrt.*zrt05;
zrt2  = zrt.*zrt;
zrt3  = zrt.*zrt2;
zrt35 = zrt3.*zrt05;
zrt4  = zrt2.*zrt2;
zrt75 = zrt4.*zrt35;
vpwat  = zpc .* exp((ztc/TK).* (cv1.*zrt + cv2.*zrt15 + cv3.*zrt3 + cv4.*zrt35 + cv5.*zrt4 + cv6.*zrt75) ); % Vapor pressure of water

zsumb = 31.998 * sss_r ./ ( 1e3 - 1.005*sss_r );
zsmh  = zsumb * 0.5;
zsmh2 = zsmh.*zsmh;
zsmh3 = zsmh.*zsmh2;
zsmh4 = zsmh2.*zsmh2;
zosmo = co0 + co1*zsmh + co2*zsmh2 + co3*zsmh3 + co4*zsmh4;
pH2O  = vpwat .* exp(-0.018 .* zosmo .* zsumb); % Vapor pressure of seawater


% Compute pCO2
Patm = pres_r /101324.9966; % conv to atmospheres for this calc
pCO2 = (Patm - pH2O) .* xco2;

% calculate fugacity coefficient
Ptot = Patm; % no hydrostatic pressure added (surface of the ocean)
B    =  -1636.75 + 12.0408*TK - 0.0327957*(TK.^2) + 0.0000316528*(TK.^3);
Del  = 57.7-0.118*TK;
xc2  = (1 - xco2*1e-6).^2;
fugcoeff = exp( Ptot.* (B + 2.*xc2.*Del)./(82.057.*TK) );

fCO2_atm = pCO2 .* fugcoeff;

clearvars -except pres_r sst_r xco2 fCO2_atm fugcoeff lon_wanted lat_wanted ywanted

%% Open dataprod sfco2
cd('zenodo_dataprod/')
files = dir('*.nc');
files(~cellfun('isempty',regexpi({files.name},'Watson'))) = [];

dataprod   = nan(360,180,length(ywanted(1):ywanted(end))*12,length(files)-1);
for m = 1:length(files)
    time  = datenum(1959,1,1) + double(ncread(files(m).name,'time'));
    yyyy = str2num(datestr(time,'yyyy'));
    lat   = double(ncread(files(m).name,'lat'));
    sfco2 = double(ncread(files(m).name,'sfco2'));
    
    sfco2 = sfco2(:,:,ismember(yyyy,ywanted(1):ywanted(end)));
    
    dataprod(:,:,ismember(sort(repmat(ywanted(1):ywanted(end),1,12)),yyyy),m) = sfco2;
end

clearvars -except pres_r sst_r xco2 fCO2_atm fugcoeff lon_wanted lat_wanted ywanted dataprod

%% Open GOBM sfco2
cd('../zenodo_model/')
files = dir('*.nc');
files(~cellfun('isempty',regexpi({files.name},'Regional'))) = [];
files(~cellfun('isempty',regexpi({files.name},'PlankTOM12'))) = [];

datamodel   = nan(360,180,length(ywanted(1):ywanted(end))*12,length(files)-1);
for m = 1:length(files)
    time  = datenum(1959,1,1) + double(ncread(files(m).name,'time'));
    yyyy = str2num(datestr(time,'yyyy'));
    lat   = double(ncread(files(m).name,'lat'));
    sfco2 = double(ncread(files(m).name,'sfco2'));
    
    sfco2 = sfco2(:,:,ismember(yyyy,ywanted(1):ywanted(end)));
    
    datamodel(:,:,ismember(sort(repmat(ywanted(1):ywanted(end),1,12)),yyyy),m) = sfco2;
end

clearvars -except pres_r sst_r xco2 fCO2_atm fugcoeff lon_wanted lat_wanted ywanted dataprod datamodel

%% read SOCAT and regridded sfco2

file = '../SOCATv2023_tracks_gridded_monthly.nc';

SOCAT = ncread(file,'fco2_ave_unwtd');

time  = ncread(file,'tmnth');
lat   = ncread(file,'ylat');
lon   = ncread(file,'xlon');
time  = round(datenum(1970,1,1) + time);

YYYY  = str2num(datestr(time,'yyyy'));
MM    = str2num(datestr(time,'mm'));

% rearrange longitude
lon = lon([181:end 1:180]);
lon(lon < 0) = lon(lon < 0) + 360; 
SOCAT = SOCAT([181:end 1:180],:,ismember(YYYY, ywanted));

SOCAT = SOCAT ./ fugcoeff; % 0.99 = fugacity coefficient

clearvars -except pres_r sst_r xco2 fCO2_atm fugcoeff lon_wanted lat_wanted ywanted dataprod datamodel SOCAT

%% Values: dataprod dfCO2

load('/Users/nicolasmayot/Documents/UEA-postdoc/Mirror/projects/4C/data/ice_extent_fronts/RECCAP2_modified.mat')
ocmask = repmat(ocmask,1,1,size(SOCAT,3));


YYYY = sort(repmat(ywanted,1,12));
MM   = repmat(1:12,1,length(ywanted));
time = datenum(YYYY,MM,1);

rwanted = [2 3];

DATA  = SOCAT;
ATM   = fCO2_atm;

result = [];

for M = 1:size(dataprod,4) % PROD  = nanmean(dataprod,4);
    PROD  = dataprod(:,:,:,M);
    PROD(~ismember(ocmask, rwanted)) = NaN;  
    ATM(~ismember(ocmask, rwanted)) = NaN;   
    DATA(~ismember(ocmask, rwanted)) = NaN;   
    
    PROD2 = PROD;
    ATM2 = ATM;
    PROD2(isnan(DATA)) = NaN;  
    ATM2(isnan(DATA)) = NaN;  
        
    alpha = 0.05;
    ydecs  = [1990:1999; 2000:2009; 2010:2019];
    
    res_int = [];
    for d = 1 % 1 = subsampled; 2 = not-subsampled
        for yd = 1:3
            ydec = ydecs(yd,:);
            
            if d == 1
                y = prctile(reshape(-(PROD2-ATM2),[],length(ywanted)),50);
            else
                y = prctile(reshape(-(PROD-ATM),[],length(ywanted)),50);
            end
        
            datain = [ywanted(ismember(ywanted, ydec))',y(ismember(ywanted, ydec))'];
            [~, ~, h, sig, ~, ~, ~, senD] = ktaub(datain, alpha, 0);
            res_int = [res_int, round(senD*10*100)/100];
        end
    end
    result = [result; res_int];
        
end

%% Values: GOBM dfCO2

load('/Users/nicolasmayot/Documents/UEA-postdoc/Mirror/projects/4C/data/ice_extent_fronts/RECCAP2_modified.mat')
ocmask = repmat(ocmask,1,1,size(SOCAT,3));


YYYY = sort(repmat(ywanted,1,12));
MM   = repmat(1:12,1,length(ywanted));
time = datenum(YYYY,MM,1);

rwanted = [2 3];

DATA  = SOCAT;
ATM   = fCO2_atm;

result = [];

for M = 1:size(datamodel,4) 
    PROD  = datamodel(:,:,:,M);
    PROD(~ismember(ocmask, rwanted)) = NaN;  
    ATM(~ismember(ocmask, rwanted)) = NaN;   
    DATA(~ismember(ocmask, rwanted)) = NaN;   
    
    PROD2 = PROD;
    ATM2 = ATM;
    PROD2(isnan(DATA)) = NaN;  
    ATM2(isnan(DATA)) = NaN;  
        
    alpha = 0.05;
    ydecs  = [1990:1999; 2000:2009; 2010:2019];
    
    res_int = [];
    for d = 1 % 1 = subsampled; 2 = not-subsampled
        for yd = 1:3
            ydec = ydecs(yd,:);
            
            if d == 1
                y = prctile(reshape(-(PROD2-ATM2),[],length(ywanted)),50);
            else
                y = prctile(reshape(-(PROD-ATM),[],length(ywanted)),50);
            end
        
            datain = [ywanted(ismember(ywanted, ydec))',y(ismember(ywanted, ydec))'];
            [~, ~, h, sig, ~, ~, ~, senD] = ktaub(datain, alpha, 0);
            res_int = [res_int, round(senD*10*100)/100];
        end
    end
    result = [result; res_int];
        
end

%% Values: SOCAT dfCO2

load('/Users/nicolasmayot/Documents/UEA-postdoc/Mirror/projects/4C/data/ice_extent_fronts/RECCAP2_modified.mat')
ocmask = repmat(ocmask,1,1,size(SOCAT,3));


YYYY = sort(repmat(ywanted,1,12));
MM   = repmat(1:12,1,length(ywanted));
time = datenum(YYYY,MM,1);

rwanted = [2 3];

DATA  = SOCAT;
ATM   = fCO2_atm;

ATM(~ismember(ocmask, rwanted)) = NaN;   
DATA(~ismember(ocmask, rwanted)) = NaN;   
        
alpha = 0.05;
ydecs  = [1990:1999; 2000:2009; 2010:2019];

result = [];
for yd = 1:3
    ydec = ydecs(yd,:);
    
    y = prctile(reshape(-(SOCAT-ATM),[],length(ywanted)),50);
    
    datain = [ywanted(ismember(ywanted, ydec))',y(ismember(ywanted, ydec))'];
    [~, ~, h, sig, ~, ~, ~, senD] = ktaub(datain, alpha, 0);
    result = [result, round(senD*10*100)/100];
end

        