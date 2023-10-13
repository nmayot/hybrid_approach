%% To produce Fig. 4.:
%           - calculate fugacity coeff, and convert xCO2 into fCO2_atm
%           - read fCO2_sea data from fCO2-products
%           - read fCO2_sea data from SOCAT

%% calculate the fugacity coefficient
%
%

clearvars

% Time and location wanted
ywanted = 1990:2019;
lat_wanted = -89.5:89.5;
lon_wanted = 0.5:359.5;
[lon_wanted,lat_wanted] = ndgrid(lon_wanted,lat_wanted);

% open SST data
sst     = double(ncread('data/sst-OISST1-2.nc','sst'));
nanval  = ncreadatt('data/sst-OISST1-2.nc','sst','missing_value');
sst(sst == nanval) = NaN;

time = ncread('data/sst-OISST1-2.nc','time'); 
time = time + datenum(1800,1,1);
lon_sst = ncread('data/sst-OISST1-2.nc','lon');
lat_sst = ncread('data/sst-OISST1-2.nc','lat');

% regrid and subset SST
[x,y]  = ndgrid(lon_sst,lat_sst);
gdata  = griddedInterpolant(x,y,sst);
sst_r  = gdata(lon_wanted,lat_wanted);
sst_r  = sst_r(:,:,ismember(str2num(datestr(time,'yyyy')),ywanted));

% regrid and subset SSS
load('data/EN4_SSS');
lon_sss = lon;
lat_sss = lat;

[x,y]  = ndgrid(lon_sss,lat_sss);
gdata  = griddedInterpolant(x,y,sss_monthly);
sss_r  = gdata(lon_wanted,lat_wanted);

% open SP data
pres  = ncread('data/surface_air_pressure-ERA5.nc','sp',[1 1 1 1],[inf inf 1 inf]); % read Atmospheric Pressure
pres = permute(pres,[1,2,4,3]);

time = double(ncread('data/surface_air_pressure-ERA5.nc','time')); 
time = (time/24) + datenum(1900,1,1);
lon_pres = ncread('data/surface_air_pressure-ERA5.nc','longitude');
lat_pres = ncread('data/surface_air_pressure-ERA5.nc','latitude');

% regrid and subset
[x,y]   = ndgrid(lon_pres,flipud(lat_pres));
gdata   = griddedInterpolant(x,y,fliplr(pres));
pres_r  = gdata(lon_wanted,lat_wanted);
pres_r  = pres_r(:,:,ismember(str2num(datestr(time,'yyyy')),ywanted));

% read atmospheric CO2 from GCB22 time series
T = readtable('data/global_co2_merged_2022.txt');
time = table2array(T(:,1));
xco2 = table2array(T(:,2));

xco2(isnan(time)) = [];
time(isnan(time)) = [];

xco2 = xco2(ismember(floor(time), ywanted));
xco2 = permute(repmat(xco2,1,360,180),[2,3,1]);

% Compute vapor pressure of seawater (atm)
TK = sst_r + 273.15;

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
files = dir('data/zenodo_dataprod/*.nc');

dataprod   = nan(360,180,length(ywanted(1):ywanted(end))*12,length(files)-1);
for m = 1:length(files)-1
    file = [files(m).folder,'/',files(m).name];
    time  = datenum(1959,1,1) + double(ncread(file,'time'));
    yyyy = str2num(datestr(time,'yyyy'));
    lat   = double(ncread(file,'lat'));
    sfco2 = double(ncread(file,'sfco2'));
    
    sfco2 = sfco2(:,:,ismember(yyyy,ywanted(1):ywanted(end)));
    
    dataprod(:,:,ismember(sort(repmat(ywanted(1):ywanted(end),1,12)),yyyy),m) = sfco2;
end

clearvars -except pres_r sst_r xco2 fCO2_atm fugcoeff lon_wanted lat_wanted ywanted dataprod

%% read SOCAT and regridded sfco2

file = 'data/SOCATv2022_tracks_gridded_monthly.nc';

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

clearvars -except pres_r sst_r xco2 fCO2_atm fugcoeff lon_wanted lat_wanted ywanted dataprod SOCAT

%% Figures: SOCAT vs PROD pCO2

% open mask
coast  = double(ncread('data/RECCAP2_region_masks_all_v20210412.nc','coast'));
load('data/RECCAP2_modified.mat')
ocmask(coast == 1) = NaN;
ocmask = repmat(ocmask,1,1,size(SOCAT,3));


YYYY = sort(repmat(ywanted,1,12));
MM   = repmat(1:12,1,length(ywanted));
time = datenum(YYYY,MM,1);

rwanted = [2 3];

DATA  = SOCAT;
ATM   = fCO2_atm;

PROD  = nanmean(dataprod,4);
DATA(~ismember(ocmask, rwanted)) = NaN;   
PROD(~ismember(ocmask, rwanted)) = NaN;  
ATM(~ismember(ocmask, rwanted)) = NaN;   
PROD(isnan(DATA)) = NaN;  
ATM(isnan(DATA)) = NaN;  


boxplot(reshape(-(DATA-ATM),[],length(ywanted)),'PlotStyle','compact','Symbol','','Whisker',0,'Positions',ywanted,'Labels',ywanted,'Colors',[.8 .9 .6])
hold on
h1 = plot(ywanted,prctile(reshape(-(DATA-ATM),[],length(ywanted)),50),'-','linewidth',2,'color',[.2 .6 .2]);
set(gca,'Xtick',1970:10:2020,'XTickLabel',1970:10:2020,'Xgrid','on','Ygrid','on','Xlim',[1969 2020],'Ylim',[-30 70],'XMinorTick','on','XMinorGrid','on','box','on','Layer','top')

PROD  = mean(dataprod,4,'omitnan');
PROD(~ismember(ocmask, rwanted)) = NaN;
PROD(isnan(DATA)) = NaN;
plot(ywanted,prctile(reshape(-(PROD-ATM),[],length(ywanted)),75),'-','color',[.5 .6 .8]);
plot(ywanted,prctile(reshape(-(PROD-ATM),[],length(ywanted)),25),'-','color',[.5 .6 .8]);
h2 = plot(ywanted,prctile(reshape(-(PROD-ATM),[],length(ywanted)),50),'-','color',[.5 .6 .8],'linewidth',2);

alpha = 0.05;
ydecs  = [2000:2009; 2010:2019];

for d = 1:2
    for yd = 1:2
        ydec = ydecs(yd,:);
        
        if d == 1
            y = prctile(reshape(-(DATA-ATM),[],length(ywanted)),50);
            col = [.2 .6 .2];
            xs = 3;
        else
            y = prctile(reshape(-(PROD-ATM),[],length(ywanted)),50);
            col = [.5 .6 .8];
            xs = 1;
        end
    
        datain = [ywanted(ismember(ywanted, ydec))',y(ismember(ywanted, ydec))'];
        [~, ~, h, sig, ~, ~, ~, senD] = ktaub(datain, alpha, 0);
        disp([h sig senD*10])
    
        if h == 1
            vv = median(datain(:,2));
            middata = datain(round(length(datain)/2),1);
            y = vv + senD*(datain(:,1)-middata);
    
            plot(datain(:,1),y,':','color',col,'linewidth',1.5)
        end
    end
end

legend([h1,h2],{'SOCAT','fCO_2-products'},'location','NorthWest')
ylabel({'\DeltafCO_2 in the North,';'at SOCAT sampling points (μatm)'})