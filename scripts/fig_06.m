% Change in fCO2-products over the last 3 GCB updates
clearvars

T = readtable('data/dataprod_GCB_19-23.xlsx','sheet','fco2-products');

yanom = 1990:2020;
yanom = ismember(T.years,yanom);

Y = [T{:,5}-T{:,end}; flipud(T{:,5}+T{:,end})]-[repmat(nanmean(T{yanom,5}),length(T.years),1);flipud(repmat(nanmean(T{yanom,5}),length(T.years),1))];
X = [T.years; flipud(T.years)];
X(isnan(Y)) = [];
Y(isnan(Y)) = [];

yyaxis left
h1 = patch(X,Y,[.75 .75 .75]);
hold on
h2 = plot(T.years,T{:,5}-repmat(nanmean(T{yanom,5}),length(T.years),1),'k-','linewidth',2);
h3 = plot(T.years,T{:,4}-repmat(nanmean(T{yanom,4}),length(T.years),1),'-','color',[.85 .30 .10],'linewidth',2);
h4 = plot(T.years,T{:,3}-repmat(nanmean(T{yanom,3}),length(T.years),1),'-','color',[.45 .45 .70],'linewidth',2);

set(gca,'Xlim',[1990 2023],'Ylim',[-1.2 1.2],'Ygrid','on','Yminorgrid','on','Xminorgrid','on','Xgrid','on','box','on','Layer','top','ycolor','k')
ylabel('CO_2 flux anomaly (PgC yr^{-1})')

%% Change in GOBMs over the last 3 GCB updates
T = readtable('data/dataprod_GCB_19-23.xlsx','sheet','gobms');

yanom = 1990:2020;
yanom = ismember(T.years,yanom);

Y = [T{:,5}-T{:,end}; flipud(T{:,5}+T{:,end})]-[repmat(nanmean(T{yanom,5}),length(T.years),1);flipud(repmat(nanmean(T{yanom,5}),length(T.years),1))];
X = [T.years; flipud(T.years)];
X(isnan(Y)) = [];
Y(isnan(Y)) = [];

yyaxis left
h5 = plot(T.years,T{:,5}-repmat(nanmean(T{yanom,5}),length(T.years),1),'k:','linewidth',2);
h6 = plot(T.years,T{:,4}-repmat(nanmean(T{yanom,4}),length(T.years),1),':','color',[.85 .30 .10],'linewidth',2);
h7 = plot(T.years,T{:,3}-repmat(nanmean(T{yanom,3}),length(T.years),1),':','color',[.45 .45 .70],'linewidth',2);

%% Number of data in SOCAT
yyaxis right

time = ncread('data/SOCATv2023_tracks_gridded_monthly.nc','tmnth');
fco2_count_nobs = ncread('data/SOCATv2023_tracks_gridded_monthly.nc','fco2_count_nobs');

time  = round(time+datenum(1970,1,1))-1;
ytime = str2num(datestr(time,'yyyy'));

fco2_count_nobs(fco2_count_nobs ~= 0) = 1;
A = permute(sum(reshape(sum(reshape(fco2_count_nobs,[],length(ytime))),12,[],length(unique(ytime)))),[3 1 2]);
A(unique(ytime) < 2005) = NaN;

bar(unique(ytime),A/1e4,'FaceColor','k')
set(gca,'Ylim',[0 6])
hold on

%%
time = ncread('data/SOCATv2022_tracks_gridded_monthly.nc','tmnth');
fco2_count_nobs = ncread('data/SOCATv2022_tracks_gridded_monthly.nc','fco2_count_nobs');

time  = round(time+datenum(1970,1,1))-1;
ytime = str2num(datestr(time,'yyyy'));

fco2_count_nobs(fco2_count_nobs ~= 0) = 1;
A = permute(sum(reshape(sum(reshape(fco2_count_nobs,[],length(ytime))),12,[],length(unique(ytime)))),[3 1 2]);
A(unique(ytime) < 2005) = NaN;

bar(unique(ytime),A/1e4,'FaceColor',[.85 .30 .10])

%%
time = ncread('data/SOCATv2021_tracks_gridded_monthly.nc','tmnth');
fco2_count_nobs = ncread('data/SOCATv2021_tracks_gridded_monthly.nc','fco2_count_nobs');

time  = round(time+datenum(1970,1,1))-1;
ytime = str2num(datestr(time,'yyyy'));

fco2_count_nobs(fco2_count_nobs ~= 0) = 1;
A = permute(sum(reshape(sum(reshape(fco2_count_nobs,[],length(ytime))),12,[],length(unique(ytime)))),[3 1 2]);
A(unique(ytime) < 2005) = NaN;

bar(unique(ytime),A/1e4,'FaceColor',[.45 .45 .70])

%%
set(gca,'ycolor','k','Ytick',[0 .5 1.75])
label_h = ylabel({'Monthly grid cells (1° x 1°)';'per year (x10^{4})'});
label_h.Position(2) = 0;
label_h.Position(2) = 1;
legend([h1 h2 h3 h4 h5 h6 h7],{'1\sigma 2023 fCO_2-products','2023 fCO_2-products','2022 fCO_2-products','2021 fCO_2-products',...
    '2023 GOBMs','2022 GOBMs','2021 GOBMs'},'location','NorthWest','NumColumns',2)


% set(gcf,'PaperPosition',[0 0 22 13.5])
% print('figSXX_last_GCB.jpeg','-djpeg','-r300')
