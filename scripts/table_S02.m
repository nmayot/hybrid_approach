%% Read regional data from NetCDF files (GCB 2023 - zenodo)
clearvars
files = dir('data/zenodo_dataprod/*.nc');
files(~cellfun('isempty',regexpi({files.name},'Watson'))) = [];
ywanted = 1990:2022; % Years wanted for times series

dataprod   = nan(length(ywanted),4,length(files));
for m = 1:length(files)
    filename = [files(m).folder,'/',files(m).name];
    time  = datenum(1959,1,1) + double(ncread(filename,'time'));
    fgco2 = double(ncread(filename,'fgco2_reg'));
    yyyy = str2num(datestr(time,'yyyy'));
    
    if m ~= 3
        fgco2 = reshape(fgco2(ismember(yyyy,ywanted),:),12,[],4);
        fgco2 = permute(mean(fgco2),[2 3 1]);
    elseif m == 3
        fgco2 = fgco2(:,ismember(unique(yyyy),ywanted))';
    end
    
    dataprod(:,:,m) = fgco2;
end

dataprod = permute(dataprod,[1,3,2]);

clearvars -except dataprod ywanted

oceanmodel = ncread('data/zenodo_model/GCB-2023_OceanModel_RegionalBreakdown_1959-2022.nc','fgco2_reg');
oceanmodel = permute(oceanmodel(:,:,[11 2:7 9:10],1),[1,3,2]);
oceanmodel = cat(3,sum(oceanmodel,3),oceanmodel);

yyyy = ncread('data/zenodo_model/GCB-2023_OceanModel_RegionalBreakdown_1959-2022.nc','time');
oceanmodel = oceanmodel(ismember(yyyy,ywanted),:,:);

ywanted = ywanted';
clearvars -except dataprod oceanmodel ywanted

%% read data from PlankTOM and calculate trend
regions = {'Global','North','Tropics','South'};
rwanted = {'Global'};
T    = readtable('data/Cflu_RSS_n_RNB.xlsx','sheet',cell2mat(rwanted));

% ywanted = data{:,1};
GOBM = oceanmodel(:,:,ismember(regions,rwanted));
PROD = dataprod(:,:,ismember(regions,rwanted));
% PERT = T{21:end,[2 4 5]};

yw  = 2010:2019;
alpha = 0.05;

results_gobm = nan(size(GOBM,2),3,1);


for m = 1:size(GOBM,2)

    y = GOBM(:,m)-nanmean(GOBM(:,m));
    datain = [ywanted(ismember(ywanted, yw)),y(ismember(ywanted, yw))];
    [~, ~, h, sig, ~, ~, ~, senD] = ktaub(datain, alpha, 0);
    [ydec2, ydec] = ts_decompo(datain(:,2));
    results_gobm(m,:,1) = [senD*10 std(detrend(ydec2)) std(detrend(ydec))];

end

results_prod = nan(size(PROD,2),3,1);

for m = 1:size(PROD,2)

    y = PROD(:,m)-nanmean(PROD(:,m));
    datain = [ywanted(ismember(ywanted, yw)),y(ismember(ywanted, yw))];
    [~, ~, h, sig, ~, ~, ~, senD] = ktaub(datain, alpha, 0);
    [ydec2, ydec] = ts_decompo(datain(:,2));
    results_prod(m,:,1) = [senD*10 std(detrend(ydec2)) std(detrend(ydec))];

end

y = T.Cflu_Opt_3rd;
y(T.group ~= 0) = NaN;

datain = [T.years(ismember(T.years,yw) & ~isnan(y)), y(ismember(T.years,yw) & ~isnan(y))];
[~, ~, h, sig, ~, ~, ~, senD] = ktaub(datain, alpha, 0);
[ydec2, ydec] = ts_decompo(datain(:,2));
results_hyb = [senD*10 std(detrend(ydec2)) std(detrend(ydec))];

y = T.Cflu_RNB0(ismember(T.years,yw));
datain = [T.years(ismember(T.years,yw)),y];
[~, ~, h, sig, ~, ~, ~, senD] = ktaub(datain, alpha, 0);
[ydec2, ydec] = ts_decompo(datain(:,2));
results_pihm = [senD*10 std(detrend(ydec2)) std(detrend(ydec))];


% results_pert = nan(size(PERT,2),3);
% 
% for m = 1:size(PERT,2)
% 
%     y = PERT(:,m)-nanmean(PERT(:,m));
%     datain = [ywanted(ismember(ywanted, yw)),y(ismember(ywanted, yw))];
%     [~, ~, h, sig, ~, ~, ~, senD] = ktaub(datain, alpha, 0);
%     [ydec2, ydec] = ts_decompo(datain(:,2));
%     results_pert(m,:) = [senD*10 std(detrend(ydec2)) std(detrend(ydec))];
% 
% end

function [ydec2, ydec] = ts_decompo(ts)
    % decadal component
    win = hann(10);
    win = win/sum(win);
    b2 = detrend(ts);
    ydec = conv(b2,win,'same');
    
    % interannual component
    ydec2 = detrend(ts)-ydec;

end