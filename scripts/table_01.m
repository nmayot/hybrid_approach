%% Table 1 : data from GCB table and hybrid approach results
clearvars
close all

file = 'data/Global_Carbon_Budget_2022v1.0.xlsx';
data = readtable(file,'sheet','Ocean Sink','Range',[31 1]);

ywanted = data{31:62,1};
GOBM = data{31:62,[5:11 13:14]};
PROD = data{31:62,19:25};

results_gobm = nan(size(GOBM,2),3);
yw  = 2000:2019;
alpha = 0.05;

for m = 1:size(GOBM,2)

    y = GOBM(:,m)-nanmean(GOBM(:,m));
    datain = [ywanted(ismember(ywanted, yw)),y(ismember(ywanted, yw))];
    [~, ~, h, sig, ~, ~, ~, senD] = ktaub(datain, alpha, 0);
    [ydec2, ydec] = ts_decompo(y(ismember(ywanted, yw)));
    results_gobm(m,:,1) = [senD*10 std(detrend(ydec2)) std(detrend(ydec))];
        
end

results_prod = nan(size(PROD,2),3);

for m = 1:size(PROD,2)

    y = PROD(:,m)-nanmean(PROD(:,m));
    datain = [ywanted(ismember(ywanted, yw)),y(ismember(ywanted, yw))];
    [~, ~, h, sig, ~, ~, ~, senD] = ktaub(datain, alpha, 0);
    [ydec2, ydec] = ts_decompo(y(ismember(ywanted, yw)));
    results_prod(m,:,1) = [senD*10 std(detrend(ydec2)) std(detrend(ydec))];

end

T    = readtable('data/Cflu_RSS_n_NoCoast.xlsx','sheet','Global');
PIHM = T.Cflu_PIHM;
ywanted = T.years;

results_pihm = nan(1,3);
y = PIHM-nanmean(PIHM);
datain = [ywanted(ismember(ywanted, yw)),y(ismember(ywanted, yw))];
[~, ~, h, sig, ~, ~, ~, senD] = ktaub(datain, alpha, 0);
[ydec2, ydec] = ts_decompo(y(ismember(ywanted, yw)));
results_pihm(1,:) = [senD*10 std(detrend(ydec2)) std(detrend(ydec))];

T    = readtable('data/Cflu_RSS_n_NoCoast.xlsx','sheet','Global');
HYB = T.Cflu_Opt_3rd;
ywanted = T.years;

results_hyb = nan(1,3);
y = HYB-nanmean(HYB);
datain = [ywanted(ismember(ywanted, yw)),y(ismember(ywanted, yw))];
[~, ~, h, sig, ~, ~, ~, senD] = ktaub(datain, alpha, 0);
[ydec2, ydec] = ts_decompo(y(ismember(ywanted, yw)));
results_hyb(1,:) = [senD*10 std(detrend(ydec2)) std(detrend(ydec))];

function [ydec2, ydec] = ts_decompo(ts)
    % decadal component
    win = hann(10);
    win = win/sum(win);
    b2 = detrend(ts);
    ydec = conv(b2,win,'same');
    
    % interannual component
    ydec2 = detrend(ts)-ydec;

end


%% STD from Uniform distribution
yw = 2000:2019;

file = 'data/Cflu_RSS_n_NoCoast.xlsx';
data = readtable(file,'sheet','Global');
YYYY = table2array(data(:,1));
hybrid = table2array(data(:,14:16));
hybrid = hybrid(ismember(YYYY,yw),:);
hybrid(:,2) = hybrid(:,1) - hybrid(:,2);
hybrid(:,3) = hybrid(:,1) + hybrid(:,3);

YYYY = yw;

subset = 10000;

% slope
results = nan(subset,2);
ts = nan(length(yw),subset);
row = find(ismember(YYYY,yw));
for r = 1:subset
    y = hybrid(row,2) + diff(hybrid(row,2:3),1,2).*rand(length(yw),1);
    ts(:,r) = y;
    [~, ~, h, sig, ~, ~, ~, sen] = ktaub([yw' y], 0.05, 0);
    results(r,:) = [sen*10 sig];
end

% A-IV and A-DV
results = nan(size(ts,2),2);
for m = 1:size(ts,2)
    % decadal component
    win = hann(10);
    win = win/sum(win);
    b2 = detrend(ts(:,m));
    ydec = conv(b2,win,'same');
    
    % interannual component
    ydec2 = detrend(ts(:,m))-ydec;

    results(m,:) = [std(ydec2) std(ydec)];
end