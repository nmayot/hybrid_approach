%% With GCB Table
clearvars
close all

file = 'data/Global_Carbon_Budget_2023v1.0.xlsx';
data = readtable(file,'sheet','Ocean Sink','Range',[31 1]);

ywanted = data{31:64,1};
GOBM = data{31:64,[5:11 13:14]};
PIHM = data{31:64,12};
PROD = data{31:64,19:25};

results_gobm = nan(size(GOBM,2),3,3);
yw  = 1990:2022;
alpha = 0.05;
TS    = nan(length(yw),9,5);
YDEC2 = nan(length(yw),9,5);
YDEC  = nan(length(yw),9,5);

for m = 1:size(GOBM,2)

    y = GOBM(:,m)-nanmean(GOBM(1:10,m));
    datain = [ywanted(ismember(ywanted, yw)),y(ismember(ywanted, yw))];
    [~, ~, h, sig, ~, ~, ~, senD] = ktaub(datain, alpha, 0);
    [ydec2, ydec] = ts_decompo(y(ismember(ywanted, yw)));
    results_gobm(m,:,1) = [senD*10 std(detrend(ydec2)) std(detrend(ydec))];
    TS(:,m,1) = datain(:,2);
    YDEC(:,m,1) = ydec;
    YDEC2(:,m,1) = ydec2;
    
    
end

results_prod = nan(size(PROD,2),3,3);

for m = 1:size(PROD,2)

    y = PROD(:,m)-nanmean(PROD(1:10,m));
    datain = [ywanted(ismember(ywanted, yw)),y(ismember(ywanted, yw))];
    [~, ~, h, sig, ~, ~, ~, senD] = ktaub(datain, alpha, 0);
    [ydec2, ydec] = ts_decompo(y(ismember(ywanted, yw)));
    results_prod(m,:,1) = [senD*10 std(detrend(ydec2)) std(detrend(ydec))];
    TS(:,m,2) = datain(:,2);
    YDEC(:,m,2) = ydec;
    YDEC2(:,m,2) = ydec2;

end

T    = readtable('../Cflu_RSS_n_RNB.xlsx','sheet','Global');
PIHM = T.Cflu_RNB0;
ywanted = T.years;

results_pihm = nan(3,2);
y = PIHM-nanmean(PIHM(21:30));
datain = [ywanted(ismember(ywanted, yw)),y(ismember(ywanted, yw))];
[~, ~, h, sig, ~, ~, ~, senD] = ktaub(datain, alpha, 0);
results_pihm(1,:) = [senD*10 std(detrend(datain(:,2)))];
[ydec2, ydec] = ts_decompo(y(ismember(ywanted, yw)));
TS(:,1,3) = datain(:,2);
YDEC(:,1,3) = ydec;
YDEC2(:,1,3) = ydec2;

T    = readtable('../Cflu_RSS_n_RNB.xlsx','sheet','Global');
PIHM = T.Cflu_Opt_3rd;
ywanted = T.years;

yPIHM = ywanted(~isnan(PIHM));
PIHM = interp1(yPIHM,PIHM(~isnan(PIHM)),ywanted);

results_hyb = nan(3,2);
y = PIHM-nanmean(PIHM(21:30));
datain = [ywanted(ismember(ywanted, yw)),y(ismember(ywanted, yw))];
[~, ~, h, sig, ~, ~, ~, senD] = ktaub(datain, alpha, 0);
results_hyb(1,:) = [senD*10 std(detrend(datain(:,2)))];
[ydec2, ydec] = ts_decompo(y(ismember(ywanted, yw)));
TS(:,1,4) = datain(:,2);
YDEC(:,1,4) = ydec;
YDEC2(:,1,4) = ydec2;

% T    = readtable('../Cflu_RSS_n_RNB.xlsx','sheet','Global');
% PIHM = T.Cflu_Opt_3rd_3y;
% ywanted = T.years;
% 
% results_hyb = nan(3,2);
% y = PIHM-nanmean(PIHM(21:30));
% datain = [ywanted(ismember(ywanted, yw)),y(ismember(ywanted, yw))];
% [~, ~, h, sig, ~, ~, ~, senD] = ktaub(datain, alpha, 0);
% results_hyb(1,:) = [senD*10 std(detrend(datain(:,2)))];
% [ydec2, ydec] = ts_decompo(y(ismember(ywanted, yw(1:end-1))));
% TS(:,1,5) = datain(:,2);
% YDEC(1:end-1,1,5) = ydec;
% YDEC2(1:end-1,1,5) = ydec2;

% Figure
subplot(3,1,1)
h1 = plot(yw,mean(TS(:,:,1),2,'omitnan'),'-','color',[.5 .5 .5],'linewidth',2);
hold on
h2 = plot(yw,mean(TS(:,:,2),2,'omitnan'),'-','color',[.5 .6 .8],'linewidth',2);
h3 = plot(yw,TS(:,1,3),'k-','linewidth',1.5);
h4 = plot(yw,TS(:,1,4),'o-','color',[.9 .1 .1],'MarkerEdgeColor',[.9 .1 .1],'MarkerFaceColor',[1 .7 .7]);
h5 = plot(yw,TS(:,1,5),'--','color',[.9 .1 .1],'linewidth',2);
set(gca,'Xtick',1990:10:2022,'Xgrid','on','Ygrid','on','Xlim',[1990 2022],'Ylim',[-1 1.5],'XMinorTick','on','XMinorGrid','on','box','on','Layer','top')
ylabel('CO_2 flux anomaly (PgC yr^{-1})')
title("(a) Ocean CO_2 sink (global annual flux)",'fontweight','normal')
legend([h2,h1,h3,h4,h5],{'pCO_2-products','GOBMs','PlankTOM12.1','Hybrid approach','Hybrid approach (3 years)'},'location','NorthWest','NumColumns',2)

subplot(3,1,2)
h1 = plot(yw,mean(YDEC2(:,:,1),2,'omitnan'),'-','color',[.5 .5 .5],'linewidth',2);
hold on
h2 = plot(yw,mean(YDEC2(:,:,2),2,'omitnan'),'-','color',[.5 .6 .8],'linewidth',2);
h3 = plot(yw,YDEC2(:,1,3),'k-','linewidth',1.5);
h4 = plot(yw,YDEC2(:,1,4),'o-','color',[.9 .1 .1],'MarkerEdgeColor',[.9 .1 .1],'MarkerFaceColor',[1 .7 .7]);
h5 = plot(yw,YDEC2(:,1,5),'--','color',[.9 .1 .1],'linewidth',2);
set(gca,'Xtick',1990:10:2022,'Xgrid','on','Ygrid','on','Xlim',[1990 2022],'Ylim',[-.4 .4],'XMinorTick','on','XMinorGrid','on','box','on','Layer','top')
ylabel('CO_2 flux anomaly (PgC yr^{-1})')
title("(b) Interannual component of the ocean CO_2 sink variability",'fontweight','normal')

subplot(3,1,3)
h1 = plot(yw,mean(YDEC(:,:,1),2,'omitnan'),'-','color',[.5 .5 .5],'linewidth',2);
hold on
h2 = plot(yw,mean(YDEC(:,:,2),2,'omitnan'),'-','color',[.5 .6 .8],'linewidth',2);
h3 = plot(yw,YDEC(:,1,3),'k-','linewidth',1.5);
h4 = plot(yw,YDEC(:,1,4),'o-','color',[.9 .1 .1],'MarkerEdgeColor',[.9 .1 .1],'MarkerFaceColor',[1 .7 .7]);
h5 = plot(yw,YDEC(:,1,5),'--','color',[.9 .1 .1],'linewidth',2);
set(gca,'Xtick',1990:10:2022,'Xgrid','on','Ygrid','on','Xlim',[1990 2022],'Ylim',[-.3 .3],'XMinorTick','on','XMinorGrid','on','box','on','Layer','top')
ylabel('CO_2 flux anomaly (PgC yr^{-1})')
title("(c) Decadal component of the ocean CO_2 sink variability",'fontweight','normal')

% set(gcf,'PaperPosition',[1 1 16 34])
% print('/Users/nicolasmayot/Documents/UEA-postdoc/Mirror/redaction/20220505 - Erik/FigS1.jpeg','-djpeg','-r300')

function [ydec2, ydec] = ts_decompo(ts)
    % decadal component
    win = hann(10);
    win = win/sum(win);
    b2 = detrend(ts);
    ydec = conv(b2,win,'same');
    
    % interannual component
    ydec2 = detrend(ts)-ydec;

end

% [mean(results_prod(:,1,1)) mean([abs(mean(results_prod(:,1,2)) - mean(results_prod(:,1,1))) abs(mean(results_prod(:,1,3)) - mean(results_prod(:,1,1)))]);
%     mean(results_gobm(:,1,1)) mean([abs(mean(results_gobm(:,1,2)) - mean(results_gobm(:,1,1))) abs(mean(results_gobm(:,1,3)) - mean(results_gobm(:,1,1)))]);
%     results_pihm(1,1) mean([abs(mean(results_pihm(2,1)) - mean(results_pihm(1,1))) abs(mean(results_pihm(3,1)) - mean(results_pihm(1,1)))]);
%     results_hyb(1,1) mean([abs(mean(results_hyb(2,1)) - mean(results_hyb(1,1))) abs(mean(results_hyb(3,1)) - mean(results_hyb(1,1)))])]

