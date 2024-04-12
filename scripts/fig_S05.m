yyyy = 1970:2022;

T = readtable('../Cflu_RSS_n_RNB.xlsx','sheet','Global'); % North, Tropics, South, latitude
X = [T.Cflu_RNB1 T.Cflu_RNA0 T.Cflu_RNB0 T.Cflu_RNB2 T.Cflu_RNB4];
subplot(3,2,1)
plot(yyyy,X,'k-')
hold on
plot(yyyy,X(:,3),'k-','LineWidth',2)
set(gca,'Xtick',1970:10:2020,'Xgrid','on','Ygrid','on','Xlim',[1972 2025],'Ylim',[-.5 4],'XMinorTick','on','XMinorGrid','on','box','on','Layer','top')
title('(a) NCEP and bact. half-saturation','fontweight','normal')
ylabel('CO_2 flux (PgC yr^{-1})')

T = readtable('../Cflu_RSS_n_RVB.xlsx','sheet','Global'); % North, Tropics, South, latitude
X = [T.Cflu_RVB0 T.Cflu_RVB1 T.Cflu_RVA0 T.Cflu_RVB2 T.Cflu_RVB3];
subplot(3,2,2)
plot(yyyy,X,'k-')
hold on
plot(yyyy,X(:,3),'k-','LineWidth',2)
set(gca,'Xtick',1970:10:2020,'Xgrid','on','Ygrid','on','Xlim',[1972 2025],'Ylim',[-.5 4],'XMinorTick','on','XMinorGrid','on','box','on','Layer','top')
title('(b) ERA5 and bact. half-saturation','fontweight','normal')
ylabel('CO_2 flux (PgC yr^{-1})')

T = readtable('../Cflu_RSS_n_RNY.xlsx','sheet','Global'); % North, Tropics, South, latitude
X = [T.Cflu_RNY1 T.Cflu_RNY2 T.Cflu_RNB0 T.Cflu_RNY4 T.Cflu_RNY3];
subplot(3,2,3)
plot(yyyy,X,'k-')
hold on
plot(yyyy,X(:,3),'k-','LineWidth',2)
set(gca,'Xtick',1970:10:2020,'Xgrid','on','Ygrid','on','Xlim',[1972 2025],'Ylim',[-.5 4],'XMinorTick','on','XMinorGrid','on','box','on','Layer','top')
title('(c) NCEP and phyto. respiration','fontweight','normal')
ylabel('CO_2 flux (PgC yr^{-1})')

T = readtable('../Cflu_RSS_n_RVY.xlsx','sheet','Global'); % North, Tropics, South, latitude
X = [T.Cflu_RVY0 T.Cflu_RVY1 T.Cflu_RVA0 T.Cflu_RVY4 T.Cflu_RVY3];
subplot(3,2,4)
plot(yyyy,X,'k-')
hold on
plot(yyyy,X(:,3),'k-','LineWidth',2)
set(gca,'Xtick',1970:10:2020,'Xgrid','on','Ygrid','on','Xlim',[1972 2025],'Ylim',[-.5 4],'XMinorTick','on','XMinorGrid','on','box','on','Layer','top')
title('(d) ERA5 and phyto. respiration','fontweight','normal')
ylabel('CO_2 flux (PgC yr^{-1})')

T = readtable('../Cflu_RSS_n_RNB_RNY.xlsx','sheet','Global'); % North, Tropics, South, latitude
X = [T.Cflu_RNY1 T.Cflu_RNY2 T.Cflu_RNB0 T.Cflu_RNY4 T.Cflu_RNY3 T.Cflu_RNB1  T.Cflu_RNA0  T.Cflu_RNB2 T.Cflu_RNB4];
subplot(3,2,5)
plot(yyyy,X,'k-')
hold on
plot(yyyy,X(:,3),'k-','LineWidth',2)
set(gca,'Xtick',1970:10:2020,'Xgrid','on','Ygrid','on','Xlim',[1972 2025],'Ylim',[-.5 4],'XMinorTick','on','XMinorGrid','on','box','on','Layer','top')
title('(e) NCEP and bact. half saturation & phyto. respiration','fontweight','normal')
ylabel('CO_2 flux (PgC yr^{-1})')

T = readtable('../Cflu_RSS_n_RVY_RVB.xlsx','sheet','Global'); % North, Tropics, South, latitude
X = [T.Cflu_RVB1 T.Cflu_RVA0 T.Cflu_RVB2 T.Cflu_RVB3 T.Cflu_RVY0 T.Cflu_RVY1 T.Cflu_RVY4 T.Cflu_RVY3];
subplot(3,2,6)
plot(yyyy,X,'k-')
hold on
plot(yyyy,X(:,2),'k-','LineWidth',2)
set(gca,'Xtick',1970:10:2020,'Xgrid','on','Ygrid','on','Xlim',[1972 2025],'Ylim',[-.5 4],'XMinorTick','on','XMinorGrid','on','box','on','Layer','top')
title('(f) ERA5 and bact. half saturation & phyto. respiration','fontweight','normal')
ylabel('CO_2 flux (PgC yr^{-1})')

% set(gcf,'PaperPosition',[0 0 30 30])
% print('six_analysis.jpeg','-djpeg','-r300')



