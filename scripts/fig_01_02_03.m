%% Fig01 - Global and regional
close all
clearvars

rwanted  = {'North','Tropics','South','Global'};
txtTitle = {'(b) North','(c) Tropics','(d) South','(a) Global'};
subr     = [2 3 4 1];
ylimr     = [-.5 1.4; -.5 1.4; -.5 2.5; -1.2 4];
ytick    = [.5 .5 .5 .5];

subp(4).num = [13:19 25:31 37:43 49:55];
subp(1).num = [9:12 21:24];
subp(2).num = [33:36 45:48];
subp(3).num = [57:60 69:72];

for r = 1:4
    subplot(6,12,subp(r).num)
    T        = readtable('data/Cflu_RSS_n_RNB.xlsx','sheet',cell2mat(rwanted(r)));

    Cflu = [T.Cflu_RNB1 T.Cflu_RNA0 T.Cflu_RNB0 T.Cflu_RNB2 T.Cflu_RNB4];
    ylim  = ylimr(r,:);
    
    if r == 4
        Cflu_opt = [T.Cflu_Opt_3rd T.Err_low_3rd T.Err_high_3rd];
        yplim = [ylim(1) ylim(1) + (ylim(2)-ylim(1))*.025];
        msize = 30;
    else
        Cflu_opt = [T.Cflu_Opt T.Err_low T.Err_high];
        yplim = [ylim(1) ylim(1) + (ylim(2)-ylim(1))*.05];
        msize = 15;
    end
    N        = T.n;
    yyyy     = T.years;
    group    = T.group;

    h1 = plot(yyyy,Cflu(:,3),'k-','linewidth',1.5);
    hold on
    h2 = plot(yyyy,Cflu,'k-');
    errorbar(yyyy, Cflu_opt(:,1),Cflu_opt(:,2),Cflu_opt(:,3),'vertical','k.','CapSize',0,'color',[.9 .1 .1],'linewidth',1);
    h3 = scatter(yyyy(group == 0), Cflu_opt(group == 0,1),msize,'filled','MarkerEdgeColor',[.9 .1 .1],'MarkerFaceColor',[1 .7 .7],'linewidth',1);
    scatter(yyyy(group == 3), Cflu_opt(group == 3,1),msize,'filled','MarkerEdgeColor',[.9 .1 .1],'MarkerFaceColor','w','linewidth',1)
    ylabel('CO_2 flux (PgC yr^{-1})')
    title(cell2mat(txtTitle(r)),'fontweight','normal')
    set(gca,'Xtick',1970:10:2020,'Xgrid','on','Ygrid','on','Xlim',[1972 2025],'Ylim',ylim,'Ytick',ylim(1):ytick(r):ylim(2),'XMinorTick','on','XMinorGrid','on','box','on','Layer','top')

    for cpt = 1:length(yyyy)
        xp = [yyyy(cpt)-.5; yyyy(cpt)-.5; yyyy(cpt)+.5; yyyy(cpt)+.5];
        yp = [yplim(1); yplim(2); yplim(2); yplim(1)];
        if group(cpt) == 1
            patch(xp,yp,[.8 .8 .8],'EdgeColor','none')
        elseif group(cpt) == 2
            patch(xp,yp,[1 .7 .7],'EdgeColor','none')
        end
    end
        
    
    if r == 1 || r == 2
        set(gca,'xticklabel','')
    elseif r == 4
        legend([h1,h2(1),h3],{'NEMO-PlankTOM12.1','Perturbed simulations', 'Hybrid approach'},'location','SouthEast')
    end

end

% set(gcf,'PaperPosition',[1 1 30 16.5])
% print('fig01_RNB.jpeg','-djpeg','-r300')

%% Fig02 - Latitudinal

clearvars
% close all

T        = readtable('data/Cflu_RSS_n_RNB.xlsx','sheet','latitude');
Cflu     = [T.Cflu_RNB1 T.Cflu_RNA0 T.Cflu_RNB0 T.Cflu_RNB2 T.Cflu_RNB4];
Cflu_opt = [T.Cflu_Opt T.Err_low T.Err_high];
N        = T.n;
yyyy     = T{:,1};
group    = T.group;
ylim     = [-0.06 0.09];
yplim    = [ylim(1) ylim(1) + (ylim(2)-ylim(1))*.025];

for cpt = 1:length(yyyy)
    xp = [yyyy(cpt)-2.5; yyyy(cpt)-2.5; yyyy(cpt)+2.5; yyyy(cpt)+2.5];
    yp = [yplim(1); yplim(2); yplim(2); yplim(1)];
    if group(cpt) == 1
        patch(xp,yp,[.8 .8 .8],'EdgeColor','none')
    elseif group(cpt) == 2
        patch(xp,yp,[1 .7 .7],'EdgeColor','none')
    end
end
    
hold on
plot([-80 90],[0 0],'k-')
h1 = plot(yyyy,Cflu(:,3),'k-','linewidth',1.5);
h2 = plot(yyyy,Cflu,'k-');
errorbar(yyyy, Cflu_opt(:,1),Cflu_opt(:,2),'vertical','k.','CapSize',0,'color',[.9 .1 .1],'linewidth',1.5)
h3 = scatter(yyyy(group == 0), Cflu_opt(group == 0,1),60,'filled','MarkerEdgeColor',[.9 .1 .1],'MarkerFaceColor',[1 .7 .7],'linewidth',1.5);
scatter(yyyy(group == 3), Cflu_opt(group == 3,1),60,'filled','MarkerEdgeColor',[.9 .1 .1],'MarkerFaceColor','w','linewidth',1.5)
set(gca,'Xtick',-80:10:90,'Xgrid','on','Ygrid','on','Xlim',[-80 90],'Ylim',ylim,'XMinorTick','on','box','on','Layer','top')
ax = gca;
ax.XAxis.MinorTickValues = -80:5:90;
ylabel('CO_2 flux (PgC yr^{-1} degree^{-1})')
xlabel('Latitude (°N)')
legend([h1,h2(1),h3],{'NEMO-PlankTOM12.1','Perturbed simulations', 'Hybrid approach'},'location','NorthEast')

% set(gcf,'PaperPosition',[0 0 20 11.75])
% print('fig02_RNB.jpeg','-djpeg','-r300')

%% Fig03 - 1st Optimized, 2nd Optimized, pCO2-products

clearvars
close all

rwanted  = {'North','Tropics','South','Global'};
txtTitle = {'(b) North','(c) Tropics','(d) South','(a) Global'};
subr     = [2 3 4 1];
ylim     = [-.8 .9; -.8 .9; -.8 .9; -1.6 1.4];
ytick    = [.4 .4 .4 .4];

yavg = [1990:1999];

subp(4).num = [13:19 25:31 37:43 49:55];
subp(1).num = [9:12 21:24];
subp(2).num = [33:36 45:48];
subp(3).num = [57:60 69:72];

for r = 1:4
    subplot(6,12,subp(r).num)
    T        = readtable('data/Cflu_RSS_n_RNB.xlsx','sheet',cell2mat(rwanted(r)));
    PIHM     = T.Cflu_RNB0;
    DATA     = T.Cflu_products;
    DATAstd  = T.Err_products;
    group    = T.group;
    yyyy     = T.years;
    if r == 4
        OPT2 = T.Cflu_Opt_3rd;
    else
        OPT2 = T.Cflu_Opt;
    end
    OPT2_grp0 =  OPT2;
    OPT2_grp0(group ~= 0) =  NaN;

    LongTerm      = mean(PIHM(ismember(yyyy,yavg)));
    LongTerm_opt  = nanmean(OPT2(ismember(yyyy,yavg) & T.group == 0));
    LongTerm_data = mean(DATA(ismember(yyyy,yavg)));
    
    y = [DATA-LongTerm_data+DATAstd;flipud(DATA-LongTerm_data-DATAstd)];
    x = [yyyy; flipud(yyyy)];
    x(isnan(y)) = [];
    y(isnan(y)) = [];
    
    patch(x,y,[.6 .8 .9],'EdgeColor','none','FaceAlpha',.4)
    hold on
    plot([1969 2025],[0 0],'k-')
    h1 = plot(yyyy,DATA-LongTerm_data,'-','color',[.5 .6 .8],'linewidth',1.5);
    h2 = plot(yyyy,PIHM-LongTerm,'k-','linewidth',1.5);

    h3 = plot(yyyy,OPT2_grp0-LongTerm_opt,'o-','color',[.9 .1 .1],'MarkerEdgeColor',[.9 .1 .1],'MarkerFaceColor',[1 .7 .7]);
    h4 = plot(yyyy(group ~= 0), OPT2(group ~= 0)-LongTerm_opt,'o','MarkerEdgeColor',[.9 .1 .1],'MarkerFaceColor','w');
    
    set(gca,'Xtick',1970:10:2020,'Xgrid','on','Ygrid','on','Xlim',[1972 2025],'Ylim',ylim(r,:),'Ytick',ylim(r,1):ytick(r):ylim(r,2),'XMinorTick','on','XMinorGrid','on','box','on','Layer','top')
    
    if r <= 3
        h3.MarkerSize = 4;
        h4.MarkerSize = 4;
    end

    if r == 1 || r == 2
        set(gca,'xticklabel','')
    elseif r == 4
        legend([h2,h3,h1],{'NEMO-PlankTOM12.1','Hybrid approach', 'fCO_2-products'},'location','SouthEast')
    end
    
    ylabel('CO_2 flux (PgC yr^{-1})')
    title(cell2mat(txtTitle(r)),'fontweight','normal')
end

% set(gcf,'PaperPosition',[1 1 30 16.5])
% print('fig03_RNB.jpeg','-djpeg','-r300')
