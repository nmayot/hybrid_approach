%% Global and regional
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

for r = 4
    for p = [0:5]
        
        if p == 0
            T = readtable('../Cflu_RSS_n_RNB.xlsx','sheet',cell2mat(rwanted(r)));
            Cflu = [T.Cflu_RNB1 T.Cflu_RNA0 T.Cflu_RNB0 T.Cflu_RNB2 T.Cflu_RNB4];
            Cflu_opt = [T.Cflu_Opt_3rd T.Err_low_3rd T.Err_high_3rd];
            group    = T.group;
            msize = 30;
            yyyy = T.years;
        elseif p == 1
            T = readtable('../Cflu_RSS_n_RNY.xlsx','sheet',cell2mat(rwanted(r)));
            Cflu = [T.Cflu_RNY1 T.Cflu_RNY2 T.Cflu_RNB0 T.Cflu_RNY4 T.Cflu_RNY3];
            Cflu_opt = [T.Cflu_Opt_3rd T.Err_low_3rd T.Err_high_3rd];
            group    = T.group;
            yyyy = T.years;
            h3 =plot(yyyy, Cflu_opt(:,1),'-','linewidth',1.5);
        
        elseif p == 2
            T = readtable('../Cflu_RSS_n_RNB_RNY.xlsx','sheet',cell2mat(rwanted(r)));
            Cflu = [[T.Cflu_RNY1 T.Cflu_RNY2 T.Cflu_RNB0 T.Cflu_RNY4 T.Cflu_RNY3 T.Cflu_RNB1  T.Cflu_RNA0  T.Cflu_RNB2 T.Cflu_RNB4]];
            Cflu_opt = [T.Cflu_Opt_2nd T.Err_low_2nd T.Err_high_2nd];
            group    = T.group_2nd;
            yyyy = T.years;
            h4 =plot(yyyy, Cflu_opt(:,1),'-','linewidth',1.5);

        elseif p == 3
            T = readtable('../Cflu_RSS_n_RVB.xlsx','sheet',cell2mat(rwanted(r)));
            Cflu = [T.Cflu_RVB0 T.Cflu_RVB1 T.Cflu_RVA0 T.Cflu_RVB2 T.Cflu_RVB3];
            Cflu_opt = [T.Cflu_Opt_3rd T.Err_low_3rd T.Err_high_3rd];
            group    = T.group;
            yyyy = T.years;
            h5=plot(yyyy, Cflu_opt(:,1),'-','linewidth',1.5);

        elseif p == 4
            T = readtable('../Cflu_RSS_n_RVY.xlsx','sheet',cell2mat(rwanted(r)));
            Cflu = [T.Cflu_RVA0 T.Cflu_RVY0 T.Cflu_RVY1 T.Cflu_RVY4 T.Cflu_RVY3];
            Cflu_opt = [T.Cflu_Opt_2nd T.Err_low_2nd T.Err_high_2nd];
            group    = T.group_2nd;
            yyyy = T.years;
            h6=plot(yyyy, Cflu_opt(:,1),'-','linewidth',1.5);

        elseif p == 5
            T = readtable('../Cflu_RSS_n_RVY_RVB.xlsx','sheet',cell2mat(rwanted(r)));
            Cflu = [T.Cflu_RVB1 T.Cflu_RVA0 T.Cflu_RVB2 T.Cflu_RVB3 T.Cflu_RVY0 T.Cflu_RVY1 T.Cflu_RVY4 T.Cflu_RVY3];
            Cflu_opt = [T.Cflu_Opt_2nd T.Err_low_2nd T.Err_high_2nd];
            group    = T.group_2nd;
            yyyy = T.years;
            h7=plot(yyyy, Cflu_opt(:,1),'-','linewidth',1.5);

        end


    
        if p == 0
            yyyy(isnan(Cflu_opt(:,1))) = [];
            Cflu_opt(isnan(Cflu_opt(:,1)),:) = [];
            
            h1 = patch([yyyy; flipud(yyyy)], [Cflu_opt(:,1)-Cflu_opt(:,2); flipud(Cflu_opt(:,1)+Cflu_opt(:,3))],[.75 .75 .75]);
            hold on
            h2 = plot(yyyy, Cflu_opt(:,1),'k-','linewidth',1.5);
            ylabel('CO_2 flux (PgC yr^{-1})')
            title(cell2mat(txtTitle(r)),'fontweight','normal')
            set(gca,'Xtick',1970:10:2020,'Xgrid','on','Ygrid','on','Xlim',[1972 2025],'XMinorTick','on','XMinorGrid','on','box','on','Layer','top')
        else
            plot(yyyy(group == 3), Cflu_opt(group == 3,1),'ko','MarkerFaceColor','w');
            hold on
        end
    end

        
    

legend([h1 h2 h3 h4 h5 h6 h7],{'uncertainty','NCEP - bact.','NCEP - phyto','NCEP - bact. + phyto','ERA - bact','ERA - phyto.','ERA - bact. + phyto'},'location','SouthEast')


end

% set(gcf,'PaperPosition',[1 1 20 10])
% print('fig01_sensitivity.jpeg','-djpeg','-r300')