%% FigS02 - validation

variable = {'AMOC','SO_SSS','SO_SI','Revelle'};
ylabels = {'AMOC at 26°N (Sv)','Southern Ocean SSS (in psu)','Southern Ocean SI (in kg m^{-3})','Surface ocean Revelle factor'};
titles = {'(a) Atlantic Meridional Overturning Circulation (AMOC)','(b) Sea Surface Salinity (SSS)','(c) Stratification index (SI)','(d) Revelle factor'};
for sub = 1:4
    T = readtable('data/GCB2023_metrics.xlsx','Sheet',cell2mat(variable(sub)));
    
    subplot(2,2,sub)
    plot(T{:,1},T{:,2:end-2},'b-')
    hold on
    plot(T{:,1},nanmean(T{:,2:end-2},2),'b-','linewidth',2)
    plot(T{:,1},T{:,end-1},'k-','linewidth',2)
    if sub ~= 4
        plot(T{:,1},T{:,end},'k-')
    end
    ylabel(cell2mat(ylabels(sub)))
    set(gca,'Xlim',[2005 2021])
    title(cell2mat(titles(sub)))
end

sub = 4;
subplot(2,2,sub)
T = readtable('data/GCB2023_metrics.xlsx','Sheet',cell2mat(variable(sub)));
plot(2013,nanmean(T{:,2:end-2}),'bo')
hold on
plot(2013,nanmean(nanmean(T{:,2:end-2})),'bo','linewidth',2)
plot(2013,nanmean(T{:,end-1}),'ko','linewidth',2)
plot(2013,nanmean(T{:,end}),'ko')

% set(gcf,'PaperPosition',[1 1 30 15])
% print('figSXX_validation.jpeg','-djpeg','-r300')