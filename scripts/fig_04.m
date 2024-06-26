clearvars
close all

PROD = [-2.97	17.20	4.86		2.93	5.83	1.90
        -1.36	14.09	9.75		-3.28	0.76	3.68
        2.64	18.35	5.70		6.71	6.90	2.52
        2.01	18.13	8.43		2.21	8.61	1.89
        -3.15	16.53	7.91		-0.18	6.94	2.70
        -4.68	15.63	6.22		0.51	2.10	3.45
        -2.58	17.02	7.40		2.17	4.48	4.65];

GOBM = [5.87	10.13	-4.29		0.31	1.60	2.43
6.97	9.11	-0.82		1.70	3.33	-0.06
11.33	23.63	7.53		1.12	2.56	-0.60
-8.90	19.93	-0.04		1.71	3.90	-0.74
0.93	14.29	6.83		1.82	3.03	2.10
-3.62	26.93	1.20		-1.38	2.78	-0.45
-21.89	12.98	5.25		-0.20	3.15	-0.66
-15.71	7.10	3.35		-2.22	3.10	1.73
-4.43	17.74	7.02		-1.14	3.99	0.64];

SOCAT = [1.91	18.76	4.72];


data = PROD(:,2);
h1 = boxchart(.5*ones(length(data),1), data, 'BoxFaceColor', [.5 .6 .8]);
hold on
scatter(.5*ones(length(data),1), data,'MarkerEdgeColor', 'k','MarkerFaceColor', [.5 .6 .8]);
h2 = scatter(1,SOCAT(2),50,.5,'o','MarkerEdgeColor','k','MarkerFaceColor',[.2 .6 .2]);
data = GOBM(:,2);
h3 = boxchart(1.5*ones(length(data),1), data, 'BoxFaceColor', [.5 .5 .5]);
scatter(1.5*ones(length(data),1), data,'MarkerEdgeColor', 'k','MarkerFaceColor', [.5 .5 .5]);

data = PROD(:,3);
boxchart(3*ones(length(data),1), data, 'BoxFaceColor', [.5 .6 .8])
scatter(3*ones(length(data),1), data,'MarkerEdgeColor', 'k','MarkerFaceColor', [.5 .6 .8]);
scatter(3.5,SOCAT(3),50,'o','MarkerEdgeColor','k','MarkerFaceColor',[.2 .6 .2])
data = GOBM(:,3);
boxchart(4*ones(length(data),1), data, 'BoxFaceColor', [.5 .5 .5])
scatter(4*ones(length(data),1), data,'MarkerEdgeColor', 'k','MarkerFaceColor', [.5 .5 .5]);

plot([5 5],[-25 30],'k-')

data = PROD(:,5);
boxchart(6*ones(length(data),1), data, 'BoxFaceColor', [.5 .6 .8])
scatter(6*ones(length(data),1), data,'MarkerEdgeColor', 'k','MarkerFaceColor', [.5 .6 .8]);
data = GOBM(:,5);
b = boxchart(7*ones(length(data),1), data, 'BoxFaceColor', [.5 .5 .5]);
b.MarkerStyle = '.';
b.MarkerColor = 'w';
scatter(7*ones(length(data),1), data,'MarkerEdgeColor', 'k','MarkerFaceColor', [.5 .5 .5]);

data = PROD(:,6);
boxchart(8.5*ones(length(data),1), data, 'BoxFaceColor', [.5 .6 .8])
scatter(8.5*ones(length(data),1), data,'MarkerEdgeColor', 'k','MarkerFaceColor', [.5 .6 .8]);
data = GOBM(:,6);
boxchart(9.5*ones(length(data),1), data, 'BoxFaceColor', [.5 .5 .5])
scatter(9.5*ones(length(data),1), data,'MarkerEdgeColor', 'k','MarkerFaceColor', [.5 .5 .5]);

legend([h1 h3 h2],{'fCO_2-products','GOBMs','SOCAT'},'location','north')
ylabel('Decadal trend in \DeltafCO_2 (μatm decade^{-1})')
set(gca,'box','on','Xtick',[1 3.5 6.5 9],'XTickLabel',{'2000s','2010s','2000s','2010s'},'Ylim',[-10 30],'Xlim',[0 10],'Ygrid','on','Yminorgrid','on','layer','bottom')

text(1.4,32,'Subsampled','Fontweight','bold')
text(6.55,32,'Not Subsampled','Fontweight','bold')

% set(gcf,'PaperPosition',[1 1 15 10])
% print('fig04_boxplot.jpeg','-djpeg','-r300')