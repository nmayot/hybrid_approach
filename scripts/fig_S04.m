clearvars

% ---- extract values from table
T = readtable(['data/Cflu_RSS_n_RNB.xlsx'],'sheet','Global');
X = [T.Cflu_RNB1 T.Cflu_RNA0 T.Cflu_RNB0 T.Cflu_RNB2 T.Cflu_RNB4];
Y = [T.MSE_RNB1  T.MSE_RNA0  T.MSE_RNB0  T.MSE_RNB2 T.MSE_RNB4];
N = T.n;

yyyy = T.years;
unconstrained = [];
uncertained   = [];
ErrLow  = nan(length(N),1);
ErrHig  = nan(length(N),1);
CflOpt  = nan(length(N),1);
RSSOpt  = nan(length(N),1);
GoodFit = nan(length(N),1);

cpt_w = 33;

for cpt = cpt_w
    
    if N(cpt) ~= 0
        x = X(cpt,:); % Cflu from model runs
        y = Y(cpt,:); % RSS from model runs
        
        p = polyfit(x,y,3); % fit a 3rd polynomial

        yhat = polyval(p,x); % predicted y, to calculate residuals
        RMSE = sqrt(mean((y - yhat).^2));  % Root Mean Squared Error

        rp2 = roots([3*p(1) 2*p(2) p(3)]); % roots 2nd derivative
        discriminant = (2*p(2))^2 - 4*(3*p(1)*p(3)); % D = b2-4ac

        [~, id] = sort([p(1)*rp2.^3 + p(2)*rp2.^2 + p(3)*rp2 + p(4)]);
    
        xmin = rp2(id == 1);
        xmax = rp2(id == 2);
        
        MSE_min = p(1)*xmin.^3 + p(2)*xmin.^2 + p(3)*xmin + p(4);
        MSE_68  = (0.468*N(cpt)/(N(cpt)-2)*sqrt(2*(2*N(cpt)-2)/(N(cpt)*(N(cpt)-4)))+N(cpt)/(N(cpt)-2))*MSE_min;
    
        ConfInt = roots([p(1) p(2) p(3) p(4)-MSE_68])-xmin;
    
        ConfInt_min = abs(max(ConfInt(ConfInt < 0)));
        ConfInt_max = abs(min(ConfInt(ConfInt > 0)));
    
        if xmin > xmax && ~isempty(find(x < xmax,1))
            uncertained = [uncertained, cpt];
        elseif xmin < xmax && ~isempty(find(x > xmax,1))
            uncertained = [uncertained, cpt];
        elseif xmin < min(x) || xmin > max(x)
            uncertained = [uncertained, cpt];
        end

        CflOpt(cpt) = xmin;
        ErrLow(cpt) = ConfInt_min;
        ErrHig(cpt) = ConfInt_max;
        RSSOpt(cpt) = MSE_min;
        GoodFit(cpt) = RMSE;

        if discriminant < 0 % 2nd derivative has no roots, no min turn point
            unconstrained = [unconstrained, cpt];
            CflOpt(cpt) = NaN;
            ErrLow(cpt) = NaN;
            ErrHig(cpt) = NaN;
            RSSOpt(cpt) = NaN;
            GoodFit(cpt)= NaN; 
        end

    end
end

% unconstrained = find(ismember(yyyy,[1973 1974 1978 1980 1981 1982]))'; % for global

% First estimation of constrained values
constrained = find(~ismember(1:length(N),[find(N == 0)',unconstrained, uncertained]));
outlier = diff(prctile(GoodFit(constrained),[25 75]))*1.5+prctile(GoodFit(constrained),75);
uncertained = unique([uncertained, find(GoodFit > outlier)']);

% Remove constrained values with bad RMSE
while ~isempty(find(GoodFit(constrained) > outlier, 1))
    constrained(GoodFit(constrained) > outlier) = [];
    outlier = diff(prctile(GoodFit(constrained),[25 75]))*1.5+prctile(GoodFit(constrained),75);
    uncertained = unique([uncertained, find(GoodFit > outlier)']);
end

ywanted = cpt_w;
for cpt = ywanted
    
    x = X(cpt,:); % Cflu from model runs
    y = Y(cpt,:); % RSS from model runs
    
    xval = [0.5:0.00001:3.25];
    p = polyfit(x,y,3); % Fits a fifth degree polynomial to exp(a) in the Least squares sense
    yval = polyval(p,xval);
    xval2 = [min(x):0.00001:max(x) max(x)];
    yval2 = polyval(p,xval2);
    
    plot(xval,yval,'k--')
    hold on
    h3 = plot(xval2,yval2,'k-','linewidth',1.5);
    plot([xval(1) xval(end)],[RSSOpt(cpt) RSSOpt(cpt)],'k:')
    plot([xval(1) xval(end)],[MSE_68 MSE_68],'k-')
    plot([CflOpt(cpt)-ErrLow(cpt) CflOpt(cpt)-ErrLow(cpt)],[1075 MSE_68],'k:')
    plot([CflOpt(cpt)+ErrHig(cpt) CflOpt(cpt)+ErrHig(cpt)],[1075 MSE_68],'k:')
    plot([CflOpt(cpt) CflOpt(cpt)],[1075 MSE_min],'k:')
    
    h2 = scatter(x([1 2 4 5]),y([1 2 4 5]),100,'MarkerEdgeColor','b','MarkerFaceColor','w');
    h1 = scatter(x(3),y(3),100,'MarkerEdgeColor','b','MarkerFaceColor','w','linewidth',1.5);

    h4 = scatter(CflOpt(cpt),RSSOpt(cpt),50,'ro','filled'); 
    plot([CflOpt(cpt)-ErrLow(cpt) CflOpt(cpt)+ErrHig(cpt)],[RSSOpt(cpt) RSSOpt(cpt)],'r-')
    
    set(gca,'Xlim',[xval(1) xval(end)],'Ylim',[1075 1400],'Xgrid','on','Ygrid','on','box','on','YTick',1000:100:1400,'XTick',.5:.5:3,'XMinorGrid','on','YMinorGrid','on')
    ylabel('MSE (between simulations and SOCAT)')
    xlabel('Annual ocean CO_2 sink (PgC yr^{-1})')
    legend([h1,h2,h3,h4],{'NEMO-PlankTOM12.1','Perturbed simulations','Fitted cubic function','Hybrid approach'})
end

text(1,1140,'MSE_{68%}','HorizontalAlignment', 'right')
text(2.75,1110,'MSE_{min}')

% set(gcf,'PaperPosition',[1 1 20 15]) % 14 11
% print('test_expl.jpg','-djpeg')