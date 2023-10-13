%% Hybrid approach with a 3rd order polynomial = Global Ocean

clearvars

% ---- extract values from table
T = readtable('data/Cflu_RSS_n_NoCoast.xlsx','sheet','Global');
X = [T.Cflu_PIHQ T.Cflu_PIHM T.Cflu_PIHR T.Cflu_PIHT];
Y = [T.MSE_PIHQ  T.MSE_PIHM  T.MSE_PIHR  T.MSE_PIHT];
N = T.n;

yyyy = T.years;
unconstrained = [];
uncertained   = [];
ErrLow  = nan(length(N),1);
ErrHig  = nan(length(N),1);
CflOpt  = nan(length(N),1);
RSSOpt  = nan(length(N),1);
GoodFit = nan(length(N),1);

for cpt = 1:length(yyyy)
    
    if N(cpt) ~= 0
        x = X(cpt,:); % Cflu from model runs
        y = Y(cpt,:); % RSS from model runs
        
        p = polyfit(x,y,3); % fit a 3rd polynomial

        yhat = polyval(p,x); % predicted y, to calculate residuals
        RMSE = sqrt(mean((y - yhat).^2));  % Root Mean Squared Error

        rp2 = roots([3*p(1) 2*p(2) p(3)]); % roots derivative
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

        if discriminant < 0 % derivative has no roots = no turn point
            unconstrained = [unconstrained, cpt];
            CflOpt(cpt) = NaN;
            ErrLow(cpt) = NaN;
            ErrHig(cpt) = NaN;
            RSSOpt(cpt) = NaN;
            GoodFit(cpt)= NaN; 
        end

    end
end

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

ywanted = uncertained; % represent uncertained years for visual check
for cpt = ywanted
    
    x = X(cpt,:); % Cflu from model runs
    y = Y(cpt,:); % RSS from model runs
    
    xval = [min(x):0.00001:max(x) max(x)];
    p = polyfit(x,y,3);
    yval = polyval(p,xval);

    subplot(ceil(length(ywanted)/4),4,find(ywanted==cpt))
    scatter(x,y) 
    hold on
    plot(xval,yval,'k-')
    scatter(CflOpt(cpt),RSSOpt(cpt),'ro','filled') 
    plot([CflOpt(cpt)-ErrLow(cpt) CflOpt(cpt)+ErrHig(cpt)],[RSSOpt(cpt) RSSOpt(cpt)],'r-')

    if GoodFit(cpt) < outlier
        title([num2str(yyyy(cpt)), ', residual = OK'])
    else
        title(num2str(yyyy(cpt)))
    end
    
end


%% Hybrid approach with a 2nd order polynomial = Regional analysis
clearvars

% ---- extract values from table
T = readtable('data/Cflu_RSS_n_NoCoast.xlsx','sheet','North'); % North, Tropics, South, latitude
X = [T.Cflu_PIHQ T.Cflu_PIHM T.Cflu_PIHR T.Cflu_PIHT];
Y = [T.MSE_PIHQ  T.MSE_PIHM  T.MSE_PIHR  T.MSE_PIHT];
N = T.n;

yyyy = T.years;
unconstrained = [];
uncertained   = [];
ErrLow  = nan(length(N),1);
ErrHig  = nan(length(N),1);
CflOpt  = nan(length(N),1);
RSSOpt  = nan(length(N),1);
GoodFit = nan(length(N),1);

for cpt = 1:length(yyyy)
    
    if N(cpt) ~= 0
        x = X(cpt,:); % Cflu from model runs
        y = Y(cpt,:); % RSS from model runs
        
        [p,~] = polyfit(x,y,2); % fit a 2nd polynomial
        
        yhat = polyval(p,x); % predicted y, to calculate residuals
        RMSE = sqrt(mean((y - yhat).^2));  % Root Mean Squared Error

        xmin  = roots([2*p(1) p(2)]); % roots 2nd derivative
        
        MSE_min = p(1)*xmin.^2 + p(2)*xmin + p(3);
        MSE_68  = (0.468*N(cpt)/(N(cpt)-2)*sqrt(2*(2*N(cpt)-2)/(N(cpt)*(N(cpt)-4)))+N(cpt)/(N(cpt)-2))*MSE_min;
    
        discriminant = (p(2)^2 - 4*p(1)*(p(3)-MSE_68)); % D = b2-4ac
        ConfInt = roots([p(1) p(2) p(3)-MSE_68])-xmin;
        
        ConfInt_min = abs(max(ConfInt(ConfInt < 0)));
        ConfInt_max = abs(min(ConfInt(ConfInt > 0)));
    
        if discriminant < 0
            unconstrained = [unconstrained, cpt];
            xmin = NaN;
            ConfInt_min = NaN;
            ConfInt_max = NaN;
            MSE_min     = NaN;
            RMSE        = NaN;
        elseif xmin < min(x) || xmin > max(x)
            uncertained = [uncertained, cpt];
        end

        CflOpt(cpt)  = xmin;
        ErrLow(cpt)  = ConfInt_min;
        ErrHig(cpt)  = ConfInt_max;
        RSSOpt(cpt)  = MSE_min;
        GoodFit(cpt) = RMSE;

    end
end

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

ywanted = uncertained;
for cpt = ywanted % represent uncertained years for visual check
    
    x = X(cpt,:); % Cflu from model runs
    y = Y(cpt,:); % RSS from model runs
    
    xval = [min(x):0.00001:max(x) max(x)];
    p = polyfit(x,y,2);
    yval = polyval(p,xval);

    subplot(ceil(length(ywanted)/4),4,find(ywanted==cpt))
    scatter(x,y) 
    hold on
    plot(xval,yval,'k-')
    scatter(CflOpt(cpt),RSSOpt(cpt),'ro','filled') 
    plot([CflOpt(cpt)-ErrLow(cpt) CflOpt(cpt)+ErrHig(cpt)],[RSSOpt(cpt) RSSOpt(cpt)],'r-')

    if GoodFit(cpt) < outlier
        title([num2str(yyyy(cpt)), ', residual = OK'])
    else
        title(num2str(yyyy(cpt)))
    end
    
end

%% Hybrid approach with a 3rd order polynomial + 3 years = Global Ocean
clearvars

% ---- extract values from table
T = readtable('data/Cflu_RSS_n_NoCoast.xlsx','sheet','Global');
X = [T.Cflu_PIHQ T.Cflu_PIHM T.Cflu_PIHR T.Cflu_PIHT];
Y = [T.MSE_PIHQ  T.MSE_PIHM  T.MSE_PIHR  T.MSE_PIHT];
N = T.n;
N_3 = sum([N(1:end-2),N(2:end-1),N(3:end)],2);
X_3 = movmean(X,3,'Endpoints','discard');
Y_3 = nan(length(N_3),size(Y,2));
for m = 1:4
    Y_3(:,m) = sum([Y(1:end-2,m),Y(2:end-1,m),Y(3:end,m)].*[N(1:end-2),N(2:end-1),N(3:end)],2,'omitnan')./N_3;
end


yyyy_3 = T.years(2:end-1);
unconstrained = [];
uncertained   = [];
ErrLow  = nan(length(N_3),1);
ErrHig  = nan(length(N_3),1);
CflOpt  = nan(length(N_3),1);
RSSOpt  = nan(length(N_3),1);
GoodFit = nan(length(N_3),1);

clear T X Y N

for cpt = 1:length(yyyy_3)
    
    if N_3(cpt) ~= 0
        x = X_3(cpt,:); % Cflu from model runs
        y = Y_3(cpt,:); % MSE from model runs
        
        p = polyfit(x,y,3); % fit a 3rd polynomial

        yhat = polyval(p,x); % predicted y, to calculate residuals
        RMSE = sqrt(mean((y - yhat).^2));  % Root Mean Squared Error

        rp2 = roots([3*p(1) 2*p(2) p(3)]); % roots 2nd derivative
        discriminant = (2*p(2))^2 - 4*(3*p(1)*p(3)); % D = b2-4ac

        [~, id] = sort([p(1)*rp2.^3 + p(2)*rp2.^2 + p(3)*rp2 + p(4)]);
    
        xmin = rp2(id == 1);
        xmax = rp2(id == 2);
        
        MSE_min = p(1)*xmin.^3 + p(2)*xmin.^2 + p(3)*xmin + p(4);
        MSE_68  = (0.468*N_3(cpt)/(N_3(cpt)-2)*sqrt(2*(2*N_3(cpt)-2)/(N_3(cpt)*(N_3(cpt)-4)))+N_3(cpt)/(N_3(cpt)-2))*MSE_min;
    
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

% First estimation of constrained values
constrained = find(~ismember(1:length(N_3),[find(N_3 == 0)',unconstrained, uncertained]));
outlier = diff(prctile(GoodFit(constrained),[25 75]))*1.5+prctile(GoodFit(constrained),75);
uncertained = unique([uncertained, find(GoodFit > outlier)']);

% Remove constrained values with bad RMSE
while ~isempty(find(GoodFit(constrained) > outlier, 1))
    constrained(GoodFit(constrained) > outlier) = [];
    outlier = diff(prctile(GoodFit(constrained),[25 75]))*1.5+prctile(GoodFit(constrained),75);
    uncertained = unique([uncertained, find(GoodFit > outlier)']);
end

uncertained(end) = [];
CflOpt(uncertained) = NaN;

ywanted = uncertained;
for cpt = ywanted % represent uncertained years for visual check
    
    x = X_3(cpt,:); % Cflu from model runs
    y = Y_3(cpt,:); % RSS from model runs
    
    xval = [min(x):0.00001:max(x) max(x)];
    p = polyfit(x,y,3);
    yval = polyval(p,xval);

    subplot(ceil(length(ywanted)/4),4,find(ywanted==cpt))
    scatter(x,y) 
    hold on
    plot(xval,yval,'k-')
    scatter(CflOpt(cpt),RSSOpt(cpt),'ro','filled') 
    plot([CflOpt(cpt)-ErrLow(cpt) CflOpt(cpt)+ErrHig(cpt)],[RSSOpt(cpt) RSSOpt(cpt)],'r-')

    if GoodFit(cpt) < outlier
        title([num2str(yyyy_3(cpt)), ', residual = OK'])
    else
        title(num2str(yyyy_3(cpt)))
    end
    
end
