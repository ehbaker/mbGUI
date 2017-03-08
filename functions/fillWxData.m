function [ready]=fillWxData(glacier)
% this function fixes the missing glacier temperature and precip data by
% using the closest city data and estimating regression functions between
% the city data and the glacier data.
%
dbstop if error

formatSecondaryWxData(glacier); %correct/reformat data from cities
PrimaryWx=readtable(['data/',glacier,'/Input/Input_',glacier,'_Daily_Weather.csv']); %Weather from the nearest Wx station
SecondaryWxNames=importdata(['data/',glacier,'/Input/Input_',glacier,'_SecondaryWxData.csv']);  %file with names of secondary Wx stations
SecondaryWx1=readtable(cell2mat(['data/',glacier,'/Input/SecondaryWxData/Output_',SecondaryWxNames(1),'data.csv'])); %
SecondaryWx2=readtable(cell2mat(['data/',glacier,'/Input/SecondaryWxData/Output_',SecondaryWxNames(2),'data.csv'])); %
outputTemperatureRegressions=['data/',glacier,'/Output/Output_',glacier,'TemperatureRegressions.csv']; % New regressions stored here
outputPrecipitationRegressions=['data/',glacier,'/Output/Output_',glacier,'PrecipitationRegressions.csv']; 
outputWx=['data/',glacier,'/Output/Output_',glacier,'filledWx.csv']; % Filled Wx stored here
 
PrimaryWx.date = datetime(PrimaryWx.date);
SecondaryWx1.date = datetime(SecondaryWx1.date);
SecondaryWx2.date = datetime(SecondaryWx2.date);
AllWx_t = outerjoin(PrimaryWx, SecondaryWx1,'Keys','date','MergeKeys',1); %
AllWx = outerjoin(AllWx_t, SecondaryWx2,'Keys','date','MergeKeys',1); %make 1 big table organized by date
AllWx.Properties.VariableNames = {'date' 'T_primary' 'P_primary' 'T_secondary1' 'P_secondary1' 'T_secondary2' 'P_secondary2'};

%%Fill short gaps based on adjacent temperatures
Tnancount(1) = sum(isnan(AllWx.T_primary));
filledT = interp1gap(AllWx.T_primary,3)';% linear interpolation of temperature through small (<3) gaps 
Tnancount(2) = sum(isnan(filledT));

%%Calculate monthly regressions and fill
figure (); hold on 
title(['Monthly temperature regressions ',glacier, 'Glacier'])
monthofyear = month(AllWx.date,'monthofyear');
MonthlyTemperatureRegressions = table;
for m = 1:12;
    %%iteratively reweighted least squares regression
    index = find(monthofyear==m);
    X = AllWx.T_primary(index);
    Y1 = AllWx.T_secondary1(index);
    Y2 = AllWx.T_secondary2(index);
    [b1 stats1] = robustfit(X,Y1);
    [b2 stats2] = robustfit(X,Y2);
    
    %%store the coefficients in a table
    MonthlyTemperatureRegressions.month(m,1) = m;
    MonthlyTemperatureRegressions.intercept1(m,1) = b1(1);
    MonthlyTemperatureRegressions.slope1(m,1) = b1(2);
    MonthlyTemperatureRegressions.rsquared1(m,1) = corr(Y1(~isnan(Y1) & ~isnan(X)),b1(1) + b1(2)*X(~isnan(Y1) & ~isnan(X)))^2;
    MonthlyTemperatureRegressions.intercept2(m,1) = b2(1);
    MonthlyTemperatureRegressions.slope2(m,1) = b2(2);
    MonthlyTemperatureRegressions.rsquared2(m,1) = corr(Y2(~isnan(Y2) & ~isnan(X)),b2(1) + b2(2)*X(~isnan(Y2) & ~isnan(X)))^2;
    
    %%plot the regression
    subplot(3,4,m)
    scatter(X,Y1, 'r+'); hold on
    scatter(X,Y2, 'b+'); hold on
    title(m)
    
    %%fill gaps
    index1 = find(monthofyear==m & isnan(filledT));
    filledT(index1) = b1(1) + AllWx.T_secondary1(index1).*b1(2);
    index2 = find(monthofyear==m & isnan(filledT));
    filledT(index2) = b2(1) + AllWx.T_secondary2(index2).*b2(2);
end
Tnancount(3) = sum(isnan(filledT));
writetable(MonthlyTemperatureRegressions,outputTemperatureRegressions)

%%Fill the remaining NaNs with the average daily temperature
dayofyear = day(AllWx.date,'dayofyear');
for d = 1:366;
    meanDailyT(d) = nanmean(AllWx.T_primary(dayofyear==d));
    index3 = find(dayofyear==d & isnan(filledT));
    filledT(index3) = meanDailyT(d);
end
Tnancount(4) = sum(isnan(filledT));
 
%%Precipitation 
filledP = AllWx.P_primary;
Pnancount(1) = sum(isnan(filledP));

%%Calculate monthly precipitation regressions and fill
figure (); hold on 
MonthlyPrecipitationRegressions = table;
for m = 1:12;
    %%iteratively reweighted least squares regression
    index = find(monthofyear==m);
    X = AllWx.P_primary(index);
    Y1 = AllWx.P_secondary1(index);
    Y2 = AllWx.P_secondary2(index);
    [b1 stats1] = robustfit(X,Y1);
    [b2 stats2] = robustfit(X,Y2);
    
   %%store the coefficients in a table
    MonthlyPrecipitationRegressions.month(m,1) = m;
    MonthlyPrecipitationRegressions.intercept1(m,1) = b1(1);
    MonthlyPrecipitationRegressions.slope1(m,1) = b1(2);
    MonthlyPrecipitationRegressions.rsquared1(m,1) = corr(Y1(~isnan(Y1) & ~isnan(X)),b1(1) + b1(2)*X(~isnan(Y1) & ~isnan(X)))^2;
    MonthlyPrecipitationRegressions.intercept2(m,1) = b2(1);
    MonthlyPrecipitationRegressions.slope2(m,1) = b2(2);
    MonthlyPrecipitationRegressions.rsquared2(m,1) = corr(Y2(~isnan(Y2) & ~isnan(X)),b2(1) + b2(2)*X(~isnan(Y2) & ~isnan(X)))^2;
    
    %%plot the regression
    subplot(3,4,m)
    scatter(X,Y1, 'r+'); hold on
    scatter(X,Y2, 'b+'); hold on
    title(m)
    
    %%fill gaps
    index1 = find(monthofyear==m & isnan(filledP));
    filledP(index1) = b1(1) + AllWx.P_secondary1(index1).*b1(2);
    index2 = find(monthofyear==m & isnan(filledP));
    filledP(index2) = b2(1) + AllWx.P_secondary2(index2).*b2(2);
end
title(['Monthly precipitation regressions ',glacier, 'Glacier'])
Pnancount(2) = sum(isnan(filledP));
writetable(MonthlyPrecipitationRegressions,outputPrecipitationRegressions)

%%Fill the remaining NaNs with the average daily temperature
dayofyear = day(AllWx.date,'dayofyear');
for d = 1:366;
    meanDailyP(d) = nanmean(AllWx.P_primary(dayofyear==d));
    index3 = find(dayofyear==d & isnan(filledP));
    filledP(index3) = meanDailyP(d);
end
Pnancount(3) = sum(isnan(filledP));

filledWx = table;
filledWx.date = AllWx.date;
filledWx.T = filledT;
filledWx.P = filledP;
writetable(filledWx,outputWx)

if Tnancount(4)==0 && Pnancount(3)==0; 
ready=1;
end
