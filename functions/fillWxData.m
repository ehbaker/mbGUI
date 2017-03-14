function [ready]=fillWxData(glacier)
% this function fixes the missing glacier temperature and precip data by
% using the closest city data and estimating regression functions between
% the city data and the glacier data.
%
dbstop if error
warning off MATLAB:table:RowsAddedExistingVars
warning off stats:statrobustfit:IterationLimit

formatSecondaryWxData(glacier); %correct/reformat data from cities
PrimaryWx=readtable(['data/',glacier,'/Input/Input_',glacier,'_Daily_Weather.csv']); %Weather from the nearest Wx station
SecondaryWxNames=importdata(['data/',glacier,'/Input/Input_',glacier,'_SecondaryWxData.csv']);  %file with names of secondary Wx stations
SecondaryWx1=readtable(cell2mat(['data/',glacier,'/Input/SecondaryWxData/Output_',SecondaryWxNames(1),'data.csv'])); %
SecondaryWx2=readtable(cell2mat(['data/',glacier,'/Input/SecondaryWxData/Output_',SecondaryWxNames(2),'data.csv'])); %
MissingPrecipitation=readtable(['data/',glacier,'/Input/Input_',glacier,'_Missing_Precipitation.csv']); % Manual Precip totals from data logger failure dates
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
MonthlyTemperatureRegressions = table;
figure (); hold on 
title(['Monthly temperature regressions ',glacier, 'Glacier'])
monthofyear = month(AllWx.date,'monthofyear');
for m = 1:12
    %%iteratively reweighted least squares regression
    index = find(monthofyear==m);
    X1 = AllWx.T_secondary1(index);
    X2 = AllWx.T_secondary2(index);
    Y = AllWx.T_primary(index);
    b1 = robustfit(X1,Y);
    b2 = robustfit(X2,Y);
    
    %%store the coefficients in a table
    MonthlyTemperatureRegressions.month(m,1) = m;
    MonthlyTemperatureRegressions.intercept1(m,1) = b1(1);
    MonthlyTemperatureRegressions.slope1(m,1) = b1(2);
    MonthlyTemperatureRegressions.rsquared1(m,1) = corr(Y(~isnan(Y) & ~isnan(X1)),b1(1) + b1(2)*X1(~isnan(Y) & ~isnan(X1)))^2;
    MonthlyTemperatureRegressions.intercept2(m,1) = b2(1);
    MonthlyTemperatureRegressions.slope2(m,1) = b2(2);
    MonthlyTemperatureRegressions.rsquared2(m,1) = corr(Y(~isnan(Y) & ~isnan(X2)),b2(1) + b2(2)*X2(~isnan(Y) & ~isnan(X2)))^2;
    
    %%plot the regression
    x = -40:10:30;
    yp1 = b1(1) + b1(2).*x;
    yp2 = b2(1) + b2(2).*x;
    subplot(3,4,m); hold on
    scatter(X1,Y, 'r+')
    plot(x,yp1,'r')
    scatter(X2,Y, 'b+')
    plot(x,yp2,'b')
    text(-28,27, month(AllWx.date(index(1)),'name'))
    text(-28,17,['r^2 = ' num2str(MonthlyTemperatureRegressions.rsquared1(m),2)], 'color', 'r')
    text(-28,7,['r^2 = ' num2str(MonthlyTemperatureRegressions.rsquared2(m),2)], 'color', 'b')
    axis([-30 30 -30 30])
    
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
meanDailyT = nan(366,1);
for d = 1:366
    meanDailyT(d) = nanmean(AllWx.T_primary(dayofyear==d));
    filledT(dayofyear==d & isnan(filledT)) = meanDailyT(d);
end
Tnancount(4) = sum(isnan(filledT));
 
%%Precipitation 
filledP = AllWx.P_primary;
Pnancount(1) = sum(isnan(filledP));

%%Calculate monthly precipitation regressions and fill
MonthlyPrecipitationRegressions = table;
figure (); hold on 
title(['Monthly precipitation regressions ',glacier, 'Glacier'])
for m = 1:12
    %%iteratively reweighted least squares regression
    index = find(monthofyear==m & AllWx.P_primary>=0.001);
    X1 = AllWx.P_secondary1(index);
    X2 = AllWx.P_secondary2(index);
    Y = AllWx.P_primary(index);
    b1 = robustfit(X1,Y,'bisquare',4.685,'off'); %removes the intercept term so that 0 precip secondary = 0 precip primary
    b2 = robustfit(X2,Y,'bisquare',4.685,'off');
    if b1==0 b1 = NaN; end
    if b2==0 b2 = NaN; end
    
   %%store the coefficients in a table
    MonthlyPrecipitationRegressions.month(m,1) = m;
    MonthlyPrecipitationRegressions.slope1(m,1) = b1;
    MonthlyPrecipitationRegressions.rsquared1(m,1) = corr(Y(~isnan(Y) & ~isnan(X1)),b1*X1(~isnan(Y) & ~isnan(X1)))^2;
    MonthlyPrecipitationRegressions.slope2(m,1) = b2;
    MonthlyPrecipitationRegressions.rsquared2(m,1) = corr(Y(~isnan(Y) & ~isnan(X2)),b2*X2(~isnan(Y) & ~isnan(X2)))^2;
    
    %%plot the regression
    x = [0,200];
    yp1 = b1.*x;
    yp2 = b2.*x;
    subplot(3,4,m); hold on
    scatter(X1,Y, 'r+') 
    plot(x,yp1,'r') 
    scatter(X2,Y, 'b+') 
    plot(x,yp2,'b') 
    text(10,220, month(AllWx.date(index(1)),'name'))
    text(10,180,['r^2 = ' num2str(MonthlyPrecipitationRegressions.rsquared1(m),2)], 'color', 'r')
    text(10,140,['r^2 = ' num2str(MonthlyPrecipitationRegressions.rsquared2(m),2)], 'color', 'b')
    axis([0 250 0 250])
    
    %%fill gaps
    index1 = find(monthofyear==m & isnan(filledP));
    filledP(index1) = AllWx.P_secondary1(index1).*b1;
    index2 = find(monthofyear==m & isnan(filledP));
    filledP(index2) = AllWx.P_secondary2(index2).*b2;
end

Pnancount(2) = sum(isnan(filledP));
writetable(MonthlyPrecipitationRegressions,outputPrecipitationRegressions)

%%Fill the remaining NaNs with the average daily precip
dayofyear = day(AllWx.date,'dayofyear');
meanDailyP = nan(366,1);
for d = 1:366
    meanDailyP(d) = nanmean(AllWx.P_primary(dayofyear==d));
    filledP(dayofyear==d & isnan(filledP)) = meanDailyP(d);
end
Pnancount(3) = sum(isnan(filledP));

%use available data to correct some of the missing values
MissingPrecipitation.start_date = datetime(MissingPrecipitation.start_date);
MissingPrecipitation.end_date = datetime(MissingPrecipitation.end_date);
for n = 1:height(MissingPrecipitation) 
    start = find(AllWx.date==MissingPrecipitation.start_date(n));
    finish = find(AllWx.date==MissingPrecipitation.end_date(n));
    estP = sum(filledP(start:finish));
    adj = MissingPrecipitation.P(n)/estP;
    filledP(start:finish) = filledP(start:finish).*adj;
end

filledWx = table;
filledWx.date = AllWx.date;
filledWx.T = filledT;
filledWx.P = filledP;
writetable(filledWx,outputWx)

if Tnancount(4)==0 && Pnancount(3)==0
ready=1;
end
