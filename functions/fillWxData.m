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
for m = 1:12;
    %%linear regression
    index = find(monthofyear==m);
    tbl1 = table(AllWx.T_primary(index),AllWx.T_secondary1(index));
    lm1 = fitlm(tbl1,'linear'); 
    tbl2 = table(AllWx.T_primary(index),AllWx.T_secondary2(index));
    lm2 = fitlm(tbl2,'linear');
    
    %%store the coefficients
    mth(m) = m;
    int1(m) = lm1.Coefficients.Estimate(1);
    slp1(m) = lm1.Coefficients.Estimate(2);
    rsq1(m) = lm1.Rsquared.Ordinary;
    int2(m) = lm2.Coefficients.Estimate(1);
    slp2(m) = lm2.Coefficients.Estimate(2);
    rsq2(m) = lm2.Rsquared.Ordinary;
    
    %%plot the regression
    subplot(3,4,m)
    scatter(tbl1.Var1,tbl1.Var2, 'r+'); hold on
    scatter(tbl2.Var1, tbl2.Var2, 'b+'); hold on
    title(m)
    
    %%fill gaps
    index1 = find(monthofyear==m & isnan(filledT));
    filledT(index1) = AllWx.T_secondary1(index1).*slp1(m) + int1(m);
    index2 = find(monthofyear==m & isnan(filledT));
    filledT(index2) = AllWx.T_secondary2(index2).*slp2(m) + int2(m);
end
Tnancount(3) = sum(isnan(filledT));
MonthlyTemperatureRegressions = table(mth,slp1,int1,rsq1,slp2,int2,rsq2,{'month','slope1','intercept1','rsquared1','slope2','intercept2','rsquared2'});
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
for m = 1:12;
    %%linear regression
    tbl1 = table(AllWx.P_primary(index),AllWx.P_secondary1(index));
    lm1 = fitlm(tbl1,'linear'); 
    tbl2 = table(AllWx.P_primary(index),AllWx.P_secondary2(index));
    lm2 = fitlm(tbl2,'linear');
    
    %%store the coefficients
    int1(m) = lm1.Coefficients.Estimate(1);
    slp1(m) = lm1.Coefficients.Estimate(2);
    rsq1(m) = lm1.Rsquared.Ordinary;
    int2(m) = lm2.Coefficients.Estimate(1);
    slp2(m) = lm2.Coefficients.Estimate(2);
    rsq2(m) = lm2.Rsquared.Ordinary;
    
    %%plot the regression
    subplot(3,4,m)
    scatter(tbl1.Var1,tbl1.Var2, 'r+'); hold on
    scatter(tbl2.Var1, tbl2.Var2, 'b+'); hold on
    title(m)
    
    %%fill gaps
    index1 = find(monthofyear==m & isnan(filledP));
    filledP(index1) = AllWx.P_secondary1(index1).*slp1(m) + int1(m);
    index2 = find(monthofyear==m & isnan(filledP));
    filledP(index2) = AllWx.P_secondary2(index2).*slp2(m) + int2(m);
end
title(['Monthly precipitation regressions ',glacier, 'Glacier'])
Pnancount(2) = sum(isnan(filledP));
MonthlyPrecipitationRegressions = table(mth,slp1,int1,rsq1,slp2,int2,rsq2,{'month','slope1','intercept1','rsquared1','slope2','intercept2','rsquared2'});
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
