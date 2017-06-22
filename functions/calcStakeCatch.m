function [ ready, mb_params ] = calcStakeCatch( glacier )
%findCatch calculates the relative catch efficiency for each stake
%   
dbstop if error

%ins
mb_data = readtable(['data/',glacier,'/Input/Input_',glacier,'_Glaciological_Data.csv']);
wx_data = readtable(['data/',glacier,'/Intermediate/',glacier,'FilledWx.csv']);
wx_meta = readtable(['data/',glacier,'/Input/Input_',glacier,'_Wx_MetaData.csv']);

%outs
out = (['data/',glacier,'/Intermediate/',glacier,'_CatchEfficiency.csv']); 


lapse=-6.5; %moist adiabatic lapse rate 

mb_data.Prior_fall_date = datetime(mb_data.Prior_fall_date); %let matlab know we are using dates
mb_data.Spring_date = datetime(mb_data.Spring_date);
mb_data.Fall_date = datetime(mb_data.Fall_date);
wx_data.date = datetime(wx_data.date);

stakes = unique(mb_data.Name);
mb_params = table;
mb_params.Year = mb_data.Year;
mb_params.Name = mb_data.Name;

%% calculate the catch efficiency at each stake
for n = 1:length(stakes)
    stakedata = mb_data(strcmp(mb_data.Name,stakes(n)),:);
    Z_offset = stakedata.Z - wx_meta.Elevation;
    T_offset = Z_offset*lapse/1000;
    for m = 1:height(stakedata)
        start = find(wx_data.date==stakedata.Prior_fall_date(m));
        finish = find(wx_data.date==stakedata.Spring_date(m));
        PrecipGauge = wx_data.P(start:finish)/1000; %precip in the gauge during the measurement period, converted to m
        temp = wx_data.T(start:finish) + T_offset(m); %estimated temperatures at the stake 
        firstaccumulation = find(PrecipGauge>=0.0002 & temp<=2); %first snowfall at a particular stake...
        SnowEstimated(m) = sum(PrecipGauge(firstaccumulation:end)); %so we include the swe added from rain on snow, but not rain at the beginning of the season
        SnowMeasured(m) = stakedata.Winter_accumulation(m);
        clear PrecipGauge temp
    end
    tempCatch = SnowMeasured'./SnowEstimated';
    medCatch = nanmedian(tempCatch);
    tempCatch(~isfinite(tempCatch)) = medCatch; 
    mb_params.Catch(strcmp(mb_data.Name,stakes(n)),1) = tempCatch;
    
    clear tempCatch medCatch SnowMeasured SnowEstimated 
end

%% calculate melt coefficients

Z_offset = mb_data.Z - wx_meta.Elevation;
T_offset = Z_offset*(lapse/1000);
r = 1.25;

for n = 1:height(mb_data) %   
    start = find(wx_data.date==mb_data.Spring_date(n));
    finish = find(wx_data.date==mb_data.Fall_date(n));
    summer_T = wx_data.T(start:finish) + T_offset(n);
    summer_P = wx_data.P(start:finish) * mb_params.Catch(n) / 1000;
    PDD_t(n) = sum(summer_T(summer_T>=0));
    summer_snow(n) = sum(summer_P(summer_T<=0)); %this is the snow that the Wx station data implies fell in mid summer and melted away...
    snow(n) = mb_data.Winter_accumulation(n) + mb_data.Pre_accumulation(n) + summer_snow(n);
    melt(n) = mb_data.Annual_balance(n) - mb_data.Pre_accumulation(n) - mb_data.Winter_accumulation(n) - summer_snow(n);
    if -melt(n)<=snow(n)
        ki(n) = NaN;
        ks(n) = melt(n)/PDD_t(n);
    else
        ki(n) = (snow(n) + melt(n) - snow(n)*r)/PDD_t(n);
        ks(n) = ki(n)/r;
    end
       
end
figure(); hold on
scatter(PDD_t, melt)
fit = fitlm(PDD_t,melt)

mb_params.k_s = ks';
mb_params.k_i = ki';
 figure(); hold on
for n = 1:length(stakes)
    scatter(mb_params.Year(strcmp(mb_params.Name,stakes(n))),mb_params.k_i(strcmp(mb_params.Name,stakes(n))),'o');
    scatter(mb_params.Year(strcmp(mb_params.Name,stakes(n))),mb_params.k_s(strcmp(mb_params.Name,stakes(n))),'+');
    
    medki(n) = nanmedian(mb_params.k_i(strcmp(mb_params.Name,stakes(n))));
    medks(n) = nanmedian(mb_params.k_s(strcmp(mb_params.Name,stakes(n))));
    stdki(n) = nanstd(mb_params.k_i(strcmp(mb_params.Name,stakes(n))));
    stdks(n) = nanstd(mb_params.k_s(strcmp(mb_params.Name,stakes(n))));
    mb_params.k_i(strcmp(mb_params.Name,stakes(n)) & ~isfinite(mb_params.k_i)) = medki(n);
    mb_params.k_s(strcmp(mb_params.Name,stakes(n)) & ~isfinite(mb_params.k_s)) = medks(n);
    
end    

ki_metric = nanmean(stdki);
ks_metric = nanmean(stdks);

figure(); hold on
years = unique(mb_data.Year);
for n = 1:length(years);
    plot(mb_data.Z(mb_data.Year == years(n)),mb_data.Annual_balance(mb_data.Year == years(n))); hold on
end

writetable(mb_params,out);

if sum(isnan(mb_params.Catch))==0
    ready = 1;   
end



