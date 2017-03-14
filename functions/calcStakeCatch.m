function [ ready,StakeCatch ] = calcStakeCatch( glacier )
%findCatch calculates the relative catch efficiency for each stake
%   
dbstop if error

mb_data = readtable(['data/',glacier,'/Input/Input_',glacier,'_Glaciological_Data.csv']);
wx_data = readtable(['data/',glacier,'/Output/Output_',glacier,'filledWx.csv']);
wx_meta = readtable(['data/',glacier,'/Input/Input_',glacier,'_Wx_MetaData.csv']);
stakecatch=['data/',glacier,'/Output/Calibrated_',glacier,'_CatchEfficiency.csv']; 
lapse=-6.5; %moist adiabatic lapse rate 

mb_data.date_1 = datetime(mb_data.date_1); %let matlab know we are using dates
mb_data.date_2 = datetime(mb_data.date_2);
wx_data.date = datetime(wx_data.date);

stakes = unique(mb_data.site_name);
mb_params = table;
mb_params.bal_year = mb_data.bal_year;
mb_params.site_name = mb_data.site_name;

%% calculate the catch efficiency at each stake
for n = 1:length(stakes)
    stakedata = mb_data(strcmp(mb_data.site_name,stakes(n)),:);
    priorFallDate = datetime(year(stakedata.date_1)-1,10,1);%initialize with october 1 of prior year
    tempdate = [NaT;stakedata.date_2(1:end-1)];
    priorFallDate(~isnat(tempdate)) = tempdate(~isnat(tempdate));%overwrite with prior fall measurement date on that stake
    Z_offset = stakedata.z - wx_meta.Elevation;
    T_offset = Z_offset*lapse/1000;
    for m = 1:height(stakedata)
        start = find(wx_data.date==priorFallDate(m));
        finish = find(wx_data.date==stakedata.date_1(m));
        PrecipGauge = wx_data.P(start:finish)/1000; %precip in the gauge during the measurement period, converted to m
        temp = wx_data.T(start:finish) + T_offset(m); %estimated temperatures at the stake 
        firstaccumulation = find(PrecipGauge>=0.0002 & temp<=2); %first snowfall at a particular stake...
        SnowEstimated(m) = sum(PrecipGauge(firstaccumulation:end)); %so we include the swe added from rain on snow, but not rain at the beginning of the season
        SnowMeasured(m) = stakedata.winter(m) - stakedata.winter_ablation(m);
        clear PrecipGauge temp
    end
    tempCatch = SnowMeasured'./SnowEstimated';
    medCatch = nanmedian(tempCatch);
    tempCatch(isnan(tempCatch) | ~isfinite(tempCatch)) = medCatch; 
    mb_params.Catch(strcmp(mb_data.site_name,stakes(n)),1) = tempCatch;
    
    clear tempCatch medCatch SnowMeasured SnowEstimated 
end

%% calculate melt coefficients

bw = mb_data.winter - mb_data.winter_ablation;
bs = mb_data.net - (mb_data.summer_accumulation + bw);

Z_offset = mb_data.z - wx_meta.Elevation;
T_offset = Z_offset*lapse/1000;
r = 1.00;

for n = 1:height(mb_data) %   NEED to add summer accumulation
    start = find(wx_data.date==mb_data.date_1(n));
    finish = find(wx_data.date==mb_data.date_2(n));
    PDD_t = sum(wx_data.T(start:finish) + T_offset(n));
    if -bs(n)<=bw(n)
        ki(n) = NaN;
        ks(n) = bs(n)/PDD_t;
    else
        ki(n) = (bs(n)*r + bs(n) + bw(n))/PDD_t;
        ks(n) = ki(n)/r;
    end
end

mb_params.k_s = ks';
mb_params.k_i = ki';

figure(); hold on
for n = 1:length(stakes)
    scatter(mb_params.bal_year(strcmp(mb_params.site_name,stakes(n))),mb_params.k_i(strcmp(mb_params.site_name,stakes(n))),'o');
    scatter(mb_params.bal_year(strcmp(mb_params.site_name,stakes(n))),mb_params.k_s(strcmp(mb_params.site_name,stakes(n))),'+');
    
    medki(n) = nanmedian(mb_params.k_i(strcmp(mb_params.site_name,stakes(n))));
    medks(n) = nanmedian(mb_params.k_s(strcmp(mb_params.site_name,stakes(n))));
    stdki(n) = nanstd(mb_params.k_i(strcmp(mb_params.site_name,stakes(n))));
    stdks(n) = nanstd(mb_params.k_s(strcmp(mb_params.site_name,stakes(n))));
    
    
end    

if sum(isnan(StakeCatch))==0
    ready = 1;
    
end



