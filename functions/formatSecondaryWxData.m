function [] = formatSecondaryWxData(glacier)
% this function converts units and replaces no data values with nans

SecondaryData=importdata(['data/',glacier,'/Input/Input_',glacier,'_SecondaryWxData.csv']);  %import selected glacier Secondary data
        
for j = 1:2;
    name = SecondaryData(j);
    mySecondary=cell2mat(['data/',glacier,'/Input/SecondaryWxData/Input_',name,'data.csv']);
    outdata =cell2mat(['data/',glacier,'/Input/SecondaryWxData/Output_',name,'data.csv']);
       
    Secondary = readtable(mySecondary);

    Secondary.Maximum_TemperatureInF(Secondary.Maximum_TemperatureInF==9999)=nan; %replace no data value with nan
    Secondary.Minimum_TemperatureInF(Secondary.Minimum_TemperatureInF==9999)=nan;
    Secondary.PrecipitationInInches(Secondary.PrecipitationInInches==9999)=nan;
    date = Secondary.Date;
    Temperature = ((Secondary.Maximum_TemperatureInF - 32)*(5/9) + (Secondary.Minimum_TemperatureInF - 32)*(5/9))/2;
    Precipitation = Secondary.PrecipitationInInches * 25.4;

    newSecondary = table(date,Temperature,Precipitation);
    writetable(newSecondary,outdata);
end
end

