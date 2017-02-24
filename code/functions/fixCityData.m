function [] = fixCityData(glacier)
% this function converts units and replaces no data values with nans

mycitydata=importdata(['../data/',glacier,'/Input/Input_',glacier,'_CityData.csv']);  %import selected glacier city data
        
for j = 1:2;
    name = mycitydata.textdata(j);
        mycity=cell2mat(['../data/',glacier,'/Input/CityData_wx/Input_',name,'data.csv']);
        outdata =cell2mat(['../data/',glacier,'/Input/CityData_wx/Output_',name,'data.csv']);
       
city = readtable(mycity);

city.Maximum_TemperatureInF(city.Maximum_TemperatureInF==9999)=nan; %replace no data value with nan
city.Minimum_TemperatureInF(city.Minimum_TemperatureInF==9999)=nan;
city.PrecipitationInInches(city.PrecipitationInInches==9999)=nan;
date = city.Date;
Temperature_C = ((city.Maximum_TemperatureInF - 32)*(5/9) + (city.Minimum_TemperatureInF - 32)*(5/9))/2;
Precipitation_mm = city.PrecipitationInInches * 25.4;

newcity = table(date,Temperature_C,Precipitation_mm);
writetable(newcity,outdata);
end
end

