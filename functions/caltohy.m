function [d,hy]=caltohy(dc,year)
% this function will convert day of year for a calendar year to day of
% hydrolic year: goes in as (day, year) comes out same format
% 
hy=NaN*ones(length(year),1);
d=NaN*ones(length(year),1);
for yi=1:length(year)
    if mod(year(yi),4)~=0  % this is not a leap year, all leapyrs divisible by 4
        if dc(yi)<274
            hy(yi)=year(yi);
            d(yi)=dc(yi)+92;
        else
            hy(yi)=year(yi)+1;
            d(yi)=dc(yi)-273;
        end
    else  % must be a leap year
        if dc(yi)<275
            hy(yi)=year(yi);
            d(yi)=dc(yi)+92;
        else
            hy(yi)=year(yi)+1;
            d(yi)=dc(yi)-274;
        end
    end
end