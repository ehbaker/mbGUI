function [d,y]=hytocal(dh,year)
% this function will convert day of year for a day of hydrolic year to 
% calendar year: goes in as (day, year) comes out same format
%
y=NaN*ones(length(year),1);
d=NaN*ones(length(year),1);
for yi=1:length(year)
    if dh(yi)<=92  % Oct1-Dec31
        y(yi)=year(yi)-1;
        if mod(y(yi),4)~=0  % this is not a leap year, all leapyrs divisible by 4
            d(yi)=dh(yi)+273;
        else % it is a leap year
            d(yi)=dh(yi)+274;
        end
    else % after jan 1
        y(yi)=year(yi);
        d(yi)=dh(yi)-92;
    end  
end  
