function [degDays]=getXYmelt(year,zee,datespringobs,datefallobs,glacier,lapse)

% this function gives degree days between the spring and fall 
% measurements in a year:
%   datespringobs - date of winter/spring visit
%   datefallobs - this is the following fall's visit
% also need:
%   zee - the altitude of the site
%   year - balance year
%   glacier - 0 1 or 2 for the 2 glaciers and the test data
% 	lapse - lapse rate

% outputs
%   degDays = over spring to fall
		precip=load(['../data/',glacier,'/Corrected_',glacier,'Precip.mat']); %location of precip data for glacier
        precip=precip.precip;
		temp=load(['../data/',glacier,'/Corrected_',glacier,'Temp.mat']);    %location of temperature data for glacier
        temp=temp.temp;
        wx_info=importdata(['../data/',glacier,'/Input_',glacier,'Wx_info.txt']);
        weatherAlt=wx_info(1);
        yrBeg = wx_info(2); %2 yrs before have a balance
        pre=wx_info(3); %Wolverine has later minimums than Gulkana
        pos=wx_info(4);

%
minT=0; % minimum temp below which dry snow falls rather than rain
%

SiteElevdiff=zee - weatherAlt;
%
% want to start the day after the spring obs
ydw=yearday(datespringobs+1);% gives [yr, julian_day]
[hdayindw,hyearw]=caltohy(ydw(2),ydw(1)); % this is the hydrologic year and day of the observation (1=oct 1, day 275)
hyearindw=hyearw-yrBeg;
%
ydf=yearday(datefallobs);% gives [yr, julian_day]
[hdayindf,hyearf]=caltohy(ydf(2),ydf(1)); % this is the hydrologic year and day of the observation (1=oct 1, day 275)
hyearindf=hyearf-yrBeg;
if mod(year,4)~=0  % balance yr is not a leap year, all leapyrs divisible by 4 leap year 
  	lengthofyearf=365;
else
    lengthofyearf=366;
end
%
if hyearindw==hyearindf % first obs and second obs in same hydro year    
    finish=hdayindf;
    hyearind2=NaN;
    mdaysL=[];
elseif hyearindw~=hyearindf % first obs before Oct1 but second isn't
    finish=lengthofyearf;
    hyearind2=hyearindf;
    mdaysL=1:hdayindf;
end
mdays0=hdayindw:finish; 
mdays=[mdays0 mdaysL];
% build constants for inversion
degDays=0;
for ii=1:length(mdays) 
    hyearindn=hyearindw; % start with winter year index and move to late year index
    if ii>length(mdays0)
        hyearindn=hyearind2;
    end
    T=temp(mdays(ii),hyearindn)+lapse*SiteElevdiff/1000;    
    if T>minT % melt
        degDays=degDays+T;
    end  
end
