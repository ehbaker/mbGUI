function [precDays]=getXYprec(year,zee,dateprevobs,datespringobs,glacier,lapse)

% this function gives days of precip and between the 2 measurements
% in a year:
%   dateprevobs - date of the prev obs in summer/fall
%   datespringobs - date of winter/spring visit
% also need:
%   year - balance yr
%   zee - the altitude of the site
%   glacier - 0 1 or 2 for the 2 glaciers and the test data
% 	lapse - lapse rate
%
% outputs
%   degDays = days of precip
		precip=load(['../data/',glacier,'/Corrected_',glacier,'Precip.mat']); %location of precip data for glacier
        precip=precip.precip;
		temp=load(['../data/',glacier,'/Corrected_',glacier,'Temp.mat']);    %location of temperature data for glacier
        temp=temp.temp;
        wx_info=importdata(['../data/',glacier,'/Input_',glacier,'Wx_info.txt']);
        weatherAlt=wx_info(1);
        yrBeg = wx_info(2); %2 yrs before have a balance
        pre=wx_info(3); %Wolverine has later minimums than Gulkana
        pos=wx_info(4);


maxT=1.7; % maximum temp above which rain falls, in between wet snow falls
%
SiteElevdiff=zee- weatherAlt;
%
% want to start the day after the prev obs
yd=yearday(dateprevobs+1);% gives [yr, julian_day]
[hdayind,hyear]=caltohy(yd(2),yd(1)); % this is the hydrologic year and day of the observation 
hyearind=hyear-yrBeg;
if mod(year-1,4)~=0  % if balance yr is not a leap year, all leapyrs divisible by 4 leap year 
  	lengthofyear=365;
else
    lengthofyear=366;
end
%
ydw=yearday(datespringobs);% gives [yr, julian_day]
[hdayindw,hyearw]=caltohy(ydw(2),ydw(1)); % this is the hydrologic year and day of the observation 
hyearindw=hyearw-yrBeg;
%      
if hyearind~=hyearindw % first obs before Oct1 but second isn't
    hyearind0=hyearind;
    mdaysE=hdayind:lengthofyear;
    hyearind1=hyearindw;
    start=1;
elseif hyearind==hyearindw % first obs and second obs in same hydro year
    hyearind0=NaN;
    mdaysE=[];
    hyearind1=hyearindw;
    start=hdayind;
end
mdays0=start:hdayindw; 
mdays=[mdaysE mdays0];
% build constants for inversion
precDays=0;
for ii=1:length(mdays) 
    hyearindn=hyearind0; % start with early year index and move to year index then late year index
    if ii>length(mdaysE)
        hyearindn=hyearind1;
    end
    T=temp(mdays(ii),hyearindn)+lapse*SiteElevdiff/1000;    
    P=precip(mdays(ii),hyearindn)/1000; 
    if T<maxT % snow
        precDays=precDays+P;
    end
end
