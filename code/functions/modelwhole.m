function [net_mod, win_mod]=modelwhole(year,zee,dateprevobs,datespringobs,datefallobs,earlysnow,site,glacier,lapse,meltrateS,meltrateI,precipratio)
% this function estimates the net between these 3 days if the net or winter is missing
%   dateprevobs - date of obs in the previous summer/fall
%   datespringobs - date of obs in winter/spring
%   datefallobs - date of obs in summer/fall
% also need:
%   year - balance year
%   zee - the altitude of the site
%   earlysnow - the snow that came before the prev obs
%   site - A,B,C or D
%   glacier - 0 1 or 2 for the 2 glaciers and the test data
% 	lapse - lapse rate
%   meltrateS,meltrateI - melt rates snow,ice in m per degree
%	precipratio - ratio for each site between gage catch and snow
%
% outputs
%   net_mod - the net amount
%   win_mod - the winter amount from min to max, not exact because need to
%   correct by if came early or late

		precip=load(['../data/',glacier,'/Corrected_',glacier,'Precip.mat']); %location of precip data for glacier
        precip=precip.precip;
		temp=load(['../data/',glacier,'/Corrected_',glacier,'Temp.mat']);    %location of temperature data for glacier
        temp=temp.temp;
        wx_info=importdata(['../data/',glacier,'/Input_',glacier,'Wx_info.csv']);
        weatherAlt=wx_info(1);
        yrBeg = wx_info(2); %2 yrs before have a balance
        pre=wx_info(3); %Wolverine has later minimums than Gulkana
        pos=wx_info(4);

%
minT=0; % minimum temp below which dry snow falls rather than rain
maxT=1.7; % maximum temp above which rain falls, in between wet snow falls
%
SiteElevdiff=zee - weatherAlt;
% 
kk=get_siteInd(site,glacier);
% want to start the day after the prev obs
yd=yearday(dateprevobs+1);% gives [yr, julian_day]
d=yd(2);
[hdayind,hyear]=caltohy(yd(2),yd(1)); % this is the hydrologic year and day of the observation (1=oct 1, day 274 or 275)
hyearind=hyear-yrBeg;
%
if mod(year-1,4)~=0  % previous balance yr a leap year, all leapyrs divisible by 4 leap year 
  	lengthofyear=365;
    oct1=274;
else
    lengthofyear=366;
    oct1=275;
end
%
ydw=yearday(datespringobs);% gives [yr, julian_day]
hdayindw=caltohy(ydw(2),ydw(1)); % this is the hydrologic year and day of the observation (1=oct 1, day 274 or 275)
%
ydf=yearday(datefallobs);% gives [yr, julian_day]
df=ydf(2);
[hdayindf,hyearf]=caltohy(ydf(2),ydf(1)); % this is the hydrologic year and day of the observation (1=oct 1, day 274 or 275)
hyearindf=hyearf-yrBeg;
%
if mod(year,4)~=0  % balance year not a leap year, all leapyrs divisible by 4 leap year 
   	lengthofyearf=365;
    oct1f=274;
else
    lengthofyearf=366;
    oct1f=275;
end       
if d < oct1 && d > 182 %in summer (after July 1) or early fall
    hyearind0=hyearind;
    mdaysE=hdayind:lengthofyear;
    hyearind1=hyearind+1;
    start1=1;
elseif d >= oct1 || d < 90% in later fall or winter(till March 31)
    hyearind0=NaN;
    mdaysE=[];
    hyearind1=hyearind;
    start1=hdayind;
end
if df < oct1f && df > 182 %in summer (after July 1) or early fall
    hyearind1=hyearindf; %is the same number as the previous hyearind1
    finish=hdayindf;
    hyearind2=NaN;
    mdaysL=[];
elseif df >= oct1f || df < 90% in later fall or winter(till March 31)
    hyearind1=hyearindf-1; %is the same number as the previous hyearind1
    finish=lengthofyearf;
    hyearind2=hyearindf;
    mdaysL=1:hdayindf;
end
mdays0=start1:finish; 
mdays=[mdaysE mdays0 mdaysL];
net= NaN*ones(length(mdays),1);
netw= NaN*ones(length(mdaysE)+length(start1:hdayindw),1);
if earlysnow>0
    totSnow=earlysnow;
else
    totSnow=0; 
end
for ii=1:length(mdays) 
    hyearindn=hyearind0; % start with early year index and move to year index then late year index
    if ii>length(mdaysE)
        hyearindn=hyearind1;
    elseif ii>length(mdaysE)+length(mdays0)
        hyearindn=hyearind2;
    end
    melt =0;
    snow =0;
    T=temp(mdays(ii),hyearindn)+lapse*SiteElevdiff/1000;    
    P=precip(mdays(ii),hyearindn)/1000*precipratio(kk); 
    if T<maxT % snow
        snow=P; 
        totSnow=totSnow+snow;
    end
    if T>minT % melt
        snmelt=meltrateS*T;
        if totSnow+snmelt>=0 % melt all as snow
            melt=meltrateS*T;
            totSnow=totSnow+melt;
        else % melt some as snow, and rest as ice, (if totSnow=0, all as ice)
            melt=-totSnow + meltrateI*(T+totSnow/meltrateS);% first term == meltrateS*(-totSnow/meltrateS)
            totSnow=0;
        end
    end         
    net(ii)=melt+snow;
    if ii<=(length(mdaysE)+length(start1:hdayindw))
        netw(ii)=melt+snow;
    end
end
stop(here)
cumnet=cumsum(net); % cummulative net over period
net_mod=cumnet(end);% net is last entry in cumnet vector
cumnetw=cumsum(netw); % cummulative net till winter
win_mod=cumnetw(end); % winter is last entry in cumnetw vector
stop(here)
