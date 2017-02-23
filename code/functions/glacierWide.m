function [glacier_net_date,glacier_win_date,min_to_glacmin,max_to_glacmax]=glacierWide(year,win_date_numM, search_max, net_date_numM, search_min,glacier,zone_weights)
% this function estimates the -ablation that occurs between the site minimums and
% the glacier-wide min, and the -snow that occurs between the site maximums
% and the glacier-wide max 
% All multipliers affect the minimum location (adjusts and zone_weights)
% but addition of constants (internal accumulation and net/win measured
% balances) doesn't affect location
% inputs
%   year
%   win_dateM - vector of dates of maximum site balance
%   search_max - matrix of excess snow from observation day (=0) at each site where
%   includes days from earliest site maximum and latest site maximum
%   rows are for each site + matlab datenum is the last row
%   net_dateM - vector of dates of the last!! siteminimum balance (so the annual balance would be from 
%   the previous minimum to this minimum
%   search_min - matrix of excess ablation from observation day (=0) at each site where
%   includes days from earliest site minimum and latest site minimum
%   rows are for each site + matlab datenum is the last row
%   glacier - 0 1 or 2 for the 2 glaciers and the test data
% outputs
%   glacier_net_date - vector of the oct 1 date, just passing it back out
%   glacier_win_date - vector of the adjustments between the datefallobsM and the annual_dateM. 
%   min_to_glacmin - vector of amount glacier min more than site min (+)
%   max_to_glacmax - vector of amount glacier max less than site max (-)

dbstop if error
    
plotver =0; %make =1 if want to plot yearly&site-ly ablation models, else make =0;
plotverw=0; %make =1 if want to plot yearly&site-ly ablation models, else make =0;
% NOTE: IF YOU TRY TO PROCESS WITH plotver=1, MORE THAN 60 OR SO YEARS AT A
% TIME, MATLAB WILL RUN OUT OF MEMORY AND CRASH
% get adjustment for misrepresenative sites
%adjust=load(myadjust); %format: year adjustNet adjustWin ----> !!!!!
%removed crazy adjustment scheme relics
numSites=length(win_date_numM);


numDaysN=length(search_min(1,:));
summDays=search_min(numSites+1,:);%last row is dates
numDaysW=length(search_max(1,:));
wintDays=search_max(numSites+1,:);%last row is dates

s_minadjd=search_min;
s_maxadjd=search_max;
indSiteMin=NaN*ones(numSites,1);
indSiteMax=NaN*ones(numSites,1);
for j=1:numSites
    indSiteMin(j)=find(summDays==net_date_numM(j)); %where min site date is 
    indSiteMax(j)=find(wintDays==win_date_numM(j)); %where max site date is    
end

% find glacierwides
b_baradj=sum(diag(zone_weights)*s_minadjd(1:numSites,:));
[netGmin,iinet]=min(b_baradj); % find the mimimum glacierwide net
glacier_net_date=summDays(iinet)*ones(numSites,1);%vector of same days
%
b_barWadj=sum(diag(zone_weights)*s_maxadjd(1:numSites,:));
[winGmax,iiwin]=max(b_barWadj); % find the maximum glacierwide win
glacier_win_date=wintDays(iiwin)*ones(numSites,1);%vector of same days
% find corrections
min_to_glacmin=NaN*ones(numSites,1);
max_to_glacmax=NaN*ones(numSites,1);
for j=1:numSites
    indSiteMin(j)=find(summDays==net_date_numM(j)); %where min site date is 
    min_to_glacmin(j)=search_min(j,iinet)-search_min(j,indSiteMin(j));%should be positive
    indSiteMax(j)=find(wintDays==win_date_numM(j)); %where max site date is    
    max_to_glacmax(j)=search_max(j,iiwin)-search_max(j,indSiteMax(j));%should be negative
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotver==1
    s2=zeros(numDaysN,1);
    for i=1:numDaysN
        s20=yearday(summDays(i));%get year, day
        s2(i)=s20(2);
    end
    if mod(s20(1),4)~=0  % if balance yr not a leap year, all leapyrs divisible by 4 leap year 
        oct1=274;
    else
        oct1=275;
    end
    figure
    for j=1:numSites
        plot(s2.',s_minadjd(j,:).','.-','color',[1,1-(j-1)/(numSites-1),0])
        hold on;
    end
    plot(s2.',b_baradj.','k.-')
    hold on;
    leg=char(ones(numSites,15));
    for j=1:numSites
        plot(s2(indSiteMin(j)),s_minadjd(j,indSiteMin(j)),'bo',s2(iinet),s_minadjd(j,iinet),'co','markersize',8,'linewidth',2)
        if j==1
            leg(j,:)=sprintf('lowest site    ');
        elseif j<numSites/2+.5
            leg(j,:)=sprintf('lower mid site ');
        elseif j==numSites/2+.5
            leg(j,:)=sprintf('middle site    ');
        elseif j==numSites
            leg(j,:)=sprintf('highest site   ');
        elseif j>numSites/2+.5
            leg(j,:)=sprintf('higher mid site');
        end
        hold on;
    end
    plot(s2(iinet),b_baradj(iinet),'co','markersize',8,'linewidth',2)
    hold off;
    set(gca,'fontsize',12);
    leg=[leg;'integrated     ';'site minimum   ';'glacier minimum'];
    legend(leg,'Location','EO')    
    grid on
    name=sprintf('Site mins vs. glacier min of %4d balance year',s20(1));
    title(name)
    xlab=sprintf('Calendar day (%3d is 9/30/%4d)',oct1-1,s20(1));
    xlabel(xlab)
    ylabel('Balance zeroed on obs. for each site (m w.e.)')    
end
if plotverw==1
    w2=zeros(numDaysW,1);
    for i=1:numDaysW
        w20=yearday(wintDays(i));%get year, day
        w2(i)=w20(2);
    end
    if mod(w20(1),4)~=0  % if balance yr not a leap year, all leapyrs divisible by 4 leap year 
        may1=121;
    else
        may1=122;
    end
    figure
    for j=1:numSites
        plot(w2.',s_maxadjd(j,:).','.-','color',[0,1-(j-1)/(numSites-1),1])
        hold on;
    end
    plot(w2.',b_barWadj.','k.-')
    hold on;
    leg=char(ones(numSites,15));
    for j=1:numSites
        plot(w2(indSiteMax(j)),s_maxadjd(j,indSiteMax(j)),'ro',w2(iiwin),s_maxadjd(j,iiwin),'mo','markersize',8,'linewidth',2)
        if j==1
            leg(j,:)=sprintf('lowest site    ');
        elseif j<numSites/2+.5
            leg(j,:)=sprintf('lower mid site ');
        elseif j==numSites/2+.5
            leg(j,:)=sprintf('middle site    ');
        elseif j==numSites
            leg(j,:)=sprintf('highest site   ');
        elseif j>numSites/2+.5
            leg(j,:)=sprintf('higher mid site');
        end
        hold on;
    end
    plot(w2(iiwin),b_barWadj(iiwin),'mo','markersize',8,'linewidth',2)
    hold off;
    set(gca,'fontsize',12);
    leg=[leg;'integrated     ';'site maximum   ';'glacier maximum'];
    legend(leg,'Location','EO')    
    grid on
    name=sprintf('Site maxs vs. glacier max of %4d balance year',w20(1));
    title(name)
    xlab=sprintf('Calendar day (%3d is 5/1/%4d)',may1,w20(1));
    xlabel(xlab)
    ylabel('Balance zeroed on obs. for each site (m w.e.)')    
end
%stop(here)
    
    
    