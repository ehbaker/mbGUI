function [win_date_numM, win_acc_adjM, win_fxd_date_num, win_fxd_date_adj, net_date_numM, net_abl_adjM, annual_date_numM, annual_abl_adjM,search_min,search_max]=ablationmodel(year,zees,datespringobsM,datefallobsM,datenet,siteM,glacier,lapse,meltrateS,meltrateI,precipratio,plot_ablation_model)
% this function estimates the ablation that occurs between the fall visit
% to the glacier and the end of the summer season and the snow that occurs 
% between the spring visit to the glacier and the end of the winter season
%   datespringobsM - vector of dates of obs in winter/spring
%   datefallobsM - vector of dates of obs in summer/fall
% also need:
%   year - balance year
%   zees - vector of the altitude of the sites
%   datenet - vector of the nets up till fall date (so know if melting ice or snow at start)
%   siteM - vector of 3 A,B,C or D
%   glacier - string stating the glacier name, which matches how it is name
%   in the data directory
% 	lapse - lapse rate
%   meltrateS,meltrateI - melt rates snow,ice in m per degree day
%	precipratio - ratio for each site between gage catch and snow
%   option - 1=compute search_min and search_max, 0=don't
%
% outputs
%   win_date_numM - vector of dates of maximum balance
%   win_acc_adjM - vector of amount of excess snow between datespringobsM
%   and the win_dateM (+)
%   win_fxd_date_num - vector of spring fixed dates: May 1st
%   win_fxd_date_adj - vector of adjustments from spring observation to our
%   spring fixed date on May 1
%   net_dateM - vector of dates from stratigraphic minimum to stratigraphic
%   minimum
%   net_abl_adjM - adjustment to the stratigraphic mass minimum(-)
%   annual_dateM - vector of the end of hydro year date
%   annual_abl_adjM - Adjustment to the end of the hydro year 
%   search_min - matrix of excess ablation from observation day (=0) at each site where
%   includes days from earliest site minimum and latest site minimum
%   rows are for each site + matlab datenum is the last row
%   search_max - matrix of excess snow from observation day (=0) at each site where
%   includes days from earliest site maximum and latest site maximum
%   rows are for each site + matlab datenum is the last row

%% General function outline
% lines ~100-500 model the daily mass balance for the time-period around
% the fall observation to find stratigraphic, and fixed date mass min. Also
% the daily mass balance is later used in glacierwide function to
% calculate the 

%lines ~550 - 700 model daily mass baance for the time-period around spring
%observations to find stratigraphic mass max
% 
dbstop if error
		% the precip and temp data are both in matrix with years as columns 
		% and day of hydrologic year as rows (starting with day 275). Make sure to correct them.
        
		precip=load(['../data/',glacier,'/Input/Corrected_',glacier,'_Precipitation.mat']); %location of precip data for glacier
        precip=precip.precip;
		temp=load(['../data/',glacier,'/Input/Corrected_',glacier,'_Temperature.mat']);    %location of temperature data for glacier
        temp=temp.temp;
        wx_info=importdata(['../data/',glacier,'/Input/Input_',glacier,'_Weather_Station_Information.csv']);% Weather station information. See template file 'GlacierWe_info.txt' in data subdirectory
        weatherAlt=wx_info(1);                                              % Altitdude of glacier weather station
        yrBeg = wx_info(2);                                                 %2 yrs before have a balance
        pre=wx_info(3);                                                     % this is the 
        pos=wx_info(4);

plotver =plot_ablation_model; %make =1 if want to plot yearly&site-ly ablation models, else make =0;
plotverw=plot_ablation_model; %make =1 if want to plot yearly&site-ly ablation models, else make =0;
% NOTE: IF YOU TRY TO PROCESS WITH plotver=1, MORE THAN 20 OR SO YEARS AT A
% TIME, MATLAB WILL RUN OUT OF MEMORY AND CRASH
minT=0; % minimum temp below which dry snow falls rather than rain
maxT=1.7; % maximum temp above which rain falls, in between wet snow falls

[mt,nt]=size(temp);
SiteElevdiff = zees - weatherAlt; 
% if bal==3 %want glac_min & glac_max
    countmax=12; %max number of loops
% else %don't want so cut out the while loop effort
%     countmax=1;
% end  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize matrices  %%%%%%%%%%%%%%%%%
net_date_numM   =NaN*ones(length(siteM(:,1)),1); % vector for stratigraphic fall dates   
net_abl_adjM    =NaN*ones(length(siteM(:,1)),1);  % vector for stratigraphic fall adjustments   
annual_date_numM=NaN*ones(length(siteM(:,1)),1);  % vector for fall fixed dates
annual_abl_adjM =NaN*ones(length(siteM(:,1)),1);  % vector for fall fixed date adjustments 
abl_adj=[];
first_net=inf;
last_net=0;
indfn0=NaN*ones(1,length(siteM(:,1))*2);
indln0=NaN*ones(1,length(siteM(:,1))*2);
indfn=NaN*ones(1,length(siteM(:,1)));
indln=NaN*ones(1,length(siteM(:,1)));
lenn=zeros(length(siteM(:,1))*2,1);
count=0;
numTimes=1;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%      FIND FALL MINIMUM       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lines ~100-550
%loop through wach site and determine the fall mass min date if not known
%find mass minimum
bad_ind=isnan(indfn);
while(sum(bad_ind)>0) && count<countmax %i.e. if a solution has not been found, keep looking for one up till countmax...
%bad_sites exist, meaning their search_min isn't the right length
%shouldn't need more than 2 times through, but put in as error stopper

    for j = find(bad_ind==1);                               %find sites with no solution
        site = siteM(j,:);                                  %get name of site missing measurement
        datefallobs = datefallobsM(j);                      %get date of missing fall observation
        kk=get_siteInd(site,glacier);                       %get the index number of that site

        if isnan(datefallobs)                               %if fall date is a nan 
            datefallobs0 = datenum(year,9,30);              % call fall date September 30th of that year could be any day technically but easier if end of hydro
        
        else                                                % if fall date is not a nan
            datefallobs0 = datefallobs;                     %use the fall date
        end
        
    yd = yearday(datefallobs0);                             % get the year and Julian day of fall observation (either real of september 30)
    d = yd(2);                                              % get the julian day of fall observation
    [hdayind,hyear]=caltohy(yd(2),yd(1));                   % get the hydrologic year and day of the observation (1=oct 1, day 274)
    hyearind = hyear-yrBeg;                                 % hydrologic year index= hydrologic year - first year wx data(2 year before first balance)     
    hyearend = datenum(year,9,30);                          % get the last day of hydrologic year. decimal date for end of hydro year (corresponding with summer/fall visit)
    
    %%%%%%%%%%%%%%Determine if it was a leap year%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if mod(year,4)~=0                                       % if balance yr not a leap year, all leapyrs divisible by 4 leap year 
        lengthofyear = 365;
        oct1 = 274;
    else                                                    %if it was a leap year
        lengthofyear = 366;
        oct1 = 275;
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% DETERMINE IF OBSERVATION WAS BEFORE THE END OF HYDRO YEAR %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if  datefallobs0 < hyearend+1 %obs is before Oct 1 find the end of ablation season
                    % for first time through while loop, make sure at least previous mins 
                    % in there for building search_min 
                    % this won't help for next mins and if datefallobs=NaN till second time
                    % through while loop
        
                    
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%% Pretty sure this section looks for issues
         %%% that are not possible the way the code is writtien
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (j~=1||count>0)&&(first_net<datefallobs-pre)
            pre=datefallobs-first_net;
        end
        if (j~=1||count>0)&&(last_net>hyearend+pos)
            pos=last_net-hyearend;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Set time frame to look at temp and precip before and after
        %%% observations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        prev=(hdayind-pre):lengthofyear;            %days last month of year
        prevYr=hyearind;                            %hydro Year
        post=1:pos;                                 %days in next year
        postYr=hyearind+1;                          %next year
        time = 1;                                   % indicates in later section that this was before end of hydro year
        prevsnow=datenet(j);                        %date of previous snow

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% FIND HOW MUCH SNOW OR MELT OCCURED BEFORE end of hyrdo year %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for ii=pre:-1:1                         % determine how much snow was there to start with
            if prevYr<1                             % you are going out of data range
                T=0;                                %no temp
                P=0;                                %no precip
                prevsnow=0;                         %no previous snow
                
            else                                % Otherwise construct arrays of temperature and precip modeled to the site elevation 
                T=temp(prev(ii),prevYr)+lapse*SiteElevdiff(j)/1000;  %%if you get an error here it is probably a date error in the mbdb file
                P=precip(prev(ii),prevYr)/1000*precipratio(kk);
            end 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% DETERMINE IF PRECIP FELL AS SNOW AT SITE %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if T<maxT                               %Precip fell as snow
                snow=P;                             % The much SWE fell on this day
                prevsnow=prevsnow-snow;             % Subtract that SWE from previous snow fall
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% DETERMINE HOW MUCH MELT HAPPEND %%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if T>minT                               % melt occured
                icmelt=meltrateI*T;                 % potentional ice melt based on temp
                if prevsnow-icmelt<=0           % if ice melt is more than prevsnow all ice melt occured
                    melt=meltrateI*T;
                else                            % melt some (T-x) as ice, and rest (x) as snow
    % then -meltrateI*(T-x)=(prevsnow had on day before) and -meltrateS*x= prevsnow on that day
                    melt=-prevsnow + meltrateS*(T+prevsnow/meltrateI);
                end
                prevsnow=prevsnow-melt;
            end        
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% IF FALL OBSERVATIONS WERE MADE AFTER THE
    %%%%%%%%%%%%%%%%% END OF THE HYDRO
    %%%%%%%%%%%%%%%%% YEAR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    elseif  datefallobs0 >= hyearend+1 % the obs was after the end of the hydro year, need to work backwards as well as forwards
        % for first time through while loop, make sure at least previous mins 
        % in there for building search_min 
        % this won't help for next mins 
        
        if (j~=1||count>0)&&(first_net<hyearend-pre)
            pre=hyearend-first_net;
        end
        if (j~=1||count>0)&&(last_net>datefallobs+pos)
            pos=last_net-datefallobs;
        end
        
        %if after end of hydro year, need to look back to previous year
        prev = (lengthofyear-pre):lengthofyear;             %these are the days of the previous year we will look at
        prevYr=hyearind-1;                                  %this is the previous year we will look at
        post =1:hdayind+pos;                                    %these are the days will will look at after the fall visit
        postYr=hyearind;                                    %this is the year we will look at
        time = 2;                                           % indicates in later section that this was before end of hydro year
        prevsnow=datenet(j);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% FIND HOW MUCH SNOW OR MELT OCCURED AFTER THE END OF THE HYDRO YR %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for ii=hdayind-1:-1:1 % do before obs date, find how much snow start with
            if postYr>nt % don't go out of data year range!
                T=0;
                P=0;
                prevsnow=0;
            else
                T=temp(post(ii),postYr)+lapse*SiteElevdiff(j)/1000;    
                P=precip(post(ii),postYr)/1000*precipratio(kk);
            end 
            if T<maxT % snow
                snow=P;
                prevsnow=prevsnow-snow;
            end
            if T>minT % melt
                icmelt=meltrateI*T;
                if prevsnow-icmelt<=0 % melt all as ice
                    melt=meltrateI*T;
                else % melt some as ice, and rest as snow
                    melt=-prevsnow + meltrateS*(T+prevsnow/meltrateI);
                end
                prevsnow=prevsnow-melt;
            end        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% FIND HOW MUCH SNOW OR MELT OCCURED BEFORE end of hyrdo year %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for ii=length(prev):-1:1 % do before obs date, find how much snow start with
            if prevYr<1 % don't go out of data year range!
                T=0;
                P=0;
                prevsnow=0;
            else
                T=temp(prev(ii),prevYr)+lapse*SiteElevdiff(j)/1000;    
                P=precip(prev(ii),prevYr)/1000*precipratio(kk);
            end 
            if T<maxT % snow
                snow=P;
                prevsnow=prevsnow-snow;
            end
            if T>minT % melt
                icmelt=meltrateI*T;
                if prevsnow-icmelt<=0 % melt all as ice
                    melt=meltrateI*T;
                else % melt some as ice, and rest as snow
                    melt=-prevsnow + meltrateS*(T+prevsnow/meltrateI);
                end
                prevsnow=prevsnow-melt;
            end        
        end
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% IF WE DO NOT KNOW THE DATE OF THE ANNUAL BALANCE
    %%%%    WE NEED TO FIND IT
    %%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if isnan(datefallobs) % we know the annual balance but not a net associated with a date
       time=3; % will be in time=1 case, but with the fake date 9/30/year
    end
    % Initialize matrices
    net= NaN*ones(length(prev),1);
    net2= NaN*ones(length(post),1);
    if prevsnow>0                           % IF we had snow at the site 
        totSnow=prevsnow;                       %this is how much snow
    else                                    % Else we didnt have snow at the site
        totSnow=0;
    end
    
    
    for ii=1:length(prev) % do before hy year end
        melt =0;
        snow =0;
        if prevYr<1 % don't go out of data year range!
            T=0;
            P=0;
            net(ii)=0;
        else
            T=temp(prev(ii),prevYr)+lapse*SiteElevdiff(j)/1000;    
            P=precip(prev(ii),prevYr)/1000*precipratio(kk);
        end 
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
    end
    for ii=1:length(post) % do after hy year end
        melt =0;
        snow =0;
        if postYr>nt % nt= years of data, don't exceed data range
            T=0;
            P=0;
            net2(ii)=0;
        else
            T=temp(post(ii),postYr)+lapse*SiteElevdiff(j)/1000;

            P=precip(post(ii),postYr)/1000*precipratio(kk);
        end 
        if T<maxT % snow
            snow=P; 
            totSnow=totSnow+snow;
        end
        if T>minT 
             snmelt=meltrateS*T;
            if totSnow+snmelt>=0 % melt all as snow
                melt=meltrateS*T; 
                totSnow=totSnow+melt;
            else % melt some as snow, and rest as ice, (if totSnow=0, all as ice)
                melt=-totSnow + meltrateI*(T+totSnow/meltrateS);
                totSnow=0;
            end
        end        
        net2(ii)=melt+snow;
    end   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch time %if fall observation was made before of after the end of the hydro year
        
        case 1              
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Observation was made before the end of the hydro year %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % first find the adjustment to hydroend
            cumnet=cumsum(net); % cummulative net since before obs to hy end
            hynetadj=cumnet(end)-cumnet(pre+1);  % the cummulative net since last obs to the hydrological year end.
            cumnet2=cumnet(end)+cumsum(net2); % cummulative net since hy end to later(+earlier) net)
            % now find the minimum
            [netmin,iinet]=min([cumnet;cumnet2]); % find the mimimum net
            netmin=netmin-cumnet(pre+1); % subtract out part before obs
            if iinet<=length(cumnet) % before Oct 1
                netminday=prev(iinet);
                netminyearind=prevYr;
            else
                netminday=iinet-length(cumnet);
                netminyearind=postYr; % we are into the next hy now
            end
            hynetadj=abs(hynetadj-netmin); %take it as correction from minimum to hynet (so can check always positive)
    % for plotting ablation model for verification     
            mel=[cumnet;cumnet2]-cumnet(pre+1)*ones(length([cumnet;cumnet2]),1);
            melday=[prev,post].';
            [mel1_date(1),mel1_date(2)]=hytocal(melday(1),prevYr+yrBeg); % day,year
            [mel2_date0(1),mel2_date(2)]=hytocal(melday(end),postYr+yrBeg);%year, day
            [net_date(1),net_date(2)]=hytocal(netminday,netminyearind+yrBeg); %year, day
            obsd=d;
            netminb=netmin;
            mel2_date(1)=mel2_date0(1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Observation made after the end of the hdryo year %%%%%%%%%%%%
        %%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        case 2              
            % first find the adjustment to hydroend
            netrev=net2(hdayind:-1:1); %reverse this one
            cumhynet=cumsum(netrev);  % cumulative from hyend to datefallobs
            hynetadj=-cumhynet(end);  %the cummulative melt since last obs to the hy end.
            cumnet=cumsum(net); % cummulative net since before obs to hy end
            cumnet2=cumnet(end)+cumsum(net2); % cummulative net since hy end to later(+earlier) net)
            % now find the minimum
            [netmin,iinet]=min([cumnet;cumnet2]); % find the mimimum net
            netmin=netmin-cumnet2(hdayind); % subtract out part before obs
            if iinet<=length(cumnet) % before Oct 1
                netminday=prev(iinet);
                netminyearind=prevYr;
            else
                netminday=iinet-length(cumnet);
                netminyearind=postYr; % we are into the next hy now
            end
            hynetadj=abs(hynetadj-netmin); %take it as correction from minimum to hynet (so can check always positive)
    % for plotting ablation model for verification              
            mel=[cumnet;cumnet2]-cumnet2(hdayind)*ones(length([cumnet;cumnet2]),1);
                %[cumnet;cumnet2]-(cumnet(pre+1)+netmin)*ones(length([cumnet;cumnet2]),1);
            melday=[prev,post].';
            [mel1_date(1),mel1_date(2)]=hytocal(melday(1),prevYr+yrBeg); %year, day
            [mel2_date0(1),mel2_date(2)]=hytocal(melday(end),postYr+yrBeg);%year, day
            [net_date(1),net_date(2)]=hytocal(netminday,netminyearind+yrBeg);%year, day
            obsd=d;
            netminb=netmin;
            mel2_date(1)=mel2_date0(1);
            if mel2_date(1)<mel1_date(1)%started into next calendar yr b/c obs really late
                mel2_date(1)=mel2_date0(1)+lengthofyear;
                if obsd<=lengthofyear+1-pos %if obs bigger, still in same cal yr but mel2_date(1) is in next
                    obsd=obsd+lengthofyear;
                end
            end
            %stop(here)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% DATE OF Annual Blanace NOT KNOWN
        %%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 3
            cumnet=cumsum(net); % cummulative net since before obs to hy end
            % this netminday is the real day our annual balance was measured on
            % this netmin is the amount have to melt from the real 9/30/year to get to the minimum
            cumnet2=cumnet(end)+cumsum(net2); % cummulative net since hy end to later(+earlier) net)
            % now find the minimum
            [netmin,iinet]=min([cumnet;cumnet2]); % find the mimimum net
            netmin=netmin-cumnet(pre+1); % subtract out part before obs
            if iinet<=length(cumnet) % before Oct 1
                netminday=prev(iinet);
                netminyearind=prevYr;
            else
                netminday=iinet-length(cumnet);
                netminyearind=postYr; % we are into the next hy now
            end
            % the real hynetadj is -netmin (positive) b/c set guess date at hy end and this is 
            % how much were off by
            % hynet happens at a point greater than our minimum where the annual balance is measured
            hynetadj=-netmin; %since measured at minimum day
    % for plotting ablation model for verification             
            mel=[cumnet;cumnet2]-(cumnet(pre+1)+netmin)*ones(length([cumnet;cumnet2]),1);
            melday=[prev,post].';
            [mel1_date(1),mel1_date(2)]=hytocal(melday(1),prevYr+yrBeg); %year, day
            [mel2_date0(1),mel2_date(2)]=hytocal(melday(end),postYr+yrBeg);%year, day
            [net_date(1),net_date(2)]=hytocal(netminday,netminyearind+yrBeg);%year, day
            obsd=net_date(1);
            netminb=0;
            mel2_date(1)=mel2_date0(1);
    end
    
    % for finding glac_min
    mel1=yearday([mel1_date(2),mel1_date(1)]);% goes into this as [year day]
    mel2=yearday([mel2_date(2),mel2_date0(1)]);% goes into this as [year day]
    abl_adj=[abl_adj;(mel1:mel2).',j*ones(length(mel),1),mel]; %#ok<AGROW>
    lenn(numTimes)=length(mel);
%     www(numTimes)=length(mel)
%     mel
    % for plotting ablation model for verification first time only     
    if plotver==1 && count==0
        figure
        plot((mel1_date(1):mel2_date(1)).',mel,'r.-')
        hold on;  
        plot(obsd,0,'ko',net_date(1),netminb,'bo',oct1-1,hynetadj+netminb,'go','markersize',8,'linewidth',2)%+netminb
        hold off;
        set(gca,'fontsize',12);
        legend('Degree Day Model','Field Observation','Stratigraphic Minimum','Fixed Date Minimum','Location','EO')    
        grid on
        name=sprintf('Site %s from Fall obs. to end of %4d balance year',site,postYr-1+yrBeg);
        title(name)
        xlab=sprintf('Calendar day (%3d is 9/30/%4d, %3d is 1/1/%4d)',oct1-1,postYr-1+yrBeg,lengthofyear+1,postYr+yrBeg);
        xlabel(xlab)
        ylabel('Balance zeroed on obs. (m w.e.)')    
    end
    %   

    [net_date(1),net_date(2)]=hytocal(netminday,netminyearind+yrBeg);           % Get calender day and year of mass min 
    net_date_num=yearday([net_date(2) net_date(1)]);                            % Get datenum of mass min
   
    
    
    if time==3                  %%%If annual balance date wasnt known to begin with
        datefallobsM(j)=net_date_num; %so if loops in while, includes the earliest and latest min dates 
    end
    net_abl_adj=netminb; %adjustment param
    annual_date_num=hyearend; %%hrylogic year end date
    annual_abl_adj=hynetadj+netminb;%%hrylogic year adjment
    % build matrix output
    net_date_numM(j)   =net_date_num;   
    net_abl_adjM(j)    =net_abl_adj;    %negative
    annual_date_numM(j)=annual_date_num;
    annual_abl_adjM(j) =annual_abl_adj; %positive
    if net_date_num<first_net
        first_net=net_date_num;
    end
    if net_date_num>last_net
        last_net=net_date_num;
    end
    %
    numTimes=numTimes+1;
    end %out of site loop
    % need to search for glacier min from earliest min day to latest min day so 
    % must check for records not containing earliest min day to latest min day
    n1=1;
   
    for j=1:(numTimes-1) %look for index in site's record
        fnj=find(abl_adj(n1:n1-1+lenn(j),1)==first_net);
        lnj=find(abl_adj(n1:n1-1+lenn(j),1)==last_net);    
        if isempty(fnj)||isempty(lnj) %might as well NaN both, have to redo
            indfn0(j)=NaN;
            indln0(j)=NaN;
        else 
            indfn0(j)=fnj+n1-1;
            indln0(j)=lnj+n1-1;
        end   
        n1=n1+lenn(j);
    end
    if count>0 %not on first run
        doSites=find(bad_ind==1);
        indfn0(doSites)=indfn0(length(siteM(:,1))+doSites);
        indln0(doSites)=indln0(length(siteM(:,1))+doSites);
    end
    indfn=indfn0(1:length(siteM(:,1)));
    indln=indln0(1:length(siteM(:,1)));
    bad_ind=isnan(indfn);
    count=count+1;
end % out of while loop
%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%           FIND SPRING MAXIMUM       %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% FIND THE MAXIMUM OF THE YEAR has a very flat maximum, so not very
% accurate at finding true max unless search a huge number of days
% Initialize matrices 
win_date_numM   =NaN*ones(length(siteM(:,1)),1);  
win_acc_adjM    =NaN*ones(length(siteM(:,1)),1);  
win_fxd_date_adj=NaN*ones(length(siteM(:,1)),1);  
win_fxd_date_num=NaN*ones(length(siteM(:,1)),1); 
acc_adj=[];
first_win=inf;
last_win=0;
indfw0=NaN*ones(1,length(siteM(:,1))*2);
indlw0=NaN*ones(1,length(siteM(:,1))*2);
indfw=NaN*ones(1,length(siteM(:,1)));
indlw=NaN*ones(1,length(siteM(:,1)));
lenw=zeros(length(siteM(:,1))*2,1);
count=0;
numTimes=1;
%
bad_ind=isnan(indfw);
while(sum(bad_ind)>0) && count<countmax   %????? was while(sum(bad_ind>0)) && count<countmax which makes no sense while 5 && 
%bad_sites exist, meaning their search_max isn't the right length
%shouldn't need more than 2 times through, but put in as error stopper
    for j=find(bad_ind==1);%do bad sites only
        site=siteM(j,:);
        datespringobs=datespringobsM(j);
        kk=get_siteInd(site,glacier);      
    %
    ydw=yearday(datespringobs);% gives [yr, julian_day]
    yindw=ydw(1)-yrBeg; % the index (or year # ???) in the precip and temp matrices
    dw=ydw(2); %julian day
    [hdayindw,hyearw]=caltohy(ydw(2),ydw(1)); % this is the hydrologic year and day of the observation (1=oct 1, day 275 cal)
    hyearindw=hyearw-yrBeg;
    %
    if mod(year,4)~=0  % if balance year is not a leap year, all leapyrs divisible by 4 leap year 
        may1=121;
    else
        may1=122;
    end
    [hy_fixed_date_day,hy_fixed_date_yr]=caltohy(may1,ydw(1));
    if hdayindw > 290 %past mid July 
        fprintf(1,'ERROR: %4d, %4d is too late to be a spring/winter visit./n',dw,yindw);
        fprintf(1,'Use 4/1/year for visit date and NaN for winter net value./n');
    end
    mdays=170:260; %start the model in mid March, end it mid June(hydrodays) 
    % Initialize matrices
    win=NaN*ones(length(mdays),1);
    if hdayindw <= 170 %need obs day to be in modelled days, use this if early winter date
        mdays=hdayindw-5:260;
    elseif hdayindw >= 260 
        mdays=170:hdayindw+5;
    end
    % for first time through while loop, make sure at least previous maxs 
    % in there for building search_max 
    % this won't help for next maxs    
    if (j~=1||count>0)&&(first_win<datespringobs-(hdayindw-mdays(1))) % ???
        mdays=(hdayindw-(datespringobs-first_win)):mdays(end); % ???
    end
    if (j~=1||count>0)&&(last_win>datespringobs+(mdays(end)-hdayindw))
        mdays=mdays(1):(hdayindw+(last_win-datespringobs));
    end
    for ii=1:length(mdays) 
        melt =0;
        snow =0;
        
        yr_length = length(temp(:,hyearindw));
        if mdays(ii) <= yr_length; 
            T=temp(mdays(ii),hyearindw)+lapse*SiteElevdiff(j)/1000;
            P=precip(mdays(ii),hyearindw)/1000*precipratio(kk);
        else% if statement added by LS 2013.10 to roll over into next hydro year for late geodetic obs...
            T=temp(mdays(ii)-yr_length,hyearindw+1)+lapse*SiteElevdiff(j)/1000;
            P=precip(mdays(ii)-yr_length,hyearindw+1)/1000*precipratio(kk);
        end
        
        if T<maxT % snow
            snow=P;
        end
        if T>minT % melt
            melt=meltrateS*T; % melt all as snow, ASSUMING not down to ice in winter!
        end        
        win(ii)=melt+snow;
    end
    cumwin=cumsum(win);
    [winmax,iiwin]=max(cumwin); %value and index of the maximum snow thickness in the model + acc - abl
    winmax = winmax-cumwin(hdayindw - mdays(1)) ;%max - snow on day of obs here our adjustment param is max snow
    fixed_date_spring_adj=cumwin(hy_fixed_date_day- mdays(1)+1)-cumwin(hdayindw - mdays(1));
    winmaxday= mdays(iiwin);
    winmaxyearind=hyearindw;
    % for plotting ablation model for verification     
    sno=cumwin-cumwin(hdayindw - mdays(1))*ones(length(cumwin),1);
    snoday=mdays.';
    [sno1_date(1),sno1_date(2)]=hytocal(snoday(1),winmaxyearind+yrBeg); %day,year
    [sno2_date(1),sno2_date(2)]=hytocal(snoday(end),winmaxyearind+yrBeg);%day,year
    [max_date(1),max_date(2)]=hytocal(winmaxday,winmaxyearind+yrBeg); %day,year
    obsd=dw;
    % for finding glac_max
    sno1=yearday([sno1_date(2),sno1_date(1)]);% goes into this as [year day]
    sno2=yearday([sno2_date(2),sno2_date(1)]);% goes into this as [year day]
    acc_adj=[acc_adj;(sno1:sno2).',j*ones(length(sno),1),sno]; %#ok<AGROW>
    lenw(numTimes)=length(sno);
    % for plotting ablation model for verification first time only    
    if plotverw==1 && count==0
        figure
        plot((sno1_date(1):sno2_date(1)).',sno,'b.-')
        hold on;  
        plot(obsd,0,'ko',max_date(1),winmax,'ro',may1,fixed_date_spring_adj,'go','markersize',8,'linewidth',2)
        hold off;
        set(gca,'fontsize',12);
        legend('Degree Day Model','Field Observation','Stratigraphic Maximum','Fixed Date Maximum','Location','EO')    
        grid on
        name=sprintf('Site %s from Spring obs. to max of %4d balance year',site,winmaxyearind+yrBeg);
        title(name)
        xlab=sprintf('Calendar day (%3d is 5/1/%4d)',may1,postYr-1+yrBeg);
        xlabel(xlab)
        ylabel('Balance zeroed on obs. (m w.e.)')    
    end
    %   
    [win_date(1),win_date(2)]=hytocal(winmaxday,winmaxyearind+yrBeg); % should spit out [day year]
    win_date_num=yearday([win_date(2) win_date(1)]); % goes into this as [year day]
    win_acc_adj=winmax; %adjustment param
    win_fxd_date_num(j)=yearday([win_date(2) may1]);
    win_fxd_date_adj(j)=fixed_date_spring_adj;
    % build matrix output
    win_date_numM(j)   =win_date_num;
    win_acc_adjM(j)    =win_acc_adj;    %positive
    if win_date_num<first_win
        first_win=win_date_num;
    end
    if win_date_num>last_win
        last_win=win_date_num;
    end
    %
    numTimes=numTimes+1;
    end %out of site loop
    % need to search for glacier max from earliest max day to latest max day so 
    % must check for records not containing earliest max day to latest max day
    w1=1;
    for j=1:(numTimes-1) %look for index in site's record
        fwj=find(acc_adj(w1:w1-1+lenw(j),1)==first_win);
        lwj=find(acc_adj(w1:w1-1+lenw(j),1)==last_win);
        if isempty(fwj)||isempty(lwj) %might as well NaN both, have to redo
            indfw0(j)=NaN;
            indlw0(j)=NaN;
        else
            indfw0(j)=fwj+w1-1;
            indlw0(j)=lwj+w1-1;
        end  
        w1=w1+lenw(j);
    end
    if count>0 %not on first run
        doSites=find(bad_ind==1);
        indfw0(doSites)=indfw0(length(siteM(:,1))+doSites);
        indlw0(doSites)=indlw0(length(siteM(:,1))+doSites);
    end
    indfw=indfw0(1:length(siteM(:,1)));
    indlw=indlw0(1:length(siteM(:,1)));
    bad_ind=isnan(indfw);
    count=count+1;
end %out of while loop
%
% stop(here)
% if option == 1 %we want stratigraphic balance
%     search_min=NaN;
%     search_max=NaN;
    
% elseif option == 2 %% we want hydro year balance
%     search_min=NaN;
%     search_max=NaN; 
    
    
% build search_min and search_max
% elseif option==3  %want glac_min & glac_max
    minDays=first_net:last_net;
    search_min=NaN*ones(length(siteM(:,1))+1,length(minDays));
    maxDays=first_win:last_win;
    search_max=NaN*ones(length(siteM(:,1))+1,length(maxDays));
    %build with a row for each site, and last row is day
    for j=1:length(siteM(:,1)) 
        search_min(j,:)=abl_adj(indfn(j):indln(j),3).';
        search_max(j,:)=acc_adj(indfw(j):indlw(j),3).';
    end
    search_min(length(siteM(:,1))+1,:)=minDays;
    search_max(length(siteM(:,1))+1,:)=maxDays;
% else %don't want glac_min & glac_max
%     search_min=NaN;
%     search_max=NaN;

% stop(here)
end

    %if year == 1994;
       % stopthisshit(here)
    %end