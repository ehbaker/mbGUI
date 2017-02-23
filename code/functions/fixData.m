function [ready]=fixData(glacier)
% this function fixes the missing glacier temperature and precip data by
% using the closest city data and estimating regression functions between
% the city data and the glacier data. The new data are stored in .csv files
% e.g use regression function Temp(at Wolverine met station)=Mf*T(seward) - Ms 
% Mf and Ms have monthly values, old from open file report 91-246 based on 
% 67-88 data correlations between seward and wolverine.
% For Wolverine open file old values are:
%       Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep    Oct    Nov    Dec   
%OLDMf=[0.795  0.740  0.952  1.195  1.308  0.885  1.076  1.106  1.496  1.480  0.994  0.560] 
%OLDMs=[4.68   4.95   5.99   7.34   8.98   5.25   7.16   7.13   10.76  8.29   5.31   5.58 ]
% We also have Prec(at Wolverine met station)=Pf*P(seward)
% Doesn't make sense to add an intercept because will get a negative precip
% possibly then; however, note that Rsquareds can be negative now because
% not doing proper regression
% inputs
%   glacier - 0 1 or 2 for the 2 glaciers and the test data
%
dbstop if error
warning off MATLAB:rankDeficientMatrix
warning off MATLAB:xlswrite:NoCOMServer
warning off MATLAB:xlswrite:AddSheet 
ready=0;
fixCityData(glacier); %correct/reformat data from cities

    	myglacier=glacier;  %Selected Glacier
        mycitydata=importdata(['../data/',glacier,'/Input_',glacier,'CityData.txt']);  %import selected glacier city data
        mycit1=mycitydata.textdata(1); %Longest city data location if second/better location does not cover entire glacier record
        yrBeg=mycitydata.data(1);  %First year of city data
        mycity2=mycitydata.textdata(2);     %Second/best city data location
        yrSplit=mycitydata.data(2);  % if your data is split between two periods, you have to end split on 12/31/yrSplit-1
        % the precip and temp data are both in matrix with years as columns
        % and day of hydrologic year as rows (starting with day 275). 
        % first year is hydrological year 1965, so start at 10/1/64
        dailyprecip=importdata(['../data/',glacier,'/Input_',glacier,'DailyPrecip.txt']); %(366,nyears) non-leap years have one too many days
        dailytemp=importdata(['../data/',glacier,'/Input_',glacier,'DailyTemp.txt']);
        missingprecip = load(['../data/',glacier,'/Input_',glacier,'MissingPrecip.txt']); % intervals of missing data plus the total precip in mm recorded at the gage during that interval
        % Missing data will be corrected by using the Seward, AK, data till 12/31/86,
        % then the Seward 19 N, AK, data on 1/1/87 and after
        myweather=['../data/',glacier,'/Output_',glacier,'TempPrecipRegressions.txt']; % New regressions stored here
        load(cell2mat(['../data/',glacier,'/CityData_wx/Input_',mycity2,'dates.mat']));
        load(cell2mat(['../data/',glacier,'/CityData_wx/Input_',mycity2,'data.mat'])); % these newcit data are Year Month Day PrecipInches MaxDailyT MinDailyT in mm and C
       
         %creates thehy and theyd
    %stop(here)
%
file = fopen(myweather,'w');
cityPT= newcit;
cityHy = thehy;
cityYd = theyd;    
precip=dailyprecip.data;
temp=dailytemp.data;
[mt,nt]=size(temp);
[mp,np]=size(precip);
if mp==mt && mp~=366 % for some reason Macs think the file is 368 long
    precip=dailyprecip.data(2:end,:);
    temp=dailytemp.data(2:end,:);
    [mt,nt]=size(temp);
    [mp,np]=size(precip);
end
if np~=nt
    fprintf(1,'ERROR: you have %d years of precip data and %d years of temp data.\n',np,nt);
    fprintf(1,'    If you are using a PC and your last year has Oct1 as NaN, you \n');
    fprintf(1,'    must make this a number (suggest 0) to get the columns loaded.\n');
end
% the old excel might save dailyprecip odd, so it thinks the header row/col
% is a .data and not a .txt, be careful
if mp~=mt || mp~=366 || mt~=366
    fprintf(1,'ERROR: you have %d days of precip data and %d days of temp data and you should have 366.\n',mp,mt);
end
% 
for i=1:nt % go through the years and deal with leapyears 
    %START WITH 1965 hydro year! I put NaNs in for all dates in 1965 to 1967 because no glacier data
    y= i + yrBeg;
    if mod(y,4)~=0  % this is not a leap year, all leapyrs divisible by 4 leap year 
        temp(152:end-1,i)=temp(153:end,i);  % get rid of the Feb 29 data (was nan)
        precip(152:end-1,i)=precip(153:end,i); %This should leave the last value duplicated in a non leap year, and it shouldnt ever be used
    end
end 
%find regression equation between temp at city and temp at glacier, so must
%build the ytemp and xtemp for ytemp= Mf*xtemp-Ms
% build at once and separate into months
%
% xtempAll structure 
%    month      1           2          . . 12
%year,day 1,1   x1,1     1  x1,3     1     x1,23     1  
%         1,2   x2,1     1  x2,3     1     x2,23     1
%         1,3   x3,1     1  x3,3     1     x3,23     1
%         .
%         .
%         1,31  x31,1    1  x31,3    1     x31,23    1
%         2,1   x32,1    1  x32,3    1     x32,23    1
%         .
%         .
%         nt,31 x31*nt,1 1  x31*nt,3 1     x31*nt,23 1
%
% so julianyear-1963,day,month = cc,d,m goes in entry x(c-1)*31+d,2*m-1
% ytemp structure is same except that goes in entry y(c-1)*31+d,m
% WILL SPLIT INTO 2 SECTIONS BECAUSE HAVE TO REGRESS FOR 2 EQUATIONS FOR
% GULKANA
% if the x or y doesn't exist, will put NaN's there
% the temp matrix is in hydrologic years.
numyrs1 = yrSplit-yrBeg; %MUST REMEMBER TO END ON Dec 31
if numyrs1==0 
    numyrs1b=1; %so has a matrix if don't split
else
    numyrs1b=numyrs1;
end
ytemp1=zeros(numyrs1b*31,12);
xtemp1=zeros(numyrs1b*31,24);
%temp has last day on 9/30/year=hyyearr, but first day on 10/1/year=hyyear-1
% so julian set of years is hydro set of years +1
numyrs2 = nt+1 -numyrs1;
ytemp2=zeros(numyrs2*31,12);
xtemp2=zeros(numyrs2*31,24);
for cc=1:nt % go through the hydro years 
    for rr=1:mt % go through the hydro days
        if ~isnan(temp(rr,cc))
            cityInd=find(cityHy(:,1)==cc+yrBeg & cityHy(:,2)==rr);
            minT=cityPT(cityInd,6);
            maxT=cityPT(cityInd,5);            
            if ~isempty(cityInd)&&~isnan(minT)&&~isnan(maxT)
                %need year,month,day index to build regression matrices
                m = cityPT(cityInd,2);
                d = cityPT(cityInd,3);
                if cityPT(cityInd,1)<yrSplit
                    c=cityPT(cityInd,1)-(yrBeg-1);
                	xtemp1((c-1)*31+d,2*m-1:2*m)= [mean([minT,maxT]),1]; 
                    ytemp1((c-1)*31+d,m)= temp(rr,cc);
                else
                    c=cityPT(cityInd,1)-(yrSplit-1);
                	xtemp2((c-1)*31+d,2*m-1:2*m)= [mean([minT,maxT]),1]; 
                    ytemp2((c-1)*31+d,m)= temp(rr,cc);
                end
            end
        end
    end
end
% do regressions
beta1=NaN*ones(12,2);
beta2=NaN*ones(12,2);
mSSE1=NaN*ones(numyrs1b*31,12);
mSSE2=NaN*ones(numyrs2*31,12);
for mm=1:12
    beta1(mm,:)=xtemp1(:,2*mm-1:2*mm)\ytemp1(:,mm);
    beta2(mm,:)=xtemp2(:,2*mm-1:2*mm)\ytemp2(:,mm);
    mSSE1(:,mm)= (ytemp1(:,mm) - xtemp1(:,2*mm-1:2*mm)*beta1(mm,:).').^2;
    mSSE2(:,mm)= (ytemp2(:,mm) - xtemp2(:,2*mm-1:2*mm)*beta2(mm,:).').^2;
end
Mf=[ beta1(:,1).'; beta2(:,1).'];
Ms=[-beta1(:,2).';-beta2(:,2).'];
% SSE is vector 1X12 because sum of matrix is a row vector containing the sum value of 
% each column.
SSE1=sum(mSSE1);
SSE2=sum(mSSE2);
% SST is vector 1X12 because mean of matrix is a row vector containing the mean value of 
% each column and sum of matrix is a row vector containing the sum value of each column.
SST1= sum((ytemp1 - ones(numyrs1b*31,1)*mean(ytemp1)).^2); 
SST2= sum((ytemp2 - ones(numyrs2*31,1)*mean(ytemp2)).^2); 
%
Rsq=ones(2,12)-[SSE1./SST1;SSE2./SST2];
%
%temperature corrections
% the matrix is in hydrologic years.
kk=0;
kkk=0;
for cc=1:nt % go through the hydro years 
    for rr=1:mt % go through the hydro days            
        if isnan(temp(rr,cc))
            if mod(cc+yrBeg,4)~=0  % this is not a leap year, all leapyrs divisible by 4 leap year 
                ly=365;
            else
                ly=366;
            end
            cityInd=find(cityHy(:,1)==cc+yrBeg & cityHy(:,2)==rr);
            minT=cityPT(cityInd,6);
            maxT=cityPT(cityInd,5);
            %need year index and month index for regression coeffs
            if cityPT(cityInd,1)<yrSplit
                yrRow=1;
            else
                yrRow=2;
            end
            month = cityPT(cityInd,2);
            if ( isempty(cityInd) || isnan(minT) || isnan(maxT) ) && (rr+1<=ly && rr-1>=1)
                Tnew=(temp(rr-1,cc)+temp(rr+1,cc))/2; % not good, but will do for now. %if no data for city, so average of before and after
                kkk=kkk+1;
            elseif (isempty(cityInd)||isnan(minT)||isnan(maxT))&&(rr+1>ly || rr-1<1)    
                Tnew=NaN;
            else
                Tnew=Mf(yrRow,month).*mean([minT,maxT])-Ms(yrRow,month);   
                kk=kk+1;
            end
            if isnan(Tnew)&&~isempty(cityInd-1)&&~isempty(cityInd+1)&&(rr+1<=ly && rr-1>=1) %average again  
                minT=mean([cityPT(cityInd-1,6),cityPT(cityInd+1,6)]);
                maxT=mean([cityPT(cityInd-1,5),cityPT(cityInd+1,5)]);
                Tnew=Mf(yrRow,month).*mean([minT,maxT])-Ms(yrRow,month);
                kk=kk+1;
            end
            if isnan(Tnew) && (rr+2<=ly && rr-2>=1) %and again
                Tnew=(temp(rr-2,cc)+temp(rr+2,cc))/2;   
                kkk=kkk+1;
            end
            if isnan(Tnew)&&~isempty(cityInd-2)&&~isempty(cityInd+2)&&(rr+2<=ly && rr-2>=1) %and again
                minT=mean([cityPT(cityInd-2,6),cityPT(cityInd+2,6)]);
                maxT=mean([cityPT(cityInd-2,5),cityPT(cityInd+2,5)]);
                Tnew=Mf(yrRow,month).*mean([minT,maxT])-Ms(yrRow,month);
                kk=kk+1;
            end
            temp(rr,cc)=Tnew;
        end
    end
end
fprintf(file,'%d of %d temps have been replaced by the average of temp days before and after.\r\n',kkk,mt*nt);
fprintf(file,'%d of %d temps have been replaced by the regression function of the city temp.\r\n',kk,mt*nt);
%
fprintf(1,'%d of %d temps have been replaced by the average of temp days before and after.\n',kkk,mt*nt);
fprintf(1,'%d of %d temps have been replaced by the regression function of the city temp.\n',kk,mt*nt);
k=0;
avgtemp=NaN*ones(mt,1);
for rr=1:mt  % go through the hydro days
%these are the daily averages over all years, if there is a missing data point, use the day average
	avgtemp(rr)=gmean(temp(rr,1:nt)); 
    for cc=1:nt % go through the hydro years 
        if isnan(temp(rr,cc))
            temp(rr,cc)=avgtemp(rr);
            k=k+1;
        end
    end
end
fprintf(file,'%d of %d temps have been replaced by average temp over all yrs of that day.\r\n',k,mt*nt);
%
fprintf(1,'%d of %d temps have been replaced by average temp over all yrs of that day.\n',k,mt*nt);
%extra check
m0=0;
for cc=1:nt % go through the hydro years    
    for rr=1:mt  % go through the hydro days
        if isnan(temp(rr,cc))
            m0=m0+1;
        end
    end
end
if m0==0     
	fprintf(1,'Temp data ready!\n')	
else
	fprintf(1,'You still have %d NaNs.\n',m)
end
% The missing file for the precips does not cover all NaNs!!!
% precip corrections is based on taking the total precip for missing periods and distributing it 
% across the missing time period by assuming the pattern of precip is the same at nearest city. 
% START WITH 1965! I put NaNs in for all dates in 1965 to 1969 because no glacier data 
% Note: 1967 to 1969 are screwed up in the original dailyprecip; CHECK IF THERE IS BETTER!
%
% precip corrections
% the matrix is in hydrologic years
[m1,n1]=size(missingprecip); %#ok<NASGU>
date1  =NaN*ones(m1,2);
date2  =NaN*ones(m1,2);
percent1=NaN*ones(m1,1);
for ff=1:m1
    date1(ff,:)=yearday(missingprecip(ff,1:3)); %Julian Days
    date2(ff,:)=yearday(missingprecip(ff,4:6));
    d1ind=find(cityYd(:,1)==date1(ff,1) & cityYd(:,2)==date1(ff,2),1); %gets indecies for Julian days in missing interval
    d2ind=find(cityYd(:,1)==date2(ff,1) & cityYd(:,2)==date2(ff,2),1);
    cityprecip=cityPT(d1ind:d2ind,4); % grab precip data
    charttotal=missingprecip(ff,7); %total recorded precip on the glacier
% take precentage of total precip that fell each day and scale it for the total precip that fell at glacier.    
    percent1(ff)=charttotal/gsum(cityprecip) ;
    newprecip=cityprecip*percent1(ff); 
% now need to put it back in the daily precip matrix
    [date1hy,date1y]=caltohy(date1(ff,2),date1(ff,1)); %#ok<NASGU> %hydrologic day and year
    [date2hy,date2y]=caltohy(date2(ff,2),date2(ff,1)); 
    date2y=date2y-yrBeg;   % if its in one hydro year, need only one - date2y is the hydro year.
    % column for the year.
    precip(date1hy:date2hy,date2y)=newprecip;
end
% now do like did for temp
%
yprecip1=zeros(numyrs1b*31,12);
xprecip1=zeros(numyrs1b*31,12);
yprecip2=zeros(numyrs2*31,12);
xprecip2=zeros(numyrs2*31,12);
for cc=1:np % go through the hydro years 
    for rr=1:mp % go through the hydro days
        if ~isnan(precip(rr,cc))
            cityInd=find(cityHy(:,1)==cc+yrBeg & cityHy(:,2)==rr);
            P=cityPT(cityInd,4);           
            if ~isempty(cityInd)&&~isnan(P)
                %need year,month,day index to build regression matrices
                m = cityPT(cityInd,2);
                d = cityPT(cityInd,3);
                if cityPT(cityInd,1)<yrSplit
                    c=cityPT(cityInd,1)-(yrBeg-1);
                	xprecip1((c-1)*31+d,m)= P; 
                    yprecip1((c-1)*31+d,m)= precip(rr,cc);
                else
                    c=cityPT(cityInd,1)-(yrSplit-1);
                	xprecip2((c-1)*31+d,m)= P; 
                    yprecip2((c-1)*31+d,m)= precip(rr,cc);
                end
            end
        end
    end
end
% do regressions
beta1=NaN*ones(12,1);
beta2=NaN*ones(12,1);
mSSE1=NaN*ones(numyrs1b*31,12);
mSSE2=NaN*ones(numyrs2*31,12);
for mm=1:12
    beta1(mm)=xprecip1(:,mm)\yprecip1(:,mm);
    beta2(mm)=xprecip2(:,mm)\yprecip2(:,mm);
    mSSE1(:,mm)= (yprecip1(:,mm) - xprecip1(:,mm)*beta1(mm).').^2;
    mSSE2(:,mm)= (yprecip2(:,mm) - xprecip2(:,mm)*beta2(mm).').^2;
end
Pf=[ beta1.'; beta2.'];
%
SSE1=sum(mSSE1);
SSE2=sum(mSSE2);
% 
SST1= sum((yprecip1 - ones(numyrs1b*31,1)*mean(yprecip1)).^2); 
SST2= sum((yprecip2 - ones(numyrs2*31,1)*mean(yprecip2)).^2);
%
RsqP=ones(2,12)-[SSE1./SST1;SSE2./SST2];
%
%precip corrections
kk=0;
for cc=1:np % go through the hydro years   
    for rr=1:mp  % go through the hydro days
        if isnan(precip(rr,cc))
            cityInd=find(cityHy(:,1)==cc+yrBeg & cityHy(:,2)==rr);
            P=cityPT(cityInd,4);            
            %need year index and month index for regression coeffs
            if cityPT(cityInd,1)<yrSplit
                yrRow=1;
            else
                yrRow=2;
            end
            month = cityPT(cityInd,2);
            if isempty(cityInd)||isnan(P) %if no data for city
                Pnew=NaN;
            else
                Pnew=Pf(yrRow,month).*P;   
                kk=kk+1;
            end
            precip(rr,cc)=Pnew;
        end
    end
end
fprintf(1,'%d of %d precips have been replaced by the regression function of the city precip.\n',kk,mp*np);
%
fprintf(file,'%d of %d precips have been replaced by the regression function of the city precip.\r\n',kk,mp*np);
%
avgprecip=NaN*ones(mp,1);
k=0;
for rr=1:mp  % go through the hydro days 
%these are the daily averages over all years, if there is still a missing data point, use the day average 
	avgprecip(rr)=gmean(precip(rr,1:np));
    for cc=1:np %go through the hydro years 
        if isnan(precip(rr,cc))
            precip(rr,cc)=avgprecip(rr);
            k=k+1;
        end
    end
end
fprintf(1,'%d of %d precips have been replaced by average precip over all yrs of that day.\n',k,mp*np);
%
fprintf(file,'%d of %d precips have been replaced by average precip over all yrs of that day.\r\n',k,mp*np);
%extra check
m2=0;
for rr=1:mp  % go through the hydro days    
    for cc=1:np %go through the hydro years
        if isnan(precip(rr,cc))
            m2=m2+1;
        end
    end
end
if m2==0
    fprintf(1,'Precip data ready!\n')	
else
	fprintf(1,'ERROR: you still have %d NaNs.\n',m2)
end
%
% print to metadata and screen and plot regressions
%
months='JanFebMarAprMayJunJulAugSepOctNovDec';
fprintf(file,'Temp(glacier_met_station) = Mf*Temp(city_met_station) - Ms\r\n');
fprintf(file,'Mf & Ms have monthly values and Rsquareds\r\n');
if yrBeg~=yrSplit %don't print out if no split
    fprintf(file,'Using the weather data from %s, calendar year %4d to %4d\r\n',cell2mat(mycit1),yrBeg,yrSplit);
    fprintf(file,'       Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep    Oct    Nov    Dec  \r\n');
    fprintf(file,' Mf=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\r\n',Mf(1,1:12));
    fprintf(file,' Ms=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\r\n',Ms(1,1:12));
    fprintf(file,'Rsq=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\r\n',Rsq(1,1:12));
end
fprintf(file,'Using the weather data from %s, calendar year %4d to %4d\r\n',cell2mat(mycity2),yrSplit,nt+yrBeg);
fprintf(file,'       Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep    Oct    Nov    Dec\r\n');
fprintf(file,' Mf=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\r\n',Mf(2,1:12));
fprintf(file,' Ms=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\r\n',Ms(2,1:12));
fprintf(file,'Rsq=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\r\n',Rsq(2,1:12));
%
fprintf(1,'Temp(glacier_met_station) = Mf*Temp(city_met_station) - Ms\n');
fprintf(1,'Mf & Ms have monthly values and Rsquareds\n');
if yrBeg~=yrSplit %don't print out if no split
    fprintf(1,'Using the city of %s, calendar year %4d to %4d\n',cell2mat(mycit1),yrBeg,yrSplit);
    fprintf(1,'       Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep    Oct    Nov    Dec  \n');
    fprintf(1,' Mf=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\n',Mf(1,1:12));
    fprintf(1,' Ms=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\n',Ms(1,1:12));
    fprintf(1,'Rsq=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\n',Rsq(1,1:12));
    figure;
    hold on
    for mm=1:12 
        x=xtemp1(:,2*mm-1);
        xx=min(x)-2:max(x)+2;
        subplot(3,4,mm), plot(x,ytemp1(:,mm),'b*',xx,Mf(1,mm)*xx-Ms(1,mm).','r--');
        name=sprintf('%s: Tg = %6.3f*Tc -%6.3f in degC',months(3*mm-2:3*mm),Mf(1,mm),Ms(1,mm));
        title(name,'fontsize',10);
        ylabel('glacier','fontsize',10);
        xlabel(cell2mat(mycit1),'fontsize',10);
        %set(gca,'fontsize',10)
    end
    hold off
end
fprintf(1,'Using the city of %s, calendar year %4d to %4d\n',cell2mat(mycity2),yrSplit,nt+yrBeg);
fprintf(1,'       Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep    Oct    Nov    Dec\n');
fprintf(1,' Mf=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\n',Mf(2,1:12));
fprintf(1,' Ms=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\n',Ms(2,1:12));
fprintf(1,'Rsq=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\n',Rsq(2,1:12));
figure;
hold on
for mm=1:12
    x=xtemp2(:,2*mm-1);
    xx=min(x)-2:max(x)+2;
    subplot(3,4,mm), plot(x,ytemp2(:,mm),'b*',xx,Mf(2,mm)*xx-Ms(2,mm),'r--');
    name=sprintf('%s: Tg = %6.3f*Tc -%6.3f in degC',months(3*mm-2:3*mm),Mf(2,mm),Ms(2,mm));
    title(name,'fontsize',10);
    ylabel('glacier','fontsize',10);
    xlabel(cell2mat(mycity2),'fontsize',10);
    %set(gca,'fontsize',10)
end
hold off
%
fprintf(file,'Prec(glacier_met_station) = Pf*Prec(city_met_station)\r\n');
fprintf(file,'Pf has monthly values and Rsquareds\r\n');
if yrBeg~=yrSplit %don't print out if no split
    fprintf(file,'Using the city of %s, calendar year %4d to %4d\r\n',cell2mat(mycit1),yrBeg,yrSplit);
    fprintf(file,'       Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep    Oct    Nov    Dec  \r\n');
    fprintf(file,' Pf=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\r\n',Pf(1,1:12));
    fprintf(file,'Rsq=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\r\n',RsqP(1,1:12));
end
fprintf(file,'Using the city of %s, calendar year %4d to %4d\r\n',cell2mat(mycity2),yrSplit,nt+yrBeg);
fprintf(file,'       Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep    Oct    Nov    Dec\r\n');
fprintf(file,' Pf=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\r\n',Pf(2,1:12));
fprintf(file,'Rsq=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\r\n',RsqP(2,1:12));
%
fprintf(1,'Prec(glacier_met_station) = Pf*Prec(city_met_station)\n');
fprintf(1,'Pf has monthly values and Rsquareds\n');
if yrBeg~=yrSplit %don't print out if no split
    fprintf(1,'Using the city of %s, calendar year %4d to %4d\n',cell2mat(mycit1),yrBeg,yrSplit);
    fprintf(1,'       Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep    Oct    Nov    Dec  \n');
    fprintf(1,' Pf=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\n',Pf(1,1:12));
    fprintf(1,'Rsq=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\n',RsqP(1,1:12));
    figure;
    hold on
    for mm=1:12 
        x=xprecip1(:,mm);
        xx=min(x):max(x)+2;
        subplot(3,4,mm), plot(x,yprecip1(:,mm),'b*',xx,Pf(1,mm)*xx.','r--');
        name=sprintf('%s: Pg = %6.3f*Pc in mm',months(3*mm-2:3*mm),Pf(1,mm));
        title(name,'fontsize',10);
        ylabel('glacier','fontsize',10);
        xlabel(cell2mat(mycit1),'fontsize',10);
        set(gca,'fontsize',10)
    end
    hold off
end
fprintf(1,'Using the city of %s, calendar year %4d to %4d\n',cell2mat(mycity2),yrSplit,nt+yrBeg);
fprintf(1,'       Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep    Oct    Nov    Dec\n');
fprintf(1,' Pf=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\n',Pf(2,1:12));
fprintf(1,'Rsq=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\n',RsqP(2,1:12));
figure;
hold on
for mm=1:12
    x=xprecip2(:,mm);
    xx=min(x):max(x)+2;
    subplot(3,4,mm), plot(x,yprecip2(:,mm),'b*',xx,Pf(2,mm)*xx,'r--');
    name=sprintf('%s: Pg = %6.3f*Pc in mm',months(3*mm-2:3*mm),Pf(2,mm));
    title(name,'fontsize',10);
    ylabel('glacier','fontsize',10);
    xlabel(cell2mat(mycity2),'fontsize',10);
    %set(gca,'fontsize',10)
end

hold off
%
fclose(file);
cd('../data')
cd(glacier)
save(['Corrected_',glacier,'Temp.mat'],'temp')
save(['Corrected_',glacier,'Precip.mat'], 'precip')

if m0==0 && m2==0
	ready=1;
    
stop(here)
end
