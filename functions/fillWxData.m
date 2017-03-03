function [ready]=fillWxData(glacier)
% this function fixes the missing glacier temperature and precip data by
% using the closest city data and estimating regression functions between
% the city data and the glacier data.
%
dbstop if error

formatSecondaryWxData(glacier); %correct/reformat data from cities
        PrimaryWx=readtable(['data/',glacier,'/Input/Input_',glacier,'_Daily_Weather.csv']); %Weather from the nearest Wx station
        SecondaryWxNames=importdata(['data/',glacier,'/Input/Input_',glacier,'_SecondaryWxData.csv']);  %file with names of secondary Wx stations
        SecondaryWx1=readtable(cell2mat(['data/',glacier,'/Input/SecondaryWxData/Output_',SecondaryWxNames(1),'data.csv'])); %
        SecondaryWx2=readtable(cell2mat(['data/',glacier,'/Input/SecondaryWxData/Output_',SecondaryWxNames(2),'data.csv'])); %
        myweather=['data/',glacier,'/Output_',glacier,'TempPrecipRegressions.txt']; % New regressions stored here
        
        file = fopen(myweather,'w');
 
PrimaryWx.date = datetime(PrimaryWx.date);
SecondaryWx1.date = datetime(SecondaryWx1.date);
SecondaryWx2.date = datetime(SecondaryWx2.date);
AllWx_t = outerjoin(PrimaryWx, SecondaryWx1,'Keys','date','MergeKeys',1); %
AllWx = outerjoin(AllWx_t, SecondaryWx2,'Keys','date','MergeKeys',1); %make 1 big table organized by date
AllWx.Properties.VariableNames = {'date' 'T_primary' 'P_primary' 'T_secondary1' 'P_secondary1' 'T_secondary2' 'P_secondary2'};

figure (); hold on
title(['Monthly temperature regressions ',glacier, 'Glacier'])
monthofyear = month(AllWx.date,'monthofyear');
for m = 1:12;
    index = find(monthofyear==m);
    tbl1 = table(AllWx.T_primary(index),AllWx.T_secondary1(index));
    lm1 = fitlm(tbl1,'linear');
    tbl2 = table(AllWx.T_primary(index),AllWx.T_secondary2(index));
    lm2 = fitlm(tbl2,'linear');
    
    mth(m) = m;
    int1(m) = lm1.Coefficients.Estimate(1);
    slp1(m) = lm1.Coefficients.Estimate(2);
    rsq1(m) = lm1.Rsquared.Ordinary;
    int2(m) = lm2.Coefficients.Estimate(1);
    slp2(m) = lm2.Coefficients.Estimate(2);
    rsq2(m) = lm2.Rsquared.Ordinary;
  
    subplot(3,4,m)
    scatter(tbl1.Var1,tbl1.Var2, 'r+'); hold on
    scatter(tbl2.Var1, tbl2.Var2, 'b+'); hold on
    title(m)
    
end
  stop(here)


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
        if isnan(primary_temp(rr,cc))
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
            oldmonth = cityPT(cityInd,2);
            if ( isempty(cityInd) || isnan(minT) || isnan(maxT) ) && (rr+1<=ly && rr-1>=1)
                Tnew=(primary_temp(rr-1,cc)+primary_temp(rr+1,cc))/2; % not good, but will do for now. %if no data for city, so average of before and after
                kkk=kkk+1;
            elseif (isempty(cityInd)||isnan(minT)||isnan(maxT))&&(rr+1>ly || rr-1<1)    
                Tnew=NaN;
            else
                Tnew=Mf(yrRow,oldmonth).*mean([minT,maxT])-Ms(yrRow,oldmonth);   
                kk=kk+1;
            end
            if isnan(Tnew)&&~isempty(cityInd-1)&&~isempty(cityInd+1)&&(rr+1<=ly && rr-1>=1) %average again  
                minT=mean([cityPT(cityInd-1,6),cityPT(cityInd+1,6)]);
                maxT=mean([cityPT(cityInd-1,5),cityPT(cityInd+1,5)]);
                Tnew=Mf(yrRow,oldmonth).*mean([minT,maxT])-Ms(yrRow,oldmonth);
                kk=kk+1;
            end
            if isnan(Tnew) && (rr+2<=ly && rr-2>=1) %and again
                Tnew=(primary_temp(rr-2,cc)+primary_temp(rr+2,cc))/2;   
                kkk=kkk+1;
            end
            if isnan(Tnew)&&~isempty(cityInd-2)&&~isempty(cityInd+2)&&(rr+2<=ly && rr-2>=1) %and again
                minT=mean([cityPT(cityInd-2,6),cityPT(cityInd+2,6)]);
                maxT=mean([cityPT(cityInd-2,5),cityPT(cityInd+2,5)]);
                Tnew=Mf(yrRow,oldmonth).*mean([minT,maxT])-Ms(yrRow,oldmonth);
                kk=kk+1;
            end
            primary_temp(rr,cc)=Tnew;
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
	avgtemp(rr)=gmean(primary_temp(rr,1:nt)); 
    for cc=1:nt % go through the hydro years 
        if isnan(primary_temp(rr,cc))
            primary_temp(rr,cc)=avgtemp(rr);
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
        if isnan(primary_temp(rr,cc))
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
    primary_precip(date1hy:date2hy,date2y)=newprecip;
end
% now do like did for temp
%
yprecip1=zeros(numyrs1b*31,12);
xprecip1=zeros(numyrs1b*31,12);
yprecip2=zeros(numyrs2*31,12);
xprecip2=zeros(numyrs2*31,12);
for cc=1:np % go through the hydro years 
    for rr=1:mp % go through the hydro days
        if ~isnan(primary_precip(rr,cc))
            cityInd=find(cityHy(:,1)==cc+yrBeg & cityHy(:,2)==rr);
            P=cityPT(cityInd,4);           
            if ~isempty(cityInd)&&~isnan(P)
                %need year,oldmonth,day index to build regression matrices
                m = cityPT(cityInd,2);
                d = cityPT(cityInd,3);
                if cityPT(cityInd,1)<yrSplit
                    c=cityPT(cityInd,1)-(yrBeg-1);
                	xprecip1((c-1)*31+d,m)= P; 
                    yprecip1((c-1)*31+d,m)= primary_precip(rr,cc);
                else
                    c=cityPT(cityInd,1)-(yrSplit-1);
                	xprecip2((c-1)*31+d,m)= P; 
                    yprecip2((c-1)*31+d,m)= primary_precip(rr,cc);
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
        if isnan(primary_precip(rr,cc))
            cityInd=find(cityHy(:,1)==cc+yrBeg & cityHy(:,2)==rr);
            P=cityPT(cityInd,4);            
            %need year index and oldmonth index for regression coeffs
            if cityPT(cityInd,1)<yrSplit
                yrRow=1;
            else
                yrRow=2;
            end
            oldmonth = cityPT(cityInd,2);
            if isempty(cityInd)||isnan(P) %if no data for city
                Pnew=NaN;
            else
                Pnew=Pf(yrRow,oldmonth).*P;   
                kk=kk+1;
            end
            primary_precip(rr,cc)=Pnew;
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
	avgprecip(rr)=gmean(primary_precip(rr,1:np));
    for cc=1:np %go through the hydro years 
        if isnan(primary_precip(rr,cc))
            primary_precip(rr,cc)=avgprecip(rr);
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
        if isnan(primary_precip(rr,cc))
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
oldmonths='JanFebMarAprMayJunJulAugSepOctNovDec';
fprintf(file,'Temp(glacier_met_station) = Mf*Temp(city_met_station) - Ms\r\n');
fprintf(file,'Mf & Ms have oldmonthly values and Rsquareds\r\n');
if yrBeg~=yrSplit %don't print out if no split
    fprintf(file,'Using the weather data from %s, calendar year %4d to %4d\r\n',cell2mat(SecondaryWx1),yrBeg,yrSplit);
    fprintf(file,'       Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep    Oct    Nov    Dec  \r\n');
    fprintf(file,' Mf=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\r\n',Mf(1,1:12));
    fprintf(file,' Ms=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\r\n',Ms(1,1:12));
    fprintf(file,'Rsq=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\r\n',Rsq(1,1:12));
end
fprintf(file,'Using the weather data from %s, calendar year %4d to %4d\r\n',cell2mat(SecondaryWx2),yrSplit,nt+yrBeg);
fprintf(file,'       Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep    Oct    Nov    Dec\r\n');
fprintf(file,' Mf=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\r\n',Mf(2,1:12));
fprintf(file,' Ms=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\r\n',Ms(2,1:12));
fprintf(file,'Rsq=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\r\n',Rsq(2,1:12));
%
fprintf(1,'Temp(glacier_met_station) = Mf*Temp(city_met_station) - Ms\n');
fprintf(1,'Mf & Ms have oldmonthly values and Rsquareds\n');
if yrBeg~=yrSplit %don't print out if no split
    fprintf(1,'Using the city of %s, calendar year %4d to %4d\n',cell2mat(SecondaryWx1),yrBeg,yrSplit);
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
        name=sprintf('%s: Tg = %6.3f*Tc -%6.3f in degC',oldmonths(3*mm-2:3*mm),Mf(1,mm),Ms(1,mm));
        title(name,'fontsize',10);
        ylabel('glacier','fontsize',10);
        xlabel(cell2mat(SecondaryWx1),'fontsize',10);
        %set(gca,'fontsize',10)
    end
    hold off
end
fprintf(1,'Using the city of %s, calendar year %4d to %4d\n',cell2mat(SecondaryWx2),yrSplit,nt+yrBeg);
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
    name=sprintf('%s: Tg = %6.3f*Tc -%6.3f in degC',oldmonths(3*mm-2:3*mm),Mf(2,mm),Ms(2,mm));
    title(name,'fontsize',10);
    ylabel('glacier','fontsize',10);
    xlabel(cell2mat(SecondaryWx2),'fontsize',10);
    %set(gca,'fontsize',10)
end
hold off
%
fprintf(file,'Prec(glacier_met_station) = Pf*Prec(city_met_station)\r\n');
fprintf(file,'Pf has oldmonthly values and Rsquareds\r\n');
if yrBeg~=yrSplit %don't print out if no split
    fprintf(file,'Using the city of %s, calendar year %4d to %4d\r\n',cell2mat(SecondaryWx1),yrBeg,yrSplit);
    fprintf(file,'       Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep    Oct    Nov    Dec  \r\n');
    fprintf(file,' Pf=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\r\n',Pf(1,1:12));
    fprintf(file,'Rsq=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\r\n',RsqP(1,1:12));
end
fprintf(file,'Using the city of %s, calendar year %4d to %4d\r\n',cell2mat(SecondaryWx2),yrSplit,nt+yrBeg);
fprintf(file,'       Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep    Oct    Nov    Dec\r\n');
fprintf(file,' Pf=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\r\n',Pf(2,1:12));
fprintf(file,'Rsq=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\r\n',RsqP(2,1:12));
%
fprintf(1,'Prec(glacier_met_station) = Pf*Prec(city_met_station)\n');
fprintf(1,'Pf has oldmonthly values and Rsquareds\n');
if yrBeg~=yrSplit %don't print out if no split
    fprintf(1,'Using the city of %s, calendar year %4d to %4d\n',cell2mat(SecondaryWx1),yrBeg,yrSplit);
    fprintf(1,'       Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep    Oct    Nov    Dec  \n');
    fprintf(1,' Pf=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\n',Pf(1,1:12));
    fprintf(1,'Rsq=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\n',RsqP(1,1:12));
    figure;
    hold on
    for mm=1:12 
        x=xprecip1(:,mm);
        xx=min(x):max(x)+2;
        subplot(3,4,mm), plot(x,yprecip1(:,mm),'b*',xx,Pf(1,mm)*xx.','r--');
        name=sprintf('%s: Pg = %6.3f*Pc in mm',oldmonths(3*mm-2:3*mm),Pf(1,mm));
        title(name,'fontsize',10);
        ylabel('glacier','fontsize',10);
        xlabel(cell2mat(SecondaryWx1),'fontsize',10);
        set(gca,'fontsize',10)
    end
    hold off
end
fprintf(1,'Using the city of %s, calendar year %4d to %4d\n',cell2mat(SecondaryWx2),yrSplit,nt+yrBeg);
fprintf(1,'       Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep    Oct    Nov    Dec\n');
fprintf(1,' Pf=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\n',Pf(2,1:12));
fprintf(1,'Rsq=[%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f]\n',RsqP(2,1:12));
figure;
hold on
for mm=1:12
    x=xprecip2(:,mm);
    xx=min(x):max(x)+2;
    subplot(3,4,mm), plot(x,yprecip2(:,mm),'b*',xx,Pf(2,mm)*xx,'r--');
    name=sprintf('%s: Pg = %6.3f*Pc in mm',oldmonths(3*mm-2:3*mm),Pf(2,mm));
    title(name,'fontsize',10);
    ylabel('glacier','fontsize',10);
    xlabel(cell2mat(SecondaryWx2),'fontsize',10);
    %set(gca,'fontsize',10)
end

hold off
%
fclose(file);
cd('data')
cd(glacier)
save(['Corrected_',glacier,'Temp.mat'],'temp')
save(['Corrected_',glacier,'Precip.mat'], 'precip')

if m0==0 && m2==0
	ready=1;
    
stop(here)
end
