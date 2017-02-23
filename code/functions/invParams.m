function [ready,usebetap]= invParams(glacier) 
% this function inverts for precip catch ratios and follows the method of 
% Shea et al 2009 to invert for snow & ice meltrates in m/K/day
% inputs 
%   glacier - 0 1 or 2 for the 2 glaciers and the test data
%   
% outputs to file and/or screen
%   ks - meltrate snow in m/degC/day
%   ki - meltrate ice in m/degC/day
%   Rsq - Rsquared for meltrates (R is correlation coeff), want close to 1
%   n - number of data inverted off for meltrates 
%   m - number of data inverted off for meltrates that are above the ELA, 
%   or melting only snow, not ice
%   betap - ratio between actual precip and gage catch at each site
%   Rsqp - Rsquared for ratios (R is correlation coeff), want close to 1
%   np - number of data inverted off for ratios
%
% outputs used for ploting
%   PDD - degree days for melt
%   bw - winter mass
%   a_bs - absolute value of summer melt
%   PrDD - precipitation days for each site MORE PROPERLY called catch
%   prec - precipitation at each site
%   name - the name of each site
%
% outputs
%   useks,useki,usebetap - ks,ki,betaABCD will compare to previous ks,ki,betaABCD
%   which are named meltrateS,meltrateI,precipratio
%
% For Wolverine previously they used
% 	meltrateS=-0.0045; % m/Cday snow
% 	meltrateI=-0.0045; % m/Cday ice
% 	precipratio(1)=1.41;  %site A or LA
% 	precipratio(2)=2.71;  %site B
% 	precipratio(3)=3.69;  %site C
% For Gulkana previously they used
% 	meltrateS=-0.0053; % m/Cday snow
% 	meltrateI=-0.0053; % m/Cday ice
% 	precipratio(1)=1.30?;  %site A
% 	precipratio(2)=1.30?;  %site B
% 	precipratio(3)=1.47?;  %site C
%   precipratio(4)=1.47?;   %site D
%   precipratio(5)=NaN;   %site L or LA
dbstop if error

disp('%%%%%%%%%%%%%%%Inverting Melt and Precip Ratios%%%%%%%%%%%%%%%%%')
warning off MATLAB:rankDeficientMatrix
warning off MATLAB:ezplotfeval:NotVectorized

        mydata =['../data/',glacier,'/Input_',glacier,'_mb.txt'];
        mymeltrate=['../data/',glacier,'/Calibrated_',glacier,'MeltRate.txt'];
        myprecipratio=['../data/',glacier,'/Calibrated_',glacier,'PrecipRatio.txt'];
        mysites=['../data/',glacier,'/Input_',glacier,'Sites.txt'];
        myparams=['../data/',glacier,'/Output_',glacier,'Parameters.txt']; 
        lapse=-6.5; %moist adiabatic lapse rate 

%
mbdb = importdata(mydata);                 %import stake data
mysites=importdata(mysites);
% reset degree day factor and we will calibrate a new coefficent  
bal_yr = str2num(cell2mat(mbdb.textdata(2:end,1)));      %get vector of all mass balance years in mb_input.txt file
site0All=mbdb.textdata(2:end,2);                         %get vector of all site names in mb_input.txt file
siteAll=char(zeros(length(site0All),100));               % create empty character matrix.if ever had a name longer than 100, would have a problem
namlen=NaN*ones(length(site0All),1);                     % create vector of nans 
for i=1:length(site0All)                    % get length of each site name
    nam0=cell2mat(site0All(i,:));                        % get the site name
    namlen(i)=length(nam0);                              % get the length of each site
    siteAll(i,1:namlen(i))=nam0;                         % get each site name
end
maxNam=max(namlen);                                      %get the maximum name length
siteAll=siteAll(:,1:maxNam);                             % put all name strings into a character vector
dates=NaN*ones(length(bal_yr),2);
for i=1:length(bal_yr)                      %get all dates
    % spring date not allowed to be NaN but fall date is
    if strcmp(mbdb.textdata(i+1,3),'NaN')||strcmp(mbdb.textdata(i+1,3),'nan')||strcmp(mbdb.textdata(i+1,3),'NAN')
        fprintf(1,'ERROR: year site %4d %s, do not NaN the spring date.\n',bal_yr(i),siteAll(i,:))
		fprintf(1,'Enter a date such as 4/1/year and put NaN for the winter net.\n')
    elseif strcmp(mbdb.textdata(i+1,4),'NaN')||strcmp(mbdb.textdata(i+1,4),'nan')||strcmp(mbdb.textdata(i+1,4),'NAN')
        dates(i,:)=[datenum(mbdb.textdata(i+1,3)), NaN];
    else  
        dates(i,:)=[datenum(mbdb.textdata(i+1,3)), datenum(mbdb.textdata(i+1,4))];
    end
end      
mbdb.data = [bal_yr dates mbdb.data];  %replace original data with properly formated data from previous steps
my_yr=(bal_yr(1):bal_yr(end)).';
ind=find(bal_yr == my_yr(end)); %check if last year on sheet is complete
if isnan(mbdb.data(ind(1),2)) %last thing to fill in is altitudes, if NaN don't use
    my_yr=my_yr(1:end-1);
end
%
ready = 0;
count = 0;
countmax=10;
ksIter=zeros(countmax,1);
kiIter=zeros(countmax,1);
RsqIter=zeros(countmax,1);
readylen=100; %random number for start


while ready~=readylen && count<countmax% have too loop because parameter finding is "slightly" circular
precipratio(1:length(mysites.data(:,1)))=2; %start with initial value of one and calibrate new ration
meltrateS=-.3;% reset degree day factor and we will calibrate a new coefficent  
meltrateI=-.3; 
% load previous data
%   THESE MAKE IT SLIGHTLY CIRCULAR, but need for computing tomax and
%   somewhat for net_date_num (date mostly off temp/precip pattern, not parameters.) 
%   Variable tomax is a small number so not a huge deal.
% precipratio=load(myprecipratio); % order site A,B,C,D
% meltrate=load(mymeltrate);


%now invert for meltrates
%
% must set an upper bound on ks, optimization is pathological and you need to
% play with this value to get a small f(x) = residual (see iteration
% print out). Depending on the max, the starting point is set, and the
% success of finding the actual global minimum depends on the starting
% point. For example, Wolverine & Gulkana 1966-2008 data do best with a  
% max around 0.01; too big and ks gets pushed up to its max with a 
% large f(x); too small and there's nothing interior the bounds and the
% optimization fails
maxbd=0.01;              
%
% From Shea et al 2009
% for bw< ks*PDD (1):abs(bs)= ki*PDD + (1-ki/ks)*bw +((ki/ks)*e1+e2)
% for bw>=ks*PDD (2):abs(bs)= ks*PDDs +e1,above ELA
% where bs is summer balance and bw is winter balance between observation dates
% PDD and PDDs is total degree days and snow melt degree days
%   (PDDi+PDDs=PDD) for PDDi ice melt degree days
% ki and ks are ice and snow meltrates
% e1 and e2 are error terms
% DO: Multiple regression for (1), dependents are PDD and bw, independent
% abs(bs)= a_bs)-- use data from below the ELA (site A)
%
% set indices
maxStak=0;
numYrs=length(my_yr);
proc_index20=NaN*ones(300,numYrs); %if ever had more than 300 stakes would have an error
prlen=NaN*ones(1,numYrs);
for i=1:numYrs % find max stakes ever have
    proc_index0 = find(bal_yr == my_yr(i)); %index of data avail for next balance year
    stak=length(proc_index0);
    if stak>maxStak
        maxStak=stak;
    end
    proc_index20(1:stak,i)=proc_index0;
    prlen(i)=stak;
end
% one year before
proc_indexA = find(bal_yr == my_yr(1)-1); %index of data avail for balance year
% collect
prlen(numYrs+1)=stak;
proc_index20=proc_index20(1:maxStak,:); %don't store all the extra stuff
% initialize
%
PDD=[];
bw=[];
a_bs=[];
for i=1:numYrs
proc_index = proc_index20(1:prlen(i),i); %index of data avail for balance year
zees = mbdb.data(proc_index,4); %index site altitudes needed for lapses
f = length(proc_index);
siteM=siteAll(proc_index,:);
datespringobsM=dates(proc_index,1);
datefallobsM=dates(proc_index,2);
[win_date_numM, win_acc_adjM, net_date_numM]=ablationmodel(my_yr(i),zees,datespringobsM,datefallobsM,mbdb.data(proc_index,6),siteM,glacier,lapse,meltrateS,meltrateI,precipratio,0);
%
for j=1:f %no temp data for 1964 or past end of 2008
    if isnan(datefallobsM(j)) % net at mass min, need to store datefallobs
    % this won't too adversely affect calculations because date mostly depends
    % on precip and temp pattern, not on meltrates and precip ratios
         datefallobsM(j)=net_date_numM(j);
    end
    if ~isnan(mbdb.data(proc_index(j),5))&&~isnan(mbdb.data(proc_index(j),6))
degDays=getXYmelt(my_yr(i),zees(j),datespringobsM(j),datefallobsM(j),glacier,lapse);
        if win_date_numM(j)>datespringobsM(j)% snow after the spring date that we missed
            tomax=win_acc_adjM(j);%circular, yes, but necessary
        else
            tomax=0;
        end
        % want from degree days between spring and fall
        PDD=[PDD; degDays]; %#ok<AGROW>		
        % winter mass is winter net + the snow after the spring date that we missed
        snow = mbdb.data(proc_index(j),5);
		bw=[bw; snow+tomax]; %#ok<AGROW>
        % abs(summer mass) is snow melt + ice melt, or starting at snow
        % height melt down to fall level correcting for early snow next yr
        abs_sprToFall=snow-(mbdb.data(proc_index(j),6)-mbdb.data(proc_index(j),8));
		a_bs=[a_bs; abs_sprToFall]; %#ok<AGROW>
    end
    if ~isnan(mbdb.data(proc_index(j),6))&&~isnan(mbdb.data(proc_index(j)+f,7))&&net_date_numM(j)>datefallobsM(j)
degDays=getXYmelt(my_yr(i),zees(j),datefallobsM(j),net_date_numM(j),glacier,lapse);
        % want from degree days between fall and min
        PDD=[PDD; degDays]; %#ok<AGROW>		
        % winter mass is net
        snow = mbdb.data(proc_index(j),6);
        if snow<0
            snow=0;
        end
		bw=[bw; snow]; %#ok<AGROW>
        % abs(summer mass) is snow melt + ice melt, or starting at net
        % height melt down to minimum level
        abs_sprToFall=abs(mbdb.data(proc_index(j)+f,7));
		a_bs=[a_bs; abs_sprToFall]; %#ok<AGROW>
    end
end
end
n=length(a_bs); % number of data inverted off
% do a free knot spline with one knot (based off ks*PDD) and two regression equations
%
coefEst = @(ks0) [PDD,bw.*(bw<ks0*PDD),(bw<ks0*PDD),ks0*ones(n,1)];
beta =@(ks0) nlinfit(coefEst(ks0),a_bs,@modelmelt,0.003);
%,statset('display','iter','robust','on')
ks = fminbnd(@(ks0) sum(a_bs - modelmelt(beta(ks0),coefEst(ks0))).^2,0,maxbd);
SST= sum((a_bs - mean(a_bs)).^2);
SSE= sum((a_bs - modelmelt(beta(ks),coefEst(ks))).^2);
Rsq=1-SSE/SST;
ki = beta(ks);
%
% m= number of data that are from above ELA
m=0;
for j= 1:n
if bw(j)>=ks*PDD(j)
    m= m+1;
end
end
%
useks=-round(ks*10000)/10000; %only keep 4 decimals, and make negative
useki=-round(ki*10000)/10000;
%
%print to file for use in mass balance calcs
file = fopen(mymeltrate,'w');
fprintf(file,'%9.4f %9.4f',useks,useki);
fclose(file);
%
%
%now invert for precipratios
%
% (1):prec= p*PrDD + e1
% where prec is snow measured at glacier between days + the melt that is modelled to happen
% and PrDD is (snow) precip catch days 
% solve for ratio between precip catch and actual precip at each site
% e1 is error term
% DO: Regression for (1), dependents are PrDD, independent prec
year=char(zeros(2*numYrs,6));
PrDD=NaN*ones(2*numYrs,300); %if ever had more than 300 stakes would have a problem
prec=NaN*ones(2*numYrs,300); %if ever had more than 300 stakes would have a problem
%
% compute for first previous year
%
% one year before
proc_index = proc_indexA; %index of data avail for balance year
if ~isempty(proc_index)
    zees0 = mbdb.data(proc_index,4); %index site altitudes needed for lapses
    f=length(proc_index);
    siteM=siteAll(proc_index,:);
    datespringobsM=dates(proc_index,1);
    datefallobsM=dates(proc_index,2);
    for j=1:f
        if isnan(datefallobsM(j)) % net at mass minimum, need to store datefallobs       
[win_date, win_acc_adj, net_date]=ablationmodel(my_yr(1)-1,zees0(j),datespringobsM(j),datefallobsM(j),mbdb.data(proc_index(j),6),siteM(j,:),glacier,lapse,useks,useki,precipratio,0);
            datefallobsM(j)=net_date;
        end
    end
else
    siteM=[];
end  
% do for all years
%
name=char(zeros(300,300*(maxNam+1)));%if ever had more than 300 stakes would have a problem
nlastend=zeros(300,1);
maxkk=0;
for i=1:numYrs
proc_index2 = proc_index20(1:prlen(i),i);       %index of data avail for balance year
zees = mbdb.data(proc_index2,4);                %site altitudes needed for lapses
f2 = length(proc_index2);
siteM2=siteAll(proc_index2,:);
dateprevobsM2=NaN*ones(f2,1);
datespringobsM2=dates(proc_index2,1);
datefallobsM2=dates(proc_index2,2);
lasFalMassOut=NaN*ones(f2,1);
year(2*i-1,1:6)=sprintf('%4dw ',my_yr(i)); %for winter snow data
year(2*i  ,1:6)=sprintf('%4d+ ',my_yr(i)); %for massoutside snow data (net to fall obs)
%
[win_date_numM2, win_acc_adjM2, net_date_numM2]=ablationmodel(my_yr(i),zees,datespringobsM2,datefallobsM2,mbdb.data(proc_index2,6),siteM2,glacier,lapse,useks,useki,precipratio,0);
mass_max=win_date_numM2;
mass_min=net_date_numM2;
for j=1:f2 
    if isnan(datefallobsM2(j)) % net at mass minimum, need to store datefallobs
    % this won't too adversely affect calculations because date mostly depends
    % on precip and temp pattern, not on meltrates and precip ratios
        datefallobsM2(j)=net_date_numM2(j);
    end
    if ~isempty(siteM)
        ilast0=1:f;
       for h=1:maxNam
            ilasta=ilast0;
            ilast0=find(siteM(ilast0,h)==siteM2(j,h)); %index of site for last balance year    
            ilast=ilasta(ilast0);
        end
    else
        ilast=[];
    end
    if ~isempty(ilast)
        lasFalMassOut(j)=mbdb.data(proc_index2(1)-f+ilast,8);
        dateprevobsM2(j)=datefallobsM(ilast);
    end
    kk=get_siteInd(siteM2(j,:),glacier);
    if kk>maxkk
        maxkk=kk;
    end
    namfound=0;
    for h=1:nlastend(kk)/(maxNam+1) % number of names have so far
        % have to take max, not sum of total, because one index might be
        % the same but not the whole name
        namfound=max(namfound,sum(siteM2(j,:)==name(kk,(h-1)*(maxNam+1)+1:h*(maxNam+1)-1)));
    end
    if namfound~=maxNam %don't have this one
       name(kk,nlastend(kk)+1:nlastend(kk)+maxNam)=siteM2(j,:);
       nlastend(kk)=nlastend(kk)+maxNam+1;
    end         
    if ~isnan(mbdb.data(proc_index2(j),5)) && ~isnan(lasFalMassOut(j)) % if winter exists
precDays=getXYprec(my_yr(i),zees(j),dateprevobsM2(j),datespringobsM2(j),glacier,lapse);
        if win_date_numM2(j)<datespringobsM2(j)%the part that melted if you went past the max
            tomax=win_acc_adjM2(j); %circular, yes, but necessary
        else
            tomax=0;
        end
        % find mass from fall date to winter date, + the part that melted if you went past the max
        snow = mbdb.data(proc_index2(j),5)-lasFalMassOut(j);
        % only want from degree days between prev fall obs and spring obs
        PrDD(2*i-1,kk)=precDays; %a column for each site
        prec(2*i-1,kk)=snow+tomax;
    end       
    if ~isnan(mbdb.data(proc_index2(j),8))&&net_date_numM2(j)<datefallobsM2(j)
    % have to use my_yr(i)+1 because my_yr(i)+1-1 is the balance year of the net_date_numM2(j) date
precDays=getXYprec(my_yr(i)+1,zees(j),net_date_numM2(j),datefallobsM2(j),glacier,lapse);
        % find mass from min date to fall date
        snow = mbdb.data(proc_index2(j),8);
        % only want from degree days between years end minimum and fall obs  
        PrDD(2*i,kk)=precDays;
        prec(2*i,kk)=snow;
    end
end
% reset params
datefallobsM = datefallobsM2;
f=f2;
siteM=siteM2;
end
% invert, will use robust because precipitation particularly
% ugly in modelling (versus temperature) and so error prone
%
maxn=max(nlastend);
PrDD=PrDD(:,1:maxkk); %cut off extra stuff
prec=prec(:,1:maxkk); %cut off extra stuff
name=name(1:maxkk,1:maxn); %cut off extra stuff
betap=NaN*ones(1,maxkk);
Rsqp=NaN*ones(1,maxkk);
np=zeros(1,maxkk);
bb=ones(1,maxkk);
bad=char(zeros(maxkk,6*2*numYrs)); %if all years were bad would be full
%stop(here)
for j=1:maxkk
    ind=find(~isnan(prec(:,j))); % where our values live

    if length(ind) == 1 
        betap(j)=nanmean(prec(:,j)./PrDD(:,j));
        Rsqp(j)=1;
        
%         bb(j)=bb(j)+1;
    elseif length(ind) >1 && length(ind) < 3
        mdl=fitlm(PrDD(:,j),prec(:,j));
        betap(j)=table2array(mdl.Coefficients(2,1));
        Rsqp(j)=mdl.Rsquared.Ordinary;
    elseif length(ind) >= 3%then we have values to invert
        np(j)=sum(~isnan(prec(:,j))); % number of data can invert off, that aren't NaNs
  
        [betap(j),stats]=robustfit(PrDD(:,j),prec(:,j),'bisquare',4.685,'off') %default but turn constant off
        w2=stats.w;
        for i=1:2*numYrs
            if w2(i)<0.5
                bad(j,6*bb(j)-5:6*bb(j))=year(i,:);
                bb(j)=bb(j)+1;
            end
        end
        SST= sum(((prec(ind,j) - mean(prec(ind,j))).*w2(ind)).^2);
        SSE= sum(((prec(ind,j) - betap(j)*PrDD(ind,j)).*w2(ind)).^2);
        Rsqp(j)=1-SSE/SST;
    else %could happen if didn't collect enough years of data for a site
fprintf(1,'ERROR: Site %s does not have enough data to give a precipratio, thus does not\n',name(j,:));
fprintf(1,'       deserve its own index in get_siteInd.m. Combine with the nearest site.\n');
        np(j)=NaN;
        betap(j)=NaN;
        bad(j,1:3)='NaN';
        bb(j)=1.5;
        Rsqp(j)=NaN;
    end
end
usebetap=round(betap*100)/100; %only keep 2 decimals
readylen=maxkk;
%
%print to file for use in mass balance calcs
fmtkk = repmat('%6.2f ', 1, maxkk);
file = fopen(myprecipratio,'w');
fprintf(file,fmtkk,usebetap);
fclose(file);
%
% while loop stuff
r = 1;
%if useks==meltrateS && useki==meltrateI 
%     for nn=1:readylen
%         if usebetap(nn)==precipratio(nn)
%             r=r+1;
%         end
%     end
%end 
ready=readylen;
count=count+1;
ksIter(count)=useks;
kiIter(count)=useki;
RsqIter(count)=round(Rsq*10000)/10000;%look at first 4 decimals
if count>2
%we've been here before, stuck in a loop. this is very likely
%to happen because we won't change the parameters very much
%each year we add data, so we are essentially starting at the
%solution, which is a very bad place to start
    ind1=find(ksIter(1:count-1)==ksIter(count));
    ind2=find(kiIter(ind1)==kiIter(count));
    ind3=find(RsqIter(ind1(ind2))==RsqIter(count),1);
    if ~isempty(ind3) 
        bestRsq=max(RsqIter(ind1(ind2(ind3)):count));%ignore part before loop
        if bestRsq==RsqIter(count)%get out, else hopefully you can next time if this isn't the max 
            ready=readylen;
        end                   
    end
end
end %exit while loop
if ready==readylen        
    % plot meltrates
    %
    figure
    x1 = PDD;
    x2 = bw;
    x3 = a_bs;
    xx1 =min(x1);
    xx2 =min(x2);
    hold on
ezmesh('s','t',@(s,t)modelmelt(-useki,[s,t*(t<-useks*s),(t<-useks*s),-useks]),[xx1-100,1500],[xx2-2,10]), view([20,20]);
    scatter3(x1,x2,x3,20,'ko','filled');
    %axis tight
    ylabel('bw (m w.e.)','fontsize',12);
    xlabel('PDD','fontsize',12);
    nameM=sprintf('ks, ki =%8.4f,%8.4f in m/degC/day',useks,useki);
    title(nameM,'fontsize',12);
    set(gca,'fontsize',12);
    grid on
    zlabel('abs(bs) (m w.e.)','fontsize',12);
    hold off
    % plot precip ratios
    %
    figure;
    maxkk=length(name(:,1));
    cols=round(maxkk/2); %will round up
    for j=1:maxkk
        ind=find(~isnan(prec(:,j))); % where our values live
        if ~isempty(ind)%then we have values to plot
subplot(2,cols,j), plot(PrDD(ind,j),prec(ind,j),'b*',PrDD(ind,j),PrDD(ind,j)*usebetap(j),'r--');
            nameP0='';
            for h=1:length(name(j,:))/(maxNam+1)
nameP0=sprintf('%s %s',nameP0,name(j,(h-1)*(maxNam+1)+1:h*(maxNam+1)-1));
            end
            nameP=sprintf('%s ratio =%5.2f',nameP0,usebetap(j));
            title(nameP,'fontsize',12);
            set(gca,'fontsize',12);
            ylabel('prec (m)','fontsize',12);
            xlabel('catch (m)','fontsize',12);
        end
    end
end
%print to metadata and screen
%
file = fopen(myparams,'w');
fprintf(file,'MELTRATES FROM PIECEWISE(LINEAR) LEAST SQUARES WITH FREE-KNOT SPLINE\r\n');
fprintf(file,'meltrateS_ks  meltrateI_ki  Rsquared  #data_regressed  #data_above_ELA\r\n');
fprintf(file,'%9.4f %13.4f %11.4f %11d %16d\n',useks,useki,Rsq,n,m);
fprintf(file,'PRECIP/CATCH RATIOS FROM ITERATIVELY REWEIGHTED LEAST SQUARES\r\n'); 
fprintf(file,'site  precip/catch_ratio  Rsquared  #data_regressed\r\n');
for j=1:maxkk
    for h=1:maxn/(maxNam+1)
        fprintf(file,' %s',name(j,(h-1)*(maxNam+1)+1:h*(maxNam+1)-1));
    end
    fprintf(file,' %14.2f %15.4f %10d\r\n',usebetap(j),Rsqp(j),np(j));
end
fprintf(file,'site years_with_weight<0.5 (for precip ratio inversion)\r\n');
for j=1:maxkk
    for h=1:maxn/(maxNam+1)
        fprintf(file,' %s',name(j,(h-1)*(maxNam+1)+1:h*(maxNam+1)-1));
    end
    fprintf(file,' [%s]\r\n',bad(j,1:6*(bb(j)-1)));
end
fclose(file);
%
fprintf(1,'MELTRATES FROM PIECEWISE(LINEAR) LEAST SQUARES WITH FREE-KNOT SPLINE\n');
fprintf(1,'meltrateS_ks  meltrateI_ki  Rsquared  #data_regressed  #data_above_ELA\n');
fprintf(1,'%9.4f %13.4f %11.4f %11d %16d\n',useks,useki,Rsq,n,m);
fprintf(1,'PRECIP/CATCH RATIOS FROM ITERATIVELY REWEIGHTED LEAST SQUARES\n'); 
fprintf(1,'site  precip/catch_ratio  Rsquared  #data_regressed\n');
for j=1:maxkk
    for h=1:maxn/(maxNam+1)
        fprintf(1,' %s',name(j,(h-1)*(maxNam+1)+1:h*(maxNam+1)-1));
    end
    fprintf(1,' %14.2f %15.4f %10d\n',usebetap(j),Rsqp(j),np(j));
end
fprintf(1,'site years_with_weight<0.5 (for precip ratio inversion)\n');
for j=1:maxkk
    for h=1:maxn/(maxNam+1)
        fprintf(1,' %s',name(j,(h-1)*(maxNam+1)+1:h*(maxNam+1)-1));
    end
    fprintf(1,' [%s]\n',bad(j,1:6*(bb(j)-1)));
end