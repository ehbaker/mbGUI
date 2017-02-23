function [ready]=fitBalGrad(glacier,balance_sites,grad)
% This function estimates the net or winter if the data is missing for a
% year or if the measurement dates are really early or late (from the
% minimum). We don't want to have to rely too much on the degree-day model
% to find the minimum amount, because our tracking of precipitation and
% melt is not that good. (Run the cummulative mass completely modelled and
% see how bad it is compared to the real data).
% So, we fit a polynomial to the balance vs altitude with a floating
% constant that changes based on the year. Then as long as we have one good
% data point per year, we can estimate the other ones.
% The balances are first from glacier-wide min days, and then we make them 
% to the site min days.
% inputs
%   glacier - 0 1 or 2 for the 2 glaciers and the test data 
%   grad - use index method==1, use gradient method==2
% outputs written to file
%   winter estimate and/or net estimate for bad year, and date supposedly
%   at 
%
dbstop if error
% switch glacier
%     case 0
        mydata = ['../data/',glacier,'/Input_',glacier,'_mb.csv'];
        myAAD = ['../data/',glacier,'/Input_',glacier,'_AAD.csv'];
        mybalgrad=['../data/',glacier,'/Output_',glacier,'BalGradReplacements.csv']; 
        mymissinb=['../data/',glacier,'/Output_',glacier,'MissinB.csv']; 
        myprecipratio =['../data/',glacier,'/Calibrated_',glacier,'PrecipRatio.csv'];
        mymeltrate = ['../data/',glacier,'/Calibrated_',glacier,'Meltrate.csv'];
        lapse=-6.5; %moist adiabatic lapse rate 

file = fopen(mybalgrad,'w');
plotver=0; %if=1, will plot each year's net data vs the line fit
plotverw=0; %if=1, will plot each year's win data vs the line fit
winlim= 60;%# days obs around max that allowed to model, use precip model=bad so fewer days
netlim= 180;%# days obs around min that allowed to model, use precip model=bad so fewer days
precipratio=load(myprecipratio); % order site A,B,C,D
meltrate=load(mymeltrate);
meltrateS=meltrate(1); 
meltrateI=meltrate(2);
mbdb = importdata(mydata);
bal_yr = str2num(cell2mat(mbdb.textdata(2:end,1)));
site0All=mbdb.textdata(2:end,2);
siteAll=char(zeros(length(site0All),100));% if ever had a name longer than 100, would have a problem
namlen=NaN*ones(length(site0All),1);
for i=1:length(site0All) % have to do this because stake names are different lengths
    nam0=cell2mat(site0All(i,:));
    namlen(i)=length(nam0);
    siteAll(i,1:namlen(i))=nam0;
end
maxNam=max(namlen);
siteAll=cellstr(siteAll(:,1:maxNam));
dates=NaN*ones(length(bal_yr),2);
for i=1:length(bal_yr)
    % spring date not allowed to be NaN but fall date is
    if strcmp(mbdb.textdata(i+1,4),'NaN')||strcmp(mbdb.textdata(i+1,4),'nan')||strcmp(mbdb.textdata(i+1,4),'NAN')
        dates(i,:)=[datenum(mbdb.textdata(i+1,3)), NaN];
    else  
        dates(i,:)=[datenum(mbdb.textdata(i+1,3)), datenum(mbdb.textdata(i+1,4))];
    end
end     
mbdb.data = [bal_yr dates mbdb.data];
my_yr=(bal_yr(1):bal_yr(end)).';
ind=find(bal_yr == my_yr(end)); %check if last year on sheet is complete
if isnan(mbdb.data(ind(1),4)) %last thing to fill in is altitudes, if NaN don't use
    my_yr=my_yr(1:end-1);
end
AADs = importdata(myAAD);
if isnan(AADs(1,2)) %some excel versions screw this up and see blanks as a line
    AADs=AADs(2:end,:);
end
bin_centers = AADs(1,2:end);
binWid0=100; %100 m for now, may change to 30 m
binWid1=(bin_centers(2)-binWid0/2-bin_centers(1))*2;
binWide=(bin_centers(end)-binWid0/2-bin_centers(end-1))*2;
bot=bin_centers(1)-binWid1/2;%glacier limits assuming in 100 m bins
top=bin_centers(end)+binWide/2;
numYrs=length(my_yr);
% initialize matrices
%
proc_index20=NaN*ones(300,numYrs); %if ever had more than 300 stakes would have an error
prlen=NaN*ones(1,numYrs+1);
for i=1:numYrs % find max stakes ever have
    proc_index0 = find(bal_yr == my_yr(i)); %index of data avail for next balance year
    stak=length(proc_index0);
    proc_index20(1:stak,i)=proc_index0;
    prlen(i)=stak;
end
maxStak=max(prlen);
proc_index20=proc_index20(1:maxStak,:);%cut excess
% choose annual modelled only if data missing(bal=1)
% surface type is conventional(surf=1)
% want glacier min so representative of actual day on glacier(callGlac=1)
% no int acc (int=0)
% don't want internal ablation, that would be circular (iabl=0)
% don't want to adjust sites, will do that later (adj=0)
% don't want to use bal grad fit to fill missing values, would be circular(miss=0)
bal=1;
surf=1;
callGlac=1;
miss=0;
[b_weqW,b_weqG,b_weqGW,b_weq,b_barW,b_bar,zone_weights,dateAMin,dateAMax]=calcBalance(glacier,balance_sites,my_yr,bal,surf,callGlac,miss,grad);
%
xyn=NaN*ones(maxStak,2*numYrs);
wgtn=zeros(maxStak,numYrs);
xyw=NaN*ones(maxStak,2*numYrs);
wgtw=zeros(maxStak,numYrs);
col3=ones(numYrs,3);
name=char(zeros(numYrs,4));
fallAlt=NaN*ones(maxStak);
kk=NaN*ones(maxStak);
maxkk=0;
badn=[];
badw=[];
%
for i=1:numYrs
proc_index2 = proc_index20(1:prlen(i),i);
siteM2=siteAll(proc_index2,:);
f2=prlen(i);
datespringobsM2=dates(proc_index2,1); 
datefallobsM2=dates(proc_index2,2); 
alt = mbdb.data(proc_index2,4);
% need altitude of surface of site 
% we collect the surface, so top of snow, altitude for the spring 
% (springAlt = minAlt+snow) and fall (fallAlt = minAlt+net ) and average the two to 
% get a altitude alt ,== ( minAlt+snow+(minAlt+net) )/2 , then 
% minalt = ( 2*alt-snow-net )/2,
%
for j=proc_index2(1):proc_index2(f2)
    abcind=j-proc_index2(1)+1;
    dw=datespringobsM2(abcind);
    d=datefallobsM2(abcind);    
    snow = b_weqW(abcind,i);
    net = b_weq(abcind,i);
    fallAlt(abcind) = ( 2*alt(abcind)-snow+net )/2;
    kk(abcind)=get_siteInd(siteM2(abcind,:),glacier);
    if kk(abcind)>maxkk
        maxkk=kk(abcind);
    end
    if (~isnan(mbdb.data(j,6))&&d<dateAMin(i)+netlim&&d>dateAMin(i)-netlim)||isnan(d)
    % net no NaN or fall date not really late or early
        xyn(abcind,2*i-1)=fallAlt(abcind);
        xyn(abcind,2*i)=b_weq(abcind,i); 
        wgtn(abcind,i)=kk(abcind);
    end
    if (~isnan(mbdb.data(j,5))&&dw<dateAMax(i)+winlim&&dw>dateAMax(i)-winlim)
    % win no NaN or spring date not really late or early
        xyw(abcind,2*i-1)=fallAlt(abcind);
        xyw(abcind,2*i)=b_weqW(abcind,i);
        wgtw(abcind,i)=kk(abcind);
    end    
end
if sum(wgtn(:,i))==0 %need at least one day per year, so cut some restrictions
    need=0;
    for j=proc_index2(1):proc_index2(f2)
        abcind=j-proc_index2(1)+1;
        kk(abcind)=get_siteInd(siteM2(abcind,:),glacier);
        if ~isnan(mbdb.data(j,6))
            if need==0
                badn=[badn,my_yr(i)]; %#ok<AGROW> %years that have to use late/early dates
                need=1;
            end
            xyn(abcind,2*i-1)=fallAlt(abcind);
            xyn(abcind,2*i)=b_weq(abcind,i); 
            wgtn(abcind,i)=kk(abcind);
        end
    end
end
if sum(wgtw(:,i))==0 %need at least one day per year, so cut some restrictions
    need=0;
    for j=proc_index2(1):proc_index2(f2)
        abcind=j-proc_index2(1)+1;
        kk(abcind)=get_siteInd(siteM2(abcind,:),glacier);
        if ~isnan(mbdb.data(j,5))
            if need==0
                badw=[badw,my_yr(i)];%#ok<AGROW> %years that have to use late/early dates
                need=1;
            end
            xyw(abcind,2*i-1)=fallAlt(abcind);
            xyw(abcind,2*i)=b_weqW(abcind,i); 
            wgtw(abcind,i)=kk(abcind);
        end
    end
end 
end
numData=length(bal_yr);
%last 3 are place holders for fall_alt, fixed_win, fixed_net
fixThese=[bal_yr,dates,mbdb.data(:,4),zeros*ones(numData,3)];
for i=1:numYrs
proc_index2 = proc_index20(1:prlen(i),i);
f2=prlen(i);
alt = mbdb.data(proc_index2,4);
for j=proc_index2(1):proc_index2(f2)
    abcind=j-proc_index2(1)+1;
    snow = b_weqW(abcind,i);
    net = b_weq(abcind,i);
    fallAlt(abcind) = ( 2*alt(abcind)-snow+net )/2;
    if isnan(xyw(abcind,2*i))&&isnan(xyn(abcind,2*i))%need to fix both
    % put glacier min/max dates in because that's the date the regression will fit
fixThese(j,:)=[my_yr(i),dateAMax(i),dateAMin(i),alt(abcind),fallAlt(abcind),NaN,NaN]; 
    elseif isnan(xyw(abcind,2*i))%need to fix winter only
fixThese(j,:)=[my_yr(i),dateAMax(i),dateAMin(i),alt(abcind),fallAlt(abcind),NaN,0]; 
    elseif isnan(xyn(abcind,2*i))%need to fix net only
fixThese(j,:)=[my_yr(i),dateAMax(i),dateAMin(i),alt(abcind),fallAlt(abcind),0,NaN]; 
    else
fixThese(j,:)=[my_yr(i),dateAMax(i),dateAMin(i),alt(abcind),fallAlt(abcind),0,0];         
    end
end
end
%
for season=1:2
if season==1
    xy=xyn;
    ylab='Annual balance (m w.e.)';
    wgt=wgtn;
    tlab='Annual';
else 
    xy=xyw;
    ylab='Winter balance (m w.e.)';
    wgt=wgtw;
    tlab='Winter';
end
%
xn1=NaN*ones(maxStak*numYrs,1);
yn =NaN*ones(maxStak*numYrs,1);
xn2=zeros(maxStak*numYrs,numYrs);
wgn=NaN*ones(maxStak*numYrs,1);
% if fit 4th order with constant term allowed to vary by year
% build matrices for inversion
% format X*beta=Y
% X=[x1A x1A^2 x1A^3 x1A^4 1 0 0 . . 0
%    x1B x1B^2 x1B^3 x1B^4 1 0 0 . . 0
%    x1C x1C^2 x1C^3 x1C^4 1 0 0 . . 0
%    x2A x2A^2 x2A^3 x2A^4 0 1 0 . . 0
%    x2B x2B^2 x2B^3 x2B^4 0 1 0 . . 0
%    x2C x2C^2 x2C^3 x2C^4 0 1 0 . . 0
%      .     .     .    .   .   
%      .     .     .    .     .
%    xNA xNA^2 xNA^3 xNA^4 0 0 0 . . 1
%    xNB xNB^2 xNB^3 xNB^4 0 0 0 . . 1
%    xNC xNC^2 xNC^3 xNC^4 0 0 0 . . 1]
% so there's one set of x's for each site, and for each year N=numYrs, then
% the second block is ones in column diagonal, for a maxStak*NXN matrix
% the beta vector will be the polynomial coefficients and a constant for
% each year.
% we might fit another order polynomial
%
for i=1:numYrs
    ind0=maxStak*(i-1)+1:maxStak*i;
    xn1(ind0)=xy(:,2*i-1)./1000;%put in km so fewer rounding errors
    yn(ind0)=xy(:,2*i);
    xn2(ind0,i)=ones(maxStak,1);
    wgn(ind0)=wgt(:,i);
end
% weight LS points by how many have at each site
%
wgn2=NaN*ones(maxStak*numYrs,1);
ordr=-1; %can fit n points perfectly with an n-1 order polynomial
for i=1:maxkk
    sit=find(wgn==i);
    meanAlt=mean(xn1(sit));
    sit2=[];
    sit1=[];
    sit0=[];
    for k=1:length(sit);
        if xn1(sit(k))>meanAlt+0.05 %50/1000, km
        %kind of in a separate alt bin then, not really same site 
            sit2=[sit2;sit(k)]; %#ok<AGROW>
        elseif xn1(sit(k))<meanAlt-0.05 %50/1000, km
            sit0=[sit0;sit(k)]; %#ok<AGROW>
        else
            sit1=[sit1;sit(k)]; %#ok<AGROW>
        end
    end
    if ~isempty(sit2)
        ordr=ordr+1;
    end
    if ~isempty(sit0)
        ordr=ordr+1;
    end
    if ~isempty(sit1)
        ordr=ordr+1;
    end
    wgn2(sit2)=1/length(sit2);
    wgn2(sit1) =1/length(sit1);
    wgn2(sit0)=1/length(sit0);
end
%ordr=3;
%
ind=find(~isnan(xn1));
X=ones(length(ind),ordr+numYrs);
for j=1:ordr
    X(:,j)=xn1(ind).^j;
end
W=diag(wgn2(ind));%cut out nan's
X(:,ordr+1:ordr+numYrs)=xn2(ind,:);
beta=(W*X)\(W*yn(ind));
SST= sum((W*yn(ind) - mean(W*yn(ind))).^2);
SSE= sum((W*yn(ind) - W*X*beta).^2);
Rsq=1-SSE/SST;

%stop here
%plotting data
%
figure;
namn=[];
test_b = [NaN,NaN];
for i=1:numYrs
    if i-1<=(numYrs-1)/3
        col3(i,:)=[1,3*(i-1)/(numYrs-1),0];
    elseif i-1<=2*(numYrs-1)/3
        col3(i,:)=[0,1,3*(i-1)/(numYrs-1)-1,];
    else
        col3(i,:)=[3*(i-1)/(numYrs-1)-2,0,1];
    end
    name(i,:)=sprintf('%4d',my_yr(i));  
    ind1=find(~isnan(xy(:,2*i)));
    if ~isempty(ind1)
        namn=[namn;name(i,:)]; %#ok<AGROW>
        plot(xy(ind1,2*i-1),xy(ind1,2*i),'.','Color', col3(i,:));
        test1 = [xy(ind1,2*i-1),xy(ind1,2*i)];
        test_b = [test_b;test1];
        hold on;
    end
end 
test_b;

for i=1:numYrs
    ind1=find(~isnan(xy(:,2*i)));
    if ~isempty(ind1)
        conn=[];
        for j=1:length(ind1)-1
            if ind1(j)+1==ind1(j+1) %only plot points next to each other
                conn=[conn,ind1(j),ind1(j+1)]; %#ok<AGROW>
            end
        end
        plot(xy(conn,2*i-1),xy(conn,2*i),'-','Color', col3(i,:));
        hold on;
    end
end     
stop(here)
hold off;
title([tlab,' balance gradients'],'fontsize',12);
set(gca,'fontsize',7.7);
legend(namn,'Location','EO');
ylabel(ylab,'fontsize',12)
xlabel('Altitude (m) terminus to head','fontsize',12)
xlim([bot,top]);
%plotting fitted polynomials
%
x=min(xn1(ind))-.05:.025:max(xn1(ind))+.05;% add 50 meters extrapolation
yhatex=@(ex) beta(1)*ex;
word=ones(1,2*(ordr-1));
for j=2:ordr
    yhatex=@(ex) yhatex(ex)+beta(j)*ex.^j;
    word(2*j-3:2*j-2)=[beta(j),j];
end
titl1=sprintf('const(yr)+%3.1fx',beta(1));
fmtj = repmat('+%3.1fx^%d', 1, ordr-1);
titl2=sprintf(fmtj,word);
titl=[titl1,titl2,', x in km'];
figure;
noyr=[];
yhat=yhatex(x);
if season==1
    ELA=NaN*ones(numYrs,1);
end
for i=1:numYrs 
    if beta(i+ordr)~=0 % if does, then no data points for that year
        plot(x*1000,yhat+beta(i+ordr),'-','Color', col3(i,:));
        hold on;
        %find ELA in interval bot to top, the zero of the net polynomial
        if season==1
            ELA0=fzero(@(ex0) yhatex(ex0)+beta(i+ordr),(bot+top)/2000);%start at average
            for ii=1:length(ELA0)%incase more than one zero or out of range
                ELAlo=top/1000;
                if ELA0(ii)<=top/1000 && ELA0(ii)>=bot/1000 && ELA0(ii)<ELAlo
                    ELAlo=ELA0(ii); %find smallest
                end
            end
            %didn't find and all of glacier positive, not likely
            if isnan(ELA0)&&(yhatex(bot/1000)+beta(i+ordr)>=0)
                ELAlo=bot/1000;
            end
            ELA(i)=ELAlo*1000; %put in meters
        end
    else
        noyr=[noyr;my_yr(i)]; %#ok<AGROW>
    end
end
hold off;
title(titl,'fontsize',12);
set(gca,'fontsize',7.7);
legend(namn,'Location','EO');
ylabel(ylab,'fontsize',12)
xlabel('Altitude (m) terminus to head','fontsize',12)
xlim([bot,top]);
%
if (plotver==1&&season==1) ||(plotverw==1&&season==2) 
dotimes=ceil(numYrs/12); %round up
for c=1:dotimes
    figure;
    for i=1:numYrs
    if i<=12*c && i>12*(c-1)
        ind1=find(~isnan(xy(:,2*i)));
        if ~isempty(ind1)
            place=mod(i,12);
            if place==0
                place=12;
            end
            subplot(3,4,place);
            plot(xy(ind1,2*i-1),xy(ind1,2*i),'.',x*1000,yhat+beta(i+ordr),'-','Color', col3(i,:));
            titl2=sprintf('Year %s fit',name(i,:));
            title(titl2,'fontsize',10);
            set(gca,'fontsize',10);
            ylabel(ylab,'fontsize',10)
            xlabel('Altitude (m)','fontsize',10)
            xlim([bot,top]);
        end
    end
    end
end



end
if season==1 % keep info
    noyr_net=noyr;
    beta_net=beta;
    ordr_net=ordr;
    Rsq_net=Rsq;
    equation_net=titl;
else
    noyr_win=noyr;
    beta_win=beta;
    ordr_win=ordr;
    Rsq_win=Rsq;
    equation_win=titl;
end
end
% format of fixThese is year spr_date fall_date average_alt fall_alt
% then 2 spots for the wint_bal and the net_bal
% these are 0 if don't need to be fixed, and NaN if do
% also, fall_alt only exists if needs to be fixed, else it's 0
%
%regression for winter except constant
yhatW=beta_win(1)*fixThese(:,5)/1000; %put x in km
for j=2:ordr_win
    yhatW=yhatW+beta_win(j)*(fixThese(:,5)/1000).^j; %put x in km
end
%regression for net, except constant
yhatN=beta_net(1)*fixThese(:,5)/1000; %put x in km
for j=2:ordr_net
    yhatN=yhatN+beta_net(j)*(fixThese(:,5)/1000).^j; %put x in km
end
indfix=[];%keep the record of which ones fixed for metaData and screen printing
forfile=char(zeros(numData,42));%initalize string
for k=1:numData
    yrind=fixThese(k,1)-my_yr(1);
    if isnan(fixThese(k,6))|| isnan(fixThese(k,7))      
        % need to get to site max/min date from glacier max/min date 
        % ablationmodel(yr,alt,springDate,fallDate,snowAtfallObs,
        % site, . . .callGlac=0 because want site maxs/mins
        indfix=[indfix,k]; %#ok<AGROW>
        %stop(here)
[win_date, win_adj, net_date, net_adj]=ablationmodel(fixThese(k,1),fixThese(k,4),fixThese(k,2),fixThese(k,3),NaN,siteAll(k),glacier,lapse,meltrateS,meltrateI,precipratio,0); 
        proc_index2 = proc_index20(1:prlen(yrind-1),yrind-1);
        ilast00=proc_index2;
        %for h=1:maxNam
            ilasta=ilast00;
            %ilast00=find(siteAll(ilast00,h)==siteAll(k,h)); %index of site for last balance year  
            ilast00=strcmp(siteAll(ilast00,1),siteAll(k,1)); %index of site for last balance year  
            ilast0=ilasta(ilast00);
        %end
        datenet=NaN;
        theAlt=fixThese(k,4);
        if ~isempty(ilast0)
            ilast= ilast0; %index of site for last balance year in whole vector
%             if maxNam>1
%                 ilast =ilast0+proc_index2(1)-1;
%             end
            datenet=b_weq(ilast00,yrind-1);
            theAlt=fixThese(ilast,4);
            springd=fixThese(ilast,2);
            falld=fixThese(ilast,3);
        elseif isempty(ilast0)&&~isempty(proc_index2)
            ilast=proc_index2(1);
            springd=fixThese(ilast,2);
            falld=fixThese(ilast,3);
        else %bad situation, means missing data on first year collected. . . have to ignore beginning correction
            springd=datenum(my_yr(yrind)-1,4,1);
            falld=NaN;
        end  
%         stop(here)
[win0_date, win0_adj,net0_date, net0_adj]=ablationmodel(my_yr(yrind-1),theAlt,springd,falld,datenet,siteAll(k),glacier,lapse,meltrateS,meltrateI,precipratio,0); 
        if isnan(fixThese(k,6)) %need to fix winter
            if isempty(find(noyr_win==fixThese(k,1),1))% if has a constant
                glacMax=yhatW(k)+beta_win(ordr_win+yrind);
                fixThese(k,2)=win_date; 
                fixThese(k,6)=glacMax+win_adj-net0_adj;
            else %can't fix
                fixThese(k,6)=0;
            end
        end
        if isnan(fixThese(k,7)) %need to fix net
            if isempty(find(noyr_net==fixThese(k,1),1))% if has a constant
                glacMin=yhatN(k)+beta_net(ordr_net+yrind);
                fixThese(k,3)=net_date; 
                fixThese(k,7)=glacMin+net_adj-net0_adj;
            else %can't fix
                fixThese(k,7)=0;
            end
        end
    end
    % keep year spr_date fall_date fixed_win fixed_net ELA for file 
    if fixThese(k,1)>my_yr(end) || fixThese(k,1)<my_yr(1)
        ELA(yrind)=NaN; %no ELA computed
    end
forfile(k,:)=sprintf('%4d %7d %7d %7.4f %7.4f %5.0f',fixThese(k,1),fixThese(k,2:3),fixThese(k,6:7),ELA(yrind));
end
spr_date=datestr(fixThese(indfix,2),2);%make in month/day/yr format
fall_date=datestr(fixThese(indfix,3),2);
forprint=char(zeros(length(indfix),55+maxNam));%initalize string
% keep year site spr_date fall_date fixed_win fixed_net for metaData andscreen
for k=1:length(indfix)
    m=sprintf('%s',cell2mat(siteAll(indfix(k),:)));
forprint(k,1:(55+length(m)))=sprintf('%4d  %s   %s   %s %10.4f %15.4f',fixThese(indfix(k),1),cell2mat(siteAll(indfix(k),:)),spr_date(k,:),fall_date(k,:),fixThese(indfix(k),6:7));
end
%
% print to screen and metadata file
%
if sum(wgtn(:,i))==0 %all sites are NaNs, this would be bad
fprintf(1,'There are no viable net balance measurements from the year %d sites.\n',my_yr(i));
fprintf(1,'This net will have to be completely modelled with the degree day model.\n');
fprintf(file,'There are no viable net balance measurements from the year %d sites.\r\n',my_yr(i));
fprintf(file,'This net will have to be completely modelled with the degree day model.\r\n');
end
if sum(wgtw(:,i))==0 %all sites are NaNs, this would be bad
fprintf(1,'There are no viable winter balance measurements from the year %d sites.\n',my_yr(i));
fprintf(1,'This winter will have to be completely modelled with the degree day model.\n');
fprintf(file,'There are no viable winter balance measurements from the year %d sites.\r\n',my_yr(i));
fprintf(file,'This winter will have to be completely modelled with the degree day model.\r\n');
end    
%
if ~isempty(badn)||~isempty(badw)
fmtbadn = repmat('%4d ', 1, length(badn));
fmtbadw = repmat('%4d ', 1, length(badw));
fprintf(1,'To get at least one point per year to locate the balance gradient curve\n');
fprintf(1,'for some years we had to use observations made on days more than +-60 from\n'); 
fprintf(1,'the min/max (so these 60 days were modelled with the degree day model.\n');
fprintf(1,'These years are:\n');
fprintf(file,'To get at least one point per year to locate the balance gradient curve\r\n');
fprintf(file,'for some years we had to use observations made on days more than +-60 from\r\n'); 
fprintf(file,'the min/max (so these 60 days were modelled with the degree day model.\r\n');
fprintf(file,'These years are:\r\n');
    if ~isempty(badn)
        fprintf(1,'[');
        fprintf(1, fmtbadn,badn);
        fprintf(1,'] for the net\r\n');
        fprintf(file,'[');
        fprintf(file, fmtbadn,badn);
        fprintf(file,'] for the net\r\n');
    end
    if ~isempty(badw)
        fprintf(1,'[');
        fprintf(1, fmtbadw,badw);
        fprintf(1,'] for the winter\r\n');
        fprintf(file,'[');
        fprintf(file, fmtbadw,badw);
        fprintf(file,'] for the winter\r\n');
    end
end
fprintf(1,'Fit the weighted LS polynomial to the net balances:\n');
fprintf(1,'%s\n',equation_net);
fprintf(1,'The constant term varies yearly. Rsquared is %6.4f\n',Rsq_net);
fprintf(1,'Fit the weighted LS polynomial to the winter balances:\n');
fprintf(1,'%s\n',equation_win);
fprintf(1,'The constant term varies yearly. Rsquared is %6.4f\n',Rsq_win);
fprintf(file,'Fit the weighted LS polynomial to the net balances:\r\n');
fprintf(file,'%s\r\n',equation_net);
fprintf(file,'The constant term varies yearly. Rsquared is %6.4f\r\n',Rsq_net);
fprintf(file,'Fit the weighted LS polynomial to the winter balances:\r\n');
fprintf(file,'%s\r\n',equation_win);
fprintf(file,'The constant term varies yearly. Rsquared is %6.4f\r\n',Rsq_win);
%
fprintf(1,'NONZERO BALANCES ARE COMPUTED FROM THE FITTED BALANCE GRAD AT\n');
fprintf(1,'GLACIER MIN/MAX DATE & ADJUSTED TO SITE MIN/MAX DATE. THESE\n');
fprintf(1,'WILL REPLACE UNACCEPTABLE VALUES IN MEASURED SITE BALANCES\n');
fprintf(1,'year site spr_date  fall_date fitted_wint_bal fitted_net_bal\n');
for i=1:length(forprint(:,1))
    fprintf(1,'%s\n', forprint(i,:).');
end
fprintf(file,'NONZERO BALANCES ARE COMPUTED FROM THE FITTED BALANCE GRAD AT\r\n');
fprintf(file,'GLACIER MIN/MAX DATE & ADJUSTED TO SITE MIN/MAX DATE. THESE\r\n');
fprintf(file,'WILL REPLACE UNACCEPTABLE VALUES IN MEASURED SITE BALANCES\r\n');
fprintf(file,'year site spr_date  fall_date fitted_wint_bal fitted_net_bal\r\n');
for i=1:length(forprint(:,1))
    fprintf(file,'%s\r\n', forprint(i,:).');
end
fclose(file);

%hell = wtf(xarrgh);
%
% print to file to use in computing balances
file=fopen(mymissinb,'w');   
for i=1:length(forfile(:,1))
    fprintf(file,'%s\r\n', forfile(i,:).');
end
ready=1;