function [zone_weights,zwFake,altF] = get_weightsFake(bin_centers1,gl_area1,gl_area0,zees0,zees1)

% This function calculates the weight for each zone of a glacier with any
% number of mass balance stakes, adding a "fake stake" located at the 
% halfway elevation between each subsequent pair of stakes. Weights are 
% then used to scale measured balances, yielding the glacier wide average, 
% or specific balance. This is off a LINEAR SPLINE, OR BALANCE GRADIENT 
% METHOD BETWEEN STAKES
%
% inputs: Area Altitude Distribution information, in 100 m bins from
% terminus up to the head of the glacier
% (1) bin_centers - a vector of AAD bin median values in meters
% (2) gl_area - area of the glacier in each of the bins in (1), in km^2
% (3) gl_area0 - previous years' area of the glacier in each of the bins in (1)
% (4) zees0 - surface altitude of each index site for a given balance year
% if it is zero, means that shouldn't use that site in balance zones
% (5) zees1 - surface altitude of each index site for a given balance year;
% this has no zeros
%
% output: zone_weights - percentage of glacier in zones represented by each
%           stake subsequently used to estimate specific balance.  
%         zw_fake - the zone weights including the fake sites 
%           form=[regular stakes;extra stakes]-- so regular stakes in same
%           place in matrix in zone_weights and zwFake
%         altF - the altitude of stakes ncluding the fake sites 
%           form=[regular stakes;extra stakes]-- so regular stakes in same
%           place in matrix in zone_weights and zwFake
% For 3 stakes: Bot zone from terminus to low stake. Low zone from low 
% stake to first fake stake. Mid-low zone from first fake to mid stake. 
% Mid-high zone from mid to second fake stake. High zone from second fake
% to high stake. Top zone high to top.
%
% use for all calculations
binWid0=100; %100 m for now, may change to 30 m
zone_weights=NaN*ones(length(zees0),1); %initialize
zwFake=zeros(2*length(zees0)-1,1); %initialize
altF=zeros(2*length(zees0)-1,1); %initialize
iextra=find(zees0==0); %stakes not using in balance 
iuse=find(zees0~=0); %stakes using
zeeFake=NaN*ones(length(iuse)-1,1); %initialize
for i=1:length(iuse)-1
    zf0 = sort(zees0(iuse)); %sorts the zees in ascending order
    zeeFake(i)=( zf0(i)+zf0(i+1) )/2;
end    
%
%initial values
gl_area = gl_area1.';
bin_centers = bin_centers1.';
binWid=binWid0*ones(length(bin_centers),1);
binWid(1)=(bin_centers(2)-binWid0/2-bin_centers(1))*2;
binWid(end)=(bin_centers(end)-binWid0/2-bin_centers(end-1))*2;
bot=bin_centers(1)-binWid(1)/2;%terminus
top=bin_centers(end)+binWid(end)/2;%top of glacier
%top elev of each bin
topBins=[bot;bin_centers+binWid/2]; %fake first bin
%parse boundary bins
numBins= length(bin_centers);

%regular weights
zees=zees0(iuse);
%  
zees = sort(zees); %sorts the zees in ascending order
numSites=length(zees);
%elevation for site so can do linear interpolation between 
elev= [bot;zees(1:end);top];
%number of zones need to deal with is one less than number of elevs
numZones=length(elev)-1;
%parse boundary bins
dist=NaN*ones(numBins+1,numZones-1);
firstB=ones(1,numZones+1);
shareRow=NaN*ones(numZones-1,1);
for i=1:numZones-1
    %distance from sites to each top, one col per bot, site, and top
    dist(:,i)= elev(i+1)-topBins;
    firstB(i+1)=find(dist(:,i)<=0,1); %index first top not below site
    if firstB(i+1)-1>0 && topBins(firstB(i+1)-1)<elev(i+1) ;
        %so bottom of bin below site, will split bin, keep info
        shareRow(i)=firstB(i+1);
    end
end
for i=1:numZones-2
    if shareRow(i)==shareRow(i+1) %should have smaller bins so divide problem bin
        s=shareRow(i)-1;
        wid=binWid(s)/3;
        add=[bin_centers(s)-binWid(s)/2+wid;bin_centers(s)+binWid(s)/2-wid];
        bin_centers=[bin_centers(1:s-1);add;bin_centers(s+1:end)];
        gl_area = [gl_area(1:s-1);gl_area(s)/2;gl_area(s)/2;gl_area(s+1:end)];
        binWid=[binWid(1:s-1);2*wid;2*wid;binWid(s+1:end)];
        topBins=[bin_centers(1)-binWid(1);bin_centers+binWid/2];
        %parse boundary bins
        numBins= length(bin_centers);
        dist=NaN*ones(numBins+1,numZones-1);
        for j=1:numZones-1
            %distance from sites to each top, one col per site and top
            dist(:,j)= elev(j+1) - topBins;
            firstB(j+1)=find(dist(:,j)<=0,1); %index first top not below site
            if firstB(j+1)-1>0 && topBins(firstB(j+1)-1)<elev(j+1) 
                shareRow(j)=firstB(j+1);
            end
        end
    end
end    
firstB(numZones+1)=numBins+1; %fake last first bin (no bin)
%partitioning ratios of split bins
distAbove=topBins(shareRow) - elev(2:end-1);
midAbove=elev(2:end-1)+distAbove/2;
fracAbove=distAbove./binWid(shareRow-1);
distBelow=elev(2:end-1)-(topBins(shareRow-1));
midBelow=elev(2:end-1)-distBelow/2;
fracBelow=distBelow./binWid(shareRow-1);
for i=1:length(shareRow)
    if fracAbove(i)~=0 && fracBelow(i)~=0 %don't need extra bin if ends in right place
        s=shareRow(i)-1;
        bin_centers=[bin_centers(1:s-1);midBelow(i);midAbove(i);bin_centers(s+1:end)];
        gl_area = [gl_area(1:s-1);gl_area(s)*fracBelow(i);gl_area(s)*fracAbove(i);gl_area(s+1:end)];
        binWid=[binWid(1:s-1);binWid(s)*fracBelow(i);binWid(s)*fracAbove(i);binWid(s+1:end)];
        shareRow=shareRow+1;%added 1 index because divided bin
        firstB(i+1:end)=firstB(i+1:end)+1;%added 1 index because divided bin
    end
end          
elevDiff=elev(2:end)-elev(1:end-1);
w=zeros(numSites,1); %compile zone weights
for i=2:numZones
    for k=firstB(i-1):firstB(i)-1 %contribution of site to lower spline
        w(i-1)=w(i-1)+( (bin_centers(k)-elev(i-1))/elevDiff(i-1) )*gl_area(k)*binWid(k);
    end
    for k=firstB(i):firstB(i+1)-1 %contribution of site to upper spline
        w(i-1)=w(i-1)+( 1 - ((bin_centers(k)-elev(i))/elevDiff(i)) )*gl_area(k)*binWid(k);
    end
end
for k=firstB(1):firstB(2)-1 %an upper for site 1, from terminus to site 1
    w(1)=w(1)+( 1 - ((bin_centers(k)-elev(1))/elevDiff(1)) )*gl_area(k)*binWid(k);
end
for k=firstB(numZones):firstB(numZones+1)-1%a lower for site last, from site last to top
    w(numSites)=w(numSites)+( (bin_centers(k)-elev(numZones))/elevDiff(numZones) )*gl_area(k)*binWid(k);
end
% calculate desired weights by scaling by total
zone_weights0 = w./sum(w);
zone_weights(iuse)=zone_weights0(1:length(iuse));
zone_weights(iextra)=0; %zero extra stakes
%
%
% fake bins
%initial values
gl_area = gl_area1.';
bin_centers = bin_centers1.';
binWid=binWid0*ones(length(bin_centers),1);
binWid(1)=(bin_centers(2)-binWid0/2-bin_centers(1))*2;
binWid(end)=(bin_centers(end)-binWid0/2-bin_centers(end-1))*2;
bot=bin_centers(1)-binWid(1)/2;%terminus
top=bin_centers(end)+binWid(end)/2;%top of glacier
%top elev of each bin
topBins=[bin_centers(1)-binWid(1)/2;bin_centers+binWid/2]; %fake first bin
%parse boundary bins
numBins= length(bin_centers);
%
zees=[zees0(iuse);zeeFake];
%  
zees = sort(zees); %sorts the zees in ascending order
numSites=length(zees);
%elevation for site so can do linear interpolation between 
elev= [bot;zees(1:end);top];
%number of zones need to deal with is one less than number of elevs
numZones=length(elev)-1;
%parse boundary bins
dist=NaN*ones(numBins+1,numZones-1);
firstB=ones(1,numZones+1);
for i=1:numZones-1
    %distance from sites to each top, one col per bot, site, and top
    dist(:,i)= elev(i+1)-topBins;
    firstB(i+1)=find(dist(:,i)<=0,1); %index first top not below site
    if firstB(i+1)-1>0 && topBins(firstB(i+1)-1)<elev(i+1) ;
        %so bottom of bin below site, will split bin, keep info
        shareRow(i)=firstB(i+1);
    end
end
for i=1:numZones-2
    if shareRow(i)==shareRow(i+1) %should have smaller bins so divide problem bin
        s=shareRow(i)-1;
        wid=binWid(s)/3;
        add=[bin_centers(s)-binWid(s)/2+wid;bin_centers(s)+binWid(s)/2-wid];
        bin_centers=[bin_centers(1:s-1);add;bin_centers(s+1:end)];
        gl_area = [gl_area(1:s-1);gl_area(s)/2;gl_area(s)/2;gl_area(s+1:end)];
        binWid=[binWid(1:s-1);2*wid;2*wid;binWid(s+1:end)];
        topBins=[bin_centers(1)-binWid(1);bin_centers+binWid/2];
        %parse boundary bins
        numBins= length(bin_centers);
        dist=NaN*ones(numBins+1,numZones-1);
        for j=1:numZones-1
            %distance from sites to each top, one col per site and top
            dist(:,j)= elev(j+1) - topBins;
            firstB(j+1)=find(dist(:,j)<=0,1); %index first top not below site
            if firstB(j+1)-1>0 && topBins(firstB(j+1)-1)<elev(j+1) 
                shareRow(j)=firstB(j+1);
            end
        end
    end
end    
firstB(numZones+1)=numBins+1; %fake last first bin (no bin)
%partitioning ratios of split bins
distAbove=topBins(shareRow) - elev(2:end-1);
midAbove=elev(2:end-1)+distAbove/2;
fracAbove=distAbove./binWid(shareRow-1);
distBelow=elev(2:end-1)-(topBins(shareRow-1));
midBelow=elev(2:end-1)-distBelow/2;
fracBelow=distBelow./binWid(shareRow-1);
for i=1:length(shareRow)
    if fracAbove(i)~=0 && fracBelow(i)~=0 %don't need extra bin if ends in right place
        s=shareRow(i)-1;
        bin_centers=[bin_centers(1:s-1);midBelow(i);midAbove(i);bin_centers(s+1:end)];
        gl_area = [gl_area(1:s-1);gl_area(s)*fracBelow(i);gl_area(s)*fracAbove(i);gl_area(s+1:end)];
        binWid=[binWid(1:s-1);binWid(s)*fracBelow(i);binWid(s)*fracAbove(i);binWid(s+1:end)];
        shareRow=shareRow+1;%added 1 index because divided bin
        firstB(i+1:end)=firstB(i+1:end)+1;%added 1 index because divided bin
    end
end        
elevDiff=elev(2:end)-elev(1:end-1);
w=zeros(numSites,1); %compile zone weights
for i=2:numZones
    for k=firstB(i-1):firstB(i)-1 %contribution of site to lower spline
        w(i-1)=w(i-1)+( (bin_centers(k)-elev(i-1))/elevDiff(i-1) )*gl_area(k)*binWid(k);
    end
    for k=firstB(i):firstB(i+1)-1 %contribution of site to upper spline
         w(i-1)=w(i-1)+( 1 - ((bin_centers(k)-elev(i))/elevDiff(i)) )*gl_area(k)*binWid(k);
    end
end
for k=firstB(1):firstB(2)-1 %an upper for site 1, from terminus to site 1
    w(1)=w(1)+( 1 - ((bin_centers(k)-elev(1))/elevDiff(1)) )*gl_area(k)*binWid(k);
end
for k=firstB(numZones):firstB(numZones+1)-1%a lower for site last, from site last to top
    w(numSites)=w(numSites)+( (bin_centers(k)-elev(numZones))/elevDiff(numZones) )*gl_area(k)*binWid(k);
end
% calculate desired weights by scaling by total
zwFake0 = w./sum(w);
zwFake(iuse)=zwFake0(2*(1:length(iuse))-1);
zwFake(iextra)=0;
if ~isempty(iextra)
    if iextra(end)==length(zees0) %cut top stake for now
        iextra=iextra(1:end-1);
    end
end
zwFake(length(zees0)+iextra)=zeros(length(iextra),1);
zwFake(length(zees0)+iuse(1:end-1))=zwFake0(2*(1:length(iuse)-1)); %cut top stake for now
altF(1:length(zees0))=zees1; %use all sites altitudes
altF(length(zees0)+iuse(1:end-1))=zeeFake; %altitudes
% form=[regular stakes;extra stakes]-- so regular stakes in same place in matrix in zw and zwFake
