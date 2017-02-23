function [] = fixCityData(glacier)
% this function fixes the city temperature and precipitation data by
% putting NaNs in the missing values
%   glacier - 0 1 or 2 for the 2 glaciers and the test data
myglacier=glacier;  %Selected Glacier
        mycitydata=importdata(['../data/',glacier,'/Input_',glacier,'CityData.txt']);  %import selected glacier city data
        mycit1=mycitydata.textdata(1); %Longest city data location if second/better location does not cover entire glacier record
        yrBeg=mycitydata.data(1);  %First year of city data
        mycity2=mycitydata.textdata(2);     %Second/best city data location
        yrSplit=mycitydata.data(2);

        mycity=cell2mat(['../data/',glacier,'/CityData_wx/Input_',mycity2,'data.txt']);
        outdata =cell2mat(['../data/',glacier,'/CityData_wx/Input_',mycity2,'data.mat']);
        outdates=cell2mat(['../data/',glacier,'/CityData_wx/Input_',mycity2,'dates.mat']); 


% AFTER DOWNLOADING .csv FILE FROM UTAH CLIMATE CENTER, YOU MUST MAKE ALL 
% CELLS HAVE FORMAT "NUMBER" OR "GENERAL"
city = importdata(mycity);

% I think Excel (2007) uses '30-Dec-1899' as day 1 but Excel (2004) 
% uses '1-Jan-1904' as day 1, Matlab uses '1-Jan-0000' as day 1
% Depending on which you save your data file in, you get a different day
% Respectively, Oct 1 1964 should be day 23651 or 22189 or 717611 
% We want day1 to be 717611
%IF YOU ADD DATA, MAKE SURE THE DATES ARE CONSECUTIVE! YOU MAY HAVE
%SWITCHED EXCEL VERSIONS BETWEEN THE UPDATES TO THE DATA
date = datenum(city.textdata);
citydata = [date city.data(1:end,1:3)];
% these data are date MaxDailyT MinDailyT PrecipInches
% make year month day precip maxT minT
for j= 2:4
for k = 1:length(date)
    if citydata(k,j)== 9999
        citydata(k,j)= nan;
    end
end
end
citydata(:,5)=(5/9)*(citydata(:,2)-32);  % convert F to C
citydata(:,6)=(5/9)*(citydata(:,3)-32);  % convert F to C
citydata(:,4)=citydata(:,4)*2.54*10;  % convert inches to mm
%
yd=yearday(date);%output year,day
[Y,MO,D]=datevec(date);
citydata(:,1:3) = [Y,MO,D];
%
%Nan missing datarows
numtotal=0;
 for si=1:length(date)-1
    numtotal=numtotal+1;
    newcit(numtotal,:)=citydata(si,:); %#ok<AGROW>
    theyd(numtotal,:)=yd(si,:); %#ok<AGROW>
    missing=date(si+1)-date(si)-1; % number of missing days
    if missing>1
        missdates=date(si)+1:date(si)+missing; % vector of missing dates
        missdatevec=datevec(missdates);  %make the year, month, day vector)
        newcit(numtotal+1:numtotal+missing,1:3)=missdatevec(:,1:3);
        if length(missdates)==2 %yearday does different thing for length 2 vector
            theyd(numtotal+1:numtotal+missing,:)=[yearday(missdates(1));yearday(missdates(2))];
        else
            theyd(numtotal+1:numtotal+missing,:)=yearday(missdates.');
        end
        nanvec=[nan nan nan];
        nanvec=repmat(nanvec,missing,1);
        newcit(numtotal+1:numtotal+missing,4:6)=nanvec;
        numtotal=numtotal+missing;
     end
 end 
numtotal=numtotal+1;
newcit(numtotal,:)=citydata(length(date),:); %#ok<NASGU>
theyd(numtotal,:)=yd(length(date),:);
[thehy1,thehy2]=caltohy(theyd(:,2),theyd(:,1)); %input and output day,year
thehy=[thehy2,thehy1]; %#ok<NASGU>
%
% store corrected data in a mat file
save(outdata,'newcit'); % these data are Year Month Day Precip MaxDailyT MinDailyT in meters and C
save(outdates,'thehy','theyd'); 


end

