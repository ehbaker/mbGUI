function [b_weqW,b_weqS,b_weqG,b_weqGW,b_weq,b_barW,b_barS,b_bar,ELA,zone_weights,dateAMin,dateAMax,meltrateS,meltrateI,Zs,winter_gradients,summer_gradients,annual_gradients]=calcBalance(glacier,balance_sites,my_yr,bal,surf,grad,plot_integration,plot_ablation_model)
% This function calculates the balance every year and 
% inputs: these are all from the GUI!
%   glacier - name of glacier you have selected in the GUI
%   my_yr - years calculate over
%   bal - time-system you want to process mb over 
%   surf - conventional =1, reference =2
%   grad - use index method==1, use gradient method==2
% 
% outputs
%   b_weqW - amount at each stake per winter, matrix 3X(length(my_yr))
%   b_weq - amount at each stake per year, matrix 3X(length(my_yr))
%   b_barW - integrated winter amound either glacierwide min or site mins
%   b_bar - integrated net amount either glacierwide min or site mins
%   zone_weights - contribution of each stake to total mass
%   dateAMin - average min day per year, col vector length(my_yr)
%   dateAMax - average max day per year, col vector length(my_yr)
%   meltrateS - meltrate for snow passed through 
%   meltrateI - meltrate for ice passed through 
dbstop if error

        mydata = ['../data/',glacier,'/Input/Input_',glacier,'_Glaciological_Data.csv'];
        myAAD = ['../data/',glacier,'/Input/Input_',glacier,'_Area_Altitude_Distribution.csv'];
        mymissinb=['../data/',glacier,'/Input/Input_',glacier,'_Missing_Glaciological_Data.csv'];
        myprecipratio = ['../data/',glacier,'/Input/Calibrated_',glacier,'_Precipitation_Ratios.csv'];
        mymeltrate = ['../data/',glacier,'/Input/Calibrated_',glacier,'_Degree_Day_Factors.csv'];
        lapse=-6.5; %moist adiabatic lapse rate 


precipratio=load(myprecipratio); % load precip ratios
meltrate=load(mymeltrate); %load melt coefficents
meltrateS=meltrate(1); %Melt rate for snow
meltrateI=meltrate(2); %melt rate for ice
miss=1;
% load and format database: text out has date strings and index_site
% Original data: bal_yr stake visit_1 visit_2 z winter net winter_ablation summer_accumulation
% Balance year net, annual are from min to min, Oct to Oct, 
% so about visit_2 of previous Bal yr to visit_2 of this Bal yr
% ie for Bal yr 2005, winter is 9/04 to 5/05, summer 5/05 to 9/05.
% Measure winter abaltion at 5/05 (visit_1) and says Bal yr 04 had less net mass 
% or previous summer wasn't done by 9/04
% Measure summer accumulation at 9/05 (visit_2) and says Bal yr 05 has less net mass 
% or next winter started before 9/05
% YOU CAN ALWAYS FIND THE FIRST MINIMUM BECAUSE YOU HIT PREVIOUS SUMMER SURFACE WHEN PIT
%
% NOTE: YOU MUST SORT DATA SO IT HAS YEARS ASCENDING, AND FOR EACH YEAR THE
% LOWEST SITE IS FIRST AND THE HIGHEST SITE IS LAST (e.g. usually L A B C)
mbdb = importdata(mydata);  %import mb data
bal_yr = str2num(cell2mat(mbdb.textdata(2:end,1)));  %get list of all mb years
site0All=mbdb.textdata(2:end,2);  %list of stake names
siteAll=char(zeros(length(site0All),100)); % need to change names from cell to string (if ever had a name longer than 100, would have a problem)
namlen=NaN*ones(length(site0All),1); 
for i=1:length(site0All) % have to do this because stake names are different lengths and cell2mat conversions
    nam0=cell2mat(site0All(i,:));
    namlen(i)=length(nam0);
    siteAll(i,1:namlen(i))=nam0;
end
maxNam=max(namlen);
siteAll=siteAll(:,1:maxNam);
dates=NaN*ones(length(bal_yr),2);

% check to make sure input mb data is formatted correctly. Spring can not
% be nans. Must have place holder of April 1st , YYYY and NaN in the
% following column
for i=1:length(bal_yr)
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
mbdb.data = [bal_yr dates mbdb.data];

%%% import AAD of glacier and strip out headers and data
AADs = importdata(myAAD);
if isnan(AADs(1,2)) %some excel versions screw this up and see blanks as a line
    AADs=AADs(2:end,:);
end
bin_centers = AADs(1,2:end);
numYrs=length(my_yr);


%
maxStak=0;
proc_index20=NaN*ones(300,numYrs); %if ever had more than 300 stakes would have an error
prlen=NaN*ones(1,numYrs+1);
for i=1:numYrs % find max stakes ever have
    proc_index0 = find(bal_yr == my_yr(i)); %index of data avail for next balance year
    stak=length(proc_index0);
    if stak>maxStak
        maxStak=stak;
    end
    proc_index20(1:stak,i)=proc_index0;
    prlen(i)=stak;
    for ii=1:stak-1
        if mbdb.data(proc_index0(ii),4)>mbdb.data(proc_index0(ii+1),4)
fprintf(1,'ERROR: year %4d needs to have sites in order lowest to highest.\n',my_yr(i));
        end
    end
end
% one year before
proc_indexA = find(bal_yr == my_yr(1)-1); %index of data avail for balance year
if ~isempty(proc_indexA) %my_yr(1) is not first year of data
    proc_index=proc_indexA;
    stak0=length(proc_index);
    if stak0>maxStak
        maxStak=stak0;
    end
else
    proc_index=proc_index20(1:prlen(1),1);
    stak0=prlen(1);
end 
% go one year past
proc_index0 = find(bal_yr == my_yr(numYrs)+1); %index of data avail for next balance year
stak=length(proc_index0);
if stak>maxStak
    maxStak=stak;
end
prlen(numYrs+1)=stak;
proc_index20=proc_index20(1:maxStak,:); %don't store all the extra stuff
%
% compute for first previous year
%
zees0 = mbdb.data(proc_index,4); %index site altitudes needed for lapses
f=stak0;
f2=prlen(1);
siteM=siteAll(proc_index,1:maxNam); %stakes for this particular year

if ~isempty(proc_indexA)    % if 
    siteM2=siteAll(proc_index(f)+1:proc_index(f)+f2,1:maxNam);
    datespringobsM=dates(proc_index,1);
    datefallobsM=dates(proc_index,2);
    sprMassOut=NaN*ones(f,1);
    if plot_ablation_model==1
        promptMessage = sprintf('Do you want to see the ablation model results for the fall before first year?');
        titleBarCaption = 'Plotting ablation model';
        choice = questdlg(promptMessage,titleBarCaption,'Yes','No','Don''t know')
                if strcmpi(choice, 'Yes')
                    plot_ablation_model1 = 1;
                elseif strcmpi(choice, 'No')
                    plot_ablation_model1=0;
                end
     else
         plot_ablation_model1=0;
    end
    [win_date_numM, win_acc_adjM,win_fxd_date_num, win_fxd_date_adj, net_date_numM, net_abl_adjM, annual_date_numM, annual_abl_adjM,search_min,search_max]=ablationmodel(my_yr(1)-1,zees0,datespringobsM,datefallobsM,mbdb.data(proc_index,6),siteM,glacier,lapse,meltrateS,meltrateI,precipratio,plot_ablation_model1);
    %
    % model something if data missing, and to get dates for balance yr
    % annual is net + (end corrections)
    % end: too early or too late so subtract extra mass that realized in next 
    % yr and add extra melt that realized in spring next yr
    % fixed-date needs beginning correction and end correction
    for j=proc_index(1):proc_index(f)
        abcind=j-proc_index(1)+1;
        inext0=1:f2;
        for h=1:maxNam
            inexta=inext0;
            inext0=find(siteM2(inext0,h)==siteM(abcind,h)); %index of site for next balance year    
            inext=inexta(inext0);
        end
        if ~isempty(inext)
            sprMassOut(abcind)=mbdb.data(proc_index(f)+inext,7);
        end
        % Error checking    
        if isnan(datefallobsM(abcind)) 
            sprMassOut(abcind)=0; %just incase they didn't make these 0
            mbdb.data(j,8)=0;
            if isnan(mbdb.data(j,6))
    fprintf(1,'ERROR: year site %4d %s, if you NaN the fall date, you must know the annual.\n',my_yr(i),siteM2(abcind,:))
    fprintf(1,'Else enter a date such as 10/1/year and put NaN for the net and mass outside.\n')
            end
        end
        if isnan(mbdb.data(j,8))||isnan(sprMassOut(abcind))% then data won't be adjusted properly
            if net_date_numM(abcind)<datefallobsM(abcind) % went late
                sprMassOut(abcind)=0;
                mbdb.data(j,8)=abs(net_abl_adjM(abcind)); % extra mass positive
            else % went early or exact right time
                sprMassOut(abcind)=net_abl_adjM(abcind);
                mbdb.data(j,8)=0;
            end
        end
    end

    % correct net if have data from next year (f indices later); corrections
    % are negative 
    adjyear=my_yr(1)-1;
else %first year of site, so measured from minimum
    siteM2=siteM;
    datespringobsM=datenum(my_yr(1)-1,4,1)*ones(f,1);
    datefallobsM=NaN*ones(f,1);
     if plot_ablation_model==1
        promptMessage = sprintf('Do you want to see the ablation model results for the fall before first year?');
        titleBarCaption = 'Plotting ablation model';
        choice = questdlg(promptMessage,titleBarCaption,'Yes','No','Don''t know')
                if strcmpi(choice, 'Yes')
                    plot_ablation_model1 = 1;
                elseif strcmpi(choice, 'No')
                    plot_ablation_model1=0;
                end
     else
         plot_ablation_model1=0;
    end
[win_date_numM, win_acc_adjM, win_fxd_date_num, win_fxd_date_adj,net_date_numM, net_abl_adjM, annual_date_numM, annual_abl_adjM,search_min,search_max]=ablationmodel(my_yr(1)-1,zees0,datespringobsM,datefallobsM,NaN*ones(f,1),siteM,glacier,lapse,meltrateS,meltrateI,precipratio,plot_ablation_model1);
    adjyear=my_yr(1);%because first year, my_yr(i)-1 doesn't exist so can't use

end

% find correction from minimum (annual or to glacier wide)
switch surf

    case 1 %is off current geometry == conventional
        AAD_index = find(AADs(:,1) == adjyear);
        gl_area = AADs(AAD_index,2:end);
        
    case 2 %is off oldest geometry == reference
        gl_area = AADs(2,2:end); %first year have data
end
             
if bal ==1 %we want only stratigraphic balances
    end_abl_adjM=net_abl_adjM;
    end_date=net_date_numM;
elseif bal==2 %is from Oct 1 to Sept 30 == fixed date
    end_abl_adjM=annual_abl_adjM; %always positive
    end_date=annual_date_numM;
elseif bal==3%is from min to min
    zw = get_weights(bin_centers,gl_area,zees0); %just need yr for adjusting data, and no adjusts for this year because no data
    [glacier_net_date,glacier_win_date,min_to_glacmin]=glacierWide(my_yr(1),win_date_numM, search_max, net_date_numM, search_min,glacier,zw);
    end_abl_adjM= min_to_glacmin; %always positive
    end_date=glacier_net_date;
end

% compute for next yrs
% Initialize matrices
dateAMin =NaN*ones(1,numYrs);
dateAMax =NaN*ones(1,numYrs);

Zs = NaN*ones(maxStak,numYrs);
b_weq =zeros(maxStak,numYrs);
b_weqW=zeros(maxStak,numYrs);
stake_zs=nan*zeros(maxStak,numYrs);
net_add=zeros(maxStak,numYrs);
win_add=zeros(maxStak,numYrs);
zone_weights=zeros(maxStak,numYrs);
altitudeF=zeros(2*maxStak-1,numYrs);
if miss==1 %then looking for a balance, so should see which sites want to include                     MAY NEED TO REVISIT! Changed this back to importdata as txt2mat was fing up AU in 2010-2011....
    useSites0 = balance_sites; %doesn't read in all text files on some computers
    useSites=char(zeros(length(useSites0),maxNam));
    sitlen=NaN*ones(length(useSites0),1);

    for i=1:length(useSites0(:,1)) % have to do this because stake names are different lengths
        sit0=cell2mat(useSites0(i,1)); %only need if using importdata, else is a mat;
        sitlen(i)=length(sit0);
        useSites(i,1:sitlen(i))=sit0;
    end
    
    missinb=load(mymissinb);%load modeled missing site values produced in 'fitbalgrad.m' in file update functions
end
keepers=ones(maxStak,numYrs);%default keep all stakes, will zero ones don't want
%stop(here)
for i=1:numYrs
%
% Find balance for yr 
proc_index2 = proc_index20(1:prlen(i),i);
if miss==1
    misind= find(missinb(:,1) == my_yr(i));
else 
    missinb=zeros(prlen(i),5);%never got defined if miss~=1
    misind=1:prlen(i);
end
zees = mbdb.data(proc_index2,4); %index site altitudes needed for AAD weights and lapses
Zs(1:length(zees),i) = zees;

f3=prlen(i+1);
siteM3=siteAll(proc_index2(f2)+1:proc_index2(f2)+f3,1:maxNam); %?? 
datespringobsM2=dates(proc_index2,1);
datefallobsM2=dates(proc_index2,2); 
sprMassOut=NaN*ones(f2,1);
prevfalMass=NaN*ones(f2,1);
beg_adj0=zeros(f2,1);
beg_adjW=zeros(f2,1);
beg_adj=zeros(f2,1);
%
for j=proc_index2(1):proc_index2(f2)
    abcind=j-proc_index2(1)+1;    
    inext0=1:f3;
    for h=1:maxNam
        inexta=inext0;
        inext0=find(siteM3(inext0,h)==siteM2(abcind,h)); %index of site for next balance year    
        inext=inexta(inext0);
    end
    if ~isempty(inext) %then site not missing at its net or has late/early obs dates
        sprMassOut(abcind)=mbdb.data(proc_index2(f2)+inext,7);
    end
    % Error checking    
    if isnan(datefallobsM2(abcind))|| missinb(misind(abcind),5)~=0
        %at true site min if date is nan or if used balance grad to
        %find net
        sprMassOut(abcind)=0; 
        mbdb.data(j,8)=0;
        if isnan(mbdb.data(j,6))&&isnan(datefallobsM2(abcind))
            fprintf(1,'ERROR: year site %4d %s, if you NaN the fall date, you must know the annual.\n',my_yr(i),siteM2(abcind,:))
            fprintf(1,'Else enter a date such as 10/1/year and put NaN for the net and mass outside.\n')
        end
    end
    ilast0=1:f;
    for h=1:maxNam
        ilasta=ilast0;
        ilast0=find(siteM(ilast0,h)==siteM2(abcind,h)); %index of site for last balance year    
        ilast=ilasta(ilast0);
    end
    if ~isempty(ilast)
        dateprev=net_date_numM(ilast);
        %so can will start whole year model with that snow
        if (proc_index2(1)-f+ilast)>1 %first year data collected, has fake proc_index so not empty
            prevfalMass(abcind)=mbdb.data(proc_index2(1)-f+ilast,6);
        end
        beg_adj(abcind)=end_abl_adjM(ilast); %always positive
        beg_adjW(abcind)=beg_adj(abcind);
    end
    if isempty(ilast)|| missinb(misind(abcind),4)~=0||missinb(misind(abcind),5)~=0
    %first year of site so measured from minimum, or missing, so from minimum
        dateW=datenum(my_yr(i)-1,4,1);
         if plot_ablation_model==1
        promptMessage = sprintf('Do you want to see the ablation model results for previous year?');
        titleBarCaption = 'Plotting ablation model';
        choice = questdlg(promptMessage,titleBarCaption,'Yes','No','Don''t know')
                if strcmpi(choice, 'Yes')
                    plot_ablation_model1 = 1;
                elseif strcmpi(choice, 'No')
                    plot_ablation_model1=0;
                end
                    else
         plot_ablation_model1=0;
    end
[win_date, win_adj,win_fxd_date_num, win_fxd_date_adj, net_date, net_adj, annual_date, annual_adj]=ablationmodel(my_yr(i)-1,zees(abcind),dateW,NaN,NaN,siteM2(abcind,:),glacier,lapse,meltrateS,meltrateI,precipratio,plot_ablation_model1);  
        % net at annual, need to store date because don't know
        dateprev=net_date;
        if bal==3
        % because the glac_min for this year wasn't calculated with this
        % site so take back to other sites glac_min (if start at that date
        % correction is the negative version to get from annual_date to
        % glac_min_date
%         if (my_yr(i)-1)==1988
%             stop(here)
%         else
         if plot_ablation_model==1
        promptMessage = sprintf('Do you want to see the ablation model results for previous year?');
        titleBarCaption = 'Plotting ablation model';
        choice = questdlg(promptMessage,titleBarCaption,'Yes','No','Don''t know')
                if strcmpi(choice, 'Yes')
                    plot_ablation_model1 = 1;
                elseif strcmpi(choice, 'No')
                    plot_ablation_model1=0;
                end
                    else
         plot_ablation_model1=0;
    end
[win_dateG, win_adjG,win_fxd_date_numG, win_fxd_date_adjG, net_dateG, net_adjG]=ablationmodel(my_yr(i)-1,zees(abcind),dateW,glacier_net_date2(1),NaN,siteM2(abcind,:),glacier,lapse,meltrateS,meltrateI,precipratio,plot_ablation_model);  
            glacMin=-net_adjG; 
%         end
        else
            glacMin=0;
        end              
        if bal==1 %stratigraphic
            end_abj=net_abl_adjM;
            end_date0=net_date_numM;
        elseif bal==2 %fixed date time-system
            end_adj =annual_adj; %always positive
            end_date0=annual_date;
                
        elseif bal==3 %combined time-system
            
            end_adj= glacMin; %always positive
            end_date0=glacier_net_date(1);
        end

    end
    if isempty(ilast)
        beg_adj(abcind)=beg_adj0(abcind);
        beg_adjW(abcind)=beg_adj0(abcind);
    end
    
    if missinb(misind(abcind),4)~=0 %missing winter
        beg_adjW(abcind)=beg_adj0(abcind);
        datespringobsM2(abcind)=missinb(misind(abcind),2);
        mbdb.data(proc_index2(abcind),5)=missinb(misind(abcind),4);
    end
    if missinb(misind(abcind),5)~=0 %missing net
        beg_adj(abcind)=beg_adj0(abcind);
        datefallobsM2(abcind)=missinb(misind(abcind),3);
        mbdb.data(proc_index2(abcind),6)=missinb(misind(abcind),5);
    end
    
    if (missinb(misind(abcind),4)==0&&isnan(mbdb.data(j,5)))||(missinb(misind(abcind),5)==0&&isnan(mbdb.data(j,6)))
    % ONLY HAPPENS IF HAVE NO POINTS FOR A YEAR, ELSE WILL USE BALANCE
    % GRADIENT
    % have to model whole year, and correct to start at minimum
    % subtract beginning adjustment(negative) to get to reported net
    % beg corr should always be positive because annual net larger if start at true 
    % minimum than at a point larger than minimum

        if isnan(datefallobsM2(abcind)) 
            fakedate=datenum(my_yr(i),9,30);
[net_mod,win_mod]=modelwhole(my_yr(i),zees(abcind),dateprev,datespringobsM2(abcind),fakedate,prevfalMass(abcind),siteM2(abcind,:),glacier,lapse,meltrateS,meltrateI,precipratio);
        %  net can't be NaN if date is Nan, but winter could be. Use fake fall date then 
        %	(won't use the resulting net_mod anyhow)
        else 
[net_mod,win_mod]=modelwhole(my_yr(i),zees(abcind),dateprev,datespringobsM2(abcind),datefallobsM2(abcind),prevfalMass(abcind),siteM2(abcind,:),glacier,lapse,meltrateS,meltrateI,precipratio);
        %model whole outputs are stratigraphic max and minimum
        
        end
        if missinb(misind(abcind),4)==0&&isnan(mbdb.data(j,5)) %started at min, just like regular year
            mbdb.data(j,5)=win_mod;
        end
        if missinb(misind(abcind),5)==0&&isnan(mbdb.data(j,6)) %started at min, just like regular year
            mbdb.data(j,6)=net_mod;
        end
    end
    if miss==1 %then looking for a balance, so should see which sites want to include
        % find sites will keep 
        ikeep0=1:length(useSites(:,1));%start with whole thing
        for h=1:maxNam
            ikeepa=ikeep0;
            ikeep0=find(useSites(ikeep0,h)==siteM2(abcind,h)); %check each index                                  %%%%%%%%%%%%%%%%Whoa
            ikeep=ikeepa(ikeep0);
        end
        if ~isempty(ikeep) %then site should be used in balance
            keepers(abcind,i)=1;%keep
        else
            keepers(abcind,i)=0;%don't keep
        end
    end
end 
% stop(here)
 if plot_ablation_model==1
        promptMessage = sprintf('Do you want to see the ablation model results for this year?');
        titleBarCaption = 'Plotting ablation model';
        choice = questdlg(promptMessage,titleBarCaption,'Yes','No','Don''t know')
                if strcmpi(choice, 'Yes')
                    plot_ablation_model1 = 1;
                elseif strcmpi(choice, 'No')
                    plot_ablation_model1=0;
                end
                    else
         plot_ablation_model1=0;
    end
[win_date_numM2, win_acc_adjM2,win_fxd_date_num2, win_fxd_date_adj2, net_date_numM2, net_abl_adjM2, annual_date_numM2, annual_abl_adjM2,search_min,search_max]=ablationmodel(my_yr(i),zees,datespringobsM2,datefallobsM2,mbdb.data(proc_index2,6),siteM2,glacier,lapse,meltrateS,meltrateI,precipratio,plot_ablation_model1);
%
% net_adj is always negative-- says should have started at this more
% negative place (the minimum)
for j=proc_index2(1):proc_index2(f2)
    abcind=j-proc_index2(1)+1;
    if isnan(mbdb.data(j,8))||isnan(sprMassOut(abcind))% then data won't be adjusted properly
        if net_date_numM2(abcind)<datefallobsM2(abcind) % went late
            sprMassOut(abcind)=0;
            mbdb.data(j,8)=abs(net_abl_adjM2(abcind));% extra mass positive
        else % went early or exact right time
            sprMassOut(abcind)=net_abl_adjM2(abcind);
            mbdb.data(j,8)=0; %exact right time

        end
    end
    if isnan(datefallobsM2(abcind)) % net at annual date, need to store datefallobs
      	datefallobsM2(abcind)=net_date_numM2(abcind);
    end
end
% correct net if have data from next year (f2 indices later); corrections
% are negative
net_abl_adjEnd=sprMassOut-mbdb.data(proc_index2,8);

% %zero zees that shouldn't be included in balanceend

if bal ==1 %we want only stratigraphic balances
    end_abl_adjM2=net_abl_adjM2;
    end_date2=net_date_numM2;
    endW_acc_adjM2= win_acc_adjM2;
    endW_date2=win_date_numM2;
elseif bal==2 %is from Oct 1 to Sept 30 == fixed date
    end_abl_adjM2=annual_abl_adjM2; %always positive
    end_date2=annual_date_numM2;
    endW_acc_adjM2= win_fxd_date_adj2;
    endW_date2=win_fxd_date_num2;
elseif bal==3%is from min to min
    zw = get_weights(bin_centers,gl_area,zees.*keepers(1:f2,i)); %just need yr for adjusting data, and no adjusts for this year because no data
    [glacier_net_date2,glacier_win_date2,min_to_glacmin2,max_to_glacmax2]=glacierWide(my_yr(i),win_date_numM2, search_max, net_date_numM2, search_min,glacier,zw);
    end_abl_adjM2= min_to_glacmin2; %always positive
    end_date2=glacier_net_date2;
    endW_acc_adjM2= max_to_glacmax2; %always negative
    endW_date2=glacier_win_date2;
end
             

% separate corrections 
% to take to balance to max (measurement+ end adjustment(positive))
% not adjusted yet because correction depends on total amount have till
% maximum (the total snow)
b_weqW(1:f2,i) = mbdb.data(proc_index2,5)+win_acc_adjM2;
stake_zs(1:f2,i)= mbdb.data(proc_index2,4);
% and corrections to take to annual or glacierwide, already adjusted
win_add(1:f2,i)=-beg_adjW+endW_acc_adjM2; %=0 if no glacierwide, NAN if annual
% % to take to balance to min (measurement+ end adjustment(negative))
% not adjusted yet because correction depends on total amount have till
% minimum (the total melt)
b_weq(1:f2,i) = mbdb.data(proc_index2,6)+net_abl_adjEnd;
% and corrections to take to annual or glacierwide, already adjusted        
net_add(1:f2,i)=-beg_adj+end_abl_adjM2; %=0 if no glacierwide
dateMin=end_date2;
dateMax=endW_date2; 
%
% zone_weights(1:f2,i) = zw;

dmin=nanmean(dateMin);%average if to site minimum, or all same date if to glacier/annual min
dateAMin(i)=round(dmin);
dmax=nanmean(dateMax);%average if to site maximum, or all same date if to glacier max
dateAMax(i)=round(dmax);

% reset params
siteM =siteM2;
siteM2=siteM3;
f =f2;
f2=f3;
glacier_net_date=dateAMin(i);
net_date_numM=net_date_numM2;
end_abl_adjM=end_abl_adjM2;

end 

%now add all together and zero out sites not keeping

b_weqW=(b_weqW+win_add).*keepers;
b_weqG =(b_weq +net_add); %Keep it all for use in the geodetic density correction
b_weqGW =(b_weqW +win_add); %Keep it all for use in the geodetic density correction
b_weq =(b_weq +net_add).*keepers;
b_weqS=b_weq-b_weqW;
stake_zs= stake_zs.*keepers;
[ b_barW,winter_gradients, b_barS,summer_gradients,annual_gradients,b_bar,ELA] = integrate_balance(b_weqW,b_weqS,b_weq,AADs,my_yr,stake_zs,grad,surf,plot_integration);
end


%stop(here)