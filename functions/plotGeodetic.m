function [errb,geod2,Badj,SEC2]=plotGeodetic(glacier,bal,surf,my_yr,option,grad,b_weqG, b_weqGW,bal_name,folder,surf3,plot_ablation_model)
% Function calculates errorbars on cummulative
% inputs
%   glacier - 0 1 or 2 for the 2 glaciers and the test data
%   bal - depending on which balance you want to plot, annual or fixed-date
%   my_yr - yrs want to plot
%   option - compute geod2==1, don't==0
%   grad - use index method==1, use gradient method==2
%
% outputs
%   errb - error bar vector for cummulative, error(yr1+yr2)=swrt(error(yr1)^2+error(yr2)^2)
%   geod2 - [yr,geodetic measurement,error bar]
%   delH - average height change on glacier, or mass balance, from Bahr, or
%     Arendt, or Luthi, scaling relationship 
%
dbstop if error
%% Input file names
        mygeod0 = ['data/',glacier,'/Input/Input_',glacier,'_Geodetics.csv']; %wolverine dem differences, reevaluated by L Sass Spring 2012
        mysec = ['data/',glacier,'/Input/Input_',glacier,'_Altimetry.csv']; %secondary geodetic check, designed to run on laser altimetry but could be from other sources 
        err0=0.2; % m/area per year for index method, estimated by Rod March (1998, a gulkana data report)
        mydata = ['data/',glacier,'/Input/Input_',glacier,'_Glaciological_Data.csv'];
        myAAD = ['data/',glacier,'/Input/Input_',glacier,'_Area_Altitude_Distribution.csv'];
        myprecipratio = ['data/',glacier,'/Input/Calibrated_',glacier,'_Precipitation_Ratios.csv'];
        mymeltrate = ['data/',glacier,'/Input/Calibrated_',glacier,'_Degree_Day_Factors.csv'];
%% Output file names
        mybal = ['data/',glacier,'/Output/',folder,'/Output_',glacier,'_',surf3,'_',bal_name,'_oldmb.csv'];
        mygeodetic=['data/',glacier,'/Output/',folder,'/Output_',glacier,'_',surf3,'_',bal_name,'_Geodetic_Calibration.csv']; % geodetic measurements stored here
        myaltimetry=['data/',glacier,'/Output/',folder,'/Output_',glacier,'_',surf3,'_',bal_name,'_Altimetry_Comparison.csv']; % Altimetry measurements stored here
        myadjseasonalbal=['data/',glacier,'/Output/',folder,'/Output_',glacier,'_',surf3,'_',bal_name,'_Calibrated_Seasonal_Balance.csv']; % corrected balance time series stored here
        myadjannualbal=['data/',glacier,'/Output/',folder,'/Output_',glacier,'_',surf3,'_',bal_name,'_Calibrated_Annual_Balance.csv']; % corrected balance time series stored here
        lapse=-6.5; %moist adiabatic lapse rate 

geod0 = importdata(mygeod0);%,'NumHeaderLines', 3, 'NumColumns', 5, 'Format', '%d %d %d %f %f');
geod0=geod0.data;
 geod0(:,4) = geod0(:,4) * 850/1000;% -- WATER EQUIVALENT UNITS ASSUMING ALL CHANGES ARE IN ICE VOLUME!!! (wtf?)
  % This should get cleaned up, it makes later corrections rather confusing. 
secondary = importdata(mysec);
secondary=secondary.data;
if isempty(secondary)
else
secondary(:,4) =  secondary(:,4) * 850/1000;%
end
%
my_yr=[my_yr(1)-1,my_yr];%before first year have a 0 cummulative
numYrs=length(my_yr);
ind=[];
for i=1:numYrs
    yes=find(geod0(:,1)==my_yr(i));
    if ~isempty(yes)
        ind=[ind,yes]; %#ok<AGROW>
    end
    if length(ind)>=2 %need at least 2 photos to have a balance check
        %geod0(ind(1),5)=0; %zero error bar on first year !!!!!!!!!!!!!!!!
        geod1=geod0(ind,:);
    else
        geod1=[]; 
    end
end

%Repeat same stuff to correct secondary check !!!!!!!!!!!!!!!!!!!!!!!!!!!!
%added by LS 2013.10
SECind=[];
for i=1:numYrs
    Uyes=find(secondary(:,1)==my_yr(i));
    if ~isempty(Uyes)
        SECind=[SECind,Uyes]; %#ok<AGROW>
    end
    if length(SECind)>=2 %need at least 2 photos to have a balance check
        secondary(SECind(1),5)=0; %zero error bar on first year
        SEC1=secondary(SECind,:);
        numSEC=length(SEC1(:,1));
    else
        SEC1=[]; 
        numSEC=0;
    end
end

%Load background data to build adjustments

if ~isempty(geod1) && option==1 % don't bother if doesn't exist
    numGeod=length(geod1(:,1));
    
    Bts=txt2mat(mybal);
    precipratio=load(myprecipratio); 
    meltrate=load(mymeltrate);
    meltrateS=meltrate(1); 
    meltrateI=meltrate(2);
    %
    % load and format data-- need altitudes and AAD
    % Original data: bal_yr stake visit_1 visit_2 z winter net spring_mass_outside fall_mass_outside
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
    siteAll=siteAll(:,1:maxNam);
    AADs = importdata(myAAD);
    if isnan(AADs(1,2)) %some excel versions screw this up and see blanks as a line
        AADs=AADs(2:end,:);
    end
    bin_centers = AADs(1,2:end);
    maxStak=0;
    proc_index00=NaN*ones(300,numGeod); %if ever had more than 300 stakes would have an error
    prlen=NaN*ones(1,numGeod);
    for i=1:numGeod % find max stakes ever have
        proc_index0 = find(bal_yr == geod1(i,1)); %index of data avail for next balance year
        stak=length(proc_index0);
        if stak>maxStak
            maxStak=stak;
        end
        proc_index00(1:stak,i)=proc_index0;
        prlen(i)=stak;
    end
    obsminadjd =zeros(maxStak,numGeod);
    endminadjd =zeros(maxStak,numGeod);
    obsmaxadjd =zeros(maxStak,numGeod);
    startmaxadjd =zeros(maxStak,numGeod);
    stake_zs=nan*zeros(maxStak,numGeod);
    zone_weights=zeros(maxStak,numGeod);
    spring_max_date =zeros(maxStak,numGeod);
    fall_min_date =zeros(maxStak,numGeod);
    gl_areaTot=zeros(numGeod,1);
     
    % calculate mass outside for Geodetics to be at minimum year end
    %
    %stakeB = [];
    for i=1:numGeod
        index(i)=find(Bts(:,1)==geod1(i,1));
        if geod1(i,1)~=my_yr(1)% if at 0 cummulative yr, skip
            % surf is off current geometry == conventional
            AAD_index = find(AADs(:,1) == geod1(i,1)); 
            gl_area = AADs(AAD_index,2:end); 
            % find ablation after date till end of year
            proc_index = proc_index00(1:prlen(i),i); %index of data avail for balance year
            zees = mbdb.data(proc_index,1); %index site altitudes needed for lapses
            stake_zs(1:length(zees),i) = zees;
            f=prlen(i);
            siteM=siteAll(proc_index,:);
            photo_date = geod1(i,1:3);
            datefallobsM=datenum(photo_date)*ones(f,1);
            datespringobsM=datenum(photo_date)*ones(f,1);
            datenet=NaN*ones(f,1); % amount up to date-- don't know
             if plot_ablation_model==1
        promptMessage = sprintf('Do you want to see the ablation model results for adjusting geodetics?');
        titleBarCaption = 'Plotting ablation model';
        choice = questdlg(promptMessage,titleBarCaption,'Yes','No','Don''t know');
                if strcmpi(choice, 'Yes')
                    plot_ablation_model1 = 1;
                elseif strcmpi(choice, 'No')
                    plot_ablation_model1=0;
                end
                    else
         plot_ablation_model1=0;
    end
[win_date_numM, win_acc_adjM,  win_fxd_date_num, win_fxd_date_adj,net_date_numM, net_abl_adjM, annual_date_numM, annual_abl_adjM,search_min,search_max]=ablationmodel(geod1(i,1),zees,datespringobsM,datefallobsM,datenet,siteM,glacier,lapse,meltrateS,meltrateI,precipratio,plot_ablation_model1);


           switch bal 
                case 1 %we want only stratigraphic balances
                    end_abl_adjM = net_abl_adjM; 
                    begin_abl_adjM = win_acc_adjM; 
                    bal2='Stratigraphic Balance'; %for plotting
                    dat='min_dat = BAL';
                    end_date=round(nanmean(net_date_numM));
                    spring_max_date(1:length(win_date_numM),i) = round(nanmean(win_date_numM)); %mass max dates for years of flights
                    fall_min_date(1:length(net_date_numM),i) = round(nanmean(net_date_numM)); % mass min dates for years of flights
                case 2 %Fixed date balance Hydro year
                    end_abl_adjM = annual_abl_adjM;
                    begin_abl_adjM = win_fxd_date_adj;
                    bal2='Fixed-Date Balance'; %for plotting
                    dat='min_dat = BAL';
                    end_date=round(nanmean(annual_abl_adjM));
                    spring_max_date(1:length(win_fxd_date_num),i) = round(nanmean(win_fxd_date_num)); %mass max dates for years of flights
                    fall_min_date(1:length(annual_date_numM),i) = round(nanmean(annual_date_numM)); % mass min dates for years of flights
                case 3 %FIxed-date system 
                    zone_weights(1:f,i) = get_weights(bin_centers,gl_area,zees);
                    [glacier_net_date,glacier_win_date,min_to_glacmin,max_to_glacmax]=glacierWide(geod1(i,1),win_date_numM, search_max, net_date_numM, search_min,glacier,zone_weights(1:f,i));
                    end_abl_adjM = min_to_glacmin;
                    begin_abl_adjM = max_to_glacmax;
                    bal2='Combine-Date System'; %for plotting
                    dat='min_dat = BAL';
                    end_date=glacier_net_date;
                    spring_max_date(1:length(glacier_win_date),i) = glacier_win_date; %mass max dates for years of flights
                    fall_min_date(1:length(glacier_win_date),i) = glacier_net_date; % mass min dates for years of flights
            end 
            
        else % geod(i,1)~=my_yr(1) at 0 cummulative yr, skip
            geod1(i,4)=0;
            geod1(i,5)=0;
        end
        gl_areaTot(i)=sum(gl_area);
        stakeB(:,i)=b_weqG(:,index(i));
        S_stakeB(:,i)=b_weqG(:,index(i)) - b_weqGW(:,index(i));
%         spring_max_date(1:length(glacier_win_date),i) = glacier_win_date;
%         fall_min_date(1:length(glacier_win_date),i) = glacier_net_date;
        for a = 1:length(spring_max_date(:,1));
            photoday(a,i) = datenum(photo_date);
        end
        
    end
    
    b_bar_model_max = (obsmaxadjd+startmaxadjd);%MMM - model max module - for adjusting the output from ablationmodel.m to incorporate both spring and fall measurements
    [m,n] = size(b_bar_model_max);
    b_bar_max_sub = S_stakeB(1:m,1:n) + b_bar_model_max;
    
    L_summer = fall_min_date-spring_max_date;
    R_after = (fall_min_date-photoday)./L_summer;
    ind_inf = find(isinf(R_after)==1);
    R_after(ind_inf) = 0;
    ind_early = find(R_after>=1);
    R_after(ind_early) = 1;
    ind_late = find(R_after<=0);
    R_after(ind_late) = 0;
    R_before = ones(size(R_after))-R_after;
   
    
    b_baradj_p_original = (obsminadjd+endminadjd);
    

    b_baradj_p = b_bar_max_sub.*R_after + b_baradj_p_original.*R_before; %hybrid, uses proportions of each model based on flight date
    
    % from scaling, delta_volume=ca*( area^gamma-previous_area^gamma )
    gamma=1.375; %from Arendt et al 2006
    ca=0.28; %from Arendt 2006 his in m^(3-2*gamma)
    delH=ca*( (gl_areaTot(1:end-1,:).*1e6).^gamma - (gl_areaTot(end,:).*1e6).^gamma)./(gl_areaTot(2:end).*1e6);
    delH=[delH;0]; 
    b_baradj_p=(obsminadjd+endminadjd);
    [m,n]=size(b_baradj_p);
    for i = 1:m;
        for j = 1:n;
            if b_baradj_p(i,j)>=0; %new snow on surface at time of photography
                b_baradj_p2(i,j)=b_baradj_p(i,j)/0.2; %probably light density
            elseif b_baradj_p(i,j)<=0 && stakeB(i,j)>=0; % snow melt
                b_baradj_p2(i,j)=b_baradj_p(i,j)/0.59; %We already accounted for the density difference between ice and water, so 0.53 assumes 450 kg.m^3, 0.59 assumes 500 
            elseif b_baradj_p(i,j)<=0 && stakeB(i,j)<=0 && b_baradj_p(i,j)<=stakeB(i,j); % snow and ice melt
                density = (0.59 + (0.47*(stakeB(i,j)/b_baradj_p(i,j))));
                b_baradj_p2(i,j)=b_baradj_p(i,j)/density; %density starts at 500 and approaches 850 with increasing proportions of icemelt
            else
                b_baradj_p2(i,j)=b_baradj_p(i,j); %all ice, so leave at 850
            end
        end
    end
      [ ~,~, ~,~,~,b_baradj,~] = integrate_balance(NaN,NaN,b_baradj_p2,AADs,geod1(:,1)',stake_zs,grad,surf,0);


    
    geod2a=[geod1(:,1),geod1(:,4)+b_baradj.',geod1(:,5)];
    geod2=[geod2a(:,1),geod2a(:,2)-geod2a(end,2)*ones(numGeod,1),geod2a(:,3)]; %normalize to first year
    year1=geod1(1,1); %year to zero 
    result=[geod2,delH,b_baradj.'].';
fprintf(1,'%s, FOR CONVENTIONAL BALANCE\n',bal2);
fprintf(1,'GEODETIC BALANCES ADJUSTED TO THE MASS MINIMUM\n');
fprintf(1,'ZEROED TO LAST AERIAL PHOTO ENDING AT %s YEAR END\n',dat);
fprintf(1,'ALSO INCLUDE VOL PREDICTION FROM VArea-SCALING LAW\n');
fprintf(1,'year geo_cumul_bal error_bar V_scal_pred ablation_adj(included)\n');
fprintf(1,'%4d %10.4f    +-%6.4f %10.4f %12.4f\n', result);
file = fopen(mygeodetic,'w');
fprintf(file,'%s FOR CONVENTIONAL BALANCE\r\n',bal2);
fprintf(file,'GEODETIC BALANCES ADJUSTED TO THE MASS MINIMUM\r\n');
fprintf(file,'ZEROED TO LAST AERIAL PHOTO ENDING AT %s YEAR END\r\n',dat);
fprintf(file,'ALSO INCLUDE VOL PREDICTION FROM VArea-SCALING LAW\r\n');
fprintf(file,'year geo_cumul_bal error_bar V_scal_pred ablation_adj(included)\r\n');
fprintf(file,'%4d %10.4f    +-%6.4f %10.4f %12.4f\r\n', result);

fclose(file); 

G = result'; 
Bmax = Bts(:,6);%(datenum([Bts(:,1),Bts(:,6:7)]) - datenum([Bts(:,1),ones(length(Bts),2)]))/365 + Bts(:,1);
Bmin = Bts(:,7);%(datenum([Bts(:,1),Bts(:,9:10)]) - datenum([Bts(:,1),ones(length(Bts),2)]))/365 + Bts(:,1);
FirstMin = (Bts(1,1)-1) + (datenum([(Bts(1,1)-1),10,1])-datenum([(Bts(1,1)-1),1,1]))/365; %assumed minimum for start of first winter balance
Ttot = Bmin - (FirstMin); %time elapsed since october 1 prev yr

B = [Bts(:,1:5), Bmax, Bmin, Ttot];
[m,n] = size(G);
for i=1:m;
    ind(i) = find(B(:,1)==G(i,1));
    V(i) = B(ind(i),5)-B(ind(1),5); %volume change from the balance time series in the row with the appropriate year
    R(i) = V(i)-G(i,2); %residual
    T(i) = B(ind(i),7)-B(ind(1),7);%Time
    W(i) = 1/(G(i,3)+1e-10);
end

mdl = LinearModel.fit(T, R, 'linear','weights', W);
adj = mdl.Coefficients{2,1};

[m,n] = size(B);
Badj = [B(:,1),nan(length(B),4)];
Badj(1,2) = B(1,2)-((B(1,6) - FirstMin) * adj); %winter
Badj(1,3) = B(1,3)-((B(1,7) - B(1,6)) * adj); %summer
Badj(1,4) = B(1,4)-((B(1,7) - FirstMin) * adj); %annual
Badj(1,5) = B(1,5)-((B(1,7) - FirstMin) * adj); %cumulative

B(1,3)-B(1,8)*adj;
for i=2:m;
    Badj(i,2) = B(i,2)-((B(i,6) - B(i-1,7)) * adj); %winter
    Badj(i,3) = B(i,3)-((B(i,7) - B(i,6)) * adj); %summer
    Badj(i,4) = B(i,4)-((B(i,7) - B(i-1,7)) * adj); %annual
    Badj(i,5) = B(i,5)-((B(i,7) - FirstMin) * adj); %cumulative
    B(i,9) = ((.2^2)*i)^(1/2);
end

Badj = [Badj,Bmax,Bmin];
temp = [Bmax;Bmin];
y = fix(temp);
%%
Date = datestr(datenum(y,0,(temp - y).*365),'yyyy/mm/dd'); % corrected
SeasonalBalance_m_we = round([Badj(:,2);Badj(:,3)],1);
AnnualBalances_m_we(:,1) = round([Badj(:,2)+Badj(:,3)],1);
CumulativeBalance_m_we = round([Badj(:,5)-Badj(:,3);Badj(:,5)],1);
mb_temp = table(Date,SeasonalBalance_m_we,CumulativeBalance_m_we);
mb_out = sortrows(mb_temp);
bw_index=find(table2array(mb_out(:,2))>0);
bs_index = find(table2array(mb_out(:,2))<0);
[~,~,bw_rank] = unique(table2array(mb_out(bw_index,2)));
[~,~,bs_rank] = unique(table2array(mb_out(bs_index,2)));

mb_out(bs_index,4)=array2table(nan,'VariableNames',{'Bw_Rank'});
mb_out(bw_index,5)=array2table(nan,'VariableNames',{'Bs_Rank'});
mb_out(bw_index,4)=array2table((max(bw_rank)-bw_rank)+1,'VariableNames',{'Bw_Rank'});
mb_out(bs_index,5)=array2table((bs_rank),'VariableNames',{'Bs_Rank'});
mb_out.Properties.VariableNames = {'Date','Seasonal_Balance_m_we','Cumulative_Balance_m_we','Bw_Rank','Bs_Rank'};
%stop(here)
writetable(mb_out,myadjseasonalbal)
for i=1:length(bw_index)
Ba_amplitude(i,1)=(table2array(mb_out(bw_index(i),2))-table2array(mb_out(bs_index(i),2)))/2;
end

[~,~,ba_rank] = unique(AnnualBalances_m_we);
[~,~,amp_rank] = unique(Ba_amplitude(:,1));
annual_mb_out(:,1:6)=[mb_out(2:2:end,1) table(AnnualBalances_m_we) mb_out(2:2:end,3) table(Ba_amplitude,'VariableNames',{'Amplitude_m_we'}) table((max(ba_rank)-ba_rank)+1,'VariableNames',{'Ba_Rank'}) table((max(amp_rank)-amp_rank)+1,'VariableNames',{'Amplitude_Rank'})];
annual_mb_out.Properties.VariableNames = {'Date', 'Annual_Balance_m_we', 'Cumulative_Balance_m_we','Amplitude', 'Ba_Rank','Amplitude_Rank'}
% stop(here)
writetable(annual_mb_out,myadjannualbal)

fprintf(1,'WEIGHTED LEAST SQUARES ADJUSTMENT TO THE BALANCE TIME SERIES\r\n');
%fprintf(1,'year geo_cumul_bal error_bar V_scal_pred ablation_adj(included)\n');
fprintf(1,'%6.4f m a^(-1)\r\n', -adj);
fprintf(1,'REFER TO calibrated_output.txt FOR ADJUSTED ANNUAL AND CUMULATIVE MASS BALANCE VALUES\r\n');

clear L_summer R_after ind_inf ind_early ind_late b_bar_model_max b_bar_max_sub m n i j 

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% start secondary geodetic module -- basically do all of the same things to the
% volume change estimates from secondary geodetic checks that we did to the volume 
% change estimates from photogrammetric DEMs, except we don't adjust the 
% time series based on those outputs. Used for UAF altimetry
if numSEC==0 || isempty(SEC1)
    SEC2=[];
else
for i=1:numSEC % find max stakes ever have
        proc_index0SEC = find(bal_yr == SEC1(i,1)); %index of data avail for next balance year
        stak=length(proc_index0SEC);
        if stak>maxStak
            maxStak=stak;
        end
        proc_index00SEC(1:stak,i)=proc_index0SEC;
        prlenSEC(i)=stak;
end
    obsminadjd =zeros(maxStak,numSEC); %initialize matrices and overwrite old values
    endminadjd =zeros(maxStak,numSEC);
    obsmaxadjd =zeros(maxStak,numSEC);
    startmaxadjd =zeros(maxStak,numSEC);
    stake_zs=nan*zeros(maxStak,numSEC);
    zone_weights =zeros(maxStak,numSEC);
    spring_max_date =zeros(maxStak,numSEC);
    fall_min_date =zeros(maxStak,numSEC);

    gl_areaTotSEC=zeros(numSEC,1);
    % calculate mass outside for Geodetics to be at minimum year end
    %
    %stakeB = [];
    for i=1:numSEC
        indexSEC(i)=find(Bts(:,1)==SEC1(i,1));
        if SEC1(i,1)~=my_yr(1)% if at 0 cummulative yr, skip
            % surf is off current geometry == conventional
            AAD_indexSEC = find(AADs(:,1) == SEC1(i,1)); 
            gl_areaSEC = AADs(AAD_indexSEC,2:end); 
            % find ablation after date till end of year
            proc_indexSEC = proc_index00SEC(1:prlenSEC(i),i); %index of data avail for balance year
            zees = mbdb.data(proc_indexSEC,1); %index site altitudes needed for lapses
            stake_zs(1:length(zees),i) = zees;
            f=prlenSEC(i);
            siteM=siteAll(proc_indexSEC,:);
            flight_date = SEC1(i,1:3);
            datefallobsM=datenum(flight_date)*ones(f,1);
            datespringobsM=datenum(flight_date)*ones(f,1);
            datenet=NaN*ones(f,1); % amount up to date-- don't know
             if plot_ablation_model==1
        promptMessage = sprintf('Do you want to see the ablation model results for adjusting the altimetry data?');
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
             
[win_date_numM, win_acc_adjM, win_fxd_date_num, win_fxd_date_adj, net_date_numM, net_abl_adjM, annual_date_numM, annual_abl_adjM,search_min,search_max]=ablationmodel(SEC1(i,1),zees,datespringobsM,datefallobsM,datenet,siteM,glacier,lapse,meltrateS,meltrateI,precipratio,plot_ablation_model1);
           switch bal 
                case 1 %we want only stratigraphic balances
                    end_abl_adjM = net_abl_adjM; 
                    begin_abl_adjM = win_acc_adjM; 
                    bal2='Stratigraphic Balance'; %for plotting
                    dat='min_dat = BAL';
                    end_date=round(nanmean(net_date_numM));
                    spring_max_date(1:length(win_date_numM),i) = round(nanmean(win_date_numM)); %mass max dates for years of flights
                    fall_min_date(1:length(net_date_numM),i) = round(nanmean(net_date_numM)); % mass min dates for years of flights
                case 2 %Fixed date balance Hydro year
                    end_abl_adjM = annual_abl_adjM;
                    begin_abl_adjM = win_fxd_date_adj;
                    bal2='Fixed-Date Balance'; %for plotting
                    dat='min_dat = BAL';
                    end_date=round(nanmean(annual_abl_adjM));
                    spring_max_date(1:length(win_fxd_date_num),i) = round(nanmean(win_fxd_date_num)); %mass max dates for years of flights
                    fall_min_date(1:length(annual_date_numM),i) = round(nanmean(annual_date_numM)); % mass min dates for years of flights
                case 3 %FIxed-date system 
                    zone_weights(1:f,i) = get_weights(bin_centers,gl_areaSEC,zees);
                    [glacier_net_date,glacier_win_date,min_to_glacmin,max_to_glacmax]=glacierWide(SEC1(i,1),win_date_numM, search_max, net_date_numM, search_min,glacier,zone_weights(1:f,i));
                    end_abl_adjM = min_to_glacmin;
                    begin_abl_adjM = max_to_glacmax;
                    bal2='Combine-Date System'; %for plotting
                    dat='min_dat = BAL';
                    end_date=glacier_net_date;
                    spring_max_date(1:length(glacier_win_date),i) = glacier_win_date; %mass max dates for years of flights
                    fall_min_date(1:length(glacier_win_date),i) = glacier_net_date; % mass min dates for years of flights
            end 

           
        else % geod(i,1)~=my_yr(1) at 0 cummulative yr, skip
            SEC1(i,4)=0;
            SEC1(i,5)=0;
        end
        gl_areaTotSEC(i)=sum(gl_areaSEC); %glacier area for each flight
        stakeBSEC(:,i)=b_weqG(:,indexSEC(i)); % pull out annual stake balances for flight years
        S_stakeBSEC(:,i)=b_weqG(:,indexSEC(i))-b_weqGW(:,indexSEC(i));% pull out stake summer balances for flight years
%         spring_max_date(1:length(glacier_win_date),i) = glacier_win_date; %mass max dates for years of flights
%         fall_min_date(1:length(glacier_win_date),i) = glacier_net_date; % mass min dates for years of flights
        for a = 1:length(spring_max_date);
            flyday(a,i) = datenum(flight_date); %get date number of each flight date
        end
    end
    
    b_bar_model_max = (obsmaxadjd+startmaxadjd);%MMM - model max module - for adjusting the output from ablationmodel.m to incorporate both spring and fall measurements
    [m,n] = size(b_bar_model_max);
    b_bar_max_sub = S_stakeBSEC(1:m,1:n) + b_bar_model_max;
    
    L_summer = fall_min_date-spring_max_date; %determine the length of ablation season
    R_after = (fall_min_date-flyday)./L_summer; % get ratio of flight date to mass minimum 
    infind = find(isinf(R_after)==1);
    R_after(infind) = 0;
    ind_early = find(R_after>=1);
    R_after(ind_early) = 1;
    ind_late = find(R_after<=0);
    R_after(ind_late) = 0;
    R_before = ones(size(R_after))-R_after;
   
    
    b_baradj_p_original = (obsminadjd+endminadjd);
    
    b_baradj_p = b_bar_max_sub.*R_after + b_baradj_p_original.*R_before; %hybrid, uses proportions of each model based on flight date
    
    [m,n]=size(b_baradj_p);
    b_baradj_p2 = [];
   
    for i = 1:m;
        for j = 1:n;
           
            if b_baradj_p(i,j)>0; %new snow on surface at time of photography
                b_baradj_p2(i,j)=b_baradj_p(i,j)/0.35; %probably light density
                %correction = 'new_snow'
            elseif b_baradj_p(i,j)<=0 && stakeBSEC(i,j)>=0; % snow melt
                b_baradj_p2(i,j)=b_baradj_p(i,j)/0.53; %high density
                %correction = 'snow'
            elseif b_baradj_p(i,j)<=0 && stakeBSEC(i,j)<=0 && b_baradj_p(i,j)<=stakeBSEC(i,j); % snow and ice melt
                density = (0.53 + (0.47*(stakeBSEC(i,j)/b_baradj_p(i,j))));
                b_baradj_p2(i,j)=b_baradj_p(i,j)/density; %high density
                %correction = 'snowandice'
            else
                b_baradj_p2(i,j)=b_baradj_p(i,j); %all ice, noa adjustment
                %correction = 'ice'
            end
            
        end
    end
    [ ~,~, ~,~,~,b_baradj,~] = integrate_balance(NaN,NaN,b_baradj_p2,AADs,SEC1(:,1)',stake_zs,grad,surf,0);
    
    
    SEC2a=[SEC1(:,1),SEC1(:,4)+b_baradj.',SEC1(:,5)];
    SEC2=[SEC2a(:,1),SEC2a(:,2)-SEC2a(1,2)*ones(numSEC,1),SEC2a(:,3)]; %normalize to first year
    SECyear1=SEC1(1,1); %year to zero 
    resultSEC=[SEC2,b_baradj.'].';
%fprintf(1,'%s, FOR CONVENTIONAL BALANCE\r\n',bal2);
fprintf(1,'UAF ALTIMETRY BALANCES ADJUSTED TO THE MASS MINIMUM\n');
fprintf(1,'ZEROED TO FIRST FLIGHT DATE ENDING AT %s YEAR END\n',dat);
fprintf(1,'year geo_cumul_bal error_bar ablation_adj(included)\n');
fprintf(1,'%4d %10.4f    +-%6.4f %12.4f\n', resultSEC);
file = fopen(myaltimetry,'w');
%fprintf(file,'%s FOR CONVENTIONAL BALANCE\r\n',bal2);
fprintf(file,'UAF ALTIMETRY BALANCES ADJUSTED TO THE MASS MINIMUM\r\n');
fprintf(file,'ZEROED TO FIRST FLIGHT ENDING AT %s YEAR END\r\n',dat);
fprintf(file,'year geo_cumul_bal error_bar ablation_adj(included)\r\n');
fprintf(file,'%4d %10.4f    +-%6.4f %12.4f\r\n', resultSEC);

fclose(file); 
end      
else
    geod2=[];
    year1=my_yr(1)-1; %year want to zero to-- so has err0 for first year
fprintf(1,'GEODETIC BALANCES NOT RELEVANT BECAUSE OF PROCESSING YEARS \r\n');
fprintf(1,'SELECTION OR BECAUSE NOT PLOTTING CONVENTIONAL BALANCE.\r\n');  

end

errb=(err0^2*abs(my_yr(2:end).'-year1*ones(numYrs-1,1))).^0.5;%error propagation over years for cummulative
end
    