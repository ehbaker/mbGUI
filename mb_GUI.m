function varargout = mb_GUI(varargin)
dbstop if error
 
addpath functions/
% MB_GUI M-file for mb_GUI.fig
%      MB_GUI, by itself, creates a new MB_GUI or raises the existing
%      singleton*.
%
%      H = MB_GUI returns the handle to a new MB_GUI or the handle to
%      the existing singleton*.
%
%      MB_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MB_GUI.M with the given input arguments.
%
%      MB_GUI('Property','Value',...) creates a new MB_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mb_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mb_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mb_GUI

% Last Modified by GUIDE v2.5 16-Feb-2017 20:04:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mb_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @mb_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
% 

%%%%%%%%%%% Setup of GUI for User  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mb_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for mb_GUI
addpath functions/ %addpath to functions
handles.output = hObject;
mydir = pwd; % get working directory
% set(handles.MYDIR_EDIT,'String',mydir)
handles.params.mydir = mydir;
% Choose default command line output for setup_inputs
handles.output = hObject;


%%%%%% setup site map %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Gulkana by default %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(handles.Site_Map);
    glacier_map=['data/Gulkana/Input/Input_Gulkana_Glaciological_Sites.jpg']; 
    locafig = imread(glacier_map);
        image(locafig);
        axis off          % Remove axis ticks and numbers
        axis image 
        
        glacier = 'Gulkana';
        mydata = ['data/',glacier,'/Input/Input_',glacier,'_Glaciological_Data.csv'];
        data=importdata(mydata);
        all_years=cell2mat(data.textdata(2:end,1));
        unique_years=unique(str2num(all_years));
        first_yr=num2str(unique_years(1,1));
        last_yr=num2str(unique_years(end-1,1));
        my_yr=[first_yr,':',last_yr];
        
        set(handles.YEAR_EDIT,'String',my_yr)
        handles.my_yr = my_yr;
use=zeros(14,1);%return all names
    [mysit] =use_sites(glacier,use);
    set(handles.SITE1,'String',mysit(1,:));
    set(handles.SITE2,'string',mysit(2,:));
    set(handles.SITE3,'string',mysit(3,:));
    set(handles.SITE4,'string',mysit(4,:));
    set(handles.SITE5,'string',mysit(5,:));
    set(handles.SITE6,'string',mysit(6,:));
    set(handles.SITE7,'string',mysit(7,:));
    set(handles.SITE8,'string',mysit(8,:));
    set(handles.SITE9,'string',mysit(9,:));
    set(handles.SITE10,'string',mysit(10,:));
    set(handles.SITE11,'string',mysit(11,:));
    set(handles.SITE12,'string',mysit(12,:));
    set(handles.SITE13,'string',mysit(13,:));
    set(handles.SITE14,'string',mysit(14,:));
    
set(handles.SURF_POPUP,'Value',1);
set(handles.BAL_POPUP,'Value',1);
% Update handles structure
guidata(hObject, handles);
%
function varargout = mb_GUI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function GRAD_POPUP_Callback(hObject, eventdata, handles)
%
function GRAD_POPUP_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in GLACIER_ID.
function GLACIER_ID_Callback(hObject, eventdata, handles)
% hObject    handle to GLACIER_ID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
glacier = get(handles.Call_Glacier,'String');
mydata = ['data/',glacier,'/Input/Input_',glacier,'_Glaciological_Data.csv'];
        data=importdata(mydata);
        all_years=cell2mat(data.textdata(2:end,1));
        unique_years=unique(str2num(all_years));
        first_yr=num2str(unique_years(1,1));
        last_yr=num2str(unique_years(end-1,1));
        my_yr=[first_yr,':',last_yr];
        
        set(handles.YEAR_EDIT,'String',my_yr)
        handles.my_yr = my_yr;
use=zeros(14,1);%return all names
    [mysit] =use_sites(glacier,use);
    set(handles.SITE1,'String',mysit(1,:));
    set(handles.SITE2,'string',mysit(2,:));
    set(handles.SITE3,'string',mysit(3,:));
    set(handles.SITE4,'string',mysit(4,:));
    set(handles.SITE5,'string',mysit(5,:));
    set(handles.SITE6,'string',mysit(6,:));
    set(handles.SITE7,'string',mysit(7,:));
    set(handles.SITE8,'string',mysit(8,:));
    set(handles.SITE9,'string',mysit(9,:));
    set(handles.SITE10,'string',mysit(10,:));
    set(handles.SITE11,'string',mysit(11,:));
    set(handles.SITE12,'string',mysit(12,:));
    set(handles.SITE13,'string',mysit(13,:));
    set(handles.SITE14,'string',mysit(14,:));
    
    axes(handles.Site_Map);
    glacier_map=['data/',glacier,'/Input/Input_',glacier,'_Glaciological_Sites.jpg']; 
    locafig = imread(glacier_map);
        image(locafig);
        axis off          % Remove axis ticks and numbers
        axis image        % Set aspect ratio to obtain square pixels
       
       
    
    guidata(hObject,handles)
    
function Call_Glacier_Callback(hObject, eventdata, handles)
% hObject    handle to Call_Glacier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Call_Glacier as text
%        str2double(get(hObject,'String')) returns contents of Call_Glacier as a double


% --- Executes during object creation, after setting all properties.
function Call_Glacier_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Call_Glacier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%
function YEAR_EDIT_Callback(hObject, eventdata, handles)
%
function YEAR_EDIT_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% 
% 
%function will write a file that says which sites will use to compute
%balances
%
%
function SITE1_Callback(hObject, eventdata, handles)
function SITE2_Callback(hObject, eventdata, handles)
function SITE3_Callback(hObject, eventdata, handles)
function SITE4_Callback(hObject, eventdata, handles)
function SITE5_Callback(hObject, eventdata, handles)
function SITE6_Callback(hObject, eventdata, handles)
function SITE7_Callback(hObject, eventdata, handles)
function SITE8_Callback(hObject, eventdata, handles)
function SITE9_Callback(hObject, eventdata, handles)
function SITE10_Callback(hObject, eventdata, handles)
function SITE11_Callback(hObject, eventdata, handles)
function SITE12_Callback(hObject, eventdata, handles)
function SITE13_Callback(hObject, eventdata, handles)
function SITE14_Callback(hObject, eventdata, handles)
function ALLSITES_Callback(hObject, eventdata, handles)
%
%
%
function SURF_POPUP_Callback(hObject, eventdata, handles)
%
function SURF_POPUP_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% 
function BAL_POPUP_Callback(hObject, eventdata, handles)
%
function BAL_POPUP_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
%
function ERRORBAR_CHECK_Callback(hObject, eventdata, handles)
%
%

function UPDATE_BUTTON_CreateFcn(hObject, eventdata, handles)
%
function UPDATE_POPUP_Callback(hObject, eventdata, handles)
%
function UPDATE_POPUP_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%       UPDATE INPUTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   Afer user has entered a glacier name  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function UPDATE_BUTTON_Callback(hObject, eventdata, handles)
% Get user input from GUI
mydir = pwd;                                            %gets present working directory
glacier = str2mat(get(handles.Call_Glacier,'String'))   %get glacier entered
 grad = get(handles.GRAD_POPUP,'Value');                %get gradients you are using
 
%

    
locaS=['data/',glacier,'/Input/Input_',glacier,'_Glaciological_Sites.jpg']; 
   

%
cd(mydir)
addpath functions

%
switch get(handles.UPDATE_POPUP,'Value') 
    case 1 %correct temp and precip data
        ready=fillWxData(glacier);
        if ready==1
            set(handles.UPDATE_BUTTON,'string','Can now invert for params');
        else
            set(handles.UPDATE_BUTTON,'string','Weather data has NaNs');
        end
        cd(mydir)
        guidata(hObject, handles);
    case 2 %calculate precip/catch ratios
        [ready,usebeta]= calcStakeCatch(glacier);
        if ready==length(usebeta)
            set(handles.UPDATE_BUTTON,'string','Can now calculate melt coefficients');
        else            
            set(handles.UPDATE_BUTTON,'string','Uh oh. Try that again');
        end
        cd(mydir)
        guidata(hObject, handles);
    case 3 %calculate melt coefficients
        [ready,usebeta]= calculateMeltCoefs(glacier);
        if ready==length(usebeta)
            set(handles.UPDATE_BUTTON,'string','Can now fit balance gradients');
        else            
            set(handles.UPDATE_BUTTON,'string','Not converged, rerun');
        end
        cd(mydir)
        guidata(hObject, handles);
    case 4 %fit balance gradient to fix missing measurements, or really late/early measurements
        use=NaN*ones(12,1); %initialize
    use(1)=get(handles.SITE1,'Value');
    use(2)=get(handles.SITE2,'Value');
    use(3)=get(handles.SITE3,'Value');
    use(4)=get(handles.SITE4,'Value');
    use(5)=get(handles.SITE5,'Value');
    use(6)=get(handles.SITE6,'Value');
    use(7)=get(handles.SITE7,'Value');
    use(8)=get(handles.SITE8,'Value');
    use(9)=get(handles.SITE9,'Value');
    use(10)=get(handles.SITE10,'Value');
    use(11)=get(handles.SITE11,'Value');
    use(12)=get(handles.SITE12,'Value');
    if get(handles.ALLSITES,'Value')==1
        use(1:end,1)==1;
    end
    
    [mysit,balance_sites,extra_sites,ready] =use_sites(glacier,use);
        ready=fitBalGrad(glacier,balance_sites,grad);
        
        if ready==1
            set(handles.UPDATE_BUTTON,'string','Ready to GO!');
        else %should never get here
            set(handles.UPDATE_BUTTON,'string','Oops! try again...');
        end
        cd(mydir)
        guidata(hObject, handles);     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%          Calculating MB         %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  afer users has defined glacier for processing %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GO_BUTTON_Callback(hObject, eventdata, handles)
% Get user input from GUI
mydir = pwd;%get(handles.MYDIR_EDIT,'string');
glacier = str2mat(get(handles.Call_Glacier,'String'));
if isempty(glacier) %if user did not enter a glacier name
    x = inputdlg('Enter glacier name',...
             'Sample', [1 50]);
glacier = str2num(x{:}); 
end
disp(['%%%%%%%%%%%%   Processing ',glacier,' Glacier Data   %%%%%%%%%%%%%%%%'])
grad = get(handles.GRAD_POPUP,'Value');
dbstop if error

    col1=[1,.5,0];
    
    use=NaN*ones(12,1); %initialize
    use(1)=get(handles.SITE1,'Value');
    use(2)=get(handles.SITE2,'Value');
    use(3)=get(handles.SITE3,'Value');
    use(4)=get(handles.SITE4,'Value');
    use(5)=get(handles.SITE5,'Value');
    use(6)=get(handles.SITE6,'Value');
    use(7)=get(handles.SITE7,'Value');
    use(8)=get(handles.SITE8,'Value');
    use(9)=get(handles.SITE9,'Value');
    use(10)=get(handles.SITE10,'Value');
    use(11)=get(handles.SITE11,'Value');
    use(12)=get(handles.SITE12,'Value');
    use(12)=get(handles.SITE13,'Value');
    use(12)=get(handles.SITE14,'Value');
    
   if get(handles.ALLSITES,'Value')==1
        use(1:end,1)=1;
   end
   if sum(use)==0
        errordlg('No sites selected')
        pause('on')
    end
    [mysit,balance_sites,extra_sites,ready] =use_sites(glacier,use);
    disp('                                                 ')
    disp(['%%%%%%%%%% You Have Selected To Use The Following %%%%%%%%%%%'])
    disp(balance_sites.')
    
    
% stope(here)

my_yr = eval(get(handles.YEAR_EDIT,'String'));
numYrs=length(my_yr);
% stop(here)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Defining Time-system used  %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch get(handles.BAL_POPUP,'Value') 
    case 1 %Stratigraphic time-system used at each site. Site mass min to site mass min. Site balance is models using P/PDD model if missing
        bal=1;
        
        ylab='Stratigraphic Year B_a (m w.e.)';
        bal2='Stratigraphic Mass Balance';
        bal_name='Stratigraphic';
    case 2 %glacier-wide from October 1 thru Sept 30. AKA the Hydrologic year
        bal=2;
        ylab='Hyrologic Year B_a (m w.e.)';
        bal2= 'Hyrologic Year Mass Balance';
        bal_name='Hyrologic'
    case 3 %Combined date system uses the the P/PDD model to find the weighted mean mass minimum date to calculate the glacier-wide balnce for one date
        bal=3;
        ylab='Fixed-Date B_a (m w.e.)';
        bal2= 'Fixed-Date Mass bBalance';
        bal_name='Fixed_Date';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Defining the surface used %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch get(handles.SURF_POPUP,'Value') 
    case 1 %is off current geometry == conventional
        surf=1;
    case 2 %is off oldest geometry == reference
        surf=2;
    case 3 % want both
        surf=[1,2];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% If you want to plot Integration %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if get(handles.PLOTINTERATION,'Value')==1
%     close all
   if length(my_yr)>1 %warn user plotting integrations for every year might take a while
       promptMessage = sprintf('You selected to plot integrations for more than one year. This might take a while. Do you still want to plot integrations?');
        titleBarCaption = 'Plotting Integrations';
        choice = questdlg(promptMessage,titleBarCaption,'Yes','No','Don''t know');
                if strcmpi(choice, 'Yes')
                    plot_integration = 1;
                elseif strcmpi(choice, 'No')
                    plot_integration=0;
                end
   else      
        plot_integration=1;
   end
else
        plot_integration=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  If you want to plot the ablation model  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if get(handles.PLOTABLATIONMODEL,'Value')==1
       % close all
       if length(my_yr)>1 %warn user plotting integrations for every year might take a while
       promptMessage = sprintf('You selected to plot the ablation model for more than one year. This probably blow up your machine. Do you still want to?');
        titleBarCaption = 'Plotting ablation model';
        choice = questdlg(promptMessage,titleBarCaption,'Yes','No','Don''t know');
                if strcmpi(choice, 'Yes')
                    plot_ablation_model = 1;
                    msgbox('Ok... Take cover!')
                elseif strcmpi(choice, 'No')
                    plot_ablation_model=0;
                end
   else      
        plot_ablation_model=1;
   end
else
        plot_ablation_model=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
cd(mydir)
addpath functions

% initialize
b_barCumu2=zeros(length(surf),numYrs);
b_bar2=zeros(length(surf),numYrs);
b_barW2=zeros(length(surf),numYrs);
b_barS2=zeros(length(surf),numYrs);
dateAMin2=zeros(length(surf),numYrs);
dateAMax2=zeros(length(surf),numYrs);
result2=zeros(7*length(surf),numYrs);
max_date=char(ones(numYrs,8*length(surf)));
min_date=char(ones(numYrs,8*length(surf)));
theNaN=char(ones(numYrs,8));
for j=1:numYrs
    theNaN(j,:)='  NaN   ';
end
for i=1:length(surf)
    if surf(i)==1 % for screen output
        surf2='CONVENTIONAL'; 
        surf3='Conventional';
        col='.-';
        col2=col1;
    elseif surf(i)==2
        surf2='REFERENCE-SURFACE';
        surf3='Reference_Surface';
        col='.-';
        col2=[.7 0 1];
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Start of glacier-wide calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%want glacier wide minimum, so callGlac=1    
% want to use bal grad fit to fill missing values (miss=1)
[b_weqW,b_weqS,b_weqG, b_weqGW,b_weq,b_barW,b_barS,b_bar,ELA,zone_weights,dateAMin,dateAMax,meltrateS,meltrateI,Zs,winter_gradients,summer_gradients,annual_gradients]=calcBalance(glacier,balance_sites,my_yr,bal,surf(i),grad,plot_integration,plot_ablation_model);
    b_barCumu2(i,:)=cumsum(b_bar);
    b_bar2(i,:)=b_bar;
    b_barW2(i,:)=b_barW;
    b_barS2(i,:)=b_barS;
    dateAMin2(i,:)=dateAMin;
    dateAMax2(i,:)=dateAMax;
    if grad==1
    for jj = 1:numYrs
        stake_index = find(b_weq(:,jj)~=0);
        for ii = 1:length(stake_index)-1
            winter_gradients(ii,jj) = round((b_weqW(stake_index(ii+1),jj)-b_weqW(stake_index(ii),jj)) / (Zs(stake_index(ii+1),jj)-Zs(stake_index(ii),jj)),4);
            summer_gradients(ii,jj) = round((b_weqS(stake_index(ii+1),jj)-b_weqS(stake_index(ii),jj)) / (Zs(stake_index(ii+1),jj)-Zs(stake_index(ii),jj)),4);
            annual_gradients(ii,jj) = round((b_weq(stake_index(ii+1),jj)-b_weq(stake_index(ii),jj)) / (Zs(stake_index(ii+1),jj)-Zs(stake_index(ii),jj)),4);
        end
        ELA(jj) = interp1(b_weq(stake_index,jj),Zs(stake_index,jj),0,'linear','extrap');
    end
    
        sizG = 2 + 3.* (size(winter_gradients,1));
        result_G(sizG*i-sizG+1:sizG*i,:)=[my_yr;winter_gradients;summer_gradients;annual_gradients;ELA];
    else
        sizG = 2 + 3.* (size(winter_gradients,1));
        result_G(sizG*i-sizG+1:sizG*i,:)=[my_yr;winter_gradients;summer_gradients;annual_gradients;ELA];  
    end
    result2(7*i-6:7*i,:)=[my_yr;cumsum(b_bar);b_bar;b_barW;b_barS;dateAMax;dateAMin];
  if get(handles.BAELA,'Value')==1
   
    ela_ba_lm=fitlm(ELA,b_bar);
    figure(1);hold on
    scatter(ELA,b_bar,'k','filled')
    lsline
    text(min(ELA)+100,-2,['r^2 = ',num2str(round(ela_ba_lm.Rsquared.Ordinary*100)/100)]) 
   
    title([glacier,' Glacier ELA B_a Regression'],'fontname','arial ','fontsize',14,'fontweight','bold')
    set(gca,'fontname','arial ','fontsize',14,'fontweight','bold','TickLength',[0.025 0.025],'linewidth',2)
    axis square
    box on
    set(gcf, 'PaperPositionMode', 'auto');
    print -depsc2 gates_epoch2.eps 
    hold off
end         
    if ~isnan(dateAMax)
        max_date(:,8*i-7:8*i)=datestr(result2(7*i-1,:).',2);%make in month/day/yr format
    else
        max_date(:,8*i-7:8*i)=theNaN;
    end
    min_date(:,8*i-7:8*i)=datestr(result2(7*i,:).',2);
end
if grad==1
    folder='Index_Method';
elseif grad==2
    folder='Gradient_Method';
end
seasonaloutput = ['data/',glacier,'/Output/',folder,'/Output_',glacier,'_',surf3,'_',bal_name,'_seasonal.csv'];
annualoutput = ['data/',glacier,'/Output/',folder,'/Output_',glacier,'_',surf3,'_',bal_name,'_annual.csv'];
myoutput = ['data/',glacier,'/Output/',folder,'/Output_',glacier,'_',surf3,'_',bal_name,'_oldmb.csv'];
gradientoutput = ['data/',glacier,'/Output/',folder,'/Output_',glacier,'_',surf3,'_',bal_name,'Gradients.csv'];

Date = datestr([max_date;min_date],'yyyy/mm/dd',1950);
SeasonalBalance_m_we = round([result2(4,:)';result2(5,:)'],1);
AnnualBalances_m_we(:,1) = round([result2(4,:)+result2(5,:)],1);
CumulativeBalance_m_we = round([result2(2,:)'-result2(5,:)';result2(2,:)'],1);

mb_temp = table(Date,SeasonalBalance_m_we,CumulativeBalance_m_we);
mb_out = sortrows(mb_temp);

annual_mb_out(:,1:3)=[mb_out(2:2:end,1) table(AnnualBalances_m_we) mb_out(2:2:end,3)];
writetable(mb_out,seasonaloutput)
writetable(annual_mb_out,annualoutput)
Bmax = (datenum(max_date) - datenum([result2(1,:)',ones(length(result2(1,:)'),2)]))/365 + result2(1,:)';
Bmax(Bmax>=2050) = Bmax(Bmax>=2050) - 100; %fix rollover year problem
Bmin = (datenum(min_date) - datenum([result2(1,:)',ones(length(result2(1,:)'),2)]))/365 + result2(1,:)';
Bmin(Bmin>=2050) = Bmin(Bmin>=2050) - 100;

file=fopen(myoutput,'w');
err = get(handles.ERRORBAR_CHECK,'Value');
if length(surf)==1    
    res_print=char(ones(numYrs,64));    
    for i=1:numYrs
res_print(i,:)=sprintf('%4d %8.1f %8.1f %8.1f %8.1f %11.3f %11.3f', result2(1,i),result2(4:5,i),result2(3,i),result2(2,i),Bmax(i),Bmin(i));
    end
    fprintf(1,'%s, FOR %s SURFACE\n',bal2,surf2);
    fprintf(1,'year B_winter B_summer B_annual B_cumulative max_date min_date(END)\n');
    for i=1:numYrs
        fprintf(1,'%s\n', res_print(i,:).');
    end
    file = fopen(myoutput,'w');%----------------------------------------------------------------------------------------------raw output file creation
    fprintf(file,'Glaciological mass balance time series\r\n');
    fprintf(file,'  \r\n');
    fprintf(file,'year  B_winter  B_summer  B_annual  B_cumulative  max_date  min_date\r\n');
    for i=1:numYrs
        fprintf(file,'%s\r\n', res_print(i,:).');
    end
    fclose(file); 
   
Year=result_G(1,:)';
Winter_Gradients_m_we_km1=result_G(2:length(winter_gradients(:,1))+1,:)';
Summer_Gradients_m_we_km1=result_G(length(winter_gradients(:,1))+2:(length(winter_gradients(:,1))*2)+1,:)';
Annual_Gradients_m_we_km1=result_G((length(winter_gradients(:,1))*2)+2:(length(winter_gradients(:,1))*3)+1,:)';
ELA=result_G(end,:)';
Mass_Maximum_Date=max_date;
Mass_Minimum=min_date;
if grad ==1
    temp = table(Year,Winter_Gradients_m_we_km1,Summer_Gradients_m_we_km1,Annual_Gradients_m_we_km1,ELA,Mass_Maximum_Date,Mass_Minimum);

writetable(temp,gradientoutput)

 
elseif grad == 2
     temp = table(Year,Winter_Gradients_m_we_km1,Summer_Gradients_m_we_km1,Annual_Gradients_m_we_km1,ELA,Mass_Maximum_Date,Mass_Minimum);

writetable(temp,gradientoutput)
    
 
end

    % Create cumulative balance plot
    axes(handles.CUMU_AXES)
    %figure
    plot([my_yr(1)-1,my_yr].',[0,b_barCumu2].',col,'color',col2)
    hold off;
    if err==1  %plot errorbars
        [errb,geod,Badj,UAF]=plotGeodetic(glacier,bal,surf,my_yr,1,grad,b_weqG,b_weqGW,bal_name,folder,surf3,plot_ablation_model);
        if ~isempty(geod) && ~isempty(UAF)
            first=geod(end,1); %first year that zero to
            ind=find(my_yr==first);
            firstUAF=UAF(1,1); %first year that zero to
            indUAF=find(my_yr==firstUAF);
            if isempty(ind)
                addme=zeros(length(geod(:,1)),1); %first photo at year 0
                addme=zeros(length(UAF(:,1)),1); %first altimetry at year 0
            else
                addme=Badj(ind,5)*ones(length(geod(:,1)),1);%need to add to geod to make first photo agree
                addmeUAF=Badj(indUAF,5)*ones(length(UAF(:,1)),1);%need to add to geod to make first photo agree
            
            end
            errorbar([my_yr(1)-1,my_yr].',[0,b_barCumu2].',[0;errb],col,'color',col2); hold on;
            errorbar([Badj(1,1)-1;Badj(:,1)],[0;Badj(:,5)],[0;errb],col,'color',[1 .1 .4],'linewidth',1);
            errorbar(geod(:,1),geod(:,2)+addme,geod(:,3),'ks','markersize',8,'linewidth',1);
            errorbar(UAF(:,1),UAF(:,2)+addmeUAF,UAF(:,3),'ko','markersize',12,'linewidth',1);
        elseif ~isempty(geod) && isempty(UAF)
            first=geod(end,1); %first year that zero to
            ind=find(my_yr==first);
           
            
            if isempty(ind)
                addme=zeros(length(geod(:,1)),1); %first photo at year 0
                
            else
                addme=Badj(ind,5)*ones(length(geod(:,1)),1);%need to add to geod to make first photo agree
                
            
            end
            errorbar([my_yr(1)-1,my_yr].',[0,b_barCumu2].',[0;errb],col,'color',col2); hold on;
            errorbar([Badj(1,1)-1;Badj(:,1)],[0;Badj(:,5)],[0;errb],col,'color',[1 .1 .4],'linewidth',1);
            errorbar(geod(:,1),geod(:,2)+addme,geod(:,3),'ks','markersize',8,'linewidth',1);
            
        else
            errorbar([my_yr(1)-1,my_yr].',[0,b_barCumu2].',[0;errb],col,'color',col2); 
        end
        
        hold off;   
        if ~isempty(geod)
            legend(surf3,'Calibrated','Geodetic','UAF altimetry','Location','SW')
        else
            legend(surf3,'Location','SW')
        end
    elseif err==1 && bal==2 
        fprintf(1,'MEASUREMENT ERROR BARS AND GEODETIC BALANCES NOT RELEVANT\r\n');
        fprintf(1,'BECAUSE BALANCE COMPELETELY MODELLED.\r\n');
        hold off;
        legend(surf3,'Location','SW')
    elseif err==1 && bal~=2
        errb=plotGeodetic(glacier,bal,my_yr,0,grad,b_weq);
        errorbar([my_yr(1)-1,my_yr].',[0,b_barCumu2].',[0;errb],col,'color',col2);     
        hold off;   
        legend(surf3,'Location','SW')
    else
        hold off;
        legend(surf3,'Location','SW')
    end
    set(handles.CUMU_AXES,'YMinorTick','on')
    xlabel('Year')
    ylabel('Cumulative balance (m w.e.)')
    grid on
    %
    % Create net/annual balance plot
    % for added excitement, plots different color every time
    axes(handles.ANN_NET_AXES)
%figure
    r= rand(1,6);

        bar(my_yr.', b_barW.',1,'FaceColor', [r(1),r(2),1])
        hold on;
        bar(my_yr.',b_barS.',1,'FaceColor', [1,r(3),r(4)])
        hold on;
        bar(my_yr.',b_bar.',1,'FaceColor', [r(5),1,r(6)]) 
        hold off; 
        legend('Winter','Summer','Annual','Location','EO')
        ylabel(ylab);
    set(handles.ANN_NET_AXES,'YMinorTick','on')
    xlabel('Year')
    xlim([my_yr(1)-.5,my_yr(numYrs)+.5])
    grid on
else % need to keep both ref and conventional bals
    res_print1=char(ones(numYrs,64));
    res_print2=char(ones(numYrs,64));
    for i=1:numYrs
res_print1(i,:)=sprintf('%4d %9.4f %8.4f %8.4f %9.4f   %s   %s', result2(1:5,i),max_date(i,1:8),min_date(i,1:8));
res_print2(i,:)=sprintf('%4d %9.4f %8.4f %8.4f %9.4f   %s   %s', result2(8:12,i),max_date(i,9:16),min_date(i,9:16));
    end
    fprintf(1,'%s, FOR CONVENTIONAL BALANCE\n',bal2);
    fprintf(1,'NETS RUN FROM min_dat PREVIOUS BALANCE YR TO min_date THIS YEAR\n');
    fprintf(1,'year  cumul_bal  yr_bal  wint_bal  summ_bal  max_dat  min_dat(END)\n');
    for i=1:numYrs
        fprintf(1,'%s\r\n', res_print1(i,:).');
    end
    fprintf(1,'%s, FOR REFERENCE-SURFACE BALANCE\n',bal2);
    fprintf(1,'NETS RUN FROM min_dat PREVIOUS BALANCE YR TO min_date THIS YEAR\n');
    fprintf(1,'year  cumul_bal  yr_bal  wint_bal  summ_bal  max_dat  min_dat(END)\n');
    for i=1:numYrs
        fprintf(1,'%s\r\n', res_print2(i,:).');
    end
    fprintf(file,'%s, FOR CONVENTIONAL BALANCE\r\n',bal2);
    fprintf(file,'NETS RUN FROM min_dat PREVIOUS BALANCE YR TO min_date THIS YEAR\r\n');
    fprintf(file,'year  cumul_bal  yr_bal  wint_bal  summ_bal  max_dat  min_dat(END)\r\n');
    for i=1:numYrs
        fprintf(file,'%s\r\n', res_print1(i,:).');
    end
    fprintf(file,'%s, FOR REFERENCE-SURFACE BALANCE\r\n',bal2);
    fprintf(file,'NETS RUN FROM min_dat PREVIOUS BALANCE YR TO min_date THIS YEAR\r\n');
    fprintf(file,'year  cumul_bal  yr_bal  wint_bal  summ_bal  max_dat  min_dat(END)\r\n');
    for i=1:numYrs
        fprintf(file,'%s\r\n', res_print2(i,:).');
    end
    fclose(file); 

    % Create cumulative balance plot
    axes(handles.CUMU_AXES)
    plot([my_yr(1)-1,my_yr].',[0,b_barCumu2(2,:)].','.-','color',[.7 0 1]);
    hold on;
    plot([my_yr(1)-1,my_yr].',[0,b_barCumu2(1,:)].','.-','color',col1);
    hold on;
    if err==1 %&& bal~=2 %plot errorbars

        [errb,geod,Badj]=plotGeodetic(glacier,bal,my_yr,1,grad,b_weq);
        if ~isempty(geod)
            first=geod(1,1); %first year that zero to
            ind=find(my_yr==first);
            if isempty(ind)
                addme=zeros(length(geod(:,1)),1); %first photo at year 0
            else
                addme=Badj(ind,3)*ones(length(geod(:,1)),1);%need to add to geod to make first photo agree        
            end
            errorbar([my_yr(1)-1,my_yr].',[0,b_barCumu2(2,:)].',[0;errb],col,'color',[.7 0 1]);hold on;
            errorbar([Badj(1,1)-1;Badj(:,1)],[0;Badj(:,3)],[0;errb],col,'color',[1 .1 .3],'linewidth',2);
            errorbar([my_yr(1)-1,my_yr].',[0,b_barCumu2(1,:)].',[0;errb],col,'color',col1);
            errorbar(geod(:,1),geod(:,2)+addme,geod(:,3),'ks','markersize',8,'linewidth',2);
            
        else
            errorbar([my_yr(1)-1,my_yr].',[0,b_barCumu2(2,:)].',[0;errb],'.','color',[.7 0 1]); hold on;
            errorbar([my_yr(1)-1,my_yr].',[0,b_barCumu2(1,:)].',[0;errb],'.','color',col1);
        end
        
   
        if ~isempty(geod)
            legend('reference-surface','adjusted','conventional','geodetic','Location','SW')
        else
            legend('reference-surface','conventional','Location','SW')
        end
    elseif err==1 && bal==2 
        fprintf(1,'MEASUREMENT ERROR BARS AND GEODETIC BALANCES NOT RELEVANT\n');
        fprintf(1,'BECAUSE BALANCE COMPELETELY MODELLED.\n');
        hold off;
        legend('reference-surface','conventional','Location','SW')
    else
        hold off;
        legend('reference-surface','conventional','Location','SW')
    end
    grid on
    set(handles.CUMU_AXES,'YMinorTick','on')
    xlabel('Year')
    ylabel('Cumulative balance (m w.e.)')
    grid on
    %
    % Create net/annual balance plot
    % for added excitement, plots different color every time
    axes(handles.ANN_NET_AXES)
    r= rand(1,12);
    p5=0.25*ones(numYrs,1);
    
        bar(my_yr.'-p5, b_barW2(1,:).',.5,'FaceColor', [r(1),r(2),1])
        hold on;
        bar(my_yr.'-p5,b_barS2(1,:).',.5,'FaceColor', [1,r(3),r(4)])
        hold on;
        bar(my_yr.'-p5,b_bar2(1,:).',.5,'FaceColor', [r(5),1,r(6)]) 
        hold on; 
        bar(my_yr.'+p5, b_barW2(2,:).',.5,'FaceColor', [r(7),r(8),1])
        hold on;
        bar(my_yr.'+p5,b_barS2(2,:).',.5,'FaceColor', [1,r(9),r(10)])
        hold on;
        bar(my_yr.'+p5,b_bar2(2,:).',.5,'FaceColor', [r(11),1,r(12)]) 
        hold off; 
        legend('con. Winter','con. Summer','con. Annual','ref. Winter','ref. Summer','ref. Annual','Location','EO')
        ylabel(ylab);

    set(handles.ANN_NET_AXES,'YMinorTick','on')
    xlabel('Year')
    xlim([my_yr(1)-.5,my_yr(numYrs)+.5])
    grid on
end




% --- Executes on button press in BAELA.
function BAELA_Callback(hObject, eventdata, handles)
% hObject    handle to BAELA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BAELA


% --- Executes on button press in PLOTABLATIONMODEL.
function PLOTABLATIONMODEL_Callback(hObject, eventdata, handles)
% hObject    handle to PLOTABLATIONMODEL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PLOTABLATIONMODEL


% --- Executes on button press in PLOTINTERATION.
function PLOTINTERATION_Callback(hObject, eventdata, handles)
% hObject    handle to PLOTINTERATION (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PLOTINTERATION
