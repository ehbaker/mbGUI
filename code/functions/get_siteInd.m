function [kk]=get_siteInd(site,glacier) %kk is the "index number" of a site???
% function gets proper index in precip ratio (ablationmodel.m, modelwhole.m) 
% and for inverting for precip ratios (invParams.m, getXYprec.m, getXYmelt.m)
% and for H0 equation(internalAcc.m) and for  in file for site
% Here is where we decide if there's enough data to warrant it's own site
% designation
%
dbstop if error;
file = ['data/',glacier,'/Input/Input_',glacier,'_Glaciological_Sites.csv'];
sites=importdata(file);
site0All=sites.textdata(:,1);                         %get vector of all site names in mb_input.txt file

site_ind=strcmp(site0All,cellstr(site));

kk=sites.data(site_ind,1);
                           
end