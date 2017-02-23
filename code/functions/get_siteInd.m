function [kk]=get_siteInd(site,glacier) %kk is the "index number" of a site???
% function gets proper index in precip ratio (ablationmodel.m, modelwhole.m) 
% and for inverting for precip ratios (invParams.m, getXYprec.m, getXYmelt.m)
% and for H0 equation(internalAcc.m) and for  in file for site
% Here is where we decide if there's enough data to warrant it's own site
% designation
% Gulkana (glacier==1)'LA' has 8 years so give an index
% Gulkana (glacier==1)'AB' should get own index if continued, else close to A?
% Gulkana (glacier==1)'AE' in rocks very close to A
% Gulkana (glacier==1)'W'and 'X' and 'V' should get own index if continued,
% else close to old C?
% Wolverine(glacier==0)'LA' was stopped after 2.5 yrs, not own index (close to A)
% Wolverine(glacier==0)'AU' should get own index if continued, maybe T (terminus pole)
% Wolverine(glacier==0)'Y'and 'Z' should get own index if continued, else close to C
%
dbstop if error;
file = ['../data/',glacier,'/Input/Input_',glacier,'_Glaciological_Sites.csv'];
sites=importdata(file);
site0All=sites.textdata(:,1);                         %get vector of all site names in mb_input.txt file

site_ind=strcmp(site0All,cellstr(site));

kk=sites.data(site_ind,1);
                           
%stop(here)
end