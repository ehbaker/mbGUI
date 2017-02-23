function [mysites,balance_sites,extra_sites,ready]=use_sites(glacier,use)
%%
% list all possible sites for a particular glacier so user can select which
% ones to use for computing balance
% leave first four spaces for index sites
% then depending on input "use" from GUI user, will print the sites using
% in the balance to a file mysites
dbstop if error
ready=0;
numnams=16; %if have more than 16 stakes ever, increase

mysites=cell(numnams,1);

                                            %NOTE!!!  if you add sites here, 
                                            %also add sites to get_siteInd.m
                                            %and revise kk distribution.
                                            %also, add synthetic data to precipratio.csv 
        
file = ['../data/',glacier,'/Input/Input_',glacier,'_Glaciological_Sites.csv'];
sites=importdata(file);
mysites(1:length(sites.textdata),1)=sites.textdata;

 
for i=1:12
    if use(i)==1   
        ready=1; %using at least one site
    end
end
if ready==1
    balance_sites=' ';
    extra_sites=' ';
    use_index=find(use==1);
    xtra_index=find(use~=1);
    balance_sites=mysites(use_index,:);
    extra_sites=[mysites(xtra_index,:)];
    
    end
    %fclose(file1); 
  
end
        
        