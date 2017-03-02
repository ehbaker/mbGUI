function [ b_barW,winter_gradient, b_barS,summer_gradient,annual_gradient,b_bar,ELA] = integrate_balance(b_weqW,b_weqS,b_weq,AADs,my_yr,stake_zs,grad,surf,plot_integration)
%This function integrates point balances over the glacier hypsometry
%Integration can be done two ways
%The index method or the gradient method
% [m,n]=find(stake_zs<1);
% for i=1:length(m)
%     stake_zs(m(i),n(i))=0;
% end
 
    
for i=1:length(my_yr)
    bin_area=[];
    this_yr_bin_centers=[];
        if surf == 1
        AAD_index = find(AADs(:,1) == my_yr(i));
            if isempty(AAD_index) %first year of data
                gl_area=NaN;
            else
                gl_area = AADs(AAD_index,2:end);
            end
        elseif surf == 2
            gl_area = AADs(2,2:end);
        end
        
bin_centers=AADs(1,2:end);     
if grad==1
   
    [zone_weights(:,i),split_elevation] = get_weights(bin_centers,gl_area,stake_zs(:,i));
 
    if plot_integration==1
    bin_area=gl_area(gl_area~=0);
    this_yr_bin_centers=bin_centers(gl_area~=0);
    mean_zone_elevation=[];
    zone_width=[];
    for k=1:length(split_elevation)
        
        if k==1
            mean_zone_elevation(k,1)=split_elevation(k,1)-(split_elevation(k,1)-(this_yr_bin_centers(1)-((this_yr_bin_centers(2)-this_yr_bin_centers(1))/2)))/2;
            zone_width(k,1)=split_elevation(k,1)-(this_yr_bin_centers(1)-(this_yr_bin_centers(2)-this_yr_bin_centers(1))/2);
        elseif k==length(split_elevation)%&& length(split_elevation)==2
             zone_width(k,1)=(split_elevation(k,1)-split_elevation(k-1,1));
             zone_width(k+1,1)=(((this_yr_bin_centers(end)+(this_yr_bin_centers(end)-this_yr_bin_centers(end-1))/2)-split_elevation(k,1)));
             mean_zone_elevation(k,1)=split_elevation(k,1)-((split_elevation(k,1)-split_elevation(k-1,1))/2);
             mean_zone_elevation(k+1,1)=split_elevation(k,1)+(((this_yr_bin_centers(end)+(this_yr_bin_centers(end)-this_yr_bin_centers(end-1))/2)-split_elevation(k,1))/2);
%         elseif k==length(split_elevation)&&length(split_elevation)>2
%             zone_width(k+1,1)=(((this_yr_bin_centers(end)+(this_yr_bin_centers(end)-this_yr_bin_centers(end-1))/2)-split_elevation(k,1)));
%             mean_zone_elevation(k+1,1)=split_elevation(k,1)+(((this_yr_bin_centers(end)+(this_yr_bin_centers(end)-this_yr_bin_centers(end-1))/2)-split_elevation(k,1))/2);
        else
            zone_width(k,1)=(split_elevation(k,1)-split_elevation(k-1,1));
            mean_zone_elevation(k,1)=mean([split_elevation(k,1) split_elevation(k-1,1)]);
        end
    end
    
        zone_area=sum(gl_area).*zone_weights(:,i);
        nan_index=[];
        zero_index=[];
        area_index=[];
        areas=[];
        nan_index=~isnan(zone_area);
        zero_index=zone_area~=0;
        area_index=find(zero_index==1 & nan_index==1);
        figure(my_yr(i));hold on 
        subplot(1,3,1)
        title('Winter Balance Integration')
        hold on
        yyaxis left
        areas=zone_area(area_index,1);
        for j=1:length(mean_zone_elevation)
            bar(mean_zone_elevation(j,1),areas(j,1),zone_width(j,1),'FaceColor',[.5 .5 .5])
        end
        ylabel('Area (km^2)')
        yyaxis right
        scatter(stake_zs(area_index,i),b_weqW(area_index,i),200,'s','filled')
        ylabel('balance (m w.e.)')
        ylim([round(min(min(b_weqW)))-1 round(max(max(b_weqW)))+1])
        xlim([bin_centers(1)-50 bin_centers(end)+50])
        axis square
        box on
        
        subplot(1,3,2)
        title('Summer Balance Integration')
        hold on
        yyaxis left
        areas=zone_area(area_index,1);
        for j=1:length(mean_zone_elevation)
            bar(mean_zone_elevation(j,1),areas(j,1),zone_width(j,1),'FaceColor',[.5 .5 .5])
        end
        ylabel('Area (km^2)')
        yyaxis right
        scatter(stake_zs(area_index,i),b_weqS(area_index,i),200,'s','filled')
        ylabel('balance (m w.e.)')
        ylim([round(min(min(b_weqS)))-1 round(max(max(b_weqS)))+1])
        xlim([bin_centers(1)-50 bin_centers(end)+50])
        axis square
        box on
        
        subplot(1,3,3)
        title('Annual Balance Integration')
        hold on
        yyaxis left
        areas=zone_area(area_index,1);
        for j=1:length(mean_zone_elevation)
            bar(mean_zone_elevation(j,1),areas(j,1),zone_width(j,1),'FaceColor',[.5 .5 .5])
        end
        ylabel('Area (km^2)')
        yyaxis right
        scatter(stake_zs(area_index,i),b_weq(area_index,i),200,'s','filled')
        ylabel('balance (m w.e.)')
        ylim([round(min(min(b_weq)))-1 round(max(max(b_weq)))+1])
        xlim([bin_centers(1)-50 bin_centers(end)+50])
        axis square
        box on
        set(gcf, 'PaperPosition', [0 0 3 6]);
        print -depsc2 gates_epoch2.eps
        
    end
end
end
    
    

if grad==1
    b_barW=nansum(b_weqW.*zone_weights);
    b_barS=nansum(b_weqS.*zone_weights);
    b_bar =nansum(b_weq .*zone_weights);
    winter_gradient=NaN*ones(1,length(b_weq(1,:)));
    summer_gradient=NaN*ones(1,length(b_weq(1,:)));
    annual_gradient=NaN*ones(1,length(b_weq(1,:)));
    ELA=nan;
    
        
elseif grad==2
   
for i=1:length(my_yr)       
    stake_index=find(stake_zs(:,i)>0);
    bin_weights(:,1)=(gl_area(1:end)./sum(gl_area(1:end)));
    % get winter balance gradient and balance
    if isnan(b_weqW)
        winter_gradient(i)=nan;
         b_barW(i)=nan;
    else
    winter_lm=fitlm(stake_zs(stake_index,i),b_weqW(stake_index,i));
    winter_gradient(i)=double(cell2mat(table2cell(winter_lm.Coefficients(2,1))))*1000; 
    bw_bin_values=predict(winter_lm,bin_centers');
    weighted_bw_bin_values(:,1)=bin_weights.*bw_bin_values;
    b_barW(i)=nansum(weighted_bw_bin_values(:,1));
    end
    % get summer balance gradient and balance
    if isnan(b_weqS)
        summer_gradient(i)=nan;
         b_barS(i)=nan;
    else
    summer_lm=fitlm(stake_zs(stake_index,i),b_weqS(stake_index,i));
    summer_gradient(i)=double(cell2mat(table2cell(summer_lm.Coefficients(2,1))))*1000;
    bs_bin_values=predict(summer_lm,bin_centers');
    weighted_bs_bin_values(:,1)=bin_weights.*bs_bin_values;
    b_barS(i)=nansum(weighted_bs_bin_values(:,1));    
    end
    % get annual gradient and balance
    if isnan(b_weq)
        annual_gradient(i)=nan;
        b_bar(i)=nan;
    else
    annual_lm=fitlm(stake_zs(stake_index,i),b_weq(stake_index,i));
    annual_gradient(i)=table2array(annual_lm.Coefficients(2,1));
    annual_gradient_yintercept(i)=table2array(annual_lm.Coefficients(1,1));
    ELA(i)=(0+abs(annual_gradient_yintercept(i)))/annual_gradient(i);
    annual_gradient(i)=annual_gradient(i)*1000;
    ba_bin_values=predict(annual_lm,bin_centers');
    weighted_ba_bin_values(:,1)=bin_weights.*ba_bin_values;
    b_bar(i)=nansum(weighted_ba_bin_values(:,1));  
    end
    if plot_integration==1
        figure(my_yr(i));hold on 
        subplot(1,3,1)
        title('Winter Balance Integration')
        hold on
        yyaxis left
        bar(bin_centers,gl_area,1,'FaceColor',[.5 .5 .5])
        ylabel('Area (km^2)')
        yyaxis right
        scatter(stake_zs(stake_index,i),b_weqW(stake_index,i),'b','fill');
        plot(bin_centers,predict(winter_lm,bin_centers'),'LineWidth',1,'color','b')
        ylabel('balance (m w.e.)')
        ylim([round(min(min(b_weqW)))-1 round(max(max(b_weqW)))+1])
        xlim([bin_centers(1)-50 bin_centers(end)+50])
        axis square
        box on
        
        subplot(1,3,2)
        title('Summer Balance Integration')
        hold on
        yyaxis left
        bar(bin_centers,gl_area,1,'FaceColor',[.5 .5 .5])
        ylabel('Area (km^2)')
        yyaxis right
        scatter(stake_zs(stake_index,i),b_weqS(stake_index,i),'r','fill');
        plot(bin_centers,predict(summer_lm,bin_centers'),'LineWidth',1,'color','b')
        ylabel('balance (m w.e.)')
        ylim([round(min(min(b_weqS)))-1 round(max(max(b_weqS)))+1])
        xlim([bin_centers(1)-50 bin_centers(end)+50])
        axis square
        box on
        
        subplot(1,3,3)
        title('Annual Balance Integration')
        hold on
        yyaxis left
        bar(bin_centers(1:end),gl_area,1,'FaceColor',[.5 .5 .5])
        ylabel('Area (km^2)')
        yyaxis right
        scatter(stake_zs(stake_index,i),b_weq(stake_index,i),'o','fill');
        plot(bin_centers,predict(annual_lm,bin_centers'),'LineWidth',1,'color','b')
        ylabel('balance (m w.e.)')
        ylim([round(min(min(b_weq)))-1 round(max(max(b_weq)))+1])
        xlim([bin_centers(1)-50 bin_centers(end)+50])
        axis square
        box on
        set(gcf, 'PaperPosition', [0 0 3 6]);
        print -depsc2 gates_epoch2.eps



    end
    
        
    end
end
    



end
