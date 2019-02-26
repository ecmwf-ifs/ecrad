figure(1)
clf
set(gcf,'defaultaxeslayer','top',...
    'defaultlinelinewidth',1,...
    'units','inches','paperposition',[0.5 0.5 25 22]);
for iscene = 1 %:2

code = ['mls_london' num2str(iscene) '_bheight'];

% Load input data
input = loadnc([code '.nc']);
bheight = input.canopy_depth(3,:);

street_area_fraction = (1.0-input.building_fraction(end,1))
norm_perimeter = 4.0.*street_area_fraction.*(1.0-street_area_fraction)./input.building_scale(end,1)
wall_area_fraction = norm_perimeter.*input.canopy_depth(end,:)

wall_scaling = street_area_fraction ./ wall_area_fraction;
wall_scaling = wall_area_fraction ./ street_area_fraction;
%wall_scaling = 1;
H = input.canopy_depth(end,1);
S = input.building_scale(end,1);

% Load output data
cases = {'','_vac'};
%	 [code '_bheight_surf_vac_out.nc']};
for icase = 1:length(cases)
  rad{icase}    = loadnc([code '_surf' cases{icase} '_out.nc'])
  output{icase} = loadnc([code cases{icase} '_out.nc']);
end

bheight_axis = [0:5:50];

id = 2;
flist = [3:5];

%figure(1)
%clf

iground = 3;
iroof = 4;
iwall = 5;
icanopy = 2;
icanopy2 = 1;
cols = {'m','b','r','k','g'};
if 0
subplot(2,2,iscene)
plot(bheight,rad{id}.flux_dn_sw_facet(iroof,:),cols{iroof})
hold on
plot(bheight,rad{id}.flux_dn_sw_facet(iground,:),cols{iground})
plot(bheight,rad{id}.flux_dn_sw_facet(iwall,:).*wall_scaling,cols{iwall})
plot(bheight,rad{id}.flux_dn_direct_sw_facet(iroof,:),[cols{iroof} '--'])
plot(bheight,rad{id}.flux_dn_direct_sw_facet(iground,:),[cols{iground} '--'])
plot(bheight,rad{id}.flux_dn_direct_sw_facet(iwall,:).*wall_scaling,[cols{iwall} '--'])

%plot(bheight,rad{id}.flux_up_sw_facet(flist,:),':')
%plot(bheight,rad{1}.absorption_sw_canopy(end,:),'k');
%plot(bheight,sum(rad{id}.flux_dn_sw_facet(flist([1 3]),:),1),'k')
%plot(bheight,sum(rad{id}.flux_dn_direct_sw_facet(flist([1 3]),:),1),'k--')
%plot(bheight,rad{id}.flux_dn_sw_facet(flist,:)-rad{id}.flux_dn_direct_sw_facet(flist,:),'-.')
xlim([0 max(bheight)])
set(gca,'xtick',bheight_axis);
xlabel('Building height (m)');
ylabel('Incoming flux (W m^{-2})');
legend(['Roof (' num2str(1-street_area_fraction,2) ')'],...
       ['Ground (' num2str(street_area_fraction,2) ')'], ...
       ['Wall'],...
       'Roof direct','Ground direct','Wall direct',1);
title(['\bf(' 'c'-1+iscene ') Scene 1: {\itH} = ' num2str(H,'%3.1f') ' m, {\itS} = ' num2str(S,'%3.1f') ' m'])
ylim([0 max(ylim)])
end

subplot(2,2,iscene+2)
%plot(bheight,rad{id}.flux_dn_lw_facet(iroof,:),cols{iroof})
%hold on
%plot(bheight,rad{id}.flux_dn_lw_facet(iground,:),cols{iground})
%plot(bheight,rad{id}.flux_dn_lw_facet(iwall,:).*wall_scaling,cols{iwall})

%plot(bheight,rad{id}.flux_up_lw_facet(iroof,:),[cols{iroof} '--'])
%plot(bheight,rad{id}.flux_up_lw_facet(iground,:),[cols{iground} '--'])
%plot(bheight,rad{id}.flux_up_lw_facet(iwall,:).*wall_scaling,[cols{iwall} '--'])

% Compute urban canopy absorption as residual

for id = 1:length(cases)
  flux_up_lw_roof = rad{id}.flux_up_lw_facet(iroof,:);
  flux_up_lw_canopy = (output{id}.flux_up_lw(end,:)-input.building_fraction(3,1).*flux_up_lw_roof) ...
		      / (1.0-input.building_fraction(3,1));
  flux_net_lw_canopy{id} = output{id}.flux_dn_lw(end,:)-flux_up_lw_canopy;
  absorption_lw_canopy{id} ...
    = flux_net_lw_canopy{id} ...
      - (rad{id}.flux_dn_lw_facet(iground,:) - rad{id}.flux_up_lw_facet(iground,:)) ...
      - (rad{id}.flux_dn_lw_facet(iwall,:) - rad{id}.flux_up_lw_facet(iwall,:)).*wall_scaling;
end

syms = {'-','--','-.',':'};
iurban = 3;

for id = 1:length(cases)

%plot(bheight,rad{id}.flux_up_lw_facet(iroof,:)-rad{id}.flux_dn_lw_facet(iroof,:),[cols{iroof} syms{id}])
%hold on
%plot(bheight,rad{id}.flux_up_lw_facet(iground,:)-rad{id}.flux_dn_lw_facet(iground,:),[cols{iground} syms{id}])
%plot(bheight,(rad{id}.flux_up_lw_facet(iwall,:)-rad{id}.flux_dn_lw_facet(iwall,:)).*wall_scaling,[cols{iwall} syms{id}])
%plot(bheight,-rad{id}.absorption_lw_canopy(iurban,:),[cols{icanopy} syms{id}])

%plot(rad{id}.flux_up_lw_facet(iroof,:)-rad{id}.flux_dn_lw_facet(iroof,:),bheight,[cols{iroof} syms{id}])
plot(rad{id}.flux_up_lw_facet(iground,:)-rad{id}.flux_dn_lw_facet(iground,:),bheight,[cols{iground} syms{id}])
hold on
plot((rad{id}.flux_up_lw_facet(iwall,:)-rad{id}.flux_dn_lw_facet(iwall,:)).*wall_scaling,bheight,[cols{iwall} syms{id}])
plot(-rad{id}.absorption_lw_canopy(iurban,:),bheight,[cols{icanopy} syms{id}])


%plot(bheight,-absorption_lw_canopy{id},[cols{icanopy2} syms{id}])
%plot(bheight,-flux_net_lw_canopy{id},[cols{icanopy} syms{id}])
%plot(bheight,rad{id}.flux_up_lw_facet(iwall,:).*wall_scaling,'g^')
%plot(bheight,rad{id}.flux_dn_lw_facet(iwall,:).*wall_scaling,'gv')

end

ylim([0 max(bheight)])
set(gca,'ytick',bheight_axis,'xtick',[-60:20:140]);
ylabel('Building height (m)');
xlabel('Net outward flux (W m^{-2})');
%legend(['Roof in (' num2str(1-street_area_fraction,2) ')'],...
%       ['Ground in (' num2str(street_area_fraction,2) ')'], ...
%       ['Wall in'],...
%       'Roof net out','Ground net out','Wall net out',1);
%legend('Roof','Ground','Wall','Canopy','Canopy (residual)');
%legend('Roof','Ground','Wall','Canopy',...
%       'Roof vacuum','Ground vacuum','Wall vacuum','Canopy vacuum');
legend('Ground','Wall','Canopy',...
       'Ground vacuum','Wall vacuum','Canopy vacuum');
%,...
%       'location','eastoutside');
%title(['\bf(' 'c'-1+iscene ') Scene 1: {\itH} = ' num2str(H,'%3.1f') ' m, {\itS} = ' num2str(S,'%3.1f') ' m'])
xlim([-60 140]);
end
