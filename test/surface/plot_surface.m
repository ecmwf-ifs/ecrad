figure(1)
clf
set(gcf,'defaultaxeslayer','top',...
    'defaultlinelinewidth',1,...
    'units','inches','paperposition',[0.5 0.5 22 22]);
for iscene = 1 %1:2

code = ['mls_london' num2str(iscene) ];

% Load input data
input = loadnc([code '_sza.nc']);
sza = acosd(input.cos_solar_zenith_angle);

street_area_fraction = (1.0-input.building_fraction(end,1))
norm_perimeter = 4.0.*street_area_fraction.*(1.0-street_area_fraction)./input.building_scale(end,1)
wall_area_fraction = norm_perimeter.*input.canopy_depth(end,1)

wall_scaling = street_area_fraction ./ wall_area_fraction
H = input.canopy_depth(end,1);
S = input.building_scale(end,1);

% Load output data
cases = {[code '_sza_surf_out.nc']};
for icase = 1:length(cases)
  rad{icase} = loadnc([cases{icase}])
end


sza_axis = [0 90];
sza_tick = [0:15:90];

id = 1;
flist = [3:5];

%figure(1)
%clf
if 0
subplot(2,2,1)
plot(sza,rad{id}.flux_dn_sw(end,:),'k')
hold on
plot(sza,rad{id}.flux_dn_direct_sw(end,:),'k--')
plot(sza,rad{id}.flux_up_sw(end,:),'r')
plot(sza,rad{id}.flux_up_sw_facet(flist,:),':')
plot(sza,rad{id}.flux_dn_sw(end,:)-rad{id}.flux_dn_direct_sw(end,:),'-.')
xlim([0 90])
hold on
%plot(sza,sum(rad{id}.flux_dn_direct_sw_facet(flist([1 3]),:),1),'g--')
end

iground = 3;
iroof = 4;
iwall = 5;
cols = {'','','r','k','g'};

subplot(2,2,iscene)
plot(sza,rad{id}.flux_dn_sw_facet(iroof,:),cols{iroof})
hold on
plot(sza,rad{id}.flux_dn_sw_facet(iground,:),cols{iground})
plot(sza,rad{id}.flux_dn_sw_facet(iwall,:).*wall_scaling,cols{iwall})
plot(sza,rad{id}.flux_dn_direct_sw_facet(iroof,:),[cols{iroof} '--'])
plot(sza,rad{id}.flux_dn_direct_sw_facet(iground,:),[cols{iground} '--'])
plot(sza,rad{id}.flux_dn_direct_sw_facet(iwall,:).*wall_scaling,[cols{iwall} '--'])

%plot(sza,rad{id}.flux_up_sw_facet(flist,:),':')
%plot(sza,rad{1}.absorption_sw_canopy(end,:),'k');
%plot(sza,sum(rad{id}.flux_dn_sw_facet(flist([1 3]),:),1),'k')
%plot(sza,sum(rad{id}.flux_dn_direct_sw_facet(flist([1 3]),:),1),'k--')
%plot(sza,rad{id}.flux_dn_sw_facet(flist,:)-rad{id}.flux_dn_direct_sw_facet(flist,:),'-.')
xlim([0 90])
set(gca,'xtick',0:15:90);
xlabel('Solar zenith angle (\circ)');
ylabel('Incoming flux (W m^{-2})');
legend(['Roof (' num2str(1-street_area_fraction,2) ')'],...
       ['Ground (' num2str(street_area_fraction,2) ')'], ...
       ['Wall (' num2str(wall_area_fraction,'%3.2f') ')'],...
       'Roof direct','Ground direct','Wall direct',1);
title(['\bf(' 'c'-1+iscene ') Scene 1: {\itH} = ' num2str(H,'%3.1f') ' m, {\itS} = ' num2str(S,'%3.1f') ' m'])
if 0
subplot(2,2,4)
plot(sza,rad{id}.flux_dn_direct_sw_facet(flist,:)./rad{id}.flux_dn_sw_facet(flist,:))

subplot(2,2,3)
plot(sza,rad{id}.flux_dn_sw(end,:)-rad{id}.flux_up_sw(end,:),'k')
hold on
plot(sza,street_area_fraction.*sum(rad{id}.flux_dn_sw_facet(flist([1 3]),:)-rad{id}.flux_up_sw_facet(flist([1 3]),:),1) ...
     + (1-street_area_fraction).*(rad{id}.flux_dn_sw_facet(flist(2),:)-rad{id}.flux_up_sw_facet(flist(2),:)),'r--')
plot(sza,rad{id}.flux_dn_sw_facet(flist,:)-rad{id}.flux_up_sw_facet(flist,:),'--')
end

end
