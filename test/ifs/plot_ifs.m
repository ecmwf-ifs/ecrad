% This matlab scripts plots some of the inputs to the radiation scheme
% in Fig. 1 and the outputs in Fig. 2

% Location of loadnc.m
path(path, '../common')

code = 'ecrad_meridian';
in = loadnc([code '.nc']);
cases = {[code '_noaer_out.nc'],
	 [code '_ecckd_tc_noaer_out.nc'],
	 [code '_default_out.nc'],
	 [code '_ecckd_tc_out.nc'],
	 [code '_expexp_out.nc'],
	 [code '_tc_out.nc'],
	 [code '_spartacus_out.nc'],
	 [code '_ecckd_spartacus_out.nc'],
	 [code '_default_out_REFERENCE.nc']};
clear out
for icase = 1:length(cases)
  out{icase} = loadnc(cases{icase});
end
case_list = [1 2 3 4 5 6 7 8 9];
leg = {'McICA no aerosols',...
       'ecCKD Tripleclouds no aerosols',...
       'McICA',...
       'ecCKD McICA',...
       'McICA Exp-Exp',...
       'Tripleclouds',...
       'SPARTACUS',...
       'ecCKD SPARTACUS',...
       'McICA REFERENCE'};

styles = {'b','b--','r','r--','g','m','c','k','k--'};

p = 0.01.*0.5.*median(in.pressure_hl(1:end-1,:)+in.pressure_hl(2:end,:),2);

figure(1)
clf
set(gcf,'defaultlinelinewidth',1);
set(gcf,'units','inches','paperposition',[0.5 0.5 15 30]);
nplot = 4;

subplot(nplot,1,1)
plot(in.lat,in.cos_solar_zenith_angle,'k')
xlim([-90 90]);
ylabel('Cosine of solar zenith angle');

subplot(nplot,1,2)
plot(in.lat,in.skin_temperature-273.15,'r');
hold on
plot(in.lat,in.temperature_hl(end-1,:)-273.15,'b');
legend('Skin temperature','Air between first two model levels');
xlim([-90 90]);
ylabel('Temperature (\circC)');
ylim([-50 60]);

subplot(nplot,1,3)
%contour(in.lat,p,in.q_ice,10.^[-5 -5],'b');
%hold on
%contour(in.lat,p,in.q_liquid,10.^[-5 -5],'r');
contourf(in.lat,p,in.cloud_fraction,[0.05:0.1:0.95],'k');
%shading flat
colormap(jet)
set(gca,'ydir','reverse');
ylim([0 1013]);
xlim([-90 90]);
ylabel('Pressure (hPa)');
text(-90, 0, [' ' 10 '  Cloud fraction'],'verticalalignment','top')

subplot(nplot,1,4)
sea_salt = 1e9.*sum(squeeze(in.aerosol_mmr(end,1:3,:)),1);
dust = 1e9.*sum(squeeze(in.aerosol_mmr(end,4:6,:)),1);
organics = 1e9.*sum(squeeze(in.aerosol_mmr(end,7:8,:)),1);
black_carbon = 1e9.*sum(squeeze(in.aerosol_mmr(end,9:10,:)),1);
sulphate = 1e9.*sum(squeeze(in.aerosol_mmr(end,11:12,:)),1);

plot(in.lat,sea_salt,'b');
hold on
plot(in.lat,dust,'r');
plot(in.lat,organics,'g');
plot(in.lat,black_carbon,'k');
plot(in.lat,sulphate,'m');
legend('Sea salt','Dust','Organics','Black carbon','Sulphate');
ylabel('Aerosol mass mixing ratio (\mug kg^{-1})');
xlabel('Latitude (\circN)');
xlim([-90 90]);
set(gca,'yscale','log');

figure(2)
clf
set(gcf,'defaultlinelinewidth',1);
set(gcf,'units','inches','paperposition',[0.5 0.5 15 30]);
nplot = 4;
subplot(nplot,1,1)
for icase = case_list
  plot(in.lat, out{icase}.flux_dn_sw(end,:),styles{icase})
  hold on
%  plot(in.lat, out{icase}.flux_dn_sw_clear(end,:),[styles{icase} '--'])
%  plot(in.lat, out{icase}.flux_dn_direct_sw(end,:),[styles{icase} '--'])
end
xlim([-90 90]);
ylabel('Surface SW down (W m^{-2})');
legend(leg{case_list})

subplot(nplot,1,2)
for icase = case_list
  cre =  out{icase}.flux_dn_sw(1,:)-out{icase}.flux_dn_sw_clear(1,:) ...
	-out{icase}.flux_up_sw(1,:)+out{icase}.flux_up_sw_clear(1,:);
  plot(in.lat, cre,styles{icase})
  hold on
end
xlim([-90 90]);
ylabel('SW cloud radiative effect (W m^{-2})');

subplot(nplot,1,3)
for icase = case_list
  cre =  out{icase}.flux_dn_lw(1,:)-out{icase}.flux_dn_lw_clear(1,:) ...
	-out{icase}.flux_up_lw(1,:)+out{icase}.flux_up_lw_clear(1,:);
  plot(in.lat, cre,styles{icase})
  hold on
end
xlim([-90 90]);
ylabel('LW cloud radiative effect (W m^{-2})');

subplot(nplot,1,4)
for icase = case_list
  plot(in.lat, out{icase}.cloud_cover_sw,styles{icase});
  hold on
end
xlim([-90 90]);
ylabel('Cloud cover');
xlabel('Latitude (\circN)');
