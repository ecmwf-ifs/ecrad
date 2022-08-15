% Location of loadnc.m
path(path, '../common')

do_plot_sw_libradtran = 1;
do_plot_sw_spartacus_extras = 0; % Need to have run "make i3rc_spartacus_extra"

% Load libRadTran benchmark
load i3rc_mls_cumulus_LIBRADTRAN

% Load ECRAD/SPARTACUS cases
sp_code = 'i3rc_mls_cumulus';
scases = {[sp_code '_ECRAD_ICA_OUT.nc'],... % Stored result of running ECRAD in ICA mode
	  [sp_code '_3reg_3d_out.nc'],...
	  [sp_code '_3reg_3d_clustering_out.nc'],...
	  [sp_code '_mcica_out.nc'],...
	  [sp_code '_tc_out.nc'],...
	  [sp_code '_3reg_1d_out.nc']};
id_ica = 1;
id_1d = 6;
id_3d = 2; % Maximum entrapment
%id_3d = 7; % Computed entrapment
ecrad_sw_list = [id_ica id_3d 4 5];
ecrad_legend = {'ECRAD ICA','ECRAD SPARTACUS','ECRAD McICA','ECRAD Tripleclouds'};
ecrad_styles = {'b--','r--','g--','c--'};

if do_plot_sw_spartacus_extras
  scases(end+1:end+6)={[sp_code '_3reg_3d_computed_out.nc'],...
		       [sp_code '_3reg_3d_minimum_out.nc'],...
		       [sp_code '_3reg_1d_computed_out.nc'],...
		       [sp_code '_3reg_1d_minimum_out.nc'],...
		       [sp_code '_3reg_3d_computedleast_out.nc'],...
		       [sp_code '_3reg_3d_fractal_out.nc']};
  ecrad_sw_list = [id_3d 7 8 id_1d 9 10 5 11 12];
  ecrad_legend = {'SPARTACUS 3D Max','SPARTACUS 3D Comp','SPARTACUS 3D Min',...
		  'SPARTACUS 1D Max','SPARTACUS 1D Comp','SPARTACUS 1D Min',...
		  'Tripleclouds','SPARTACUS 3D Comp least','SPARTACUS 3D Fract'};
  ecrad_styles = {'r--','r-','r-.','b--','bo','bx','c--','m--','m-.'};
  id_3d_min = 8;
end

% Factor for conversion to heating rates in K day-1
ff = 24.*3600.*(9.81./1004);
for icase = 1:length(scases)
  sp{icase} = loadnc([scases{icase}]);
  sp{icase}.hr_sw = ff.*diff(sp{icase}.flux_up_sw-sp{icase}.flux_dn_sw)./diff(sp{icase}.pressure_hl);
  sp{icase}.hr_lw = ff.*diff(sp{icase}.flux_up_lw-sp{icase}.flux_dn_lw)./diff(sp{icase}.pressure_hl);
end

% Load input data
sp_input = loadnc([sp_code '_sza.nc']);
sp_sza = acosd(sp_input.cos_solar_zenith_angle);


sza_axis = [0 90];
sza_tick = [0:15:90];

% Plot shortwave 3D effect as shown in Fig. 4 of Hogan et al. (2016)
figure(1)
clf
set(gcf,'defaultlinelinewidth',1,'paperposition',[0.25 2.5 21 16]);

subplot(2,2,1);
if do_plot_sw_libradtran
  plot(sza,up_toa_1D,'b');
  hold on
  errorbar(sza,up_toa_3D,up_toa_std_3D,'r');
end

for ii = 1:length(ecrad_sw_list)
  iecrad = ecrad_sw_list(ii);
  plot(sp_sza,sp{iecrad}.flux_up_sw(1,:),ecrad_styles{ii});
  hold on
end

plot(sp_sza,sp{1}.flux_up_sw_clear(1,:),'k');

xlim(sza_axis);
ylabel('TOA upwelling flux (W m^{-2})');
text(0,1.02,' \bf(a)','verticalalignment','bottom','units','normalized');
ylim([0 300]);
grid on
set(gca,'xtick',sza_tick);
xlabel('Solar zenith angle (\circ)')
ylim([0 250]);
if do_plot_sw_libradtran
  leg = {'libRadtran ICA (DISORT)',...
	 'libRadtran 3D (MYSTIC)'};
  leg([1:length(ecrad_legend)]+2) = ecrad_legend;
else
  leg = ecrad_legend;
end
leg{length(leg)+1} = 'Clear sky';
hh=legend(leg,3);
set(hh,'fontsize',7)

subplot(2,2,2);
plot(sza,dn_direct_surf_1D./cosd(sza),'b');
hold on
errorbar(sza,dn_direct_surf_3D./cosd(sza), dn_direct_surf_std_3D./cosd(sza),'r');

for ii = 1:length(ecrad_sw_list)
  iecrad = ecrad_sw_list(ii);
  plot(sp_sza,sp{iecrad}.flux_dn_direct_sw(end,:)./cosd(sp_sza'),ecrad_styles{ii});
end

plot(sp_sza,sp{id_3d}.flux_dn_direct_sw_clear(end,:)./cosd(sp_sza'),'k');

xlim(sza_axis);
ylim([0 1100]);
grid on
set(gca,'xtick',sza_tick);
xlabel('Solar zenith angle (\circ)')
ylabel('Direct flux in direction of sun (W m^{-2})');
h=legend(leg,3);
set(h,'fontsize',7)
text(0,1.02,' \bf(b)','verticalalignment','bottom','units','normalized');

up1dinterp = interp1(sza,up_toa_1D,sp_sza)';

subplot(2,2,3)
errorbar(sza,...
	 100.*(up_toa_3D-up_toa_1D)./(up_toa_1D-sw_up_clear(end,:)),...
	 100.*(up_toa_std_3D)./(up_toa_1D-sw_up_clear(end,:)),...
	 'r');
hold on
effect_3d = 100.*(sp{id_3d}.flux_up_sw(1,:)-sp{id_ica}.flux_up_sw(1,:)) ...
	    ./(sp{id_ica}.flux_up_sw(1,:)-sp{id_3d}.flux_up_sw_clear(1,:));
effect_3d(end) = NaN;
plot(sp_sza,effect_3d,'r--');

if do_plot_sw_spartacus_extras
  effect_edge = 100.*(sp{id_3d_min}.flux_up_sw(1,:)-sp{id_ica}.flux_up_sw(1,:)) ...
	      ./(sp{id_ica}.flux_up_sw(1,:)-sp{id_3d_min}.flux_up_sw_clear(1,:));
  effect_edge(end) = NaN;
  plot(sp_sza,effect_edge,'r-.');
end

xlim(sza_axis);
ylim([-40 140]);
set(gca,'xtick',sza_tick);
ylabel('3D change to cloud radiative effect (%)');
xlabel('Solar zenith angle (\circ)')
text(0,1.02,' \bf(d)','verticalalignment','bottom','units','normalized');
grid on
if do_plot_sw_spartacus_extras
  h=legend('libRadtran','SPARTACUS','SPARTACUS edge only',2);
else
  h=legend('libRadtran','SPARTACUS',2);
end
set(h,'fontsize',8)

subplot(2,2,4)
errorbar(sza,-(up_toa_3D-up_toa_1D),up_toa_std_3D,'r');
hold on
effect_3d = sp{id_3d}.flux_up_sw(1,:)-sp{id_ica}.flux_up_sw(1,:);
effect_3d(end) = NaN;
plot(sp_sza,-effect_3d,'r--');

if do_plot_sw_spartacus_extras
  effect_edge = sp{id_3d_min}.flux_up_sw(1,:)-sp{id_ica}.flux_up_sw(1,:);
  effect_edge(end) = NaN;
  plot(sp_sza,-effect_edge,'r-.');
end

xlim(sza_axis);
ylim([-25 30]);
set(gca,'xtick',sza_tick);
ylabel('3D change to cloud radiative effect (W m^{-2})');
xlabel('Solar zenith angle (\circ)')
text(0,1.02,' \bf(c)','verticalalignment','bottom','units','normalized');
grid on
if do_plot_sw_spartacus_extras
  h=legend('libRadtran','SPARTACUS','SPARTACUS edge only',2);
else
  h=legend('libRadtran','SPARTACUS',2);
end

set(h,'fontsize',8)

% Write selected cloud radiative forcings to the terminal window
isza = [1 6];
sza(isza)

crf_mystic_1d_toa = -(up_toa_1D(isza)-sw_up_clear(end,isza))
crf_mystic_3deffect = up_toa_1D(isza)-up_toa_3D(isza)

crf_spartacus_1d_toa = interp1(sp_sza,-sp{id_ica}.flux_up_sw(1,:) + sp{id_ica}.flux_up_sw_clear(1,:), sza(isza))
crf_spartacus_3deffect = interp1(sp_sza,-effect_3d(:), sza(isza))

crf_mystic_3deffect_percentage = 100.* crf_mystic_3deffect ./ crf_mystic_1d_toa
crf_spartacus_3deffect_percentage = 100.* crf_spartacus_3deffect ./ crf_spartacus_1d_toa


% Plot the 3D impact on shortwave heating rates as in Fig. 5 of Hogan et al. (2016)
figure(2)
clf
set(gcf,'units','inches','paperposition',[0.25 0.25 10 9],'defaultlinelinewidth',1);

sza_map = [1 9 16 24 31 38 41 45];
sza_map = {1, [8 9], 16, [23 24], 31, [38 39], 41, 45};

for isza = [1 5]; %[2 6]
  sza(isza)
  hr_sp_clear = ff.*diff(sp{id_3d}.flux_up_sw_clear-sp{id_3d}.flux_dn_sw_clear)./diff(sp{id_3d}.pressure_hl);
  ii = sza_map{isza};
  z_mid_sp = (sp_input.height_hl(1:end-1,1) + sp_input.height_hl(2:end,1))./2;

  plot(dat1D{isza}.hr,dat1D{isza}.z_mid,'b')
  axis([1 4.5 0 3]);
  hold on

  plot(dat3D{isza}.hr,dat3D{isza}.z_mid,'r')
  if isza == 5
    errs=nan.*dat3D{isza}.hr_std;
    ind = [2 4 6 16 26 36 43:2:49];
    ind = [1:5 8:5:38 42:50]; 
    errs(ind) = dat3D{isza}.hr_std(ind);
    herrorbar(dat3D{isza}.hr,dat3D{isza}.z_mid,errs,'r')
  end
  
  SPARTACUS_SZA = sp_sza(ii);
  plot(mean(sp{id_ica}.hr_sw(:,ii),2),z_mid_sp./1000,'b--');
  plot(mean(sp{id_3d}.hr_sw(:,ii),2),z_mid_sp./1000,'r--');
  plot(mean(hr_sp_clear(:,ii),2),z_mid_sp./1000,'k-');
end

axis([1 4.5 0 3]);
xlabel('Shortwave heating rate (K d^{-1})');
ylabel('Height (km)');
text(1.25,3,['\theta_0=' num2str(sza(isza)) '\circ'],'verticalalignment','bottom')
text(2.4,3,'\theta_0=0\circ','verticalalignment','bottom')
h=legend(leg{[1 2 3 4 end]},4);
set(h,'fontsize',7);

% Plot longwave 3D effect on heating rates as shown in Fig. 6 of Hogan
% et al. (2016)
p = interp1(sp_input.height_hl(:,1),sp_input.pressure_hl(:,1),z.*1000);
z_mid = 0.5.*(z(1:end-1)+z(2:end));

hr_clear = ff.*diff(lw_up_clear-lw_dn_clear)./diff(p);
hr_1d = ff.*diff(d{1}.lw_up-d{1}.lw_dn)./diff(p);
hr_3d = ff.*diff(d{2}.lw_up-d{2}.lw_dn)./diff(p);
hr_sp_3d = ff.*diff(sp{id_3d}.flux_up_lw-sp{id_3d}.flux_dn_lw)./diff(sp{id_3d}.pressure_hl);
hr_sp_1d = ff.*diff(sp{id_1d}.flux_up_lw-sp{id_1d}.flux_dn_lw)./diff(sp{id_1d}.pressure_hl);
hr_sp_clear = ff.*diff(sp{id_1d}.flux_up_lw_clear-sp{id_1d}.flux_dn_lw_clear)./diff(sp{id_1d}.pressure_hl);
hr_sp_3d_clustering = ff.*diff(sp{3}.flux_up_lw-sp{3}.flux_dn_lw)./diff(sp{3}.pressure_hl);

cols = 'brmg';
cols_clear = 'km';

figure(3)
clf
set(gcf,'units','inches','paperposition',[0.25 0.25 10 9],'defaultlinelinewidth',1);

plot(smooth1D(hr_clear),z_mid,'k');
hold on
plot(smooth1D(hr_1d),z_mid,cols(1));
plot(smooth1D(hr_3d),z_mid,cols(2));
plot(smooth1D(hr_sp_clear(:,1)),z_mid_sp./1000,'k--');
plot(smooth1D(hr_sp_1d(:,1)),z_mid_sp./1000,[cols(1) '--']);
plot(smooth1D(hr_sp_3d_clustering(:,1)),z_mid_sp./1000,[cols(2) '--']);
h=plot(smooth1D(hr_sp_3d(:,1)),z_mid_sp./1000,[cols(3) '--']);
if cols(3) == 'm'
  set(h,'color',[1 0.65 0]);
end

ylim([0 3]);

xlabel('Longwave heating rate (K d^{-1})');
ylabel('Height (km)');
xlim([-5 -1]);
h=legend('libRadtran clear','libRadtran ICA (DISORT)','libRadtran 3D (MYSTIC)',...
	 'SPARTACUS clear','SPARTACUS 1D','SPARTACUS 3D clust', ...
	 'SPARTACUS 3D no clust',3);
set(h,'fontsize',6);
