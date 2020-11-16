function stats = evaluate_ckd_lw_fluxes(ref_file, ckd_file_in, model_name, scenario_title)
% Evaluate the fluxes and heating rates from a CKD model with
% line-by-line benchmark calculations for a CKDMIP scenario

ickd_stats = 1; % Which CKD model do we show the stats for (0=none)

if iscell(ckd_file_in)
  ckd_file = ckd_file_in;
else
  ckd_file = {ckd_file_in};
end

if iscell(model_name)
  titles = model_name;
else
  titles = {model_name}
end

combined_titles = '';
for ie = 1:length(titles)
  combined_titles = [combined_titles titles{ie} '_'];
end
ickd = 1:length(ckd_file);

col_ref = 'k';
cols = {'r','b','g','m','c'};
alpha= [0.2 0.3 0.25 0.25 0.1 0.1];

dflux_axis = [-4 4];
%dfluxscatter_axis = [-8 8];
dfluxscatter_axis = [-4 4];
hr_axis = [-1 1];

if length(ckd_file) > 1
  do_legend = 1;
else
  do_legend = 0;
end

legfontsize = 7;
error_save = [];


[ref,attr_ref] = loadnc(ref_file);
if ~exist('scenario_title','var')
  scenario_title = attr_ref.global.scenario;
end

hr_ref = calc_hr(ref,'lw')';

for ie = 1:length(ckd_file)
  ckd{ie} = loadnc(ckd_file{ie});
  hr_ckd{ie} = calc_hr(ckd{ie},'lw')';
end

% Pressure in hPa at full levels and half levels
p_fl_hPa = 0.01.*0.5.*(ref.pressure_hl(1:end-1,:) + ref.pressure_hl(2:end,:));
p_hl_hPa = 0.01.*ref.pressure_hl;

iprof = 1:size(ref.pressure_hl,2);

xloc = 0.025;yloc = 0.95;

if nargout == 0
  do_plot = 1;
else
  do_plot = 0;
end

if do_plot
clf
set(gcf,'paperposition',[0.5 0.5 27 25]);
set(gcf,'defaultaxeslayer','top','defaultlinelinewidth',0.5,'defaulttextfontsize',10);

% Plot reference flux and heating-rate profiles
subplot(3,3,1)
semilogy(ref.flux_up_lw(:,iprof), p_hl_hPa(:,iprof), col_ref);
set(gca,'ydir','reverse');
ylim([0.01 1000]);
ylabel('Pressure (hPa)');
xlabel('Upwelling longwave flux (W m^{-2})');
title('Reference profiles')
text(xloc, yloc, '\bf(a)','units','normalized');

subplot(3,3,4)
semilogy(ref.flux_dn_lw(:,iprof), p_hl_hPa, col_ref);
set(gca,'ydir','reverse');
ylim([0.01 1000]);
ylabel('Pressure (hPa)');
xlabel('Downwelling longwave flux (W m^{-2})');
text(xloc, yloc, '\bf(d)','units','normalized');

subplot(3,3,7)
semilogy(hr_ref(:,iprof), p_fl_hPa(:,iprof), col_ref);
set(gca,'ydir','reverse');
ylim([0.01 1000]);
ylabel('Pressure (hPa)');
xlabel('Heating rate (K d^{-1})');
xlim([-25 5]);
text(xloc, yloc, '\bf(g)','units','normalized');

% Plot random errors and biases in flux and heating-rate profiles
subplot(3,3,8)
end

ileg = 1;
leg = {};
for ie = ickd
  err = calc_hr_error(p_hl_hPa(:,iprof), hr_ckd{ie}(:,iprof),hr_ref(:,iprof), [0.02 4]);
  err4 = calc_hr_error(p_hl_hPa(:,iprof), hr_ckd{ie}(:,iprof),hr_ref(:,iprof), [4 1100]);
  leg{ileg} = [titles{ie} ' (RMSE=' num2str(err,'%0.3f') ' K d^{-1})'];
  ileg = ileg+1;
  error_save(1,ie) = err;
  error_save(2,ie) = err4;

  if do_plot
    plot(-1,-1,cols{ie},'linewidth',1.5);
    hold on
  end
end


if do_plot
for ie = ickd
  hr_error = hr_ckd{ie}(:,iprof) - hr_ref(:,iprof);
  hr_bias = mean(hr_error');
  hr_ci = std(hr_error').*1.96;
  hr_errmin = min(hr_error,[],2);
  hr_errmax = max(hr_error,[],2);
  pp = mean(p_fl_hPa');
  h=fill([hr_bias+hr_ci flip(hr_bias-hr_ci)],...
	 [pp flip(pp)],'r');
  set(h,'facecolor',cols{ie}(1),'edgecolor','none','facealpha',alpha(ie));
  hold on
  plot(hr_bias,pp,cols{ie},'linewidth',1.5);
%  plot(hr_errmin,pp,[cols{ie}(1) '--'])
%  plot(hr_errmax,pp,[cols{ie}(1) '--'])
end
plot(hr_axis,[0.02 0.02],'k:');
plot(hr_axis,[4 4],'k:');
set(gca,'ydir','reverse');
set(gca,'yscale','log');
ylim([0.01 1000]);
ylabel('Pressure (hPa)');
xlabel('Heating rate error (K d^{-1})');
xlim(hr_axis);
plot([0 0],[0.01 1000],'k:');
%if do_legend
%  set(legend(leg,'location','south'),'fontsize',legfontsize','box','on');
  %end
text(xloc, yloc, '\bf(h)','units','normalized');

subplot(3,3,2)
for ie = ickd
  plot(-1,-1,cols{ie},'linewidth',1.5);
  hold on
end
for ie = ickd
  flux_error = ckd{ie}.flux_up_lw(:,iprof) - ref.flux_up_lw(:,iprof);
  flux_bias = mean(flux_error');
  flux_ci = std(flux_error').*1.96;
  flux_errmin = min(flux_error,[],2);
  flux_errmax = max(flux_error,[],2);
  pp = mean(p_hl_hPa');
  h=fill([flux_bias+flux_ci flip(flux_bias-flux_ci)],...
	 [pp flip(pp)],'r');
  set(h,'facecolor',cols{ie}(1),'edgecolor','none','facealpha',alpha(ie));
  hold on
  plot(flux_bias,pp,cols{ie},'linewidth',1.5);
end
set(gca,'ydir','reverse','yscale','log');
ylim([0.01 1000]);
plot([0 0],[0.01 1000],'k:');
ylabel('Pressure (hPa)');
xlabel('Upwelling flux error (W m^{-2})');
xlim(dflux_axis);
title('Errors in profiles')
if do_legend
  set(legend(titles{ickd},'location','best'),'fontsize',legfontsize','box','on');
end
text(xloc, yloc, '\bf(b)','units','normalized');

subplot(3,3,5)
for ie = ickd
  flux_error = ckd{ie}.flux_dn_lw(:,iprof) - ref.flux_dn_lw(:,iprof);
  flux_bias = mean(flux_error');
  flux_ci = std(flux_error').*1.96;
  flux_errmin = min(flux_error,[],2);
  flux_errmax = max(flux_error,[],2);
  pp = mean(p_hl_hPa');
  h=fill([flux_bias+flux_ci flip(flux_bias-flux_ci)],...
	 [pp flip(pp)],'r');
  set(h,'facecolor',cols{ie}(1),'edgecolor','none','facealpha',alpha(ie));
  hold on
  plot(flux_bias,pp,cols{ie},'linewidth',1.5);
end
set(gca,'ydir','reverse','yscale','log');
ylim([0.01 1000]);
plot([0 0],[0.01 1000],'k:');
ylabel('Pressure (hPa)');
xlabel('Downwelling flux error (W m^{-2})');
xlim(dflux_axis);
text(xloc, yloc, '\bf(e)','units','normalized');

subplot(3,3,3)
ileg = 1;
leg = {};
end

for ie = ickd
  mydiff = ckd{ie}.flux_up_lw(1,iprof)-ref.flux_up_lw(1,iprof);
  err = sqrt(mean(mydiff.^2));
  leg{ileg} = [titles{ie} ' (RMSE=' num2str(err,'%0.2f') ' W m^{-2})'];
  error_save(3,ie) = err;
  error_save(4,ie) = mean(mydiff);
  %plot(-1,-1,cols{ie},'linewidth',1.5);
  ileg = ileg+1;
  if do_plot
    plot(ref.flux_up_lw(1,iprof), mydiff, ...
	 [cols{ie}(1) 'o']);
    hold on
  end
end
if do_plot
if do_legend
  set(legend(leg,'location','best'),'fontsize',legfontsize','box','on')
end
xlabel('Reference TOA upwelling (W m^{-2})');
ylabel('TOA upwelling error (W m^{-2})');
xlim([0 400]);
plot(xlim,[0 0],'k:','linewidth',0.5);
ylim(dfluxscatter_axis);
title('Errors at surface and TOA');
text(xloc, yloc, '\bf(c)','units','normalized');

subplot(3,3,6)
end
ileg = 1;
leg = {};
for ie = ickd
  mydiff = ckd{ie}.flux_dn_lw(end,iprof)-ref.flux_dn_lw(end,iprof);
  err = sqrt(mean(mydiff.^2));
  leg{ileg} = [titles{ie} ' (RMSE=' num2str(err,'%0.2f') ' W m^{-2})'];
  error_save(5,ie) = err;
  error_save(6,ie) = mean(mydiff);
  ileg = ileg+1;
  if do_plot
    plot(ref.flux_dn_lw(end,iprof), mydiff, ...
	 [cols{ie}(1) 'o']);
    hold on
  end
end
if do_plot
if do_legend
  set(legend(leg,'location','best'),'fontsize',legfontsize')
end
xlabel('Reference surface downwelling (W m^{-2})');
ylabel('Surface downwelling error (W m^{-2})');
plot(xlim,[0 0],'k:','linewidth',0.5);
ylim(dfluxscatter_axis);
text(xloc, yloc, '\bf(f)','units','normalized');

h=subplot(3,3,9);
set(h,'visible','off')

xstart = -0.15;
text(xstart,0.95,['\bfScenario: ' scenario_title],'fontsize',12);

if ickd_stats > 0
  for ie = ickd_stats
    text(xstart,0.85,['\bfCKD model: ' titles{ie}],'fontsize',12);
    ystart = 0.75; %-0.4.*(ie-1);
    %text(xstart,ystart,['\bfRMS errors in ' titles{ie} ':'],'units','normalized');
    text(xstart,ystart,['Bias TOA upwelling: ' num2str(error_save(4,ie),'%0.2f') ' W m^{-2}'],'units','normalized');
    text(xstart,ystart-0.1,['Bias surface downwelling: ' num2str(error_save(6,ie),'%0.2f') ' W m^{-2}'],'units','normalized');
    text(xstart,ystart-0.2,['RMSE TOA upwelling: ' num2str(error_save(3,ie),'%0.2f') ' W m^{-2}'],'units','normalized');
    text(xstart,ystart-0.3,['RMSE surface downwelling: ' num2str(error_save(5,ie),'%0.2f') ' W m^{-2}'],'units','normalized');
    text(xstart,ystart-0.4,['RMSE heating rate (0.02-4 hPa):  ' num2str(error_save(1,ie),'%0.3f') ' K d^{-1}'],'units','normalized');
    text(xstart,ystart-0.5,['RMSE heating rate (4-1100 hPa):  ' num2str(error_save(2,ie),'%0.3f') ' K d^{-1}'],'units','normalized');
    
  end
end

drawnow
end
%print('-dpng','-painters','-r100',[combined_titles 'evaluation1_fluxes_' scenario '.png']);
%end

if nargin > 0
  stats.toa_up_bias = error_save(4,:);
  stats.surf_dn_bias= error_save(6,:);
  stats.toa_up_rmse = error_save(3,:);
  stats.surf_dn_rmse= error_save(5,:);
  stats.heating_rate_low_rmse=error_save(1,:);
  stats.heating_rate_high_rmse=error_save(2,:);
end

