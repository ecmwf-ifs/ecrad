function evaluate_ckd_sw_fluxes(ref_file, ckd_file_in, model_name, scenario_title,specdef)
% Evaluate the shortwave fluxes and heating rates from a CKD model
% with line-by-line benchmark calculations for the four main CKDMIP
% scenarios

ickd_stats = 1; % Which CKD model do we show the stats for (0=none)
nsza = 5;isza_list = 1:nsza;
%nsza = 1;isza_list = 5;


do_diurnal_average = 0;

if do_diurnal_average
  hr_scaling_factor = 0.5;
  warning(['Heating rate error (only) has been scaled by 0.5 to represent a diurnal average']);
else
  hr_scaling_factor = 1;
end

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

if nargin > 4 & ~iscell(ckd_file_in)
  ckd_file = {ckd_file_in,ckd_file_in}
  do_fix_ssi = [0 1];
  ickd_stats = 2;
  titles = {model_name,[model_name ' (fix SSI)']};
else
  do_fix_ssi = zeros(size(ckd_file));
end


combined_titles = '';
for ie = 1:length(titles)
  combined_titles = [combined_titles titles{ie} '_'];
end
ickd = 1:length(ckd_file);

col_ref = 'k';
cols = {'r','b','g','m','c'};
alpha= [0.2 0.3 0.25 0.25 0.1 0.1];

ref_cols = {'y','g','c','b','k'};
%ref_shades = {[0.8 0.8 0.8],[0.6 0.6 0.6],[0.4 0.4 0.4],[0.2 0.2 0.2],[0 0 0]};
%ref_styles = {'-','-','-','-','-'};


dflux_axis = [-5 5];
hr_axis = [-3 3];
%hr_axis = [-15 15];

% Debugging ecCKD...
%dflux_axis = [-20 20];
%hr_axis = [-5 5];


if length(ckd_file) > 1
  cols = {'b','r','g','m','c'};
  do_legend = 1;
else
  do_legend = 0;
end
legfontsize = 7;
rmse_save = [];

  [ref_orig, attr_ref] = loadnc(ref_file);
  ref = flatten_sza(ref_orig,isza_list);
  hr_ref = calc_hr(ref,'sw')';
  if ~exist('scenario_title','var')
    scenario_title = attr_ref.global.scenario;
  end

  for ie = 1:length(ckd_file)
    ckd_tmp = loadnc(ckd_file{ie});
    if do_fix_ssi(ie)
      ckd_tmp = correct_ssi(ref_orig, ckd_tmp, specdef);
    end
    % The Makefile in this directory uses nco tools to concatenate the
    % 5 calculations each with a different solar zenith angle, but
    % this leads to the 50 columns and 5 solar zenith angles being
    % interlaced differently than in the reference dataset. This
    % function permutes them to the same convention.
    ckd{ie} = permute_sza(ckd_tmp,[0.1 0.3 0.5 0.7 0.9]);
    hr_ckd{ie} = calc_hr(ckd{ie},'sw')';
  end

  % Pressure in hPa at full levels and half levels
  p_fl_hPa = 0.01.*0.5.*(ref.pressure_hl(1:end-1,:) + ref.pressure_hl(2:end,:));
  p_hl_hPa = 0.01.*ref.pressure_hl;

  iprof = 1:size(ref.pressure_hl,2);
  iprof_orig = 1:size(ref_orig.pressure_hl,2);

  xloc = 0.025;yloc = 0.95;

  clf
  set(gcf,'paperposition',[0.5 0.5 27 25]);
  set(gcf,'defaultaxeslayer','top','defaultlinelinewidth',0.5,'defaulttextfontsize',10);

  % Plot reference flux and heating-rate profiles
  subplot(3,3,1)
  %semilogy(ref.flux_up_sw(:,iprof), p_hl_hPa(:,iprof), col_ref);
  for isza = 1:nsza
    semilogy(squeeze(ref_orig.flux_up_sw(:,isza,iprof_orig)), p_hl_hPa(:,iprof_orig), ref_cols{isza});
%...
%	     ref_styles{isza}, 'color', ref_shades{isza});
    hold on
  end
  set(gca,'ydir','reverse');
  ylim([0.01 1000]);
  ylabel('Pressure (hPa)');
  xlabel('Upwelling shortwave flux (W m^{-2})');
  title('Reference profiles')
  text(xloc, yloc, '\bf(a)','units','normalized')

  subplot(3,3,4)
  for isza = 1:nsza
    semilogy(-1,-1, ref_cols{isza});
    hold on
  end
  %semilogy(ref.flux_dn_sw(:,iprof), p_hl_hPa, col_ref);
  for isza = 1:nsza
    semilogy(squeeze(ref_orig.flux_dn_sw(:,isza,iprof_orig)), p_hl_hPa(:,iprof_orig), ref_cols{isza});
    hold on
  end
  set(gca,'ydir','reverse');
  ylim([0.01 1000]);
  ylabel('Pressure (hPa)');
  xlabel('Downwelling shortwave flux (W m^{-2})');
  legend('\mu_0 = 0.1',...
	 '\mu_0 = 0.3',...
	 '\mu_0 = 0.5',...
	 '\mu_0 = 0.7',...
	 '\mu_0 = 0.9','location','northeast');
  text(xloc, yloc, '\bf(d)','units','normalized')

  subplot(3,3,7)
  %semilogy(hr_ref(:,iprof), p_fl_hPa(:,iprof), col_ref);
  for isza = 1:nsza
    semilogy(-1,-1, ref_cols{isza});
    hold on
  end
  for isza = 1:nsza
    semilogy(hr_ref(:,isza:nsza:end), p_fl_hPa(:,iprof_orig), ref_cols{isza});
    hold on
  end
  text(xloc, yloc, '\bf(g)','units','normalized')
  if 0
     % Extra legend...
    legend('\mu_0 = 0.1',...
	   '\mu_0 = 0.3',...
	   '\mu_0 = 0.5',...
	   '\mu_0 = 0.7',...
	   '\mu_0 = 0.9','location','southeast');
  end

  set(gca,'ydir','reverse');
  ylim([0.01 1000]);
  ylabel('Pressure (hPa)');
  xlabel('Heating rate (K d^{-1})');
  xlim([0 40]);

  % Plot random errors and biases in flux and heating-rate profiles
  subplot(3,3,8)

  ileg = 1;
  leg = {};
  for ie = ickd
    err = calc_hr_error(p_hl_hPa(:,iprof), hr_ckd{ie}(:,iprof),hr_ref(:,iprof), [0.02 4]);
    err4 = calc_hr_error(p_hl_hPa(:,iprof), hr_ckd{ie}(:,iprof),hr_ref(:,iprof), [4 1100]);
    leg{ileg} = [titles{ie} ' (RMSE=' num2str(err,'%0.3f') ' K d^{-1})'];
    plot(-1,-1,cols{ie},'linewidth',1.5);
    ileg = ileg+1;
    rmse_save(1,ie) = err;
    rmse_save(2,ie) = err4;
    hold on
  end

  for ie = ickd
    hr_error = hr_scaling_factor .* (hr_ckd{ie}(:,iprof) - hr_ref(:,iprof));
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
  text(xloc, yloc, '\bf(h)','units','normalized')

  subplot(3,3,2)
  for ie = ickd
    plot(-1,-1,cols{ie},'linewidth',1.5);
    hold on
  end
  for ie = ickd
    flux_error = ckd{ie}.flux_up_sw(:,iprof) - ref.flux_up_sw(:,iprof);
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
    set(legend(titles{ickd},'location','northeast'),'fontsize',legfontsize','box','on');
  end
  text(xloc, yloc, '\bf(b)','units','normalized')  

  subplot(3,3,5)
  for ie = ickd
    flux_error = ckd{ie}.flux_dn_sw(:,iprof) - ref.flux_dn_sw(:,iprof);
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
  text(xloc, yloc, '\bf(e)','units','normalized')

  subplot(3,3,3)
  ileg = 1;
  leg = {};
  symsurf = {'+','o'};
  for ie = ickd
    mydiff = ckd{ie}.flux_up_sw(1,iprof)-ref.flux_up_sw(1,iprof);
    err = sqrt(mean(mydiff.^2));
    leg{ileg} = [titles{ie} ' (RMSE=' num2str(err,'%0.2f') ' W m^{-2})'];
    rmse_save(3,ie) = err;
    rmse_save(4,ie) = mean(mydiff);
    %plot(-1,-1,cols{ie},'linewidth',1.5);
    ileg = ileg+1;
    plot(ref.flux_up_sw(1,iprof), mydiff, ...
	 [cols{ie}(1) symsurf{ie}]);
    hold on
  end
  %if do_legend
  %  set(legend(leg,'location','best'),'fontsize',legfontsize','box','on')
  %end
  xlabel('Reference TOA upwelling (W m^{-2})');
  ylabel('TOA upwelling error (W m^{-2})');
  xlim([0 250]);
  plot(xlim,[0 0],'k:','linewidth',0.5);
  ylim(dflux_axis);
  title('Errors at surface and TOA');
  text(xloc, yloc, '\bf(c)','units','normalized')

  subplot(3,3,6)
  ileg = 1;
  leg = {};
  for ie = ickd
    mydiff = ckd{ie}.flux_dn_sw(end,iprof)-ref.flux_dn_sw(end,iprof);
    err = sqrt(mean(mydiff.^2));
    leg{ileg} = [titles{ie} ' (RMSE=' num2str(err,'%0.2f') ' W m^{-2})'];
    rmse_save(5,ie) = err;
    rmse_save(6,ie) = mean(mydiff);
    ileg = ileg+1;
    plot(ref.flux_dn_sw(end,iprof), mydiff, ...
	 [cols{ie}(1) symsurf{ie}]);
    hold on
  end
  %if do_legend
  %  set(legend(leg,'location','best'),'fontsize',legfontsize')
  %end
  xlabel('Reference surface downwelling (W m^{-2})');
  ylabel('Surface downwelling error (W m^{-2})');
  plot(xlim,[0 0],'k:','linewidth',0.5);
  ylim(dflux_axis);
  xlim([0 1200])
  text(xloc, yloc, '\bf(f)','units','normalized')

  h=subplot(3,3,9);
  set(h,'visible','off')

  xstart = -0.15;
  text(xstart,0.95,['\bfScenario: ' scenario_title],'fontsize',12);

  if ickd_stats > 0
  for ie = ickd_stats
    text(xstart,0.85,['\bfCKD model: ' titles{ie}],'fontsize',12);
    ystart = 0.75; %-0.4.*(ie-1);
    %text(xstart,ystart,['\bfRMS errors in ' titles{ie} ':'],'units','normalized');
    text(xstart,ystart,['Bias TOA upwelling: ' num2str(rmse_save(4,ie),'%0.2f') ' W m^{-2}'],'units','normalized');
    text(xstart,ystart-0.1,['Bias surface downwelling: ' num2str(rmse_save(6,ie),'%0.2f') ' W m^{-2}'],'units','normalized');
    text(xstart,ystart-0.2,['RMSE TOA upwelling: ' num2str(rmse_save(3,ie),'%0.2f') ' W m^{-2}'],'units','normalized');
    text(xstart,ystart-0.3,['RMSE surface downwelling: ' num2str(rmse_save(5,ie),'%0.2f') ' W m^{-2}'],'units','normalized');
    text(xstart,ystart-0.4,['RMSE heating rate (0.02-4 hPa):  ' num2str(rmse_save(1,ie),'%0.3f') ' K d^{-1}'],'units','normalized');
    text(xstart,ystart-0.5,['RMSE heating rate (4-1100 hPa):  ' num2str(rmse_save(2,ie),'%0.3f') ' K d^{-1}'],'units','normalized');

  end
  end

  drawnow
%  print('-dpng','-painters','-r100',[combined_titles 'evaluation1_fluxes_' scenario '.png']);




