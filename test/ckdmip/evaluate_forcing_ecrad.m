% Evaluate the forcing calculated by ecRad for the CKDMIP scenarios
% against line-by-line

DOMAIN='sw';
DATASET='evaluation1';
VAR=['flux_up_' DOMAIN];
REFDIR='lbl_fluxes';
EXPDIR='fluxes';
MODELS={'rrtmg','ecckd'};
LEGEND={'LBLRTM','RRTMG','ecCKD'};

gases={'co2','ch4','n2o'};
labels={'CO2 concentration (ppmv)','CH4 concentration (ppbv)','N2O concentration (ppbv)'}
concs={[180 280 415 560 1120 2240],[350 700 1200 1921 2600 3500],[190 270 332 405 540]};
ipresent = [3 4 3];
if ~exist('ref','var')
  for igas = 1:length(gases)
    conc = concs{igas};
    for iconc = 1:length(conc)
      if iconc == ipresent(igas)
	scenario_str = 'present';
      else
	scenario_str = [gases{igas} '-' num2str(conc(iconc))];
      end
      ref{igas}{iconc} = loadnc([REFDIR '/ckdmip_' DATASET '_' DOMAIN '_fluxes_' scenario_str '.nc']);
      data = ref{igas}{iconc}.(VAR)(1,:,:);
      ref_mean_toa_up{igas}(iconc) = mean(data(:));
    end
  end
end

if ~exist('expt','var')
  for igas = 1:length(gases)
    conc = concs{igas};
    for iconc = 1:length(conc)
      if iconc == ipresent(igas)
	scenario_str = 'present';
      else
	scenario_str = [gases{igas} '-' num2str(conc(iconc))];
      end
      for iexpt = 1:length(MODELS)
	expt{iexpt,igas}{iconc} = loadnc([EXPDIR '/ecrad-' MODELS{iexpt} '_' DATASET '_' DOMAIN '_fluxes_' scenario_str '.nc']);
	data = expt{iexpt,igas}{iconc}.(VAR)(1,:,:);
	expt_mean_toa_up{iexpt,igas}(iconc) = mean(data(:));
      end
    end
  end
end

clf
set(gcf,'defaultlinelinewidth',1,'paperposition',[0.5 0.5 30 10]);
stys = {'b--','r-.'};
for igas = 1:length(gases)
  subplot(1,length(gases),igas)
  conc = concs{igas};
  plot(conc, -ref_mean_toa_up{igas}+ref_mean_toa_up{igas}(ipresent(igas)),'k');
  hold on
  for iexpt = 1:length(MODELS)
    plot(conc, -expt_mean_toa_up{iexpt,igas}+expt_mean_toa_up{iexpt,igas}(ipresent(igas)),stys{iexpt});
  end
  if igas == 1
    set(gca,'xscale','log');
    legend(LEGEND,'location','northwest');
  end
  xlabel(labels{igas})
  ylabel('Instantaneous radiative forcing (W/m2)')
  grid on
end

