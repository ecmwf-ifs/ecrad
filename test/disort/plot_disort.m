ids={'disort_ref','disort_2stream','disort_elsasser','ecrad'};
titles={'DISORT-Reference','DISORT-Legendre-2s','DISORT-Elsasser-2s','ecRad-Elsasser-2s'};
iexpts=[2:4];
%ids={'disort_ref','disort_4stream','disort_4stream_jac'};
%titles={'DISORT-Reference','DISORT-Legendre-4s','DISORT-Jacobi-4s'};
%iexpts=[2:3];
icont=1;
iall=[icont iexpts];
base='mls_liquid';name='MLS'
%base='tro_cirrus';name='TRO';
in=loadnc([base '.nc']);
stys = {'k','b','r--','g-.'};
do_plot = 1;

for ii = 1:length(ids)
  filename=[base '_' ids{ii} '_out.nc'];
  d{ii}=loadnc(filename);
  hr{ii}=calc_hr(d{ii},'lw');
  hr_clear{ii}=calc_hr(d{ii},'lw_clear');
  hr_s{ii}=squeeze(calc_spectral_hr(d{ii},'lw'));
  hr_clear_s{ii}=squeeze(calc_spectral_hr(d{ii},'lw_clear'));
end

pmid = 0.5.*(d{icont}.pressure_hl(2:end)+d{icont}.pressure_hl(1:end-1))./100;

paxis = [0.1 max(d{icont}.pressure_hl)./100];
ptick=10.^[-2:3];

figure(1)
set(gcf,'paperposition',[0.5 0.5 25 15],'defaultlinelinewidth',1,'defaultaxeslayer','top');
clf
subplot(2,3,1)
for ii = 1:length(iall)
  semilogy(d{iall(ii)}.flux_dn_lw,...
	   d{iall(ii)}.pressure_hl./100,stys{iall(ii)});
  hold on
  set(gca,'ydir','reverse');
end
ylim(paxis);
set(gca,'ytick',ptick);

legend(titles{iall},'location','northeast')
xlabel(['Flux down (' Wm2 ')']);
ylabel('Pressure (hPa)');
title([name ' cloudy'])

subplot(2,3,2)
for ii = 1:length(iall)
  semilogy(d{iall(ii)}.flux_up_lw,...
	   d{iall(ii)}.pressure_hl./100,stys{iall(ii)});
  hold on
  set(gca,'ydir','reverse');
end
ylim(paxis);
set(gca,'ytick',ptick);
xlabel(['Flux up (' Wm2 ')']);

subplot(2,3,3)
for ii = 1:length(iall)
  semilogy(hr{iall(ii)},...
	   pmid,stys{iall(ii)});
  hold on
  set(gca,'ydir','reverse');
end
ylim(paxis);
set(gca,'ytick',ptick);
xlabel('Heating rate (K/day)')

subplot(2,3,4)
for ii = 1:length(iexpts)
  semilogy(d{iexpts(ii)}.flux_dn_lw-d{icont}.flux_dn_lw,...
	   d{iexpts(ii)}.pressure_hl./100,stys{iexpts(ii)});
  hold on
  set(gca,'ydir','reverse');
end
plot([0 0],paxis,'k:','linewidth',0.5);
ylim(paxis);
set(gca,'ytick',ptick);
xlabel(['Error in flux down (' Wm2 ')']);
ylabel('Pressure (hPa)');

subplot(2,3,5)
for ii = 1:length(iexpts)
  semilogy(d{iexpts(ii)}.flux_up_lw-d{icont}.flux_up_lw,...
	   d{iexpts(ii)}.pressure_hl./100,stys{iexpts(ii)});
  hold on
  set(gca,'ydir','reverse');
end
plot([0 0],paxis,'k:','linewidth',0.5);
ylim(paxis);
set(gca,'ytick',ptick);
xlabel(['Error in flux up (' Wm2 ')']);

subplot(2,3,6)
for ii = 1:length(iexpts)
  semilogy(hr{iexpts(ii)}-hr{icont},...
	   pmid,stys{iexpts(ii)});
  hold on
  set(gca,'ydir','reverse');
end
plot([0 0],paxis,'k:','linewidth',0.5);
ylim(paxis);
set(gca,'ytick',ptick);
xlabel('Error in heating rate (K/day)');
if do_plot
  print_pdf([base '_cloudy.pdf'],'-painters')
end

figure(2)
set(gcf,'paperposition',[0.5 0.5 25 15],'defaultlinelinewidth',1,'defaultaxeslayer','top');

clf
subplot(2,3,1)
for ii = 1:length(iall)
  semilogy(d{iall(ii)}.flux_dn_lw_clear,...
	   d{iall(ii)}.pressure_hl./100,stys{iall(ii)});
  hold on
  set(gca,'ydir','reverse');
end
ylim(paxis);
set(gca,'ytick',ptick);
legend(titles{iall},'location','northeast')
xlabel(['Flux down (' Wm2 ')']);
ylabel('Pressure (hPa)');
title([name ' clear'])

subplot(2,3,2)
for ii = 1:length(iall)
  semilogy(d{iall(ii)}.flux_up_lw_clear,...
	   d{iall(ii)}.pressure_hl./100,stys{iall(ii)});
  hold on
  set(gca,'ydir','reverse');
end
ylim(paxis);
set(gca,'ytick',ptick);
xlabel(['Flux up (' Wm2 ')']);

subplot(2,3,3)
for ii = 1:length(iall)
  semilogy(hr_clear{iall(ii)},...
	   pmid,stys{iall(ii)});
  hold on
  set(gca,'ydir','reverse');
end
ylim(paxis);
set(gca,'ytick',ptick);
xlabel('Heating rate (K/day)')

subplot(2,3,4)
for ii = 1:length(iexpts)
  semilogy(d{iexpts(ii)}.flux_dn_lw_clear-d{icont}.flux_dn_lw_clear,...
	   d{iexpts(ii)}.pressure_hl./100,stys{iexpts(ii)});
  hold on
  set(gca,'ydir','reverse');
end
plot([0 0],paxis,'k:','linewidth',0.5);
ylim(paxis);
set(gca,'ytick',ptick);
xlabel(['Error in flux down (' Wm2 ')']);
ylabel('Pressure (hPa)');

subplot(2,3,5)
for ii = 1:length(iexpts)
  semilogy(d{iexpts(ii)}.flux_up_lw_clear-d{icont}.flux_up_lw_clear,...
	   d{iexpts(ii)}.pressure_hl./100,stys{iexpts(ii)});
  hold on
  set(gca,'ydir','reverse');
end
plot([0 0],paxis,'k:','linewidth',0.5);
ylim(paxis);
set(gca,'ytick',ptick);
xlabel(['Error in flux up (' Wm2 ')']);

subplot(2,3,6)
for ii = 1:length(iexpts)
  semilogy(hr_clear{iexpts(ii)}-hr_clear{icont},...
	   pmid,stys{iexpts(ii)});
  hold on
  set(gca,'ydir','reverse');
end
plot([0 0],paxis,'k:','linewidth',0.5);
ylim(paxis);
set(gca,'ytick',ptick);
xlabel('Error in heating rate (K/day)');
if do_plot
  print_pdf([base '_clear.pdf'],'-painters')
end

figure(4)
set(gcf,'paperposition',[0.5 0.5 25 15],'defaultlinelinewidth',1,'defaultaxeslayer','top');
clf
subplot(2,3,1)
% Ensure legend is correct
for ii = 1:length(iall)
  semilogy(-1,-1,stys{iall(ii)})
  hold on
end
for ii = 1:length(iall)
  semilogy(d{iall(ii)}.spectral_flux_dn_lw_clear,...
	   d{iall(ii)}.pressure_hl./100,stys{iall(ii)});
  hold on
  set(gca,'ydir','reverse');
end
ylim(paxis);
set(gca,'ytick',ptick);
legend(titles{iall},'location','northeast')
xlabel(['Flux down (' Wm2 ')']);
ylabel('Pressure (hPa)');
title([name ' clear (spectral)'])

subplot(2,3,2)
for ii = 1:length(iall)
  semilogy(d{iall(ii)}.spectral_flux_up_lw_clear,...
	   d{iall(ii)}.pressure_hl./100,stys{iall(ii)});
  hold on
  set(gca,'ydir','reverse');
end
ylim(paxis);
set(gca,'ytick',ptick);
xlabel(['Flux up (' Wm2 ')']);

subplot(2,3,3)
for ii = 1:length(iall)
  semilogy(hr_clear_s{iall(ii)},...
	   pmid,stys{iall(ii)});
  hold on
  set(gca,'ydir','reverse');
end
ylim(paxis);
set(gca,'ytick',ptick);
xlabel('Heating rate (K/day)')

subplot(2,3,4)
for ii = 1:length(iexpts)
  semilogy(d{iexpts(ii)}.spectral_flux_dn_lw_clear-d{icont}.spectral_flux_dn_lw_clear,...
	   d{iexpts(ii)}.pressure_hl./100,stys{iexpts(ii)});
  hold on
  set(gca,'ydir','reverse');
end
plot([0 0],paxis,'k:','linewidth',0.5);
ylim(paxis);
set(gca,'ytick',ptick);

xlabel(['Error in flux down (' Wm2 ')']);
ylabel('Pressure (hPa)');

subplot(2,3,5)
for ii = 1:length(iexpts)
  semilogy(d{iexpts(ii)}.spectral_flux_up_lw_clear-d{icont}.spectral_flux_up_lw_clear,...
	   d{iexpts(ii)}.pressure_hl./100,stys{iexpts(ii)});
  hold on
  set(gca,'ydir','reverse');
end
plot([0 0],paxis,'k:','linewidth',0.5);
ylim(paxis);
set(gca,'ytick',ptick);
xlabel(['Error in flux up (' Wm2 ')']);

subplot(2,3,6)
for ii = 1:length(iexpts)
  semilogy(hr_clear_s{iexpts(ii)}-hr_clear_s{icont},...
	   pmid,stys{iexpts(ii)});
  hold on
  set(gca,'ydir','reverse');
end
plot([0 0],paxis,'k:','linewidth',0.5);
ylim(paxis);
set(gca,'ytick',ptick);
xlabel('Error in heating rate (K/day)');
if do_plot
  print_pdf([base '_clear_spectral.pdf'],'-painters')
end
