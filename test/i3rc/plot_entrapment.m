% This script plots horizontal migration distances predicted by
% SPARTACUS's "explicit entrapment" option, but it requires (1) the
% whole code to have been compiled with "make
% PRINT_ENTRAPMENT_DATA=1", and (2) "make i3rc_print_entrapment" to
% have been run in this directory. This produces the fort.10[12] files
% containing the migration distances.

xdata = load('fort.101');
fdata = load('fort.102');
in = loadnc('i3rc_mls_cumulus_sza.nc');
out= loadnc('i3rc_mls_cumulus_mono_entr_out.nc');

z=in.height_hl(1:end-1,1)./1000;
zh=in.height_hl(:,1)./1000;
[nlev,ncol] = size(in.cloud_fraction);

x_region_direct  = flip(reshape(xdata(:,3:5),nlev,ncol,3),1);
x_region_diffuse = flip(reshape(xdata(:,6:8),nlev,ncol,3),1);
n_scat_region = flip(reshape(xdata(:,9:11),nlev,ncol,3),1);

x_diffuse_scaling = sqrt(2.0 * (n_scat_region+exp(-n_scat_region)-1.0) + 1.0e-16) ./ (n_scat_region+1.0e-8);

tot_region_diffuse = sqrt(2.0).*x_region_diffuse.*x_diffuse_scaling;
tot_region_direct = sqrt(x_region_direct.^2 + (x_region_diffuse.*x_diffuse_scaling).^2);

flux_region_direct = reshape(fdata(:,3:5),nlev,ncol,3);
flux_region_dn_diffuse = reshape(fdata(:,6:8),nlev,ncol,3);

flux_direct = squeeze(sum(flux_region_direct,3));
flux_dn_diffuse = squeeze(sum(flux_region_dn_diffuse,3));

ind = find(flux_region_dn_diffuse(:,1,1) < 1.0e-8);
flux_region_dn_diffuse(ind,:,1) = NaN;

xdirect = squeeze(sum(flux_region_direct.*tot_region_direct,3) ...
		  ./sum(flux_region_direct,3));
xdiffuse = squeeze(sum(flux_region_dn_diffuse.*tot_region_diffuse,3) ...
		   ./sum(flux_region_dn_diffuse,3));

sza = acosd(in.cos_solar_zenith_angle);

isza = [1 16 31 41];
isza = [1 23 36];

xmax = 12;
zmax = 4;

figure(2)
clf
for ii = 1:length(isza)
subplot(2,2,ii)
%plot(squeeze(x_region_direct(:,isza(ii),:)),z);
%plot(squeeze(x_region_diffuse(:,isza(ii),:)),z,'--');

plot(xdirect(:,isza(ii))./1000,z,'r')
hold on
plot(xdiffuse(:,isza(ii))./1000,z,'b');
ylim([0 zmax]);
xlim([0 xmax])
xlabel('Horizontal distance travelled (km)');
ylabel('Height (km)');
title(['(' 'a'-1+ii ') SZA=' num2str(sza(isza(ii))) '\circ'])
if ii == 1
   legend('Direct','Diffuse','location','southeast')
end
grid on
end

subplot(2,2,4)
zz = 0.5.*(zh(1:end-1)+zh(2:end));
plot(in.cloud_fraction, zz,'k');
ylim([0 zmax])
xlabel('Cloud fraction');
ylabel('Height (km)');
grid on
