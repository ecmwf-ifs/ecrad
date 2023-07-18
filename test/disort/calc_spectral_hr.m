function hr = calc_spectral_hr(data, band, icol, prefix)

if nargin < 4
  prefix = 'spectral';
end

% If we select "direct_sw" for the band then heating rate is computed
% assuming a zero scattering atmosphere using the direct downwelling
% flux only
do_direct = 0;
if length(band) > 6
  if strcmp(band(1:6),'direct')
    do_direct = 1;
  end
end

eval(['flux_dn = data.' prefix '_flux_dn_' band ';']);

if ~do_direct
  eval(['flux_up = data.' prefix '_flux_up_' band ';']);
else
  flux_up = zeros(size(flux_dn));
end


%if isfield(data,[prefix '_flux_dn_direct_' band])
%  disp('Adding direct to down (assumed to be diffuse only)');
%  eval(['flux_dn = flux_dn + data.' prefix '_flux_dn_direct_' band ';']);
%end

flux_net = flux_dn-flux_up;
g=9.81;
scaling = 3600.*24;
nband = size(flux_net,1);
for iband = 1:nband
  %size(diff(squeeze(flux_net(iband,:,:))))
  %size(diff(data.pressure_hl)./1004)
  if size(flux_net,3) == 1
    hr_all(iband,:,:) = -scaling.*(diff(squeeze(flux_net(iband,:,:))').*g./diff(data.pressure_hl)./1004)';
  else
    hr_all(iband,:,:) = -scaling.*(diff(squeeze(flux_net(iband,:,:))).*g./diff(data.pressure_hl)./1004)';
  end
end

if nargin > 2
  if ~isempty(icol)
    hr = hr_all(:,icol,:);
  else
    hr = hr_all;
  end
else
  hr = hr_all;
end
