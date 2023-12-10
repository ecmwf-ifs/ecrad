function hr = calc_hr(data, band, icol)

% If we select "direct_sw" for the band then heating rate is computed
% assuming a zero scattering atmosphere using the direct downwelling
% flux only
do_direct = 0
if length(band) > 6
  if strcmp(band(1:6),'direct')
    do_direct = 1;
  end
end

eval(['flux_dn = data.flux_dn_' band ';']);
if ~do_direct
  eval(['flux_up = data.flux_up_' band ';']);
else
  flux_up = zeros(size(flux_dn));
end


flux_net = flux_dn-flux_up;
g=9.81;
scaling = 3600.*24;
hr_all = -scaling.*(diff(flux_net).*g./diff(data.pressure_hl)./1004)';
if nargin > 2
  hr = hr_all(icol,:);
else
  hr = hr_all;
end
