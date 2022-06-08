function hr = calc_hr(data, band, icol)

eval(['flux_up = data.flux_up_' band ';']);
eval(['flux_dn = data.flux_dn_' band ';']);
flux_net = flux_dn-flux_up;
g=9.81;
scaling = 3600.*24;
hr_all = -scaling.*(diff(flux_net).*g./diff(data.pressure_hl)./1004)';
if nargin > 2
  hr = hr_all(icol,:);
else
  hr = hr_all;
end
