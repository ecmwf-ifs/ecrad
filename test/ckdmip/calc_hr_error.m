% Compute RMS error in heating rate, in K d-1, weighting by the
% cube-root of pressure, where pressure arguments are in hPa
function err = calc_hr_error(pressure_hl, hr, hr_ref, pressure_range);

  % If final argument not provided, consider all pressure ranges
  if nargin < 4
    pressure_range = [0 Inf];
  end

  mypow = (1./3);
  
  pressure_fl = 0.5.*(pressure_hl(1:end-1,:) + pressure_hl(2:end,:));
  weight = (pressure_hl(2:end,:)).^mypow - (pressure_hl(1:end-1,:)).^mypow;

  index = find(pressure_fl < pressure_range(1) | pressure_fl >= pressure_range(2));
  if ~isempty(index)
    weight(index) = 0;
  end
  nprof = size(weight,2);
  for ii = 1:nprof
    weight(:,ii) = weight(:,ii) ./ sum(weight(:,ii));
  end
  err = sqrt(sum(weight(:) .* (hr(:)-hr_ref(:)).^2)./nprof);
