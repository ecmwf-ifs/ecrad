function out = flatten_sza(in, isza)
% Flatten the solar-zenith-angle and column dimensions

if nargin < 2
  isza = 1:length(in.mu0);
end
nsza = length(isza);

ncol = size(in.pressure_hl,2);
nz = size(in.pressure_hl,1);

out.pressure_hl = in.pressure_hl(:,ceil([1:nsza*ncol]./nsza));
if isfield(out,'temperature_hl')
  out.temperature_hl = in.temperature_hl(:,ceil([1:nsza*ncol]./nsza));
end
out.mu0 = repmat(in.mu0(isza), ncol, 1);
%out.reference_surface_mole_fraction = in.reference_surface_mole_fraction;
%out.mole_fraction_fl = in.mole_fraction_fl(:,:,ceil([1:nsza*ncol]./nsza));
out.flux_up_sw = reshape(in.flux_up_sw(:,isza,:),[nz nsza*ncol]);
out.flux_dn_sw = reshape(in.flux_dn_sw(:,isza,:),[nz nsza*ncol]);
out.flux_dn_direct_sw = reshape(in.flux_dn_direct_sw(:,isza,:),[nz nsza*ncol]);
if isfield(in, 'band_flux_up_sw')
  nband = size(in.band_flux_up_sw,1);
  out.band_wavenumber1_sw = in.band_wavenumber1_sw;
  out.band_wavenumber2_sw = in.band_wavenumber2_sw;
  out.band_flux_up_sw = reshape(in.band_flux_up_sw(:,:,isza,:),[nband nz nsza*ncol]);
  out.band_flux_dn_sw = reshape(in.band_flux_dn_sw(:,:,isza,:),[nband nz nsza*ncol]);
  out.band_flux_dn_direct_sw = reshape(in.band_flux_dn_direct_sw(:,:,isza,:),[nband nz nsza*ncol]);
end
if isfield(in, 'spectral_flux_up_sw')
  nspectral = size(in.spectral_flux_up_sw,1);
  out.spectral_flux_up_sw = reshape(in.spectral_flux_up_sw(:,:,isza,:),[nspectral nz nsza*ncol]);
  out.spectral_flux_dn_sw = reshape(in.spectral_flux_dn_sw(:,:,isza,:),[nspectral nz nsza*ncol]);
  out.spectral_flux_dn_direct_sw = reshape(in.spectral_flux_dn_direct_sw(:,:,isza,:),[nspectral nz nsza*ncol]);
end


