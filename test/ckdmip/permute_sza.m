function out = permute_sza(in, sza)
% Permute the implicit solar-zenith-angle and column dimensions

nsza = length(sza);
ncol = size(in.pressure_hl,2)./nsza;
nz = size(in.pressure_hl,1);

vars = {'flux_dn_direct_sw','flux_dn_sw','flux_up_sw','pressure_hl'};

out = in;

for ivar = 1:length(vars)
  out.(vars{ivar}) = reshape(permute(reshape(in.(vars{ivar}),[nz ncol nsza]), [1 3 2]), [nz ncol*nsza]);
end
