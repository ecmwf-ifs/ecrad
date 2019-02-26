function [Ds] = smooth1D(D)
Ds=D;
s = length(Ds);
Ds(2:s-1) = 0.25 .* (D(1:s-2)+D(3:s)) + 0.5 .* D(2:s-1);
