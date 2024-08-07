CDF       
      band_lw       band_sw       coeff               comment      CThis file provides a parameterization of ice particle scattering in the longwave and shortwave RRTM bands,
using the parameterization of Baran et al. (J. Climate, 2016, 29, 5299-5316). Effective radius is not used.
If ice mass mixing ratio is qi (in kg/kg), temperature is T (in K) and the 1-based 5-element array of
coefficients is p, then:
mass extinction coefficient (m2/kg) = qi*p[1]/T^4,
single scattering albedo = p[2] + p[3]*qi*T, and
asymmetry factor = p[4] + p[5]*qi*T,
where mass extinction coefficient is the total extinction cross section per unit mass of cloudy air.          wavenumber1_lw                  	long_name         (Lower bound wavenumber for longwave band   units         cm-1      @  �   wavenumber2_lw                  	long_name         (Upper bound wavenumber for longwave band   units         cm-1      @  �   wavenumber1_sw                 	long_name         )Lower bound wavenumber for shortwave band      units         cm-1      8  ,   wavenumber2_sw                 	long_name         )Upper bound wavenumber for shortwave band      units         cm-1      8  d   coeff_lw                   	long_name         (Longwave droplet scattering coefficients     @  �   coeff_sw                  	long_name         )Shortwave droplet scattering coefficients          �A   C�  C�  D� D/  DM  Du  D�  D�� D�� D�  D�  E  E� E� E"� C�  C�  D� D/  DM  Du  D�  D�� D�� D�  D�  E  E� E� E"� EK  E"� EK  Ez  E�P E�� E�0 E� E�� FH� Fz  F�� F� Gp DM  EK  Ez  E�P E�� E�0 E� E�� FH� Fz  F�� F� Gp GCP E"� R��?!G���  ?6��?($R"�@?E����<?U��>��R!j?I�<0�|?ba|>;��R"�?�<o4�?bn�>1&�R!j?I�<0�|?ba|>;��RF?ی<-��?o�=�v`Rlv?&ff�hs?p�e=�J�Rlv?&ff�hs?p�e=�J�RS)?���6_�?m(�>�cR /�?���6_�?nߤ>��RA�?ԕ���?hr�>4m�RA�?ԕ���?hr�>4m�RA�?ԕ���?hr�>4m�RA�?ԕ���?hr�>4m�RA�?ԕ���?hr�>4m�RA�?ԕ���?hr�>4m�R�?@  �p�|?i�^=��R�?@  �p�|?i�^=��R�?{P����+?R�>
W�R�?{P����+?R�>
W�R�?{P����+?R�>
W�R�?{P����+?R�>
W�R�?{P����+?R�>
W�R��?�ɻK)_?J)�=?�[R��?�ɻK)_?J)�=?�[R��?�r��7�?I7L=1�3R��?�r��7�?I7L=1�3R��?�  '�/�?A�7=�,R��?�  '�/�?A�7=�,R�?@  �p�|?i�^=��