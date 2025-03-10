CDF       
      band_lw       band_sw       coeff      	   	coeff_gen               comment      This file provides a parameterization of ice particle scattering in the longwave and shortwave RRTM bands,
using the data from Baran et al. (J. Climate, 2016, 29, 5299-5316). Effective radius is not used.
If ice mass mixing ratio is qi (in kg/kg), temperature is T (in K), the 1-based 9-element array of
band-specific coefficients is p (either coeff_lw or coeff_sw), the 1-based 5-element array of general
general coefficients is s (from coeff_gen), qi_mod=(qi*exp(s[1]*(T-s[2]))), then:
  mass extinction coefficient (m2/kg) = qi * (p[1] + p[2]/(1+p[3]*qi_mod^s[3])),
  single scattering albedo = p[4] + p[5]/(1+p[6]*qi_mod^s[4]), and
  asymmetry factor = p[7] + p[8]/(1+p[9]*qi^s[5]),
where mass extinction coefficient is the total extinction cross section per unit mass of cloudy air.          wavenumber1_lw                  	long_name         (Lower bound wavenumber for longwave band   units         cm-1      @  �   wavenumber2_lw                  	long_name         (Upper bound wavenumber for longwave band   units         cm-1      @     wavenumber1_sw                 	long_name         )Lower bound wavenumber for shortwave band      units         cm-1      8  D   wavenumber2_sw                 	long_name         )Upper bound wavenumber for shortwave band      units         cm-1      8  |   	coeff_gen                  	long_name         General coefficients        �   coeff_lw                   	long_name         Longwave coefficients        @  �   coeff_sw                  	long_name         Shortwave coefficients       �  
A   C�  C�  D� D/  DM  Du  D�  D�� D�� D�  D�  E  E� E� E"� C�  C�  D� D/  DM  Du  D�  D�� D�� D�  D�  E  E� E� E"� EK  E"� EK  Ez  E�P E�� E�0 E� E�� FH� Fz  F�� F� Gp DM  EK  Ez  E�P E�� E�0 E� E�� FH� Fz  F�� F� Gp GCP E"� =���Ci&f?*��>���>���BP  ��  G�` ?
=q��z�CH  ?k��
=C4  BT  �p  H/� ?�>���A�  ?p�׾��B�  BL  @   D/  ?�>8Q�B�  ?s33�W
=B�  BJ��@9��D�  ?�D<��
C\  ?p�׾L��C   BJ��?�ffD�  ?\)��C  ?n{�.{C�  BL  �P  G�@ ?
=q�aG�CH  ?u��C\  BL  �`  H/� ?�\>B�\B4  ?z�H��Q�B�  BN  ��33H/� ?�>�\)B�  ?z�H��G�B�  BJ��?���Dz  ?�>�=qB�  ?xQ��B�  BJ��?�ffE@ ?�>L��C  ?xQ��C   BJ��?ٙ�E;� ?�>.{C  ?xQ��C\  BK33@9��F;� ?�>�{B�  ?xQ�\)B�  BK33@L��G@ ?�>�  Bp  ?xQ���B�  BJ��@L��GCP ?
=q>uB�  ?xQ�\)B�  BJ��@���Gj` ?
=q>��RBp  ?xQ�#�
BH  BK33@���Gj` ?��>�=qBH  ?s33�.{B�  BK33@�  G�@ ?��>�B  ?p�׽�B�  BK33�333H/� ?��>B�\A@  ?s�ϽuA@  BK33@&ffG�� ?B�\>k�A  ?h�ý�Q�A   BK33?L��H/� ?!G�>���A�  ?p�׾\)A�  BJ��@ffH/� ?G�>W
=@�  ?^�R��Q�A�  BJ��@�33H� ?Q�>.{@�  ?\(����
A�  BJ��@s33H@ ?w��=o@�  ?O\)�49XB�  BK33?���H� ?~��;�u@�  ?MO߽Y�C\  BK{=�G�E�  ?�N9���@�  ?L�D�Y�Cp  BJ�>�E�@ ?��8Q�@�  ?K��P�`Cp  BK  >�EZ� ?�      ?�  ?J=q�D��Cp  BK{>�E�@ ?�      ?�  ?H�9�D��Cp  BK  >B�\F� ?�      ?�  ?E�˽H�9Cp  BK33@   F;� ?�>��B�  ?xQ��B�  