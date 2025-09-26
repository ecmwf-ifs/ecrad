! radiation_gauss_lognormal.F90 - Coefficients for Gauss Quadrature with lognormal weighting
!
! (C) Copyright 2021- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!

module radiation_gauss_lognormal

  use parkind1, only : jpim, jprb

  implicit none

  public

  ! The look-up table covers fractional standard deviations in the
  ! range 0 to 4

  ! Number of fractional standard deviations in look-up table
  integer(jpim), parameter :: nfsd = 21

  ! Spacing of fractional standard deviations in look-up table
  real(jprb), parameter :: dfsd = 0.2_jprb

  ! Weights and nodes for each value of fractional standard deviation,
  ! and each order from 2 to 7, assuming a gamma distribution
  real(jprb), parameter :: weights2(2,nfsd) = reshape([ &
   &          0.5_jprb,          0.5_jprb, &
   &          0.5_jprb,          0.5_jprb, &
   &          0.5_jprb,          0.5_jprb, &
   &          0.5_jprb,          0.5_jprb, &
   &          0.5_jprb,          0.5_jprb, &
   &          0.5_jprb,          0.5_jprb, &
   &          0.5_jprb,          0.5_jprb, &
   &          0.5_jprb,          0.5_jprb, &
   &          0.5_jprb,          0.5_jprb, &
   &          0.5_jprb,          0.5_jprb, &
   &          0.5_jprb,          0.5_jprb, &
   &          0.5_jprb,          0.5_jprb, &
   &          0.5_jprb,          0.5_jprb, &
   &          0.5_jprb,          0.5_jprb, &
   &          0.5_jprb,          0.5_jprb, &
   &          0.5_jprb,          0.5_jprb, &
   &          0.5_jprb,          0.5_jprb, &
   &          0.5_jprb,          0.5_jprb, &
   &          0.5_jprb,          0.5_jprb, &
   &          0.5_jprb,          0.5_jprb, &
   &          0.5_jprb,          0.5_jprb  &
   & ], [2,nfsd])
  real(jprb), parameter :: nodes2(2,nfsd) = reshape([ &
   &    1.0000000_jprb,    1.0000000_jprb, &
   &    0.8044049_jprb,    1.1953414_jprb, &
   &    0.6316226_jprb,    1.3648481_jprb, &
   &    0.4925027_jprb,    1.4929749_jprb, &
   &    0.3864725_jprb,    1.5777477_jprb, &
   &    0.3075467_jprb,    1.6257691_jprb, &
   &    0.2489619_jprb,    1.6461800_jprb, &
   &    0.2050874_jprb,    1.6472868_jprb, &
   &    0.1717498_jprb,    1.6355126_jprb, &
   &    0.1459952_jprb,    1.6154575_jprb, &
   &    0.1257631_jprb,    1.5902914_jprb, &
   &    0.1096132_jprb,    1.5621550_jprb, &
   &    0.0965293_jprb,    1.5324781_jprb, &
   &    0.0857844_jprb,    1.5022072_jprb, &
   &    0.0768512_jprb,    1.4719637_jprb, &
   &    0.0693409_jprb,    1.4421493_jprb, &
   &    0.0629631_jprb,    1.4130182_jprb, &
   &    0.0574972_jprb,    1.3847247_jprb, &
   &    0.0527741_jprb,    1.3573559_jprb, &
   &    0.0486620_jprb,    1.3309536_jprb, &
   &    0.0450572_jprb,    1.3055292_jprb  &
   & ], [2,nfsd])
  real(jprb), parameter :: weights3(3,nfsd) = reshape([ &
   &    0.1666667_jprb,    0.6666667_jprb,    0.1666667_jprb, &
   &    0.1666667_jprb,    0.6666667_jprb,    0.1666667_jprb, &
   &    0.1666667_jprb,    0.6666667_jprb,    0.1666667_jprb, &
   &    0.1666667_jprb,    0.6666667_jprb,    0.1666667_jprb, &
   &    0.1666667_jprb,    0.6666667_jprb,    0.1666667_jprb, &
   &    0.1666667_jprb,    0.6666667_jprb,    0.1666667_jprb, &
   &    0.1666667_jprb,    0.6666667_jprb,    0.1666667_jprb, &
   &    0.1666667_jprb,    0.6666667_jprb,    0.1666667_jprb, &
   &    0.1666667_jprb,    0.6666667_jprb,    0.1666667_jprb, &
   &    0.1666667_jprb,    0.6666667_jprb,    0.1666667_jprb, &
   &    0.1666667_jprb,    0.6666667_jprb,    0.1666667_jprb, &
   &    0.1666667_jprb,    0.6666667_jprb,    0.1666667_jprb, &
   &    0.1666667_jprb,    0.6666667_jprb,    0.1666667_jprb, &
   &    0.1666667_jprb,    0.6666667_jprb,    0.1666667_jprb, &
   &    0.1666667_jprb,    0.6666667_jprb,    0.1666667_jprb, &
   &    0.1666667_jprb,    0.6666667_jprb,    0.1666667_jprb, &
   &    0.1666667_jprb,    0.6666667_jprb,    0.1666667_jprb, &
   &    0.1666667_jprb,    0.6666667_jprb,    0.1666667_jprb, &
   &    0.1666667_jprb,    0.6666667_jprb,    0.1666667_jprb, &
   &    0.1666667_jprb,    0.6666667_jprb,    0.1666667_jprb, &
   &    0.1666667_jprb,    0.6666667_jprb,    0.1666667_jprb  &
   & ], [3,nfsd])
  real(jprb), parameter :: nodes3(3,nfsd) = reshape([ &
   &    1.0000000_jprb,    1.0000000_jprb,    1.0000000_jprb, &
   &    0.6958442_jprb,    0.9805807_jprb,    1.3818301_jprb, &
   &    0.4764045_jprb,    0.9284767_jprb,    1.8095317_jprb, &
   &    0.3281820_jprb,    0.8574929_jprb,    2.2405075_jprb, &
   &    0.2309440_jprb,    0.7808688_jprb,    2.6402766_jprb, &
   &    0.1671940_jprb,    0.7071068_jprb,    2.9905380_jprb, &
   &    0.1246999_jprb,    0.6401844_jprb,    3.2865802_jprb, &
   &    0.0956639_jprb,    0.5812382_jprb,    3.5315066_jprb, &
   &    0.0752742_jprb,    0.5299989_jprb,    3.7316774_jprb, &
   &    0.0605653_jprb,    0.4856429_jprb,    3.8941288_jprb, &
   &    0.0496845_jprb,    0.4472136_jprb,    4.0253967_jprb, &
   &    0.0414496_jprb,    0.4138029_jprb,    4.1311087_jprb, &
   &    0.0350881_jprb,    0.3846154_jprb,    4.2159346_jprb, &
   &    0.0300830_jprb,    0.3589791_jprb,    4.2836805_jprb, &
   &    0.0260805_jprb,    0.3363364_jprb,    4.3374240_jprb, &
   &    0.0228329_jprb,    0.3162278_jprb,    4.3796479_jprb, &
   &    0.0201634_jprb,    0.2982750_jprb,    4.4123576_jprb, &
   &    0.0179434_jprb,    0.2821663_jprb,    4.4371774_jprb, &
   &    0.0160777_jprb,    0.2676439_jprb,    4.4554283_jprb, &
   &    0.0144951_jprb,    0.2544933_jprb,    4.4681898_jprb, &
   &    0.0131410_jprb,    0.2425356_jprb,    4.4763480_jprb  &
   & ], [3,nfsd])
  real(jprb), parameter :: weights4(4,nfsd) = reshape([ &
   &   0.04587585_jprb,    0.4541241_jprb,    0.4541241_jprb,   0.04587585_jprb, &
   &   0.04587585_jprb,    0.4541241_jprb,    0.4541241_jprb,   0.04587585_jprb, &
   &   0.04587585_jprb,    0.4541241_jprb,    0.4541241_jprb,   0.04587585_jprb, &
   &   0.04587585_jprb,    0.4541241_jprb,    0.4541241_jprb,   0.04587585_jprb, &
   &   0.04587585_jprb,    0.4541241_jprb,    0.4541241_jprb,   0.04587585_jprb, &
   &   0.04587585_jprb,    0.4541241_jprb,    0.4541241_jprb,   0.04587585_jprb, &
   &   0.04587585_jprb,    0.4541241_jprb,    0.4541241_jprb,   0.04587585_jprb, &
   &   0.04587585_jprb,    0.4541241_jprb,    0.4541241_jprb,   0.04587585_jprb, &
   &   0.04587585_jprb,    0.4541241_jprb,    0.4541241_jprb,   0.04587585_jprb, &
   &   0.04587585_jprb,    0.4541241_jprb,    0.4541241_jprb,   0.04587585_jprb, &
   &   0.04587585_jprb,    0.4541241_jprb,    0.4541241_jprb,   0.04587585_jprb, &
   &   0.04587585_jprb,    0.4541241_jprb,    0.4541241_jprb,   0.04587585_jprb, &
   &   0.04587585_jprb,    0.4541241_jprb,    0.4541241_jprb,   0.04587585_jprb, &
   &   0.04587585_jprb,    0.4541241_jprb,    0.4541241_jprb,   0.04587585_jprb, &
   &   0.04587585_jprb,    0.4541241_jprb,    0.4541241_jprb,   0.04587585_jprb, &
   &   0.04587585_jprb,    0.4541241_jprb,    0.4541241_jprb,   0.04587585_jprb, &
   &   0.04587585_jprb,    0.4541241_jprb,    0.4541241_jprb,   0.04587585_jprb, &
   &   0.04587585_jprb,    0.4541241_jprb,    0.4541241_jprb,   0.04587585_jprb, &
   &   0.04587585_jprb,    0.4541241_jprb,    0.4541241_jprb,   0.04587585_jprb, &
   &   0.04587585_jprb,    0.4541241_jprb,    0.4541241_jprb,   0.04587585_jprb, &
   &   0.04587585_jprb,    0.4541241_jprb,    0.4541241_jprb,   0.04587585_jprb  &
   & ], [4,nfsd])
  real(jprb), parameter :: nodes4(4,nfsd) = reshape([ &
   &    1.0000000_jprb,    1.0000000_jprb,    1.0000000_jprb,    1.0000000_jprb, &
   &    0.6175947_jprb,    0.8465801_jprb,    1.1357915_jprb,    1.5569086_jprb, &
   &    0.3777395_jprb,    0.6976387_jprb,    1.2356954_jprb,    2.2821785_jprb, &
   &    0.2349916_jprb,    0.5682630_jprb,    1.2939328_jprb,    3.1290228_jprb, &
   &    0.1511851_jprb,    0.4633810_jprb,    1.3158850_jprb,    4.0331748_jprb, &
   &    0.1012562_jprb,    0.3812504_jprb,    1.3114741_jprb,    4.9379671_jprb, &
   &    0.0705981_jprb,    0.3176670_jprb,    1.2901435_jprb,    5.8051974_jprb, &
   &    0.0510775_jprb,    0.2683358_jprb,    1.2590116_jprb,    6.6142149_jprb, &
   &    0.0381823_jprb,    0.2297067_jprb,    1.2228591_jprb,    7.3567843_jprb, &
   &    0.0293632_jprb,    0.1990801_jprb,    1.1846940_jprb,    8.0321359_jprb, &
   &    0.0231389_jprb,    0.1744702_jprb,    1.1463273_jprb,    8.6434648_jprb, &
   &    0.0186208_jprb,    0.1544301_jprb,    1.1088050_jprb,    9.1957780_jprb, &
   &    0.0152588_jprb,    0.1379036_jprb,    1.0726985_jprb,    9.6946770_jprb, &
   &    0.0127015_jprb,    0.1241136_jprb,    1.0382905_jprb,   10.1457168_jprb, &
   &    0.0107183_jprb,    0.1124822_jprb,    1.0056899_jprb,   10.5540963_jprb, &
   &    0.0091537_jprb,    0.1025742_jprb,    0.9749037_jprb,   10.9245322_jprb, &
   &    0.0079004_jprb,    0.0940584_jprb,    0.9458806_jprb,   11.2612294_jprb, &
   &    0.0068827_jprb,    0.0866789_jprb,    0.9185377_jprb,   11.5679000_jprb, &
   &    0.0060461_jprb,    0.0802364_jprb,    0.8927770_jprb,   11.8478039_jprb, &
   &    0.0053510_jprb,    0.0745736_jprb,    0.8684951_jprb,   12.1037971_jprb, &
   &    0.0047675_jprb,    0.0695652_jprb,    0.8455889_jprb,   12.3383806_jprb  &
   & ], [4,nfsd])
  real(jprb), parameter :: weights5(5,nfsd) = reshape([ &
   &   0.01125741_jprb,    0.2220759_jprb,    0.5333333_jprb,    0.2220759_jprb,   0.01125741_jprb, &
   &   0.01125741_jprb,    0.2220759_jprb,    0.5333333_jprb,    0.2220759_jprb,   0.01125741_jprb, &
   &   0.01125741_jprb,    0.2220759_jprb,    0.5333333_jprb,    0.2220759_jprb,   0.01125741_jprb, &
   &   0.01125741_jprb,    0.2220759_jprb,    0.5333333_jprb,    0.2220759_jprb,   0.01125741_jprb, &
   &   0.01125741_jprb,    0.2220759_jprb,    0.5333333_jprb,    0.2220759_jprb,   0.01125741_jprb, &
   &   0.01125741_jprb,    0.2220759_jprb,    0.5333333_jprb,    0.2220759_jprb,   0.01125741_jprb, &
   &   0.01125741_jprb,    0.2220759_jprb,    0.5333333_jprb,    0.2220759_jprb,   0.01125741_jprb, &
   &   0.01125741_jprb,    0.2220759_jprb,    0.5333333_jprb,    0.2220759_jprb,   0.01125741_jprb, &
   &   0.01125741_jprb,    0.2220759_jprb,    0.5333333_jprb,    0.2220759_jprb,   0.01125741_jprb, &
   &   0.01125741_jprb,    0.2220759_jprb,    0.5333333_jprb,    0.2220759_jprb,   0.01125741_jprb, &
   &   0.01125741_jprb,    0.2220759_jprb,    0.5333333_jprb,    0.2220759_jprb,   0.01125741_jprb, &
   &   0.01125741_jprb,    0.2220759_jprb,    0.5333333_jprb,    0.2220759_jprb,   0.01125741_jprb, &
   &   0.01125741_jprb,    0.2220759_jprb,    0.5333333_jprb,    0.2220759_jprb,   0.01125741_jprb, &
   &   0.01125741_jprb,    0.2220759_jprb,    0.5333333_jprb,    0.2220759_jprb,   0.01125741_jprb, &
   &   0.01125741_jprb,    0.2220759_jprb,    0.5333333_jprb,    0.2220759_jprb,   0.01125741_jprb, &
   &   0.01125741_jprb,    0.2220759_jprb,    0.5333333_jprb,    0.2220759_jprb,   0.01125741_jprb, &
   &   0.01125741_jprb,    0.2220759_jprb,    0.5333333_jprb,    0.2220759_jprb,   0.01125741_jprb, &
   &   0.01125741_jprb,    0.2220759_jprb,    0.5333333_jprb,    0.2220759_jprb,   0.01125741_jprb, &
   &   0.01125741_jprb,    0.2220759_jprb,    0.5333333_jprb,    0.2220759_jprb,   0.01125741_jprb, &
   &   0.01125741_jprb,    0.2220759_jprb,    0.5333333_jprb,    0.2220759_jprb,   0.01125741_jprb, &
   &   0.01125741_jprb,    0.2220759_jprb,    0.5333333_jprb,    0.2220759_jprb,   0.01125741_jprb  &
   & ], [5,nfsd])
  real(jprb), parameter :: nodes5(5,nfsd) = reshape([ &
   &    1.0000000_jprb,    1.0000000_jprb,    1.0000000_jprb,    1.0000000_jprb,    1.0000000_jprb, &
   &    0.5568769_jprb,    0.7497005_jprb,    0.9805807_jprb,    1.2825635_jprb,    1.7266624_jprb, &
   &    0.3088601_jprb,    0.5507528_jprb,    0.9284767_jprb,    1.5652558_jprb,    2.7911306_jprb, &
   &    0.1758772_jprb,    0.4043580_jprb,    0.8574929_jprb,    1.8184237_jprb,    4.1807253_jprb, &
   &    0.1046862_jprb,    0.3009464_jprb,    0.7808688_jprb,    2.0261283_jprb,    5.8246075_jprb, &
   &    0.0655359_jprb,    0.2287315_jprb,    0.7071068_jprb,    2.1859691_jprb,    7.6294058_jprb, &
   &    0.0430978_jprb,    0.1779364_jprb,    0.6401844_jprb,    2.3032724_jprb,    9.5094384_jprb, &
   &    0.0296359_jprb,    0.1415952_jprb,    0.5812382_jprb,    2.3859418_jprb,   11.3996210_jprb, &
   &    0.0211902_jprb,    0.1150430_jprb,    0.5299989_jprb,    2.4416861_jprb,   13.2561007_jprb, &
   &    0.0156690_jprb,    0.0952159_jprb,    0.4856429_jprb,    2.4769929_jprb,   15.0519902_jprb, &
   &    0.0119243_jprb,    0.0800971_jprb,    0.4472136_jprb,    2.4969689_jprb,   16.7724346_jprb, &
   &    0.0093008_jprb,    0.0683426_jprb,    0.4138029_jprb,    2.5055077_jprb,   18.4105552_jprb, &
   &    0.0074096_jprb,    0.0590408_jprb,    0.3846154_jprb,    2.5055363_jprb,   19.9645221_jprb, &
   &    0.0060118_jprb,    0.0515620_jprb,    0.3589791_jprb,    2.4992441_jprb,   21.4355751_jprb, &
   &    0.0049557_jprb,    0.0454622_jprb,    0.3363364_jprb,    2.4882688_jprb,   22.8267305_jprb, &
   &    0.0041422_jprb,    0.0404230_jprb,    0.3162278_jprb,    2.4738395_jprb,   24.1419565_jprb, &
   &    0.0035047_jprb,    0.0362117_jprb,    0.2982750_jprb,    2.4568816_jprb,   25.3856565_jprb, &
   &    0.0029974_jprb,    0.0326557_jprb,    0.2821663_jprb,    2.4380956_jprb,   26.5623512_jprb, &
   &    0.0025882_jprb,    0.0296248_jprb,    0.2676439_jprb,    2.4180135_jprb,   27.6764891_jprb, &
   &    0.0022541_jprb,    0.0270195_jprb,    0.2544933_jprb,    2.3970416_jprb,   28.7323374_jprb, &
   &    0.0019783_jprb,    0.0247627_jprb,    0.2425356_jprb,    2.3754904_jprb,   29.7339243_jprb  &
   & ], [5,nfsd])
  real(jprb), parameter :: weights6(6,nfsd) = reshape([ &
   &  0.002555784_jprb,   0.08861575_jprb,    0.4088285_jprb,    0.4088285_jprb,   0.08861575_jprb,  0.002555784_jprb, &
   &  0.002555784_jprb,   0.08861575_jprb,    0.4088285_jprb,    0.4088285_jprb,   0.08861575_jprb,  0.002555784_jprb, &
   &  0.002555784_jprb,   0.08861575_jprb,    0.4088285_jprb,    0.4088285_jprb,   0.08861575_jprb,  0.002555784_jprb, &
   &  0.002555784_jprb,   0.08861575_jprb,    0.4088285_jprb,    0.4088285_jprb,   0.08861575_jprb,  0.002555784_jprb, &
   &  0.002555784_jprb,   0.08861575_jprb,    0.4088285_jprb,    0.4088285_jprb,   0.08861575_jprb,  0.002555784_jprb, &
   &  0.002555784_jprb,   0.08861575_jprb,    0.4088285_jprb,    0.4088285_jprb,   0.08861575_jprb,  0.002555784_jprb, &
   &  0.002555784_jprb,   0.08861575_jprb,    0.4088285_jprb,    0.4088285_jprb,   0.08861575_jprb,  0.002555784_jprb, &
   &  0.002555784_jprb,   0.08861575_jprb,    0.4088285_jprb,    0.4088285_jprb,   0.08861575_jprb,  0.002555784_jprb, &
   &  0.002555784_jprb,   0.08861575_jprb,    0.4088285_jprb,    0.4088285_jprb,   0.08861575_jprb,  0.002555784_jprb, &
   &  0.002555784_jprb,   0.08861575_jprb,    0.4088285_jprb,    0.4088285_jprb,   0.08861575_jprb,  0.002555784_jprb, &
   &  0.002555784_jprb,   0.08861575_jprb,    0.4088285_jprb,    0.4088285_jprb,   0.08861575_jprb,  0.002555784_jprb, &
   &  0.002555784_jprb,   0.08861575_jprb,    0.4088285_jprb,    0.4088285_jprb,   0.08861575_jprb,  0.002555784_jprb, &
   &  0.002555784_jprb,   0.08861575_jprb,    0.4088285_jprb,    0.4088285_jprb,   0.08861575_jprb,  0.002555784_jprb, &
   &  0.002555784_jprb,   0.08861575_jprb,    0.4088285_jprb,    0.4088285_jprb,   0.08861575_jprb,  0.002555784_jprb, &
   &  0.002555784_jprb,   0.08861575_jprb,    0.4088285_jprb,    0.4088285_jprb,   0.08861575_jprb,  0.002555784_jprb, &
   &  0.002555784_jprb,   0.08861575_jprb,    0.4088285_jprb,    0.4088285_jprb,   0.08861575_jprb,  0.002555784_jprb, &
   &  0.002555784_jprb,   0.08861575_jprb,    0.4088285_jprb,    0.4088285_jprb,   0.08861575_jprb,  0.002555784_jprb, &
   &  0.002555784_jprb,   0.08861575_jprb,    0.4088285_jprb,    0.4088285_jprb,   0.08861575_jprb,  0.002555784_jprb, &
   &  0.002555784_jprb,   0.08861575_jprb,    0.4088285_jprb,    0.4088285_jprb,   0.08861575_jprb,  0.002555784_jprb, &
   &  0.002555784_jprb,   0.08861575_jprb,    0.4088285_jprb,    0.4088285_jprb,   0.08861575_jprb,  0.002555784_jprb, &
   &  0.002555784_jprb,   0.08861575_jprb,    0.4088285_jprb,    0.4088285_jprb,   0.08861575_jprb,  0.002555784_jprb  &
   & ], [6,nfsd])
  real(jprb), parameter :: nodes6(6,nfsd) = reshape([ &
   &    1.0000000_jprb,    1.0000000_jprb,    1.0000000_jprb,    1.0000000_jprb,    1.0000000_jprb,    1.0000000_jprb, &
   &    0.5076548_jprb,    0.6745248_jprb,    0.8678431_jprb,    1.1079634_jprb,    1.4255051_jprb,    1.8940794_jprb, &
   &    0.2579755_jprb,    0.4484218_jprb,    0.7321293_jprb,    1.1774819_jprb,    1.9224511_jprb,    3.3416701_jprb, &
   &    0.1357302_jprb,    0.3007984_jprb,    0.6091357_jprb,    1.2071105_jprb,    2.4444748_jprb,    5.4173226_jprb, &
   &    0.0753619_jprb,    0.2067814_jprb,    0.5060567_jprb,    1.2049165_jprb,    2.9487963_jprb,    8.0910343_jprb, &
   &    0.0444140_jprb,    0.1466927_jprb,    0.4231555_jprb,    1.1815987_jprb,    3.4084870_jprb,   11.2576992_jprb, &
   &    0.0277196_jprb,    0.1075023_jprb,    0.3575602_jprb,    1.1462016_jprb,    3.8123484_jprb,   14.7850635_jprb, &
   &    0.0182142_jprb,    0.0812199_jprb,    0.3057362_jprb,    1.1049980_jprb,    4.1595444_jprb,   18.5480510_jprb, &
   &    0.0125157_jprb,    0.0630597_jprb,    0.2645281_jprb,    1.0618867_jprb,    4.4544879_jprb,   22.4438000_jprb, &
   &    0.0089357_jprb,    0.0501427_jprb,    0.2314246_jprb,    1.0191184_jprb,    4.7035549_jprb,   26.3941634_jprb, &
   &    0.0065914_jprb,    0.0407053_jprb,    0.2045188_jprb,    0.9779051_jprb,    4.9133600_jprb,   30.3426556_jprb, &
   &    0.0049995_jprb,    0.0336411_jprb,    0.1823882_jprb,    0.9388376_jprb,    5.0899893_jprb,   34.2498976_jprb, &
   &    0.0038837_jprb,    0.0282375_jprb,    0.1639744_jprb,    0.9021470_jprb,    5.2387378_jprb,   38.0893140_jprb, &
   &    0.0030797_jprb,    0.0240239_jprb,    0.1484865_jprb,    0.8678634_jprb,    5.3640838_jprb,   41.8436314_jprb, &
   &    0.0024861_jprb,    0.0206814_jprb,    0.1353284_jprb,    0.8359086_jprb,    5.4697617_jprb,   45.5022142_jprb, &
   &    0.0020384_jprb,    0.0179893_jprb,    0.1240464_jprb,    0.8061502_jprb,    5.5588663_jprb,   49.0591086_jprb, &
   &    0.0016943_jprb,    0.0157914_jprb,    0.1142911_jprb,    0.7784329_jprb,    5.6339572_jprb,   52.5116294_jprb, &
   &    0.0014253_jprb,    0.0139750_jprb,    0.1057909_jprb,    0.7525962_jprb,    5.6971538_jprb,   55.8593513_jprb, &
   &    0.0012120_jprb,    0.0124575_jprb,    0.0983320_jprb,    0.7284833_jprb,    5.7502148_jprb,   59.1033900_jprb, &
   &    0.0010405_jprb,    0.0111771_jprb,    0.0917447_jprb,    0.7059465_jprb,    5.7946053_jprb,   62.2458923_jprb, &
   &    0.0009010_jprb,    0.0100871_jprb,    0.0858927_jprb,    0.6848491_jprb,    5.8315501_jprb,   65.2896748_jprb  &
   & ], [6,nfsd])
  real(jprb), parameter :: weights7(7,nfsd) = reshape([ &
   & 0.000548269_jprb,   0.0307571_jprb,    0.240123_jprb,    0.457143_jprb,    0.240123_jprb,   0.0307571_jprb, 0.000548269_jprb, &
   & 0.000548269_jprb,   0.0307571_jprb,    0.240123_jprb,    0.457143_jprb,    0.240123_jprb,   0.0307571_jprb, 0.000548269_jprb, &
   & 0.000548269_jprb,   0.0307571_jprb,    0.240123_jprb,    0.457143_jprb,    0.240123_jprb,   0.0307571_jprb, 0.000548269_jprb, &
   & 0.000548269_jprb,   0.0307571_jprb,    0.240123_jprb,    0.457143_jprb,    0.240123_jprb,   0.0307571_jprb, 0.000548269_jprb, &
   & 0.000548269_jprb,   0.0307571_jprb,    0.240123_jprb,    0.457143_jprb,    0.240123_jprb,   0.0307571_jprb, 0.000548269_jprb, &
   & 0.000548269_jprb,   0.0307571_jprb,    0.240123_jprb,    0.457143_jprb,    0.240123_jprb,   0.0307571_jprb, 0.000548269_jprb, &
   & 0.000548269_jprb,   0.0307571_jprb,    0.240123_jprb,    0.457143_jprb,    0.240123_jprb,   0.0307571_jprb, 0.000548269_jprb, &
   & 0.000548269_jprb,   0.0307571_jprb,    0.240123_jprb,    0.457143_jprb,    0.240123_jprb,   0.0307571_jprb, 0.000548269_jprb, &
   & 0.000548269_jprb,   0.0307571_jprb,    0.240123_jprb,    0.457143_jprb,    0.240123_jprb,   0.0307571_jprb, 0.000548269_jprb, &
   & 0.000548269_jprb,   0.0307571_jprb,    0.240123_jprb,    0.457143_jprb,    0.240123_jprb,   0.0307571_jprb, 0.000548269_jprb, &
   & 0.000548269_jprb,   0.0307571_jprb,    0.240123_jprb,    0.457143_jprb,    0.240123_jprb,   0.0307571_jprb, 0.000548269_jprb, &
   & 0.000548269_jprb,   0.0307571_jprb,    0.240123_jprb,    0.457143_jprb,    0.240123_jprb,   0.0307571_jprb, 0.000548269_jprb, &
   & 0.000548269_jprb,   0.0307571_jprb,    0.240123_jprb,    0.457143_jprb,    0.240123_jprb,   0.0307571_jprb, 0.000548269_jprb, &
   & 0.000548269_jprb,   0.0307571_jprb,    0.240123_jprb,    0.457143_jprb,    0.240123_jprb,   0.0307571_jprb, 0.000548269_jprb, &
   & 0.000548269_jprb,   0.0307571_jprb,    0.240123_jprb,    0.457143_jprb,    0.240123_jprb,   0.0307571_jprb, 0.000548269_jprb, &
   & 0.000548269_jprb,   0.0307571_jprb,    0.240123_jprb,    0.457143_jprb,    0.240123_jprb,   0.0307571_jprb, 0.000548269_jprb, &
   & 0.000548269_jprb,   0.0307571_jprb,    0.240123_jprb,    0.457143_jprb,    0.240123_jprb,   0.0307571_jprb, 0.000548269_jprb, &
   & 0.000548269_jprb,   0.0307571_jprb,    0.240123_jprb,    0.457143_jprb,    0.240123_jprb,   0.0307571_jprb, 0.000548269_jprb, &
   & 0.000548269_jprb,   0.0307571_jprb,    0.240123_jprb,    0.457143_jprb,    0.240123_jprb,   0.0307571_jprb, 0.000548269_jprb, &
   & 0.000548269_jprb,   0.0307571_jprb,    0.240123_jprb,    0.457143_jprb,    0.240123_jprb,   0.0307571_jprb, 0.000548269_jprb, &
   & 0.000548269_jprb,   0.0307571_jprb,    0.240123_jprb,    0.457143_jprb,    0.240123_jprb,   0.0307571_jprb, 0.000548269_jprb  &
   & ], [7,nfsd])
  real(jprb), parameter :: nodes7(7,nfsd) = reshape([ &
   &    1.000000_jprb,    1.000000_jprb,    1.000000_jprb,    1.000000_jprb,    1.000000_jprb,    1.000000_jprb,    1.000000_jprb, &
   &    0.466566_jprb,    0.613651_jprb,    0.780179_jprb,    0.980581_jprb,    1.232458_jprb,    1.566914_jprb,    2.060884_jprb, &
   &    0.218913_jprb,    0.373062_jprb,    0.595146_jprb,    0.928477_jprb,    1.448500_jprb,    2.310795_jprb,    3.937944_jprb, &
   &    0.107162_jprb,    0.230814_jprb,    0.452089_jprb,    0.857493_jprb,    1.626435_jprb,    3.185651_jprb,    6.861495_jprb, &
   &    0.055843_jprb,    0.147785_jprb,    0.346700_jprb,    0.780869_jprb,    1.758742_jprb,    4.125981_jprb,   10.919065_jprb, &
   &    0.031148_jprb,    0.098566_jprb,    0.270447_jprb,    0.707107_jprb,    1.848791_jprb,    5.072749_jprb,   16.052623_jprb, &
   &    0.018534_jprb,    0.068474_jprb,    0.215179_jprb,    0.640184_jprb,    1.904628_jprb,    5.985275_jprb,   22.112169_jprb, &
   &    0.011684_jprb,    0.049385_jprb,    0.174616_jprb,    0.581238_jprb,    1.934749_jprb,    6.840877_jprb,   28.914093_jprb, &
   &    0.007743_jprb,    0.036816_jprb,    0.144322_jprb,    0.529999_jprb,    1.946331_jprb,    7.629869_jprb,   36.279476_jprb, &
   &    0.005354_jprb,    0.028244_jprb,    0.121267_jprb,    0.485643_jprb,    1.944876_jprb,    8.350540_jprb,   44.052025_jprb, &
   &    0.003839_jprb,    0.022209_jprb,    0.103391_jprb,    0.447214_jprb,    1.934408_jprb,    9.005521_jprb,   52.103118_jprb, &
   &    0.002838_jprb,    0.017838_jprb,    0.089286_jprb,    0.413803_jprb,    1.917809_jprb,    9.599518_jprb,   60.330491_jprb, &
   &    0.002155_jprb,    0.014592_jprb,    0.077976_jprb,    0.384615_jprb,    1.897118_jprb,   10.138001_jprb,   68.654529_jprb, &
   &    0.001673_jprb,    0.012127_jprb,    0.068773_jprb,    0.358979_jprb,    1.873776_jprb,   10.626505_jprb,   77.014083_jprb, &
   &    0.001325_jprb,    0.010219_jprb,    0.061187_jprb,    0.336336_jprb,    1.848800_jprb,   11.070272_jprb,   85.362651_jprb, &
   &    0.001068_jprb,    0.008715_jprb,    0.054857_jprb,    0.316228_jprb,    1.822909_jprb,   11.474101_jprb,   93.665172_jprb, &
   &    0.000873_jprb,    0.007513_jprb,    0.049520_jprb,    0.298275_jprb,    1.796608_jprb,   11.842298_jprb,  101.895452_jprb, &
   &    0.000724_jprb,    0.006537_jprb,    0.044975_jprb,    0.282166_jprb,    1.770254_jprb,   12.178687_jprb,  110.034137_jprb, &
   &    0.000607_jprb,    0.005737_jprb,    0.041072_jprb,    0.267644_jprb,    1.744097_jprb,   12.486641_jprb,  118.067144_jprb, &
   &    0.000514_jprb,    0.005072_jprb,    0.037692_jprb,    0.254493_jprb,    1.718308_jprb,   12.769127_jprb,  125.984451_jprb, &
   &    0.000440_jprb,    0.004515_jprb,    0.034745_jprb,    0.242536_jprb,    1.693007_jprb,   13.028753_jprb,  133.779162_jprb  &
   & ], [7,nfsd])

contains

  !---------------------------------------------------------------------
  ! Compute the weights and nodes for integrating across a lognormal
  ! distribution given the fractional standard deviation (=standard
  ! deviation divided by the mean) assuming a mean of unity, using
  ! Gauss-Hermite Quadrature. The user can then multiply the nodes
  ! output from this routine by the mean of the distribution.  If the
  ! quadrature order "norder" is less than 1, the weights and nodes
  ! will not be written, while if it is larger than 7, only 7 nodes
  ! will be written and the remainder will be assigned zero weight.
  subroutine calc_gauss_lognormal(ng, norder, fsd, weights, nodes)

    ! Number of independent sets of weights required (typically this
    ! corresponds to spectral interval or atmospheric column)
    integer(jpim), intent(in)  :: ng
    ! Number of nodes to be used for the quadrature
    integer(jpim), intent(in)  :: norder
    ! Fractional standard deviation
    real(jprb),    intent(in)  :: fsd(ng)

    ! Weighting of each node
    real(jprb),    intent(out) :: weights(ng,norder)
    ! Value of each node assuming a mean of one
    real(jprb),    intent(out) :: nodes(ng,norder)

    ! Index for look-up table in fractional standard deviation
    integer(jpim) :: ifsd(ng)
    ! Weights for look-up table
    real(jprb)    :: w1(ng), w2(ng)

    ! Loop index
    integer(jpim) :: jg

    ! Compute look-up table indices as real numbers
    w1 = 1.0_jprb + max(0.0_jprb, min(fsd / dfsd, nfsd-1.0001_jprb))
    ! Indices of the first point in the interpolation
    ifsd = floor(w1)
    ! Weight of the second point in the interpolation
    w2 = w1 - ifsd
    ! Weight of the first point in the interpolation
    w1 = 1.0_jprb - w2

    ! According to the order of the quadrature, we interpolate from
    ! different parameter arrays weights[2-7] and nodes[2-7]
    select case (norder)
      case (2)
        do jg = 1,ng
          weights(jg,:) = w1(jg) * weights2(:,ifsd(jg)) + w2(jg) * weights2(:,ifsd(jg)+1)
          nodes(jg,:)   = w1(jg) * nodes2(:,ifsd(jg))   + w2(jg) * nodes2(:,ifsd(jg)+1)
        end do
      case (3)
        do jg = 1,ng
          weights(jg,:) = w1(jg) * weights3(:,ifsd(jg)) + w2(jg) * weights3(:,ifsd(jg)+1)
          nodes(jg,:)   = w1(jg) * nodes3(:,ifsd(jg))   + w2(jg) * nodes3(:,ifsd(jg)+1)
        end do
      case (4)
        do jg = 1,ng
          weights(jg,:) = w1(jg) * weights4(:,ifsd(jg)) + w2(jg) * weights4(:,ifsd(jg)+1)
          nodes(jg,:)   = w1(jg) * nodes4(:,ifsd(jg))   + w2(jg) * nodes4(:,ifsd(jg)+1)
        end do
      case (5)
        do jg = 1,ng
          weights(jg,:) = w1(jg) * weights5(:,ifsd(jg)) + w2(jg) * weights5(:,ifsd(jg)+1)
          nodes(jg,:)   = w1(jg) * nodes5(:,ifsd(jg))   + w2(jg) * nodes5(:,ifsd(jg)+1)
        end do
      case (6)
        do jg = 1,ng
          weights(jg,:) = w1(jg) * weights6(:,ifsd(jg)) + w2(jg) * weights6(:,ifsd(jg)+1)
          nodes(jg,:)   = w1(jg) * nodes6(:,ifsd(jg))   + w2(jg) * nodes6(:,ifsd(jg)+1)
        end do
      case (7:)
        do jg = 1,ng
          weights(jg,1:7) = w1(jg) * weights7(:,ifsd(jg)) + w2(jg) * weights7(:,ifsd(jg)+1)
          nodes(jg,1:7)   = w1(jg) * nodes7(:,ifsd(jg))   + w2(jg) * nodes7(:,ifsd(jg)+1)
        end do
        weights(:,8:) = 0.0_jprb
        nodes(:,8:)   = 1.0_jprb
      case (1)
        weights = 1.0_jprb
        nodes   = 1.0_jprb
      case default
        ! norder <= 0: do nothing
    end select

  end subroutine calc_gauss_lognormal

end module radiation_gauss_lognormal
