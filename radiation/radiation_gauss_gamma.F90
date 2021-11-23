! radiation_gauss_gamma.F90 - Coefficients for Gauss Quadrature with gamma weighting
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

module radiation_gauss_gamma

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
   &    0.8011502_jprb,    1.1988813_jprb, &
   &    0.6108430_jprb,    1.3897531_jprb, &
   &    0.4372079_jprb,    1.5650577_jprb, &
   &    0.2888331_jprb,    1.7176263_jprb, &
   &    0.1727553_jprb,    1.8412462_jprb, &
   &    0.0917505_jprb,    1.9312177_jprb, &
   &    0.0425524_jprb,    1.9847447_jprb, &
   &    0.0170291_jprb,    2.0011383_jprb, &
   &    0.0058401_jprb,    1.9816189_jprb, &
   &    0.0017113_jprb,    1.9286890_jprb, &
   &    0.0004282_jprb,    1.8462520_jprb, &
   &    0.0000916_jprb,    1.7390198_jprb, &
   &    0.0000167_jprb,    1.6125971_jprb, &
   &    0.0000026_jprb,    1.4722552_jprb, &
   &    0.0000004_jprb,    1.3235979_jprb, &
   &    0.0000000_jprb,    1.1718219_jprb, &
   &    0.0000000_jprb,    1.0216426_jprb, &
   &    0.0000000_jprb,    0.8771323_jprb, &
   &    0.0000000_jprb,    0.7414203_jprb, &
   &    0.0000000_jprb,    0.6170549_jprb  &
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
   &    0.6814394_jprb,    0.9867433_jprb,    1.3722117_jprb, &
   &    0.4248339_jprb,    0.9472399_jprb,    1.7867542_jprb, &
   &    0.2346598_jprb,    0.8829374_jprb,    2.2337751_jprb, &
   &    0.1104546_jprb,    0.7966385_jprb,    2.7019865_jprb, &
   &    0.0425238_jprb,    0.6931569_jprb,    3.1793568_jprb, &
   &    0.0129424_jprb,    0.5792789_jprb,    3.6548119_jprb, &
   &    0.0030503_jprb,    0.4630838_jprb,    4.1185323_jprb, &
   &    0.0005518_jprb,    0.3527570_jprb,    4.5613796_jprb, &
   &    0.0000764_jprb,    0.2551593_jprb,    4.9778419_jprb, &
   &    0.0000081_jprb,    0.1747017_jprb,    5.3607725_jprb, &
   &    0.0000007_jprb,    0.1129345_jprb,    5.7075865_jprb, &
   &    0.0000000_jprb,    0.0687933_jprb,    6.0134069_jprb, &
   &    0.0000000_jprb,    0.0394370_jprb,    6.2777727_jprb, &
   &    0.0000000_jprb,    0.0212605_jprb,    6.4985616_jprb, &
   &    0.0000000_jprb,    0.0107762_jprb,    6.6755274_jprb, &
   &    0.0000000_jprb,    0.0051357_jprb,    6.8085711_jprb, &
   &    0.0000000_jprb,    0.0023022_jprb,    6.8984824_jprb, &
   &    0.0000000_jprb,    0.0009710_jprb,    6.9470503_jprb, &
   &    0.0000000_jprb,    0.0003855_jprb,    6.9554890_jprb, &
   &    0.0000000_jprb,    0.0001441_jprb,    6.9253733_jprb  &
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
   &    0.5927012_jprb,    0.8465862_jprb,    1.1416107_jprb,    1.5265063_jprb, &
   &    0.3054363_jprb,    0.6881427_jprb,    1.2647692_jprb,    2.1629128_jprb, &
   &    0.1292354_jprb,    0.5319233_jprb,    1.3643437_jprb,    2.9001567_jprb, &
   &    0.0419627_jprb,    0.3863907_jprb,    1.4352144_jprb,    3.7260156_jprb, &
   &    0.0098352_jprb,    0.2601423_jprb,    1.4738697_jprb,    4.6281102_jprb, &
   &    0.0016016_jprb,    0.1600257_jprb,    1.4786155_jprb,    5.5904412_jprb, &
   &    0.0001784_jprb,    0.0887642_jprb,    1.4500311_jprb,    6.6002152_jprb, &
   &    0.0000136_jprb,    0.0439287_jprb,    1.3905235_jprb,    7.6482125_jprb, &
   &    0.0000007_jprb,    0.0192542_jprb,    1.3043413_jprb,    8.7184701_jprb, &
   &    0.0000000_jprb,    0.0074432_jprb,    1.1969975_jprb,    9.8060001_jprb, &
   &    0.0000000_jprb,    0.0025336_jprb,    1.0748114_jprb,   10.9000539_jprb, &
   &    0.0000000_jprb,    0.0007592_jprb,    0.9442652_jprb,   11.9917174_jprb, &
   &    0.0000000_jprb,    0.0002004_jprb,    0.8114326_jprb,   13.0777795_jprb, &
   &    0.0000000_jprb,    0.0000466_jprb,    0.6819954_jprb,   14.1475366_jprb, &
   &    0.0000000_jprb,    0.0000096_jprb,    0.5603933_jprb,   15.2012100_jprb, &
   &    0.0000000_jprb,    0.0000017_jprb,    0.4501054_jprb,   16.2296065_jprb, &
   &    0.0000000_jprb,    0.0000003_jprb,    0.3532535_jprb,   17.2344830_jprb, &
   &    0.0000000_jprb,    0.0000000_jprb,    0.2708529_jprb,   18.2043706_jprb, &
   &    0.0000000_jprb,    0.0000000_jprb,    0.2028404_jprb,   19.1448734_jprb, &
   &    0.0000000_jprb,    0.0000000_jprb,    0.1483481_jprb,   20.0477800_jprb  &
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
   &    0.5226172_jprb,    0.7415683_jprb,    0.9867433_jprb,    1.2812262_jprb,    1.6685114_jprb, &
   &    0.2236963_jprb,    0.5145558_jprb,    0.9472399_jprb,    1.5746715_jprb,    2.5257058_jprb, &
   &    0.0716935_jprb,    0.3273229_jprb,    0.8829374_jprb,    1.8712092_jprb,    3.5661068_jprb, &
   &    0.0156017_jprb,    0.1859842_jprb,    0.7966385_jprb,    2.1612939_jprb,    4.7777287_jprb, &
   &    0.0021408_jprb,    0.0916871_jprb,    0.6931569_jprb,    2.4352272_jprb,    6.1487340_jprb, &
   &    0.0001791_jprb,    0.0381751_jprb,    0.5792789_jprb,    2.6847897_jprb,    7.6647237_jprb, &
   &    0.0000091_jprb,    0.0131562_jprb,    0.4630838_jprb,    2.9038421_jprb,    9.3067221_jprb, &
   &    0.0000003_jprb,    0.0037096_jprb,    0.3527570_jprb,    3.0866813_jprb,   11.0674629_jprb, &
   &    0.0000000_jprb,    0.0008517_jprb,    0.2551593_jprb,    3.2307859_jprb,   12.9270164_jprb, &
   &    0.0000000_jprb,    0.0001591_jprb,    0.1747017_jprb,    3.3344466_jprb,   14.8798783_jprb, &
   &    0.0000000_jprb,    0.0000242_jprb,    0.1129345_jprb,    3.3968674_jprb,   16.9106920_jprb, &
   &    0.0000000_jprb,    0.0000030_jprb,    0.0687933_jprb,    3.4196745_jprb,   19.0108320_jprb, &
   &    0.0000000_jprb,    0.0000003_jprb,    0.0394370_jprb,    3.4041109_jprb,   21.1765462_jprb, &
   &    0.0000000_jprb,    0.0000000_jprb,    0.0212605_jprb,    3.3533316_jprb,   23.3883436_jprb, &
   &    0.0000000_jprb,    0.0000000_jprb,    0.0107762_jprb,    3.2702962_jprb,   25.6488585_jprb, &
   &    0.0000000_jprb,    0.0000000_jprb,    0.0051357_jprb,    3.1586362_jprb,   27.9520814_jprb, &
   &    0.0000000_jprb,    0.0000000_jprb,    0.0023022_jprb,    3.0228078_jprb,   30.2803352_jprb, &
   &    0.0000000_jprb,    0.0000000_jprb,    0.0009710_jprb,    2.8669688_jprb,   32.6427096_jprb, &
   &    0.0000000_jprb,    0.0000000_jprb,    0.0003855_jprb,    2.6950136_jprb,   35.0229974_jprb, &
   &    0.0000000_jprb,    0.0000000_jprb,    0.0001441_jprb,    2.5115262_jprb,   37.4183027_jprb  &
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
   &    0.4649384_jprb,    0.6575435_jprb,    0.8692894_jprb,    1.1144298_jprb,    1.4107885_jprb,    1.8035396_jprb, &
   &    0.1655820_jprb,    0.3909163_jprb,    0.7279844_jprb,    1.2068560_jprb,    1.8806884_jprb,    2.8794333_jprb, &
   &    0.0397520_jprb,    0.2024317_jprb,    0.5827610_jprb,    1.2732663_jprb,    2.3970312_jprb,    4.2336833_jprb, &
   &    0.0056642_jprb,    0.0872529_jprb,    0.4415205_jprb,    1.3097355_jprb,    2.9496924_jprb,    5.8522408_jprb, &
   &    0.0004434_jprb,    0.0298761_jprb,    0.3129573_jprb,    1.3141606_jprb,    3.5262524_jprb,    7.7240713_jprb, &
   &    0.0000186_jprb,    0.0078393_jprb,    0.2050423_jprb,    1.2863942_jprb,    4.1137350_jprb,    9.8355955_jprb, &
   &    0.0000004_jprb,    0.0015452_jprb,    0.1227669_jprb,    1.2284878_jprb,    4.7011426_jprb,   12.1717083_jprb, &
   &    0.0000000_jprb,    0.0002271_jprb,    0.0665317_jprb,    1.1445065_jprb,    5.2795616_jprb,   14.7142934_jprb, &
   &    0.0000000_jprb,    0.0000249_jprb,    0.0324005_jprb,    1.0400598_jprb,    5.8404992_jprb,   17.4528483_jprb, &
   &    0.0000000_jprb,    0.0000020_jprb,    0.0141167_jprb,    0.9216928_jprb,    6.3769440_jprb,   20.3726358_jprb, &
   &    0.0000000_jprb,    0.0000001_jprb,    0.0054905_jprb,    0.7962454_jprb,    6.8829759_jprb,   23.4626696_jprb, &
   &    0.0000000_jprb,    0.0000000_jprb,    0.0019050_jprb,    0.6703773_jprb,    7.3549780_jprb,   26.7209257_jprb, &
   &    0.0000000_jprb,    0.0000000_jprb,    0.0005898_jprb,    0.5497090_jprb,    7.7882857_jprb,   30.1173978_jprb, &
   &    0.0000000_jprb,    0.0000000_jprb,    0.0001630_jprb,    0.4388950_jprb,    8.1816260_jprb,   33.6627997_jprb, &
   &    0.0000000_jprb,    0.0000000_jprb,    0.0000403_jprb,    0.3409924_jprb,    8.5315371_jprb,   37.3332192_jprb, &
   &    0.0000000_jprb,    0.0000000_jprb,    0.0000089_jprb,    0.2576965_jprb,    8.8385926_jprb,   41.1332983_jprb, &
   &    0.0000000_jprb,    0.0000000_jprb,    0.0000018_jprb,    0.1893815_jprb,    9.0998980_jprb,   45.0446523_jprb, &
   &    0.0000000_jprb,    0.0000000_jprb,    0.0000003_jprb,    0.1352946_jprb,    9.3169600_jprb,   49.0554379_jprb, &
   &    0.0000000_jprb,    0.0000000_jprb,    0.0000000_jprb,    0.0939516_jprb,    9.4901347_jprb,   53.1901355_jprb, &
   &    0.0000000_jprb,    0.0000000_jprb,    0.0000000_jprb,    0.0634077_jprb,    9.6189996_jprb,   57.4001190_jprb  &
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
   &    0.416720_jprb,    0.588438_jprb,    0.774948_jprb,    0.986743_jprb,    1.234125_jprb,    1.535165_jprb,    1.932037_jprb, &
   &    0.123400_jprb,    0.299824_jprb,    0.567600_jprb,    0.947240_jprb,    1.468230_jprb,    2.183559_jprb,    3.228274_jprb, &
   &    0.021947_jprb,    0.124855_jprb,    0.386665_jprb,    0.882937_jprb,    1.693574_jprb,    2.938839_jprb,    4.898805_jprb, &
   &    0.002012_jprb,    0.039637_jprb,    0.239983_jprb,    0.796638_jprb,    1.902153_jprb,    3.786872_jprb,    6.940638_jprb, &
   &    0.000088_jprb,    0.009013_jprb,    0.132581_jprb,    0.693157_jprb,    2.086198_jprb,    4.714939_jprb,    9.335492_jprb, &
   &    0.000002_jprb,    0.001413_jprb,    0.063738_jprb,    0.579279_jprb,    2.239306_jprb,    5.708285_jprb,   12.079228_jprb, &
   &    0.000000_jprb,    0.000150_jprb,    0.026176_jprb,    0.463084_jprb,    2.357211_jprb,    6.752719_jprb,   15.143095_jprb, &
   &    0.000000_jprb,    0.000011_jprb,    0.009072_jprb,    0.352757_jprb,    2.436885_jprb,    7.839836_jprb,   18.529595_jprb, &
   &    0.000000_jprb,    0.000001_jprb,    0.002637_jprb,    0.255159_jprb,    2.477649_jprb,    8.953793_jprb,   22.212150_jprb, &
   &    0.000000_jprb,    0.000000_jprb,    0.000642_jprb,    0.174702_jprb,    2.480283_jprb,   10.088058_jprb,   26.174841_jprb, &
   &    0.000000_jprb,    0.000000_jprb,    0.000131_jprb,    0.112935_jprb,    2.446814_jprb,   11.230058_jprb,   30.414806_jprb, &
   &    0.000000_jprb,    0.000000_jprb,    0.000022_jprb,    0.068793_jprb,    2.380739_jprb,   12.377779_jprb,   34.913505_jprb, &
   &    0.000000_jprb,    0.000000_jprb,    0.000003_jprb,    0.039437_jprb,    2.285681_jprb,   13.518538_jprb,   39.658049_jprb, &
   &    0.000000_jprb,    0.000000_jprb,    0.000000_jprb,    0.021260_jprb,    2.166583_jprb,   14.651661_jprb,   44.637921_jprb, &
   &    0.000000_jprb,    0.000000_jprb,    0.000000_jprb,    0.010776_jprb,    2.027992_jprb,   15.766711_jprb,   49.868737_jprb, &
   &    0.000000_jprb,    0.000000_jprb,    0.000000_jprb,    0.005136_jprb,    1.875078_jprb,   16.862527_jprb,   55.312810_jprb, &
   &    0.000000_jprb,    0.000000_jprb,    0.000000_jprb,    0.002302_jprb,    1.712769_jprb,   17.933645_jprb,   60.973715_jprb, &
   &    0.000000_jprb,    0.000000_jprb,    0.000000_jprb,    0.000971_jprb,    1.545684_jprb,   18.975474_jprb,   66.836273_jprb, &
   &    0.000000_jprb,    0.000000_jprb,    0.000000_jprb,    0.000385_jprb,    1.378294_jprb,   19.986007_jprb,   72.891393_jprb, &
   &    0.000000_jprb,    0.000000_jprb,    0.000000_jprb,    0.000144_jprb,    1.214433_jprb,   20.964243_jprb,   79.138372_jprb  &
   & ], [7,nfsd])


contains

  !---------------------------------------------------------------------
  ! Compute the weights and nodes for integrating across a gamma
  ! distribution given the fractional standard deviation (=standard
  ! deviation divided by the mean) assuming a mean of unity, using
  ! Gauss-Hermite Quadrature. The user can then multiply the nodes
  ! output from this routine by the mean of the distribution.  If the
  ! quadrature order "norder" is less than 1, the weights and nodes
  ! will not be written, while if it is larger than 7, only 7 nodes
  ! will be written and the remainder will be assigned zero weight.
  subroutine calc_gauss_gamma(ng, norder, fsd, weights, nodes)

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

  end subroutine calc_gauss_gamma

end module radiation_gauss_gamma
