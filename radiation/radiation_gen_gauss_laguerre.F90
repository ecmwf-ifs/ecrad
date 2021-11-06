! radiation_gen_gauss_laguerre.F90 - Coefficients for Generalized Gauss-Laguerre Quadrature
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

module radiation_gen_gauss_laguerre

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
   &    0.5497519_jprb,    0.4502481_jprb, &
   &    0.5980581_jprb,    0.4019419_jprb, &
   &    0.6856953_jprb,    0.3143047_jprb, &
   &    0.7572479_jprb,    0.2427521_jprb, &
   &    0.8123475_jprb,    0.1876525_jprb, &
   &    0.8535534_jprb,    0.1464466_jprb, &
   &    0.8841106_jprb,    0.1158894_jprb, &
   &    0.9068667_jprb,   0.09313326_jprb, &
   &    0.9239992_jprb,   0.07600085_jprb, &
   &    0.9370786_jprb,   0.06292136_jprb, &
   &    0.9472136_jprb,    0.0527864_jprb, &
   &    0.9551832_jprb,   0.04481676_jprb, &
   &    0.9615385_jprb,   0.03846154_jprb, &
   &    0.9666728_jprb,    0.0333272_jprb, &
   &     0.970871_jprb,   0.02912904_jprb, &
   &    0.9743416_jprb,   0.02565835_jprb, &
   &      0.97724_jprb,   0.02276001_jprb, &
   &    0.9796828_jprb,   0.02031725_jprb, &
   &     0.981759_jprb,   0.01824105_jprb, &
   &    0.9835373_jprb,   0.01646273_jprb, &
   &    0.9850713_jprb,   0.01492875_jprb  &
   & ], [2,nfsd])
  real(jprb), parameter :: nodes2(2,nfsd) = reshape([ &
   &    1.0000000_jprb,    1.0000000_jprb, &
   &    0.8360392_jprb,    1.2439608_jprb, &
   &    0.7291868_jprb,    1.5908132_jprb, &
   &    0.6602858_jprb,    2.0597142_jprb, &
   &    0.6155001_jprb,    2.6644999_jprb, &
   &    0.5857864_jprb,    3.4142136_jprb, &
   &    0.5655401_jprb,    4.3144599_jprb, &
   &    0.5513489_jprb,    5.3686511_jprb, &
   &    0.5411260_jprb,    6.5788740_jprb, &
   &    0.5335731_jprb,    7.9464269_jprb, &
   &    0.5278640_jprb,    9.4721360_jprb, &
   &    0.5234598_jprb,   11.1565402_jprb, &
   &    0.5200000_jprb,   13.0000000_jprb, &
   &    0.5172381_jprb,   15.0027619_jprb, &
   &    0.5150015_jprb,   17.1649985_jprb, &
   &    0.5131670_jprb,   19.4868330_jprb, &
   &    0.5116450_jprb,   21.9683550_jprb, &
   &    0.5103693_jprb,   24.6096307_jprb, &
   &    0.5092900_jprb,   27.4107100_jprb, &
   &    0.5083691_jprb,   30.3716309_jprb, &
   &    0.5075775_jprb,   33.4924225_jprb  &
   & ], [2,nfsd])
  real(jprb), parameter :: weights3(3,nfsd) = reshape([ &
   &    0.2193311_jprb,    0.6571823_jprb,    0.1234866_jprb, &
   &    0.2795311_jprb,    0.6303667_jprb,   0.09010218_jprb, &
   &    0.4095007_jprb,    0.5426731_jprb,   0.04782615_jprb, &
   &    0.5316424_jprb,    0.4417012_jprb,   0.02665638_jprb, &
   &    0.6327083_jprb,    0.3512629_jprb,   0.01602879_jprb, &
   &     0.711093_jprb,    0.2785177_jprb,   0.01038926_jprb, &
   &    0.7703438_jprb,    0.2224788_jprb,  0.007177441_jprb, &
   &    0.8149387_jprb,    0.1798414_jprb,  0.005219891_jprb, &
   &     0.848723_jprb,    0.1473227_jprb,  0.003954283_jprb, &
   &    0.8746152_jprb,    0.1222903_jprb,  0.003094526_jprb, &
   &    0.8947291_jprb,    0.1027849_jprb,  0.002485996_jprb, &
   &    0.9105729_jprb,   0.08738674_jprb,  0.002040344_jprb, &
   &    0.9232225_jprb,   0.07507291_jprb,  0.001704543_jprb, &
   &    0.9334511_jprb,   0.06510351_jprb,  0.001445363_jprb, &
   &    0.9418202_jprb,   0.05693864_jprb,  0.001241192_jprb, &
   &    0.9487425_jprb,   0.05017996_jprb,  0.001077508_jprb, &
   &    0.9545256_jprb,   0.04453011_jprb, 0.0009442703_jprb, &
   &    0.9594012_jprb,   0.03976441_jprb, 0.0008343665_jprb, &
   &    0.9635463_jprb,   0.03571109_jprb, 0.0007426437_jprb, &
   &    0.9670973_jprb,   0.03223737_jprb, 0.0006652965_jprb, &
   &     0.970161_jprb,    0.0292395_jprb,  0.000599466_jprb  &
   & ], [3,nfsd])
  real(jprb), parameter :: nodes3(3,nfsd) = reshape([ &
   &    1.0000000_jprb,    1.0000000_jprb,    1.0000000_jprb, &
   &    0.7341575_jprb,    1.0531846_jprb,    1.4526579_jprb, &
   &    0.5839422_jprb,    1.2113066_jprb,    2.1647512_jprb, &
   &    0.4982403_jprb,    1.4717656_jprb,    3.1899942_jprb, &
   &    0.4474056_jprb,    1.8329197_jprb,    4.5596747_jprb, &
   &    0.4157746_jprb,    2.2942804_jprb,    6.2899451_jprb, &
   &    0.3951560_jprb,    2.8559266_jprb,    8.3889174_jprb, &
   &    0.3811439_jprb,    3.5180946_jprb,   10.8607615_jprb, &
   &    0.3712698_jprb,    4.2810143_jprb,   13.7077159_jprb, &
   &    0.3640907_jprb,    5.1448675_jprb,   16.9310418_jprb, &
   &    0.3587286_jprb,    6.1097882_jprb,   20.5314832_jprb, &
   &    0.3546296_jprb,    7.1758724_jprb,   24.5094980_jprb, &
   &    0.3514323_jprb,    8.3431889_jprb,   28.8653788_jprb, &
   &    0.3488942_jprb,    9.6117870_jprb,   33.5993188_jprb, &
   &    0.3468480_jprb,   10.9817028_jprb,   38.7114492_jprb, &
   &    0.3451758_jprb,   12.4529625_jprb,   44.2018617_jprb, &
   &    0.3437927_jprb,   14.0255857_jprb,   50.0706217_jprb, &
   &    0.3426362_jprb,   15.6995872_jprb,   56.3177766_jprb, &
   &    0.3416598_jprb,   17.4749783_jprb,   62.9433618_jprb, &
   &    0.3408283_jprb,   19.3517677_jprb,   69.9474040_jprb, &
   &    0.3401145_jprb,   21.3299622_jprb,   77.3299233_jprb  &
   & ], [3,nfsd])
  real(jprb), parameter :: weights4(4,nfsd) = reshape([ &
   &   0.07703671_jprb,    0.5247076_jprb,    0.3721074_jprb,   0.02614831_jprb, &
   &    0.1212432_jprb,     0.572072_jprb,    0.2919959_jprb,   0.01468886_jprb, &
   &    0.2418224_jprb,    0.5848138_jprb,    0.1684241_jprb,  0.004939671_jprb, &
   &    0.3779293_jprb,    0.5227364_jprb,   0.09735132_jprb,  0.001982909_jprb, &
   &    0.5019118_jprb,    0.4376179_jprb,   0.05951096_jprb,  0.000959309_jprb, &
   &    0.6031541_jprb,    0.3574187_jprb,   0.03888791_jprb, 0.0005392947_jprb, &
   &    0.6819533_jprb,    0.2907232_jprb,   0.02698461_jprb, 0.0003388365_jprb, &
   &    0.7423189_jprb,    0.2377742_jprb,   0.01967602_jprb,   0.00023088_jprb, &
   &    0.7885715_jprb,    0.1963316_jprb,   0.01492989_jprb, 0.0001670202_jprb, &
   &    0.8242907_jprb,    0.1638863_jprb,   0.01169658_jprb, 0.0001263908_jprb, &
   &    0.8521875_jprb,    0.1383098_jprb,  0.009403678_jprb, 9.902142e-05_jprb, &
   &    0.8742477_jprb,    0.1179504_jprb,  0.007722218_jprb,  7.97275e-05_jprb, &
   &    0.8919119_jprb,    0.1015685_jprb,  0.006453969_jprb, 6.561805e-05_jprb, &
   &    0.9062273_jprb,   0.08824336_jprb,  0.005474372_jprb, 5.498501e-05_jprb, &
   &    0.9179608_jprb,   0.07729024_jprb,  0.004702241_jprb, 4.676906e-05_jprb, &
   &    0.9276796_jprb,   0.06819722_jprb,  0.004082938_jprb, 4.028583e-05_jprb, &
   &    0.9358081_jprb,    0.0605782_jprb,  0.003578649_jprb, 3.507743e-05_jprb, &
   &    0.9426675_jprb,   0.05413915_jprb,   0.00316255_jprb, 3.082818e-05_jprb, &
   &    0.9485036_jprb,    0.0486539_jprb,    0.0028152_jprb, 2.731471e-05_jprb, &
   &    0.9535067_jprb,   0.04394675_jprb,  0.002522229_jprb,  2.43753e-05_jprb, &
   &    0.9578255_jprb,   0.03987981_jprb,  0.002272837_jprb, 2.189054e-05_jprb  &
   & ], [4,nfsd])
  real(jprb), parameter :: nodes4(4,nfsd) = reshape([ &
   &    1.0000000_jprb,    1.0000000_jprb,    1.0000000_jprb,    1.0000000_jprb, &
   &    0.6605266_jprb,    0.9314195_jprb,    1.2428620_jprb,    1.6451919_jprb, &
   &    0.4899950_jprb,    0.9950343_jprb,    1.6985784_jprb,    2.7363924_jprb, &
   &    0.4013801_jprb,    1.1608195_jprb,    2.3918190_jprb,    4.3659814_jprb, &
   &    0.3520317_jprb,    1.4133146_jprb,    3.3356554_jprb,    6.5789983_jprb, &
   &    0.3225477_jprb,    1.7457611_jprb,    4.5366203_jprb,    9.3950709_jprb, &
   &    0.3038289_jprb,    2.1553306_jprb,    5.9979856_jprb,   12.8228549_jprb, &
   &    0.2913302_jprb,    2.6408189_jprb,    7.7214427_jprb,   16.8664082_jprb, &
   &    0.2826292_jprb,    3.2016910_jprb,    9.7079054_jprb,   21.5277743_jprb, &
   &    0.2763578_jprb,    3.8376970_jprb,   11.9578917_jprb,   26.8080535_jprb, &
   &    0.2717034_jprb,    4.5487141_jprb,   14.4717083_jprb,   32.7078742_jprb, &
   &    0.2681625_jprb,    5.3346788_jprb,   17.2495447_jprb,   39.2276140_jprb, &
   &    0.2654107_jprb,    6.1955573_jprb,   20.2915224_jprb,   46.3675096_jprb, &
   &    0.2632326_jprb,    7.1313304_jprb,   23.5977218_jprb,   54.1277153_jprb, &
   &    0.2614807_jprb,    8.1419871_jprb,   27.1681974_jprb,   62.5083348_jprb, &
   &    0.2600517_jprb,    9.2275209_jprb,   31.0029876_jprb,   71.5094399_jprb, &
   &    0.2588715_jprb,   10.3879277_jprb,   35.1021196_jprb,   81.1310813_jprb, &
   &    0.2578860_jprb,   11.6232050_jprb,   39.4656133_jprb,   91.3732957_jprb, &
   &    0.2570549_jprb,   12.9333513_jprb,   44.0934835_jprb,  102.2361104_jprb, &
   &    0.2563477_jprb,   14.3183654_jprb,   48.9857413_jprb,  113.7195456_jprb, &
   &    0.2557411_jprb,   15.7782468_jprb,   54.1423952_jprb,  125.8236169_jprb  &
   & ], [4,nfsd])
  real(jprb), parameter :: weights5(5,nfsd) = reshape([ &
   &   0.02543686_jprb,    0.3181873_jprb,    0.5094245_jprb,    0.1422811_jprb,  0.004670252_jprb, &
   &   0.05149768_jprb,    0.4122326_jprb,    0.4478745_jprb,   0.08646434_jprb,  0.001930875_jprb, &
   &    0.1455323_jprb,     0.527301_jprb,     0.295269_jprb,   0.03150115_jprb, 0.0003965826_jprb, &
   &    0.2758596_jprb,    0.5294526_jprb,    0.1815356_jprb,   0.01303585_jprb, 0.0001163199_jprb, &
   &     0.407794_jprb,    0.4714102_jprb,    0.1143626_jprb,   0.00638657_jprb, 4.662396e-05_jprb, &
   &    0.5217556_jprb,    0.3986668_jprb,   0.07594245_jprb,  0.003611759_jprb, 2.336997e-05_jprb, &
   &    0.6133317_jprb,    0.3311849_jprb,   0.05319331_jprb,  0.002276443_jprb, 1.365113e-05_jprb, &
   &    0.6848595_jprb,    0.2745629_jprb,   0.03901475_jprb,  0.001554024_jprb, 8.862296e-06_jprb, &
   &    0.7403545_jprb,    0.2287943_jprb,   0.02971952_jprb,  0.001125509_jprb, 6.200098e-06_jprb, &
   &    0.7835761_jprb,    0.1922204_jprb,   0.02334652_jprb, 0.0008523843_jprb, 4.580796e-06_jprb, &
   &    0.8175351_jprb,    0.1629865_jprb,   0.01880666_jprb, 0.0006681716_jprb, 3.525969e-06_jprb, &
   &    0.8445069_jprb,    0.1394857_jprb,   0.01546642_jprb, 0.0005381966_jprb, 2.801214e-06_jprb, &
   &    0.8661751_jprb,    0.1204388_jprb,   0.01294073_jprb, 0.0004430846_jprb, 2.281729e-06_jprb, &
   &      0.88378_jprb,    0.1048606_jprb,   0.01098613_jprb, 0.0003713712_jprb, 1.896433e-06_jprb, &
   &    0.8982385_jprb,   0.09200084_jprb,  0.009443137_jprb, 0.0003159378_jprb, 1.602533e-06_jprb, &
   &    0.9102335_jprb,   0.08128889_jprb,  0.008204036_jprb, 0.0002721815_jprb, 1.373048e-06_jprb, &
   &    0.9202788_jprb,    0.0722889_jprb,  0.007194048_jprb, 0.0002370203_jprb, 1.190291e-06_jprb, &
   &    0.9287648_jprb,   0.06466582_jprb,  0.006360003_jprb, 0.0002083283_jprb, 1.042276e-06_jprb, &
   &    0.9359913_jprb,   0.05815994_jprb,   0.00566328_jprb, 0.0001846002_jprb, 9.206446e-07_jprb, &
   &    0.9421909_jprb,   0.05256828_jprb,  0.005075292_jprb, 0.0001647463_jprb, 8.194216e-07_jprb, &
   &     0.947546_jprb,   0.04773082_jprb,  0.004574518_jprb, 0.0001479611_jprb, 7.342414e-07_jprb  &
   & ], [5,nfsd])
  real(jprb), parameter :: nodes5(5,nfsd) = reshape([ &
   &    1.0000000_jprb,    1.0000000_jprb,    1.0000000_jprb,    1.0000000_jprb,    1.0000000_jprb, &
   &    0.6033591_jprb,    0.8421437_jprb,    1.1062717_jprb,    1.4198698_jprb,    1.8283556_jprb, &
   &    0.4233187_jprb,    0.8497619_jprb,    1.4218739_jprb,    2.1956486_jprb,    3.3093969_jprb, &
   &    0.3365111_jprb,    0.9627847_jprb,    1.9426537_jprb,    3.3839070_jprb,    5.5741435_jprb, &
   &    0.2903659_jprb,    1.1542686_jprb,    2.6674181_jprb,    5.0103596_jprb,    8.6775877_jprb, &
   &    0.2635603_jprb,    1.4134031_jprb,    3.5964258_jprb,    7.0858100_jprb,   12.6408008_jprb, &
   &    0.2468376_jprb,    1.7357930_jprb,    4.7301416_jprb,    9.6149437_jprb,   17.4722841_jprb, &
   &    0.2357983_jprb,    2.1195430_jprb,    6.0689197_jprb,   12.5999399_jprb,   23.1757991_jprb, &
   &    0.2281725_jprb,    2.5637758_jprb,    7.6129909_jprb,   16.0418916_jprb,   29.7531691_jprb, &
   &    0.2227059_jprb,    3.0680563_jprb,    9.3625010_jprb,   19.9413860_jprb,   37.2053509_jprb, &
   &    0.2186650_jprb,    3.6321540_jprb,   11.3175423_jprb,   24.2987577_jprb,   45.5328810_jprb, &
   &    0.2155998_jprb,    4.2559402_jprb,   13.4781745_jprb,   29.1142077_jprb,   54.7360778_jprb, &
   &    0.2132232_jprb,    4.9393388_jprb,   15.8444371_jprb,   34.3878621_jprb,   64.8151388_jprb, &
   &    0.2113453_jprb,    5.6823030_jprb,   18.4163570_jprb,   40.1198028_jprb,   75.7701920_jprb, &
   &    0.2098371_jprb,    6.4848029_jprb,   21.1939525_jprb,   46.3100849_jprb,   87.6013226_jprb, &
   &    0.2086083_jprb,    7.3468188_jprb,   24.1772370_jprb,   52.9587466_jprb,  100.3085894_jprb, &
   &    0.2075944_jprb,    8.2683371_jprb,   27.3662198_jprb,   60.0658147_jprb,  113.8920341_jprb, &
   &    0.2067483_jprb,    9.2493483_jprb,   30.7609080_jprb,   67.6313090_jprb,  128.3516864_jprb, &
   &    0.2060353_jprb,   10.2898458_jprb,   34.3613067_jprb,   75.6552437_jprb,  143.6875684_jprb, &
   &    0.2054290_jprb,   11.3898246_jprb,   38.1674200_jprb,   84.1376296_jprb,  159.8996968_jprb, &
   &    0.2049091_jprb,   12.5492811_jprb,   42.1792508_jprb,   93.0784751_jprb,  176.9880839_jprb  &
   & ], [5,nfsd])
  real(jprb), parameter :: weights6(6,nfsd) = reshape([ &
   &  0.008136579_jprb,    0.1652795_jprb,    0.4778739_jprb,    0.3058383_jprb,   0.04213606_jprb, 0.0007357211_jprb, &
   &    0.0219008_jprb,    0.2638498_jprb,    0.4883225_jprb,    0.2067111_jprb,   0.01899979_jprb,  0.000216044_jprb, &
   &   0.08996503_jprb,     0.438117_jprb,    0.3834913_jprb,   0.08416662_jprb,  0.004233512_jprb, 2.658794e-05_jprb, &
   &    0.2068295_jprb,    0.5001251_jprb,      0.25534_jprb,   0.03642204_jprb,   0.00127754_jprb,  5.80398e-06_jprb, &
   &    0.3385675_jprb,    0.4754932_jprb,    0.1672159_jprb,   0.01820345_jprb, 0.0005179473_jprb, 1.973028e-06_jprb, &
   &    0.4589647_jprb,    0.4170008_jprb,    0.1133734_jprb,    0.0103992_jprb, 0.0002610172_jprb, 8.985479e-07_jprb, &
   &     0.558841_jprb,    0.3540278_jprb,   0.08038572_jprb,   0.00659209_jprb, 0.0001529091_jprb, 4.947338e-07_jprb, &
   &    0.6383703_jprb,    0.2976003_jprb,   0.05941366_jprb,  0.004515983_jprb, 9.943862e-05_jprb, 3.090502e-07_jprb, &
   &    0.7008443_jprb,    0.2503168_jprb,   0.04549079_jprb,  0.003278278_jprb, 6.964388e-05_jprb, 2.106086e-07_jprb, &
   &    0.7499126_jprb,    0.2116853_jprb,   0.03586381_jprb,  0.002486702_jprb, 5.149291e-05_jprb, 1.527318e-07_jprb, &
   &     0.788695_jprb,    0.1803489_jprb,    0.0289648_jprb,  0.001951519_jprb, 3.965637e-05_jprb, 1.159684e-07_jprb, &
   &    0.8196319_jprb,    0.1548968_jprb,   0.02386646_jprb,  0.001573237_jprb, 3.151723e-05_jprb, 9.118917e-08_jprb, &
   &    0.8445671_jprb,    0.1341125_jprb,   0.01999861_jprb,  0.001296046_jprb, 2.567984e-05_jprb, 7.369193e-08_jprb, &
   &    0.8648774_jprb,    0.1170167_jprb,   0.01699763_jprb,  0.001086827_jprb, 2.134832e-05_jprb, 6.086776e-08_jprb, &
   &    0.8815908_jprb,    0.1028423_jprb,   0.01462381_jprb, 0.0009249699_jprb, 1.804307e-05_jprb,   5.1179e-08_jprb, &
   &    0.8954787_jprb,   0.09099429_jprb,    0.0127144_jprb,  0.000797122_jprb, 1.546148e-05_jprb, 4.367294e-08_jprb, &
   &    0.9071241_jprb,    0.0810121_jprb,   0.01115598_jprb, 0.0006943309_jprb, 1.340508e-05_jprb, 3.773405e-08_jprb, &
   &    0.9169723_jprb,   0.07253786_jprb,  0.009867624_jprb, 0.0006104137_jprb, 1.173926e-05_jprb, 3.295026e-08_jprb, &
   &    0.9253662_jprb,   0.06529196_jprb,  0.008790403_jprb, 0.0005409886_jprb, 1.037016e-05_jprb, 2.903723e-08_jprb, &
   &    0.9325728_jprb,   0.05905448_jprb,  0.007880594_jprb, 0.0004828797_jprb, 9.230606e-06_jprb, 2.579353e-08_jprb, &
   &    0.9388016_jprb,   0.05365112_jprb,  0.007105222_jprb, 0.0004337391_jprb, 8.271551e-06_jprb, 2.307311e-08_jprb  &
   & ], [6,nfsd])
  real(jprb), parameter :: nodes6(6,nfsd) = reshape([ &
   &    1.0000000_jprb,    1.0000000_jprb,    1.0000000_jprb,    1.0000000_jprb,    1.0000000_jprb,    1.0000000_jprb, &
   &    0.5570279_jprb,    0.7721714_jprb,    1.0048965_jprb,    1.2709300_jprb,    1.5895551_jprb,    2.0054191_jprb, &
   &    0.3732017_jprb,    0.7438383_jprb,    1.2308476_jprb,    1.8650413_jprb,    2.7022203_jprb,    3.8848508_jprb, &
   &    0.2898885_jprb,    0.8241764_jprb,    1.6444915_jprb,    2.8073827_jprb,    4.4277269_jprb,    6.8063340_jprb, &
   &    0.2471600_jprb,    0.9770085_jprb,    2.2331265_jprb,    4.1091197_jprb,    6.8017265_jprb,   10.8318588_jprb, &
   &    0.2228466_jprb,    1.1889321_jprb,    2.9927363_jprb,    5.7751436_jprb,    9.8374674_jprb,   15.9828740_jprb, &
   &    0.2078670_jprb,    1.4547098_jprb,    3.9220572_jprb,    7.8076626_jprb,   13.5402137_jprb,   22.2674896_jprb, &
   &    0.1980570_jprb,    1.7721147_jprb,    5.0206657_jprb,   10.2077326_jprb,   17.9122761_jprb,   29.6891539_jprb, &
   &    0.1913166_jprb,    2.1401158_jprb,    6.2884090_jprb,   12.9758927_jprb,   22.9547693_jprb,   38.2494966_jprb, &
   &    0.1865027_jprb,    2.5581976_jprb,    7.7252281_jprb,   16.1124362_jprb,   28.6682767_jprb,   47.9493587_jprb, &
   &    0.1829539_jprb,    3.0260840_jprb,    9.3310989_jprb,   19.6175317_jprb,   35.0531250_jprb,   58.7892067_jprb, &
   &    0.1802674_jprb,    3.5436183_jprb,   11.1060114_jprb,   23.4912808_jprb,   42.1095074_jprb,   70.7693146_jprb, &
   &    0.1781876_jprb,    4.1107069_jprb,   13.0499615_jprb,   27.7337478_jprb,   49.8375441_jprb,   83.8898520_jprb, &
   &    0.1765463_jprb,    4.7272914_jprb,   15.1629474_jprb,   32.3449746_jprb,   58.2373126_jprb,   98.1509277_jprb, &
   &    0.1752293_jprb,    5.3933338_jprb,   17.4449687_jprb,   37.3249893_jprb,   67.3088646_jprb,  113.5526144_jprb, &
   &    0.1741571_jprb,    6.1088087_jprb,   19.8960251_jprb,   42.6738114_jprb,   77.0522359_jprb,  130.0949617_jprb, &
   &    0.1732730_jprb,    6.8736987_jprb,   22.5161169_jprb,   48.3914550_jprb,   87.4674515_jprb,  147.7780049_jprb, &
   &    0.1725356_jprb,    7.6879915_jprb,   25.3052443_jprb,   54.4779299_jprb,   98.5545296_jprb,  166.6017691_jprb, &
   &    0.1719145_jprb,    8.5516779_jprb,   28.2634074_jprb,   60.9332437_jprb,  110.3134835_jprb,  186.5662730_jprb, &
   &    0.1713864_jprb,    9.4647516_jprb,   31.3906065_jprb,   67.7574018_jprb,  122.7443233_jprb,  207.6715304_jprb, &
   &    0.1709339_jprb,   10.4272075_jprb,   34.6868417_jprb,   74.9504086_jprb,  135.8470565_jprb,  229.9175518_jprb  &
   & ], [6,nfsd])
  real(jprb), parameter :: weights7(7,nfsd) = reshape([ &
   &  0.00256206_jprb,    0.077959_jprb,    0.358848_jprb,    0.416072_jprb,    0.134087_jprb,   0.0103669_jprb, 0.000104967_jprb, &
   &  0.00942248_jprb,    0.158067_jprb,    0.439903_jprb,    0.322218_jprb,   0.0670167_jprb,  0.00335092_jprb,  2.1269e-05_jprb, &
   &   0.0571872_jprb,     0.34974_jprb,    0.424796_jprb,    0.151368_jprb,   0.0164588_jprb, 0.000448909_jprb, 1.55472e-06_jprb, &
   &    0.158874_jprb,    0.456563_jprb,    0.309893_jprb,   0.0694258_jprb,   0.0051434_jprb, 0.000100796_jprb, 2.57333e-07_jprb, &
   &    0.286301_jprb,    0.463914_jprb,    0.212006_jprb,   0.0356283_jprb,  0.00211598_jprb, 3.46518e-05_jprb, 7.55845e-08_jprb, &
   &    0.409319_jprb,    0.421831_jprb,    0.147126_jprb,   0.0206335_jprb,  0.00107401_jprb, 1.58655e-05_jprb, 3.17032e-08_jprb, &
   &     0.51455_jprb,    0.365881_jprb,    0.105746_jprb,   0.0131825_jprb, 0.000631683_jprb, 8.76083e-06_jprb, 1.66087e-08_jprb, &
   &    0.599908_jprb,    0.311771_jprb,    0.078829_jprb,   0.0090749_jprb, 0.000411785_jprb, 5.48227e-06_jprb, 1.00488e-08_jprb, &
   &    0.667763_jprb,    0.264634_jprb,   0.0607013_jprb,  0.00660896_jprb, 0.000288857_jprb, 3.74022e-06_jprb, 6.70153e-09_jprb, &
   &    0.721489_jprb,    0.225224_jprb,   0.0480461_jprb,  0.00502437_jprb, 0.000213804_jprb, 2.71447e-06_jprb, 4.78642e-09_jprb, &
   &    0.764194_jprb,    0.192774_jprb,   0.0389155_jprb,  0.00394939_jprb, 0.000164785_jprb, 2.06221e-06_jprb, 3.59413e-09_jprb, &
   &    0.798402_jprb,    0.166143_jprb,   0.0321347_jprb,  0.00318766_jprb, 0.000131039_jprb, 1.62223e-06_jprb, 2.80268e-09_jprb, &
   &     0.82606_jprb,    0.144232_jprb,   0.0269712_jprb,  0.00262843_jprb, 0.000106816_jprb, 1.31136e-06_jprb, 2.25042e-09_jprb, &
   &    0.848643_jprb,    0.126108_jprb,   0.0229535_jprb,   0.0022057_jprb, 8.88286e-05_jprb, 1.08341e-06_jprb, 1.84946e-09_jprb, &
   &    0.867261_jprb,    0.111016_jprb,   0.0197681_jprb,  0.00187828_jprb, 7.50959e-05_jprb, 9.11127e-07_jprb, 1.54882e-09_jprb, &
   &    0.882756_jprb,   0.0983585_jprb,   0.0172013_jprb,  0.00161941_jprb, 6.43652e-05_jprb, 7.77617e-07_jprb, 1.31737e-09_jprb, &
   &    0.895764_jprb,   0.0876651_jprb,   0.0151032_jprb,  0.00141112_jprb, 5.58144e-05_jprb, 6.71956e-07_jprb, 1.13517e-09_jprb, &
   &    0.906776_jprb,   0.0785668_jprb,   0.0133666_jprb,  0.00124096_jprb, 4.88857e-05_jprb, 5.86828e-07_jprb, 9.89052e-10_jprb, &
   &     0.91617_jprb,   0.0707731_jprb,    0.011913_jprb,  0.00110011_jprb, 4.31897e-05_jprb, 5.17183e-07_jprb, 8.69964e-10_jprb, &
   &    0.924241_jprb,   0.0640538_jprb,   0.0106843_jprb,  0.00098216_jprb, 3.84477e-05_jprb, 4.59443e-07_jprb, 7.71554e-10_jprb, &
   &    0.931221_jprb,   0.0582255_jprb,  0.00963641_jprb, 0.000882377_jprb,  3.4456e-05_jprb, 4.11012e-07_jprb,  6.8924e-10_jprb  &
   & ], [7,nfsd])
  real(jprb), parameter :: nodes7(7,nfsd) = reshape([ &
   &    1.000000_jprb,    1.000000_jprb,    1.000000_jprb,    1.000000_jprb,    1.000000_jprb,    1.000000_jprb,    1.000000_jprb, &
   &    0.518376_jprb,    0.715068_jprb,    0.924753_jprb,    1.159276_jprb,    1.429789_jprb,    1.754517_jprb,    2.178222_jprb, &
   &    0.334007_jprb,    0.662557_jprb,    1.088684_jprb,    1.631998_jprb,    2.322462_jprb,    3.217274_jprb,    4.463018_jprb, &
   &    0.254708_jprb,    0.721236_jprb,    1.429353_jprb,    2.413188_jprb,    3.733503_jprb,    5.510764_jprb,    8.057248_jprb, &
   &    0.215182_jprb,    0.847617_jprb,    1.924744_jprb,    3.502406_jprb,    5.684117_jprb,    8.678206_jprb,   13.027729_jprb, &
   &    0.193044_jprb,    1.026665_jprb,    2.567877_jprb,    4.900353_jprb,    8.182153_jprb,   12.734180_jprb,   19.395728_jprb, &
   &    0.179531_jprb,    1.252709_jprb,    3.356390_jprb,    6.607593_jprb,   11.230746_jprb,   17.684133_jprb,   27.168899_jprb, &
   &    0.170735_jprb,    1.523375_jprb,    4.289390_jprb,    8.624461_jprb,   14.831275_jprb,   23.530363_jprb,   36.350401_jprb, &
   &    0.164714_jprb,    1.837571_jprb,    5.366498_jprb,   10.951151_jprb,   18.984409_jprb,   30.273953_jprb,   46.941704_jprb, &
   &    0.160426_jprb,    2.194750_jprb,    6.587536_jprb,   13.587774_jprb,   23.690499_jprb,   37.915461_jprb,   58.943554_jprb, &
   &    0.157271_jprb,    2.594620_jprb,    7.952414_jprb,   16.534396_jprb,   28.949741_jprb,   46.455195_jprb,   72.356364_jprb, &
   &    0.154886_jprb,    3.037011_jprb,    9.461082_jprb,   19.791058_jprb,   34.762252_jprb,   55.893336_jprb,   87.180375_jprb, &
   &    0.153041_jprb,    3.521825_jprb,   11.113510_jprb,   23.357788_jprb,   41.128105_jprb,   66.229996_jprb,  103.415734_jprb, &
   &    0.151587_jprb,    4.048998_jprb,   12.909682_jprb,   27.234603_jprb,   48.047347_jprb,   77.465247_jprb,  121.062536_jprb, &
   &    0.150421_jprb,    4.618490_jprb,   14.849585_jprb,   31.421516_jprb,   55.520007_jprb,   89.599137_jprb,  140.120844_jprb, &
   &    0.149472_jprb,    5.230272_jprb,   16.933213_jprb,   35.918533_jprb,   63.546109_jprb,  102.631698_jprb,  160.590703_jprb, &
   &    0.148690_jprb,    5.884325_jprb,   19.160560_jprb,   40.725662_jprb,   72.125667_jprb,  116.562954_jprb,  182.472141_jprb, &
   &    0.148038_jprb,    6.580637_jprb,   21.531622_jprb,   45.842907_jprb,   81.258693_jprb,  131.392922_jprb,  205.765181_jprb, &
   &    0.147489_jprb,    7.319197_jprb,   24.046398_jprb,   51.270271_jprb,   90.945193_jprb,  147.121613_jprb,  230.469839_jprb, &
   &    0.147023_jprb,    8.099997_jprb,   26.704884_jprb,   57.007756_jprb,  101.185176_jprb,  163.749037_jprb,  256.586128_jprb, &
   &    0.146623_jprb,    8.923033_jprb,   29.507080_jprb,   63.055365_jprb,  111.978644_jprb,  181.275201_jprb,  284.114055_jprb  &
   & ], [7,nfsd])

contains

  !---------------------------------------------------------------------
  ! Compute the weights and nodes for integrating across a gamma
  ! distribution given the fractional standard deviation (=standard
  ! deviation divided by the mean) assuming a mean of unity, using
  ! Generalized Gauss-Laguerre Quadrature. The user can then multiply
  ! the nodes output from this routine by the mean of the
  ! distribution.  If the quadrature order "norder" is less than 1,
  ! the weights and nodes will not be written, while if it is larger
  ! than 7, only 7 nodes will be written and the remainder will be
  ! assigned zero weight.
  subroutine calc_gen_gauss_laguerre(ng, norder, fsd, weights, nodes)

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

  end subroutine calc_gen_gauss_laguerre

end module radiation_gen_gauss_laguerre
