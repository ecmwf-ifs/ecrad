! radiation_gauss_logbeta.F90 - Coefficients for Gauss Quadrature with logbeta weighting
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

module radiation_gauss_logbeta

  use parkind1, only : jpim, jprb

  implicit none

  public

  ! The look-up table covers fractional standard deviations in the
  ! range 0 to 6

  ! Number of fractional standard deviations (FSD) in look-up table
  integer(jpim), parameter :: nfsd = 31

  ! Maximum number of quadrature points
  integer(jpim), parameter :: nmaxorder = 7

  ! Spacing of fractional standard deviations in look-up table
  real(jprb), parameter :: dfsd = 0.2_jprb

  ! Weights and nodes for a beta distribution of the form
  ! x^alpha*(1-x)^beta between -1 and 1, where alpha=30 and beta=12,
  ! for quadrature orders between 2 and 7
  real(jprb), parameter :: weights_lut(nmaxorder,2:nmaxorder) = reshape([ &
   & 0.43517283_jprb, 0.56482717_jprb, 0.00000000_jprb, 0.00000000_jprb, 0.00000000_jprb, 0.00000000_jprb, 0.00000000_jprb, &
   & 0.11940635_jprb, 0.63669247_jprb, 0.24390118_jprb, 0.00000000_jprb, 0.00000000_jprb, 0.00000000_jprb, 0.00000000_jprb, &
   & 0.02693743_jprb, 0.34888686_jprb, 0.52691917_jprb, 0.09725654_jprb, 0.00000000_jprb, 0.00000000_jprb, 0.00000000_jprb, &
   & 0.00551856_jprb, 0.13598842_jprb, 0.47733558_jprb, 0.34314991_jprb, 0.03800754_jprb, 0.00000000_jprb, 0.00000000_jprb, &
   & 0.00107838_jprb, 0.04376495_jprb, 0.28116568_jprb, 0.46112581_jprb, 0.19794849_jprb, 0.01491668_jprb, 0.00000000_jprb, &
   & 0.00020666_jprb, 0.01254801_jprb, 0.12844383_jprb, 0.37888647_jprb, 0.36704848_jprb, 0.10692111_jprb, 0.00594546_jprb  &
   & ], [nmaxorder,nmaxorder-1])
  real(jprb), parameter :: nodes_lut(nmaxorder,2:nmaxorder) = reshape([ &
   & 0.62705994_jprb, 0.76424440_jprb, 0.00000000_jprb, 0.00000000_jprb, 0.00000000_jprb, 0.00000000_jprb, 0.00000000_jprb, &
   & 0.56776739_jprb, 0.69294656_jprb, 0.80178604_jprb, 0.00000000_jprb, 0.00000000_jprb, 0.00000000_jprb, 0.00000000_jprb, &
   & 0.51880779_jprb, 0.63574091_jprb, 0.73666653_jprb, 0.82878477_jprb, 0.00000000_jprb, 0.00000000_jprb, 0.00000000_jprb, &
   & 0.47703412_jprb, 0.58702052_jprb, 0.68263959_jprb, 0.76919611_jprb, 0.84949428_jprb, 0.00000000_jprb, 0.00000000_jprb, &
   & 0.44072304_jprb, 0.54448008_jprb, 0.63552501_jprb, 0.71848247_jprb, 0.79477408_jprb, 0.86601532_jprb, 0.00000000_jprb, &
   & 0.40877045_jprb, 0.50681503_jprb, 0.59362297_jprb, 0.67341172_jprb, 0.74724683_jprb, 0.81557800_jprb, 0.87955502_jprb  &
   & ], [nmaxorder,nmaxorder-1])

  ! The scale_factor and add_offset are functions of FSD from 0 to 6,
  ! and the nodes for a particular FSD are found from
  ! node*scale_factor(FSD)+add_offset(FSD).
  real(jprb), parameter :: scale_factor(nfsd) = [ &
   &    0.0000000_jprb, &
   &    2.9615086_jprb, &
   &    5.9444930_jprb, &
   &    8.8261148_jprb, &
   &   11.5294927_jprb, &
   &   14.0236521_jprb, &
   &   16.3073088_jprb, &
   &   18.3954084_jprb, &
   &   20.3074288_jprb, &
   &   22.0648828_jprb, &
   &   23.6868571_jprb, &
   &   25.1907103_jprb, &
   &   26.5907377_jprb, &
   &   27.8994722_jprb, &
   &   29.1280000_jprb, &
   &   30.2851173_jprb, &
   &   31.3787409_jprb, &
   &   32.4151820_jprb, &
   &   33.4002873_jprb, &
   &   34.3386342_jprb, &
   &   35.2350406_jprb, &
   &   36.0927689_jprb, &
   &   36.9153560_jprb, &
   &   37.7055509_jprb, &
   &   38.4659337_jprb, &
   &   39.1989111_jprb, &
   &   39.9058551_jprb, &
   &   40.5894558_jprb, &
   &   41.2505629_jprb, &
   &   41.8910471_jprb, &
   &   42.5122083_jprb  &
   & ]
  real(jprb), parameter :: add_offset(nfsd) = [ &
   &    0.0000000_jprb, &
   &   -2.1357891_jprb, &
   &   -4.3266978_jprb, &
   &   -6.4786490_jprb, &
   &   -8.5272931_jprb, &
   &  -10.4414400_jprb, &
   &  -12.2131841_jprb, &
   &  -13.8483628_jprb, &
   &  -15.3577937_jprb, &
   &  -16.7549531_jprb, &
   &  -18.0523755_jprb, &
   &  -19.2618770_jprb, &
   &  -20.3933775_jprb, &
   &  -21.4557417_jprb, &
   &  -22.4569424_jprb, &
   &  -23.4033531_jprb, &
   &  -24.3007813_jprb, &
   &  -25.1538782_jprb, &
   &  -25.9669979_jprb, &
   &  -26.7435555_jprb, &
   &  -27.4871947_jprb, &
   &  -28.2003714_jprb, &
   &  -28.8857818_jprb, &
   &  -29.5455215_jprb, &
   &  -30.1815691_jprb, &
   &  -30.7957805_jprb, &
   &  -31.3891969_jprb, &
   &  -31.9639235_jprb, &
   &  -32.5206012_jprb, &
   &  -33.0606960_jprb, &
   &  -33.5852225_jprb  &
   & ]


contains
  !---------------------------------------------------------------------
  ! Compute the weights and nodes for integrating across a logbeta
  ! distribution given the fractional standard deviation (=standard
  ! deviation divided by the mean) assuming a mean of unity, using
  ! Gauss-Jacobi Quadrature. The user can then multiply the nodes
  ! output from this routine by the mean of the distribution.  If the
  ! quadrature order "norder" is less than 1, the weights and nodes
  ! will not be written, while if it is larger than 7, only 7 nodes
  ! will be written and the remainder will be assigned zero weight.
  subroutine calc_gauss_logbeta(ng, norder, fsd, weights, nodes)

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

    ! Capped quadrature order
    integer(jpim) :: myorder

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

    myorder = min(nmaxorder, norder)

    if (myorder > 1) then
      do jg = 1,ng
        weights(jg,1:myorder) = weights_lut(1:myorder,myorder)
        nodes(jg,1:myorder)   = nodes_lut(1:myorder,myorder)
        nodes(jg,1:myorder)   = exp(nodes(jg,1:myorder)*(w1(jg)*scale_factor(ifsd(jg)) &
        &                                               +w2(jg)*scale_factor(ifsd(jg)+1)) &
        &                     + w1(jg)*add_offset(ifsd(jg)) + w2(jg)*add_offset(ifsd(jg)+1))
        ! Optionally provide a local normalization to ensure the mean
        ! optical depth is exactly the requested value
        !nodes(jg,1:myorder) = nodes(jg,1:myorder) &
        !     &                / sum(weights(jg,1:myorder)*nodes(jg,1:myorder))
      end do
      weights(:,myorder+1:) = 0.0_jprb
      nodes(:,myorder+1:)  = 1.0_jprb
    else
      weights(:,1)  = 1.0_jprb
      weights(:,2:) = 0.0_jprb
      nodes         = 1.0_jprb
    end if

  end subroutine calc_gauss_logbeta

end module radiation_gauss_logbeta
