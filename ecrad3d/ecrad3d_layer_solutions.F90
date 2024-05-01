module ecrad3d_layer_solutions

  public

contains
  
  !---------------------------------------------------------------------
  ! Compute the shortwave reflectance and transmittance to diffuse
  ! radiation using the Meador & Weaver formulas, as well as the
  ! "direct" reflection and transmission, which really means the rate
  ! of transfer of direct solar radiation (into a plane perpendicular
  ! to the direct beam) into diffuse upward and downward streams at
  ! the top and bottom of the layer, respectively.  Finally,
  ! trans_dir_dir is the transmittance of the atmosphere to direct
  ! radiation with no scattering. This version incorporates the
  ! calculation of the gamma terms.
  subroutine calc_ref_trans_sw(ncol, nlay, nmaxspec, nspec, &
       &      mu0, od, ssa, &
       &      asymmetry, ref_diff, trans_diff, &
       &      ref_dir, trans_dir_diff, trans_dir_dir)

    use parkind1, only : jprb
    
#ifdef DO_DR_HOOK_TWO_STREAM
    use yomhook, only : lhook, dr_hook, jphook
#endif

    implicit none
    
    integer, intent(in) :: ncol, nlay, nmaxspec, nspec

    ! Cosine of solar zenith angle
    real(jprb), intent(in) :: mu0(ncol)

    ! Optical depth and single scattering albedo
    real(jprb), intent(in), dimension(ncol,nmaxspec,nlay) :: od, ssa, asymmetry

    ! The direct reflectance and transmittance, i.e. the fraction of
    ! incoming direct solar radiation incident at the top of a layer
    ! that is either reflected back (ref_dir) or scattered but
    ! transmitted through the layer to the base (trans_dir_diff)
    real(jprb), intent(out), dimension(ncol,nspec,nlay) :: ref_dir, trans_dir_diff

    ! The diffuse reflectance and transmittance, i.e. the fraction of
    ! diffuse radiation incident on a layer from either top or bottom
    ! that is reflected back or transmitted through
    real(jprb), intent(out), dimension(ncol,nspec,nlay) :: ref_diff, trans_diff

    ! Transmittance of the direct been with no scattering
    real(jprb), intent(out), dimension(ncol,nspec,nlay) :: trans_dir_dir

    ! The three transfer coefficients from the two-stream
    ! differentiatial equations 
    real(jprb), dimension(:), allocatable :: gamma1, gamma2, gamma3, gamma4 
    real(jprb), dimension(:), allocatable :: alpha1, alpha2, k_exponent
    real(jprb), dimension(:), allocatable :: exponential ! = exp(-k_exponent*od)
    
    real(jprb) :: reftrans_factor, factor
    real(jprb) :: exponential2 ! = exp(-2*k_exponent*od)
    real(jprb) :: k_mu0, k_gamma3, k_gamma4
    real(jprb) :: k_2_exponential, one_minus_kmu0_sqr
    integer    :: jc, js, jl

#ifdef DO_DR_HOOK_TWO_STREAM
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('ecrad3d_layer_solutions:calc_ref_trans_sw',0,hook_handle)
#endif

    !$OMP PARALLEL PRIVATE(gamma1, gamma2, gamma3, gamma4, alpha1, alpha2, &
    !$OMP&                 k_exponent, exponential, reftrans_factor, factor, &
    !$OMP&                 exponential2, k_mu0, k_gamma3, k_gamma4, &
    !$OMP&                 k_2_exponential, one_minus_kmu0_sqr, jc, js)
    allocate(gamma1(ncol))
    allocate(gamma2(ncol))
    allocate(gamma3(ncol))
    allocate(gamma4(ncol))
    allocate(alpha1(ncol))
    allocate(alpha2(ncol))
    allocate(k_exponent(ncol))
    allocate(exponential(ncol))
    
    !$OMP DO SCHEDULE(DYNAMIC)
    do jl = 1,nlay
      do js = 1,nspec
        ! GCC 9.3 strange error: intermediate values of ~ -8000 cause
        ! a FPE when vectorizing exp(), but not in non-vectorized
        ! loop, nor with larger negative values!
        trans_dir_dir(:,js,jl) = max(-max(od(:,js,jl) * (1.0_jprb/mu0(:)), 0.0_jprb),-1000.0_jprb)
        trans_dir_dir(:,js,jl) = exp(trans_dir_dir(:,js,jl))
        do jc = 1,ncol
          ! Zdunkowski "PIFM" (Zdunkowski et al., 1980; Contributions
          ! to Atmospheric Physics 53, 147-66)
          factor = 0.75_jprb*asymmetry(jc,js,jl)
          gamma1(jc) = 2.0_jprb  - ssa(jc,js,jl) * (1.25_jprb + factor)
          gamma2(jc) = ssa(jc,js,jl) * (0.75_jprb - factor)
          gamma3(jc) = 0.5_jprb  - mu0(jc)*factor
          gamma4(jc) = 1.0_jprb - gamma3(jc)
          alpha1(jc) = gamma1(jc)*gamma4(jc) + gamma2(jc)*gamma3(jc) ! Eq. 16
          alpha2(jc) = gamma1(jc)*gamma3(jc) + gamma2(jc)*gamma4(jc) ! Eq. 17
          ! The following line crashes inexplicably with gfortran 8.5.0 in
          ! single precision - try a later version. Note that the minimum
          ! value is needed to produce correct results for single
          ! scattering albedos very close to or equal to one.
#ifdef PARKIND1_SINGLE
          k_exponent(jc) = sqrt(max((gamma1(jc) - gamma2(jc)) * (gamma1(jc) + gamma2(jc)), &
               &       1.0e-6_jprb)) ! Eq 18
#else
          k_exponent(jc) = sqrt(max((gamma1(jc) - gamma2(jc)) * (gamma1(jc) + gamma2(jc)), &
               &       1.0e-12_jprb)) ! Eq 18
#endif
        end do

        exponential = exp(-k_exponent*od(:,js,jl))

        do jc = 1,ncol
          k_mu0 = k_exponent(jc)*mu0(jc)
          one_minus_kmu0_sqr = 1.0_jprb - k_mu0*k_mu0
          k_gamma3 = k_exponent(jc)*gamma3(jc)
          k_gamma4 = k_exponent(jc)*gamma4(jc)
          exponential2 = exponential(jc)*exponential(jc)
          k_2_exponential = 2.0_jprb * k_exponent(jc) * exponential(jc)
          reftrans_factor = 1.0_jprb / (k_exponent(jc) + gamma1(jc) &
               &                 + (k_exponent(jc) - gamma1(jc))*exponential2)
        
          ! Meador & Weaver (1980) Eq. 25
          ref_diff(jc,js,jl) = gamma2(jc) * (1.0_jprb - exponential2) * reftrans_factor

          ! Meador & Weaver (1980) Eq. 26, with security (which is
          ! sometimes needed, but apparently not on ref_diff)
          trans_diff(jc,js,jl) = max(0.0_jprb, min(k_2_exponential * reftrans_factor, &
               &                             1.0_jprb-ref_diff(jc,js,jl)))

          ! Here we need mu0 even though it wasn't in Meador and Weaver
          ! because we are assuming the incoming direct flux is defined to
          ! be the flux into a plane perpendicular to the direction of the
          ! sun, not into a horizontal plane
          reftrans_factor = mu0(jc) * ssa(jc,js,jl) * reftrans_factor &
               &  / merge(one_minus_kmu0_sqr, epsilon(1.0_jprb), abs(one_minus_kmu0_sqr) > epsilon(1.0_jprb))
      
          ! Meador & Weaver (1980) Eq. 14, multiplying top & bottom by
          ! exp(-k_exponent*od) in case of very high optical depths
          ref_dir(jc,js,jl) = reftrans_factor &
               &  * ( (1.0_jprb - k_mu0) * (alpha2(jc) + k_gamma3) &
               &     -(1.0_jprb + k_mu0) * (alpha2(jc) - k_gamma3)*exponential2 &
               &     -k_2_exponential*(gamma3(jc) - alpha2(jc)*mu0(jc))*trans_dir_dir(jc,js,jl) )
        
          ! Meador & Weaver (1980) Eq. 15, multiplying top & bottom by
          ! exp(-k_exponent*od), minus the 1*exp(-od/mu0) term
          ! representing direct unscattered transmittance.
          trans_dir_diff(jc,js,jl) = reftrans_factor &
               & * ( k_2_exponential*(gamma4(jc) + alpha1(jc)*mu0(jc)) &
               &   - trans_dir_dir(jc,js,jl) &
               &   * ( (1.0_jprb + k_mu0) * (alpha1(jc) + k_gamma4) &
               &      -(1.0_jprb - k_mu0) * (alpha1(jc) - k_gamma4) * exponential2) )

          ! Final check that ref_dir + trans_dir_diff <= 1
          ref_dir(jc,js,jl)        = max(0.0_jprb, min(ref_dir(jc,js,jl), &
               &                   mu0(jc)*(1.0_jprb-trans_dir_dir(jc,js,jl))))
          trans_dir_diff(jc,js,jl) = max(0.0_jprb, min(trans_dir_diff(jc,js,jl), &
               &                   mu0(jc)*(1.0_jprb-trans_dir_dir(jc,js,jl))-ref_dir(jc,js,jl)))
        end do
      end do
    end do

    deallocate(gamma1, gamma2, gamma3, gamma4, alpha1, alpha2, &
         &     k_exponent, exponential)
    
    !$OMP END PARALLEL
    
#ifdef DO_DR_HOOK_TWO_STREAM
    if (lhook) call dr_hook('ecrad3d_two_stream:calc_ref_trans_sw',1,hook_handle)
#endif
 
  end subroutine calc_ref_trans_sw

end module ecrad3d_layer_solutions
