module RadCheck

  use dvr_spline_mod
  use dvr_diag_mod
  use DVRData
  use OrbData, only : orb

  use omp_lib

  implicit none

contains

  subroutine radial_check()

    integer                    :: i, j, a, b, c, d, l, l_val, norb
    real(idp)                  :: full_r_max, curr_r, curr_r_prime, integral_val, dr
    real(idp)                  :: integral_val_1,  integral_val_2
 
    !norb = orb%nSpatialOrbs
    norb = 2

    allocate(two_e_rad_int_num(norb,norb,norb,norb, 2*para%l+1))

    dr = 0.0001d0  

    do l = 1, 2*para%l+1
      l_val=l-1
 
      open(13,file="numerical_integrals_l"//trim(int2str(l_val))//".dat",      &
      &    form="formatted", action="write", recl=100000)

      do a = 1, norb
      do b = 1, norb
        do c = 1, norb
        do d = 1, norb
          integral_val = zero
          do i = 1, 2001
            curr_r = (i-1) * dr
            if (abs(dvr_primitive_val_esteban(a, para, grid, curr_r)) < 1d-16) cycle
            if (abs(dvr_primitive_val_esteban(c, para, grid, curr_r)) < 1d-16) cycle
            do j = 1, 2001
              curr_r_prime = (j-1) * dr
              if (curr_r + curr_r_prime < 1d-16) cycle
              integral_val = integral_val +                                        &
              &              dvr_primitive_val_esteban(a, para, grid, curr_r) *            &
              &              dvr_primitive_val_esteban(c, para, grid, curr_r) *            &
              &              dvr_primitive_val_esteban(b, para, grid, curr_r_prime) *      &
              &              dvr_primitive_val_esteban(d, para, grid, curr_r_prime) *      &
              &              ( min(curr_r, curr_r_prime)**(l_val) /                &
              &                max(curr_r, curr_r_prime)**(l_val+1) ) * dr * dr
              !write(*,*) integral_val, dvr_primitive_val(a, para, grid, curr_r),   &
              !& dvr_primitive_val(b, para, grid, curr_r_prime), min(curr_r, curr_r_prime),&
              !& max(curr_r, curr_r_prime)
            end do
          end do
          if (abs(integral_val) < 1d-16) cycle
          two_e_rad_int_num(a,b,c,d,l) = integral_val
          write(13,'(4I4, ES25.17)') a, b, c, d, integral_val
!         write(*,*) 'Writing integral', a, b, c, d
        end do
        end do
      end do
      end do
      close(13)

    end do
    
    l_val= 0

    open(14,file="GLL_numerical_integrals_l"//trim(int2str(l_val))//".dat",  &
    &    form="formatted", action="write", recl=100000)

    open(15,file="GLL_numerical_integrals_normalized_l"//trim(int2str(l_val))//".dat",  &
    &    form="formatted", action="write", recl=100000)

    do a = 1, 4
    do b = 1,4
      do c = 1, 4
      do d = 1, 4
        integral_val = zero
        integral_val_1 = zero
        do i = 1, para%nr-2
          curr_r = grid%r(i)
          do j = 1, para%nr-2
            curr_r_prime = grid%r(j) 
            if (curr_r + curr_r_prime < 1d-16) cycle
            integral_val = integral_val +                                        &
            &              dvr_primitive_val_esteban(a, para, grid, curr_r) *            &
            &              dvr_primitive_val_esteban(c, para, grid, curr_r) *            &
            &              dvr_primitive_val_esteban(b, para, grid, curr_r_prime) *      &
            &              dvr_primitive_val_esteban(d, para, grid, curr_r_prime) *      &
            &              ( min(curr_r, curr_r_prime)**(l_val) /                &
            &                max(curr_r, curr_r_prime)**(l_val+1) ) *            &
            &              grid%weights(i) * grid%weights(j)
            integral_val_1 = integral_val_1 +                                        &
            &              dvr_primitive_val_esteban(a, para, grid, curr_r) *            &
            &              dvr_primitive_val_esteban(c, para, grid, curr_r) *            &
            &              dvr_primitive_val_esteban(b, para, grid, curr_r_prime) *      &
            &              dvr_primitive_val_esteban(d, para, grid, curr_r_prime) *      &
            &              ( min(curr_r, curr_r_prime)**(l_val) /                &
            &                max(curr_r, curr_r_prime)**(l_val+1) )
            !write(*,*) integral_val, dvr_primitive_val(a, para, grid, curr_r),   &
            !& dvr_primitive_val(b, para, grid, curr_r_prime), min(curr_r, curr_r_prime),&
            !& max(curr_r, curr_r_prime)
          end do
        end do
!       if (abs(integral_val) < 1d-16) integral_val = zero !truncate num. zeroes
        if (abs(integral_val) < 1d-16) cycle
        write(14,'(4I4, ES25.17)') a, b, c, d, integral_val
        write(15,'(4I4, 2ES25.17)') a, b, c, d, integral_val_1
!       write(*,*) 'Writing GLL integral', a, b, c, d
      end do
      end do
    end do
    end do

    close(13)
    close(14)
    close(15)


! Calculate the overlap between basis functions

  ! open(15,file="basis_overlap_l"//trim(int2str(l_val))//".dat",      &
  ! &    form="formatted", action="write", recl=100000)

  ! do a = 1, 10
  !   do b = 1, 10
  !     integral_val_1 = zero
  !     integral_val_2 = zero
  !     do i = 1, para%nr-2
  !       curr_r = grid%r(i)
! !       do j = 1, para%nr-2
! !         curr_r_prime = grid%r(j) 
! !         if (curr_r + curr_r_prime < 1d-16) cycle
  !         integral_val_1 = integral_val_1 +                                        &
  !         &              dvr_primitive_val_esteban(a, para, grid, curr_r) *    &
  !         &              dvr_primitive_val_esteban(b, para, grid, curr_r)*    &
  !         &              grid%weights(i)
  !         integral_val_2 = integral_val_2 +                                        &
  !         &              dvr_primitive_val_esteban(a, para, grid, curr_r) *    &
  !         &              dvr_primitive_val_esteban(b, para, grid, curr_r)
  !      end do
  !      write(15,'(2I4, ES25.17)') a, b, integral_val_1
  !    end do
  !  end do
  !  close(15)

  end subroutine radial_check

  subroutine check_with_analytic()

    integer :: norb, l_max, lp_max
    integer :: a, b, c, d, la, lb, lc, ld, lac, lbd, l, l_val, i, j, n_points

    real(dp) :: dr, curr_r, curr_rp, fact_1, fact_2, val_a, val_b, val_c, val_c_xc, val_d, val_d_xc, nrm, tol
    real(dp), allocatable :: an(:,:,:), int_val_dr(:), int_val_xc(:)
    real(dp), allocatable :: hydrogen_wf(:,:,:), TwoEInts(:,:,:,:)

    
    n_points = 1001
    dr = 0.1d0

    norb = orb%n_max

    l_max = para%l + 1
    lp_max = 2*para%l + 1

    tol = 1e-12

    call init_values(hydrogen_wf, norb, l_max, n_points, dr)

    allocate(two_e_rad_int_ana_dr(norb,norb,norb,norb,l_max,l_max,lp_max))
    allocate(two_e_rad_int_ana_xc(norb,norb,norb,norb,l_max,l_max,lp_max))

    allocate(int_val_dr(lp_max))
    allocate(int_val_xc(lp_max))

    two_e_rad_int_ana_dr = zero
    two_e_rad_int_ana_xc = zero


!   do a = 1,1 
!   do la = l_max, l_max
!   do i = 1, n_points
!     write(89,'(3i4,f16.8)') a, la, i, hydrogen_wf(i, a, la)
!   end do
!   end do
!   end do

    do a = 1, norb
      do c = 1, norb
!       if(c.ne.a) cycle
        do lac = 1, l_max
          do b = 1, norb
            do d = 1, norb
!             if(d.ne.b) cycle
              do lbd = 1, l_max
                int_val_dr = zero
                int_val_xc = zero
                do i = 1, n_points
                  curr_r = (i-1)*dr
                  val_a = hydrogen_wf(i, a, lac)
                  val_c = hydrogen_wf(i, c, lac)
                  val_c_xc = hydrogen_wf(i, d, lbd)
                  !val_c_xc = 1.0d0
!                 write(92,'(3i4,f16.8)')  a, b, i, val_a
                  do j = 1, n_points
                    curr_rp = (j-1)*dr
                    if (curr_r + curr_rp < 1d-16) cycle
                    fact_1 = dr*dr*(curr_r**2)*(curr_rp**2)
                    val_b = hydrogen_wf(j, b, lbd)
                    val_d = hydrogen_wf(j, d, lbd)
                    val_d_xc = hydrogen_wf(j, c, lac)
                    !val_d_xc = 1.0d0
!                   if (i.eq.1) write(92,'(3i4,f16.8)')  a, b, j, val_b

                    do l = 1, lp_max
                      l_val = l-1
                      fact_2 = (min(curr_r,curr_rp)**l_val)/(max(curr_r,curr_rp)**(l_val+1)) 
                      int_val_dr(l)=int_val_dr(l)+ fact_1*fact_2*val_a*val_b*val_c*val_d
                     !val_a = 1.0d0
                     !val_b = 1.0d0
                     !val_c_xc = 1.0d0
                     !val_d_xc = 1.0d0
                      int_val_xc(l)=int_val_xc(l)+ fact_1*fact_2*val_a*val_b*val_c_xc*val_d_xc
                    end do

                  end do
                end do   
                two_e_rad_int_ana_dr(a,b,c,d,lac,lbd,:) = int_val_dr(:)
                two_e_rad_int_ana_xc(a,b,c,d,lac,lbd,:) = int_val_xc(:)
              end do
            end do 
          end do
        end do
      end do
    end do

    do l = 1, lp_max
    do a = 1, norb
      do c = 1, norb
!       if(c.ne.a) cycle
        do lac = 1, l_max
          do b = 1, norb
            do d = 1, norb
!             if (d.ne.b) cycle
              do lbd = 1, l_max
                  write(91,'(7i3,f16.8)') a, b, c, d, lac, lbd, l, two_e_rad_int_ana_dr(a,b,c,d,lac,lbd,l)
                  write(92,'(7i3,f16.8)') a, b, c, d, lac, lbd, l, two_e_rad_int_ana_xc(a,b,c,d,lac,lbd,l)
!                 write(91,'(2i4,f16.8)') a, b, two_e_rad_int_ana_dr(a,b,c,d,lac,lbd,l)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    
    call CombineOrbInts_H(two_e_rad_int_ana_dr, two_e_rad_int_ana_xc, TwoEInts)

    call WriteInts_H(TwoEInts, tol)

    deallocate(int_val_dr, int_val_xc)
    deallocate(two_e_rad_int_ana_dr, two_e_rad_int_ana_xc)
    deallocate(hydrogen_wf)

  end subroutine check_with_analytic


  subroutine init_values(hydrogen_wf, norb, l_max, n_points, dr_p)

    real(dp), intent(in) :: dr_p
    integer, intent(in) :: n_points
    real(dp), allocatable, intent(inout) :: hydrogen_wf(:,:,:)

    integer :: norb, l_max
    integer :: a, la, k, i, n
    real(dp), allocatable :: an(:)
    real(dp) :: dr, drho, val, nrm, curr_r, curr_rho, norm

    allocate(hydrogen_wf(n_points,norb,l_max))

    do a = 1, norb
      do la = 1, l_max
!     do la = l_max, l_max
        n = a+la-1
        if (allocated(an)) deallocate(an)
        allocate(an(0:n-la))
        an = zero
        an(0) = one
        do k = 0, n-la-1
          an(k+1) = an(k)*(real(k+la-n,kind=idp) /  &
                       &             real((k + 1) * (k + 2*la ),kind=idp))
        end do
        norm = zero
        do i = 1, 10001
          val = zero
          dr = 0.01_idp
          drho = (two / real(n, kind=idp)) * dr 
          curr_r = (i-1) * dr
          curr_rho = (i-1) * drho
          do k = 0, n-la
            val = val + an(k) * curr_rho**k *exp(-curr_rho/two)
          end do
          val = curr_rho**(la-1) * val
          norm = norm + (val**2 * dr * curr_r**2)
        end do
        do i  = 1, n_points
          curr_r = (i-1)*dr_p
          hydrogen_wf(i,a,la) = radial_hydrogen_ef(curr_r, n, la, an, norm)
        end do
      end do
    end do
    
    deallocate(an)

  end subroutine init_values


  !! @description: Returns value of radial hydrogen eigenfunction at given point 
  !! @param: r_val Value of radial coordinate 
  !! @param: n     Principal quantum number
  !! @param: l     Orbital angular momentum quantum number
  !! @param: an    Coefficient array (should be unallocated on first call) 
  real(idp) function radial_hydrogen_ef(r_val, n, l, an, norm)
  
    real(idp), intent(in) :: r_val
    integer, intent(in)   :: n
    integer, intent(in)   :: l
    real(idp), allocatable, intent(in) :: an(:)
    real(idp), intent(in) :: norm

!   real(idp), allocatable :: an(:)
  
    integer :: i, k
    real(idp) :: val, dr, curr_r, drho, curr_rho
    
!   integer, save :: prev_n = -1, prev_l = -1
!   real(idp), save :: norm = - one


!   if ((n .ne. prev_n) .or. (l .ne. prev_l)) then
!     prev_n = n
!     prev_l = l
!     norm = zero
!     if (allocated(an)) deallocate(an)
!     allocate(an(0:n-l-1))
!     an(0) = one
!     do k = 0, n-l-2
!       an(k+1) = an(k) * (real(k + l + 1 - n,kind=idp) /                      &
!       &                  real((k + 1) * (k + 2*l + 2),kind=idp))
!     end do
!     do i = 1, 10001
!       val = zero
!       dr = 0.01_idp
!       drho = (two / real(n, kind=idp)) * dr 
!       curr_r = (i-1) * dr
!       curr_rho = (i-1) * drho
!       do k = 0, n-l-1
!         val = val + an(k) * curr_rho**k *exp(-curr_rho/two)
!       end do
!       val = curr_rho**l * val
!       norm = norm + (val**2 * dr * curr_r**2)
!     end do
!   end if

    curr_rho = (two / real(n, kind=idp)) * r_val
    val = zero
    do k = 0, n-l
      val = val + an(k) * curr_rho**k *exp(-curr_rho/two)
    end do
    val = curr_rho**(l-1) * val

    radial_hydrogen_ef = val / sqrt(norm)
  
  end function radial_hydrogen_ef

  real(idp) function dvr_primitive_val(ind, para, grid, r)

    integer,                    intent(in) :: ind
    type(para_t),               intent(in) :: para
    type(grid_t),               intent(in) :: grid 
    real(idp),                  intent(in) :: r

    integer :: m_val
    integer :: i_val

    i_val = ind / para%nl
    m_val = mod(ind, para%nl)

    !write(*,*) ind, i_val, m_val

    if (m_val == 0) then
      !Bridge function
      !write(*,*) 'BRIDGE'
      dvr_primitive_val = (dvr_product(r, i_val - 1, para%nl - 1, para, grid) +&
      &                    dvr_product(r, i_val, 0, para, grid)) /         &
      &                    sqrt(grid%weights(ind) + grid%weights(ind+1))
    else
      !write(*,*) 'STUFF'
      dvr_primitive_val = (dvr_product(r, i_val, m_val, para, grid) /          &
      &                   sqrt(grid%weights(ind)))
    end if

  end function dvr_primitive_val
  
  !! @description: Calculates DVR eigenfunction for arbitrary point on real axis 
  !! @param: ind   Index of DVR eigenfunction
  !! @param: para  Parameter data
  !! @param: grid  Spatial grid
  !! @param: r     Evaluation point
  real(idp) function dvr_primitive_val_esteban(ind, para, grid, r)

    integer,                    intent(in) :: ind
    type(para_t),               intent(in) :: para
    type(grid_t),               intent(in) :: grid 
    real(idp),                  intent(in) :: r

    integer :: m_val
    integer :: i_val
    integer :: i
    real(idp) :: r_min, r_max
    real(idp) :: xi, xi_ref, xi_curr
    real(idp) :: numerator, denominator
    
    i_val = ind / para%nl
    m_val = mod(ind, para%nl)
    if (i_val == 0) then
      r_min = zero
    else
      r_min = grid%r(i_val*para%nl)
    end if
    r_max = grid%r(i_val*para%nl+para%nl)
    
    if (r < r_min) then
      dvr_primitive_val_esteban = zero
      return 
    end if
    if (r > r_max) then
      dvr_primitive_val_esteban = zero
      return
    end if

    xi     = two * ((r           - r_min) / (r_max - r_min)) - one!map to [-1,1]
    xi_ref = two * ((grid%r(ind) - r_min) / (r_max - r_min)) - one!map to [-1,1]

    if (abs(xi) > one) then
      dvr_primitive_val_esteban = zero
      return 
    end if

    if (abs(xi - xi_ref) < 1d-16) then
      dvr_primitive_val_esteban = one
      return
    end if

    numerator = DerivLegPol(xi,para%nl) * (one - xi**2)
    denominator = real(para%nl*(para%nl+1),kind=8)                             &
    &             *LegPol(xi_ref,para%nl)*(xi-xi_ref)
    dvr_primitive_val_esteban = - numerator / (denominator ) 

  end function dvr_primitive_val_esteban
   
  real(idp) function LegPol(point, nl)

    implicit none
    real(idp), intent(in)    :: point
    integer,   intent(in)    :: nl 
   
    real(idp) :: Leg(0:nl)
    real(idp) :: fac0, fac1 
    integer   :: j 

    Leg(0) = one 
    Leg(1) = point 
    do j = 1, nl-1 
      fac1 = (two*real(j,kind = 8) + one) / (real(j,kind=8) + one) 
      fac0 = real(j,kind = 8)  / (real(j,kind=8) + one) 
      Leg(j+1) = point *  Leg(j) * fac1 - Leg(j-1) * fac0  
    enddo
    LegPol = Leg(nl)

  end function LegPol 
   
  real(idp) function DerivLegPol(point, nl)

    implicit none
    real(idp), intent(in)    :: point
    integer,   intent(in)    :: nl 
  
    real(idp) :: DLeg(0:nl)
    real(idp) :: fac0, fac1 
    integer   :: j 
   
    DLeg(0) = zero 
    DLeg(1) = one 
    do j = 1, nl-1 
      fac1 = (real(2*j,kind = 8) + one) / real(j,kind=8) 
      fac0 = real(j+1,kind = 8)  / (real(j,kind=8) ) 
      DLeg(j+1) = point *  DLeg(j) * fac1 - DLeg(j-1) * fac0  
    enddo
    DerivLegPol = DLeg(nl)

  end function DerivLegPol 

  real(idp) function dvr_product(r, i_val, m_val, para, grid)

    real(idp),                  intent(in) :: r
    integer,                    intent(in) :: i_val
    integer,                    intent(in) :: m_val
    type(para_t),               intent(in) :: para 
    type(grid_t),               intent(in) :: grid 

    integer :: k, ind
    real(idp) :: r_min, r_max

    ind = i_val * para%nl + m_val
    r_min = grid%r(i_val*para%nl)
    r_max = grid%r(i_val*para%nl+para%nl-1)

    if (r < r_min) then
      dvr_product = zero
      return 
    end if
    if (r > r_max) then
      dvr_product = zero
      return
    end if
    
    !write(*,*) i_val*para%nl, i_val*para%nl+para%nl
    !write(*,*) r, r_min, r_max

    dvr_product = one
    do k = 0, para%nl - 1
      if (k == m_val) cycle
      !write(*,*) r-grid%r(i_val*para%nl+k), grid%r(ind), grid%r(i_val*para%nl+k)
      !write(*,*) i_val*para%nl+k, ind, i_val, m_val
      dvr_product = dvr_product * ((r - grid%r(i_val*para%nl+k)) /             &
      &                            (grid%r(ind) - grid%r(i_val*para%nl+k)))
    end do

  end function dvr_product

  !! @description: Legendre function
  !! @param: j  Positive integer
  !! @param: x  Argument with |x|<=1
  real(idp) function legendre(j,x)

    integer   :: j
    real(idp) :: x

    real(idp) :: pmmp1, pmm, fact, somx2, pll
    integer   :: m, ji, ll, kk

    m = 0

    ji = j
    if (m < 0 .or. m > ji .or. abs(x) > one) then
      write(*,*) 'Bad arguments'
    end if
    pmm = one
    if (m > 0) then
      somx2 = sqrt((one - x) * (one + x))
      fact = one
      do kk = 1,m
        pmm = -pmm * fact * somx2
        fact = fact + two
      end do
    end if
    if (ji == m) then
      legendre = pmm
    else
      pmmp1 = x * (two * real(m,idp) + one) * pmm
      if (ji == m + 1) then
        legendre = pmmp1
      else
        do ll = m+2,ji
          pll = (x * (two * real(ll,idp) - one) * pmmp1 - &
          &     real(ll + m - 1,idp) * pmm) / real(ll - m,idp)
          pmm = pmmp1
          pmmp1 = pll
        end do
        legendre = pll
      end if
    end if

  end function legendre
  
  !! @description: Derivative of Legendre function
  !! @param: j  Positive integer
  !! @param: x  Argument with |x|<=1
  real(idp) function legendre_deriv(j, x)

    integer   :: j
    real(idp) :: x

    if (j == 0) then
      legendre_deriv = 0
      return
    end if

    legendre_deriv = (real(j, kind=idp) / (x**2 - 1)) *                        &
    &                (x * legendre(j,x) - legendre(j-1,x))


  end function legendre_deriv

  subroutine CombineOrbInts_H(TwoERadOrbInts_dr, TwoERadOrbInts_xc, TwoEInts)

    use DVRData, only : integrals_ang, sph_harm, para
    use OrbData, only : orb, SpatialOrbInd

    real(dp), allocatable, intent(in)  :: TwoERadOrbInts_dr(:,:,:,:,:,:,:)
    real(dp), allocatable, intent(in)  :: TwoERadOrbInts_xc(:,:,:,:,:,:,:)
    real(dp), allocatable, intent(inout) :: TwoEInts(:,:,:,:)

    integer  :: n_l, n_mp, l1, l2, l3, l4, m1, m2, m3, m4, k1, k2, k3, k4, l, indx
    integer  :: la, lb, lc, ld, ma, mb, mc, md, n_m1, n_m2, n_m3, n_m4, error
    integer  :: m1_init, m2_init, m3_init, m4_init, lma, lmb, lmc, lmd
    integer  :: klm_1, klm_2, klm_3, klm_4
    real(dp) :: start, finish, int_value_dr, int_value_xc


    call cpu_time(start)

    allocate(TwoEInts(orb%nSpatialOrbs,orb%nSpatialOrbs, orb%nSpatialOrbs, &
    & orb%nSpatialOrbs), stat=error)
    call allocerror(error)

    TwoEInts = zero

    n_l  = sph_harm%n_l
    n_mp = sph_harm%n_mp

    indx = 0
    do l1 = 1, n_l
      la = l1 - 1
      n_m1 = 2*l1 - 1
      m1_init = -1*l1

        do l2 = 1, n_l
          lb = l2 - 1
          n_m2 = 2*l2 - 1
          m2_init = -1*l2

            do m1 = 1, n_m1
              ma = m1_init + m1
              lma = (l1-1)**2 + m1
              do m3 = 1, n_m1
                mc = m1_init + m3
                lmc = (l1-1)**2 + m3

                do m2 = 1, n_m2
                  mb = m2_init + m2
                  lmb = (l2-1)**2 + m2
                  do m4 = 1, n_m2
                    md = m2_init + m4
                    lmd = (l2-1)**2 + m4

                    do k1 = 1, orb%n_max - l1 + 1
                      do k2 = 1, orb%n_max -l2 + 1  
                        do k3 = 1, orb%n_max - l1 + 1
                          do k4 = 1, orb%n_max -l2 + 1  

                            indx = indx + 1
                            int_value_dr = 0.0d0
                            int_value_xc = 0.0d0

                            klm_1 = SpatialOrbInd(k1,l1,m1)
                            klm_2 = SpatialOrbInd(k2,l2,m2)
                            klm_3 = SpatialOrbInd(k3,l1,m3)
                            klm_4 = SpatialOrbInd(k4,l2,m4)

                            do l = 1, 2*para%l + 1

                              int_value_dr = int_value_dr + (integrals_ang(l, lma, lmb, lmc, lmd)* &
                              &  TwoERadOrbInts_dr(k1,k2,k3,k4,l1,l2,l))

                              write(85,'(7i3,3f13.8)') k1,k2,k3,k4,l1,l2,l,integrals_ang(l, lma, lmb, lmc, lmd), TwoERadOrbInts_dr(k1,k2,k3,k4,l1,l2,l)!, int_value_dr

                              if (l1.ne.l2.or.m1.ne.m2) &
                              & int_value_xc = int_value_xc + (integrals_ang(l, lma, lmb, lmd, lmc)* &
                              &  TwoERadOrbInts_xc(k1,k2,k3,k4,l1,l2,l))

                              if (l1.ne.l2.or.m1.ne.m2) &
                              & write(87,'(11i3,3f13.8)') k1,k2,k3,k4,l1,l2,l, m1, m2, m3, m4, integrals_ang(l, lma, lmb, lmd, lmc), TwoERadOrbInts_xc(k1,k2,k3,k4,l1,l2,l), int_value_xc

!                             if (l1.ne.l2.or.m1.ne.m2) &
!                             & write(90,'(5i3,3f13.8)') klm_1, klm_2, klm_3, klm_4, l,integrals_ang(l, lma, lmb, lmd, lmc), TwoERadOrbInts_xc(k1,k2,k3,k4,l1,l2,l), int_value_xc

!                             if (klm_1.eq.1.and.klm_2.eq.1.and.klm_3.eq.1) then
!                               write(77,'(7I4,X,2F15.10)') k1, k2, lma, lmb, lmc, lmd, klm_4, integrals_ang(l, lma, lmb, lmc, lmd), TwoERadOrbInts(k1,k2,k3,k4,l1,l2,l)
!                             end if

                            end do
!                             write(85,'(6i5,f15.10)') k1,k2,k3,k4,l1,l2, int_value_dr
!                             write(86,'(6i5,f15.10)') k1,k2,k3,k4,l1,l2, int_value_xc
!                             if (klm_1.eq.1.and.klm_2.eq.1.and.klm_3.eq.1) then
!                               write(77,*) ''
!                             end if
                            TwoEInts(klm_1, klm_2, klm_3, klm_4) = int_value_dr
!                           write(77, '(6I5,X,f15.10)') k1, k2, lma, lmb, lmc, lmd, int_value
!                           write(78, '(4I5,X,f15.10)') klm_1, klm_2, klm_3, klm_4, TwoEInts(klm_1, klm_2, klm_3, klm_4)
                            if (l1.ne.l2.or.m1.ne.m2) then 
                              TwoEInts(klm_1, klm_2, klm_4, klm_3) = int_value_xc
                              TwoEInts(klm_4, klm_2, klm_1, klm_3) = int_value_xc
                              TwoEInts(klm_1, klm_3, klm_4, klm_2) = int_value_xc
                              TwoEInts(klm_4, klm_3, klm_1, klm_2) = int_value_xc
!                             write(79, '(4I5,X,f15.10)') klm_1, klm_2, klm_4, klm_3, TwoEInts(klm_1, klm_2, klm_4, klm_3)
                            end if
!                           write(79, '(4I5,X,f15.10)') lma, lmb, lmc, lmd, int_value
                          end do
                        end do
                      end do
                    end do
                  
                  end do
                end do  
              end do
            end do

        end do
    end do  

!   write(iout, *) 'Index:', indx

!   do l1 = 1, orb%nSpatialOrbs
!     do l2 = 1, orb%nSpatialOrbs
!       do l3 = 1, orb%nSpatialOrbs
!         do l4 = 1, orb%nSpatialOrbs
!!          write(80,'(4i5, f20.12)')  l1, l2, l3, l4, TwoEInts(l1,l2,l3,l4)
!           if (abs(TwoEInts(l1,l2,l3,l4)).gt.0.0d0) then
!!            write(81,*)  l1, l2, l3, l4, TwoEInts(l1,l2,l3,l4)
!           end if
!         end do
!       end do
!     end do
!   end do

    call cpu_time(finish)

    write(iout,'(X,a,f10.5,X,a)') 'Time taken for combining integrals = ', finish-start, 'seconds.'
  end subroutine CombineOrbInts_H

  subroutine WriteInts_H(TwoEInts, tol)
    
    use DVRData, only : para

    real(dp), intent(in) :: tol
    real(dp), allocatable, intent(in) :: TwoEInts(:,:,:,:)

    integer :: i, j, k, l, f_int, norbs, ij, kl, ijkl

    f_int = 15
    open(f_int, file='FCIDUMP_H', status='unknown', form="formatted")

    norbs = orb%nSpatialOrbs

    write(f_int,1001) norbs, nint(para%z), 0

1001 format(' &FCI NORB=', i5, ',NELEC=', i3, ',MS2=', i3,',')

    write(f_int,1002, advance='no') 

1002 format('  ORBSYM=')
 
    do i = 1, norbs

      write(f_int, '(i1)', advance='no') 1

      if (i.lt.norbs) then
        write(f_int,1003, advance='no')
      else
        write(f_int,1004) 
      endif

    end do
    
1003 format(',')
1004 format(',')

    write(f_int,*) ' ISYM=1,'
    write(f_int, *) ' &END'

    write(82,'(f16.8)') TwoEInts

    ij  = 0
    ijkl = 0
    do i = 1, norbs
      do j = 1, i
        kl = 0
!       do k = 1, norbs
        do k = 1, i
          do l = 1, k
            if (ij.ge.kl) then
              if (abs(TwoEInts(i,k,j,l)).gt.tol) &
              & write(f_int, 1005) TwoEInts(i,k,j,l), i, j, k, l
              ijkl = ijkl + 1
            end if
            kl = kl + 1
          end do
        end do
        ij = ij + 1
      end do
    end do

 !  do i = 1, norbs
 !    do j = 1, i
 !      if (abs(OneEInts(i,j)).gt.tol) &
!!      if (i.eq.j) &
 !      & write(f_int, 1005) OneEInts(i,j), i, j, 0, 0
 !    end do 
 !  end do

 !  write(f_int, 1005) 0.0_dp, 0, 0, 0, 0

1005 format(f20.16,x,5i4)


    close(f_int)
  end subroutine WriteInts_H

end module RadCheck

