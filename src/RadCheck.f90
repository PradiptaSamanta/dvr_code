module RadCheck

  use dvr_spline_mod
  use dvr_diag_mod
  use DVRData

  use omp_lib

  implicit none

contains

  subroutine radial_check()

    integer                    :: i, j, a, b, c, d, l, l_val
    real(idp)                  :: full_r_max, curr_r, curr_r_prime, integral_val, dr
    real(idp)                  :: integral_val_1,  integral_val_2
 
    l_val=0
 
    open(13,file="numerical_integrals_l"//trim(int2str(l_val))//".dat",      &
    &    form="formatted", action="write", recl=100000)

    open(14,file="GLL_numerical_integrals_l"//trim(int2str(l_val))//".dat",  &
    &    form="formatted", action="write", recl=100000)

    open(15,file="GLL_numerical_integrals_normalized_l"//trim(int2str(l_val))//".dat",  &
    &    form="formatted", action="write", recl=100000)

    dr = 0.0001d0  
    do a = 1, 2
    do b = 1, 2
      do c = 1, 2
      do d = 1, 2
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
        write(13,'(4I4, ES25.17)') a, b, c, d, integral_val
!       write(*,*) 'Writing integral', a, b, c, d
      end do
      end do
    end do
    end do
    
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

    open(15,file="basis_overlap_l"//trim(int2str(l_val))//".dat",      &
    &    form="formatted", action="write", recl=100000)

    do a = 1, 10
      do b = 1, 10
        integral_val_1 = zero
        integral_val_2 = zero
        do i = 1, para%nr-2
          curr_r = grid%r(i)
!         do j = 1, para%nr-2
!           curr_r_prime = grid%r(j) 
!           if (curr_r + curr_r_prime < 1d-16) cycle
            integral_val_1 = integral_val_1 +                                        &
            &              dvr_primitive_val_esteban(a, para, grid, curr_r) *    &
            &              dvr_primitive_val_esteban(b, para, grid, curr_r)*    &
            &              grid%weights(i)
            integral_val_2 = integral_val_2 +                                        &
            &              dvr_primitive_val_esteban(a, para, grid, curr_r) *    &
            &              dvr_primitive_val_esteban(b, para, grid, curr_r)
         end do
         write(15,'(2I4, ES25.17)') a, b, integral_val_1
       end do
     end do
     close(15)

  end subroutine radial_check

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

end module RadCheck

