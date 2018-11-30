module radial_mod

  use dvr_spline_mod
  use dvr_diag_mod

  implicit none

  public

contains

  ! Symmetric matrix factorization
  !! @description: See BLAS
  !! [dsytrf](http://www.netlib.org/blas/dsytrf.f) documentation.
  !!
  !! @param: n      Dimension of `a`, `x`, `y`
  !! @param: a      Matrix
  !! @param: ipiv   Auxiliary integer array
  !! @param: work   Auxiliary double precision array
  subroutine wrap_dsytrf(a, n, ipiv, work)

    real(idp),     intent(in)    :: a(n,n)
    integer,       intent(in)    :: n
    integer,       intent(in)    :: ipiv(n)
    real(idp),     intent(in)    :: work(n)

    integer :: error

    call dsytrf('U', n, a, n, ipiv, work, n, error)

    if (error .ne. 0) then
      write(*,*) 'Error in Matrix inversion routine dsytri, error value: ',    &
      &          trim(int2str(error))
    end if

  end subroutine wrap_dsytrf

  ! Symmetric matrix inverse
  !! @description: See BLAS
  !! [dsytri](http://www.netlib.org/blas/dsytri.f) documentation.
  !!
  !! @param: n      Dimension of `a`, `x`, `y`
  !! @param: a      Matrix
  !! @param: ipiv   Auxiliary integer array
  !! @param: work   Auxiliary double precision array
  subroutine wrap_dsytri(a, n, ipiv, work)

    real(idp),     intent(in)    :: a(n,n)
    integer,       intent(in)    :: n
    integer,       intent(in)    :: ipiv(n)
    real(idp),     intent(in)    :: work(n)

    integer :: error

    call dsytri('U', n, a, n, ipiv, work, error)

    if (error .ne. 0) then
      write(*,*) 'Error in Matrix inversion routine dsytri, error value: ',    &
      &          trim(int2str(error))
    end if

  end subroutine wrap_dsytri
  
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

  !! @description: Calculates DVR eigenfunction for arbitrary point on real axis 
  !! @param: ind   Index of DVR eigenfunction
  !! @param: para  Parameter data
  !! @param: grid  Spatial grid
  !! @param: r     Evaluation point
  real(idp) function dvr_primitive_val(ind, para, grid, r)

    integer,                    intent(in) :: ind
    type(para_t),               intent(in) :: para
    type(grid_t),               intent(in) :: grid 
    real(idp),                  intent(in) :: r

    integer :: m_val
    integer :: i_val
    integer :: i
    real(idp) :: r_min, r_max
    real(idp) :: xi, xi_ref, xi_curr

    i_val = ind / para%nl
    m_val = mod(ind, para%nl)
    if (i_val == 0) then
      r_min = zero
    else
      r_min = grid%r(i_val*para%nl)
    end if
    r_max = grid%r(i_val*para%nl+para%nl)
    
    if (r < r_min) then
      dvr_primitive_val = zero
      return 
    end if
    if (r > r_max) then
      dvr_primitive_val = zero
      return
    end if

    xi     = two * ((r           - r_min) / (r_max - r_min)) - one!map to [-1,1]
    xi_ref = two * ((grid%r(ind) - r_min) / (r_max - r_min)) - one!map to [-1,1]

    if (abs(xi) > one) then
      dvr_primitive_val = zero
      return 
    end if

    dvr_primitive_val = - ( legendre_deriv(para%nl,xi) * (one - xi**2)) /      & 
    &             (real(para%nl*(para%nl+1),kind=8)                            &
    &             *legendre(para%nl,xi_ref)*(xi-xi_ref))

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

end module radial_mod
