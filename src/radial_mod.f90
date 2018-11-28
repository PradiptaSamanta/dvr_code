module radial_mod

  use dvr_spline_mod
  use dvr_diag_mod, only: int2str, allocerror, mat_banded_to_full
  use DVRData, only: grid_t, para_t

  implicit none

! public

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

! real(idp) function dvr_primitive_val(ind, para, grid, r)

!   integer,                    intent(in) :: ind
!   type(para_t),               intent(in) :: para
!   type(grid_t),               intent(in) :: grid 
!   real(idp),                  intent(in) :: r

!   integer :: m_val
!   integer :: i_val

!   i_val = ind / para%nl
!   m_val = mod(ind, para%nl)

!   !write(*,*) ind, i_val, m_val

!   if (m_val == 0) then
!     !Bridge function
!     !write(*,*) 'BRIDGE'
!     dvr_primitive_val = (dvr_product(r, i_val - 1, para%nl - 1, para, grid) +&
!     &                    dvr_product(r, i_val, 0, para, grid)) /         &
!     &                    sqrt(grid%weights(ind) + grid%weights(ind+1))
!   else
!     !write(*,*) 'STUFF'
!     dvr_primitive_val = (dvr_product(r, i_val, m_val, para, grid) /          &
!     &                   sqrt(grid%weights(ind)))
!   end if

! end function dvr_primitive_val
! 
! real(idp) function dvr_product(r, i_val, m_val, para, grid)

!   real(idp),                  intent(in) :: r
!   integer,                    intent(in) :: i_val
!   integer,                    intent(in) :: m_val
!   type(para_t),               intent(in) :: para 
!   type(grid_t),               intent(in) :: grid 

!   integer :: k, ind
!   real(idp) :: r_min, r_max

!   ind = i_val * para%nl + m_val
!   r_min = grid%r(i_val*para%nl)
!   r_max = grid%r(i_val*para%nl+para%nl-1)

!   if (r < r_min) then
!     dvr_product = zero
!     return 
!   end if
!   if (r > r_max) then
!     dvr_product = zero
!     return
!   end if
!   
!   !write(*,*) i_val*para%nl, i_val*para%nl+para%nl
!   !write(*,*) r, r_min, r_max

!   dvr_product = one
!   write(*,*) ind, i_val, m_val 
!   do k = 0, para%nl - 1
!     if (k == m_val) cycle
!     write(*,*) r, k, grid%r(i_val*para%nl+k) 
!     !write(*,*) r-grid%r(i_val*para%nl+k), grid%r(ind), grid%r(i_val*para%nl+k)
!     !write(*,*) i_val*para%nl+k, ind, i_val, m_val
!     dvr_product = dvr_product * ((r - grid%r(i_val*para%nl+k)) /             &
!     &                            (grid%r(ind) - grid%r(i_val*para%nl+k)))
!   end do

! end function dvr_product

end module radial_mod
