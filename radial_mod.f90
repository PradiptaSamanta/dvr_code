module radial_mod

  use dvr_spline_mod
  use dvr_diag_mod

  implicit none

  public

contains

  !! @description: Convert banded storage matrix to full matrix, assuming the
  !! layout indicated by the combination of `mode`, `kl`, and `ku`.
  !!
  !! * `mode='g'` (general matrix) or `mode='t'` (triangular matrix):
  !!   convert directly, assuming `full(i,j)` is stored in
  !!   `packed(ku+1+i-j,j)` for `max(1,j-ku) <= i <= min(m, j+kl)
  !! * `mode='h'` (Hermitian matrix) or `mode='s'` (symmetric matrix):
  !!   If `kl=0`, assume information in `packed` describes the upper triangle
  !!   of the matrix. In this case `full(i,j)` is stored in
  !!   `packed(kd+1+i-j,j)` for `max(1,j-ku) <= i <= j` and `full(j,i)` is
  !!   either the complex conjugate of `full(i,j)` (`'h'`) or identical to
  !!   `full(i,j)` (`'s'`).
  !!   If `ku=0`, assume information in `packed` describes the lower triangle of
  !!   the matrix. `full(i,j)` is stored in `packed(1+i-j,j)` for
  !!   `j <= i <= min(n, j+kd)`. Again `full(j,i)` is again completed as
  !!   Hermitian or symmetric
  !!
  !! @param: full              Storage for full matrix. If already allocated,
  !!                           must be of size $n \times n$.
  !! @param: banded            Banded data, for input.
  !! @param: n                 Dimension of full matrix.
  !! @param: kl                Number of sub-diagonals. Must be 0 for Hermitian
  !!                           matrices
  !! @param: ku                Number of super-diagonals
  subroutine mat_banded_to_full(full, banded, n, kl, ku)

    real(idp), allocatable,    intent(inout) :: full(:,:)
    real(idp),                 intent(in)    :: banded(1+kl+ku,n)
    integer,                   intent(in)    :: n
    integer,                   intent(in)    :: kl
    integer,                   intent(in)    :: ku

    integer :: i, j, error, kd

    if (allocated(full)) then
      if ((lbound(full, 1) /= 1) .or.  (ubound(full, 1) /= n) .or.             &
      &   (lbound(full, 2) /= 1) .or.  (ubound(full, 2) /= n)) then
        write(*,*) 'ERROR: Full matrix is allocated to the wrong size.'
      end if
    else
      allocate(full(n,n), stat=error)
      call allocerror(error)
    end if

    full = czero

    ! See formula in http://www.netlib.org/lapack/lug/node124.html
    do j = 1, n
      do i = max(1, j-ku), min(n, j+kl)
        full(i,j) = full(i,j) + banded(ku+1+i-j,j)
      end do
    end do

  end subroutine mat_banded_to_full

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

  character(len=converted_l) function int2str(i, format)

    integer,                    intent(in) :: i
    character(len=*), optional, intent(in) :: format

    if (present(format)) then
      write(int2str, format) i
    else
      write(int2str, '(I25)') i
    end if
    int2str = adjustl(int2str)

  end function int2str

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
    write(*,*) ind, i_val, m_val 
    do k = 0, para%nl - 1
      if (k == m_val) cycle
      write(*,*) r, k, grid%r(i_val*para%nl+k) 
      !write(*,*) r-grid%r(i_val*para%nl+k), grid%r(ind), grid%r(i_val*para%nl+k)
      !write(*,*) i_val*para%nl+k, ind, i_val, m_val
      dvr_product = dvr_product * ((r - grid%r(i_val*para%nl+k)) /             &
      &                            (grid%r(ind) - grid%r(i_val*para%nl+k)))
    end do

  end function dvr_product

end module radial_mod
