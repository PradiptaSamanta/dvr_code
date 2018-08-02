module dvr_diag_mod

  use constants
  use DVRData, only: grid
  use dvr_spline_mod
  use util_mod
  
  implicit none

! public
  
contains

  !! @description: Construct the explicit `nl+1` $\times$ `nr` matrix (because
  !!               of banded storage!) that is the Hamiltonian for a single
  !!               surface, consisting of the potential for that surface
  !!               (ignoring imaginary potentials, hence 'real'), and the
  !!               kinetic operator.
  !! @param: matrix       On output, explicit Hamiltonian matrix
  !! @param: grid         Coordinate Grid based on Gauss-Legendre-Lobatto points
  !! @param: Tkin_carinal Work array 
  subroutine get_real_surf_matrix_cardinal(matrix, grid, pot, Tkin_cardinal)

    use DVRData, only : grid_t

    real(idp), allocatable, intent(inout) :: matrix(:,:)
    type(grid_t),           intent(in)    :: grid
    real(idp), intent(in)    :: pot(:) 
    real(idp), allocatable, intent(in)    :: Tkin_cardinal(:,:)

    integer                 :: nl, ku
    integer                 :: i, ir, error, j
    complex(idp)            :: cpulse
    real(idp)               :: rpulse

    ku = grid%nl

    if (allocated(matrix)) deallocate(matrix)
    allocate(matrix(ku+1, size(grid%r)), stat=error)
    call allocerror(error)
    ! kinetic operator
    do i = 1,size(Tkin_cardinal(:,1))
      do j = 1,size(Tkin_cardinal(1,:))
        matrix(i,j) = Tkin_cardinal(i,j)
      end do
    end do

    ku = grid%nl

    do ir = 1, size(pot)
      matrix(ku + 1,ir) = matrix(ku+1,ir) + pot(ir)
    end do

  end subroutine get_real_surf_matrix_cardinal
  
  !! @description: This subroutine returns the converged approximations to
  !! to the problem $\hat{A}_{red}z = \lambda z$ via a call to the `Arpack`
  !! routines `dsaupd` and `dseupd`.
  !! @param: matrix      Matrix that should be diagonalized, will be replaced
  !!                     by matrix of `nev` eigenvectors of dimension `n`
  !! @param: formt       Format of the input matrix to be diagonzalied. 'formt'
  !!                     *   `formt=full`: Full storage matrix format
  !!                     *   `format=banded`: Matrix stored in `Lapack` band
  !!                          storage format.
  !!                     could be one of the following:
  !! @param: n           Dimension of the eigenproblem
  !! @param: nev         Number of eigenvalues to be computed
  !! @param: which       Specifies which eigenvalues to be
  !!                     computed. `which` could be one of the following:
  !!                     *   `which=LA`: get the `nev` eigenvalues of largest
  !!                          amplitude.
  !!                     *   `which=SA': get the `nev` eigenvalues of smallest
  !!                          amplitude
  !!                     *   `which=BE`: Compute `nev` eigenvalues, half from
  !!                          each end of the spectrum.  When `nev` is odd,
  !!                          compute one more from the high end than from the
  !!                          low end.
  !! @param: eigenvals   Array of eigenvalues of the matrix containing the `nev`
  !!                     desired eigenvalues
  !! @param: rvec        Specifies wheather the associated `nev` eigenvector
  !!                     will be calculated or not.
  subroutine diag_arpack_real_sym_matrix(matrix, formt, n, nev, which,         &
  & eigenvals, rvec)

    real(idp), allocatable, intent(inout) :: matrix(:,:)
    character(len=*),       intent(in)    :: formt
    integer,                intent(in)    :: n
    integer,                intent(in)    :: nev
    character(len=2),       intent(in)    :: which
    real(idp),              intent(inout) :: eigenvals(:)
    logical,                intent(in)    :: rvec

    integer                :: ldv, mode, maxitr, lworkl
    integer                :: i, j, ido, error, iter
    real(idp)              :: sigma, tol
    character(len=1)       :: bmat
    integer                :: ncv, info
    integer                :: iparam(11)
    logical,   allocatable :: selct(:)
    real(idp), allocatable :: resid(:)
    real(idp), allocatable :: v(:,:)
    real(idp), allocatable :: workd(:)
    real(idp), allocatable :: workl(:)
    integer                :: ipntr(14)
    real(idp), allocatable :: d(:)

    ido    = 0
    info   = 0
    sigma  = zero
    tol    = zero
    ldv    = n
    ncv    = n
    lworkl = ncv**2 + 8*ncv
    maxitr = 900
    mode   = 1
    bmat   = 'I'
    iparam(1) = 1 !
    iparam(3) = maxitr
    iparam(4) = 1
    iparam(7) = mode

    allocate(v(ldv,ncv),stat=error)
    call allocerror(error)
    allocate(resid(n),stat=error)
    call allocerror(error)
    allocate(workd(3*n),stat=error)
    call allocerror(error)
    allocate(workl(lworkl),stat=error)
    call allocerror(error)
    allocate(selct(ncv),stat=error)
    call allocerror(error)
    allocate(d(nev),stat=error)
    call allocerror(error)

    iter = 0

    select case (formt)

      case ('full')

        do
          iter = iter + 1
          call dsaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,       &
          &           iparam, ipntr, workd, workl, lworkl, info)
          select case (abs(ido))
            case (1)
              ! compute  Y = OP * X  where
              ! IPNTR(1) is the pointer into WORKD for X,
              ! IPNTR(2) is the pointer into WORKD for Y.
              ! This is for the initialization phase to force the
              ! starting vector into the range of OP.
              call wrap_dsymv(n, one, matrix, workd(ipntr(1):ipntr(1)+n-1),    &
              &               zero, workd(ipntr(2):ipntr(2)+n-1))
              if (ido > 0) then       !no need to do this loop B=1
                ! compute Z = B * X for B = unity
                ! IPNTR(1) is the pointer into WORKD for X,
                ! IPNTR(3) is the pointer into WORKD for Z.
                do i = 1, n ! simply copy
                  workd(ipntr(3)+i-1) = workd(ipntr(1)+i-1)
                end do
              end if
            case (99)
              exit
            case default
              write(*,*) "ERROR: Unknown value for ido: ", ido
              stop
          end select
        end do

      case ('banded')

        do
          iter = iter + 1
          call dsaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,       &
          &           iparam, ipntr, workd, workl, lworkl, info)
          select case (abs(ido))
            case (1)
              ! compute  Y = OP * X  where
              ! IPNTR(1) is the pointer into WORKD for X,
              ! IPNTR(2) is the pointer into WORKD for Y.
              ! This is for the initialization phase to force the
              ! starting vector into the range of OP.
              call wrap_dsbmv(size(matrix(1,:)),size(matrix(:,1))-1, one,      &
              &               matrix, workd(ipntr(1):ipntr(1)+n-1), zero,      &
              &               workd(ipntr(2):ipntr(2)+n-1))
              if (ido > 0) then       !no need to do this loop B=1
                ! compute Z = B * X for B = unity
                ! IPNTR(1) is the pointer into WORKD for X,
                ! IPNTR(3) is the pointer into WORKD for Z.
                do i = 1, n ! simply copy
                  workd(ipntr(3)+i-1) = workd(ipntr(1)+i-1)
                end do
              end if
            case (99)
              exit
            case default
              write(*,*) "ERROR: Unknown value for ido: ", ido
              stop 
          end select
        end do

      case default

        write(*,*) "ERROR: Unknown format"
        stop

    end select ! formt

    if (allocated(matrix)) deallocate(matrix)
    allocate(matrix(n,nev), stat=error)
    call allocerror(error)

    do j = 1, nev
      do i = 1, n
        matrix(i,j) = v(i,j)
      end do
    end do

    call dseupd(rvec, 'A', selct, d, matrix, ldv, sigma, bmat, n, which, nev,  &
    &           tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl,  &
    &           info)

    if (.not. rvec) then
      do i = 1,nev
        eigenvals(i) = d(nev+1-i)
      end do
    else
      do i = 1,nev
        eigenvals(i) = d(i)
      end do
    end if

    deallocate(selct, d, v, workd)
    if (allocated(resid)) deallocate(resid)
    if (allocated(workl)) deallocate(workl)

  end subroutine diag_arpack_real_sym_matrix

  ! Symmetric matrix-vector multiplication
  !! @description: See BLAS
  !! [dsymv](http://www.netlib.org/blas/dsymv.f) documentation.
  !!
  !!  $\vec{y} = \alpha \hat{a} \vec{x} + \beta \vec{y}$
  !!
  !! @param: n      Dimension of `a`, `x`, `y`
  !! @param: alpha  Scalar for `x`
  !! @param: a      Matrix
  !! @param: x      Vector
  !! @param: beta   Scalar for `y`
  !! @param: y      Vector
  subroutine wrap_dsymv(n, alpha, a, x, beta, y)

    integer,       intent(in)    :: n
    real(idp),     intent(in)    :: alpha
    real(idp),     intent(in)    :: a(n,n)
    real(idp),     intent(in)    :: x(n)
    real(idp),     intent(in)    :: beta
    real(idp),     intent(inout) :: y(n)

    call dsymv('U', n, alpha, a, n, x, 1, beta, y, 1)

  end subroutine wrap_dsymv

  ! Symmetric banded matrix-vector multiplication
  !! @description: See BLAS
  !! [dsbmv](http://www.netlib.org/blas/dsbmv.f) documentation.
  !!
  !!  $\vec{y} = \alpha \hat{a} \vec{x} + \beta \vec{y}$
  !!
  !! @param: n      Dimension of `a`, `x`, `y`
  !! @param: ku     Number of super-diagonals
  !! @param: alpha  Scalar for `x`
  !! @param: a      Matrix (banded storage)
  !! @param: x      Vector
  !! @param: beta   Scalar for `y`
  !! @param: y      Vector
  subroutine wrap_dsbmv(n, ku, alpha, a, x, beta, y)

    integer,       intent(in)    :: n
    integer,       intent(in)    :: ku
    real(idp),     intent(in)    :: alpha
    real(idp),     intent(in)    :: a(ku+1,n)
    real(idp),     intent(in)    :: x(n)
    real(idp),     intent(in)    :: beta
    real(idp),     intent(inout) :: y(n)

    call dsbmv('U', n, ku, alpha, a, ku+1, x, 1, beta, y, 1)

  end subroutine wrap_dsbmv

end module dvr_diag_mod
