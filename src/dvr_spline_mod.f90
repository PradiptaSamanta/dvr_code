module dvr_spline_mod

  use constants
  use DVRData

  implicit none

  public

  contains

  !! @description: Initialize the given spline for interpolation
  !!               of $y(x)$
  !! @param: spline  Spline structure
  !! @param: x       Array of $x$-coordinate values
  !! @param: y       Array of function values
  subroutine init_spline(spline, x, y)

    type(spline_t), intent(inout) :: spline
    real(idp),      intent(in)    :: x(:)
    real(idp),      intent(in)    :: y(:)

    integer :: ir, nr, error

    nr = size(x)
    if (size(y) /= nr) then
      write(*,*) "ERROR: x and y have to be of the same size."
    end if
    if (allocated(spline%cm)) then
      if (size(spline%cm) /= nr) deallocate(spline%cm)
    end if
    if (allocated(spline%x)) then
      if (size(spline%x) /= nr) deallocate(spline%x)
    end if
    if (allocated(spline%y)) then
      if (size(spline%y) /= nr) deallocate(spline%y)
    end if
    if (.not. allocated(spline%cm)) then
      allocate(spline%cm(nr), stat=error)
    end if
    if (.not. allocated(spline%x)) then
      allocate(spline%x(nr), stat=error)
    end if
    if (.not. allocated(spline%y)) then
      allocate(spline%y(nr), stat=error)
    end if
    do ir=1, nr
      spline%x(ir) = x(ir)
      spline%y(ir) = y(ir)
    end do
    spline%n = nr
    call spli(x, y, spline%cm)

  end subroutine init_spline


  !! @description: Initialize the given splines for the
  !!               interpolation of the complex function $y(x)$
  !! @param: rspline  Spline structure for the real part of $y$
  !! @param: ispline  Spline structure for the imaginary part of $y$
  !! @param: x        Array of $x$-coordinate values
  !! @param: y        Array of complex function values
  subroutine init_cspline(rspline, ispline, x, y)

    type(spline_t), intent(inout) :: rspline
    type(spline_t), intent(inout) :: ispline
    real(idp),      intent(in)    :: x(:)
    complex(idp),   intent(in)    :: y(:)

    real(idp), allocatable :: ry(:), iy(:)
    integer :: i

    allocate(ry(size(y)))
    allocate(iy(size(y)))
    do i=1, size(y)
      ry(i) = real(y(i), idp)
      iy(i) = aimag(y(i))
    end do
    call init_spline(rspline, x, ry)
    call init_spline(ispline, x, iy)
    deallocate(ry, iy)

  end subroutine init_cspline


  !! @description: Interpolate a function on $x$ using the given
  !!               `spline` structure. The interpolation result
  !!               is saved in `y`.
  !! @param: spline  Spline structure
  !! @param: x       Interpolation point
  !! @param: y       Interpolated value
  subroutine spline_value(spline, x, y)

    type(spline_t), intent(in)  :: spline
    real(idp),      intent(in)  :: x
    real(idp),      intent(out) :: y

    y = spl(spline%n, spline%x, spline%y, spline%cm, x)

  end subroutine spline_value


  !! @description: Interpolate a function of $x$ using the given
  !!               spline structure. The interpolation result is
  !!               saved in `y`.
  !! @param: rspline  Spline structure for real part
  !! @param: ispline  Spline structure for imaginary part
  !! @param: x        Interpolation point
  !! @param: y        Interpolated complex value
  subroutine cspline_value(rspline, ispline, x, y)

    type(spline_t), intent(in)  :: rspline
    type(spline_t), intent(in)  :: ispline
    real(idp),      intent(in)  :: x
    complex(idp),   intent(out) :: y

    y =         spl(rspline%n, rspline%x, rspline%y, rspline%cm, x) &
    &    + ci * spl(ispline%n, ispline%x, ispline%y, ispline%cm, x)

  end subroutine cspline_value


  !! @description: calculate the derivative $dy = f^\prime(x)$ as
  !!               $(f(x + \gamma) - f(x)) / \gamma$
  !!               where $\gamma$ is a small constant and $f(x)$ is a function
  !!               that is accessible through a spline interpolation
  !! @param: spline  Spline interpolation routine that can be used
  !!                 to calculate $f(x)$
  !! @param: x       A specific position $x$
  !! @param: dy      The derivative $f^\prime(x)$ that is to be calculated
  subroutine getderivative(spline, x, dy)

    type (spline_t),    intent(in)  :: spline
    real (idp),         intent(in)  :: x
    real (idp),         intent(out) :: dy

    real (idp) :: y, yy, gamma

    gamma = 3.0d-6

    y = spl(spline%n, spline%x, spline%y, spline%cm, x)
    yy = spl(spline%n, spline%x, spline%y, spline%cm, x+gamma)
    dy = (yy-y)/gamma

  end subroutine getderivative


  !-----------------------------------------------------------------------------
  !
  ! The routines spli and spl are implemented according to the conventional
  ! algorithm to programing natural Cubic Spline. The formalism is based upon
  ! the application of so called c0, c1 and c2 continuity conditions at data
  ! points. The third order interpolation polynomial commonly denoted as Si(x)
  ! (spline) in the interval "i", xi <= x <= x(i+1), is rewritten as a Taylor
  ! expansion such as:
  !
  !          Si(x) = ai + bi(x-xi) + ci(x-xi)**2 + di(x-xi)**3
  !
  ! This form is more suitable for computation purposes. For more details on the
  ! general algorithm, see:
  !
  !http://www.google.de/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&ved=
  !0CDgQFjAB&url=http%3A%2F%2Fgeophysics-old.tau.ac.il%2Fheb%2Fstudent%2Fexams%
  !2Ffilelist.asp%3Faction%3Ddownload%26dir%3D%255C%25EB%25EC%25EC%25E9%255C%
  !25F9%25E9%25E8%25E5%25FA%2B%25EE%25FA%25EE%25E8%25E9%25E5%25FA%2B%25E1%25E2%
  !25E9%25E0%25E5%25F4%25E9%25E6%25E9%25F7%25E4%255C%25F1%25E9%25EB%25E5%25EE%
  !25E9%25ED%2B%25E5%25E3%25F4%25E9%2B%25F0%25E5%25F1%25E7%25E0%25E5%25FA%
  !26file3Dcubic_spline_p2.pdf&ei=GtLAVP_ENojAOfGRgcAO&usg=
  !AFQjCNG6D2gTqFe6YRd5As2typw4kYPBEg&bvm=bv.83829542,d.ZWU
  !
  !-----------------------------------------------------------------------------
  !! @description: the routine spli initializes the spline calculation by
  !!               computing the interpolation help coefficients in terms of the
  !!               second derivatives of the spline S"(xi), usually written as
  !!               "zi" in the literature. They are calculated here recurrently,
  !!               after have been analytically solved for over Gaussian
  !!               elimination (see reference above).
  !! @param: x   Array of $x$-coordinate values
  !! @param: y   Array of function values
  !! @param: cm  Interpolation help coefficients "zi"
  subroutine spli(x, y, cm)

    real(idp), intent(in)    :: x(:)
    real(idp), intent(in)    :: y(:)
    real(idp), intent(inout) :: cm(:)

    integer :: i, n, error
    real(idp), allocatable :: dx(:), dy(:), diff_forward(:)
    real(idp), allocatable :: u(:), v(:)

    ! Note that the consistency check as for the matching of the dimension of
    ! the input arguments x and y has been already conducted in the routine
    ! init_spline in which the call to the present routine has been initiated.
    n = size(x)

    ! dx(i) = x(i+1) - x(i) : length of the interval "i", building an array
    !                         consisting of n-1 components because dx(n) cannot
    !                         be defined according to the discretization scheme
    allocate (dx(n-1), stat=error)

    ! dy(i) = y(i+1) - y(i) : variation of the function to be interpolated
    !                         at data points "i", evaluated conform to the
    !                         forward differentiation scheme, building an array
    !                         containing n-1 components because dy(n) is not
    !                         defined
    allocate (dy(n-1), stat=error)

    ! diff_forward(i) = dy(i)/dx(i) : forward differentiation, building
    !                                 an array of n-1 components
    allocate (diff_forward(n-1), stat=error)

    ! Help array of dimension n-2 obtained from the Gaussian elimination, not
    ! defined at the boundary data points, i=1 and i=n
    ! i = 2        : u(2) = 2*(dx(2) + dx(1)),
    ! i = 3,.., n-1: u(i) = 2*(dx(i) + dx(i-1)) - dx(i-1)**2/u(i-1)
    allocate (u(2:n-1), stat=error)

    ! Help array of dimension n-2 obtained from the Gaussian elimination, not
    ! defined at the boundary data points, i=1 and i=n
    ! i = 2        : v(2) = 6*(diff_forward(2) - diff_forward(1)),
    ! i = 3,.., n-1: v(i) = 6*(diff_forward(i) - diff_forward(i-1)) -
    !                       - dx(i-1)*(v(i-1)/u(i-1))
    allocate (v(2:n-1), stat=error)

    ! Calculating dx, dy, diff_forward using here the ability of fortran in
    ! performing compact array operation, once the dimension of the involved
    ! array matches. We have checked it!
    dx(1:n-1) = x(2:n) - x(1:n-1)
    dy(1:n-1) = y(2:n) - y(1:n-1)

    diff_forward(1:n-1) = dy(1:n-1)/dx(1:n-1) ! Or simply: diff_forward = dy/dx
                                              ! will also work

    ! Calculating the help array u(n-2) and v(n-2)
    u(2) = 2.0d0*(dx(2) + dx(1))
    v(2) = 6.0d0*(diff_forward(2) - diff_forward(1))

    do i = 3, n-1
      u(i) = 2.0d0*(dx(i) + dx(i-1)) - (dx(i-1))**2/u(i-1)
      v(i) = 6.0d0*(diff_forward(i) - diff_forward(i-1)) - &
      &      dx(i-1)*(v(i-1)/u(i-1))
    end do

    ! Calculating the interpolation help coefficients, "cm" <---> "zi", defined
    ! as the second derivative of the spline S"(xi) at the data points. Note
    ! that the dummy argument "cm" has been already properly allocated from the
    ! calling routine init_spline

    ! Upper boundary condition for the natural cubic spline, zn=0
    cm(n) = 0.0d0

    do i = n-1, 2, -1
      cm(i) = (v(i) - dx(i)*cm(i+1))/u(i)
    end do

    ! Lower boundary condition for the natural cubic spline, z1=0
    cm(1) = 0.0d0

    deallocate(dx, dy, diff_forward, u, v)

  end subroutine spli


  !! @description: the spl function calculates interpolated function values at
  !!                given interpolation points, where the arguments n, x, y, and
  !!                zi <--> cm are variables within the defined spline type
  !!                spline_t
  !! @param: n   Array dimension
  !! @param: x   Array of $x$-coordinate values
  !! @param: y   Array of function values
  !! @param: zi  Array of interpolation help coefficients, zi = S"(xi)
  !! @param: xp  Interpolation point
  real(idp) function spl(n, x, y, zi, xp) result(yp)

    integer,   intent(in) :: n
    real(idp), intent(in) :: x(n)
    real(idp), intent(in) :: y(n)
    real(idp), intent(in) :: zi(n)
    real(idp), intent(in) :: xp

    integer   :: i, caseindex
    real(idp) :: ai, bi, cci, di, dxi, dyi

    ! Case indexing
    if (xp < x(1)) then
      caseindex = 1
    elseif (xp > x(n)) then
      caseindex = 2
    else
      caseindex = 3
    end if

    ! Case selecting
    select case(caseindex)
    case(1)
      ! Interpolation point within the lower outer range
      dxi  = x(2) - x(1)
      dyi  = y(2) - y(1)
      yp = y(1) + (dyi/dxi - zi(2)*(dxi/6.0d0))*(xp - x(1)) ! Asymptotic
    case(2)
      ! Interpolation point within the upper outer range
      dxi = x(n) - x(n-1)
      dyi = y(n) - y(n-1)
      yp = y(n) + (dyi/dxi - zi(n-1)*(dxi/3.0d0))*(xp - x(n)) ! Asymptotic
    case(3)
      ! Interpolation point within the given $x$-coordinate data range
      i=1
      do while(xp > x(i+1))
        i = i + 1           ! Locating the interall within which xp lies
      end do
      dxi = x(i+1) - x(i)
      dyi = y(i+1) - y(i)
      ai  = y(i)
      bi  = -(dxi/6.0d0)*zi(i+1) - (dxi/3.0d0)*zi(i) + dyi/dxi
      cci = zi(i)/2.0d0
      di  = (zi(i+1) - zi(i))/(6.0d0*dxi)

      yp = ai + (xp - x(i))*(bi + (xp - x(i))*(cci + (xp - x(i))*di)) ! Result

    end select

    return

  end function spl

end module
