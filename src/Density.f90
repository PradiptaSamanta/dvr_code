module Density

  use constants
  use util_mod

  implicit none

  contains

  subroutine CalcDensity(EigVecs)

    use DensityData
    use DVRData, only : para, grid
  
    real(dp), allocatable, intent(in) :: EigVecs(:,:,:)

    real(dp), allocatable :: MOCoeff(:,:)
    integer, allocatable :: OrbInd_cntr(:,:,:), OrbInd_prim(:,:,:)
    integer :: n_l, n_cntr, i, j, k, error, tot_prim, n_prim, start_prim, tot_cntr
    real(dp) :: start, finish

    call cpu_time(start)

    ! upper limit for the number of l quantum number
    n_l = para%l + 1
    ! upper limit for the n quantum number in the contracted basis in which the RDMs will be read
    n_cntr = sum(n_orb_den)
    if (n_cntr.lt.0) call stop_all('CalcDenstiy','Number of orbitals are put wrong in the input')

    ! First, the number of n quantum number which will be involved in calculating the photonization yield (n_prim) is calculated.
    ! Photonization yied is calculated by integrating electron density outside a certian spherical region. That is why
    ! the initial grid points are not needed while calculating the yield
    if (abs(r_val).gt.1e-8) then
      call CalcStartOrb(grid%r, para%ng, r_val, n_prim, start_prim)
    else 
      start_prim = 1
      n_prim = para%ng
    end if

    ! The array that gives the orbital index in the primitive basis is allocated here
    allocate(OrbInd_cntr(n_cntr, n_l,2*n_l-1), stat=error)
    call allocerror(error)
    OrbInd_cntr = 0
    allocate(OrbInd_prim(n_prim, n_l,2*n_l-1), stat=error)
    call allocerror(error)
    OrbInd_prim = 0

    ! Get the total number of orbitals involved in the contracted basis (tot_cntr) and 
    ! an array that gives the index of a contraced orbital for given, n, l and m
    tot_cntr = 0
    do i = 1, n_cntr
      do j = 1, min(n_l,i)
        do k = 1, 2*j-1
          tot_cntr = tot_cntr + 1
          OrbInd_cntr(i,j,k) = tot_cntr
        end do
      end do
    end do

    ! Get the total number of orbitals involved in the primitive basis (tot_prim) and
    ! an array that gives the index of a primitive orbital for given, n, l and m
    tot_prim = 0
    do i = 1, n_prim
      !do j = 1, min(n_l,i)
      do j = 1, n_l
        do k = 1, 2*j-1
          tot_prim = tot_prim + 1
          OrbInd_prim(i,j,k) = tot_prim
        end do
      end do
    end do

    ! First read one and two electron RDMs obtained from a real-time (FCIQMC) simulation
    call ReadOneRDM(DensOrb1e, file_1rdm, tot_cntr)
    call ReadTwoRDM(DensOrb2e, file_2rdm, tot_cntr)

    ! Transform one and two electron RDMs from the MO orbitals basis to the basis of the primitive orbitals
    call GetOrbCoeff(EigVecs, MOCoeff, tot_cntr, tot_prim, start_prim, n_prim, n_cntr, n_l, OrbInd_cntr, OrbInd_prim)
    call TransformDens1e(DensOrb1e, PrimDens1e, MOCoeff, tot_cntr, tot_prim)
    call TransformDens2e(DensOrb2e, PrimDens2e, MOCoeff, tot_cntr, tot_prim)

    call Get1eYield(PrimDens1e, Yield1e, tot_prim)
    call Get2eYield(PrimDens2e, Yield2e, tot_prim)

    deallocate(OrbInd_cntr, OrbInd_prim)
    deallocate(DensOrb1e, DensOrb2e)

    call cpu_time(finish)
    write(iout,'(a, f10.5)') ' Time taken to calculate the Density: ', finish - start

  end subroutine CalcDensity

  ! The following subroutine cut the total array of grid points at 'r_val'
  ! It returns the number of grid points after 'r_val' (n_prim) and 
  ! the position of this particular grid point in the whole array (start_prim)
  subroutine CalcStartOrb(grid_r, ng, r_val, n_prim, start_prim)

    real(dp), allocatable, intent(in) :: grid_r(:)
    integer, intent(in) :: ng
    real(dp), intent(in) :: r_val
    integer, intent(out) :: n_prim, start_prim

    real(dp) :: r_diff, sign_diff
    integer :: i
    if (size(grid_r).ne.ng) call stop_all('CalcStartOrb', 'Somehow number of grid points does not match')
    
    do i = 1, ng
      r_diff = r_val - grid_r(i)
      sign_diff = sign(1.0d0, r_diff)
      if (sign_diff.eq.-1.0d0) then
        n_prim = ng - i + 1
        start_prim = i
        exit
      else 
        cycle
      end if
    end do

  end subroutine CalcStartOrb

  subroutine ReadOneRDM(Dens1e, file_1, tot_orb)

    complex(idp), allocatable, intent(inout) :: Dens1e(:)
    character(len=32), intent(in) :: file_1
    integer, intent(in) :: tot_orb

    integer :: n_elements, i, j, error, ij
    real(dp) :: val
    character(len=32) :: filename
    logical :: file_exists

    n_elements = tot_orb*(tot_orb + 1)/2

    allocate(Dens1e(n_elements))

    Dens1e = czero

    ! First read the real part of the one body reduced density matrix
    filename = trim(file_1)//"1"

    write(iout, *)  'Reading the real part of the 1-body RDM from the file ', trim(filename)
    inquire(file=trim(filename), exist=file_exists)

    if (file_exists) then
      open(12, file=trim(filename), form="formatted",&
      &    action="read", recl=100000)
    else
      call stop_all('ReadOneRDM', 'File for OneRDM is not present for reading')
    end if

    do
      read(12, *, iostat=error) i, j, val
      if  (error < 0) exit 
      if (i.le.j) then
        ij = i*(i-1)/2 + j 
        Dens1e(ij) = cmplx(val, zero)
        !write(81, '(2i5, f25.17)') i, j, val
      else
        cycle
      endif
    end do

    close(12)

    ! Read the imaginary part of the one body reduced density matrix
    filename = trim(file_1)//"2"

    write(iout, *)  'Reading the imaginary part of the 1-body RDM from the file ', trim(filename)
    inquire(file=trim(filename), exist=file_exists)

    if (file_exists) then
      open(12, file=trim(filename), form="formatted",&
      &    action="read", recl=100000)
    else
      call stop_all('ReadOneRDM', 'File for OneRDM is not present for reading')
    end if

    do
      read(12, *, iostat=error) i, j, val
      if  (error < 0) exit 
      if (i.le.j) then
        ij = i*(i-1)/2 + j 
        Dens1e(ij) = cmplx(real(Dens1e(ij)), val)
        !write(81, '(2i5, f25.17)') i, j, val
      else
        cycle
      endif
    end do

    close(12)

  end subroutine ReadOneRDM

  subroutine ReadTwoRDM(Dens2e, file_1, tot_orb)

    complex(idp), allocatable, intent(inout) :: Dens2e(:,:)
    character(len=32), intent(in) :: file_1
    integer, intent(in) :: tot_orb

    integer :: n_elements, i, j, k, l, error, ij, kl
    real(dp) :: val
    real(dp) :: check(tot_orb), trace
    character(len=32) :: filename
    logical :: file_exists

    n_elements = tot_orb*(tot_orb + 1)/2

   !allocate(Dens2e(n_elements, n_elements), stat=error)
   !call allocerror(error)
    allocate(Dens2e(n_elements, n_elements))

    Dens2e = czero

    ! First read the real part of the two body reduced density matrix
    filename = trim(file_1)//"1"

    write(iout, *)  'Reading the real part of the 2-body RDM from the file ', trim(filename)
    inquire(file=trim(filename), exist=file_exists)

    if (file_exists) then
      open(12, file=trim(filename), form="formatted",&
      &    action="read", recl=100000)
    else
      call stop_all('ReadTwoRDM', 'File for TwoRDM is not present for reading')
    end if

    check = 0.0d0
    trace = 0.0d0

    do
      read(12, *, iostat=error) i, j, k, l, val
      if  (error < 0) exit 
      if  (i < 0) exit 
      if ((i.eq.k).and.(j.eq.l)) check(i) = check(i) + val
      if ((i.eq.k).and.(j.eq.l)) trace = trace + val
      if (i.le.j.and.k.le.l) then
        ij = i*(i-1)/2 + j 
        kl = k*(l-1)/2 + l
        Dens2e(ij,kl) = cmplx(val, zero)
      else
        cycle
      endif
    end do

    close(12)

    !do i = 1, tot_orb
    !  write(77, '(2i6, g25.17)') i, i, check(i)
    !end do

    write(iout, '(a, g25.17)') ' Trace of the 1-RDM read from the file:', sum(check)
    write(iout, '(a, g25.17)') ' Trace of the 2-RDM read from the file:', trace

    ! Read the imaginary part of the one body reduced density matrix
    filename = trim(file_1)//"2"

    write(iout, *)  'Reading the imaginary part of the 2-body RDM from the file ', trim(filename)
    inquire(file=trim(filename), exist=file_exists)

    if (file_exists) then
      open(12, file=trim(filename), form="formatted",&
      &    action="read", recl=100000)
    else
      call stop_all('ReadTwoRDM', 'File for TwoRDM is not present for reading')
    end if

    do
      read(12, *, iostat=error) i, j, k, l, val
      if  (error < 0) exit 
      if  (i < 0) exit 
      if (i.le.j.and.k.le.l) then
        ij = i*(i-1)/2 + j 
        kl = k*(l-1)/2 + l
        Dens2e(ij,kl) = cmplx(real(Dens2e(ij,kl)), val)
        !write(84, '(4i5, f25.17)') i, j, k, l, val
      else
        cycle
      endif
    end do

    close(12)

  end subroutine ReadTwoRDM

  subroutine GetOrbCoeff(EigVecs, MOCoeff, tot_orb, tot_prim, ni, ng, n_cntr, n_l, OrbInd_cntr, OrbInd_prim)

    real(dp), allocatable, intent(in) :: EigVecs(:,:,:)
    real(dp), allocatable, intent(out) :: MOCoeff(:,:)
    integer, allocatable, intent(in) :: OrbInd_cntr(:,:,:), OrbInd_prim(:,:,:)
    integer, intent(in) :: tot_orb, tot_prim, ni, ng, n_cntr, n_l

    integer :: error, indx_1, indx_2, n, i, l, m, i_p
    real(dp) :: val

    allocate(MOCoeff(tot_orb,tot_prim), stat=error)
    call allocerror(error)

    indx_1 = 0
    indx_2 = 0
    do l = 1, n_l
      do n = l, n_cntr
        do i = 1, ng
        i_p = i + ni - 1
          val = EigVecs(i_p,n,l)
          do m = 1, 2*l-1
            indx_1 = OrbInd_cntr(n,l,m)
            indx_2 = OrbInd_prim(i,l,m)
            MOCoeff(indx_1, indx_2) = val
            !write(77,'(6i4,f15.8)') n, l, m, i, indx_1, indx_2, MOCoeff(indx_1,indx_2)
          end do
        end do
      end do
    end do

  end subroutine GetOrbCoeff

  subroutine TransformDens1e(Dens1, Dens2, MOCoeff, tot_orb, tot_prim)

    complex(idp), allocatable, intent(in) :: Dens1(:)
    complex(idp), allocatable, intent(out) :: Dens2(:)
    real(dp), allocatable, intent(in) :: MOCoeff(:,:)
    integer, intent(in) :: tot_orb, tot_prim

    integer :: error, p, q, pp, pq, k

    allocate(Dens2(tot_prim), stat=error)
    call allocerror(error)
    
    Dens2 = czero

    do k = 1, tot_prim
      do p = 1, tot_orb
        pp = p*(p-1)/2 + p
        do q = 1, p-1

          pq = p*(p-1)/2 + q
          Dens2(k) = Dens2(k) + 2 * MOCoeff(p,k) * MOCoeff(q,k) * real(Dens1(pq))
         !Dens2(k) = Dens2(k) + MOCoeff(p,k) * MOCoeff(q,k) * Dens1(pq) + &
         !                      MOCoeff(p,k) * MOCoeff(q,k) * dconjg(Dens1(pq))
          
        end do
        Dens2(k) = Dens2(k) + MOCoeff(p,k) * MOCoeff(p,k) * real(Dens1(pp))
        !Dens2(k) = Dens2(k) + MOCoeff(p,k) * MOCoeff(p,k) * Dens1(pp)
      end do
      write(78,'(i5,2f25.17)') k, Dens2(k)
    end do

  end subroutine TransformDens1e

  subroutine Get1eYield(Dens, Yield, tot_prim)
    complex(idp), allocatable, intent(in) :: Dens(:)
    real(dp), intent(out) :: Yield
    integer, intent(in) :: tot_prim

    integer :: i

    Yield = 0.0d0

    do i = 1, tot_prim
      Yield = Yield + Dens(i)
    end do

    write(iout, '(a, f25.17)') ' One electron photoionization yield: ', Yield

  end subroutine Get1eYield

  subroutine TransformDens2e(Dens1, Dens2, MOCoeff, tot_orb, tot_prim)

    complex(idp), allocatable, intent(in) :: Dens1(:,:)
    complex(idp), allocatable, intent(out) :: Dens2(:)
    real(dp), allocatable, intent(in) :: MOCoeff(:,:)
    integer, intent(in) :: tot_orb, tot_prim

    integer :: error, p, q, pq, k, l,  kl, r, s, rs
    integer :: dim_dens
    real :: value

    dim_dens = tot_prim*(tot_prim+1)/2

    allocate(Dens2(dim_dens), stat=error)
    call allocerror(error)
    
    Dens2 = czero

    do k = 1, tot_prim
      do l = 1, k
        kl = k*(k-1)/2 + l
        value = 0.0d0
        do r = 1, tot_orb
          do s = 1, r
            do p = 1, r
              do q = 1, p

                pq = p*(p-1)/2 + q
                rs = r*(r-1)/2 + s
                if (abs(Dens1(pq,rs)).gt.1e-12) then
                  if (p.eq.q.and.r.eq.s) then
                    if (q.eq.r) then
                      value = value + real(Dens1(pq,rs)) * MOCoeff(p,k)**2 * MOCoeff(p,l)**2
                      cycle
                    else
                      value = value + 2.0d0 * real(Dens1(pq,rs)) * MOCoeff(p,k) * MOCoeff(p,l) * MOCoeff(r,k) * MOCoeff(r,l)
                    end if
                  elseif (p.eq.r.and.q.eq.s) then
                      value = value + real(Dens1(pq,rs)) * ( (MOCoeff(p,k)**2 * MOCoeff(q,l)**2) + (MOCoeff(p,l)**2 * MOCoeff(q,k)**2) )
                  elseif (p.eq.s.and.q.eq.r) then
                      value = value + 2.0d0 * real(Dens1(pq,rs)) * MOCoeff(p,k) * MOCoeff(p,l) * MOCoeff(q,k) * MOCoeff(q,l)
                  else
                      value = value + 2.0d0 * real(Dens1(pq,rs)) * ( (MOCoeff(p,k)*MOCoeff(q,l)*MOCoeff(r,k)*MOCoeff(s,l)) + &
                                                                     (MOCoeff(p,l)*MOCoeff(q,k)*MOCoeff(r,l)*MOCoeff(s,k)) )
                  endif
                end if
            
              end do
            end do
          end do
        end do
        Dens2(kl) = cmplx(value, zero)

        if (value.gt.1e-12) write(79,'(2i5,2f25.17)') k, l, Dens2(kl)
      end do
    end do

  end subroutine TransformDens2e

  subroutine Get2eYield(Dens, Yield, tot_prim)
    complex(idp), allocatable, intent(in) :: Dens(:)
    real(dp), intent(out) :: Yield
    integer, intent(in) :: tot_prim

    integer :: i, j, ij

    Yield = 0.0d0

    do i = 1, tot_prim
      do j = 1, i
        ij = i*(i-1)/2 + j
        if (abs(Dens(ij)).gt.1e-12) then
          Yield = Yield + Dens(ij) 
        end if
      end do
    end do

    write(iout, '(a, f25.17)') ' Two electron photoionization yield: ', Yield

  end subroutine Get2eYield

end module Density
