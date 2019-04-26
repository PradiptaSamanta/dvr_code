module DensityPrim

  use constants
  use util_mod

  implicit none

  contains

  subroutine CalcDensityPrim(EigVecs)

    use DensityData
    use DVRData, only : para, grid
    use OrbData, only : break_inn
    use OrbInts, only : SetUpEigVec, SetUpEigVecBreak
  
    real(dp), allocatable, intent(in) :: EigVecs(:,:,:)

    real(dp), allocatable :: MOCoeff(:,:)
    complex(idp), allocatable :: CheckDens(:)
    real(dp), allocatable :: EigVecs_mod(:,:,:)
    integer, allocatable :: OrbInd_cntr(:,:,:), OrbInd_prim(:,:,:), n_m(:), get_l(:), get_m(:)
    integer :: n_l, n_cntr, i, j, k, error, tot_prim, n_prim, start_prim, tot_cntr
    integer :: n_m1, n_m2, ni
    real(dp) :: start, finish, value_r
    logical :: tOldForm

    tOldForm = .true.

    call cpu_time(start)

    ! upper limit for the number of l quantum number
    n_l = para%l + 1

    allocate(n_m(n_l))

    n_m1 = 2*para%ml_max + 1
    do i = 1, n_l
      n_m2 = 2*i - 1
      n_m(i) = min(n_m1, n_m2)
    end do

    ! First, the number of n quantum number which will be involved in calculating the photonization yield (n_prim) is calculated.
    ! Photonization yied is calculated by integrating electron density outside a certian spherical region. That is why
    ! the initial grid points are not needed while calculating the yield
    start_prim = 1
    n_prim = para%ng

!   if (abs(r_val).gt.1e-8) then
!     call CalcStartOrb(grid%r, para%ng, r_val, n_prim, start_prim)
!   end if

    write(iout, *) 'Here: ', n_prim
    ! The array that gives the orbital index in the primitive basis is allocated here
    allocate(OrbInd_prim(n_prim, n_l,2*n_l-1), stat=error)
    call allocerror(error)
    OrbInd_prim = 0

    ! Get the total number of orbitals involved in the primitive basis (tot_prim) and
    ! an array that gives the index of a primitive orbital for given, n, l and m

    tot_prim = 0

    do i = 1, n_prim
      do j = 1, n_l
        do k = 1, n_m(j)
          tot_prim = tot_prim + 1
          OrbInd_prim(i,j,k) = tot_prim
        end do
      end do
    end do

    allocate(get_l(tot_prim))
    allocate(get_m(tot_prim))

    do i = 1, n_prim
      do j = 1, n_l
        do k = 1, n_m(j)
          ni = OrbInd_prim(i,j,k)
          get_l(ni) = j
          get_m(ni) = k
        end do
      end do
    end do

    ! First read one and two electron RDMs obtained from a real-time (FCIQMC) simulation
   !if (tAvRDM) then
   !  call ReadAvOneRDM(DensOrb1e, file_1rdm, tot_cntr, nReadRDMs, MoCoeff, tot_cntr, tot_prim)
   !  !call ReadAvTwoRDM(DensOrb2e, file_2rdm, tot_cntr, nReadRDMs, tBinaryRDM)
   !else
      call ReadOneRDM(DensOrb1e, file_1rdm, tot_prim)
      !call ReadTwoRDM(DensOrb2e, file_2rdm, tot_cntr)
   !end if


    ! Transform one and two electron RDMs from the MO orbitals basis to the basis of the primitive orbitals

    call Get1eYieldPrim(DensOrb1e, Yield1e, n_prim, n_l, n_m,  OrbInd_prim, grid%r, r_val)

    !call Get2eYield(PrimDens2e, Yield2e, tot_prim)

    deallocate(OrbInd_prim)
    deallocate(DensOrb1e, get_l)
    !deallocate(DensOrb1e, DensOrb2e, get_l)

    call cpu_time(finish)
    write(iout,'(a, f10.5)') ' Time taken to calculate the Density: ', finish - start

  end subroutine CalcDensityPrim

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
    write(iout, *) 'debugging:', tot_orb, n_elements
    flush(6)

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
      if (i.ge.j) then
        ij = i*(i-1)/2 + j 
        Dens1e(ij) = dcmplx(val, zero)
        !write(80, '(2i5, g25.17)') i, j, real(Dens1e(ij))
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
      if (i.ge.j) then
        ij = i*(i-1)/2 + j 
        Dens1e(ij) = dcmplx(real(Dens1e(ij)), val)
        !write(81, '(2i5, f25.17)') i, j, val
      else
        cycle
      endif
    end do

    close(12)

  end subroutine ReadOneRDM

  subroutine ReadAvOneRDM(Dens1e, file_1, tot_orb, nFiles, MOCoeff, tot_cntr, tot_prim)

    complex(idp), allocatable, intent(inout) :: Dens1e(:)
    real(dp), allocatable, intent(in) :: MOCoeff(:,:)
    character(len=32), intent(in) :: file_1
    integer, intent(in) :: tot_orb, nFiles
    integer, intent(in) :: tot_cntr, tot_prim

    complex(idp), allocatable :: DensTemp(:,:)
    complex(idp), allocatable :: Dens2(:)
    real(dp) :: Yield
    complex(idp) :: csum
    integer :: n_elements, i, j, error, ij, iFile, i_f
    real(dp) :: val, start, finish
    character(len=32) :: filename
    logical :: file_exists

    call cpu_time(start)

    n_elements = tot_orb*(tot_orb + 1)/2

    allocate(Dens1e(n_elements), DensTemp(n_elements, nFiles))

    Dens1e = czero
    DensTemp = czero

    do iFile = 1, nFiles
    ! First read the real part of the one body reduced density matrix
      i_f = 2*iFile - 1
      filename = trim(file_1)//int2str(i_f)
 
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
        if (i.ge.j) then
          ij = i*(i-1)/2 + j 
          DensTemp(ij, iFile) = dcmplx(val, zero)
          !write(81, '(2i5, f25.17)') i, j, val
        else
          cycle
        endif
      end do
 
      close(12)
 
      ! Read the imaginary part of the one body reduced density matrix
      i_f = 2*iFile
      filename = trim(file_1)//int2str(i_f)
 
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
        if (i.ge.j) then
          ij = i*(i-1)/2 + j 
          DensTemp(ij, iFile) = dcmplx(real(DensTemp(ij, iFile)), val)
          !write(81, '(2i5, f25.17)') i, j, val
        else
          cycle
        endif
      end do
 
      close(12)

    end do

   !do i = 1, nFiles
   !  Dens1e(:) = DensTemp(:,i)
   !  call TransformDens1e(Dens1e, Dens2, MOCoeff, tot_cntr, tot_prim)
   !  call Get1eYield(Dens2, Yield, tot_prim)
   !  write(81, *) i, Yield
   !end do

    Dens1e = 0.0d0
    do i = 1, n_elements
        csum = sum(DensTemp(i,:))
        if (abs(csum).gt.1e-12) Dens1e(i) = csum/nFiles
    end do

    deallocate(DensTemp)

    write(iout, *) 'One electron RDM is averaged out'
    call cpu_time(finish)
    write(iout, *) 'Time taken ... ', finish-start

  end subroutine ReadAvOneRDM

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
      if (i.ge.j.and.k.ge.l) then
        ij = i*(i-1)/2 + j 
        kl = k*(k-1)/2 + l
        Dens2e(ij,kl) = dcmplx(val, zero)
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
      if (i.ge.j.and.k.ge.l) then
        ij = i*(i-1)/2 + j 
        kl = k*(k-1)/2 + l
        Dens2e(ij,kl) = dcmplx(real(Dens2e(ij,kl)), val)
        write(84, '(6i5, 2f25.17)') i, j, k, l, ij, kl, real(Dens2e(ij,kl)), aimag(Dens2e(ij,kl)) 
      else
        cycle
      endif
    end do

    close(12)

  end subroutine ReadTwoRDM

  subroutine ReadAvTwoRDM(Dens2e, file_1, tot_orb, nFiles, tBinaryRDM)

    complex(idp), allocatable, intent(inout) :: Dens2e(:,:)
    character(len=32), intent(in) :: file_1
    integer, intent(in) :: tot_orb, nFiles
    logical, intent(in) :: tBinaryRDM

    complex(idp), allocatable :: DensTemp(:,:,:)
    complex(idp) :: csum
    integer :: n_elements, i, j, k, l, error, ij, kl, iFile, i_f
    real(dp) :: val
    real(dp) :: check(tot_orb), trace, start, finish, t1, t2
    character(len=32) :: filename
    logical :: file_exists

    call cpu_time(start)

    n_elements = tot_orb*(tot_orb + 1)/2

   !allocate(Dens2e(n_elements, n_elements), stat=error)
   !call allocerror(error)
    allocate(Dens2e(n_elements, n_elements))
    allocate(DensTemp(n_elements, n_elements, nFiles))

    Dens2e = czero
    DensTemp = czero


    do iFile = 1, nFiles
      ! First read the real part of the two body reduced density matrix
      i_f = 2*iFile - 1
      filename = trim(file_1)//int2str(i_f)
 
      write(iout, *)  'Reading the real part of the 2-body RDM from the file ', trim(filename)
      inquire(file=trim(filename), exist=file_exists)
 
      check = 0.0d0
      trace = 0.0d0
 
      call cpu_time(t1)

      if (tBinaryRDM) then
        if (file_exists) then
          open(12, file=trim(filename), form="unformatted",&
          &    action="read", recl=100000)
        else
          call stop_all('ReadAvTwoRDM', 'File for TwoRDM is not present for reading')
        end if
 
        do
          read(12, iostat=error) i, j, k, l, val
          if  (error < 0) exit 
          if  (i < 0) exit 
          !if ((i.eq.k).and.(j.eq.l)) check(i) = check(i) + val
          !if ((i.eq.k).and.(j.eq.l)) trace = trace + val
          if (i.ge.j.and.k.ge.l) then
            ij = i*(i-1)/2 + j 
            kl = k*(k-1)/2 + l
            DensTemp(ij,kl,iFile) = dcmplx(val, zero)
          else
            cycle
          endif
        end do
 
      else 

        if (file_exists) then
          open(12, file=trim(filename), form="formatted",&
          &    action="read", recl=100000)
        else
          call stop_all('ReadAvTwoRDM', 'File for TwoRDM is not present for reading')
        end if
 
        do
          read(12, *, iostat=error) i, j, k, l, val
          if  (error < 0) exit 
          if  (i < 0) exit 
          !if ((i.eq.k).and.(j.eq.l)) check(i) = check(i) + val
          !if ((i.eq.k).and.(j.eq.l)) trace = trace + val
          if (i.ge.j.and.k.ge.l) then
            ij = i*(i-1)/2 + j 
            kl = k*(k-1)/2 + l
            DensTemp(ij,kl,iFile) = dcmplx(val, zero)
          else
            cycle
          endif
        end do
 

      end if

      close(12)
      call cpu_time(t2)
      write(iout, *) 'Time 1 ... ', t2-t1
 
      !do i = 1, tot_orb
      !  write(77, '(2i6, g25.17)') i, i, check(i)
      !end do
 
     !write(iout, '(a, g25.17)') ' Trace of the 1-RDM read from the file:', sum(check)
     !write(iout, '(a, g25.17)') ' Trace of the 2-RDM read from the file:', trace
 
      ! Read the imaginary part of the one body reduced density matrix
      i_f = 2*iFile
      filename = trim(file_1)//int2str(i_f)
 
      write(iout, *)  'Reading the imaginary part of the 2-body RDM from the file ', trim(filename)
      inquire(file=trim(filename), exist=file_exists)
 
      if(tBinaryRDM) then 
        if (file_exists) then
          open(12, file=trim(filename), form="unformatted",&
          &    action="read", recl=100000)
        else
          call stop_all('ReadAvTwoRDM', 'File for TwoRDM is not present for reading')
        end if
  
        call cpu_time(t1)
        do
          read(12, iostat=error) i, j, k, l, val
          !read(12, *) i, j, k, l, val
          if  (error < 0) exit 
          if  (i < 0) exit 
          if (i.ge.j.and.k.ge.l) then
            ij = i*(i-1)/2 + j 
            kl = k*(k-1)/2 + l
            DensTemp(ij,kl,iFile) = dcmplx(real(DensTemp(ij,kl,iFile)), val)
            !write(84, '(4i5, f25.17)') i, j, k, l, val
          else
            cycle
          endif
        end do
      else
        if (file_exists) then
          open(12, file=trim(filename), form="formatted",&
          &    action="read", recl=100000)
        else
          call stop_all('ReadAvTwoRDM', 'File for TwoRDM is not present for reading')
        end if
  
        call cpu_time(t1)
        do
          read(12, *, iostat=error) i, j, k, l, val
          !read(12, *) i, j, k, l, val
          if  (error < 0) exit 
          if  (i < 0) exit 
          if (i.ge.j.and.k.ge.l) then
            ij = i*(i-1)/2 + j 
            kl = k*(k-1)/2 + l
            DensTemp(ij,kl,iFile) = dcmplx(real(DensTemp(ij,kl,iFile)), val)
            !write(84, '(4i5, f25.17)') i, j, k, l, val
          else
            cycle
          endif
        end do
      end if
 
      close(12)
      call cpu_time(t2)
      write(iout, *) 'Time 2 ... ', t2-t1
    end do

    call cpu_time(t1)
    do i = 1, n_elements
      do j = 1, i
        csum = sum(DensTemp(i,j,:))
        if (abs(csum).gt.1e-12) then 
          Dens2e(i,j) = csum/nFiles
          Dens2e(j,i) = dconjg(csum)/nFiles
        end if
      end do
    end do
    call cpu_time(t2)
    write(iout, *) 'Time 3 ... ', t2-t1

    deallocate(DensTemp)

    write(iout, *) 'Two electron RDM is averaged out'
    call cpu_time(finish)
    write(iout, *) 'Time taken ... ', finish-start

  end subroutine ReadAvTwoRDM

  subroutine Get1eYieldPrim(Dens, Yield, n_prim, n_l, n_m, OrbInd, grid_r, r_val)
    complex(idp), allocatable, intent(in) :: Dens(:)
    real(dp), intent(inout) :: Yield
    integer, intent(in) :: n_prim, n_l
    real(dp), allocatable, intent(in) :: grid_r(:)
    integer, intent(in), allocatable :: n_m(:), OrbInd(:,:,:)
    real(dp), intent(in) :: r_val

    integer :: i, p, q, pp, pq, l1, m1, l2, m2
    real(dp) :: Yield_check

    Yield = 0.0d0
    Yield_check = 0.0d0

    do i = 1, n_prim
      if (grid_r(i).le.r_val) cycle
      do l1 = 1, n_l
        do m1 = 1, n_m(l1)
          p = OrbInd(i, l1, m1)
          pp = p*(p-1)/2 + p
          do l2 = 1, n_l
            do m2 = 1, n_m(l2)
              q = OrbInd(i, l2, m2)

              if (q.ge.p) cycle
              pq = p*(p-1)/2 + q
              if (abs(Dens(pq)).gt.1e-12) then
                write(78, '(5i4,e25.17)') p, q, i, l2, m2, real(Dens(pq))
                Yield = Yield + 2.0*real(Dens(pq))
              end if
            end do
          end do
          Yield = Yield + real(Dens(pp))
          Yield_check = Yield_check + real(Dens(pp))
        end do
      end do
    end do

    write(iout, '(a, f25.17)') ' One electron photoionization yield: ', Yield

    write(iout, '(a, f25.17)') ' Contribution from the diagonal parts: ', Yield_check
    write(iout, '(a, f25.17)') ' Contribution from the off-diagonal parts: ', Yield - Yield_check


  end subroutine Get1eYieldPrim

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

end module DensityPrim
