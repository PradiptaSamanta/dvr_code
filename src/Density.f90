module Density

  use constants
  use util_mod

  implicit none

  contains

  subroutine CalcDensity(EigVecs)

    use DensityData
    use DVRData, only : para, grid
    use OrbData, only : break_inn
    use OrbInts, only : SetUpEigVec, SetUpEigVecBreak
  
    real(dp), allocatable, intent(in) :: EigVecs(:,:,:)

    real(dp), allocatable :: MOCoeff(:,:)
    real(dp), allocatable :: EigVecs_mod(:,:,:)
    integer, allocatable :: OrbInd_cntr(:,:,:), OrbInd_prim(:,:,:), n_m(:), get_l(:), get_m(:)
    integer :: n_l, n_cntr, i, j, k, error, tot_prim, n_prim, start_prim, tot_cntr
    integer :: n_m1, n_m2, ni
    real(dp) :: start, finish

    call cpu_time(start)

    ! upper limit for the number of l quantum number
    n_l = para%l + 1

    allocate(n_m(n_l))

    n_m1 = 2*para%ml_max + 1
    do i = 1, n_l
      n_m2 = 2*i - 1
      n_m(i) = min(n_m1, n_m2)
    end do

    ! upper limit for the n quantum number in the contracted basis in which the RDMs will be read
    n_cntr = sum(n_orb_den)
    if (n_cntr.lt.0) call stop_all('CalcDenstiy','Number of orbitals are put wrong in the input')

    ! First, the number of n quantum number which will be involved in calculating the photonization yield (n_prim) is calculated.
    ! Photonization yied is calculated by integrating electron density outside a certian spherical region. That is why
    ! the initial grid points are not needed while calculating the yield
    start_prim = 1
    n_prim = para%ng

    if (abs(r_val).gt.1e-8) then
      call CalcStartOrb(grid%r, para%ng, r_val, n_prim, start_prim)
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
        do k = 1, n_m(j)
          tot_cntr = tot_cntr + 1
          OrbInd_cntr(i,j,k) = tot_cntr
        end do
      end do
    end do

    allocate(get_l(tot_cntr))
    allocate(get_m(tot_cntr))

    do i = 1, n_cntr
      do j = 1, min(n_l,i)
        do k = 1, n_m(j)
          ni = OrbInd_cntr(i,j,k)
          get_l(ni) = j
          get_m(ni) = k
        end do
      end do
    end do

    ! Get the total number of orbitals involved in the primitive basis (tot_prim) and
    ! an array that gives the index of a primitive orbital for given, n, l and m
    tot_prim = 0
    !do i = start_prim, para%ng
    do i = 1, n_prim
      !do j = 1, min(n_l,i)
      do j = 1, n_l
        do k = 1, n_m(j)
          tot_prim = tot_prim + 1
          OrbInd_prim(i,j,k) = tot_prim
        end do
      end do
    end do

    if (para%split_grid) then
      call SetUpEigVec(EigVecs, EigVecs_mod)
      call GetOrbCoeff(EigVecs_mod, MOCoeff, tot_cntr, tot_prim, start_prim, n_prim, n_cntr, n_l, n_m, OrbInd_cntr, OrbInd_prim)
    elseif (break_inn) then
        call SetUpEigVecBreak(EigVecs, EigVecs_mod) ! ** Not fully checked
        call GetOrbCoeff(EigVecs_mod, MOCoeff, tot_cntr, tot_prim, start_prim, n_prim, n_cntr, n_l, n_m, OrbInd_cntr, OrbInd_prim)
    else
      call GetOrbCoeff(EigVecs, MOCoeff, tot_cntr, tot_prim, start_prim, n_prim, n_cntr, n_l, n_m, OrbInd_cntr, OrbInd_prim)
    endif

    ! First read one and two electron RDMs obtained from a real-time (FCIQMC) simulation
    if (tAvRDM) then
      call ReadAvOneRDM(DensOrb1e, file_1rdm, tot_cntr, nReadRDMs, MoCoeff, tot_cntr, tot_prim)
      !call ReadAvTwoRDM(DensOrb2e, file_2rdm, tot_cntr, nReadRDMs, tBinaryRDM)
    else
      call ReadOneRDM(DensOrb1e, file_1rdm, tot_cntr)
      !call ReadTwoRDM(DensOrb2e, file_2rdm, tot_cntr)
    end if


    ! Transform one and two electron RDMs from the MO orbitals basis to the basis of the primitive orbitals

    if (tOneRDMDiag) then
      call TransformDens1e(DensOrb1e, PrimDens1e, MOCoeff, tot_cntr, tot_prim)
      call Get1eYield(PrimDens1e, Yield1e, tot_prim)
    else
      call TransformDens1e_alt(DensOrb1e, PrimDens1e, MOCoeff, tot_cntr, n_prim, get_l, get_m, OrbInd_prim)
      call Get1eYield(PrimDens1e, Yield1e, n_prim)
    end if

    !!call TransformDens2e(DensOrb2e, PrimDens2e, MOCoeff, tot_cntr, tot_prim)
    !call TransformDens2e_alt(DensOrb2e, PrimDens2e, MOCoeff, tot_cntr, tot_prim, OrbInd_prim, n_prim, n_l, n_m)

    !call Get2eYield(PrimDens2e, Yield2e, tot_prim)

    deallocate(OrbInd_cntr, OrbInd_prim)
    !deallocate(DensOrb1e, DensOrb2e, get_l)
    deallocate(DensOrb1e, get_l)

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

    do i = 1, nFiles
      Dens1e(:) = DensTemp(:,i)
      call TransformDens1e(Dens1e, Dens2, MOCoeff, tot_cntr, tot_prim)
      call Get1eYield(Dens2, Yield, tot_prim)
      write(81, *) i, Yield
    end do

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

  subroutine GetOrbCoeff(EigVecs, MOCoeff, tot_orb, tot_prim, ni, ng, n_cntr, n_l, n_m, OrbInd_cntr, OrbInd_prim)

    real(dp), allocatable, intent(in) :: EigVecs(:,:,:)
    real(dp), allocatable, intent(out) :: MOCoeff(:,:)
    integer, allocatable, intent(in) :: n_m(:), OrbInd_cntr(:,:,:), OrbInd_prim(:,:,:)
    integer, intent(in) :: tot_orb, tot_prim, ni, ng, n_cntr, n_l

    integer :: error, indx_1, indx_2, n, i, l, m, i_p
    real(dp) :: val

    allocate(MOCoeff(tot_orb,tot_prim), stat=error)
    call allocerror(error)

    MOCoeff = 0.0d0
    indx_1 = 0
    indx_2 = 0
    do n = 1, n_cntr
      do l = 1, min(n_l, n) 
        do i = 1, ng
          i_p = i + ni - 1
          val = EigVecs(i_p,n-l+1,l)
          do m = 1, n_m(l)
            indx_1 = OrbInd_cntr(n,l,m)
            indx_2 = OrbInd_prim(i,l,m)
            MOCoeff(indx_1, indx_2) = val
            !if (abs(val).gt.1e-12) write(77,'(6i4,f20.16)') n, l, m, i_p, indx_1, indx_2, MOCoeff(indx_1,indx_2)
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
    real(dp) :: start, finish, norm

    call cpu_time(start)

    write(iout, *) 'Transforming 1-RDM from the contracted basis to the primitive one'

    allocate(Dens2(tot_prim), stat=error)
    call allocerror(error)
 
    Dens2 = czero
    norm = zero

    do k = 1, tot_prim
      do p = 1, tot_orb
        if (abs(MOCoeff(p,k)).lt.1e-12) cycle
        pp = p*(p-1)/2 + p
        do q = 1, p-1

          pq = p*(p-1)/2 + q
          if (abs(Dens1(pq)).gt.1e-12) then
          Dens2(k) = Dens2(k) + 2 * MOCoeff(p,k) * MOCoeff(q,k) * real(Dens1(pq))
         !Dens2(k) = Dens2(k) + MOCoeff(p,k) * MOCoeff(q,k) * Dens1(pq) + &
         !                      MOCoeff(p,k) * MOCoeff(q,k) * dconjg(Dens1(pq))
          
          end if
        end do
        Dens2(k) = Dens2(k) + MOCoeff(p,k) * MOCoeff(p,k) * real(Dens1(pp))
      end do
      norm = norm + Dens2(k)
    end do

    call cpu_time(finish)
    write(iout, *) 'Time taken ... ', finish-start

  end subroutine TransformDens1e

  subroutine TransformDens1e_alt(Dens1, Dens2, MOCoeff, tot_orb, n_prim, get_l, get_m, OrbInd)

    complex(idp), allocatable, intent(in) :: Dens1(:)
    complex(idp), allocatable, intent(out) :: Dens2(:)
    real(dp), allocatable, intent(in) :: MOCoeff(:,:)
    integer, intent(in) :: tot_orb, n_prim
    integer, allocatable, intent(in) :: get_l(:), get_m(:), OrbInd(:,:,:)

    integer :: error, p, q, pp, pq, k
    integer :: k1, k2, l1, l2, m1, m2
    real(dp) :: start, finish, norm

    call cpu_time(start)

    write(iout, *) 'Transforming 1-RDM from the contracted basis to the primitive one'

    allocate(Dens2(n_prim), stat=error)
    call allocerror(error)
 
    Dens2 = czero
    norm = zero

    do k = 1, n_prim
      do p = 1, tot_orb
        if (abs(MOCoeff(p,k)).lt.1e-12) cycle
        pp = p*(p-1)/2 + p
        l1 = get_l(p)
        m1 = get_m(p)
        k1 = OrbInd(k, l1, m1)
        do q = 1, p-1

          pq = p*(p-1)/2 + q
          if (abs(Dens1(pq)).gt.1e-12) then

          l2 = get_l(q)
          m2 = get_m(q)
          k2 = OrbInd(k, l2, m2)
!         if (p.eq.28.and.q.eq.1) then
!           write(78, '(2i5, 4f20.12)') k1, k2, MOCoeff(p,k1), MOCoeff(q,k2), real(Dens1(pq)), 2*MOCoeff(p,k1) * MOCoeff(q,k2) * real(Dens1(pq))   
!         end if

          Dens2(k) = Dens2(k) + 2 * MOCoeff(p,k1) * MOCoeff(q,k1) * real(Dens1(pq))
          
          end if
        end do
        Dens2(k) = Dens2(k) + MOCoeff(p,k1) * MOCoeff(p,k1) * real(Dens1(pp))
      end do
    end do

    call cpu_time(finish)
    write(iout, *) 'Time taken ... ', finish-start

  end subroutine TransformDens1e_alt

  subroutine Get1eYield(Dens, Yield, tot_prim)
    complex(idp), allocatable, intent(in) :: Dens(:)
    real(dp), intent(inout) :: Yield
    integer, intent(in) :: tot_prim

    integer :: i

    Yield = 0.0d0

    do i = 1, tot_prim
      Yield = Yield + real(Dens(i))
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
    real(dp) :: value, start, finish

    call cpu_time(start)

    write(iout, *) 'Transforming 2-RDM from the contracted basis to the primitive one'

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
        Dens2(kl) = dcmplx(value, zero)

        !if (value.gt.1e-12) write(79,'(2i5,2f25.17)') k, l, Dens2(kl)
      end do
    end do

    call cpu_time(finish)
    write(iout, *) 'Time taken ... ', finish-start

  end subroutine TransformDens2e

  !subroutine TransformDens2e_alt(Dens1, Dens2, MOCoeff, tot_orb, tot_prim, OrbInd, n_prim, n_l, n_m, get_l)
  subroutine TransformDens2e_alt(Dens1, Dens2, MOCoeff, tot_orb, tot_prim, OrbInd, n_prim, n_l, n_m)

    complex(idp), allocatable, intent(in) :: Dens1(:,:)
    complex(idp), allocatable, intent(out) :: Dens2(:)
    real(dp), allocatable, intent(in) :: MOCoeff(:,:)
    integer, intent(in) :: tot_orb, tot_prim, n_prim, n_l
    integer, allocatable, intent(in) :: OrbInd(:,:,:), n_m(:)

    integer :: error, p, q, pq, k, l,  kl, r, s, rs, l1, l2, n1, n2, m1, m2
    integer :: dim_dens
    real(dp) :: value, start, finish
    real(dp), allocatable :: coeff_1(:), coeff_2(:)

    call cpu_time(start)

    write(iout, *) 'Transforming 2-RDM from the contracted basis to the primitive one'

    dim_dens = tot_prim*(tot_prim+1)/2

    allocate(Dens2(dim_dens), stat=error)
    call allocerror(error)
    allocate(coeff_1(tot_orb), stat=error)
    call allocerror(error)
    allocate(coeff_2(tot_orb), stat=error)
    call allocerror(error)
    
    Dens2 = czero

    do l1 = 1, n_l
      do l2 = 1, n_l
        do n1 = 1, n_prim
          do n2 = 1, n_prim
            do m1 = 1, n_m(l1)
              do m2 = 1, n_m(l2)
                k = OrbInd(n1, l1, m1)
                l = OrbInd(n2, l2, m2)
                if (l.gt.k) cycle
                kl = k*(k-1)/2 + l
                value = 0.0d0
                coeff_1(:) = MOCoeff(:,k) 
                coeff_2(:) = MOCoeff(:,l) 
                do r = 1, tot_orb
                  !if (get_l(r).ne.l1) cycle
                  if (abs(coeff_1(r)).lt.1e-12) cycle
                  do s = 1, r
                    !if (get_l(s).ne.l2) cycle
                    if (abs(coeff_2(s)).lt.1e-12) cycle
                    do p = 1, r
                      !if (get_l(p).ne.l1) cycle
                      if (abs(coeff_1(p)).lt.1e-12) cycle
                      do q = 1, p
                        !if (get_l(q).ne.l2) cycle
                        if (abs(coeff_2(q)).lt.1e-12) cycle
 
                        pq = p*(p-1)/2 + q
                        rs = r*(r-1)/2 + s
                        if (abs(Dens1(pq,rs)).gt.1e-12) then
                          if (p.eq.q.and.r.eq.s) then
                            if (q.eq.r) then
                         !    value = value + real(Dens1(pq,rs)) * coeff_1(p)**2 * coeff_2(p)**2
                              value = value + real(Dens1(pq,rs)) * MOCoeff(p,k)**2 * MOCoeff(p,l)**2
                              cycle
                            else
                         !    value = value + 2.0d0 * real(Dens1(pq,rs)) * coeff_1(p) * coeff_2(p) * coeff_1(r) * coeff_2(r)
                              value = value + 2.0d0 * real(Dens1(pq,rs)) * MOCoeff(p,k) * MOCoeff(p,l) * MOCoeff(r,k) * MOCoeff(r,l)
                            end if
                          elseif (p.eq.r.and.q.eq.s) then
                         !    value = value + real(Dens1(pq,rs)) * ( (coeff_1(p)**2 * coeff_2(q)**2) + (coeff_2(p)**2 * coeff_1(q)**2) )
                              value = value + real(Dens1(pq,rs)) * ( (MOCoeff(p,k)**2 * MOCoeff(q,l)**2) + (MOCoeff(p,l)**2 * MOCoeff(q,k)**2) )
                          elseif (p.eq.s.and.q.eq.r) then
                         !    value = value + 2.0d0 * real(Dens1(pq,rs)) * coeff_1(p) * coeff_2(p) * coeff_1(q) * coeff_2(q)
                              value = value + 2.0d0 * real(Dens1(pq,rs)) * MOCoeff(p,k) * MOCoeff(p,l) * MOCoeff(q,k) * MOCoeff(q,l)
                          else
                         !    value = value + 2.0d0 * real(Dens1(pq,rs)) * ( (coeff_1(p)*coeff_2(q)*coeff_1(r)*coeff_2(s)) + &
                         !                                                   (coeff_2(p)*coeff_1(q)*coeff_2(r)*coeff_1(s)) )
                              value = value + 2.0d0 * real(Dens1(pq,rs)) * ( (MOCoeff(p,k)*MOCoeff(q,l)*MOCoeff(r,k)*MOCoeff(s,l)) + &
                                                                             (MOCoeff(p,l)*MOCoeff(q,k)*MOCoeff(r,l)*MOCoeff(s,k)) )
                          endif
                        end if
                    
                      end do
                    end do
                  end do
                end do
                Dens2(kl) = dcmplx(value, zero)
              end do
            end do
          end do
        end do

        !if (value.gt.1e-12) write(79,'(2i5,2f25.17)') k, l, Dens2(kl)
      end do
    end do

    call cpu_time(finish)
    write(iout, *) 'Time taken ... ', finish-start
    deallocate(coeff_1, coeff_2)

  end subroutine TransformDens2e_alt

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
