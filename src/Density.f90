module Density

  use constants
  use util_mod

  implicit none

  contains

  subroutine CalcDensity(EigVecs)

    use DensityData
    use DVRData, only : para
  
    real(dp), allocatable, intent(in) :: EigVecs(:,:,:)

    real(dp), allocatable :: MOCoeff(:,:)
    integer, allocatable :: OrbInd(:,:,:)
    integer :: n_l, n_nqn, i, j, k, error, tot_prim
    real(dp) :: start, finish

    call cpu_time(start)

    n_l = para%l + 1
    n_nqn = sum(n_orb_den)
    if (n_nqn.lt.0) call stop_all('CalcDenstiy','Number of orbitals are put wrong in the input')

    allocate(OrbInd(para%ng, n_l,2*n_l-1), stat=error)
    call allocerror(error)

    tot_prim = 0
    do i = 1, para%ng
      do j = 1, min(n_l,i)
        do k = 1, 2*j-1
          tot_prim = tot_prim + 1
          OrbInd(i,j,k) = tot_prim
          if (i.eq.n_nqn) tot_orb = tot_prim
        end do
      end do
    end do

    ! First read one and two electron RDMs obtained from a real-time simulation
    call ReadOneRDM(DensOrb1e, file_1rdm, tot_orb)
    call ReadTwoRDM(DensOrb2e, file_2rdm, tot_orb)

    ! Transform one and two electron RDMs from the MO orbitals basis to the basis of the primitive orbitals
    call GetOrbCoeff(EigVecs, MOCoeff, tot_orb, tot_prim, para%ng, n_nqn, n_l, OrbInd)
    call TransformDens1e(DensOrb1e, PrimDens1e, MOCoeff, tot_orb, tot_prim)
    call TransformDens2e(DensOrb2e, PrimDens2e, MOCoeff, tot_orb, tot_prim)

    call Get1eYield(PrimDens1e, Yield1e, tot_prim)
    call Get2eYield(PrimDens2e, Yield2e, tot_prim)

    deallocate(DensOrb1e, DensOrb2e)

    call cpu_time(finish)
    write(iout,'(a, f10.5)') ' Time taken to calculate the Density: ', finish - start

  end subroutine CalcDensity

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
        write(81, '(2i5, f25.17)') i, j, val
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

  subroutine GetOrbCoeff(EigVecs, MOCoeff, tot_orb, tot_prim, ng, n_nqn, n_l, OrbInd)

    real(dp), allocatable, intent(in) :: EigVecs(:,:,:)
    real(dp), allocatable, intent(out) :: MOCoeff(:,:)
    integer, allocatable, intent(in) :: OrbInd(:,:,:)
    integer, intent(in) :: tot_orb, tot_prim, ng, n_nqn, n_l

    integer :: error, indx_1, indx_2, n, i, l_val, l, m
    real(dp) :: val

    allocate(MOCoeff(tot_orb,tot_prim), stat=error)
    call allocerror(error)

    indx_1 = 0
    indx_2 = 0
    do n = 1, n_nqn
      do i = 1, ng
        l_val =  min(i,n,n_l)
        do l = 1, l_val
          val = EigVecs(i,n,l)
          do m = 1, 2*l-1
            indx_1 = OrbInd(n,l,m)
            indx_2 = OrbInd(i,l,m)
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
