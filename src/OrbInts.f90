module OrbInts

  use constants
  use util_mod, only : stop_all, allocerror
  use ReadInput, only : n_max, two_e_int, nfrz, shift_int, reduce_int, red_start, red_num, n_shift_out_orb
  use OrbData
  use CombineInts, only : CombineOrbInts, Calc2eRadOrbInts

  implicit none

  contains

  subroutine SetOrbData()

    use DVRData, only : para, sph_harm
    use util_mod, only : allocerror

    integer :: i, j, k, l_val, indx, error

    if (para%split_grid) then
      if (n_max(2) == -1) n_max(2) = n_max(1)
      orb%n_inner = n_max(1)
      orb%n_outer = n_max(2)
      orb%n_max = sum(n_max)
      orb%shift_int = shift_int
      orb%n_shift_out_orb = n_shift_out_orb
    else 
      orb%n_max = n_max(1)
    end if

    orb%reduce_orb = reduce_int
    if (orb%reduce_orb) then
      orb%break = red_start
      orb%n_red = red_num
    end if
    
    allocate(SpatialOrbInd(orb%n_max,para%l+1,2*para%l+1),stat=error)
    call allocerror(error)

    SpatialOrbInd = 0

    indx = 0
    do i = 1, orb%n_max
      l_val = min(i,para%l+1)
!     write(*,*) 'limit for l_val:', l_val
      do j = 1, l_val
        !do k = 1, 2*j-1
        do k = 1, sph_harm%n_m(j)
          indx = indx + 1
          SpatialOrbInd(i-j+1,j,k) = indx
!         write(iout, *) i-j+1, j, k, SpatialOrbInd(i-j+1, j, k)
        end do
      end do
    end do

    orb%nSpatialOrbs = indx

    nFrozen = nfrz

    write(iout, *) '**********'
    write(iout, *) 'Setting up these parameters for the orbitals:'
    write(iout, '(X,A,3X, I6)') 'orb%n_max     =', orb%n_max
    if (para%split_grid) then
      write(iout, '(X,A,3X, I6)') 'orb%n_inner     =', orb%n_inner
      write(iout, '(X,A,3X, I6)') 'orb%n_outer     =', orb%n_outer
      write(iout, '(X,A,3X, I6)') 'orb%n_shift_out_orb     =', orb%n_shift_out_orb
    end if
    write(iout, '(X,A,3X, I6)') 'orb%nSpatialOrbs     =', indx
    write(iout, *) '***********' 

!   do i =  1, n_max
!     do j = 1, para%l+1
!       do k = 1, 2*para%l+1
!         write(76,'(4I5)') i, j, k, SpatialOrbInd(i,j,k)
!         if (SpatialOrbInd(i,j,k).ne.0) write(77,'(4I5)') i, j, k, SpatialOrbInd(i,j,k)
!       end do
!     end do
!   end do  
    
!   write(iout,*) 'Total number of spatial orbital', nSpatialOrbs

  end subroutine SetOrbData

! ***********************
  subroutine GetOrbInts()

    use DVRData, only : para, eigen_vecs
    use OrbData, only : break_inn, break_orb

    real(dp) :: tol

    real(dp), allocatable :: EigVecs(:,:,:)

    tol = 1e-12

    call SetOrbData()

    if (para%split_grid) then
        call SetUpEigVec(eigen_vecs, EigVecs)
    elseif (orb%reduce_orb) then
        call SetUpEigVecNorm(eigen_vecs, EigVecs) ! ** Not fully checked
    elseif (break_inn) then
        call SetUpEigVecBreak(eigen_vecs, EigVecs) ! ** Not fully checked
    end if

    ! Calculate the one-electron part of the Hamiltonian in the contracted basis
    if (para%split_grid) then
      call Calc1eOrbInts(EigVecs)
    else if (orb%reduce_orb) then
      call Calc1eOrbInts(EigVecs)
    else if (break_inn) then
      call Calc1eOrbInts(EigVecs)
    else
      call Calc1eOrbInts(eigen_vecs)
    end if

    ! Calculate the two-electron part of the Hamiltonian in the contracted basis
    ! First do the transformation from the FE-DVR basis to the contracted basis
    ! Then, Combined the radial and angular parts to get the total 2-e integrals
    if (para%split_grid) then
      call Calc2eRadOrbInts(EigVecs)
    else if (orb%reduce_orb) then
      call Calc2eRadOrbInts(EigVecs)
    else if (break_inn) then
      call Calc2eRadOrbInts(EigVecs)
    else
      call Calc2eRadOrbInts(eigen_vecs)
    end if
    call CombineOrbInts()
 
    ! Write down the integrals into FCIDUMP
    file_int = 'FCIDUMP'

    call WriteInts(tol)

  end subroutine GetOrbInts
! ***********************

  subroutine SetUpEigVec(VecsOld, EigVecs)

    use DVRData, only : para, grid

    real(dp), allocatable, intent(in) :: VecsOld(:,:,:)
    real(dp), allocatable, intent(inout) :: EigVecs(:,:,:)

    integer :: len_2, i, j, l, n_shift, n_l, len_1, j_p
    integer, allocatable :: len_mid(:)

    n_l = para%l + 1
    allocate(len_mid(n_l))
!    allocate(len_1(n_l))

    len_1 = para%m1*para%nl
    len_2 = para%m2*para%nl - 1

    n_shift = orb%n_shift_out_orb
    !n_shift = 0

    if (orb%shift_int) then
      len_mid = para%m1*para%nl + orb%n_inner + n_shift
      do l = 1, n_l
        len_mid(l) = len_mid(l) - l + 1
      end do
    else
      len_mid = para%m1*para%nl + n_shift
    end if


    write(iout, *) 'The integrals are now calculated for orbitals taken from two separated regions.'
    write(iout, '(a, i3)') ' Orbitals from the inner region: 1 -', orb%n_inner
    write(iout, '(a,i3,a,i3)') ' Orbitals from the outer region: ', len_mid(1)+1, ' -', len_mid(1)+ orb%n_outer
    if (.not.orb%shift_int) write(iout, *) '***Orbitals are not shifted in the outer region***'

!   do i = 1, size(grid%r)
!     write(78, '(11f10.6)') (VecsOld(i,j,3), j=1, size(VecsOld(1,:,2)))
!   end do

    allocate(EigVecs(size(grid%r),orb%n_max,n_l))
    EigVecs = 0.0d0
    do l = 1, n_l
      do j = 1, orb%n_inner - l + 1
        do i = 1, len_1
          EigVecs(i,j,l) = VecsOld(i,j,l)
        end do
      end do
    end do

    if (orb%n_outer == 0) return

    do l = 1, n_l
      if ((len_mid(l)+orb%n_outer+l-1).gt.size(grid%r))  &
      &  call stop_all('SetUpEigVec','Not enough orbitals in the outer region')
      do j = 1, orb%n_outer + l - 1
        j_p = j + orb%n_inner - l + 1
        do i = 1, len_2
          EigVecs(i+len_1,j_p,l) = VecsOld(i+len_1,j+len_mid(l),l)
          !EigVecs(i+len_1,j+orb%n_inner,:) = VecsOld(i+len_1,j+len_mid,:)
        end do
      end do
    end do

!   do i = 1, size(grid%r)
!     write(79, '(11f10.6)') (EigVecs(i,j,3), j=1, size(EigVecs(1,:,2)))
!   end do

  end subroutine SetUpEigVec

  subroutine SetUpEigVecBreak(VecsOld, EigVecs)

    use DVRData, only : para, grid

    real(dp), allocatable, intent(in) :: VecsOld(:,:,:)
    real(dp), allocatable, intent(inout) :: EigVecs(:,:,:)

    integer :: len_2, i, j, l, n_shift, n_l, len_1, j_p, tot_orb
    integer, allocatable :: len_mid(:)

    n_l = para%l + 1
    allocate(len_mid(n_l))

    len_mid = break_orb(1) + break_orb(2)
    do l = 1, n_l
      len_mid(l) = len_mid(l) - l + 1
    end do

    tot_orb = break_orb(1) + break_orb(3)

    write(iout, *) 'The integrals are now calculated after orbitals being distributed over two sets'
    write(iout, '(a, i3)') ' Orbitals in the first set: 1 -', break_orb(1)
    write(iout, '(a,i3,a,i3)') ' Orbitals in the second set: ', len_mid(1)+1, ' -', len_mid(1)+ break_orb(3)

!   do i = 1, size(grid%r)
!     write(78, '(11f10.6)') (VecsOld(i,j,3), j=1, size(VecsOld(1,:,2)))
!   end do

    allocate(EigVecs(size(grid%r), tot_orb, n_l))
    EigVecs = 0.0d0
    do l = 1, n_l
      do j = 1, break_orb(1) - l + 1
        do i = 1, size(grid%r)
          EigVecs(i,j,l) = VecsOld(i,j,l)
        end do
      end do
    end do

    do l = 1, n_l
      if ((len_mid(l)+break_orb(3)+l-1).gt.size(grid%r))  &
      &  call stop_all('SetUpEigVecBreak','Not enough orbitals in the outer region')
      do j = 1, break_orb(3) + l - 1
        j_p = j + break_orb(1) - l + 1
        do i = 1, size(grid%r)
          EigVecs(i,j_p,l) = VecsOld(i,j+len_mid(l),l)
          !EigVecs(i+len_1,j+orb%n_inner,:) = VecsOld(i+len_1,j+len_mid,:)
        end do
      end do
    end do

!   do i = 1, size(grid%r)
!     write(79, '(11f10.6)') (EigVecs(i,j,3), j=1, size(EigVecs(1,:,2)))
!   end do

  end subroutine SetUpEigVecBreak 

  subroutine SetUpEigVecNorm(VecsOld, EigVecs)

    use DVRData, only : para, grid

    real(dp), allocatable, intent(in) :: VecsOld(:,:,:)
    real(dp), allocatable, intent(inout) :: EigVecs(:,:,:)

    integer :: len_2, i, j, l, n_l, len_1, j_p
    integer, allocatable :: len_mid(:)

    n_l = para%l + 1
    allocate(len_mid(n_l))
!    allocate(len_1(n_l))

    len_1 = para%m1*para%nl
    len_2 = para%m2*para%nl - 1

    len_1 = orb%break
    len_2 = orb%n_max - orb%break - orb%n_red


    len_mid = orb%break + orb%n_red
    do l = 1, n_l
      len_mid(l) = len_mid(l) - l + 1
    end do


    write(iout, *) 'The integrals are now calculated for orbitals after leaving out some orbitals from middle'
    write(iout, '(a, i3)') ' Orbitals from the first half: 1 -', len_1
    write(iout, '(a,i3,a,i3)') ' Orbitals from the second half: ', len_mid(1)+1, ' -', len_mid(1)+len_2

    do i = 1, size(grid%r)
      write(78, '(11f10.6)') (VecsOld(i,j,3), j=1, size(VecsOld(1,:,2)))
    end do

    allocate(EigVecs(size(grid%r),orb%n_max-orb%n_red,n_l))
    EigVecs = 0.0d0
    do l = 1, n_l
      do j = 1, len_1 - l + 1
        do i = 1, size(grid%r)
          EigVecs(i,j,l) = VecsOld(i,j,l)
        end do
      end do
    end do

    do l = 1, n_l
      do j = 1, len_2 + l - 1
        j_p = j + len_1 - l + 1
        do i = 1, size(grid%r)
          EigVecs(i,j_p,l) = VecsOld(i,j+len_mid(l),l)
          !EigVecs(i+len_1,j+orb%n_inner,:) = VecsOld(i+len_1,j+len_mid,:)
        end do
      end do
    end do

    do i = 1, size(grid%r)
      write(79, '(11f10.6)') (EigVecs(i,j,3), j=1, size(EigVecs(1,:,2)))
    end do

  end subroutine SetUpEigVecNorm

  subroutine Calc1eOrbInts(EigVecs)

    use DVRData, only : one_e_rad_int, para, sph_harm

    real(dp), allocatable, intent(in) :: EigVecs(:,:,:)
    integer  :: i, j, l, l_val, m, n, error, ml, ind_1, ind_2
    real(dp) :: int_value, start, finish

    real(dp), allocatable :: inter_int(:,:)

    call cpu_time(start)

    allocate(inter_int(para%ng,orb%n_max), stat=error)
    call allocerror(error)

    allocate(OneEInts(orb%nSpatialOrbs,orb%nSpatialOrbs), stat=error)
    call allocerror(error)

    OneEInts = zero

    ! The transformation of 1e integrals is done in two steps
    ! i and j are basis indices, while m and n are orbital indices
    ! h_{ij} -> h_{mn} using transformation matrix b_{in}
    ! First step: h_{in} = \sum_j h_{ij} b_{in} 
    ! Second step: h_{mn} = \sum_i b^*_{im} h_{in} 

    ! The calculation is done over a loop of l

    do l = 1, para%l + 1
!   do l = 1, 1
      l_val = l - 1

      ! The first step is done here
      do n = 1, orb%n_max - (l-1)
        do i = 1, para%ng

          int_value = 0.0d0
          do j = 1, para%ng
             int_value = int_value + one_e_rad_int(i,j,l)*EigVecs(j,n,l)
!            write(77,'(3i4,3f20.10)') i, n, l, one_e_rad_int(i,j,l),EigVecs(j,n,l),int_value
          end do 
          inter_int(i,n) = int_value
!         write(iout,*) i, n, l, int_value
        end do 
      end do ! end loop over n_max

!     write(iout,*) 'The second step'
      ! The second step is done here, still it is inside the loop of l
      do m = 1, orb%n_max - (l-1)
        do n = 1, orb%n_max - (l-1)

          int_value = 0.0d0
          do i = 1, para%ng
            int_value = int_value + EigVecs(i,m,l) * inter_int(i,n)
!            write(78,'(3i4,3f20.10)') m, n, l, inter_int(i,n), EigVecs(i,m,l), int_value
          end do

!         write(iout,*) m, n, l, int_value

          do ml = 1, sph_harm%n_m(l)
!           ind_1 = SpatialOrbInd(m+l-1, l, ml)
            ind_1 = SpatialOrbInd(m, l, ml)
            if (ind_1.eq.0) write(iout,*) 'ind_1', m, l, ml, ind_1
!           ind_2 = SpatialOrbInd(n+l-1, l, ml)
            ind_2 = SpatialOrbInd(n, l, ml)
            if (ind_2.eq.0) write(iout,*) 'ind_2', n, l, ml, ind_2
            OneEInts(ind_1, ind_2) = int_value
!           write(iout,'(5i4,f20.10)') m, n, l, ind_1, ind_2, int_value
          end do

        end do ! end loop over n_max
      end do ! end loop over n_max

    end do ! end loop over para%l

    call cpu_time(finish)
    write(iout,'(X,a,f10.5,X,a)') 'Time taken for 1e transformation = ', finish-start, 'seconds.'

  end subroutine Calc1eOrbInts

  subroutine WriteInts(tol)
    
    use DVRData, only : para

    real(dp), intent(in) :: tol
    integer :: i, j, k, l, f_int, norbs, ij, kl, ijkl, i_n, j_n, k_n, l_n, nelec
    real(dp) :: h_core, int_value
    !integer :: i_p, j_p, k_p, l_p
    !integer, allocatable :: temp_pos(:)

    f_int = 15
    open(f_int, file=file_int, status='unknown', form="formatted")

    norbs = orb%nSpatialOrbs - nFrozen

    nelec = nint(para%z) - 2*nFrozen

    write(f_int,1001) norbs, nelec, 0

1001 format(' &FCI NORB=', i5, ',NELEC=', i3, ',MS2=', i3,',')

    write(f_int,1002, advance='no') 

1002 format('  ORBSYM=')
 
    do i = 1, norbs

      write(f_int, '(i1)', advance='no') 1

      if (i.lt.norbs) then
        write(f_int,1003, advance='no')
      else
        write(f_int,1004) 
      endif

    end do
    
1003 format(',')
1004 format(',')

    write(f_int,*) ' ISYM=1,'
    write(f_int, *) ' &END'

    !allocate(temp_pos(norbs))
    !temp_pos = (/1,2,3,4,5,6,7,8,9,11,12,13,15,16,17,18,20,21,22,26,27,28/)
    !temp_pos = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25,26,27,28,29/)
    ij  = 0
    ijkl = 0
    do i = 1, norbs
      i_n = i + nFrozen
      !i_p = temp_pos(i)
      do j = 1, i
        j_n = j + nFrozen
        !j_p = temp_pos(j)
        kl = 0
!       do k = 1, norbs
        do k = 1, i
          k_n = k + nFrozen
          !k_p = temp_pos(k)
          do l = 1, k
            l_n = l + nFrozen
            !l_p = temp_pos(l)
            if (ij.ge.kl) then
              if (abs(TwoEInts(i_n,k_n,j_n,l_n)).gt.tol) &
              & write(f_int, 1005) TwoEInts(i_n,k_n,j_n,l_n), i, j, k, l
              !& write(f_int, 1005) TwoEInts(i_n,k_n,j_n,l_n), i_p, j_p, k_p, l_p
              ijkl = ijkl + 1
            end if
            kl = kl + 1
          end do
        end do
        ij = ij + 1
      end do
    end do

!   ij  = 0
!   ijkl = 0
!   do i = 1, norbs
!     do j = 1, norbs
!       kl = 0
!       do k = 1, norbs
!         do l = 1, norbs
!             if (abs(TwoEInts(i,k,j,l)).gt.tol) &
!             & write(f_int, 1005) TwoEInts(i,k,j,l), i, j, k, l
!         end do
!       end do
!     end do
!   end do

    do i = 1, norbs
      i_n = i + nFrozen
      !i_p = temp_pos(i)
      do j = 1, i
        j_n = j + nFrozen
        !j_p = temp_pos(j)
        int_value = OneEInts(i_n,j_n)
!       if (nFrozen.gt.0.and.i.eq.j) then
        if (nFrozen.gt.0) then
          do k = 1, nfrozen
            int_value = int_value + 2.0d0*TwoEInts(i_n,k,j_n,k) - TwoEInts(i_n,k,k,j_n) 
          end do
        end if
        if (abs(int_value).gt.tol) &
!       if (i.eq.j) &
        & write(f_int, 1005) int_value, i, j, 0, 0
        !& write(f_int, 1005) int_value, i_p, j_p, 0, 0
      end do 
    end do

    h_core = 0.0_dp
    if (nFrozen.gt.0) then
      do i = 1, nFrozen
        h_core = h_core + 2.0d0*OneEInts(i,i)
        do j = 1, nFrozen
          h_core = h_core + 2.0d0*TwoEInts(i,j,i,j) - TwoEInts(i,j,j,i)
        end do
      end do  
    end if
    write(f_int, 1005) h_core, 0, 0, 0, 0

1005 format(f20.16,x,5i5)


    close(f_int)
  end subroutine WriteInts

end module Orbints
