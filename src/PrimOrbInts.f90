module PrimOrbInts

  use constants
  use util_mod, only : stop_all, allocerror
  use ReadInput, only : n_max, two_e_int, nfrz
  use OrbData

  implicit none

  contains

  subroutine GetPrimOrbInts()

    real(dp) :: tol

    tol = 1e-12

    call SetPrimOrbData()

    call Calc1ePrimOrbInts()

    call Calc2ePrimOrbInts()

    file_int = 'FCIDUMP'

    call WritePrimInts(tol)

  end subroutine GetPrimOrbInts

  subroutine SetPrimOrbData()

    use DVRData, only : para
    integer :: i, j, k, l_val, i_o, indx, error

    orb%n_max = para%ng

    allocate(SpatialOrbInd(n_max,para%l+1,2*para%l+1),stat=error)
    call allocerror(error)

    SpatialOrbInd = 0

    indx = 0
    l_val = para%l+1
!   write(*,*) 'limit for l_val:', l_val

    do i = 1, n_max
      do j = 1, l_val
        do k = 1, 2*j-1
          indx = indx + 1
          SpatialOrbInd(i,j,k) = indx
!         write(iout, *) i, j, k, SpatialOrbInd(i, j, k)
        end do
      end do
    end do
      
    orb%nSpatialOrbs = indx

    write(iout, *) '**********'
    write(iout, *) 'Setting up these parameters for the orbitals:'
    write(iout, '(X,A,3X, I6)') 'orb%n_max     =', n_max
    write(iout, '(X,A,3X, I6)') 'orb%nSpatialOrbs     =', indx
    write(iout, *) '***********' 

  end subroutine SetPrimOrbData

  subroutine Calc1ePrimOrbInts()

    use DVRData, only : one_e_rad_int, para, grid, eigen_vecs

    integer  :: l, m, n, error, ml, ind_1, ind_2
    real(dp) :: start, finish

    call cpu_time(start)

    allocate(OneEInts(orb%nSpatialOrbs,orb%nSpatialOrbs), stat=error)
    call allocerror(error)

    OneEInts = zero

    ! The 1-e integrals for the primitive orbitals are calculated here
    ! The 1-e integrals only has the radial part in it and has the same 
    ! value for all m quantum numbers

    do l = 1, para%l + 1

      do n = 1, orb%n_max
        do m = 1, orb%n_max
          do ml = 1, 2*l-1
            ind_1 = SpatialOrbInd(m, l, ml)
            ind_2 = SpatialOrbInd(n, l, ml)
            OneEInts(ind_1, ind_2) = one_e_rad_int(m,n,l)
!           write(77,'(4i4,f20.10)') m, n, l, ml, one_e_rad_int(m,n,l)
!           write(78,'(2i4,f20.10)') ind_1, ind_2, OneEInts(ind_1,ind_2)
          end do 
!         write(iout,*) i, n, l, int_value
        end do 
      end do ! end loop over n_max

    end do ! end loop over para%l

    call cpu_time(finish)
    write(iout,'(X,a,f10.5,X,a)') 'Time taken for 1e transformation = ', finish-start, 'seconds.'

  end subroutine Calc1ePrimOrbInts

  subroutine Calc2ePrimOrbInts()

    use DVRData, only : two_e_rad_int, para, grid, integrals_ang, sph_harm
    use OrbData, only : TwoEInts, orb, SpatialOrbInd

    integer  :: n_l, n_mp, l1, l2, l3, l4, m1, m2, m3, m4, k1, k2, l
    integer  :: klm_1, klm_2, klm_3, klm_4, error
    integer  :: n_m1, n_m2, n_m3, n_m4, lm1, lm2, lm3, lm4, n_orb
    real(dp) :: start, finish, int_value

    call cpu_time(start)

    n_l = para%l + 1
    n_mp = sph_harm%n_mp
    n_orb = orb%nSpatialOrbs

    allocate(TwoEInts(n_orb, n_orb, n_orb, n_orb), stat= error)
    call allocerror(error)

    TwoEints = zero

    do l1 = 1, n_l
      n_m1 = 2*l1 - 1
      do l2 = 1,  n_l
        n_m2 = 2*l2 - 1
        do l3 = 1,  n_l
          n_m3 = 2*l3 - 1
          do l4 = 1,  n_l
            n_m4 = 2*l4 - 1
            do m1 = 1, n_m1
              lm1 = (l1-1)**2 + m1
              do m2 = 1, n_m2
                lm2 = (l2-1)**2 + m2
                do m3 = 1, n_m3
                  lm3 = (l3-1)**2 + m3
                  do m4 = 1, n_m4
                    lm4 = (l4-1)**2 + m4
                    do k1 = 1, orb%n_max
                      do k2 = 1, orb%n_max
                        klm_1 = SpatialOrbInd(k1,l1,m1)
                        klm_2 = SpatialOrbInd(k2,l2,m2)
                        klm_3 = SpatialOrbInd(k1,l3,m3)
                        klm_4 = SpatialOrbInd(k2,l4,m4)
                      
                        int_value = 0.0d0
                        do l = 1, 2*para%l + 1
                          int_value = int_value + (integrals_ang(l, lm1, lm2, lm3, lm4)* &
                          &  two_e_rad_int(k1,k2,l))
                        end do

                        TwoEInts(klm_1, klm_2, klm_3, klm_4) = int_value

                      end do
                     end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine Calc2ePrimOrbInts

  subroutine WritePrimInts(tol)
    
    use DVRData, only : para

    real(dp), intent(in) :: tol
    integer :: i, j, k, l, f_int, norbs, ij, kl, ijkl, i_n, j_n, k_n, l_n, nelec
    real(dp) :: h_core, int_value

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

    ij  = 0
    ijkl = 0
    do i = 1, norbs
      i_n = i + nFrozen
      do j = 1, i
        j_n = j + nFrozen
        kl = 0
!       do k = 1, norbs
        do k = 1, i
          k_n = k + nFrozen
          do l = 1, k
            l_n = l + nFrozen
            if (ij.ge.kl) then
              if (abs(TwoEInts(i_n,k_n,j_n,l_n)).gt.tol) &
              & write(f_int, 1005) TwoEInts(i_n,k_n,j_n,l_n), i, j, k, l
              ijkl = ijkl + 1
            end if
            kl = kl + 1
          end do
        end do
        ij = ij + 1
      end do
    end do

    do i = 1, norbs
      i_n = i + nFrozen
      do j = 1, i
        j_n = j + nFrozen
        int_value = OneEInts(i_n,j_n)
        if (nFrozen.gt.0) then
          do k = 1, nfrozen
            int_value = int_value + 2.0d0*TwoEInts(i_n,k,j_n,k) - TwoEInts(i_n,k,k,j_n) 
          end do
        end if
        if (abs(int_value).gt.tol) &
        & write(f_int, 1005) int_value, i, j, 0, 0
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
  end subroutine WritePrimInts

end module PrimOrbInts

