module CombineInts

  use constants
  use util_mod
  use DVRData, only : combined_two_e_int, para, grid

  implicit none

  contains

  subroutine CombineIntegrals() 

  real(idp), allocatable     :: rad_tp_elements(:,:)
  real(idp), allocatable     :: ang_tp_elements(:,:)
  integer                    :: i, j, a, b, c, d, l, l_ang_max, k, k_max
  integer                    :: a_rad, a_ang, b_rad, b_ang 
  integer                    :: c_rad, c_ang, d_rad, d_ang 
  integer                    :: running_index_tot
  integer                    :: running_index_rad, running_index_ang
  integer                    :: n_lm, n_n, n_all

  k_max = 2*para%l    ! Multipole order
  
  call read_ascii_table(rad_tp_elements,                                       &
  &                     'Results/twoparticle_rad_elements_l0.dat')

  call read_ascii_table(ang_tp_elements,                                       &
  &                     'Results/ang_element_l0.dat')

  n_n = nint(sqrt(real(size(rad_tp_elements(:,1)),kind=idp)))
  n_lm = nint(sqrt(sqrt(real(size(ang_tp_elements(:,1)),kind=idp))))
  n_all = n_n * n_lm 

  write(*,*) n_n, n_lm, n_all
  allocate(combined_two_e_int(n_all**4))
  combined_two_e_int = zero
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!! Four-Index Matrix Elements !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! <a b | c d> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  do k = 0, k_max
    
    write(iout,*) 'Combining two electron radial and angular integrals for l:', k
    if (allocated(rad_tp_elements)) deallocate(rad_tp_elements)
    if (allocated(ang_tp_elements)) deallocate(ang_tp_elements)

    write(*,*) 'Get here 1'
    call read_ascii_table(rad_tp_elements,                                     &
    &            'Results/twoparticle_rad_elements_l'//trim(int2str(k))//'.dat')
    call read_ascii_table(ang_tp_elements,                                     &
    &                         'Results/ang_element_l'//trim(int2str(k))//'.dat')
    write(*,*) 'Get here 2'
    running_index_tot = 0
    running_index_rad = 0
    do a_rad = 1, n_n
    do b_rad = 1, n_n
      ! Already at this stage due to Kronecker Delta in radial part
      ! \delta_ac \delta_bd 
      running_index_rad = running_index_rad + 1
      do c_rad = 1, n_n
      do d_rad = 1, n_n
      running_index_ang = 0
        do a_ang = 1, n_lm
        do b_ang = 1, n_lm
          do c_ang = 1, n_lm
          do d_ang = 1, n_lm
            if (a_rad .ne. c_rad) cycle ! Radial Kronecker Delta \delta_ac
            if (b_rad .ne. d_rad) cycle ! Radial Kronecker Delta \delta_bd
            running_index_ang = running_index_ang + 1
            running_index_tot = running_index_tot + 1
            write(76,*) running_index_rad, running_index_ang, running_index_tot
            flush(6)
            combined_two_e_int(running_index_tot) =                          &
            & combined_two_e_int(running_index_tot) +                        &
            & rad_tp_elements(running_index_rad,3) *                           &
            & ang_tp_elements(running_index_ang,5)
          end do
          end do
        end do
        end do
      end do
      end do
    end do
    end do
    
    write(*,*) 'Get here 3'
  end do

  open(11, file="twoparticle_combined_elements.dat",                           &
  &    form="formatted", action="write")
  write(11,*) '# Primitive Combined Matrix Elements for the Four-index '//     &
  &           'Integrals for Multipole Order up to k = '//trim(int2str(k_max))
  running_index_tot = 0
  do a = 1, n_all
    do b = 1, n_all
      do c = 1, n_all
        do d = 1, n_all
          running_index_tot = running_index_tot + 1
          write(11, '(4I8,ES25.17)') a, b, c, d,                               &
          &                              combined_two_e_int(running_index_tot)
        end do
      end do
    end do
  end do
  close(11)

  end subroutine CombineIntegrals

  subroutine CombineOrbInts_old()

    use DVRData, only : integrals_ang, sph_harm
    use OrbData, only : TwoERadOrbInts, TwoEInts, orb, SpatialOrbInd

    integer  :: n_l, n_mp, l1, l2, l3, l4, m1, m2, m3, m4, k1, k2, k3, k4, l, indx
    integer  :: la, lb, lc, ld, ma, mb, mc, md, n_m1, n_m2, n_m3, n_m4, error
    integer  :: m1_init, m2_init, m3_init, m4_init, lma, lmb, lmc, lmd
    integer  :: klm_1, klm_2, klm_3, klm_4
    real(dp) :: start, finish, int_value_dr, int_value_xc


    call cpu_time(start)

    allocate(TwoEInts(orb%nSpatialOrbs,orb%nSpatialOrbs, orb%nSpatialOrbs, &
    & orb%nSpatialOrbs), stat=error)
    call allocerror(error)

    TwoEInts = zero

    n_l  = sph_harm%n_l
    n_mp = sph_harm%n_mp

    indx = 0
    do l1 = 1, n_l
      la = l1 - 1
      n_m1 = 2*l1 - 1
      m1_init = -1*l1

      do l3 = 1, n_l
        lc = l3 - 1
        n_m3 = 2*l3 - 1
        m3_init = -1*l3

        do l2 = 1, n_l
          lb = l2 - 1
          n_m2 = 2*l2 - 1
          m2_init = -1*l2

          do l4 = 1, n_l
            ld = l4 - 1
            n_m4 = 2*l4 - 1
            m4_init = -1*l4


            do m1 = 1, n_m1
              ma = m1_init + m1
              lma = (l1-1)**2 + m1
              !lma = (l1-1)**2 + int(la) + int(ma) + 1
              do m3 = 1, n_m3
                mc = m1_init + m3
                lmc = (l3-1)**2 + m3
                !lmc = (l3-1)**2 + int(lc) + int(mc) + 1

                do m2 = 1, n_m2
                  mb = m2_init + m2
                  lmb = (l2-1)**2 + m2
                  !lmb = (l2-1)**2 + int(lb) + int(mb) + 1
                  do m4 = 1, n_m4
                    md = m2_init + m4
                    lmd = (l4-1)**2 + m4
                    !lmd = (l4-1)**2 + int(ld) + int(md) + 1

                    do k1 = 1, orb%n_max - l1 + 1
                      do k2 = 1, orb%n_max - l2 + 1  
                        do k3 = 1, orb%n_max - l3 + 1
                          do k4 = 1, orb%n_max - l4 + 1  
                            indx = indx + 1
                            int_value_dr = 0.0d0
                            int_value_xc = 0.0d0

                            klm_1 = SpatialOrbInd(k1,l1,m1)
                            klm_2 = SpatialOrbInd(k2,l2,m2)
                            klm_3 = SpatialOrbInd(k3,l3,m3)
                            klm_4 = SpatialOrbInd(k4,l4,m4)

                            if (klm_1.eq.0.or.klm_2.eq.0.or.klm_3.eq.0.or.klm_4.eq.0) cycle

                            do l = 1, 2*para%l + 1

                              int_value_dr = int_value_dr + (integrals_ang(l, lma, lmb, lmc, lmd)* &
                              &  TwoERadOrbInts(k1,k2,k3,k4,l1,l2,l3,l4,l))

!                             if (l1.eq.l3.and.l2.eq.l4) then
!                               write(85,'(7i3,3f13.8)') k1,k2,k3,k4,l1,l2,l,integrals_ang(l, lma, lmb, lmc, lmd), TwoERadOrbInts(k1,k2,k3,k4,l1,l2,l1,l2,l), int_value_dr
!                             end if
                              if (lma.eq.lmd.and.lmb.eq.lmc.and.l1.ne.l2.and.m1.ne.m2) then
!                             if (lma.eq.lmd.and.lmb.eq.lmc.and.l1.eq.1.and.l2.eq.2) then
!                               write(87,'(11i3,3f13.8)') k1,k2,k4,k3,l1,l2,l, m1, m2, m4, m3, integrals_ang(l, lma, lmb, lmc, lmd), TwoERadOrbInts(k1,k2,k3,k4,l1,l2,l3,l4,l), int_value_dr
!                               write(89,'(5i3,3f13.8)') klm_1, klm_2, klm_4, klm_3,l,integrals_ang(l,lma,lmb,lmc,lmd), TwoERadOrbInts(k1,k2,k3,k4,l1,l2,l3,l4,l), int_value_dr
                              end if

!                             if (l1.ne.l2.or.m1.ne.m2) &
!                             &  int_value_xc = int_value_xc + (integrals_ang(l, lma, lmb, lmd, lmc)* &
!                             &  TwoERadOrbInts(k1,k2,k4,k3,l1,l2,l2,l1,l))

 !                            if (klm_1.eq.1.and.klm_2.eq.1.and.klm_3.eq.1) then
 !                              write(77,'(7I4,X,2F15.10)') k1, k2, lma, lmb, lmc, lmd, klm_4, integrals_ang(l, lma, lmb, lmc, lmd), TwoERadOrbInts(k1,k2,k3,k4,l1,l2,l)
 !                            end if
                            end do

 !                            write(83,'(6i5,f15.10)') k1,k2,k3,k4,l1,l2, int_value_dr
 !                            write(84,'(6i5,f15.10)') k1,k2,k3,k4,l1,l2, int_value_xc

!                           if (lma.eq.lmd.and.lmb.eq.lmc.and.l1.ne.l2.and.m1.ne.m2) then
!                           if (lma.eq.lmd.and.lmb.eq.lmc) then
!                           TwoEInts(klm_1, klm_2, klm_4, klm_3) = int_value_dr
!                           else 
                            TwoEInts(klm_1, klm_2, klm_3, klm_4) = int_value_dr
!                           end if

 !                          write(78, '(4I5,X,f15.10)') klm_1, klm_2, klm_3, klm_4, TwoEInts(klm_1, klm_2, klm_3, klm_4)

!                           if (l1.ne.l2.or.m1.ne.m2) then 
!                             TwoEInts(klm_1, klm_2, klm_4, klm_3) = int_value_xc
!                             TwoEInts(klm_4, klm_2, klm_1, klm_3) = int_value_xc
!                             TwoEInts(klm_1, klm_3, klm_4, klm_2) = int_value_xc
!                             TwoEInts(klm_4, klm_3, klm_1, klm_2) = int_value_xc
!!                            write(79, '(4I5,X,f15.10)') klm_1, klm_2, klm_4, klm_3, TwoEInts(klm_1, klm_2, klm_4, klm_3)
!                           end if

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
      end do  
    end do  

!   write(iout, *) 'Index:', indx

!   do l1 = 1, orb%nSpatialOrbs
!     do l2 = 1, orb%nSpatialOrbs
!       do l3 = 1, orb%nSpatialOrbs
!         do l4 = 1, orb%nSpatialOrbs
!!          write(80,'(4i5, f20.12)')  l1, l2, l3, l4, TwoEInts(l1,l2,l3,l4)
!           if (abs(TwoEInts(l1,l2,l3,l4)).gt.0.0d0) then
!!            write(81,*)  l1, l2, l3, l4, TwoEInts(l1,l2,l3,l4)
!           end if
!         end do
!       end do
!     end do
!   end do

    call cpu_time(finish)

    write(iout,'(X,a,f10.5,X,a)') 'Time taken for combining integrals = ', finish-start, 'seconds.'

  end subroutine CombineOrbInts_old

  subroutine CombineOrbInts()

    use DVRData, only : integrals_ang, sph_harm
    use OrbData, only : TwoERadOrbInts, TwoERadOrbInts_dr, TwoERadOrbInts_xc, TwoEInts, orb, SpatialOrbInd

    integer  :: n_l, n_mp, l1, l2, l3, l4, m1, m2, m3, m4, k1, k2, k3, k4, l, indx
    integer  :: la, lb, lc, ld, ma, mb, mc, md, n_m1, n_m2, n_m3, n_m4, error
    integer  :: m1_init, m2_init, m3_init, m4_init, lma, lmb, lmc, lmd
    integer  :: klm_1, klm_2, klm_3, klm_4
    real(dp) :: start, finish, int_value_dr, int_value_xc


    call cpu_time(start)

    allocate(TwoEInts(orb%nSpatialOrbs,orb%nSpatialOrbs, orb%nSpatialOrbs, &
    & orb%nSpatialOrbs), stat=error)
    call allocerror(error)

    TwoEInts = zero

    n_l  = sph_harm%n_l
    n_mp = sph_harm%n_mp

    indx = 0
    do l1 = 1, n_l
      la = l1 - 1
      n_m1 = 2*l1 - 1
      m1_init = -1*l1

!     do l3 = 1, n_l
!       lc = l3 - 1
!       n_m3 = 2*l3 - 1
!       m3_init = -1*l3

        do l2 = 1, n_l
          lb = l2 - 1
          n_m2 = 2*l2 - 1
          m2_init = -1*l2

!         do l4 = 1, n_l
!           ld = l4 - 1
!           n_m4 = 2*l4 - 1
!           m4_init = -1*l4


            do m1 = 1, n_m1
              ma = m1_init + m1
              lma = (l1-1)**2 + m1
              !lma = (l1-1)**2 + int(la) + int(ma) + 1
              do m3 = 1, n_m1
                mc = m1_init + m3
                lmc = (l1-1)**2 + m3
                !lmc = (l3-1)**2 + int(lc) + int(mc) + 1

                do m2 = 1, n_m2
                  mb = m2_init + m2
                  lmb = (l2-1)**2 + m2
                  !lmb = (l2-1)**2 + int(lb) + int(mb) + 1
                  do m4 = 1, n_m2
                    md = m2_init + m4
                    lmd = (l2-1)**2 + m4
                    !lmd = (l4-1)**2 + int(ld) + int(md) + 1

                    do k1 = 1, orb%n_max - l1 + 1
                      do k2 = 1, orb%n_max -l2 + 1  
                        do k3 = 1, orb%n_max - l1 + 1
                          do k4 = 1, orb%n_max -l2 + 1  

                            indx = indx + 1
                            int_value_dr = 0.0d0
                            int_value_xc = 0.0d0

                            klm_1 = SpatialOrbInd(k1,l1,m1)
                            klm_2 = SpatialOrbInd(k2,l2,m2)
                            klm_3 = SpatialOrbInd(k3,l1,m3)
                            klm_4 = SpatialOrbInd(k4,l2,m4)

                            do l = 1, 2*para%l + 1

                              int_value_dr = int_value_dr + (integrals_ang(l, lma, lmb, lmc, lmd)* &
!                             &  TwoERadOrbInts(k1,k2,k3,k4,l1,l2,l1,l2,l))
                              &  TwoERadOrbInts_dr(k1,k2,k3,k4,l1,l2,l))

!                             write(86,'(7i3,3f13.8)') k1,k2,k3,k4,l1,l2,l,integrals_ang(l, lma, lmb, lmc, lmd), TwoERadOrbInts_dr(k1,k2,k3,k4,l1,l2,l)!, int_value_dr

                              if (l1.ne.l2.or.m1.ne.m2) &
                              & int_value_xc = int_value_xc + (integrals_ang(l, lma, lmb, lmd, lmc)* &
!                             &  TwoERadOrbInts(k1,k2,k4,k3,l1,l2,l2,l1,l))
                              &  TwoERadOrbInts_xc(k1,k2,k3,k4,l1,l2,l))

!                             if (l1.ne.l2.or.m1.ne.m2) &
!                             & write(88,'(11i3,3f13.8)') k1,k2,k3,k4,l1,l2,l, m1, m2, m3, m4, integrals_ang(l, lma, lmb, lmd, lmc), TwoERadOrbInts_xc(k1,k2,k3,k4,l1,l2,l), int_value_xc

!                             if (l1.ne.l2.or.m1.ne.m2) &
!                             & write(90,'(5i3,3f13.8)') klm_1, klm_2, klm_3, klm_4, l,integrals_ang(l, lma, lmb, lmd, lmc), TwoERadOrbInts_xc(k1,k2,k3,k4,l1,l2,l), int_value_xc

!                             if (klm_1.eq.1.and.klm_2.eq.1.and.klm_3.eq.1) then
!                               write(77,'(7I4,X,2F15.10)') k1, k2, lma, lmb, lmc, lmd, klm_4, integrals_ang(l, lma, lmb, lmc, lmd), TwoERadOrbInts(k1,k2,k3,k4,l1,l2,l)
!                             end if

                            end do
!                             write(85,'(6i5,f15.10)') k1,k2,k3,k4,l1,l2, int_value_dr
!                             write(86,'(6i5,f15.10)') k1,k2,k3,k4,l1,l2, int_value_xc
!                             if (klm_1.eq.1.and.klm_2.eq.1.and.klm_3.eq.1) then
!                               write(77,*) ''
!                             end if
                            TwoEInts(klm_1, klm_2, klm_3, klm_4) = int_value_dr
!                           write(77, '(6I5,X,f15.10)') k1, k2, lma, lmb, lmc, lmd, int_value
!                           write(78, '(4I5,X,f15.10)') klm_1, klm_2, klm_3, klm_4, TwoEInts(klm_1, klm_2, klm_3, klm_4)
                            if (l1.ne.l2.or.m1.ne.m2) then 
                              TwoEInts(klm_1, klm_2, klm_4, klm_3) = int_value_xc
                              TwoEInts(klm_4, klm_2, klm_1, klm_3) = int_value_xc
                              TwoEInts(klm_1, klm_3, klm_4, klm_2) = int_value_xc
                              TwoEInts(klm_4, klm_3, klm_1, klm_2) = int_value_xc
!                             write(79, '(4I5,X,f15.10)') klm_1, klm_2, klm_4, klm_3, TwoEInts(klm_1, klm_2, klm_4, klm_3)
                            end if
!                           write(79, '(4I5,X,f15.10)') lma, lmb, lmc, lmd, int_value
                          end do
                        end do
                      end do
                    end do
                  
                  end do
                end do  
              end do
            end do

!         end do
        end do
!     end do  
    end do  

!   write(iout, *) 'Index:', indx

!   do l1 = 1, orb%nSpatialOrbs
!     do l2 = 1, orb%nSpatialOrbs
!       do l3 = 1, orb%nSpatialOrbs
!         do l4 = 1, orb%nSpatialOrbs
!!          write(80,'(4i5, f20.12)')  l1, l2, l3, l4, TwoEInts(l1,l2,l3,l4)
!           if (abs(TwoEInts(l1,l2,l3,l4)).gt.0.0d0) then
!!            write(81,*)  l1, l2, l3, l4, TwoEInts(l1,l2,l3,l4)
!           end if
!         end do
!       end do
!     end do
!   end do

    call cpu_time(finish)

    write(iout,'(X,a,f10.5,X,a)') 'Time taken for combining integrals = ', finish-start, 'seconds.'
  end subroutine CombineOrbInts

end module CombineInts
