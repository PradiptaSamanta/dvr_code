module CombineInts_alt

  use constants
  use util_mod

  implicit none

  contains

  subroutine CombineOrbInts_alt()

    use DVRData, only : integrals_ang, sph_harm, para
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

!   do k1 = 1, orb%n_max
!     do k2 = 1, k1 + 
!       do k3 = 1, k1
!         do k4 = 1, k2

    do l1 = 1, n_l 
!   do l1 = 1, min(n_l,(orb%n_max-k1+1))
      la = l1 - 1
      n_m1 = 2*l1 - 1
      m1_init = -1*l1

      do l3 = 1, n_l
!     do l3 = 1, min(n_l,(orb%n_max-k3+1))
        lc = l3 - 1
        n_m3 = 2*l3 - 1
        m3_init = -1*l3
!       if (l3.ne.l1) cycle

        do l2 = 1, n_l
!       do l2 = 1, min(n_l,(orb%n_max-k2+1))
          lb = l2 - 1
          n_m2 = 2*l2 - 1
          m2_init = -1*l2

           do l4 = 1, n_l
!         do l4 = 1, min(n_l,(orb%n_max-k4+1))
            ld = l4 - 1
            n_m4 = 2*l4 - 1
            m4_init = -1*l4
!           if (l4.ne.l2) cycle


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
!                       do k3 = 1, orb%n_max - l3 + 1
!                         do k4 = 1, orb%n_max - l4 + 1  
                        do k3 = 1, k1
                          do k4 = 1, k2
                            indx = indx + 1
                            int_value_dr = 0.0d0
                            int_value_xc = 0.0d0

                            klm_1 = SpatialOrbInd(k1,l1,m1)
                            klm_2 = SpatialOrbInd(k2,l2,m2)
                            klm_3 = SpatialOrbInd(k3,l3,m3)
                            klm_4 = SpatialOrbInd(k4,l4,m4)

                            write(85,'(4i4)') klm_1, klm_2, klm_3, klm_4

                            if (klm_1.eq.0.or.klm_2.eq.0.or.klm_3.eq.0.or.klm_4.eq.0) cycle

                            do l = 1, 2*para%l + 1

                              int_value_dr = int_value_dr + (integrals_ang(l, lma, lmb, lmc, lmd)* &
                              &  TwoERadOrbInts(k1,k2,k3,k4,l1,l2,l3,l4,l))

!                             if (l1.eq.l3.and.l2.eq.l4) then
!                               write(85,'(7i3,3f13.8)') k1,k2,k3,k4,l1,l2,l,integrals_ang(l, lma, lmb, lmc, lmd), TwoERadOrbInts(k1,k2,k3,k4,l1,l2,l1,l2,l), int_value_dr
!                             end if
                              if (lma.eq.lmd.and.lmb.eq.lmc.and.l1.ne.l2.and.m1.ne.m2) then
!                             if (lma.eq.lmd.and.lmb.eq.lmc.and.l1.eq.1.and.l2.eq.2) then
                                write(87,'(11i3,3f13.8)') k1,k2,k4,k3,l1,l2,l, m1, m2, m4, m3, integrals_ang(l, lma, lmb, lmc, lmd), TwoERadOrbInts(k1,k2,k3,k4,l1,l2,l3,l4,l), int_value_dr
!                               write(89,'(5i3,3f13.8)') klm_1, klm_2, klm_4, klm_3,l,integrals_ang(l,lma,lmb,lmc,lmd), TwoERadOrbInts(k1,k2,k3,k4,l1,l2,l3,l4,l), int_value_dr
                              end if
!                             if (l1.ne.l2.or.m1.ne.m2) then
                                if (abs(integrals_ang(l, lma, lmb, lmc, lmd)).gt.1e-12) &
                                & write(89,'(5i3,3f13.8)') klm_1, klm_2, klm_4, klm_3, l,integrals_ang(l, lma, lmb, lmc, lmd), TwoERadOrbInts(k1,k2,k3,k4,l1,l2,l3,l4,l), int_value_dr
!                             end if

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

!                           if (abs(integrals_ang(1, lma, lmb, lmc, lmd)).gt.1e-12) then
                                if (abs(TwoEInts(klm_1, klm_2, klm_3, klm_4)).gt.1e-12) then
                                write(*,'(a,4i4,2f13.8)') '1 overwriting it for:,', klm_1, klm_2, klm_3, klm_4, TwoEInts(klm_1, klm_2, klm_3, klm_4), int_value_dr
                                end if
                                TwoEInts(klm_1, klm_2, klm_3, klm_4) = int_value_dr
                                if (abs(TwoEInts(klm_1, klm_2, klm_3, klm_4)).gt.1e-12) write(91,'(4i3,f13.8)') klm_1, klm_2, klm_3, klm_4,  TwoEInts(klm_1, klm_2, klm_3, klm_4)
!                           elseif(abs(integrals_ang(2, lma, lmb, lmc, lmd)).gt.1e-12) then
!                               if (abs(TwoEInts(klm_1, klm_2, klm_4, klm_3)).gt.1e-12) then
!                               write(*,'(a,4i4,2f13.8)') '2 overwriting it for:,', klm_1, klm_2, klm_3, klm_4, TwoEInts(klm_1, klm_2, klm_4, klm_3), int_value_dr
!                               end if
!                               TwoEInts(klm_1, klm_2, klm_4, klm_3) = int_value_dr
!                           end if

!                           end if

                            write(78, '(4I5,X,f15.10)') klm_1, klm_2, klm_3, klm_4, TwoEInts(klm_1, klm_2, klm_3, klm_4)

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

    do l1 = 1, n_l
      n_m1 = 2*l1 - 1
      do l2 = 1, n_l
        n_m2 = 2*l2 - 1
        do l3 = 1, n_l
          n_m3 = 2*l3 - 1
!         if (l1.ne.l3) cycle
          do l4 = 1, n_l
            n_m4 = 2*l4 - 1
!           if (l2.ne.l4) cycle

        do m1 = 1, n_m1
          lma = (l1-1)**2 + m1
          do m3 = 1, n_m3
            lmc = (l3-1)**2 + m3

            do m2 = 1, n_m2
              lmb = (l2-1)**2 + m2
              do m4 = 1, n_m4
                lmd = (l4-1)**2 + m4

                do k1 = 1, orb%n_max - l1 + 1
                  do k2 = 1, orb%n_max -l2 + 1  
!                  do k3 = 1, orb%n_max -l3 + 1  
!                     do k4 = 1, orb%n_max - l4 + 1
                   do k3 = 1, k1
                      do k4 = 1, k1

                        klm_1 = SpatialOrbInd(k1,l1,m1)
                        klm_2 = SpatialOrbInd(k2,l2,m2)
                        klm_3 = SpatialOrbInd(k3,l3,m3)
                        klm_4 = SpatialOrbInd(k4,l4,m4)

                        int_value_dr = 0.0d0
                        do l = 1, 2*para%l + 1

!                         if (l1.ne.l2.or.m1.ne.m2) then
!!                          int_value_dr = int_value_dr + (integrals_ang(l, lma, lmb, lmd, lmc)* &
!!                        & TwoERadOrbInts(k1,k2,k4,k3,l1,l2,l4,l3,l))
!!                          if (abs(integrals_ang(l, lma, lmb, lmd, lmc)).gt.1e-12) &
!!                          & write(89,'(5i3,3f13.8)') klm_1, klm_2, klm_3, klm_4, l,integrals_ang(l, lma, lmb, lmd, lmc), TwoERadOrbInts(k1,k2,k4,k3,l1,l2,l4,l3,l), int_value_dr
!                           int_value_dr = int_value_dr + (integrals_ang(l, lma, lmb, lmc, lmd)* &
!                         & TwoERadOrbInts(k1,k2,k3,k4,l1,l2,l3,l4,l))
!                           if (abs(integrals_ang(l, lma, lmb, lmc, lmd)).gt.1e-12) &
!                           & write(89,'(5i3,3f13.8)') klm_1, klm_2, klm_4, klm_3, l,integrals_ang(l, lma, lmb, lmc, lmd), TwoERadOrbInts(k1,k2,k3,k4,l1,l2,l3,l4,l), int_value_dr
!                         end if

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

  end subroutine CombineOrbInts_alt

  subroutine Calc2eRadOrbInts_alt()

    use DVRData, only : two_e_rad_int, para, grid, eigen_vecs
    use OrbData, only : orb, TwoERadOrbInts

    integer  :: i, j, l, l1, l2, l3 ,l4, n_l, l_val, m, n, mp, np, error, ml, ind_1, ind_2
    real(dp) :: int_value, start, finish, int_value_xc, int_value_dr
    integer  :: count_i, count_j, count_l
    logical  :: split

    real(dp), allocatable :: inter_int(:,:,:,:,:)

    call cpu_time(start)

    n_l = para%l+1

    allocate(TwoERadOrbInts(orb%n_max,orb%n_max, orb%n_max, orb%n_max, n_l, n_l, n_l, n_l, 2*para%l+1), stat=error)
    call allocerror(error)
    allocate(inter_int(para%ng,orb%n_max,orb%n_max,n_l,n_l), stat=error)
    call allocerror(error)


    inter_int = zero
    TwoERadOrbInts = zero

    ! The transformation of 1e integrals is done in two steps
    ! i and j are basis indices, while m and n are orbital indices
    ! h_{ij} -> h_{mn} using transformation matrix b_{in}
    ! First step: h_{in} = \sum_j h_{ij} b_{in} 
    ! Second step: h_{mn} = \sum_i b^*_{im} h_{in} 

    ! The calculation is done over a loop of l

!   write(*,*) 'Size 1', size(eigen_vecs(:,1,1))
!   write(*,*) 'Size 2', size(eigen_vecs(1,:,1))
!   write(*,*) 'Size 3', size(eigen_vecs(1,1,:))

!   flush(6)

!   split = .false.
    split = .true.

    do l = 1, 2*para%l + 1
      l_val = l - 1

      if (split) then

      ! The first step is done here
      do m = 1, orb%n_max
        do n = 1, orb%n_max
          do l1 = 1, n_l
            do l2 = 1, n_l
              do i = 1, para%ng
       
                int_value = 0.0d0
                do j = 1, para%ng
                   int_value = int_value + two_e_rad_int(i,j,l)*eigen_vecs(j,n,l2)*eigen_vecs(j,m,l1)
                end do 
                inter_int(i,m,n,l1,l2) = int_value
!               if (l1.eq.l2.and.abs(int_value).gt.1e-12) write(87,'(4i3,f13.8)') i,m,n,l1,int_value
 !              write(76,'(3I5,X,F20.10)') i, n, l, eigen_vecs(i,n,l)
              end do
            end do 
          end do 
        end do
      end do ! end loop over n_max

      ! The second step is done here, still it is inside the loop of l

      do m = 1, orb%n_max
        do n = 1, orb%n_max
          do mp = 1, orb%n_max
            do np = 1, orb%n_max
              do l1 = 1, n_l
                  do l3 = 1, n_l
                do l2 = 1, n_l
                    do l4 = 1, n_l

                      int_value = 0.0d0
                      int_value_xc = 0.0d0
                      do i = 1, para%ng
                        int_value = int_value + eigen_vecs(i,mp,l1)*eigen_vecs(i,np,l3)*inter_int(i,m,n,l2,l4)
                        int_value_xc = int_value_xc + eigen_vecs(i,mp,l1)*eigen_vecs(i,n,l4)*inter_int(i,m,np,l2,l3)
                      end do
          
                      !if(l1.eq.l3.and.l2.eq.l4) write(82,'(7i4,f16.8)') m, n , mp, np, l1, l2, l, int_value_xc
                      !if (l1.eq.l3.and.l2.eq.l4) write(82,'(7i4,f16.8)') m, n , mp, np, l1, l2, l, int_value
          
                      TwoERadOrbInts(mp, m, np, n, l1, l2, l3, l4, l) = int_value
                    end do
                  end do
                end do
              end do ! end loop over n_max
            end do ! end loop over n_max
          end do
        end do ! end loop over n_max
      end do ! end loop over n_max

  !   else 

  !     do m = 1, orb%n_max
  !       do n = 1, orb%n_max
  !         do mp = 1, orb%n_max
  !           do np = 1, orb%n_max
  !             int_value = 0.0d0
  !             do i = 1, para%ng
  !               do j = 1, para%ng
  !                 int_value = int_value + (eigen_vecs(i,mp,l)*eigen_vecs(i,np,l))*two_e_rad_int(i,j,l)*(eigen_vecs(j,m,l)*eigen_vecs(j,n,l))
! !        write(80,'(4i4, 4f16.8)') m, n, i, j, eigen_vecs(i,m,l), eigen_vecs(j,n,l),two_e_rad_int(i,j,l), int_value
  !               end do 
  !             end do 
  !             TwoERadOrbInts(mp, m, np, n, l) = int_value
  !           write(78,'(5i4,f16.8)') mp, np, m, n, l, int_value
! !     write(82,'(2I4,ES25.17)') m, n, int_value
  !           end do ! end loop over n_max
  !         end do ! end loop over n_max
  !       end do ! end loop over n_max
  !     end do ! end loop over n_max

      end if

 !    do m = 1, orb%n_max
 !      do n = 1, orb%n_max
 !        do mp = 1, orb%n_max
 !          do np = 1, orb%n_max
 !            do l1 = 1, n_l
 !              do l2 = 1, n_l

 !                    write(82,'(7i4,f16.8)') m, n , mp, np, l1, l2, l, TwoERadOrbInts(mp, m, n, np, l1, l2, l2, l1, l)

 !              end do
 !            end do ! end loop over n_max
 !          end do ! end loop over n_max
 !        end do
 !      end do ! end loop over n_max
 !    end do ! end loop over n_max
    end do ! end loop over para%l

    call cpu_time(finish)
    write(iout,'(X,a,f10.5,X,a)') 'Time taken for 2e transformation = ', finish-start, 'seconds.'


  end subroutine Calc2eRadOrbInts_alt

end module CombineInts_alt
