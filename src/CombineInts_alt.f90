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

    do l1 = 1, n_l 
      la = l1 - 1
      n_m1 = 2*l1 - 1

      do l3 = 1, n_l
        lc = l3 - 1
        n_m3 = 2*l3 - 1

        do l2 = 1, n_l
          lb = l2 - 1
          n_m2 = 2*l2 - 1

           do l4 = 1, n_l
            ld = l4 - 1
            n_m4 = 2*l4 - 1


            do m1 = 1, n_m1
              lma = (l1-1)**2 + m1
              do m3 = 1, n_m3
                lmc = (l3-1)**2 + m3

                do m2 = 1, n_m2
                  lmb = (l2-1)**2 + m2
                  do m4 = 1, n_m4
                    lmd = (l4-1)**2 + m4

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

                            end do

                                TwoEInts(klm_1, klm_2, klm_3, klm_4) = int_value_dr

!                           write(78, '(4I5,X,f15.10)') klm_1, klm_2, klm_3, klm_4, TwoEInts(klm_1, klm_2, klm_3, klm_4)

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

    end do ! end loop over para%l

    call cpu_time(finish)
    write(iout,'(X,a,f10.5,X,a)') 'Time taken for 2e transformation = ', finish-start, 'seconds.'


  end subroutine Calc2eRadOrbInts_alt

end module CombineInts_alt
