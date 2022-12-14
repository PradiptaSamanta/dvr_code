module CombineInts

  use constants
  use util_mod

  implicit none

  contains

  subroutine CombineOrbInts()

    use DVRData, only : integrals_ang, sph_harm, para
    use OrbData, only : TwoERadOrbInts, TwoEInts, orb, SpatialOrbInd

    integer  :: n_l, n_mp, l1, l2, l3, l4, m1, m2, m3, m4, k1, k2, k3, k4, l, indx
    integer  :: n_m1, n_m2, n_m3, n_m4, error
    integer  :: lma, lmb, lmc, lmd
    integer  :: klm_1, klm_2, klm_3, klm_4
    integer, allocatable :: l_interm(:)
    real(dp) :: start, finish, int_value_dr, int_value_xc


    call cpu_time(start)

    allocate(TwoEInts(orb%nSpatialOrbs,orb%nSpatialOrbs, orb%nSpatialOrbs, &
    & orb%nSpatialOrbs), stat=error)
    call allocerror(error)

    TwoEInts = zero

    allocate(l_interm(sph_harm%n_l), stat= error)
    call allocerror(error)

    n_l  = sph_harm%n_l
    n_mp = sph_harm%n_mp

    l_interm(1) = 0
    do l = 2, sph_harm%n_l
      l_interm(l) = sum(sph_harm%n_m(1:l-1))
    end do

    indx = 0

    do l1 = 1, n_l 
      n_m1 = sph_harm%n_m(l1)
      !n_m1 = 2*l1 - 1

      do l3 = 1, n_l
        n_m3 = sph_harm%n_m(l3)
        !n_m3 = 2*l3 - 1

        do l2 = 1, n_l
          n_m2 = sph_harm%n_m(l2)
          !n_m2 = 2*l2 - 1

          do l4 = 1, n_l
            n_m4 = sph_harm%n_m(l4)
            !n_m4 = 2*l4 - 1


            do m1 = 1, n_m1
              !lma = (l1-1)**2 + m1
              lma = l_interm(l1) + m1
              do m3 = 1, n_m3
                !lmc = (l3-1)**2 + m3
                lmc = l_interm(l3) + m3

                do m2 = 1, n_m2
                  !lmb = (l2-1)**2 + m2
                  lmb = l_interm(l2) + m2
                  do m4 = 1, n_m4
                    !lmd = (l4-1)**2 + m4
                    lmd = l_interm(l4) + m4

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

  end subroutine CombineOrbInts

  subroutine Calc2eRadOrbInts(EigVecs)

    use DVRData, only : two_e_rad_int, para
    use OrbData, only : orb, TwoERadOrbInts

    real(dp), allocatable, intent(in) :: EigVecs(:,:,:)
    integer  :: i, j, l, l1, l2, l3 ,l4, n_l, l_val, m, n, mp, np, error
    real(dp) :: int_value, start, finish
    real(dp) :: time_1, time_3

    real(dp), allocatable :: inter_int_1(:,:), inter_int_2(:), inter_int_3(:)
    real(dp), allocatable :: vec_1(:)

    call cpu_time(start)

    n_l = para%l+1

    allocate(TwoERadOrbInts(orb%n_max,orb%n_max, orb%n_max, orb%n_max, n_l, n_l, n_l, n_l, 2*para%l+1), stat=error)
    call allocerror(error)
    allocate(inter_int_1(para%ng,para%ng), stat=error)
    call allocerror(error)
    allocate(inter_int_2(para%ng), stat=error)
    call allocerror(error)
    allocate(inter_int_3(para%ng), stat=error)
    call allocerror(error)
    allocate(vec_1(para%ng), stat=error)
    call allocerror(error)

    inter_int_1 = zero
    inter_int_2 = zero
    inter_int_3 = zero
    TwoERadOrbInts = zero

    ! The transformation of 1e integrals is done in two steps
    ! i and j are basis indices, while m and n are orbital indices
    ! h_{ij} -> h_{mn} using transformation matrix b_{in}
    ! First step: h_{in} = \sum_j h_{ij} b_{in} 
    ! Second step: h_{mn} = \sum_i b^*_{im} h_{in} 

    ! The calculation is done over a loop of l

    do l = 1, 2*para%l + 1
      l_val = l - 1

      call cpu_time(time_1)

      ! The first step is done here

      do l2 = 1, n_l
        do m = 1, orb%n_max
          vec_1(:) = EigVecs(:,m,l2) 
          do j = 1, para%ng
            do i = 1, para%ng
               inter_int_1(i,j) = two_e_rad_int(i,j,l)*vec_1(i)
            end do 
          end do

          do l4 = 1, n_l
            do n = 1, orb%n_max
              vec_1(:) = EigVecs(:,n,l4) 
              do j = 1, para%ng
      
                int_value = 0.0d0
                do i = 1, para%ng
                  int_value = int_value + inter_int_1(i,j)*vec_1(i)
                end do 
                inter_int_2(j) = int_value
              end do

      ! The second step is done here, still it is inside the loop of l

              do l1 = 1, n_l
                do mp = 1, orb%n_max
                  vec_1(:) = EigVecs(:,mp,l1)

                  do i = 1, para%ng
                    inter_int_3(i) = vec_1(i)*inter_int_2(i)
                  end do

                  do l3 = 1, n_l
                    do np = 1, orb%n_max
                      vec_1(:) = EigVecs(:,np,l3)
                      int_value = 0.0d0

                      do i = 1, para%ng
                        int_value = int_value + vec_1(i)*inter_int_3(i)
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

      call cpu_time(time_3)

!     write(iout,*) 'Time taken in the first step = ', time_3 - time_1, 'seconds.'


    end do ! end loop over para%l

!     do m = 1, orb%n_max
!       do l1 = 1, n_l
!         do i = 1, para%ng
!     
!           int_value = 0.0d0
!           do j = 1, para%ng
!              int_value = int_value + two_e_rad_int(i,j,l)*EigVecs(j,m,l1)
!           end do 
!           inter_int_1(i,j,m,l1) = int_value
!         end do
!       end do 
!     end do ! end loop over n_max

    ! deallocate the large arrays after they are used
    deallocate(inter_int_1, inter_int_2, inter_int_3, vec_1)

    call cpu_time(finish)
    write(iout,'(X,a,f10.5,X,a)') 'Time taken for 2e transformation = ', finish-start, 'seconds.'

  end subroutine Calc2eRadOrbInts

  subroutine Calc2eRadOrbInts_alt_1()

    use DVRData, only : two_e_rad_int, para, eigen_vecs
    use OrbData, only : orb, TwoERadOrbInts

    integer  :: i, j, l, l1, l2, l3 ,l4, n_l, l_val, m, n, mp, np, error
    real(dp) :: int_value, start, finish, int_value_xc, val
    real(dp) :: time_1, time_2, time_3
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

      call cpu_time(time_1)

      ! The first step is done here
      do m = 1, orb%n_max
        do n = 1, orb%n_max
          do l1 = 1, n_l
            do l2 = 1, n_l
              do i = 1, para%ng
       
                int_value = 0.0d0
                do j = 1, para%ng
                  val = eigen_vecs(j,n,l2)*eigen_vecs(j,m,l1)  
!                 int_value = int_value + two_e_rad_int(i,j,l)*eigen_vecs(j,n,l2)*eigen_vecs(j,m,l1)
                  int_value = int_value + two_e_rad_int(i,j,l)*val
                end do 
                inter_int(i,m,n,l1,l2) = int_value
              end do
            end do 
          end do 
        end do
      end do ! end loop over n_max

      call cpu_time(time_2)

!     write(iout,*) 'Time taken in the first step = ', time_2 - time_1, 'seconds.'

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
!                       int_value_xc = int_value_xc + eigen_vecs(i,mp,l1)*eigen_vecs(i,n,l4)*inter_int(i,m,np,l2,l3)
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

      call cpu_time(time_3)

      write(iout,*) 'Time taken in the second step = ', time_3-time_2, 'seconds.'
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

  end subroutine Calc2eRadOrbInts_alt_1

end module CombineInts
