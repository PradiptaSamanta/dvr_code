module DVRIntAng

! use dvr_spline_mod
  use DVRData
  use dvr_diag_mod, only : allocerror
  use angular_utils
  use util_mod
  use constants
  
  implicit none

  contains
 
  subroutine GetAngularElements()

    integer   :: i, n_m1, n_m2
    logical   :: all_int, tSymInts
    complex(idp), allocatable  :: integrals_ang_prim(:,:,:,:,:)
    integer, allocatable :: l_interm(:)

    character(len=32)  :: file_name_syntx
    integer   :: opt

    sph_harm%n_l   = para%l + 1
    sph_harm%n_mp  = 2*para%l + 1

    ! The following three arrays will be required further in the subroutine.
    ! Although these arrays might look redundant in normal cases, they are very useful 
    ! (and therefore attempted) when 'm_l' values for 'l' are restricted to reduce the
    ! total number of orbitals.

    ! The array sph_harm%n_m stores the number 'm_l' quantum number for each 'l'+1
    ! e.g. n_m(1) = 1, n_m(2) = 3, n_m(3) = 5, while using the normal degeneracy for 'l'
    allocate(sph_harm%n_m(sph_harm%n_l))
    ! sph_harm%m_init stores 'min(m_l) - 1' for each 'l'+1
    ! e.g. m_init(1) = -1, m_init(2) = -2, m_init(3) = -3, while using the normal degeneracy for 'l'
    allocate(sph_harm%m_init(sph_harm%n_l))
    ! array to keep track the sum of the angular functions
    allocate(l_interm(sph_harm%n_l))

    n_m1 = 2*para%ml_max + 1
    do i = 1, sph_harm%n_l
      n_m2 = 2*i - 1
      sph_harm%n_m(i) = min(n_m1, n_m2)
      sph_harm%m_init(i) = -1*min(para%ml_max, i-1) - 1
    end do

    l_interm(1) = 0
    do i = 2, sph_harm%n_l
      n_m2 = 2*i - 1
      l_interm(i) = sum(sph_harm%n_m(1:i-1))
    end do

    ! Write the parameters corresponding to the spherical harmonics

    write(iout, *) '**********'
    write(iout, *) 'Setting up these parameters for the spherical harmonics:'
    write(iout, '(X,A,3X, I6)') 'sph_harm%n_l         =', sph_harm%n_l
    write(iout, '(X,A,3X, I6)') 'sph_harm%n_mp        =', sph_harm%n_mp
    write(iout, *) '***********' 


    ! If the eigenfunctions L_z operator being used, 2e integrals is no longer symmetric. 
    ! To recover this symmetry, the real eigenfunctions are used
    tSymInts = .true.

    all_int = .true.

    if (tSymInts) then
     
      ! The integrals are calculated in two steps. In the first step, the integrals are calculated 
      ! for the eigenfunctions of L_z. In the second step, these integrals are combined to form the
      ! the desired set of integrals

      ! Allocate the integrals to be calculated in the first step
      call allocate_int_ang(integrals_ang_prim, sph_harm)
      ! Allocate the final integrals 
      call allocate_int_ang(integrals_ang, sph_harm)
 
      ! First step
      write(iout, *) 'calculating the angular integrals'
!     call calc_int_angular_primitive(integrals_ang_prim, sph_harm) ! this subroutine is moved into the file DVRIntAng_adt.f90
      call calc_int_angular_main(integrals_ang_prim, sph_harm, l_interm)
 
      if (debug.gt.5) then
        file_name_syntx = 'ang_element_prim'
        write(iout, *) 'Writing down the angular integrals'
        call write_int_angular_real(integrals_ang_prim, sph_harm, all_int, file_name_syntx, l_interm)
      end if

      ! Second step

      ! Two alternative ways are used here to calculate the integrals in the combined basis

      opt = 1
      if ( opt.eq.1) then
        call calc_int_angular_combined(integrals_ang, integrals_ang_prim, sph_harm, l_interm)
!     elseif ( opt.eq.2 ) then
!       call calc_int_angular_combined_2(integrals_ang, integrals_ang_prim, sph_harm) ! this subroutine is moved into the file DVRIntAng_adt.f90
      else
        call stop_all('GetAngIntegral','Not an available option to calculate angular part of the 2e integrals')
      end if
 
      if (debug.gt.5) then
        file_name_syntx = 'ang_element_final'
        write(iout, *) 'Writing down the angular integrals'
        call write_int_angular_imaginary(integrals_ang, sph_harm, all_int, file_name_syntx, l_interm)
      end if

    else

      ! The integrals are calculated following the first step as mentioned above

      call allocate_int_ang(integrals_ang, sph_harm)
      
      write(iout, *) 'calculating the angular integrals'
      call calc_int_angular_main(integrals_ang, sph_harm, l_interm)

      if (debug.gt.4) then
        file_name_syntx = 'ang_element_prim'
        write(iout, *) 'Writing down the angular integrals'
        call write_int_angular_real(integrals_ang, sph_harm, all_int, file_name_syntx, l_interm)
      end if

    endif

    deallocate(l_interm)
    
  end subroutine GetAngularElements

  subroutine allocate_int_ang(integrals, sph_harm)

    type(sph_harm_t),  intent(in)           :: sph_harm
    complex(idp), allocatable,  intent(inout)  :: integrals(:,:,:,:,:)
 
    integer    :: n_l, n_mp, error, dim_l

    n_l  = sph_harm%n_l
    n_mp = sph_harm%n_mp

    dim_l = para%dim_l

    allocate(integrals(n_mp, dim_l, dim_l, dim_l, dim_l), stat=error)
    call allocerror(error)

    integrals = 0.0d0

  end subroutine allocate_int_ang

  subroutine calc_int_angular_main(integrals, sph_harm, l_interm)      
      
    type(sph_harm_t),     intent(in)       :: sph_harm
    complex(idp),         intent(inout)    :: integrals(:,:,:,:,:)
    integer, allocatable, intent(in)       :: l_interm(:)
 
    integer    :: n_l, n_mp, error, dim_l, dim_mp, k, q, q_init
    integer    :: l1, l2, l3, l4, l13, l24, n_m1, n_m2, n_m3, n_m4
    integer    :: m1, m2, m3, m4, m1_init, m2_init, m3_init, m4_init
    integer    :: lma, lmb, lmc, lmd, lmk, m_abq, m_sign

    real(idp)  :: lk, mq, la, lb, lc, ld, ma, mb, mc, md
    real(idp)  :: pre_fact_prod 
    real(idp)  :: w_symb_abcd, w_symb_abcd_q
    real(idp)  :: int_value

    real(idp), allocatable :: pre_fact_ac(:,:), pre_fact_bd(:,:)
    real(idp), allocatable :: w_symb_ac(:,:), w_symb_bd(:,:)
    real(idp), allocatable :: w_symb_ac_q(:,:,:), w_symb_bd_q(:,:,:)

!   real(idp), allocatable  ::  w_symb_ac_q, w_symb_bd_q


    n_l  = sph_harm%n_l
    n_mp = sph_harm%n_mp

    dim_l = para%dim_l
    dim_mp = n_mp**2

    allocate(pre_fact_ac(n_l,n_l), pre_fact_bd(n_l,n_l))
    allocate(w_symb_ac(n_l,n_l), w_symb_bd(n_l,n_l))
    allocate(w_symb_ac_q(dim_mp, dim_l, dim_l), w_symb_bd_q(dim_mp, dim_l, dim_l))

    do k = 1, n_mp
      lk = dfloat(k-1)
      q_init = -1*k

      do l1 = 1, n_l
        la = dfloat(l1 - 1)
        n_m1 = sph_harm%n_m(l1)
        m1_init = sph_harm%m_init(l1)
        !n_m1 = 2*l1 - 1
        !m1_init = -1*l1
        do l3 = 1, n_l
          lc = dfloat(l3 - 1)
          n_m3 = sph_harm%n_m(l3)
          m3_init = sph_harm%m_init(l3)
          !n_m3 = 2*l3 - 1
          !m3_init = -1*l3

          pre_fact_ac(l1,l3) = sqrt((2.0d0*la)+1.0d0)*sqrt((2.0d0*lc)+1.0d0)

          w_symb_ac(l1,l3) = wigner3j(lk, la, lc, 0.0d0, 0.0d0, 0.0d0)

          do m1 = 1, n_m1
            ma = dfloat(m1_init + m1)
            lma = l_interm(l1) + m1
            !lma = (l1-1)**2 + int(la) + int(ma) + 1
            do m3 = 1, n_m3
              mc = dfloat(m3_init + m3)
              lmc = l_interm(l3) + m3
              !lmc = (l3-1)**2 + int(lc) + int(mc) + 1

              do q = 1, 2*k-1
!               if (lma.eq.1.and.lmb.eq.1.and.lmc.eq.2.and.lmd.eq.4) then
                mq = dfloat(q_init + q)

!               lmk = (lk-1)**2 + q  !! CHECK IF IT IS A BUG
                lmk = (k-1)**2 + q
                m_abq = int(ma + mc + mq)
                m_sign = (-1)**m_abq
!           write(*,*) 'k, lk, q, mq, m_abq:',  k, int(lk), q, int(mq), m_sign

                w_symb_ac_q(lmk, lma, lmc) =                                 &
                &     wigner3j(lk, la, lc, mq, -1.0d0*ma, mc)

              end do
            end do
          enddo  
        enddo
      enddo

      do l2 = 1, n_l
        lb = dfloat(l2 - 1)
        n_m2 = sph_harm%n_m(l2)
        m2_init = sph_harm%m_init(l2)
        !n_m2 = 2*l2 - 1
        !m2_init = -1*l2
        do l4 = 1, n_l
          ld = dfloat(l4 - 1)
          n_m4 = sph_harm%n_m(l4)
          m4_init = sph_harm%m_init(l4)
          !n_m4 = 2*l4 - 1
          !m4_init = -1*l4

          pre_fact_bd(l2,l4) = sqrt((2.0d0*lb)+1.0d0)*sqrt((2.0d0*ld)+1.0d0)

          w_symb_bd(l2, l4) = wigner3j(lk, lb, ld, 0.0d0, 0.0d0, 0.0d0)

!         pre_fact_prod = pre_fact_ac * pre_fact_bd
!         w_symb_abcd = w_symb_ac * w_symb_bd

          do m2 = 1, n_m2
            mb = dfloat(m2_init + m2)
            lmb = l_interm(l2) + m2
            !lmb = (l2-1)**2 + int(lb) + int(mb) + 1
            do m4 = 1, n_m4
              md = dfloat(m4_init + m4)
              lmd = l_interm(l4) + m4
              !lmd = (l4-1)**2 + int(ld) + int(md) + 1

              int_value = 0.0d0
              do q = 1, 2*k-1
!               if (lma.eq.1.and.lmb.eq.1.and.lmc.eq.2.and.lmd.eq.4) then
                mq = dfloat(q_init + q) 
!               lmk = (lk-1)**2 + q  !! CHECK IF IT IS A BUG
                lmk = (k-1)**2 + q

                m_abq = int(mb + md + mq)
                m_sign = (-1)**m_abq
!           write(*,*) 'k, lk, q, mq, m_abq:',  k, int(lk), q, int(mq), m_sign

                w_symb_bd_q(lmk, lmb, lmd) =                                          &
                        &     wigner3j(lk, lb, ld, -1.0d0*mq, -1.0d0*mb,md)
!           write(*,*) lk, lb, ld, w_symb_bd_q




!             end if
              end do ! end loop for q

            end do ! end loop for m4
          end do ! end loop for m2   
        end do  ! end loop for l4
      end do  ! end loop for l2

      do l1 = 1, n_l
        n_m1 = sph_harm%n_m(l1)
        m1_init = sph_harm%m_init(l1)
        !n_m1 = 2*l1 - 1
        !m1_init = -1*l1
        do l2 = 1, n_l 
          n_m2 = sph_harm%n_m(l2)
          m2_init = sph_harm%m_init(l2)
          !n_m2 = 2*l2 - 1
          !m2_init = -1*l2
          do l3 = 1, n_l 
            n_m3 = sph_harm%n_m(l3)
            !n_m3 = 2*l3 - 1
            do l4 = 1, n_l 
              n_m4 = sph_harm%n_m(l4)
              !n_m4 = 2*l4 - 1
              pre_fact_prod = pre_fact_ac(l1, l3) * pre_fact_bd(l2, l4)
              w_symb_abcd = w_symb_ac(l1,l3)*w_symb_bd(l2,l4)
              do m1 = 1, n_m1
                ma = m1_init + m1
                lma = l_interm(l1) + m1
                !lma = (l1-1)**2 + m1
                do m2 = 1, n_m2
                  mb = m2_init + m2
                  lmb = l_interm(l2) + m2
                  !lmb = (l2-1)**2 + m2
                  do m3 = 1, n_m3
                    lmc = l_interm(l3) + m3
                    !lmc = (l3-1)**2 + m3
                    do m4 = 1, n_m4
!                   m4 = (-l1 -l2 +l3 +l4) + (m1 +m2 -m3)  ! following the condition m1-m3 = m4-m2 transformed into the index we are using here for m
!                   if (0.lt.m4.and.m4.le.n_m4) then
                      lmd = l_interm(l4) + m4
                      !lmd = (l4-1)**2 + m4
                      int_value = 0.0d0
                      do q = 1, 2*k-1
!               if (lma.eq.1.and.lmb.eq.1.and.lmc.eq.2.and.lmd.eq.4) then
                        mq = dfloat(q_init + q) 
!                       lmk = (lk-1)**2 + q  !! ANOTHER PLACE TO CHECK THE BUG
                        lmk = (k-1)**2 + q
   
                        m_abq = int(ma + mb + mq)
                        m_sign = (-1)**m_abq
!                       if (lma.eq.1.and.lmb.eq.2.and.lmc.eq.1.and.lmd.eq.2) &
!                       &   write(79,'(3i4)') int(mq), m_abq, m_sign
                        w_symb_abcd_q = w_symb_ac_q(lmk, lma, lmc)        &
                        &              *w_symb_bd_q(lmk, lmb, lmd)
                        int_value = int_value +                                &
                  &     m_sign * pre_fact_prod  * w_symb_abcd * w_symb_abcd_q
                        integrals(k, lma, lmb, lmc, lmd) = int_value
!                       if (lma.eq.1.and.lmb.eq.4.and.lmc.eq.1.and.lmd.eq.4) &
!                       &   write(80,'(i4,4f15.8)') m_sign, pre_fact_prod, w_symb_abcd, w_symb_abcd_q, real(int_value)
                      end do
                    end do
!                   end if
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
!               write(*,*) q, int_value, w_symb_ac_q, w_symb_bd_q
    end do  ! end loop for k 

  end subroutine calc_int_angular_main

  subroutine calc_int_angular_combined(integrals, integrals_prim, sph_harm, l_interm)

    complex(idp),      intent(inout)        :: integrals(:,:,:,:,:)
    complex(idp),      intent(in)           :: integrals_prim(:,:,:,:,:)
    type(sph_harm_t),  intent(in)           :: sph_harm
    integer, allocatable, intent(in)        :: l_interm(:)

    integer :: n_l, n_mp, dim_l, dim_mp
    integer :: k, l1, l2, l3, l4, m1, m2, m3, m4, l10, l20, l30, l40
    integer :: lm1, lm2, lm3, lm4, lm1_p, lm2_p, lm3_p, lm4_p
    integer :: n_m1, n_m2, n_m3, n_m4, ma, mb, mc, md
    integer, allocatable :: n_m(:), l_finish(:), pos_m_zero(:)
    complex(idp) :: int_1, int_2, int_3, int_4, int_5, int_6, int_7, int_8
    complex(idp) :: int_9, int_10, int_11, int_12, int_13, int_14, int_15, int_16
    complex(idp) :: prim_fac_1, prim_fac_2
    complex(idp) :: fac_1, fac_2, fac_3, fac_4, fac_5, fac_6, fac_7, fac_8
    complex(idp) :: fac_9, fac_10, fac_11, fac_12, fac_13, fac_14, fac_15, fac_16
    real(idp), parameter :: val_1 = (1.0d0*sqrt_half)
    real(idp), parameter :: val_2 = (-1.0d0*sqrt_half)
    real(idp)    :: msign_1, msign_2, msign_3, msign_4, msign_5, msign_6, msign_7, msign_8
    real(idp)    :: msign_9, msign_10, msign_11, msign_12, msign_13, msign_14, msign_15, msign_16

    n_l  = sph_harm%n_l
    n_mp = sph_harm%n_mp

    dim_l = para%dim_l
    dim_mp = n_mp**2

    prim_fac_1 = (0.0d0, val_1)
    prim_fac_2 = (val_1,0.0d0)
    !prim_fac_2 = (val_1, 0.0d0)

    ! The following three arrays are needed to store numbers which will be required
    ! during the further loops. These are stored as separated arrays, instead of calculating them
    ! inside the loops, because calculations of these numbers are more involved when not all 'm_l' are used for 'l'. 

    ! n_m stores the number of 'm_l' quantum number =< 0 for each 'l'+1 
    ! e.g. n_m(1) = 1, only 0; n_m(2) = 2, -1 and 0; n_m(3) = -2, -1, 0, while using the normal degeneracy for 'l'
    allocate(n_m(n_l))
    ! l_finish stores the number of spherical harmonics exist upto each 'l'+1
    ! e.g. l_finish(1) = 1; l_finish(2) = 4; l_finish(3) = 9, while using the normal degeneracy for 'l'
    allocate(l_finish(n_l))
    ! pos_m_zero stores the position of the m_l=0 in the distribution of all 'm_l' for a particular 'l'+1
    ! e.g. pos_m_zero(1)=1; pos_m_zero(2)=2; pos_m_zero=3; while using the normal degeneracy for 'l' 
    allocate(pos_m_zero(n_l))

    ! These new arrays are calculated here
    do k = 1, n_l
      n_m(k) = min(k,para%ml_max+1)
      l_finish(k) = sum(sph_harm%n_m(1:k))
      pos_m_zero(k) = min(k,para%ml_max+1)
    end do
    
!   write(iout, *) n_m
!   write(iout, *) l_finish
!   write(iout, *) pos_m_zero

    do k = 1, n_mp
      do l1 = 1, n_l
        n_m1 = n_m(l1) ! Note the use of the 'n_m' array
        l10 = pos_m_zero(l1) ! Note the use of the 'pos_m_zero' array
        do l3 = 1, n_l
          n_m3 = n_m(l3)
          l30 = pos_m_zero(l3)
          do l2 = 1, n_l
            n_m2 = n_m(l2)
            l20 = pos_m_zero(l2)
            do l4 = 1, n_l
              n_m4 = n_m(l4)
              l40 = pos_m_zero(l4)
              do m1 = 1, n_m1
                lm1 = l_interm(l1) + m1
                if (m1.eq.l10) then
 
                  do m3 = 1, n_m3
                    lm3 = l_interm(l3) + m3
                    if (m3.eq.l30) then
                      do m2 = 1, n_m2
                        lm2 = l_interm(l2) + m2
                        if (m2.eq.l20) then
                          do m4 = 1, n_m4
                 
                            lm4 = l_interm(l4) + m4
                            int_1 = integrals_prim(k, lm1, lm2, lm3, lm4)
                 
                            if (m4.eq.l40) then
                              integrals(k, lm1, lm2, lm3, lm4) = int_1
                            else
                              lm4_p = l_finish(l4) - m4 + 1 ! Note the use of the 'l_finish' array
                              md = sph_harm%m_init(l4) + m4
                              msign_1 = m_one**md
                              int_2 = msign_1 * integrals_prim(k, lm1, lm2, lm3, lm4_p)
                              fac_1 = prim_fac_1
                              fac_2 = prim_fac_2
                              integrals(k, lm1, lm2, lm3, lm4  ) = fac_1 * (int_1 - int_2)
                              integrals(k, lm1, lm2, lm3, lm4_p) = fac_2 * (int_1 + int_2)
                            end if
                          end do
                        else
!                         lm2_p =  (l2-1)**2 + m2+l2
                          lm2_p = l_finish(l2) - m2 + 1
                          !mb = -1*l2 + m2
                          mb = sph_harm%m_init(l2) + m2
                          msign_1 = m_one**mb
                          do m4 = 1, n_m4
                            lm4 = l_interm(l4) + m4
                            int_1 =           integrals_prim(k, lm1, lm2  , lm3, lm4)
                            int_2 = msign_1 * integrals_prim(k, lm1, lm2_p, lm3, lm4)
                            if (m4.eq.l40) then
                              fac_1 = dconjg(prim_fac_1)
                              fac_2 = dconjg(prim_fac_2)
                              integrals(k, lm1, lm2  , lm3, lm4) = fac_1 * (int_1 - int_2)
                              integrals(k, lm1, lm2_p, lm3, lm4) = fac_2 * (int_1 + int_2)
                            else
                              lm4_p = l_finish(l4) - m4 + 1
                              !md = -1*l4 + m4
                              md = sph_harm%m_init(l4) + m4
                              msign_2 = m_one**(      md )
                              msign_3 = m_one**( mb + md )
                              int_3 = msign_2 * integrals_prim(k, lm1, lm2  , lm3, lm4_p)
                              int_4 = msign_3 * integrals_prim(k, lm1, lm2_p, lm3, lm4_p)
                              fac_1 = dconjg(prim_fac_1)*prim_fac_1
                              fac_2 = dconjg(prim_fac_2)*prim_fac_1
                              fac_3 = dconjg(prim_fac_1)*prim_fac_2
                              fac_4 = dconjg(prim_fac_2)*prim_fac_2
                              integrals(k, lm1, lm2  , lm3, lm4  ) = fac_1 * (int_1 - int_2 - int_3 + int_4)
                              integrals(k, lm1, lm2_p, lm3, lm4  ) = fac_2 * (int_1 + int_2 - int_3 - int_4)
                              integrals(k, lm1, lm2  , lm3, lm4_p) = fac_3 * (int_1 - int_2 + int_3 - int_4)
                              integrals(k, lm1, lm2_p, lm3, lm4_p) = fac_4 * (int_1 + int_2 + int_3 + int_4)
                            end if
                          end do
                        end if
                      end do
                    else
                      lm3_p = l_finish(l3) - m3 + 1
                      !mc = -1*l3 + m3
                      mc = sph_harm%m_init(l3) + m3
                      msign_1 = m_one**mc
                      do m2 = 1, n_m2
                        lm2 = l_interm(l2) + m2
                        if (m2.eq.l20) then
                          do m4 = 1, n_m4
                 
                            lm4 = l_interm(l4) + m4
                            int_1 =           integrals_prim(k, lm1, lm2, lm3  , lm4)
                            int_2 = msign_1 * integrals_prim(k, lm1, lm2, lm3_p, lm4)
                 
                            if (m4.eq.l40) then
                              fac_1 = prim_fac_1
                              fac_2 = prim_fac_2
                              integrals(k, lm1, lm2, lm3  , lm4) = fac_1 * (int_1 - int_2)
                              integrals(k, lm1, lm2, lm3_p, lm4) = fac_2 * (int_1 + int_2)
                            else
                              lm4_p = l_finish(l4) - m4 + 1
                              !md = -1*l4 + m4
                              md = sph_harm%m_init(l4) + m4
                              msign_2 = m_one**(      md )
                              msign_3 = m_one**( mc + md )
                              int_3 = msign_2 * integrals_prim(k, lm1, lm2, lm3  , lm4_p)
                              int_4 = msign_3 * integrals_prim(k, lm1, lm2, lm3_p, lm4_p)
                              fac_1 = prim_fac_1*prim_fac_1
                              fac_2 = prim_fac_2*prim_fac_1
                              fac_3 = prim_fac_2*prim_fac_2
                              integrals(k, lm1, lm2, lm3  , lm4  ) = fac_1 * (int_1 - int_2 - int_3 + int_4)
                              integrals(k, lm1, lm2, lm3_p, lm4  ) = fac_2 * (int_1 + int_2 - int_3 - int_4)
                              integrals(k, lm1, lm2, lm3  , lm4_p) = fac_2 * (int_1 - int_2 + int_3 - int_4)
                              integrals(k, lm1, lm2, lm3_p, lm4_p) = fac_3 * (int_1 + int_2 + int_3 + int_4)
!                             write(91,'(4i4,3x,5f15.8)') lm1, lm2, lm3_p, lm4_p, real(int_1), real(int_2), real(int_3), real(int_4), real(integrals(k, lm1, lm2, lm3_p, lm4_p))
                            end if
                          end do
                        else
                          lm2_p = l_finish(l2) - m2 + 1
                          !mb = -1*l2 + m2
                          mb = sph_harm%m_init(l2) + m2
                          msign_2 = m_one**( mb      )
                          msign_3 = m_one**( mb + mc )
                          do m4 = 1, n_m4
                            lm4 = l_interm(l4) + m4
                            int_1 =           integrals_prim(k, lm1, lm2  , lm3  , lm4)
                            int_2 = msign_2 * integrals_prim(k, lm1, lm2_p, lm3  , lm4)
                            int_3 = msign_1 * integrals_prim(k, lm1, lm2  , lm3_p, lm4)
                            int_4 = msign_3 * integrals_prim(k, lm1, lm2_p, lm3_p, lm4)
                            if (m4.eq.l40) then
                              fac_1 = dconjg(prim_fac_1)*prim_fac_1
                              fac_2 = dconjg(prim_fac_2)*prim_fac_1
                              fac_3 = dconjg(prim_fac_1)*prim_fac_2
                              fac_4 = dconjg(prim_fac_2)*prim_fac_2
                              integrals(k, lm1, lm2  , lm3  , lm4) = fac_1*(int_1 - int_2 - int_3 + int_4)
                              integrals(k, lm1, lm2_p, lm3  , lm4) = fac_2*(int_1 + int_2 - int_3 - int_4)
                              integrals(k, lm1, lm2  , lm3_p, lm4) = fac_3*(int_1 - int_2 + int_3 - int_4)
                              integrals(k, lm1, lm2_p, lm3_p, lm4) = fac_4*(int_1 + int_2 + int_3 + int_4)
                            else
                              lm4_p = l_finish(l4) - m4 + 1
                              !md = -1*l4 + m4
                              md = sph_harm%m_init(l4) + m4
                              msign_4 = m_one**( md           )
                              msign_5 = m_one**( mb + md      )
                              msign_6 = m_one**( mc + md      )
                              msign_7 = m_one**( mb + mc + md )
                              int_5 = msign_4 * integrals_prim(k, lm1, lm2  , lm3  , lm4_p)
                              int_6 = msign_5 * integrals_prim(k, lm1, lm2_p, lm3  , lm4_p)
                              int_7 = msign_6 * integrals_prim(k, lm1, lm2  , lm3_p, lm4_p)
                              int_8 = msign_7 * integrals_prim(k, lm1, lm2_p, lm3_p, lm4_p)
                              fac_1 = dconjg(prim_fac_1)*prim_fac_1*prim_fac_1
                              fac_2 = dconjg(prim_fac_2)*prim_fac_1*prim_fac_1
                              fac_3 = dconjg(prim_fac_1)*prim_fac_2*prim_fac_1
                              fac_4 = dconjg(prim_fac_2)*prim_fac_2*prim_fac_1
                              fac_5 = dconjg(prim_fac_1)*prim_fac_1*prim_fac_2
                              fac_6 = dconjg(prim_fac_2)*prim_fac_1*prim_fac_2
                              fac_7 = dconjg(prim_fac_1)*prim_fac_2*prim_fac_2
                              fac_8 = dconjg(prim_fac_2)*prim_fac_2*prim_fac_2
                              integrals(k, lm1, lm2,   lm3,   lm4  ) = fac_1*(int_1 - int_2 - int_3 + int_4 - int_5 + int_6 + int_7 - int_8)
                              integrals(k, lm1, lm2_p, lm3,   lm4  ) = fac_2*(int_1 + int_2 - int_3 - int_4 - int_5 - int_6 + int_7 + int_8)
                              integrals(k, lm1, lm2,   lm3_p, lm4  ) = fac_3*(int_1 - int_2 + int_3 - int_4 - int_5 + int_6 - int_7 + int_8)
                              integrals(k, lm1, lm2_p, lm3_p, lm4  ) = fac_4*(int_1 + int_2 + int_3 + int_4 - int_5 - int_6 - int_7 - int_8)
                              integrals(k, lm1, lm2,   lm3,   lm4_p) = fac_5*(int_1 - int_2 - int_3 + int_4 + int_5 - int_6 - int_7 + int_8)
                              integrals(k, lm1, lm2_p, lm3,   lm4_p) = fac_6*(int_1 + int_2 - int_3 - int_4 + int_5 + int_6 - int_7 - int_8)
                              integrals(k, lm1, lm2,   lm3_p, lm4_p) = fac_7*(int_1 - int_2 + int_3 - int_4 + int_5 - int_6 + int_7 - int_8)
                              integrals(k, lm1, lm2_p, lm3_p, lm4_p) = fac_8*(int_1 + int_2 + int_3 + int_4 + int_5 + int_6 + int_7 + int_8)
                            end if
                          end do
                        end if
                      end do
                    end if
                  end do
                else
                  lm1_p = l_finish(l1) - m1 + 1
                  !ma = -1*l1 + m1
                  ma = sph_harm%m_init(l1) + m1
                  msign_1 = m_one**ma
                  do m3 = 1, n_m3
                    lm3 = l_interm(l3) + m3
                    if (m3.eq.l30) then
                      do m2 = 1, n_m2
                        lm2 = l_interm(l2) + m2
                        if (m2.eq.l20) then
                          do m4 = 1, n_m4
                  
                            lm4 = l_interm(l4) + m4
                            int_1 =            integrals_prim(k, lm1  , lm2, lm3, lm4)
                            int_2 = msign_1  * integrals_prim(k, lm1_p, lm2, lm3, lm4)
                  
                            if (m4.eq.l40) then
                              fac_1 =  dconjg(prim_fac_1)
                              fac_2 =  dconjg(prim_fac_2)
                              integrals(k, lm1  , lm2, lm3, lm4) = fac_1 * (int_1 - int_2)
                              integrals(k, lm1_p, lm2, lm3, lm4) = fac_2 * (int_1 + int_2)
                            else
                              lm4_p = l_finish(l4) - m4 + 1
                              !md = -1*l4 + m4
                              md = sph_harm%m_init(l4) + m4
                              msign_2 = m_one**( md      )
                              msign_3 = m_one**( ma + md )
                              int_3 = msign_2 * integrals_prim(k, lm1  , lm2, lm3, lm4_p)
                              int_4 = msign_3 * integrals_prim(k, lm1_p, lm2, lm3, lm4_p)
                              fac_1 = dconjg(prim_fac_1)*prim_fac_1
                              fac_2 = dconjg(prim_fac_2)*prim_fac_1
                              fac_3 = dconjg(prim_fac_1)*prim_fac_2
                              fac_4 = dconjg(prim_fac_2)*prim_fac_2
                              integrals(k, lm1  , lm2, lm3, lm4  ) = fac_1 * (int_1 - int_2 - int_3 + int_4)
                              integrals(k, lm1_p, lm2, lm3, lm4  ) = fac_2 * (int_1 + int_2 - int_3 - int_4)
                              integrals(k, lm1  , lm2, lm3, lm4_p) = fac_3 * (int_1 - int_2 + int_3 - int_4)
                              integrals(k, lm1_p, lm2, lm3, lm4_p) = fac_4 * (int_1 + int_2 + int_3 + int_4)
                            end if
                          end do
                        else
                          lm2_p = l_finish(l2) - m2 + 1
                          !mb = -1*l2 + m2
                          mb = sph_harm%m_init(l2) + m2
                          msign_2 = m_one**(      mb )
                          msign_3 = m_one**( ma + mb )
                          do m4 = 1, n_m4
                            lm4 = l_interm(l4) + m4
                            int_1 =           integrals_prim(k, lm1  , lm2  , lm3, lm4)
                            int_2 = msign_1 * integrals_prim(k, lm1_p, lm2  , lm3, lm4)
                            int_3 = msign_2 * integrals_prim(k, lm1  , lm2_p, lm3, lm4)
                            int_4 = msign_3 * integrals_prim(k, lm1_p, lm2_p, lm3, lm4)
                            if (m4.eq.l40) then
                              fac_1 = dconjg(prim_fac_1)*dconjg(prim_fac_1)
                              fac_2 = dconjg(prim_fac_2)*dconjg(prim_fac_1)
                              fac_3 = dconjg(prim_fac_1)*dconjg(prim_fac_2)
                              fac_4 = dconjg(prim_fac_2)*dconjg(prim_fac_2)
                              integrals(k, lm1  , lm2  , lm3, lm4) = fac_1 * (int_1 - int_2 - int_3 + int_4)
                              integrals(k, lm1_p, lm2  , lm3, lm4) = fac_2 * (int_1 + int_2 - int_3 - int_4)
                              integrals(k, lm1  , lm2_p, lm3, lm4) = fac_3 * (int_1 - int_2 + int_3 - int_4)
                              integrals(k, lm1_p, lm2_p, lm3, lm4) = fac_4 * (int_1 + int_2 + int_3 + int_4)
                            else
                              lm4_p = l_finish(l4) - m4 + 1
                              !md = -1*l4 + m4
                              md = sph_harm%m_init(l4) + m4
                              msign_4 = m_one**(           md )
                              msign_5 = m_one**( ma      + md )
                              msign_6 = m_one**(      mb + md )
                              msign_7 = m_one**( ma + mb + md )
                              int_5 = msign_4 * integrals_prim(k, lm1  , lm2  , lm3, lm4_p)
                              int_6 = msign_5 * integrals_prim(k, lm1_p, lm2  , lm3, lm4_p)
                              int_7 = msign_6 * integrals_prim(k, lm1  , lm2_p, lm3, lm4_p)
                              int_8 = msign_7 * integrals_prim(k, lm1_p, lm2_p, lm3, lm4_p)
                              fac_1 = dconjg(prim_fac_1)*dconjg(prim_fac_1)*prim_fac_1
                              fac_2 = dconjg(prim_fac_2)*dconjg(prim_fac_1)*prim_fac_1
                              fac_3 = dconjg(prim_fac_1)*dconjg(prim_fac_2)*prim_fac_1
                              fac_4 = dconjg(prim_fac_2)*dconjg(prim_fac_2)*prim_fac_1
                              fac_5 = dconjg(prim_fac_1)*dconjg(prim_fac_1)*prim_fac_2
                              fac_6 = dconjg(prim_fac_2)*dconjg(prim_fac_1)*prim_fac_2
                              fac_7 = dconjg(prim_fac_1)*dconjg(prim_fac_2)*prim_fac_2
                              fac_8 = dconjg(prim_fac_2)*dconjg(prim_fac_2)*prim_fac_2
                              integrals(k, lm1  , lm2  , lm3, lm4  ) = fac_1 * (int_1 - int_2 - int_3 + int_4 - int_5 + int_6 + int_7 - int_8)
                              integrals(k, lm1_p, lm2  , lm3, lm4  ) = fac_2 * (int_1 + int_2 - int_3 - int_4 - int_5 - int_6 + int_7 + int_8) 
                              integrals(k, lm1  , lm2_p, lm3, lm4  ) = fac_3 * (int_1 - int_2 + int_3 - int_4 - int_5 + int_6 - int_7 + int_8)
                              integrals(k, lm1_p, lm2_p, lm3, lm4  ) = fac_4 * (int_1 + int_2 + int_3 + int_4 - int_5 - int_6 - int_7 - int_8)
                              integrals(k, lm1  , lm2  , lm3, lm4_p) = fac_5 * (int_1 - int_2 - int_3 + int_4 + int_5 - int_6 - int_7 + int_8)
                              integrals(k, lm1_p, lm2  , lm3, lm4_p) = fac_6 * (int_1 + int_2 - int_3 - int_4 + int_5 + int_6 - int_7 - int_8)
                              integrals(k, lm1  , lm2_p, lm3, lm4_p) = fac_7 * (int_1 - int_2 + int_3 - int_4 + int_5 - int_6 + int_7 - int_8)
                              integrals(k, lm1_p, lm2_p, lm3, lm4_p) = fac_8 * (int_1 + int_2 + int_3 + int_4 + int_5 + int_6 + int_7 + int_8)
                            end if
                          end do
                        end if
                      end do
                    else
                      lm3_p = l_finish(l3) - m3 + 1
                      !mc = -1*l3 + m3
                      mc = sph_harm%m_init(l3) + m3
                      msign_2 = m_one**(      mc )
                      msign_3 = m_one**( ma + mc )
                      do m2 = 1, n_m2
                        lm2 = l_interm(l2) + m2
                        if (m2.eq.l20) then
                          do m4 = 1, n_m4
                  
                            lm4 = l_interm(l4) + m4
                            int_1 = integrals_prim(k, lm1  , lm2, lm3  , lm4)
                            int_2 = msign_1 * integrals_prim(k, lm1_p, lm2, lm3  , lm4)
                            int_3 = msign_2 * integrals_prim(k, lm1  , lm2, lm3_p, lm4)
                            int_4 = msign_3 * integrals_prim(k, lm1_p, lm2, lm3_p, lm4)
                  
                            if (m4.eq.l40) then
                              fac_1 = dconjg(prim_fac_1)*prim_fac_1
                              fac_2 = dconjg(prim_fac_2)*prim_fac_1
                              fac_3 = dconjg(prim_fac_1)*prim_fac_2
                              fac_4 = dconjg(prim_fac_2)*prim_fac_2
                              integrals(k, lm1  , lm2, lm3  , lm4) = fac_1 * (int_1 - int_2 - int_3 + int_4)
                              integrals(k, lm1_p, lm2, lm3  , lm4) = fac_2 * (int_1 + int_2 - int_3 - int_4) 
                              integrals(k, lm1  , lm2, lm3_p, lm4) = fac_3 * (int_1 - int_2 + int_3 - int_4)
                              integrals(k, lm1_p, lm2, lm3_p, lm4) = fac_4 * (int_1 + int_2 + int_3 + int_4)
!                             write(92,'(4i4,3x,5f15.8)') lm1_p, lm2, lm3_p, lm4, real(int_1), real(int_2), real(int_3), real(int_4), real(integrals(k, lm1_p, lm2, lm3_p, lm4))
!                             write(93,'(4i4,3x,5f15.8)') lm1, lm2, lm3, lm4, real(int_1), real(int_2), real(int_3), real(int_4), real(integrals(k, lm1, lm2, lm3, lm4))
                            else
                              lm4_p = l_finish(l4) - m4 + 1
                              !md = -1*l4 + m4
                              md = sph_harm%m_init(l4) + m4
                              msign_4 = m_one**(           md )
                              msign_5 = m_one**( ma      + md )
                              msign_6 = m_one**(      mc + md )
                              msign_7 = m_one**( ma + mc + md )
                              int_5 = msign_4 * integrals_prim(k, lm1  , lm2, lm3  , lm4_p)
                              int_6 = msign_5 * integrals_prim(k, lm1_p, lm2, lm3  , lm4_p)
                              int_7 = msign_6 * integrals_prim(k, lm1  , lm2, lm3_p, lm4_p)
                              int_8 = msign_7 * integrals_prim(k, lm1_p, lm2, lm3_p, lm4_p)
                              fac_1 = dconjg(prim_fac_1)*prim_fac_1*prim_fac_1
                              fac_2 = dconjg(prim_fac_2)*prim_fac_1*prim_fac_1
                              fac_3 = dconjg(prim_fac_1)*prim_fac_2*prim_fac_1
                              fac_4 = dconjg(prim_fac_2)*prim_fac_2*prim_fac_1
                              fac_5 = dconjg(prim_fac_1)*prim_fac_1*prim_fac_2
                              fac_6 = dconjg(prim_fac_2)*prim_fac_1*prim_fac_2
                              fac_7 = dconjg(prim_fac_1)*prim_fac_2*prim_fac_2
                              fac_8 = dconjg(prim_fac_2)*prim_fac_2*prim_fac_2
                              integrals(k, lm1  , lm2, lm3  , lm4  ) = fac_1 * (int_1 - int_2 - int_3 + int_4 - int_5 + int_6 + int_7 - int_8)
                              integrals(k, lm1_p, lm2, lm3  , lm4  ) = fac_2 * (int_1 + int_2 - int_3 - int_4 - int_5 - int_6 + int_7 + int_8)
                              integrals(k, lm1  , lm2, lm3_p, lm4  ) = fac_3 * (int_1 - int_2 + int_3 - int_4 - int_5 + int_6 - int_7 + int_8)
                              integrals(k, lm1_p, lm2, lm3_p, lm4  ) = fac_4 * (int_1 + int_2 + int_3 + int_4 - int_5 - int_6 - int_7 - int_8)
                              integrals(k, lm1  , lm2, lm3  , lm4_p) = fac_5 * (int_1 - int_2 - int_3 + int_4 + int_5 - int_6 - int_7 + int_8)
                              integrals(k, lm1_p, lm2, lm3  , lm4_p) = fac_6 * (int_1 + int_2 - int_3 - int_4 + int_5 + int_6 - int_7 - int_8)
                              integrals(k, lm1  , lm2, lm3_p, lm4_p) = fac_7 * (int_1 - int_2 + int_3 - int_4 + int_5 - int_6 + int_7 - int_8)
                              integrals(k, lm1_p, lm2, lm3_p, lm4_p) = fac_8 * (int_1 + int_2 + int_3 + int_4 + int_5 + int_6 + int_7 + int_8)
                            end if
                          end do
                        else
                          lm2_p = l_finish(l2) - m2 + 1
                          !mb = -1*l2 + m2
                          mb = sph_harm%m_init(l2) + m2
                          ! msign_2 = m_one**( mc    )  This is defined above
                          ! msign_3 = m_one**( ma + mc ) This is defined above
                          msign_4 = m_one**(      mb      )
                          msign_5 = m_one**( ma + mb      )
                          msign_6 = m_one**(      mb + mc )
                          msign_7 = m_one**( ma + mb + mc )
                          do m4 = 1, n_m4
                            lm4 = l_interm(l4) + m4
                            int_1 = integrals_prim(k, lm1  , lm2  , lm3  , lm4)
                            int_2 = msign_1 * integrals_prim(k, lm1_p, lm2  , lm3  , lm4)
                            int_3 = msign_4 * integrals_prim(k, lm1  , lm2_p, lm3  , lm4)
                            int_4 = msign_5 * integrals_prim(k, lm1_p, lm2_p, lm3  , lm4)
                            int_5 = msign_2 * integrals_prim(k, lm1  , lm2  , lm3_p, lm4)
                            int_6 = msign_3 * integrals_prim(k, lm1_p, lm2  , lm3_p, lm4)
                            int_7 = msign_6 * integrals_prim(k, lm1  , lm2_p, lm3_p, lm4)
                            int_8 = msign_7 * integrals_prim(k, lm1_p, lm2_p, lm3_p, lm4)
                            if (m4.eq.l40) then
                              fac_1 = dconjg(prim_fac_1)*dconjg(prim_fac_1)*prim_fac_1
                              fac_2 = dconjg(prim_fac_2)*dconjg(prim_fac_1)*prim_fac_1
                              fac_3 = dconjg(prim_fac_1)*dconjg(prim_fac_2)*prim_fac_1
                              fac_4 = dconjg(prim_fac_2)*dconjg(prim_fac_2)*prim_fac_1
                              fac_5 = dconjg(prim_fac_1)*dconjg(prim_fac_1)*prim_fac_2
                              fac_6 = dconjg(prim_fac_2)*dconjg(prim_fac_1)*prim_fac_2
                              fac_7 = dconjg(prim_fac_1)*dconjg(prim_fac_2)*prim_fac_2
                              fac_8 = dconjg(prim_fac_2)*dconjg(prim_fac_2)*prim_fac_2
                              integrals(k, lm1  , lm2  , lm3  , lm4) = fac_1 * (int_1 - int_2 - int_3 + int_4 - int_5 + int_6 + int_7 - int_8)
                              integrals(k, lm1_p, lm2  , lm3  , lm4) = fac_2 * (int_1 + int_2 - int_3 - int_4 - int_5 - int_6 + int_7 + int_8)
                              integrals(k, lm1  , lm2_p, lm3  , lm4) = fac_3 * (int_1 - int_2 + int_3 - int_4 - int_5 + int_6 - int_7 + int_8)
                              integrals(k, lm1_p, lm2_p, lm3  , lm4) = fac_4 * (int_1 + int_2 + int_3 + int_4 - int_5 - int_6 - int_7 - int_8)
                              integrals(k, lm1  , lm2  , lm3_p, lm4) = fac_5 * (int_1 - int_2 - int_3 + int_4 + int_5 - int_6 - int_7 + int_8)
                              integrals(k, lm1_p, lm2  , lm3_p, lm4) = fac_6 * (int_1 + int_2 - int_3 - int_4 + int_5 + int_6 - int_7 - int_8)
                              integrals(k, lm1  , lm2_p, lm3_p, lm4) = fac_7 * (int_1 - int_2 + int_3 - int_4 + int_5 - int_6 + int_7 - int_8)
                              integrals(k, lm1_p, lm2_p, lm3_p, lm4) = fac_8 * (int_1 + int_2 + int_3 + int_4 + int_5 + int_6 + int_7 + int_8)
                            else
                              lm4_p = l_finish(l4) - m4 + 1
                              !md = -1*l4 + m4
                              md = sph_harm%m_init(l4) + m4
                              msign_8  = m_one**(                md )
                              msign_9  = m_one**( ma           + md )
                              msign_10 = m_one**(      mb      + md )
                              msign_11 = m_one**( ma + mb      + md )
                              msign_12 = m_one**(           mc + md )
                              msign_13 = m_one**( ma      + mc + md )
                              msign_14 = m_one**(      mb + mc + md )
                              msign_15 = m_one**( ma + mb + mc + md )
                              int_9  = msign_8  * integrals_prim(k, lm1  , lm2  , lm3  , lm4_p)
                              int_10 = msign_9  * integrals_prim(k, lm1_p, lm2  , lm3  , lm4_p)
                              int_11 = msign_10 * integrals_prim(k, lm1  , lm2_p, lm3  , lm4_p)
                              int_12 = msign_11 * integrals_prim(k, lm1_p, lm2_p, lm3  , lm4_p)
                              int_13 = msign_12 * integrals_prim(k, lm1  , lm2  , lm3_p, lm4_p)
                              int_14 = msign_13 * integrals_prim(k, lm1_p, lm2  , lm3_p, lm4_p)
                              int_15 = msign_14 * integrals_prim(k, lm1  , lm2_p, lm3_p, lm4_p)
                              int_16 = msign_15 * integrals_prim(k, lm1_p, lm2_p, lm3_p, lm4_p)
                              fac_1  = dconjg(prim_fac_1)*dconjg(prim_fac_1)*prim_fac_1*prim_fac_1
                              fac_2  = dconjg(prim_fac_2)*dconjg(prim_fac_1)*prim_fac_1*prim_fac_1
                              fac_3  = dconjg(prim_fac_1)*dconjg(prim_fac_2)*prim_fac_1*prim_fac_1
                              fac_4  = dconjg(prim_fac_2)*dconjg(prim_fac_2)*prim_fac_1*prim_fac_1
                              fac_5  = dconjg(prim_fac_1)*dconjg(prim_fac_1)*prim_fac_2*prim_fac_1
                              fac_6  = dconjg(prim_fac_2)*dconjg(prim_fac_1)*prim_fac_2*prim_fac_1
                              fac_7  = dconjg(prim_fac_1)*dconjg(prim_fac_2)*prim_fac_2*prim_fac_1
                              fac_8  = dconjg(prim_fac_2)*dconjg(prim_fac_2)*prim_fac_2*prim_fac_1
                              fac_9  = dconjg(prim_fac_1)*dconjg(prim_fac_1)*prim_fac_1*prim_fac_2
                              fac_10 = dconjg(prim_fac_2)*dconjg(prim_fac_1)*prim_fac_1*prim_fac_2
                              fac_11 = dconjg(prim_fac_1)*dconjg(prim_fac_2)*prim_fac_1*prim_fac_2
                              fac_12 = dconjg(prim_fac_2)*dconjg(prim_fac_2)*prim_fac_1*prim_fac_2
                              fac_13 = dconjg(prim_fac_1)*dconjg(prim_fac_1)*prim_fac_2*prim_fac_2
                              fac_14 = dconjg(prim_fac_2)*dconjg(prim_fac_1)*prim_fac_2*prim_fac_2
                              fac_15 = dconjg(prim_fac_1)*dconjg(prim_fac_2)*prim_fac_2*prim_fac_2
                              fac_16 = dconjg(prim_fac_2)*dconjg(prim_fac_2)*prim_fac_2*prim_fac_2

                              integrals(k, lm1  , lm2  , lm3  , lm4  ) = fac_1  * (int_1 - int_2 - int_3 + int_4 - int_5 + int_6 + int_7 - int_8 &
                              &                             - int_9 + int_10+ int_11 - int_12 + int_13 - int_14 - int_15 + int_16)
                              integrals(k, lm1_p, lm2  , lm3  , lm4  ) = fac_2  * (int_1 + int_2 - int_3 - int_4 - int_5 - int_6 + int_7 + int_8 &
                              &                             - int_9 - int_10+ int_11 + int_12 + int_13 + int_14 - int_15 - int_16)
                              integrals(k, lm1  , lm2_p, lm3  , lm4  ) = fac_3  * (int_1 - int_2 + int_3 - int_4 - int_5 + int_6 - int_7 + int_8 &
                              &                             - int_9 + int_10- int_11 + int_12 + int_13 - int_14 + int_15 - int_16)
                              integrals(k, lm1_p, lm2_p, lm3  , lm4  ) = fac_4  * (int_1 + int_2 + int_3 + int_4 - int_5 - int_6 - int_7 - int_8 &
                              &                             - int_9 - int_10- int_11 - int_12 + int_13 + int_14 + int_15 + int_16)
                              integrals(k, lm1  , lm2  , lm3_p, lm4  ) = fac_5  * (int_1 - int_2 - int_3 + int_4 + int_5 - int_6 - int_7 + int_8 &
                              &                             - int_9 + int_10+ int_11 - int_12 - int_13 + int_14 + int_15 - int_16)
                              integrals(k, lm1_p, lm2  , lm3_p, lm4  ) = fac_6  * (int_1 + int_2 - int_3 - int_4 + int_5 + int_6 - int_7 - int_8 &
                              &                             - int_9 - int_10+ int_11 + int_12 - int_13 - int_14 + int_15 + int_16)
                              integrals(k, lm1  , lm2_p, lm3_p, lm4  ) = fac_7  * (int_1 - int_2 + int_3 - int_4 + int_5 - int_6 + int_7 - int_8 &
                              &                             - int_9 + int_10- int_11 + int_12 - int_13 + int_14 - int_15 + int_16)
                              integrals(k, lm1_p, lm2_p, lm3_p, lm4  ) = fac_8  * (int_1 + int_2 + int_3 + int_4 + int_5 + int_6 + int_7 + int_8 &
                              &                             - int_9 - int_10- int_11 - int_12 - int_13 - int_14 - int_15 - int_16)

                              integrals(k, lm1  , lm2  , lm3  , lm4_p) = fac_9  * (int_1 - int_2 - int_3 + int_4 - int_5 + int_6 + int_7 - int_8 &
                              &                             + int_9 - int_10- int_11 + int_12 - int_13 + int_14 + int_15 - int_16)
                              integrals(k, lm1_p, lm2  , lm3  , lm4_p) = fac_10 * (int_1 + int_2 - int_3 - int_4 - int_5 - int_6 + int_7 + int_8 &
                              &                             + int_9 + int_10- int_11 - int_12 - int_13 - int_14 + int_15 + int_16)
                              integrals(k, lm1  , lm2_p, lm3  , lm4_p) = fac_11 * (int_1 - int_2 + int_3 - int_4 - int_5 + int_6 - int_7 + int_8 &
                              &                             + int_9 - int_10+ int_11 - int_12 - int_13 + int_14 - int_15 + int_16)
                              integrals(k, lm1_p, lm2_p, lm3  , lm4_p) = fac_12 * (int_1 + int_2 + int_3 + int_4 - int_5 - int_6 - int_7 - int_8 &
                              &                             + int_9 + int_10+ int_11 + int_12 - int_13 - int_14 - int_15 - int_16)
                              integrals(k, lm1  , lm2  , lm3_p, lm4_p) = fac_13 * (int_1 - int_2 - int_3 + int_4 + int_5 - int_6 - int_7 + int_8 &
                              &                             + int_9 - int_10- int_11 + int_12 + int_13 - int_14 - int_15 + int_16)
                              integrals(k, lm1_p, lm2  , lm3_p, lm4_p) = fac_14 * (int_1 + int_2 - int_3 - int_4 + int_5 + int_6 - int_7 - int_8 &
                              &                             + int_9 + int_10- int_11 - int_12 + int_13 + int_14 - int_15 - int_16)
                              integrals(k, lm1  , lm2_p, lm3_p, lm4_p) = fac_15 * (int_1 - int_2 + int_3 - int_4 + int_5 - int_6 + int_7 - int_8 &
                              &                             + int_9 - int_10+ int_11 - int_12 + int_13 - int_14 + int_15 - int_16)
                              integrals(k, lm1_p, lm2_p, lm3_p, lm4_p) = fac_16 * (int_1 + int_2 + int_3 + int_4 + int_5 + int_6 + int_7 + int_8 &
                              &                             + int_9 + int_10 + int_11 + int_12 + int_13 + int_14 + int_15 + int_16)
                            end if
                          end do
                        end if
                      end do
                    end if
                  end do
 
                end if
              end do
            end do
          end do
        end do
      end do
    end do
    
    deallocate(n_m, l_finish, pos_m_zero)

  end subroutine calc_int_angular_combined

end module DVRIntAng
