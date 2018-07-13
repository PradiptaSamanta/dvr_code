module angular_mod 

  use dvr_spline_mod
  
  implicit none

  public
  
  !! Here is the data_type to containing the parameters for spherical harmonics
  type sph_harm_t
      integer           :: n_l  ! number of L quantum number
      integer           :: n_mp ! orders of multipole included

  end type sph_harm_t

contains
  
  !! @description: Evaluate Wigner 3-j symbol
  !! @param: J1 J1 angular momentum
  !! @param: J2 J2 angular momentum
  !! @param: J3 J3 angular momentum
  !! @param: M1 Projection of J1 angular momentum
  !! @param: M2 Projection of J2 angular momentum
  !! @param: M3 Projection of J3 angular momentum
  real(idp) function wigner3j(J1,J2,J3,M1,M2,M3)

    real(idp), intent(in) :: J1
    real(idp), intent(in) :: J2
    real(idp), intent(in) :: J3
    real(idp), intent(in) :: M1
    real(idp), intent(in) :: M2
    real(idp), intent(in) :: M3

    integer   :: i, k, intexp
    real(idp) :: C, sumk, term
    real(idp), dimension(0:99) :: fact

    ! Compute table of factorials
    fact(0) = one
    do i = 1, 99
      fact(i) = i * fact(i-1)
    end do

    ! Check for invalid input
    if (isfrac(J1+J2+J3) .or. isfrac(J1+M1) .or. isfrac(J2+M2) .or. &
    &   isfrac(J3-M3) .or. isfrac(-J1+J3-M2) .or. isfrac(-J2+J3+M1)) then
      write(*,'(A)') "=================== ERROR ========================"
      write(*,'(A)') "Message: Invalid input of Wigner 3j symbol"
      write(*,'(A)') "=================================================="
      stop
    end if

    ! Compute Clebsch-Gordan coefficient C
    if ( (nint(J3) < abs(nint(J1-J2))) .or.  &
    &  (nint(J3) > nint(J1+J2)) .or.       &
    &  (abs(nint(M1+M2+M3)) > 0) .or.       &
    &  (abs(nint(M1)) > nint(J1)) .or.       &
    &  (abs(nint(M2)) > nint(J2)) .or.       &
    &  (abs(nint(M3)) > nint(J3))) then
      C = zero
    else
      C = sqrt((J3+J3+one)/fact(nint(J1+J2+J3+one)))
      C = C * sqrt(fact(nint(J1+J2-J3))*fact(nint(J2+J3-J1))* &
      &   fact(nint(J3+J1-J2)))
      C = C * sqrt(fact(nint(J1+M1))*fact(nint(J1-M1))*fact(nint(J2+M2))* &
      &   fact(nint(J2-M2))*fact(nint(J3-M3))*fact(nint(J3+M3)))
      sumk = zero
      do k = 0, 99
        if (nint(J1+J2-J3-K) < 0) cycle
        if (nint(J3-J1-M2+K) < 0) cycle
        if (nint(J3-J2+M1+K) < 0) cycle
        if (nint(J1-M1-K)    < 0) cycle
        if (nint(J2+M2-K)    < 0) cycle
        term = fact(nint(J1+J2-J3-k))*fact(nint(J3-J1-M2+k))* &
        &      fact(nint(J3-J2+M1+k))*fact(nint(J1-M1-k))*    &
        &      fact(nint(J2+M2-k))*fact(k)
        if (mod(k,2) == 1) term = -term
        sumk = sumk + one/term
      end do
      C = C * sumk
    end if

    ! calculate 3j symbol from Clebsch-Gordan coefficient
    ! Note: Nagfor treats expressions like (-1)^n with real n as illegal
    !       (because it is) and will throw a floating invalid operation error.
    !       So in order to evaluate the wigner3j symbol we first have to
    !       convert the exponent to an integer expression.
    intexp = nint(J1-J2-M3)
    wigner3j = (-one) ** intexp / sqrt(two * J3 + one) * C

  end function wigner3j
  
  !! @description: Check if argument is fractional
  !! @param: x Argument
  logical function isfrac(x)

    real(idp), intent(in) :: x

    real(idp), parameter  :: eps = 1.0d-8

      if ((abs(x)-int(abs(x))) > eps) then
         isfrac = .true.
      else
         isfrac = .false.
      end if

  end function isfrac

  
  !! @description: Report an allocation error. Print out given the given
  !!               `module_name`, the given `routine_name`, a description of the
  !!               error that occured, and optionally an additional message
  !! @param: i              Error code (`stat` in allocate)
  !! @param: error_message  Additional error message to print
  subroutine allocerror(i, error_message)

    integer,                    intent(in) :: i
    character(len=*), optional, intent(in) :: error_message

    integer :: proc_id, error
    logical :: mpi_is_initialized

    if (i > 0) then
      write(*,'("")')
      write(*,'("================= ALLOCERROR =====================")')
      write(*,'("ERROR: Could not allocate memory")')
      stop
    end if

  end subroutine allocerror

  subroutine allocate_int_ang(integrals, sph_harm)

    type(sph_harm_t),  intent(in)           :: sph_harm
    real(idp), allocatable,  intent(inout)  :: integrals(:,:,:,:,:)
 
    integer    :: n_l, n_mp, error, dim_l

    n_l  = sph_harm%n_l
    n_mp = sph_harm%n_mp

    dim_l = n_l**2

    allocate(integrals(n_mp, dim_l, dim_l, dim_l, dim_l), stat=error)
    call allocerror(error)

    integrals = 0.0d0

  end subroutine allocate_int_ang

  subroutine calc_int_angular(integrals, sph_harm)      
      
    type(sph_harm_t),  intent(in)           :: sph_harm
    real(idp),         intent(inout)        :: integrals(:,:,:,:,:)
 
    integer    :: n_l, n_mp, error, dim_l, k, q, q_init
    integer    :: l1, l2, l3, l4, l13, l24, n_m1, n_m2, n_m3, n_m4
    integer    :: m1, m2, m3, m4, m1_init, m2_init, m3_init, m4_init
    integer    :: lma, lmb, lmc, lmd, m_abq, m_sign

    real(idp)  :: lk, mq, la, lb, lc, ld, ma, mb, mc, md
    real(idp)  :: pre_fact_ac, pre_fact_bd, pre_fact_prod 
    real(idp)  :: w_symb_ac, w_symb_bd, w_symb_ac_q, w_symb_bd_q, w_symb_abcd
    real(idp)  :: int_value

!   real(idp), allocatable  ::  w_symb_ac_q, w_symb_bd_q

    n_l  = sph_harm%n_l
    n_mp = sph_harm%n_mp

    dim_l = n_l**2

    do k = 1, n_mp
      lk = dfloat(k-1)
      q_init = -1*k
      do l1 = 1, n_l
        la = dfloat(l1 - 1)
        n_m1 = 2*l1 - 1
        m1_init = -1*l1
        do l3 = 1, n_l
          lc = dfloat(l3 - 1)
          n_m3 = 2*l3 - 1
          m3_init = -1*l3

          pre_fact_ac = sqrt((2.0d0*la)+1.0d0)*sqrt((2.0d0*lc)+1.0d0)

          w_symb_ac = wigner3j(lk, la, lc, 0.0d0, 0.0d0, 0.0d0)

          do l2 = 1, n_l
            lb = dfloat(l2 - 1)
            n_m2 = 2*l2 - 1
            m2_init = -1*l2
            do l4 = 1, n_l
              ld = dfloat(l4 - 1)
              n_m4 = 2*l4 - 1
              m4_init = -1*l4

              pre_fact_bd = sqrt((2.0d0*lb)+1.0d0)*sqrt((2.0d0*ld)+1.0d0)

              w_symb_bd = wigner3j(lk, lb, ld, 0.0d0, 0.0d0, 0.0d0)

              pre_fact_prod = pre_fact_ac * pre_fact_bd
              w_symb_abcd = w_symb_ac * w_symb_bd

              do m1 = 1, n_m1
                ma = dfloat(m1_init + m1)
                lma = (l1-1)**2 + m1
                !lma = (l1-1)**2 + int(la) + int(ma) + 1
                do m3 = 1, n_m3
                  mc = dfloat(m3_init + m3)
                  lmc = (l3-1)**2 + m3
                  !lmc = (l3-1)**2 + int(lc) + int(mc) + 1

                  do m2 = 1, n_m2
                    mb = dfloat(m2_init + m2)
                    lmb = (l2-1)**2 + m2
                    !lmb = (l2-1)**2 + int(lb) + int(mb) + 1
                    do m4 = 1, n_m4
                      md = dfloat(m4_init + m4)
                      lmd = (l4-1)**2 + m4
                      !lmd = (l4-1)**2 + int(ld) + int(md) + 1

                      int_value = 0.0d0
                      do q = 1, 2*k-1
!                       if (lma.eq.1.and.lmb.eq.1.and.lmc.eq.2.and.lmd.eq.4) then
                        mq = dfloat(q_init + q) 

                        m_abq = int(ma + mb + mq)
                        m_sign = (-1)**m_abq
!                   write(*,*) 'k, lk, q, mq, m_abq:',  k, int(lk), q, int(mq), m_sign

                        w_symb_ac_q =                                          &
                        &     wigner3j(lk, la, lc, mq, -1.0d0*ma, mc)

                        w_symb_bd_q =                                          &
                        &     wigner3j(lk, lb, ld, -1.0d0*mq, -1.0d0*mb,md)
!                   write(*,*) lk, lb, ld, w_symb_bd_q

                        int_value = int_value +                                &
                     &   m_sign !* pre_fact_prod !* w_symb_abcd * w_symb_ac_q*   &
                     !&   w_symb_bd_q 
                        integrals(k, lma, lmb, lmc, lmd) = int_value
!                        write(*,*) q, int_value, w_symb_ac_q, w_symb_bd_q

!                     end if
                      end do ! end loop for q

                    end do ! end loop for m4
                  end do ! end loop for m2   
                end do  ! end loop for m3
              end do  ! end loop for m1

            end do  ! end loop for l4
          end do  ! end loop for l2

        end do  ! end loop for l3
      end do  ! end loop for l1
    end do  ! end loop for k 

  end subroutine calc_int_angular

  subroutine calc_int_angular_main(integrals, sph_harm)      
      
    type(sph_harm_t),  intent(in)           :: sph_harm
    real(idp),         intent(inout)        :: integrals(:,:,:,:,:)
 
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

    dim_l = n_l**2
    dim_mp = n_mp**2

    allocate(pre_fact_ac(n_l,n_l), pre_fact_bd(n_l,n_l))
    allocate(w_symb_ac(n_l,n_l), w_symb_bd(n_l,n_l))
    allocate(w_symb_ac_q(dim_mp, dim_l, dim_l), w_symb_bd_q(dim_mp, dim_l, dim_l))

    do k = 1, n_mp
      lk = dfloat(k-1)
      q_init = -1*k

      do l1 = 1, n_l
        la = dfloat(l1 - 1)
        n_m1 = 2*l1 - 1
        m1_init = -1*l1
        do l3 = 1, n_l
          lc = dfloat(l3 - 1)
          n_m3 = 2*l3 - 1
          m3_init = -1*l3

          pre_fact_ac(l1,l3) = sqrt((2.0d0*la)+1.0d0)*sqrt((2.0d0*lc)+1.0d0)

          w_symb_ac(l1,l3) = wigner3j(lk, la, lc, 0.0d0, 0.0d0, 0.0d0)

          do m1 = 1, n_m1
            ma = dfloat(m1_init + m1)
            lma = (l1-1)**2 + m1
            !lma = (l1-1)**2 + int(la) + int(ma) + 1
            do m3 = 1, n_m3
              mc = dfloat(m3_init + m3)
              lmc = (l3-1)**2 + m3
              !lmc = (l3-1)**2 + int(lc) + int(mc) + 1

              do q = 1, 2*k-1
!               if (lma.eq.1.and.lmb.eq.1.and.lmc.eq.2.and.lmd.eq.4) then
                mq = dfloat(q_init + q)

                lmk = (lk-1)**2 + q
                m_abq = int(ma + mb + mq)
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
        n_m2 = 2*l2 - 1
        m2_init = -1*l2
        do l4 = 1, n_l
          ld = dfloat(l4 - 1)
          n_m4 = 2*l4 - 1
          m4_init = -1*l4

          pre_fact_bd(l2,l4) = sqrt((2.0d0*lb)+1.0d0)*sqrt((2.0d0*ld)+1.0d0)

          w_symb_bd(l2, l4) = wigner3j(lk, lb, ld, 0.0d0, 0.0d0, 0.0d0)

!         pre_fact_prod = pre_fact_ac * pre_fact_bd
!         w_symb_abcd = w_symb_ac * w_symb_bd

          do m2 = 1, n_m2
            mb = dfloat(m2_init + m2)
            lmb = (l2-1)**2 + m2
            !lmb = (l2-1)**2 + int(lb) + int(mb) + 1
            do m4 = 1, n_m4
              md = dfloat(m4_init + m4)
              lmd = (l4-1)**2 + m4
              !lmd = (l4-1)**2 + int(ld) + int(md) + 1

              int_value = 0.0d0
              do q = 1, 2*k-1
!               if (lma.eq.1.and.lmb.eq.1.and.lmc.eq.2.and.lmd.eq.4) then
                mq = dfloat(q_init + q) 
                lmk = (lk-1)**2 + q

                m_abq = int(ma + mb + mq)
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
        n_m1 = 2*l1 - 1
        m1_init = -1*l1
        do l2 = 1, n_l 
          n_m2 = 2*l2 - 1
          m2_init = -1*l2
          do l3 = 1, n_l 
            n_m3 = 2*l3 - 1
            do l4 = 1, n_l 
              n_m4 = 2*l4 - 1
              pre_fact_prod = pre_fact_ac(l1, l3) * pre_fact_bd(l2, l4)
              w_symb_abcd = w_symb_ac(l1,l3)*w_symb_bd(l2,l4)
              do m1 = 1, n_m1
                ma = m1_init + m1
                lma = (l1-1)**2 + m1
                do m2 = 1, n_m2
                  mb = m2_init + m2
                  lmb = (l2-1)**2 + m2
                  do m3 = 1, n_m3
                    lmc = (l3-1)**2 + m3
                    m4 = (-l1 -l2 +l3 +l4) + (m1 +m2 -m3)  ! following the condition m1-m3 = m4-m2 transformed into the index we are using here for m
                    if (0.lt.m4.and.m4.le.n_m4) then
                      lmd = (l4-1)**2 + m4
                      int_value = 0.0d0
                      do q = 1, 2*k-1
!               if (lma.eq.1.and.lmb.eq.1.and.lmc.eq.2.and.lmd.eq.4) then
                        mq = dfloat(q_init + q) 
                        lmk = (lk-1)**2 + q
   
                        m_abq = int(ma + mb + mq)
                        m_sign = (-1)**m_abq
                        w_symb_abcd_q = w_symb_ac_q(lmk, lma, lmc)        &
                        &              *w_symb_bd_q(lmk, lmb, lmd)
                        int_value = int_value +                                &
                  &     m_sign * pre_fact_prod  * w_symb_abcd * w_symb_abcd_q
                        integrals(k, lma, lmb, lmc, lmd) = int_value
                      end do
                    end if
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

  subroutine write_int_angular(integrals, sph_harm, all_int)
      
    type(sph_harm_t),  intent(in)           :: sph_harm
    real(idp),         intent(in)           :: integrals(:,:,:,:,:)
    logical,           intent(in)           :: all_int

    integer    :: n_l, n_mp, error, dim_l, k
    integer    :: l1, l2, l3, l4, m1, m2, m3, m4, lm1, lm2, lm3, lm4
    integer    :: n_m1, n_m2, n_m3, n_m4

    n_l  = sph_harm%n_l
    n_mp = sph_harm%n_mp

    dim_l = n_l**2

    do k = 1, n_mp
      open(11, file="ang_element_l"//trim(int2str(k-1))//".dat",                          &
  &    form="formatted", action="write")
      
      if (all_int) then
        do l1 = 1, n_l 
          n_m1 = 2*l1 - 1
          do m1 = 1, n_m1
            lm1 = (l1-1)**2 + m1
            do l2 = 1, n_l 
              n_m2 = 2*l2 - 1
              do m2 = 1, n_m2
                lm2 = (l2-1)**2 + m2
                do l3 = 1, n_l 
                  n_m3 = 2*l3 - 1
                  do m3 = 1, n_m3
                    lm3 = (l3-1)**2 + m3
                    do l4 = 1, n_l 
                      n_m4 = 2*l4 - 1
                      do m4 = 1, n_m4
                        lm4 = (l4-1)**2 + m4
                        if (abs(integrals(k, lm1, lm2, lm3, lm4)).gt.1e-12)      &
                        &   write(11,'(4I3, ES25.17)') lm1, lm2, lm3, lm4,       &
                        &   integrals(k, lm1, lm2, lm3, lm4)
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      else
        do l1 = 1, n_l 
          n_m1 = 2*l1 - 1
          do m1 = 1, n_m1
            lm1 = (l1-1)**2 + m1
            do l2 = 1, n_l 
              n_m2 = 2*l2 - 1
              do m2 = 1, n_m2
                lm2 = (l2-1)**2 + m2
                do l3 = 1, n_l 
                  n_m3 = 2*l3 - 1
                  do m3 = 1, n_m3
                    lm3 = (l3-1)**2 + m3
                    do l4 = 1, n_l 
                      n_m4 = 2*l4 - 1
                      m4 = (-l1 -l2 +l3 +l4) + (m1 +m2 -m3)  ! following the condition m1-m3 = m4-m2 transformed into the index we are using here for m
                      if (0.lt.m4.and.m4.le.n_m4) then
                        lm4 = (l4-1)**2 + m4
                         if (abs(integrals(k, lm1, lm2, lm3, lm4)).gt.1e-12)      &
                        &  write(11,'(4I3, ES25.17)') lm1, lm2, lm3, lm4,       &
                        &   integrals(k, lm1, lm2, lm3, lm4)
                      end if
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      endif
      close(11)
    end do


  end subroutine write_int_angular

  character(len=converted_l) function int2str(i, format)                         
                                                                                 
    integer,                    intent(in) :: i                                  
    character(len=*), optional, intent(in) :: format                             
                                                                                 
    if (present(format)) then                                                    
      write(int2str, format) i                                                   
    else                                                                         
      write(int2str, '(I25)') i                                                  
    end if                                                                       
    int2str = adjustl(int2str)                                                   
                                                                                 
  end function int2str
  
end module angular_mod 
