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
    if ( (J3 < abs(J1-J2)) .or.  &
    &  (J3 > (J1+J2)) .or.       &
    &  (abs(M1) > J1) .or.       &
    &  (abs(M2) > J2) .or.       &
    &  (abs(M3) > J3)) then
      C = zero
    else
      C = sqrt((J3+J3+one)/fact(nint(J1+J2+J3+one)))
      C = C * sqrt(fact(nint(J1+J2-J3))*fact(nint(J2+J3-J1))* &
      &   fact(nint(J3+J1-J2)))
      C = C * sqrt(fact(nint(J1+M1))*fact(nint(J1-M1))*fact(nint(J2+M2))* &
      &   fact(nint(J2-M2))*fact(nint(J3-M3))*fact(nint(J3+M3)))
      sumk = zero
      do k = 0, 99
        if (J1+J2-J3-K < zero) cycle
        if (J3-J1-M2+K < zero) cycle
        if (J3-J2+M1+K < zero) cycle
        if (J1-M1-K    < zero) cycle
        if (J2+M2-K    < zero) cycle
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
    write(*,*) "Debug:", dim_l

    allocate(integrals(n_mp, dim_l, dim_l, dim_l, dim_l), stat=error)
    call allocerror(error)

    integrals = 0.0d0

  end subroutine allocate_int_ang

  subroutine calc_int_angular(integrals, sph_harm)      
      
    type(sph_harm_t),  intent(in)           :: sph_harm
    real(idp),         intent(inout)        :: integrals(:,:,:,:,:)
 
    integer    :: n_l, n_mp, error, dim_l, k, q, q_init
    integer    :: l1, l2, l3, l4, l13, l24, n_m1, n_m2, m1, m2, m1_init, m2_init

    real(idp)  :: lk, mq, la, lb, lc, ld, ma, mb
    real(idp)  :: pre_fact_ac, pre_fact_bd 
    real(idp)  :: w_symb_ac, w_symb_bd, w_symb_ac_q, w_symb_bd_q 

!   real(idp), allocatable  ::  w_symb_ac_q, w_symb_bd_q

    n_l  = sph_harm%n_l
    n_mp = sph_harm%n_mp

    dim_l = n_l**2

    do k = 1, n_mp
      lk = dfloat(k)
      q_init = -1*(k+1)
      do l1 = 1, n_l
        la = dfloat(l1 - 1)
        do l3 = 1, n_l
          lc = dfloat(l3 - 1)

          pre_fact_ac = sqrt((2.0d0*la)+1.0d0)*sqrt((2.0d0*lc)+1.0d0)

          l13 = min(l1,l3) - 1
          n_m1 = 2*l13 + 1
          m1_init = -1*( l13 + 1)

          w_symb_ac = wigner3j(lk, la, lc, 0.0d0, 0.0d0, 0.0d0)

          do l2 = 1, n_l
            lb = dfloat(l2 - 1)
            do l4 = 1, n_l
              ld = dfloat(l4 - 1)

              pre_fact_bd = sqrt((2.0d0*lb)+1.0d0)*sqrt((2.0d0*ld)+1.0d0)
              l24 = min(lb,ld)
              n_m2 = 2*l24 + 1
              m2_init = -1*( l24 + 1)

              w_symb_bd = wigner3j(lk, lb, ld, 0.0d0, 0.0d0, 0.0d0)

              do m1 = 1, n_m1
                ma = dfloat(m1_init + m1)
                write(*,*) 'la, lc, l13, n_m1', l1-1, l3-1, int(ma)

                do m2 = 1, n_m2
                  ma = dfloat(m2_init + m2)
                  write(*,*) 'lb, lc, l24, n_m2', l2-1, l4-1, int(mb)
                 

                  do q = 1, k
                    mq = dfloat(q_init + q) 
!                   w_symb_ac_q = wigner3j(lk, la, lc, mq, -1.0d0*ma, mc)
!                   w_symb_bd_q = wigner3j(lk, lb, ld, -1.0d0*mq, -1.0d0*mb, md)


                  end do
                end do  ! end loop for l4
              end do  ! end loop for l3
            end do  ! end loop for m2

          end do  ! end loop for l2
        end do  ! end loop for m1
      end do  ! end loop for l1
      end do  ! end loop for k 

  end subroutine calc_int_angular

end module angular_mod 
