module angular_utils

  use DVRData
  use dvr_diag_mod, only : allocerror
  use constants

  
  implicit none

  contains
 
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
  
  subroutine write_int_angular_real(integrals, sph_harm, all_int, file_name_syntx, l_interm)

    type(sph_harm_t),  intent(in)           :: sph_harm
    complex(idp),      intent(in)           :: integrals(:,:,:,:,:)
    logical,           intent(in)           :: all_int
    character(len=32), intent(in)           :: file_name_syntx
    integer,           intent(in)           :: l_interm(:)

    integer    :: n_l, n_mp, error, dim_l, k
    integer    :: l1, l2, l3, l4, m1, m2, m3, m4, lm1, lm2, lm3, lm4
    integer    :: n_m1, n_m2, n_m3, n_m4

    n_l  = sph_harm%n_l
    n_mp = sph_harm%n_mp

    dim_l = n_l**2

    do k = 1, n_mp
      open(11, file=trim(file_name_syntx)//"_l"//trim(int2str(k-1))//".dat",                          &
  &    form="formatted", action="write")
      
      if (all_int) then
        do l1 = 1, n_l 
          n_m1 = sph_harm%n_m(l1)
          do m1 = 1, n_m1
            lm1 = l_interm(l1) + m1
            do l2 = 1, n_l 
              n_m2 = sph_harm%n_m(l2)
              do m2 = 1, n_m2
                lm2 = l_interm(l2) + m2
                do l3 = 1, n_l 
                  n_m3 = sph_harm%n_m(l3)
                  do m3 = 1, n_m3
                    lm3 = l_interm(l3) + m3
                    do l4 = 1, n_l 
                      n_m4 = sph_harm%n_m(l4)
                      do m4 = 1, n_m4
                        lm4 = l_interm(l4) + m4
                        if (abs(integrals(k, lm1, lm2, lm3, lm4)).gt.1e-12)      &
                        &   write(11,'(4I3, ES20.12)') lm1, lm2, lm3, lm4,       &
                        &   real(integrals(k, lm1, lm2, lm3, lm4))
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
          n_m1 = sph_harm%n_m(l1)
          do m1 = 1, n_m1
            lm1 = l_interm(l1) + m1
            do l2 = 1, n_l 
              n_m2 = sph_harm%n_m(l2)
              do m2 = 1, n_m2
                lm2 = l_interm(l2) + m2
                do l3 = 1, n_l 
                  n_m3 = sph_harm%n_m(l3)
                  do m3 = 1, n_m3
                    lm3 = l_interm(l3) + m3
                    do l4 = 1, n_l 
                      n_m4 = sph_harm%n_m(l4)
                      m4 = (-l1 -l2 +l3 +l4) + (m1 +m2 -m3)  ! following the condition m1-m3 = m4-m2 transformed into the index we are using here for m
                      if (0.lt.m4.and.m4.le.n_m4) then
                        lm4 = l_interm(l4) + m4
                         if (abs(integrals(k, lm1, lm2, lm3, lm4)).gt.1e-12)      &
                        &  write(11,'(4I3, ES20.12)') lm1, lm2, lm3, lm4,       &
                        &   real(integrals(k, lm1, lm2, lm3, lm4))
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


  end subroutine write_int_angular_real

  subroutine write_int_angular_imaginary(integrals, sph_harm, all_int, file_name_syntx, l_interm)
      
    type(sph_harm_t),  intent(in)           :: sph_harm
    complex(idp),      intent(in)           :: integrals(:,:,:,:,:)
    logical,           intent(in)           :: all_int
    character(len=32), intent(in)           :: file_name_syntx
    integer,           intent(in)           :: l_interm(:)

    integer    :: n_l, n_mp, error, dim_l, k
    integer    :: l1, l2, l3, l4, m1, m2, m3, m4, lm1, lm2, lm3, lm4
    integer    :: n_m1, n_m2, n_m3, n_m4

    n_l  = sph_harm%n_l
    n_mp = sph_harm%n_mp

    dim_l = n_l**2

    do k = 1, n_mp
      open(11, file=trim(file_name_syntx)//"_l"//trim(int2str(k-1))//".dat",                          &
  &    form="formatted", action="write")
      
      if (all_int) then
        do l1 = 1, n_l 
          n_m1 = sph_harm%n_m(l1)
          do m1 = 1, n_m1
            lm1 = l_interm(l1) + m1
            do l2 = 1, n_l 
              n_m2 = sph_harm%n_m(l2)
              do m2 = 1, n_m2
                lm2 = l_interm(l2) + m2
                do l3 = 1, n_l 
                  n_m3 = sph_harm%n_m(l3)
                  do m3 = 1, n_m3
                    lm3 = l_interm(l3) + m3
                    do l4 = 1, n_l 
                      n_m4 = sph_harm%n_m(l4)
                      do m4 = 1, n_m4
                        lm4 = l_interm(l4) + m4
                        if (abs(integrals(k, lm1, lm2, lm3, lm4)).gt.1e-12)      &
                        &   write(11,'(4I3, 2ES20.12)') lm1, lm2, lm3, lm4,       &
                        &   real(integrals(k, lm1, lm2, lm3, lm4)),   &
                        &   aimag(integrals(k, lm1, lm2, lm3, lm4))   
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
          n_m1 = sph_harm%n_m(l1)
          do m1 = 1, n_m1
            lm1 = l_interm(l1) + m1
            do l2 = 1, n_l 
              n_m2 = sph_harm%n_m(l2)
              do m2 = 1, n_m2
                lm2 = l_interm(l2) + m2
                do l3 = 1, n_l 
                  n_m3 = sph_harm%n_m(l3)
                  do m3 = 1, n_m3
                    lm3 = l_interm(l3) + m3
                    do l4 = 1, n_l 
                      n_m4 = sph_harm%n_m(l4)
                      m4 = (-l1 -l2 +l3 +l4) + (m1 +m2 -m3)  ! following the condition m1-m3 = m4-m2 transformed into the index we are using here for m
                      if (0.lt.m4.and.m4.le.n_m4) then
                        lm4 = l_interm(l4) + m4
                         if (abs(integrals(k, lm1, lm2, lm3, lm4)).gt.1e-12)      &
                        & write(11,'(4I3, ES20.12)') lm1, lm2, lm3, lm4,       &
                        &   real(integrals(k, lm1, lm2, lm3, lm4)),  &
                        &   aimag(integrals(k, lm1, lm2, lm3, lm4))   
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

  end subroutine write_int_angular_imaginary

  integer function get_l(lm, n_l)
  
    integer, intent(in) :: lm, n_l

    integer :: i 

    do i = 1, n_l

      if (lm.gt.((i-1)**2).and.lm.le.(i**2)) then
        get_l = i
        exit
      end if
    end do

  end function get_l

  subroutine initiate_basis_matrix(matrix)

    complex(idp), intent(inout) :: matrix(:,:)

    real(idp), parameter :: val_1 = (1.0d0*sqrt_half)
    real(idp), parameter :: val_2 = (-1.0d0*sqrt_half)

    matrix = zero

!   matrix(1,1) = (0.0d0, val_2)
!   matrix(1,2) = (0.0d0, val_1)
!   matrix(2,1) = (val_1, 0.0d0)
!   matrix(2,2) = (val_1, 0.0d0)
    matrix(1,1) = (1.0d0, 0.0d0)
    matrix(2,2) = (1.0d0, 0.0d0)

    write(*,'(a, 2f15.8)') '(1,1)', real(matrix(1,1)), aimag(matrix(1,1))
    write(*,'(a, 2f15.8)') '(1,2)', real(matrix(1,2)), aimag(matrix(1,2))
    write(*,'(a, 2f15.8)') '(2,1)', real(matrix(2,1)), aimag(matrix(2,1))
    write(*,'(a, 2f15.8)') '(2,2)', real(matrix(2,2)), aimag(matrix(2,2))

  end subroutine initiate_basis_matrix

  subroutine GetTMat(TMat, basis, l)

    complex(idp), allocatable, intent(inout) :: TMat(:,:)
    complex(idp), intent(in) :: basis(2,2)
    integer,  intent(in) :: l

    integer :: n_l, m, mp

    n_l = 2*l-1
    if (allocated(TMat)) deallocate(TMat) 
    allocate(TMat(n_l, n_l))
    TMat = zero

    do m = 1, n_l
      if (m.eq.l) then
        TMat(m,m) = 1.0d0
      else 
        mp = 2*l - m
        TMat(m,1) = basis(1,1)
        TMat(1,m) = basis(1,2)
        TMat(mp,1) = basis(2,1)
        TMat(1,mp) = basis(2,2)
      end if
    end do

  end subroutine GetTMat

end module angular_utils
