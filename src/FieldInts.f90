module FieldIntegrals

  use DVRData
  use FieldData
  use dvr_diag_mod, only: allocerror
  use constants
  use angular_utils
  use util_mod
  use OrbData

  implicit none

  contains

  ! Here is the main subroutine to calculate the integrals for any dipole field
  ! applied externally
  subroutine  CalcFieldInts()

    integer :: iField 
    real(idp) :: tol, start, finish

    call cpu_time(start)

    tol = 1e-12
    
    write(iout,*) '***********'
    write(iout, '(a)') ' Calculating field integrals...'
    write(iout, '(a,i4)') ' Total number of Fields:', nFields
    write(iout, '(a,2x,3a4)') ' Fields:', FieldComp(:)
    ! allocating all the arrays and setting up other parameters related to the 
    ! integrals corresponding to the field
    call SetupFieldInts()

    ! Integrals are calculated first in the primitive FE-DVR basis which required both 
    ! the radial and angular contributions
    call DVRFieldInts()

    ! Integrals are then calculated in the basis of eigenvectors obtained from 
    ! solving Schroedinger equation
    do iField = 1, nFields
      call CombineFieldInts(iField)
      call WriteFieldInts(iField, tol)
    end do

    call cpu_time(finish)

    write(iout,'(X,a,f10.5,X,a)') 'Time taken for calculating field integrals = ', finish-start, 'seconds.'
    write(iout,*) '***********'

  end subroutine CalcFieldInts

  subroutine SetupFieldInts()

    integer :: n_l, dim_lm, error

    n_l = para%l + 1
    dim_lm = n_l**2

    ! allocate the field integrals calculated in the primitive FE-DVR basis
    allocate(PrimFieldInts(dim_lm,dim_lm,para%ng,nFields), stat=error)
    call allocerror(error)

    ! allocate the field integrals calculated in the eigenbasis
    allocate(FieldInts(orb%nSpatialOrbs,orb%nSpatialOrbs,nFields), stat=error)
    call allocerror(error)

  end subroutine SetupFieldInts

  subroutine DVRFieldInts()

    complex(idp), allocatable :: AngPart(:,:,:)

    integer :: n_l, dim_lm, i, error

    ! The pure angular part is required to calculate the 
    ! field integrals. This angular parts are nothing but integrals involving
    ! three spherical harmonics where one of the components are coming 
    ! from the external field which has l=1 and m_l=-1,0,+1. The other two 
    ! components of these spherical harmonics are coming from the orbitals

    n_l = para%l + 1
    dim_lm = n_l**2

    allocate(AngPart(dim_lm,dim_lm,3), stat=error)
    call allocerror(error)

    ! Calculate the pure angular part of the field integrals. 
    ! A combination of these angular parts will be multiplied 
    ! with the radial contribution to give the field integrals
    call GetAngParts(AngPart)

    do i = 1, nFields
      call CalcPrimFieldInts(AngPart, i)
    end do

  end subroutine DVRFieldInts

  subroutine GetAngParts(AngPart)

    complex(idp), allocatable, intent(inout) :: AngPart(:,:,:)

    integer :: la, lb, lk, ma, mb, mk, lmk, lma, lmb, n_ma, n_mb, n_l
    integer :: ma_init, mb_init, la_p, lb_p, lk_p
    real(idp) :: mk_flt, ma_flt, mb_flt, lk_flt, la_flt, lb_flt
    real(idp), allocatable ::  pre_fact(:,:)

    n_l = para%l + 1

    allocate(pre_fact(n_l,n_l))

    lk = 1
    lk_flt = float(lk)
    lk_p = 2*lk - 1

    do lb = 1, n_l
      lb_flt = float(lb-1)
      lb_p = 2*lb - 1
      do la = 1, n_l
        la_flt = float(la-1)
        la_p = 2*la - 1

        pre_fact(la,lb) = wigner3j(la_flt, lb_flt, lk_flt, zero, zero, zero)*sqrt(lk_p*la_p*lb_p/FourPi)

      end do
    end do

    do mk = 1, 3
      lmk = mk
      mk_flt = -2.0d0 + float(mk)
      do lb = 1, n_l
        lb_flt = float(lb-1)
        n_mb = 2*lb-1
        mb_init = -1*lb
        do mb = 1, n_mb 
          mb_flt = dfloat(mb_init + mb)
          lmb = (lb-1)**2 + mb
          do la = 1, n_l
            la_flt = float(la-1)
            n_ma = 2*la-1
            ma_init = -1*la
            do ma = 1, n_ma 
              ma_flt = dfloat(ma_init + ma)
              lma = (la-1)**2 + ma

!             AngPart(lma,lmb, lmk) = wigner3j(la_flt, lk_flt, lb_flt, -1.0d0*ma_flt, mk_flt, mb_flt)
              AngPart(lma,lmb, lmk) = wigner3j(la_flt, lb_flt, lk_flt, -1.0d0*ma_flt, mb_flt, mk_flt)*pre_fact(la,lb)

!             if (abs(AngPart(lma,lmb,lmk)).gt.1e-12) write(78,'(3i4,2f15.8)') lma, lmb, lmk, AngPart(lma,lmb,lmk)
!             write(79,'(6i4,f15.8)') la, lb, lk, int(-1.0d0*ma_flt), int(mb_flt), int(mk_flt), wigner3j(la_flt, lk_flt, lb_flt, -1.0d0*ma_flt, mk_flt, mb_flt)
     
            end do 
          end do
        end do
      end do
    end do

  end subroutine GetAngParts

  subroutine CalcPrimFieldInts(AngPart, iField )

    complex(idp), allocatable, intent(in) :: AngPart(:,:,:)
    integer, intent(in)  :: iField

    complex(idp), pointer :: IntPoint(:,:,:)
    integer :: i, n_l, la, lb, lma, lmb, ma, mb, lmk, lmk_p, lma_p, lmb_p
    real(idp) :: RadPart, pre_fact
    real(idp), parameter :: val_1 = (1.0d0*sqrt_half)
    real(idp), parameter :: val_2 = (-1.0d0*sqrt_half)
    complex(idp) :: int_1, int_2, int_3, int_4, int_5, int_6, int_7, int_8
    complex(idp) :: fact_1, fact_2, fact_3, fact_4, fact_5

    fact_1 = (val_1,0.0d0)
    fact_2 = (val_1, 0.0d0)

    IntPoint =>  PrimFieldInts(:,:,:,iField)

    n_l = para%l + 1

    select case (trim(FieldComp(iField)))
    case('Z')

      lmk = 2
      pre_fact = sqrt(4*pi/3)
      do i= 1, para%ng
        RadPart = pre_fact*grid%r(i)
        do la=1,n_l
          do lb=1,n_l
            do ma = 1,la
              lma = (la-1)**2 + ma
              if (ma.eq.la) then
                do mb = 1, lb
                  lmb = (lb-1)**2 + mb
                  int_1 = AngPart(lma,lmb,lmk)
                  if (mb.eq.lb) then
                    IntPoint(lma,lmb,i) = RadPart*int_1
                  else
                    lmb_p = lb**2 - mb + 1
                    int_2 = AngPart(lma,lmb_p,lmk)
                    IntPoint(lma,lmb,i) = fact_1*RadPart*(int_1 - int_2)
                    IntPoint(lma,lmb_p,i) = fact_2*RadPart*(int_1 + int_2)
                  end if
                end do
              else
                lma_p = la**2 - ma + 1
                do mb = 1, lb
                  lmb = (lb-1)**2 + mb
                  int_1 = AngPart(lma,lmb,lmk)
                  int_2 = AngPart(lma_p,lmb,lmk)
                  if (mb.eq.lb) then
                    IntPoint(lma,lmb,i) = fact_1*RadPart*(int_1 - int_2)
                    IntPoint(lma_p,lmb,i) = fact_2*RadPart*(int_1 + int_2)
                  else
                    fact_3 = fact_1*fact_1
                    fact_4 = fact_2*fact_1
                    fact_5 = fact_2*fact_2
                    lmb_p = lb**2 - mb + 1
                    int_3 = AngPart(lma,lmb_p,lmk)
                    int_4 = AngPart(lma_p,lmb_p,lmk)
                    IntPoint(lma,lmb,i) = fact_3*RadPart*(int_1 - int_2 - int_3 + int_4)
                    IntPoint(lma_p,lmb,i) = fact_4*RadPart*(int_1 + int_2 - int_3 - int_4)
                    IntPoint(lma,lmb_p,i) = fact_4*RadPart*(int_1 - int_2 + int_3 - int_4)
                    IntPoint(lma_p,lmb_p,i) = fact_5*RadPart*(int_1 + int_2 + int_3 + int_4)
!                   write(81, '(2i4, 6f15.8)') lma, lmb, real(int_1), real(int_2), real(int_3), real(int_4), real(IntPoint(lma_p,lmb,i)), real(IntPoint(lma,lmb_p,i))
                  end if
                end do
              end if
            end do
          end do
        end do
      end do

    case('X')

      lmk = 1
      lmk_p = 3
      pre_fact = fact_2*sqrt(4*pi/3)
      do i= 1, para%ng
        RadPart = pre_fact*grid%r(i)
        do la=1,n_l
          do lb=1,n_l
            do ma = 1,la
              lma = (la-1)**2 + ma
              if (ma.eq.la) then
                do mb = 1, lb
                  lmb = (lb-1)**2 + mb
                  int_1 = AngPart(lma,lmb,lmk)
                  int_2 = AngPart(lma,lmb,lmk_p)
                  if (mb.eq.lb) then
                    IntPoint(lma,lmb,i) = RadPart*(int_1 + int_2)
                  else
                    lmb_p = lb**2 - mb + 1
                    int_3 = AngPart(lma,lmb_p,lmk)
                    int_4 = AngPart(lma,lmb_p,lmk_p)
                    IntPoint(lma,lmb,i) = fact_1*RadPart*(int_1 + int_2 - int_3 - int_4)
                    IntPoint(lma,lmb_p,i) = fact_2*RadPart*(int_1 + int_2 + int_3 + int_4)
                  end if
                end do
              else
                lma_p = la**2 - ma + 1
                do mb = 1, lb
                  lmb = (lb-1)**2 + mb
                  int_1 = AngPart(lma,lmb,lmk)
                  int_2 = AngPart(lma_p,lmb,lmk)
                  int_3 = AngPart(lma,lmb,lmk_p)
                  int_4 = AngPart(lma_p,lmb,lmk_p)
                  if (mb.eq.lb) then
                    IntPoint(lma,lmb,i) = fact_1*RadPart*(int_1 - int_2 + int_3 - int_4)
                    IntPoint(lma_p,lmb,i) = fact_2*RadPart*(int_1 + int_2 + int_3 + int_4)
                  else
                    lmb_p = lb**2 - mb + 1
                    int_5 = AngPart(lma,lmb_p,lmk)
                    int_6 = AngPart(lma_p,lmb_p,lmk)
                    int_7 = AngPart(lma,lmb_p,lmk_p)
                    int_8 = AngPart(lma_p,lmb_p,lmk_p)
                    fact_3 = fact_1*fact_1
                    fact_4 = fact_2*fact_1
                    fact_5 = fact_2*fact_2
                    IntPoint(lma,lmb,i) = fact_3*RadPart*(int_1 - int_2 + int_3 - int_4 - int_5 + int_6 - int_7 + int_8)
                    IntPoint(lma_p,lmb,i) = fact_4*RadPart*(int_1 + int_2 + int_3 + int_4 - int_5 - int_6 - int_7 - int_8)
                     
                    IntPoint(lma,lmb_p,i) = fact_4*RadPart*(int_1 - int_2 + int_3 - int_4 + int_5 - int_6 + int_7 - int_8)
                    IntPoint(lma_p,lmb_p,i) = fact_5*RadPart*(int_1 + int_2 + int_3 + int_4 + int_5 + int_6 + int_7 + int_8)
                  end if
                end do
              end if
            end do
          end do
        end do
      end do

    case('Y')

      lmk = 1
      lmk_p = 3
      pre_fact = fact_1*sqrt(4*pi/3)
      do i= 1, para%ng
        RadPart = pre_fact*grid%r(i)
        do la=1,n_l
          do lb=1,n_l
            do ma = 1,la
              lma = (la-1)**2 + ma
              if (ma.eq.la) then
                do mb = 1, lb
                  lmb = (lb-1)**2 + mb
                  int_1 = AngPart(lma,lmb,lmk)
                  int_2 = AngPart(lma,lmb,lmk_p)
                  if (mb.eq.lb) then
                    IntPoint(lma,lmb,i) = RadPart*(int_1 - int_2)
                  else
                    lmb_p = lb**2 - mb + 1
                    int_3 = AngPart(lma,lmb_p,lmk)
                    int_4 = AngPart(lma,lmb_p,lmk_p)
                    IntPoint(lma,lmb,i) = fact_1*RadPart*(int_1 - int_2 - int_3 + int_4)
                    IntPoint(lma,lmb_p,i) = fact_2*RadPart*(int_1 - int_2 + int_3 - int_4)
                  end if
                end do
              else
                lma_p = la**2 - ma + 1
                do mb = 1, lb
                  lmb = (lb-1)**2 + mb
                  int_1 = AngPart(lma,lmb,lmk)
                  int_2 = AngPart(lma_p,lmb,lmk)
                  int_3 = AngPart(lma,lmb,lmk_p)
                  int_4 = AngPart(lma_p,lmb,lmk_p)
                  if (mb.eq.lb) then
                    IntPoint(lma,lmb,i) = fact_1*RadPart*(int_1 - int_2 - int_3 + int_4)
                    IntPoint(lma_p,lmb,i) = fact_2*RadPart*(int_1 + int_2 - int_3 - int_4)
                  else
                    lmb_p = lb**2 - mb + 1
                    int_5 = AngPart(lma,lmb_p,lmk)
                    int_6 = AngPart(lma_p,lmb_p,lmk)
                    int_7 = AngPart(lma,lmb_p,lmk_p)
                    int_8 = AngPart(lma_p,lmb_p,lmk_p)
                    fact_3 = fact_1*fact_1
                    fact_4 = fact_2*fact_1
                    fact_5 = fact_2*fact_2

                    IntPoint(lma,lmb,i) = fact_3*RadPart*(int_1 - int_2 - int_3 + int_4 - int_5 + int_6 + int_7 - int_8)
                    IntPoint(lma_p,lmb,i) = fact_4*RadPart*(int_1 + int_2 - int_3 - int_4 - int_5 - int_6 + int_7 + int_8)
                     
                    IntPoint(lma,lmb_p,i) = fact_4*RadPart*(int_1 - int_2 - int_3 + int_4 + int_5 - int_6 - int_7 + int_8)
                    IntPoint(lma_p,lmb_p,i) = fact_5*RadPart*(int_1 + int_2 - int_3 - int_4 + int_5 + int_6 - int_7 - int_8)
                  end if
                end do
              end if
            end do
          end do
        end do
      end do

    case default 

      call stop_all('CalcPrimFieldInts', 'The perturbation '//trim(FieldComp(iField))//' is not among the available ones')

    end select

  end subroutine CalcPrimFieldInts

  subroutine CombineFieldInts(iField)

    integer, intent(in) :: iField
    real(idp), allocatable ::  inter_int(:)

    integer :: error, na, nb, la, lb, ma, mb, lma, lmb, n_l,  i
    integer :: ind_1, ind_2
    real(idp) ::  int_val

    complex(idp), pointer :: PointInts(:,:), PrimPointInts(:,:,:)

    FieldInts = zero

    !allocate(inter_int(para%ng,orb%n_max), stat=error)
    allocate(inter_int(para%ng), stat=error)
    call allocerror(error)

    n_l = para%l + 1

    PointInts => FieldInts(:,:,iField)
    PrimPointInts => PrimFieldInts(:,:,:,iField)

    if (trim(FieldComp(iField)).eq.'Z') then
      do i= 1, para%ng
        do la=1,n_l
          do ma = 1,2*la-1
            lma = (la-1)**2 + ma
            do lb=1,n_l
              do mb = 1, 2*lb-1
                lmb = (lb-1)**2 + mb
!               write(78,'(5i4,f15.8)') la, ma, lb, mb, i, real(PrimFieldInts(lma,lmb,i,iField))
              end do
            end do
          end do
        end do
      end do
    end if

    do la = 1, n_l
      do ma = 1, 2*la-1
        lma = (la-1)**2 + ma
        do lb = 1, n_l
          do mb = 1, 2*lb-1
            lmb = (lb-1)**2 + mb
            do na = 1, orb%n_max - (la-1)
              do i = 1, para%ng
                inter_int(i) = eigen_vecs(i,na,la)*PrimPointInts(lma,lmb,i)
                if(abs(PrimPointInts(lma,lmb,i)).gt.1e-12) then
!                   write(77,'(5i4, 3f15.8)') la, ma, lb, mb, i, eigen_vecs(i,na,la), real(PrimPointInts(lma,lmb,i)), inter_int(i)
                end if
              end do
              do nb = 1, orb%n_max - (lb-1) 
                int_val = 0.0d0
                do i = 1, para%ng
                  int_val = int_val + inter_int(i)*eigen_vecs(i,nb,lb)
                  if (abs(inter_int(i).gt.1e-12)) then
!                   write(78,'(7i4, 3f15.8)') na, la, ma, nb, lb, mb, i, inter_int(i), eigen_vecs(i, nb, lb), int_val
                  end if
                end do
                ind_1 = SpatialOrbInd(na, la, ma) 
                ind_2 = SpatialOrbInd(nb, lb, mb)
                if (ind_1.eq.0) write(iout,*) 'ind_1 is 0 ', na, la, ma
                if (ind_2.eq.0) write(iout,*) 'ind_2 is 0 ', nb, lb, mb
                PointInts(ind_1, ind_2) = int_val
              end do
            end do  
          end do
        end do
      end do
    end do

  end subroutine CombineFieldInts

  subroutine WriteFieldInts(iField, tol)
    use DVRData, only : para

    integer, intent(in) :: iField
    real(dp), intent(in) :: tol

    integer :: i, j, f_int, norbs, i_n, j_n, nelec
    real(idp) :: core
    complex(idp) :: int_value
    complex(idp), pointer :: PointInts(:,:)
    character(len=32) :: file_int

    PointInts => FieldInts(:,:,iField)

    file_int =  'DIP'//FieldComp(iField)
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

    do i = 1, norbs
      i_n = i + nFrozen
      do j = 1, i
!     do j = 1, norbs
        j_n = j + nFrozen
        int_value = PointInts(i_n,j_n)
        if (abs(int_value).gt.tol) &
        & write(f_int, 1005) real(int_value), i, j, 0, 0
      end do 
    end do

    core = 0.0_dp
    write(f_int, 1005) core, 0, 0, 0, 0

    1005 format(f20.16,x,5i4)

    close(f_int)
  end subroutine WriteFieldInts

end module FieldIntegrals


