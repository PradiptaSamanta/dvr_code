module FieldIntegrals

  use DVRData
  use FieldData
  use dvr_diag_mod, only: allocerror
  use constants
  use angular_utils
  use util_mod
  use OrbData
  use OrbInts, only : SetUpEigVec, SetUpEigVecBreak

  implicit none

  contains

  ! Here is the main subroutine to calculate the integrals for any dipole field
  ! applied externally
  subroutine  CalcFieldInts()

    integer :: iField 
    real(idp) :: tol, start, finish
    real(idp), allocatable :: EigVecs(:,:,:)
    integer :: n_l, l
    integer, allocatable :: n_m(:), ml_interm(:), m_init(:)

    call cpu_time(start)

    tol = 1e-12
    
    n_l = para%l + 1
    allocate(n_m(n_l), ml_interm(n_l), m_init(n_l))

    n_m = sph_harm%n_m
    m_init = sph_harm%m_init

    ml_interm(1) = 0
    do l = 2, n_l
      ml_interm(l) = sum(n_m(1:l-1))
    end do

    write(iout,*) '***********'
    write(iout, '(a)') ' Calculating field integrals...'
    write(iout, '(a,i4)') ' Total number of Fields:', nFields
    write(iout, '(a,2x,3a4)') ' Fields:', FieldComp(:)
    ! allocating all the arrays and setting up other parameters related to the 
    ! integrals corresponding to the field
    call SetupFieldInts()

    ! Integrals are calculated first in the primitive FE-DVR basis which required both 
    ! the radial and angular contributions
    call DVRFieldInts(n_m, ml_interm, m_init)

    ! Integrals are then calculated in the basis of eigenvectors obtained from 
    ! solving Schroedinger equation
    ! First the eigenvectors are made block diagonal when the total radial grid are split
    if (para%split_grid) then
      call SetUpEigVec(eigen_vecs, EigVecs)
    elseif (break_inn) then
        call SetUpEigVecBreak(eigen_vecs, EigVecs) ! ** Not fully checked
    end if

    do iField = 1, nFields
!     call CombineFieldInts(iField)
      if (para%split_grid.or.break_inn) then
        call ConvertFieldInts(iField, EigVecs, n_m, ml_interm)
      else 
        call ConvertFieldInts(iField, eigen_vecs, n_m, ml_interm)
      end if

      call WriteFieldInts(iField, tol)
    end do

    call cpu_time(finish)

    write(iout,'(X,a,f10.5,X,a)') 'Time taken for calculating field integrals = ', finish-start, 'seconds.'
    write(iout,*) '***********'

    deallocate(n_m, ml_interm, m_init)

  end subroutine CalcFieldInts

  subroutine SetupFieldInts()

    integer :: n_l, dim_lm, error

    n_l = para%l + 1
    dim_lm = para%dim_l

    ! allocate the field integrals calculated in the primitive FE-DVR basis
    allocate(PrimFieldInts(dim_lm,dim_lm,para%ng,nFields), stat=error)
    call allocerror(error)

    ! allocate the field integrals calculated in the eigenbasis
    allocate(FieldInts(orb%nSpatialOrbs,orb%nSpatialOrbs,nFields), stat=error)
    call allocerror(error)

  end subroutine SetupFieldInts

  subroutine DVRFieldInts(n_m, ml_interm, m_init)

    complex(idp), allocatable :: AngPart(:,:,:)
    integer, allocatable, intent(in) ::  n_m(:), ml_interm(:), m_init(:)

    integer :: n_l, dim_lm, i, error

    ! The pure angular part is required to calculate the 
    ! field integrals. This angular parts are nothing but integrals involving
    ! three spherical harmonics where one of the components are coming 
    ! from the external field which has l=1 and m_l=-1,0,+1. The other two 
    ! components of these spherical harmonics are coming from the orbitals

    n_l = para%l + 1
    dim_lm = para%dim_l

    allocate(AngPart(dim_lm,dim_lm,3), stat=error)
    call allocerror(error)

    ! Calculate the pure angular part of the field integrals. 
    ! A combination of these angular parts will be multiplied 
    ! with the radial contribution to give the field integrals
    call GetAngParts(AngPart, n_m, ml_interm, m_init)

    do i = 1, nFields
!     call CalcAngRealBasis(AngPart,i)
      call CalcPrimFieldInts(AngPart, i, ml_interm, m_init)
    end do

  end subroutine DVRFieldInts

  subroutine GetAngParts(AngPart, n_m, ml_interm, m_init)

    complex(idp), allocatable, intent(inout) :: AngPart(:,:,:)
    integer, allocatable, intent(in) ::  n_m(:), ml_interm(:), m_init(:)

    integer :: la, lb, lk, ma, mb, mk, lmk, lma, lmb, n_l, lmk_p
    integer :: ma_init, mb_init, la_p, lb_p, lk_p
    real(idp) :: mk_flt, ma_flt, mb_flt, lk_flt, la_flt, lb_flt, part_1, part_2
    real(idp), allocatable ::  pre_fact(:,:)
    logical :: real_func
    real(idp), parameter :: val_1 = (1.0d0*sqrt_half)
    real(idp), parameter :: val_2 = (1.0d0*sqrt_half)
    complex(idp) :: fact_1, fact_2

    fact_1 = (val_1,0.0d0)
    !fact_2 = (val_1, 0.0d0)
    fact_2 = ( 0.0d0, val_1)

!   real_func = .true.
    real_func = .false.

    n_l = para%l + 1

    allocate(pre_fact(n_l,n_l))

    lk = 1
    lk_flt = float(lk)
    lk_p = 2*lk + 1

    do lb = 1, n_l
      lb_flt = float(lb-1)
      lb_p = 2*lb - 1
      do la = 1, n_l
        la_flt = float(la-1)
        la_p = 2*la - 1

        pre_fact(la,lb) = wigner3j(la_flt, lb_flt, lk_flt, zero, zero, zero)*sqrt(lk_p*la_p*lb_p/FourPi)

      end do
    end do

    if (real_func) then

      do mk = 1,2
        lmk = mk
        mk_flt = -2.0d0 + float(mk)
        do lb = 1, n_l
          lb_flt = float(lb-1)
          !n_mb = 2*lb-1
          mb_init = m_init(lb)
          !do mb = 1, n_mb 
          do mb = 1, n_m(lb)
            mb_flt = dfloat(mb_init + mb)
            lmb = ml_interm(lb) + mb
            do la = 1, n_l
              la_flt = float(la-1)
              ma_init = m_init(la)
              do ma = 1, n_m(la)
                ma_flt = dfloat(ma_init + ma)
                lma = ml_interm(la) + ma
                if (mk.eq.2) then
!                 AngPart(lma,lmb, lmk) = wigner3j(la_flt, lk_flt, lb_flt, -1.0d0*ma_flt, mk_flt, mb_flt)*pre_fact(la,lb)
                  AngPart(lma,lmb, lmk) = m_one**(ma_flt) * wigner3j(la_flt, lb_flt, lk_flt, -1.0d0*ma_flt, mb_flt, mk_flt)*pre_fact(la,lb)
                else
                  lmk_p = 3
                  part_1 = wigner3j(la_flt, lb_flt, lk_flt, -1.0d0*ma_flt, mb_flt, mk_flt)
                  part_2 = wigner3j(la_flt, lb_flt, lk_flt, -1.0d0*ma_flt, mb_flt, -1.0d0*mk_flt)
                  AngPart(lma,lmb, lmk  ) = m_one**(ma_flt) *pre_fact(la,lb)*fact_2*(part_1 + part_2)
                  AngPart(lma,lmb, lmk_p) = m_one**(ma_flt) *pre_fact(la,lb)*fact_1*(part_1 - part_2)
                end if
              end do
            end do
          end do
        end do
      end do

    else

      do mk = 1, 3
        lmk = mk
        mk_flt = -2.0d0 + float(mk)
        do lb = 1, n_l
          lb_flt = float(lb-1)
          mb_init = m_init(lb)
          do mb = 1, n_m(lb)
            mb_flt = dfloat(mb_init + mb)
            lmb = ml_interm(lb) + mb
            do la = 1, n_l
              la_flt = float(la-1)
              ma_init = m_init(la)
              do ma = 1, n_m(la)
                ma_flt = dfloat(ma_init + ma)
                lma = ml_interm(la) + ma
     
                AngPart(lma,lmb, lmk) = m_one**(ma_flt) * wigner3j(la_flt, lk_flt, lb_flt, -1.0d0*ma_flt, 1.0d0*mk_flt, mb_flt)*pre_fact(la,lb)
!               AngPart(lma,lmb, lmk) = m_one**ma_flt * wigner3j(la_flt, lb_flt, lk_flt, -1.0d0*ma_flt, mb_flt, mk_flt)*pre_fact(la,lb)
     
                if (abs(AngPart(lma,lmb,lmk)).gt.1e-12) write(78,'(3i4,2f15.8)') lma, lmb, lmk, AngPart(lma,lmb,lmk)
!               write(79,'(6i4,f15.8)') la, lb, lk, int(-1.0d0*ma_flt), int(mb_flt), int(mk_flt), wigner3j(la_flt, lk_flt, lb_flt, -1.0d0*ma_flt, mk_flt, mb_flt)
       
              end do 
            end do
          end do
        end do
      end do

    end if

    deallocate(pre_fact)

  end subroutine GetAngParts

  subroutine CalcAngRealBasis(AngPart, iField)

    complex(idp), allocatable, intent(in) :: AngPart(:,:,:)
    integer, intent(in)  :: iField

    complex(idp), allocatable :: IntPoint(:,:)
    integer :: n_l, la, lb, lma, lmb, ma, mb, m1, m2, lmk, lmk_p, lma_p, lmb_p, n_l2
    complex(idp) :: RadPart, pre_fact
    real(idp), parameter :: val_1 = (1.0d0*sqrt_half)
    real(idp), parameter :: val_2 = (1.0d0*sqrt_half)
    complex(idp) :: int_1, int_2, int_3, int_4, int_5, int_6, int_7, int_8
    complex(idp) :: prim_fac_1, prim_fac_2, fac_1, fac_2, fac_3, fac_4
    real(idp) :: msign_1, msign_2, msign_3

    prim_fac_1= dcmplx(0.0d0, val_1)
    prim_fac_2= dcmplx(val_1, 0.0d0)

    n_l = para%l + 1

    n_l2 = n_l*n_l

    allocate(IntPoint(n_l2,n_l2))
    IntPoint = czero

    select case (trim(FieldComp(iField)))
    case('Z')

      lmk = 2
      pre_fact = sqrt(4*pi/3)
      RadPart = pre_fact
      do la=1,n_l
        do lb=1,n_l
          do ma = 1,la
            lma = (la-1)**2 + ma
            if (ma.eq.la) then
              do mb = 1, lb
                lmb = (lb-1)**2 + mb
                int_1 = AngPart(lma,lmb,lmk)
                if (mb.eq.lb) then
                  IntPoint(lma,lmb) = RadPart*int_1
                else
                  lmb_p = lb**2 - mb + 1
                  m2 = -1*lb + mb
                  msign_1 = m_one**m2
                  int_2 = msign_1 * AngPart(lma,lmb_p,lmk)
                  fac_1 = prim_fac_1
                  fac_2 = prim_fac_2
                  IntPoint(lma,lmb  ) = fac_1*RadPart*(int_1 - int_2)
                  IntPoint(lma,lmb_p) = fac_2*RadPart*(int_1 + int_2)
                end if
              end do
            else
              lma_p = la**2 - ma + 1
              m1 = -1*la + ma
              msign_1 = m_one**m1
              do mb = 1, lb
                lmb = (lb-1)**2 + mb
                int_1 =           AngPart(lma,lmb,lmk)
                int_2 = msign_1 * AngPart(lma_p,lmb,lmk)
                if (mb.eq.lb) then
                  fac_1 = dconjg(prim_fac_1)
                  fac_2 = dconjg(prim_fac_2)
                  IntPoint(lma  ,lmb) = fac_1*RadPart*(int_1 - int_2)
                  IntPoint(lma_p,lmb) = fac_2*RadPart*(int_1 + int_2)
                else
                  lmb_p = lb**2 - mb + 1
                  m2 = -1*lb + mb
                  msign_2 = m_one**(      m2 )
                  msign_3 = m_one**( m1 + m2 )
                  fac_1 = dconjg(prim_fac_1)*prim_fac_1
                  fac_2 = dconjg(prim_fac_2)*prim_fac_1
                  fac_3 = dconjg(prim_fac_1)*prim_fac_2
                  fac_4 = dconjg(prim_fac_2)*prim_fac_2
                  int_3 = msign_2 * AngPart(lma  ,lmb_p,lmk)
                  int_4 = msign_3 * AngPart(lma_p,lmb_p,lmk)
                  IntPoint(lma  ,lmb  ) = fac_1*RadPart*(int_1 - int_2 - int_3 + int_4)
                  IntPoint(lma_p,lmb  ) = fac_2*RadPart*(int_1 + int_2 - int_3 - int_4)
                  IntPoint(lma  ,lmb_p) = fac_3*RadPart*(int_1 - int_2 + int_3 - int_4)
                  IntPoint(lma_p,lmb_p) = fac_4*RadPart*(int_1 + int_2 + int_3 + int_4)
!                 write(81, '(2i4, 6f15.8)') lma, lmb, real(int_1), real(int_2), real(int_3), real(int_4), real(IntPoint(lma_p,lmb,i)), real(IntPoint(lma,lmb_p,i))
                end if
              end do
            end if
          end do
        end do
      end do

    case('Y')

      lmk = 1
      lmk_p = 3
      pre_fact = prim_fac_1*sqrt(4*pi/3)
      RadPart = pre_fact
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
                  IntPoint(lma,lmb) = RadPart*(int_1 + int_2)
                else
                  lmb_p = lb**2 - mb + 1
                  m2 = -1*lb + mb
                  msign_1 = m_one**m2
                  fac_1 = prim_fac_1
                  fac_2 = prim_fac_2
                  int_3 = msign_1 * AngPart(lma,lmb_p,lmk  )
                  int_4 = msign_1 * AngPart(lma,lmb_p,lmk_p)
                  IntPoint(lma,lmb  ) = fac_1*RadPart*(int_1 + int_2 - int_3 - int_4)
                  IntPoint(lma,lmb_p) = fac_2*RadPart*(int_1 + int_2 + int_3 + int_4)
                end if
              end do
            else
              lma_p = la**2 - ma + 1
              m1 = -1*la + ma
              msign_1 = m_one**m1
              do mb = 1, lb
                lmb = (lb-1)**2 + mb
                int_1 =           AngPart(lma  ,lmb,lmk  )
                int_2 = msign_1 * AngPart(lma_p,lmb,lmk  )
                int_3 =           AngPart(lma  ,lmb,lmk_p)
                int_4 = msign_1 * AngPart(lma_p,lmb,lmk_p)
                if (mb.eq.lb) then
                  fac_1 = dconjg(prim_fac_1)
                  fac_2 = dconjg(prim_fac_2)
                  IntPoint(lma  ,lmb) = fac_1*RadPart*(int_1 - int_2 + int_3 - int_4)
                  IntPoint(lma_p,lmb) = fac_2*RadPart*(int_1 + int_2 + int_3 + int_4)
                else
                  lmb_p = lb**2 - mb + 1
                  m2 = -1*lb + mb
                  msign_2 = m_one**(      m2 )
                  msign_3 = m_one**( m1 + m2 )
                  int_5 = msign_2 * AngPart(lma  ,lmb_p,lmk  )
                  int_6 = msign_3 * AngPart(lma_p,lmb_p,lmk  )
                  int_7 = msign_2 * AngPart(lma  ,lmb_p,lmk_p)
                  int_8 = msign_3 * AngPart(lma_p,lmb_p,lmk_p)
                  fac_1 = dconjg(prim_fac_1)*prim_fac_1
                  fac_2 = dconjg(prim_fac_2)*prim_fac_1
                  fac_3 = dconjg(prim_fac_1)*prim_fac_2
                  fac_4 = dconjg(prim_fac_2)*prim_fac_2
                  IntPoint(lma  ,lmb  ) = fac_1*RadPart*(int_1 - int_2 + int_3 - int_4 - int_5 + int_6 - int_7 + int_8)
                  IntPoint(lma_p,lmb  ) = fac_2*RadPart*(int_1 + int_2 + int_3 + int_4 - int_5 - int_6 - int_7 - int_8)
                  IntPoint(lma  ,lmb_p) = fac_3*RadPart*(int_1 - int_2 + int_3 - int_4 + int_5 - int_6 + int_7 - int_8)
                  IntPoint(lma_p,lmb_p) = fac_4*RadPart*(int_1 + int_2 + int_3 + int_4 + int_5 + int_6 + int_7 + int_8)
!                 if (abs(IntPoint(lma_p,lmb_p)).gt.0.0d0) then
!                     write(95,'(8f15.8)') real(int_1), real(int_2), real(int_3), real(int_4), real(int_5), real(int_6), real(int_7), real(int_8)
!                 end if
                end if
              end do
            end if
          end do
        end do
      end do

    case('X')

      lmk = 1
      lmk_p = 3
      pre_fact = prim_fac_2*sqrt(4*pi/3)
      RadPart = pre_fact
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
                  IntPoint(lma,lmb) = RadPart*(int_1 - int_2)
                else
                  lmb_p = lb**2 - mb + 1
                  m2 = -1*lb + mb
                  msign_1 = m_one**m2
                  fac_1 = prim_fac_1
                  fac_2 = prim_fac_2
                  int_3 = msign_1 * AngPart(lma,lmb_p,lmk)
                  int_4 = msign_1 * AngPart(lma,lmb_p,lmk_p)
                  IntPoint(lma,lmb  ) = fac_1*RadPart*(int_1 - int_2 - int_3 + int_4)
                  IntPoint(lma,lmb_p) = fac_2*RadPart*(int_1 - int_2 + int_3 - int_4)
                end if
              end do
            else
              lma_p = la**2 - ma + 1
              m1 = -1*la + ma
              msign_1 = m_one**m1
              do mb = 1, lb
                lmb = (lb-1)**2 + mb
                int_1 =           AngPart(lma  ,lmb,lmk  )
                int_2 = msign_1 * AngPart(lma_p,lmb,lmk  )
                int_3 =           AngPart(lma  ,lmb,lmk_p)
                int_4 = msign_1 * AngPart(lma_p,lmb,lmk_p)
                if (mb.eq.lb) then
                  fac_1 = dconjg(prim_fac_1)
                  fac_2 = dconjg(prim_fac_2)
                  IntPoint(lma  ,lmb) = fac_1*RadPart*(int_1 - int_2 - int_3 + int_4)
                  IntPoint(lma_p,lmb) = fac_2*RadPart*(int_1 + int_2 - int_3 - int_4)
                else
                  lmb_p = lb**2 - mb + 1
                  m2 = -1*lb + mb
                  msign_2 = m_one**(      m2 )
                  msign_3 = m_one**( m1 + m2 )
                  int_5 = msign_2 * AngPart(lma  ,lmb_p,lmk  )
                  int_6 = msign_3 * AngPart(lma_p,lmb_p,lmk  )
                  int_7 = msign_2 * AngPart(lma  ,lmb_p,lmk_p)
                  int_8 = msign_3 * AngPart(lma_p,lmb_p,lmk_p)

                  fac_1 = dconjg(prim_fac_1)*prim_fac_1
                  fac_2 = dconjg(prim_fac_2)*prim_fac_1
                  fac_3 = dconjg(prim_fac_1)*prim_fac_2
                  fac_4 = dconjg(prim_fac_2)*prim_fac_2
                  IntPoint(lma  ,lmb  ) = fac_1*RadPart*(int_1 - int_2 - int_3 + int_4 - int_5 + int_6 + int_7 - int_8)
                  IntPoint(lma_p,lmb  ) = fac_2*RadPart*(int_1 + int_2 - int_3 - int_4 - int_5 - int_6 + int_7 + int_8)
                  IntPoint(lma  ,lmb_p) = fac_3*RadPart*(int_1 - int_2 - int_3 + int_4 + int_5 - int_6 - int_7 + int_8)
                  IntPoint(lma_p,lmb_p) = fac_4*RadPart*(int_1 + int_2 - int_3 - int_4 + int_5 + int_6 - int_7 - int_8)
                end if
              end do
            end if
          end do
        end do
      end do

    case default 

      call stop_all('CalcAngRealBasis', 'The perturbation '//trim(FieldComp(iField))//' is not among the available ones')

    end select

    do la=1,n_l
      do lb=1,n_l
        do ma = 1,2*la-1
          lma = (la-1)**2 + ma
          do mb = 1, 2*lb-1
            lmb = (lb-1)**2 + mb
            if(abs(IntPoint(lma,lmb)).gt.1e-12) &
          &  write(79,'(3i5,2f15.8)') lma, lmb, iField, IntPoint(lma,lmb)
          end do
        end do
      end do
    end do  
  end subroutine CalcAngRealBasis

  subroutine CalcPrimFieldInts(AngPart, iField, ml_interm, m_init)

    complex(idp), allocatable, intent(in) :: AngPart(:,:,:)
    integer, intent(in)  :: iField
    integer, allocatable, intent(in) ::  ml_interm(:), m_init(:)

    complex(idp), pointer :: IntPoint(:,:,:)
    integer :: i, n_l, la, lb, lma, lmb, ma, mb, m1, m2, lmk, lmk_p, lma_p, lmb_p
    integer :: l10, l20
    complex(idp) :: RadPart, pre_fact
    real(idp), parameter :: val_1 = (1.0d0*sqrt_half)
    real(idp), parameter :: val_2 = (1.0d0*sqrt_half)
    complex(idp) :: int_1, int_2, int_3, int_4, int_5, int_6, int_7, int_8
    complex(idp) :: prim_fac_1, prim_fac_2, fac_1, fac_2, fac_3, fac_4
    real(idp) :: msign_1, msign_2, msign_3
    integer, allocatable :: n_m(:), l_finish(:), pos_m_zero(:)

    prim_fac_1= dcmplx(0.0d0, val_1)
    prim_fac_2= dcmplx(val_1, 0.0d0)

    IntPoint =>  PrimFieldInts(:,:,:,iField)

    n_l = para%l + 1

    ! n_m, l_finish and pos_m_zero are used here in the same way that they have been used in 
    ! the 'calc_int_angular_combined' subroutine. More details about them can be found there
    allocate(n_m(n_l), l_finish(n_l), pos_m_zero(n_l))

    do i = 1, n_l
      n_m(i) = min(i,para%ml_max+1)
      l_finish(i) = sum(sph_harm%n_m(1:i))
      pos_m_zero(i) = min(i,para%ml_max+1)
    end do
    
    select case (trim(FieldComp(iField)))
    case('Z')

      lmk = 2
      pre_fact = sqrt(4*pi/3)
      !pre_fact = 1.0d0
      do i= 1, para%ng
        RadPart = pre_fact*grid%r(i)
        do la=1,n_l
          l10 = pos_m_zero(la)
          do lb=1,n_l
            l20 = pos_m_zero(lb)
            do ma = 1, n_m(la)
              lma = ml_interm(la) + ma
              if (ma.eq.l10) then
                do mb = 1, n_m(lb)
                  lmb = ml_interm(lb) + mb
                  int_1 = AngPart(lma,lmb,lmk)
                  if (mb.eq.l20) then
                    IntPoint(lma,lmb,i) = RadPart*int_1
                  else
                    lmb_p = l_finish(lb) - mb + 1
                    m2 = m_init(lb) + mb
                    msign_1 = m_one**m2
                    int_2 = msign_1 * AngPart(lma,lmb_p,lmk)
                    fac_1 = prim_fac_1
                    fac_2 = prim_fac_2
                    IntPoint(lma,lmb  ,i) = fac_1*RadPart*(int_1 - int_2)
                    IntPoint(lma,lmb_p,i) = fac_2*RadPart*(int_1 + int_2)
                  end if
                end do
              else
                lma_p = l_finish(la) - ma + 1
                m1 = m_init(la) + ma
                msign_1 = m_one**m1
                do mb = 1, n_m(lb)
                  lmb = ml_interm(lb) + mb
                  int_1 =           AngPart(lma,lmb,lmk)
                  int_2 = msign_1 * AngPart(lma_p,lmb,lmk)
                  if (mb.eq.l20) then
                    fac_1 = dconjg(prim_fac_1)
                    fac_2 = dconjg(prim_fac_2)
                    IntPoint(lma  ,lmb,i) = fac_1*RadPart*(int_1 - int_2)
                    IntPoint(lma_p,lmb,i) = fac_2*RadPart*(int_1 + int_2)
                  else
                    lmb_p = l_finish(lb) - mb + 1
                    m2 = m_init(lb) + mb
                    msign_2 = m_one**(      m2 )
                    msign_3 = m_one**( m1 + m2 )
                    fac_1 = dconjg(prim_fac_1)*prim_fac_1
                    fac_2 = dconjg(prim_fac_2)*prim_fac_1
                    fac_3 = dconjg(prim_fac_1)*prim_fac_2
                    fac_4 = dconjg(prim_fac_2)*prim_fac_2
                    int_3 = msign_2 * AngPart(lma  ,lmb_p,lmk)
                    int_4 = msign_3 * AngPart(lma_p,lmb_p,lmk)
                    IntPoint(lma  ,lmb  ,i) = fac_1*RadPart*(int_1 - int_2 - int_3 + int_4)
                    IntPoint(lma_p,lmb  ,i) = fac_2*RadPart*(int_1 + int_2 - int_3 - int_4)
                    IntPoint(lma  ,lmb_p,i) = fac_3*RadPart*(int_1 - int_2 + int_3 - int_4)
                    IntPoint(lma_p,lmb_p,i) = fac_4*RadPart*(int_1 + int_2 + int_3 + int_4)
!                   write(81, '(2i4, 6f15.8)') lma, lmb, real(int_1), real(int_2), real(int_3), real(int_4), real(IntPoint(lma_p,lmb,i)), real(IntPoint(lma,lmb_p,i))
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
      pre_fact = prim_fac_1*sqrt(4*pi/3)
      do i= 1, para%ng
        RadPart = pre_fact*grid%r(i)
        do la=1,n_l
          l10 = pos_m_zero(la)
          do lb=1,n_l
            l20 = pos_m_zero(lb)
            do ma = 1, n_m(la)
              lma = ml_interm(la) + ma
              if (ma.eq.l10) then
                do mb = 1, n_m(lb)
                  lmb = ml_interm(lb) + mb
                  int_1 = AngPart(lma,lmb,lmk)
                  int_2 = AngPart(lma,lmb,lmk_p)
                  if (mb.eq.l20) then
                    IntPoint(lma,lmb,i) = RadPart*(int_1 + int_2)
                  else
                    lmb_p = l_finish(lb) - mb + 1
                    m2 = m_init(lb) + mb
                    msign_1 = m_one**m2
                    fac_1 = prim_fac_1
                    fac_2 = prim_fac_2
                    int_3 = msign_1 * AngPart(lma,lmb_p,lmk  )
                    int_4 = msign_1 * AngPart(lma,lmb_p,lmk_p)
                    IntPoint(lma,lmb  ,i) = fac_1*RadPart*(int_1 + int_2 - int_3 - int_4)
                    IntPoint(lma,lmb_p,i) = fac_2*RadPart*(int_1 + int_2 + int_3 + int_4)
                  end if
                end do
              else
                lma_p = l_finish(la) - ma + 1
                m1 = m_init(la) + ma
                msign_1 = m_one**m1
                do mb = 1, n_m(lb)
                  lmb = ml_interm(lb) + mb
                  int_1 =           AngPart(lma  ,lmb,lmk  )
                  int_2 = msign_1 * AngPart(lma_p,lmb,lmk  )
                  int_3 =           AngPart(lma  ,lmb,lmk_p)
                  int_4 = msign_1 * AngPart(lma_p,lmb,lmk_p)
                  if (mb.eq.l20) then
                    fac_1 = dconjg(prim_fac_1)
                    fac_2 = dconjg(prim_fac_2)
                    IntPoint(lma  ,lmb,i) = fac_1*RadPart*(int_1 - int_2 + int_3 - int_4)
                    IntPoint(lma_p,lmb,i) = fac_2*RadPart*(int_1 + int_2 + int_3 + int_4)
                  else
                    lmb_p = l_finish(lb) - mb + 1
                    m2 = m_init(lb) + mb
                    msign_2 = m_one**(      m2 )
                    msign_3 = m_one**( m1 + m2 )
                    int_5 = msign_2 * AngPart(lma  ,lmb_p,lmk  )
                    int_6 = msign_3 * AngPart(lma_p,lmb_p,lmk  )
                    int_7 = msign_2 * AngPart(lma  ,lmb_p,lmk_p)
                    int_8 = msign_3 * AngPart(lma_p,lmb_p,lmk_p)
                    fac_1 = dconjg(prim_fac_1)*prim_fac_1
                    fac_2 = dconjg(prim_fac_2)*prim_fac_1
                    fac_3 = dconjg(prim_fac_1)*prim_fac_2
                    fac_4 = dconjg(prim_fac_2)*prim_fac_2
                    IntPoint(lma  ,lmb  ,i) = fac_1*RadPart*(int_1 - int_2 + int_3 - int_4 - int_5 + int_6 - int_7 + int_8)
                    IntPoint(lma_p,lmb  ,i) = fac_2*RadPart*(int_1 + int_2 + int_3 + int_4 - int_5 - int_6 - int_7 - int_8)
                    IntPoint(lma  ,lmb_p,i) = fac_3*RadPart*(int_1 - int_2 + int_3 - int_4 + int_5 - int_6 + int_7 - int_8)
                    IntPoint(lma_p,lmb_p,i) = fac_4*RadPart*(int_1 + int_2 + int_3 + int_4 + int_5 + int_6 + int_7 + int_8)
                    if (abs(IntPoint(lma_p,lmb_p,i)).gt.0.0d0) then
                        write(95,'(8f15.8)') real(int_1), real(int_2), real(int_3), real(int_4), real(int_5), real(int_6), real(int_7), real(int_8)
                    end if
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
      pre_fact = prim_fac_2*sqrt(4*pi/3)
      !pre_fact = prim_fac_2
      do i= 1, para%ng
        RadPart = pre_fact*grid%r(i)
        do la=1,n_l
          l10 = pos_m_zero(la)
          do lb=1,n_l
            l20 = pos_m_zero(lb)
            do ma = 1, n_m(la)
              lma = ml_interm(la) + ma
              if (ma.eq.l10) then
                do mb = 1, n_m(lb)
                  lmb = ml_interm(lb) + mb
                  int_1 = AngPart(lma,lmb,lmk)
                  int_2 = AngPart(lma,lmb,lmk_p)
                  if (mb.eq.l20) then
                    IntPoint(lma,lmb,i) = RadPart*(int_1 - int_2)
                  else
                    lmb_p = l_finish(lb) - mb + 1
                    m2 = m_init(lb) + mb
                    msign_1 = m_one**m2
                    fac_1 = prim_fac_1
                    fac_2 = prim_fac_2
                    int_3 = msign_1 * AngPart(lma,lmb_p,lmk)
                    int_4 = msign_1 * AngPart(lma,lmb_p,lmk_p)
                    IntPoint(lma,lmb  ,i) = fac_1*RadPart*(int_1 - int_2 - int_3 + int_4)
                    IntPoint(lma,lmb_p,i) = fac_2*RadPart*(int_1 - int_2 + int_3 - int_4)
                  end if
                end do
              else
                lma_p = l_finish(la) - ma + 1
                m1 = m_init(la) + ma
                msign_1 = m_one**m1
                do mb = 1, n_m(lb)
                  lmb = ml_interm(lb) + mb
                  int_1 =           AngPart(lma  ,lmb,lmk  )
                  int_2 = msign_1 * AngPart(lma_p,lmb,lmk  )
                  int_3 =           AngPart(lma  ,lmb,lmk_p)
                  int_4 = msign_1 * AngPart(lma_p,lmb,lmk_p)
                  if (mb.eq.l20) then
                    fac_1 = dconjg(prim_fac_1)
                    fac_2 = dconjg(prim_fac_2)
                    IntPoint(lma  ,lmb,i) = fac_1*RadPart*(int_1 - int_2 - int_3 + int_4)
                    IntPoint(lma_p,lmb,i) = fac_2*RadPart*(int_1 + int_2 - int_3 - int_4)
                  else
                    lmb_p = l_finish(lb) - mb + 1
                    m2 = m_init(lb) + mb
                    msign_2 = m_one**(      m2 )
                    msign_3 = m_one**( m1 + m2 )
                    int_5 = msign_2 * AngPart(lma  ,lmb_p,lmk  )
                    int_6 = msign_3 * AngPart(lma_p,lmb_p,lmk  )
                    int_7 = msign_2 * AngPart(lma  ,lmb_p,lmk_p)
                    int_8 = msign_3 * AngPart(lma_p,lmb_p,lmk_p)

                    fac_1 = dconjg(prim_fac_1)*prim_fac_1
                    fac_2 = dconjg(prim_fac_2)*prim_fac_1
                    fac_3 = dconjg(prim_fac_1)*prim_fac_2
                    fac_4 = dconjg(prim_fac_2)*prim_fac_2
                    IntPoint(lma  ,lmb  ,i) = fac_1*RadPart*(int_1 - int_2 - int_3 + int_4 - int_5 + int_6 + int_7 - int_8)
                    IntPoint(lma_p,lmb  ,i) = fac_2*RadPart*(int_1 + int_2 - int_3 - int_4 - int_5 - int_6 + int_7 + int_8)
                    IntPoint(lma  ,lmb_p,i) = fac_3*RadPart*(int_1 - int_2 - int_3 + int_4 + int_5 - int_6 - int_7 + int_8)
                    IntPoint(lma_p,lmb_p,i) = fac_4*RadPart*(int_1 + int_2 - int_3 - int_4 + int_5 + int_6 - int_7 - int_8)
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

    deallocate(l_finish, pos_m_zero)

  end subroutine CalcPrimFieldInts

  subroutine CombineFieldInts(iField)

    integer, intent(in) :: iField
    complex(idp), allocatable ::  inter_int(:)

    integer :: error, na, nb, la, lb, ma, mb, lma, lmb, n_l,  i
    integer :: ind_1, ind_2
    complex(idp) ::  int_val

    complex(idp), pointer :: PointInts(:,:), PrimPointInts(:,:,:)

    FieldInts = zero

    !allocate(inter_int(para%ng,orb%n_max), stat=error)
    allocate(inter_int(para%ng), stat=error)
    call allocerror(error)

    inter_int = czero
    n_l = para%l + 1

    PointInts => FieldInts(:,:,iField)
    PrimPointInts => PrimFieldInts(:,:,:,iField)

!   if (trim(FieldComp(iField)).eq.'Z') then
!     do i= 1, para%ng
!       do la=1,n_l
!         do ma = 1,2*la-1
!           lma = (la-1)**2 + ma
!           do lb=1,n_l
!             do mb = 1, 2*lb-1
!               lmb = (lb-1)**2 + mb
!               write(78,'(5i4,f15.8)') la, ma, lb, mb, i, real(PrimFieldInts(lma,lmb,i,iField))
!             end do
!           end do
!         end do
!       end do
!     end do
!   end if

    do la = 1, n_l
      do ma = 1, 2*la-1
        lma = (la-1)**2 + ma
        do lb = 1, n_l
          do mb = 1, 2*lb-1
            lmb = (lb-1)**2 + mb
            do na = 1, orb%n_max - (la-1)
              inter_int = czero
              do i = 1, para%ng
                inter_int(i) = inter_int(i) + eigen_vecs(i,na,la)*PrimPointInts(lma,lmb,i)
!               if(abs(PrimPointInts(lma,lmb,i)).gt.1e-12) then
!                   write(77,'(5i4, 3f15.8)') la, ma, lb, mb, i, eigen_vecs(i,na,la), real(PrimPointInts(lma,lmb,i)), inter_int(i)
!               end if
              end do
              do nb = 1, orb%n_max - (lb-1) 
                int_val = 0.0d0
                do i = 1, para%ng
                  int_val = int_val + inter_int(i)*eigen_vecs(i,nb,lb)
!                 if (abs(inter_int(i).gt.1e-12)) then
!                   write(78,'(7i4, 3f15.8)') na, la, ma, nb, lb, mb, i, inter_int(i), eigen_vecs(i, nb, lb), int_val
!                 end if
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

  subroutine ConvertFieldInts(iField, EigVecs, n_m, ml_interm)

    integer, intent(in) :: iField
    real(idp), intent(in) :: EigVecs(:,:,:)
    complex(idp), allocatable ::  inter_int(:,:,:,:)
    integer, allocatable, intent(in) ::  n_m(:), ml_interm(:)

    integer :: error, na, nb, la, lb, ma, mb, lma, lmb, n_l,  i
    integer :: ind_1, ind_2, n_lm
    complex(idp) ::  int_val

    complex(idp), pointer :: PointInts(:,:), PrimPointInts(:,:,:)

    FieldInts = zero

    !allocate(inter_int(para%ng,orb%n_max), stat=error)

    inter_int = czero
    n_l = para%l + 1
    n_lm = para%dim_l

    write(iout,*) 'Converting the integrals for ', trim(FieldComp(iField)), ' from the FE-DVR basis to the contracted basis'

    allocate(inter_int(n_lm, n_lm, orb%n_max, para%ng), stat=error)
    call allocerror(error)
    inter_int = czero

    PointInts => FieldInts(:,:,iField)
    PrimPointInts => PrimFieldInts(:,:,:,iField)

!   if (trim(FieldComp(iField)).eq.'Z') then
!     do i= 1, para%ng
!       do la=1,n_l
!         do ma = 1,2*la-1
!           lma = (la-1)**2 + ma
!           do lb=1,n_l
!             do mb = 1, 2*lb-1
!               lmb = (lb-1)**2 + mb
!               write(78,'(5i4,f15.8)') la, ma, lb, mb, i, real(PrimFieldInts(lma,lmb,i,iField))
!             end do
!           end do
!         end do
!       end do
!     end do
!   end if

    do lb = 1, n_l
      do mb = 1, n_m(lb)
        lmb = ml_interm(lb) + mb
        do la = 1, n_l
          do ma = 1, n_m(la)
            lma = ml_interm(la) + ma
            do na = 1, orb%n_max - (la-1)
              do i = 1, para%ng
                inter_int(lma,lmb, na, i) = inter_int(lma,lmb, na, i) + EigVecs(i,na,la)*PrimPointInts(lma,lmb,i)
!               if(abs(PrimPointInts(lma,lmb,i)).gt.1e-12) then
!                   write(77,'(5i4, 3f15.8)') la, ma, lb, mb, i, EigVecs(i,na,la), real(PrimPointInts(lma,lmb,i)), inter_int(i)
!               end if
              end do
            end do
          end do
        end do
      end do
    end do


    do la = 1, n_l
      do ma = 1, n_m(la)
        lma = ml_interm(la) + ma
        do lb = 1, n_l
          do mb = 1, n_m(lb)
            lmb = ml_interm(lb) + mb
            do na = 1, orb%n_max - (la-1) 
              do nb = 1, orb%n_max - (lb-1) 
                int_val = 0.0d0
                do i = 1, para%ng
                  int_val = int_val + inter_int(lma,lmb, na, i)*EigVecs(i,nb,lb)
!                 if (abs(EigVecs(i,na,la)*PrimPointInts(lma,lmb,i)*EigVecs(i,nb,lb)).gt.1e-12) then
!                   write(84,'(7i4, 3f15.8)') na, la, ma, nb, lb, mb, i, real(PrimPointInts(lma,lmb,i)),real(EigVecs(i,na,la)*PrimPointInts(lma,lmb,i)*EigVecs(i,nb,lb)), real(int_val)
!                 end if
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

!   do la = 1, n_l
!     do ma = 1, 2*la-1
!       lma = (la-1)**2 + ma
!       do lb = 1, n_l
!         do mb = 1, 2*lb-1
!           lmb = (lb-1)**2 + mb
!           do na = 1, orb%n_max - (la-1) 
!             do nb = 1, orb%n_max - (lb-1) 
!               int_val = 0.0d0
!               do i = 1, para%ng
!                 int_val = int_val + EigVecs(i,na,la)*PrimPointInts(lma,lmb,i)*EigVecs(i,nb,lb)
!                 if (abs(EigVecs(i,na,la)*PrimPointInts(lma,lmb,i)*EigVecs(i,nb,lb)).gt.1e-12) then
!                   write(84,'(7i4, 3f15.8)') na, la, ma, nb, lb, mb, i, real(PrimPointInts(lma,lmb,i)),real(EigVecs(i,na,la)*PrimPointInts(lma,lmb,i)*EigVecs(i,nb,lb)), real(int_val)
!                 end if
!               end do
!               ind_1 = SpatialOrbInd(na, la, ma) 
!               ind_2 = SpatialOrbInd(nb, lb, mb)
!               if (ind_1.eq.0) write(iout,*) 'ind_1 is 0 ', na, la, ma
!               if (ind_2.eq.0) write(iout,*) 'ind_2 is 0 ', nb, lb, mb
!               PointInts(ind_1, ind_2) = int_val
!             end do
!           end do  
!         end do
!       end do
!     end do
!   end do

  end subroutine ConvertFieldInts

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
          write(f_int, 1005) real(int_value), i, j, 0, 0
!       if (abs(int_value).gt.tol) &
!       write(f_int, 1006) i, j, real(int_value), aimag(int_value)
      end do 
    end do

    core = 0.0_dp
    if (nFrozen.gt.0) then
      do i = 1, nFrozen
        core = core + 2.0d0*PointInts(i,i)
      end do  
    end if
    write(f_int, 1005) core, 0, 0, 0, 0

    1005 format(f20.16,x,5i4)
    1006 format(2i4, x, 2f20.16)

    close(f_int)
  end subroutine WriteFieldInts

end module FieldIntegrals


