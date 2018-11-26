   ! The total one-electron integral for an external dipole field is calculated here numerically using Cartesian grid points. 
   ! While the spherical harmonics and the radial part of the wavefunctions are calculated separately in these grid points and then later multiplied in order to obtain the integrals

   program NumericalInts
   implicit none
   integer i,j,k,l,m,ngrid,n,mm,ic,igrid, n_l, nao, orb, j_val
   integer :: m_val, l_val, tao, il, iao, lm
   double precision check,r,r12,rm,rn,pi,check6d,rval
   double precision, allocatable ::  wt(:),x(:),y(:),z(:), coords(:,:)
   double precision ::  wt_val, x_val, y_val, z_val, val_1, val_2, start, finish, fac_1, fac_2
   double precision, allocatable :: coords2(:,:)
   complex*16, allocatable :: zpsi(:,:),zpsi_ang(:,:),zpsi_rad(:,:), zoverlap(:,:)
   complex*16, allocatable :: dipx(:,:), dipy(:,:), dipz(:,:)
   complex*16, allocatable :: inter_ints_1(:,:,:,:), inter_ints_2(:,:,:,:)
   complex*16 :: zi,zcheck1,zcheck2,zcheck, zval, zval_1
   complex*16 :: transf_l_0, transf_l_1(3,3), transf_l_2(5,6)
   complex*16  :: get_val
   character*80 grid
   character*32 :: file_name_syntx 
   character*8 :: int2str
   logical :: real_func
   
   ! First define some common constants

   zi = dcmplx(0.0d0,1.0d0)
   pi = 3.141592653589d0


   ! Integrals can be obtained for spherical harmonics in both of their complex and real forms.
   ! The real spherical harmonics are obtained by combining the original descriptions of spherical
   ! harmonics which has the complex form
   real_func = .false.

   !file_name_syntx = 'angular_element'
    
   ! Total number of 'l' quantum number to be used to define the spherical harmonics
   n_l = 2
   ! Total number of spherical harmonics
   nao = n_l**2

   ! Sph. Harm. for l=0 does not depend on x, y and z
   transf_l_0 = 1.d0/2.d0/sqrt(pi)

   ! Initiating the transformation matrices for l=1,2
   transf_l_1 = dcmplx(0.0d0,0.0d0)

   transf_l_2 = dcmplx(0.0d0,0.0d0)

   ! transf_l_1 is contracted with the column vector (x,z,y) to produce the 3 spher. harm. for l=1
   ! transf_l_2 is contracted with the column vector (xy,yz,zx,x^2,y^2,z^2) to produce the 5 spher. harm. for l=2
   if (real_func) then
     write(*,*) 'Using real Spherical Harmonics'
     ! Elements of transf_l_1 to obtain real sph. harm. 
     val_1 = 1.0d0/2.d0*sqrt(3.d0/pi)
     transf_l_1(3,1) = val_1
     transf_l_1(2,2) = val_1
     transf_l_1(1,3) = val_1

     ! Elements of transf_l_2 to obtain real sph. harm. 
     val_1 = 1.0d0/2.d0*sqrt(15.d0/pi)
     val_2 = 1.0d0/4.d0*sqrt(5.d0/pi)
     transf_l_2(1,1) = val_1
     transf_l_2(2,2) = val_1
     !transf_l_2(2,3) = val_1
     !transf_l_2(2,2) = zi*val_1
     transf_l_2(3,4) = -1.0d0*val_2
     transf_l_2(3,5) = -1.0d0*val_2
     transf_l_2(3,6) =  2.0d0*val_2
     transf_l_2(4,3) = val_1
     !transf_l_2(4,2) = val_1
     !transf_l_2(4,3) = dconjg(zi)*val_1
     transf_l_2(5,4) = 0.5d0*val_1
     transf_l_2(5,5) = -0.5d0*val_1
   else 
     ! Elements of transf_l_1 to obtain complex sph. harm. 
     write(*,*) 'Using imaginary Spherical Harmonics'
     transf_l_1(1,1) = 1.0d0/2.d0*sqrt(3.d0/(2.d0*pi))
     transf_l_1(1,3) = -1.0d0*zi/2.d0*sqrt(3.d0/(2.d0*pi))
     transf_l_1(2,2) = 1.0d0/2.d0*sqrt(3.d0/(pi))
     transf_l_1(3,1) = -1.0d0/2.d0*sqrt(3.d0/(2.d0*pi))
     transf_l_1(3,3) = -1.0d0*zi/2.d0*sqrt(3.d0/(2.d0*pi))

     ! Elements of transf_l_2 to obtain complex sph. harm. 
     val_1 = 1.0d0/4.d0*sqrt(15.d0/2.0d0/pi)
     transf_l_2(1,1) = dconjg(zi)*2.0d0*val_1
     transf_l_2(1,4) = val_1
     transf_l_2(1,5) = -1.0d0*val_1
     transf_l_2(5,1) = zi*2.0d0*val_1
     transf_l_2(5,4) = val_1
     transf_l_2(5,5) = -1.0d0*val_1
     zval_1 = 1.0d0/2.d0*sqrt(15.d0/2.0d0/pi)
     !zval_1 = zi*1.0d0/2.d0*sqrt(15.d0/2.0d0/pi)
     transf_l_2(2,2) = dconjg(zi)*zval_1
     transf_l_2(2,3) = -1.0d0*zval_1
     zval_1 = 1.0d0/2.d0*sqrt(15.d0/2.0d0/pi)
     !zval_1 = dconjg(zi)*1.0d0/2.d0*sqrt(15.d0/2.0d0/pi)
     transf_l_2(4,2) = dconjg(zi)*zval_1
     transf_l_2(4,3) =  1.0d0*zval_1
     val_1 = 1.0d0/4.d0*sqrt(5.d0/pi)
     transf_l_2(3,4) = -1.0d0*val_1
     transf_l_2(3,5) = -1.0d0*val_1
     transf_l_2(3,6) = 2.0d0*val_1
   end if

   ! calculating the total number of orbitals involved
   tao = 0 
   do i = 1, 2
     do j = 1,i 
       do k = 1,2*j-1
       tao = tao + 1
       end do
     end do
   end do

   write(*,*) 'total number of orbitals', tao

   write(grid,'(a,i0)') '3dgrid'
   open(81,file=trim(grid))
   ! read the total number of 3d grid points
   read(81,*) ngrid
   write(*,*) 'Total number point:', ngrid
   
   allocate(zpsi_rad(3,ngrid), zpsi(tao,ngrid))
   allocate(wt(ngrid), x(ngrid), y(ngrid), z(ngrid), coords(3,ngrid), coords2(6,ngrid))

   do m = 1,ngrid
     read(81,*) x_val,y_val,z_val,wt_val
     r = dsqrt(x_val**2+y_val**2+z_val**2)

     x(m) = x_val
     y(m) = y_val
     z(m) = z_val
     wt(m) = wt_val

     ! coords contains the grid points for the column vector (x,z,y)
     coords(1,m) = x(m)/r
     coords(2,m) = z(m)/r
     coords(3,m) = y(m)/r

     ! Calculating the radial part of the wave function for the Hydrogen atom for n=0 and n=1
     ! and storing it in the array zpsi_rad
     zpsi_rad(1,m) = 2.0d0*exp(-r)
     zpsi_rad(2,m) = (2.0d0-r)*exp(-r/2.0d0)/sqrt(8.0d0)
     zpsi_rad(3,m) = r*exp(-r/2.0d0)/sqrt(24.0d0)

!    coords2(1,m) = x(m)*y(m)
!    coords2(2,m) = y(m)*z(m)
!    coords2(3,m) = z(m)*x(m)
!    coords2(4,m) = x(m)*x(m)
!    coords2(5,m) = y(m)*y(m)
!    coords2(6,m) = z(m)*z(m)
   enddo

   
! 3d check

   check  = 0.d0
   zcheck1 = 0.d0
   zcheck2 = 0.d0
   
   ! zpsi_ang stores the angular part of the wave function
   ! zoverlap stores the overlap between different orbitals
   allocate(zpsi_ang(nao, ngrid), zoverlap(tao,tao))
   allocate(inter_ints_1(ngrid,nao,nao,nao), inter_ints_2(ngrid,nao,nao,nao))
   allocate(dipx(tao,tao), dipy(tao,tao), dipz(tao,tao))

   zoverlap = 0.d0
   
   zpsi_ang = dcmplx(0.d0, 0.d0)
   ! Get zpsi_ang for l=0
   zpsi_ang(1,:) = 1.d0/2.d0/sqrt(pi)

   do m = 1,ngrid

     i = 2
     ! Get zpsi_ang for l=1
     do j = 1, 2*i-1
       orb = (i-1)**2 + j
       do k = 1,3
         zpsi_ang(orb,m) = zpsi_ang(orb,m) + transf_l_1(j,k)*coords(k,m)
       end do
     end do
    
     ! Multiply the radial and angular part of the wave function to form the full orbitals
     iao = 0
     il = 0
     do i = 1, 2
       do l = 1, i
       il = il + 1
         do lm = (l-1)**2 + 1, l**2
           iao = iao+1
           zpsi(iao, m) =  zpsi_rad(il,m)*zpsi_ang(lm, m)
         end do
       end do
     end do   
   end do
   
   ! Calculate the overlap between orbitals
   do m = 1,ngrid
     do i = 1,tao
      do j = 1,tao
        zoverlap(i,j) = zoverlap(i,j) + wt(m)*dconjg(zpsi(i,m))*zpsi(j,m)
      enddo
     enddo
   enddo
   
   
   do i = 1,tao
    do j = 1,tao
      if (abs(zoverlap(i,j)).gt.1e-12) write(6,*) i,j,zoverlap(i,j)
!     write(6,*) i,j,zoverlap(i,j)
    enddo
   enddo

   ! Here calculating the integrals for all three components of the dipole operators. 
   dipz = 0.0d0
   dipy = 0.0d0
   dipx = 0.0d0
   fac_1 = sqrt(4.0d0*pi/3.0d0)
!  fac_2 = sqrt(2.0d0*pi/3.0d0)
!  fac_1 = 1.0d0
!  fac_2 = 1.0d0
   do m = 1, ngrid
     r = sqrt(x(m)**2 + y(m)**2 + z(m)**2)
     do i =  1, tao
       do j =  1, tao
          dipz(i,j) = dipz(i,j) + fac_1*wt(m)*dconjg(zpsi(i,m))*r*zpsi_ang(3,m)*zpsi(j,m)
          dipx(i,j) = dipx(i,j) + fac_1*wt(m)*dconjg(zpsi(i,m))*r*zpsi_ang(4,m)*zpsi(j,m)
          dipy(i,j) = dipy(i,j) + fac_1*wt(m)*dconjg(zpsi(i,m))*r*zpsi_ang(2,m)*zpsi(j,m)
       end do
     end do
   end do

   
   write(*,*) 'DIPX:'
   do i = 1,tao
     do j = 1,i
       if (abs(dipx(i,j)).gt.1e-12) write(6,*) i,j,dipx(i,j)
     enddo
   enddo

   write(*,*) 'DIPY:'
   do i = 1,tao
     do j = 1,i
       if (abs(dipy(i,j)).gt.1e-12) write(6,*) i,j,dipy(i,j)
     enddo
   enddo

   write(*,*) 'DIPZ:'
   do i = 1,tao
     do j = 1,i
       if (abs(dipz(i,j)).gt.1e-12) write(6,*) i,j,dipz(i,j)
     enddo
   enddo

   deallocate(wt,x,y,z,zpsi,zoverlap,zpsi_rad, zpsi_ang, coords)

   end program NumericalInts

   pure function get_val(x,y,z,l,m) result(zval)

   implicit none

   double precision, intent(in) :: x,y,z
   integer, intent(in) :: l, m
   complex*16 zval, zi
   double precision :: pi

   pi = 3.141592653589d0

   zi = dcmplx(0.0,1.0)

   zval = dcmplx(0.0d0,0.0d0)

   if (l.eq.1) then
     zval = 1.d0/2.d0/sqrt(pi)
   elseif (l.eq.2) then
     if (m.eq.-1) then
       zval = (x - zi*y)/2.d0*sqrt(3/(2*pi))
     elseif(m.eq.0) then
       zval = z/2.d0*sqrt(3/(pi))
     else
       zval = -(x+zi*y)/2.d0*sqrt(3/(2*pi))
     end if
   end if
     
   end function get_val

!      do i = 2, n_l
!        do j = 1, 2*i-1
!          j_val = -1*i + j 
!          orb = (i-1)**2 + j 
!          zpsi(m,orb) =  get_val(x(m),y(m),z(m),i,j_val)
!        end do
!      end do

!      zpsi(m,1) = 1.d0/2.d0/sqrt(pi)
!      zpsi(m,2) = (x(m) - zi*y(m))/2.d0*sqrt(3/(2*pi))
!      zpsi(m,3) = z(m)/2.d0*sqrt(3/(pi))
!      zpsi(m,4) =  -(x(m)+zi*y(m))/2.d0*sqrt(3/(2*pi))

  character(len=8) function int2str(i, format)

    integer,                    intent(in) :: i
    character(len=*), optional, intent(in) :: format

    if (present(format)) then
      write(int2str, format) i
    else
      write(int2str, '(I7)') i
    end if
    int2str = adjustl(int2str)

  end function int2str
