module allocarray
contains


!************************************************************
! ALLOCATE_3
!************************************************************
subroutine alloc_3( num_atom, num_data_total, &
                    num_atom_x, num_atom_y, num_atom_z, &
                    num_data_x, num_data_y, num_data_z, &
                    natom_x_total, natom_xy_total, natom_xz_total, natom_xyz_total, &
                    natom_y_total, natom_yx_total, natom_yz_total, natom_yxz_total, &
                    natom_z_total, natom_zx_total, natom_zy_total, natom_zxy_total, &
                    natom_x_sum,   natom_xy_sum,  natom_xz_sum, natom_xyz_sum, &
                    natom_y_sum,   natom_yx_sum,  natom_yz_sum, natom_yxz_sum, &
                    natom_z_sum,   natom_zx_sum,  natom_zy_sum, natom_zxy_sum, &
                    tag2_x_x,  g2_x_x,  g2_deriv_x_x,  g2_x_x_maxmin,  num_g2_x_x, &
                    tag2_x_y,  g2_x_y,  g2_deriv_x_y,  g2_x_y_maxmin,  num_g2_x_y, &
                    tag2_x_z,  g2_x_z,  g2_deriv_x_z,  g2_x_z_maxmin,  num_g2_x_z, &
                    tag5_x_xx, g5_x_xx, g5_deriv_x_xx, g5_x_xx_maxmin, num_g5_x_xx, &
                    tag5_x_xy, g5_x_xy, g5_deriv_x_xy, g5_x_xy_maxmin, num_g5_x_xy, &
                    tag5_x_xz, g5_x_xz, g5_deriv_x_xz, g5_x_xz_maxmin, num_g5_x_xz, &
                    tag5_x_yy, g5_x_yy, g5_deriv_x_yy, g5_x_yy_maxmin, num_g5_x_yy, &
                    tag5_x_yz, g5_x_yz, g5_deriv_x_yz, g5_x_yz_maxmin, num_g5_x_yz, &
                    tag5_x_zz, g5_x_zz, g5_deriv_x_zz, g5_x_zz_maxmin, num_g5_x_zz, &
                    tag2_y_x,  g2_y_x,  g2_deriv_y_x,  g2_y_x_maxmin,  num_g2_y_x, &
                    tag2_y_y,  g2_y_y,  g2_deriv_y_y,  g2_y_y_maxmin,  num_g2_y_y, &
                    tag2_y_z,  g2_y_z,  g2_deriv_y_z,  g2_y_z_maxmin,  num_g2_y_z, &
                    tag5_y_xx, g5_y_xx, g5_deriv_y_xx, g5_y_xx_maxmin, num_g5_y_xx, &
                    tag5_y_xy, g5_y_xy, g5_deriv_y_xy, g5_y_xy_maxmin, num_g5_y_xy, &
                    tag5_y_xz, g5_y_xz, g5_deriv_y_xz, g5_y_xz_maxmin, num_g5_y_xz, &
                    tag5_y_yy, g5_y_yy, g5_deriv_y_yy, g5_y_yy_maxmin, num_g5_y_yy, &
                    tag5_y_yz, g5_y_yz, g5_deriv_y_yz, g5_y_yz_maxmin, num_g5_y_yz, &
                    tag5_y_zz, g5_y_zz, g5_deriv_y_zz, g5_y_zz_maxmin, num_g5_y_zz, &
                    tag2_z_x,  g2_z_x,  g2_deriv_z_x,  g2_z_x_maxmin,  num_g2_z_x, &
                    tag2_z_y,  g2_z_y,  g2_deriv_z_y,  g2_z_y_maxmin,  num_g2_z_y, &
                    tag2_z_z,  g2_z_z,  g2_deriv_z_z,  g2_z_z_maxmin,  num_g2_z_z, &
                    tag5_z_xx, g5_z_xx, g5_deriv_z_xx, g5_z_xx_maxmin, num_g5_z_xx, &
                    tag5_z_xy, g5_z_xy, g5_deriv_z_xy, g5_z_xy_maxmin, num_g5_z_xy, &
                    tag5_z_xz, g5_z_xz, g5_deriv_z_xz, g5_z_xz_maxmin, num_g5_z_xz, &
                    tag5_z_yy, g5_z_yy, g5_deriv_z_yy, g5_z_yy_maxmin, num_g5_z_yy, &
                    tag5_z_yz, g5_z_yz, g5_deriv_z_yz, g5_z_yz_maxmin, num_g5_z_yz, &
                    tag5_z_zz, g5_z_zz, g5_deriv_z_zz, g5_z_zz_maxmin, num_g5_z_zz )
                        
implicit none
character(len=120) filename
character(len=6) tag2_x_x,  tag2_x_y,  tag2_x_z
character(len=7) tag5_x_xx, tag5_x_xy, tag5_x_xz, &
                 tag5_x_yy, tag5_x_yz, &
                 tag5_x_zz
character(len=6) tag2_y_x,  tag2_y_y,  tag2_y_z
character(len=7) tag5_y_xx, tag5_y_xy, tag5_y_xz, &
                 tag5_y_yy, tag5_y_yz, &
                 tag5_y_zz
character(len=6) tag2_z_x,  tag2_z_y,  tag2_z_z
character(len=7) tag5_z_xx, tag5_z_xy, tag5_z_xz, &
                 tag5_z_yy, tag5_z_yz, &
                 tag5_z_zz
integer num_data_total, num_data_x, num_data_y, num_data_z
integer natom_x_total, natom_xy_total, natom_xz_total, natom_xyz_total
integer natom_y_total, natom_yx_total, natom_yz_total, natom_yxz_total
integer natom_z_total, natom_zx_total, natom_zy_total, natom_zxy_total
integer natom_x_sum,   natom_xy_sum,  natom_xz_sum,  natom_xyz_sum
integer natom_y_sum,   natom_yx_sum,  natom_yz_sum,  natom_yxz_sum
integer natom_z_sum,   natom_zx_sum,  natom_zy_sum,  natom_zxy_sum
integer,allocatable,dimension(:) :: &
  num_atom, num_atom_x, num_atom_y, num_atom_z
integer num_g2_x_x,  num_g2_x_y,  num_g2_x_z, &
        num_g5_x_xx, num_g5_x_xy, num_g5_x_xz, &
        num_g5_x_yy, num_g5_x_yz, &
        num_g5_x_zz
integer num_g2_y_x,  num_g2_y_y,  num_g2_y_z, &
        num_g5_y_xx, num_g5_y_xy, num_g5_y_xz, &
        num_g5_y_yy, num_g5_y_yz, &
        num_g5_y_zz
integer num_g2_z_x,  num_g2_z_y,  num_g2_z_z, &
        num_g5_z_xx, num_g5_z_xy, num_g5_z_xz, &
        num_g5_z_yy, num_g5_z_yz, &
        num_g5_z_zz
double precision,allocatable,dimension(:,:) :: &
  g2_x_x,  g2_x_y,  g2_x_z, &
  g5_x_xx, g5_x_xy, g5_x_xz, &
  g5_x_yy, g5_x_yz, &
  g5_x_zz
double precision,allocatable,dimension(:,:) :: &
  g2_y_x,  g2_y_y,  g2_y_z, &
  g5_y_xx, g5_y_xy, g5_y_xz, &
  g5_y_yy, g5_y_yz, &
  g5_y_zz
double precision,allocatable,dimension(:,:) :: &
  g2_z_x,  g2_z_y,  g2_z_z, &
  g5_z_xx, g5_z_xy, g5_z_xz, &
  g5_z_yy, g5_z_yz, &
  g5_z_zz
double precision,allocatable,dimension(:,:) :: &
  g2_x_x_maxmin,  g2_x_y_maxmin,  g2_x_z_maxmin, &
  g5_x_xx_maxmin, g5_x_xy_maxmin, g5_x_xz_maxmin, &
  g5_x_yy_maxmin, g5_x_yz_maxmin, &
  g5_x_zz_maxmin
double precision,allocatable,dimension(:,:) :: &
  g2_y_x_maxmin,  g2_y_y_maxmin,  g2_y_z_maxmin, &
  g5_y_xx_maxmin, g5_y_xy_maxmin, g5_y_xz_maxmin, &
  g5_y_yy_maxmin, g5_y_yz_maxmin, &
  g5_y_zz_maxmin
double precision,allocatable,dimension(:,:) :: &
  g2_z_x_maxmin,  g2_z_y_maxmin,  g2_z_z_maxmin, &
  g5_z_xx_maxmin, g5_z_xy_maxmin, g5_z_xz_maxmin, &
  g5_z_yy_maxmin, g5_z_yz_maxmin, &
  g5_z_zz_maxmin
double precision,allocatable,dimension(:,:,:) :: &
  g2_deriv_x_x,  g2_deriv_x_y,  g2_deriv_x_z, &
  g5_deriv_x_xx, g5_deriv_x_xy, g5_deriv_x_xz, &
  g5_deriv_x_yy, g5_deriv_x_yz, &
  g5_deriv_x_zz
double precision,allocatable,dimension(:,:,:) :: &
  g2_deriv_y_x,  g2_deriv_y_y,  g2_deriv_y_z, &
  g5_deriv_y_xx, g5_deriv_y_xy, g5_deriv_y_xz, &
  g5_deriv_y_yy, g5_deriv_y_yz, &
  g5_deriv_y_zz
double precision,allocatable,dimension(:,:,:) :: &
  g2_deriv_z_x,  g2_deriv_z_y,  g2_deriv_z_z, &
  g5_deriv_z_xx, g5_deriv_z_xy, g5_deriv_z_xz, &
  g5_deriv_z_yy, g5_deriv_z_yz, &
  g5_deriv_z_zz


 allocate( num_atom( num_data_total ) )

!X
 allocate( num_atom_x( num_data_x ) )

 allocate( g2_x_x(  natom_x_total,   num_g2_x_x ) )
 allocate( g2_x_y(  natom_xy_total,  num_g2_x_y ) )
 allocate( g2_x_z(  natom_xz_total,  num_g2_x_z ) )
 allocate( g5_x_xx( natom_x_total,   num_g5_x_xx ) )
 allocate( g5_x_xy( natom_xy_total,  num_g5_x_xy ) )
 allocate( g5_x_xz( natom_xz_total,  num_g5_x_xz ) )
 allocate( g5_x_yy( natom_xy_total,  num_g5_x_yy ) )
 allocate( g5_x_yz( natom_xyz_total, num_g5_x_yz ) )
 allocate( g5_x_zz( natom_xz_total,  num_g5_x_zz ) )

 allocate( g2_x_x_maxmin(  num_g2_x_x, 2 ) )
 allocate( g2_x_y_maxmin(  num_g2_x_y, 2 ) )
 allocate( g2_x_z_maxmin(  num_g2_x_z, 2 ) )
 allocate( g5_x_xx_maxmin( num_g5_x_xx, 2 ) )
 allocate( g5_x_xy_maxmin( num_g5_x_xy, 2 ) )
 allocate( g5_x_xz_maxmin( num_g5_x_xz, 2 ) )
 allocate( g5_x_yy_maxmin( num_g5_x_yy, 2 ) )
 allocate( g5_x_yz_maxmin( num_g5_x_yz, 2 ) )
 allocate( g5_x_zz_maxmin( num_g5_x_zz, 2 ) )

 allocate( g2_deriv_x_x(  natom_x_sum,   3, num_g2_x_x ) )
 allocate( g2_deriv_x_y(  natom_xy_sum,  3, num_g2_x_y ) )
 allocate( g2_deriv_x_z(  natom_xz_sum,  3, num_g2_x_z ) )
 allocate( g5_deriv_x_xx( natom_x_sum,   3, num_g5_x_xx ) )
 allocate( g5_deriv_x_xy( natom_xy_sum,  3, num_g5_x_xy ) )
 allocate( g5_deriv_x_xz( natom_xz_sum,  3, num_g5_x_xz ) )
 allocate( g5_deriv_x_yy( natom_xy_sum,  3, num_g5_x_yy ) )
 allocate( g5_deriv_x_yz( natom_xyz_sum, 3, num_g5_x_yz ) )
 allocate( g5_deriv_x_zz( natom_xz_sum,  3, num_g5_x_zz ) )

 filename = './data_maxmin/'//trim(adjustl(tag2_x_x))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g2_x_x_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag2_x_y))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g2_x_y_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag2_x_z))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g2_x_z_maxmin
 close(20)

 filename = './data_maxmin/'//trim(adjustl(tag5_x_xx))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_x_xx_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag5_x_xy))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_x_xy_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag5_x_xz))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_x_xz_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag5_x_yy))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_x_yy_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag5_x_yz))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_x_yz_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag5_x_zz))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_x_zz_maxmin
 close(20)

!Y
 allocate( num_atom_y( num_data_y ) )

 allocate( g2_y_x(  natom_yx_total,  num_g2_y_x ) )
 allocate( g2_y_y(  natom_y_total,   num_g2_y_y ) )
 allocate( g2_y_z(  natom_yz_total,  num_g2_y_z ) )
 allocate( g5_y_xx( natom_yx_total,  num_g5_y_xx ) )
 allocate( g5_y_xy( natom_yx_total,  num_g5_y_xy ) )
 allocate( g5_y_xz( natom_yxz_total, num_g5_y_xz ) )
 allocate( g5_y_yy( natom_y_total,   num_g5_y_yy ) )
 allocate( g5_y_yz( natom_yz_total,  num_g5_y_yz ) )
 allocate( g5_y_zz( natom_yz_total,  num_g5_y_zz ) )

 allocate( g2_y_x_maxmin(  num_g2_y_x, 2 ) )
 allocate( g2_y_y_maxmin(  num_g2_y_y, 2 ) )
 allocate( g2_y_z_maxmin(  num_g2_y_z, 2 ) )
 allocate( g5_y_xx_maxmin( num_g5_y_xx, 2 ) )
 allocate( g5_y_xy_maxmin( num_g5_y_xy, 2 ) )
 allocate( g5_y_xz_maxmin( num_g5_y_xz, 2 ) )
 allocate( g5_y_yy_maxmin( num_g5_y_yy, 2 ) )
 allocate( g5_y_yz_maxmin( num_g5_y_yz, 2 ) )
 allocate( g5_y_zz_maxmin( num_g5_y_zz, 2 ) )

 allocate( g2_deriv_y_x(  natom_yx_sum,  3, num_g2_y_x ) )
 allocate( g2_deriv_y_y(  natom_y_sum,   3, num_g2_y_y ) )
 allocate( g2_deriv_y_z(  natom_yz_sum,  3, num_g2_y_z ) )
 allocate( g5_deriv_y_xx( natom_yx_sum,  3, num_g5_y_xx ) )
 allocate( g5_deriv_y_xy( natom_yx_sum,  3, num_g5_y_xy ) )
 allocate( g5_deriv_y_xz( natom_yxz_sum, 3, num_g5_y_xz ) )
 allocate( g5_deriv_y_yy( natom_y_sum,   3, num_g5_y_yy ) )
 allocate( g5_deriv_y_yz( natom_yz_sum,  3, num_g5_y_yz ) )
 allocate( g5_deriv_y_zz( natom_yz_sum,  3, num_g5_y_zz ) )

 filename = './data_maxmin/'//trim(adjustl(tag2_y_x))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g2_y_x_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag2_y_y))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g2_y_y_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag2_y_z))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g2_y_z_maxmin
 close(20)

 filename = './data_maxmin/'//trim(adjustl(tag5_y_xx))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_y_xx_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag5_y_xy))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_y_xy_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag5_y_xz))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_y_xz_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag5_y_yy))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_y_yy_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag5_y_yz))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_y_yz_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag5_y_zz))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_y_zz_maxmin
 close(20)

!Z
 allocate( num_atom_z( num_data_z ) )

 allocate( g2_z_x(  natom_zx_total,  num_g2_z_x ) )
 allocate( g2_z_y(  natom_zy_total,  num_g2_z_y ) )
 allocate( g2_z_z(  natom_z_total,   num_g2_z_z ) )
 allocate( g5_z_xx( natom_zx_total,  num_g5_z_xx ) )
 allocate( g5_z_xy( natom_zxy_total, num_g5_z_xy ) )
 allocate( g5_z_xz( natom_zx_total,  num_g5_z_xz ) )
 allocate( g5_z_yy( natom_zy_total,  num_g5_z_yy ) )
 allocate( g5_z_yz( natom_zy_total,  num_g5_z_yz ) )
 allocate( g5_z_zz( natom_z_total,   num_g5_z_zz ) )

 allocate( g2_z_x_maxmin(  num_g2_z_x, 2 ) )
 allocate( g2_z_y_maxmin(  num_g2_z_y, 2 ) )
 allocate( g2_z_z_maxmin(  num_g2_z_z, 2 ) )
 allocate( g5_z_xx_maxmin( num_g5_z_xx, 2 ) )
 allocate( g5_z_xy_maxmin( num_g5_z_xy, 2 ) )
 allocate( g5_z_xz_maxmin( num_g5_z_xz, 2 ) )
 allocate( g5_z_yy_maxmin( num_g5_z_yy, 2 ) )
 allocate( g5_z_yz_maxmin( num_g5_z_yz, 2 ) )
 allocate( g5_z_zz_maxmin( num_g5_z_zz, 2 ) )

 allocate( g2_deriv_z_x(  natom_zx_sum,  3, num_g2_z_x ) )
 allocate( g2_deriv_z_y(  natom_zy_sum,  3, num_g2_z_y ) )
 allocate( g2_deriv_z_z(  natom_z_sum,   3, num_g2_z_z ) )
 allocate( g5_deriv_z_xx( natom_zx_sum,  3, num_g5_z_xx ) )
 allocate( g5_deriv_z_xy( natom_zxy_sum, 3, num_g5_z_xy ) )
 allocate( g5_deriv_z_xz( natom_zx_sum,  3, num_g5_z_xz ) )
 allocate( g5_deriv_z_yy( natom_zy_sum,  3, num_g5_z_yy ) )
 allocate( g5_deriv_z_yz( natom_zy_sum,  3, num_g5_z_yz ) )
 allocate( g5_deriv_z_zz( natom_z_sum,   3, num_g5_z_zz ) )

 filename = './data_maxmin/'//trim(adjustl(tag2_z_x))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g2_z_x_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag2_z_y))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g2_z_y_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag2_z_z))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g2_z_z_maxmin
 close(20)

 filename = './data_maxmin/'//trim(adjustl(tag5_z_xx))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_z_xx_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag5_z_xy))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_z_xy_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag5_z_xz))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_z_xz_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag5_z_yy))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_z_yy_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag5_z_yz))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_z_yz_maxmin
 close(20)
 filename = './data_maxmin/'//trim(adjustl(tag5_z_zz))//'_maxmin.dat'
 open(unit=20,file=filename,action="read",form="unformatted")
 read(20) g5_z_zz_maxmin
 close(20)


end subroutine alloc_3


end module
