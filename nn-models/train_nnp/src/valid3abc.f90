subroutine valid3abc( &
  io_x, io_y, io_z, &
  natom, natom_x, natom_y, natom_z, &
  natom_x_total,   counter_natom_x,   natom_x_sum,   counter_deriv_x, &
  natom_xy_total,  counter_natom_xy,  natom_xy_sum,  counter_deriv_xy, &
  natom_xz_total,  counter_natom_xz,  natom_xz_sum,  counter_deriv_xz, &
  natom_xyz_total, counter_natom_xyz, natom_xyz_sum, counter_deriv_xyz, &
  natom_y_total,   counter_natom_y,   natom_y_sum,   counter_deriv_y, &
  natom_yx_total,  counter_natom_yx,  natom_yx_sum,  counter_deriv_yx, &
  natom_yz_total,  counter_natom_yz,  natom_yz_sum,  counter_deriv_yz, &
  natom_yxz_total, counter_natom_yxz, natom_yxz_sum, counter_deriv_yxz, &
  natom_z_total,   counter_natom_z,   natom_z_sum,   counter_deriv_z, &
  natom_zx_total,  counter_natom_zx,  natom_zx_sum,  counter_deriv_zx, &
  natom_zy_total,  counter_natom_zy,  natom_zy_sum,  counter_deriv_zy, &
  natom_zxy_total, counter_natom_zxy, natom_zxy_sum, counter_deriv_zxy, &
  g2_x_x_total,  g2_deriv_x_x_total,  num_g2_x_x,  g2_x_x_maxmin, &
  g2_x_y_total,  g2_deriv_x_y_total,  num_g2_x_y,  g2_x_y_maxmin, &
  g2_x_z_total,  g2_deriv_x_z_total,  num_g2_x_z,  g2_x_z_maxmin, &
  g5_x_xx_total, g5_deriv_x_xx_total, num_g5_x_xx, g5_x_xx_maxmin, &
  g5_x_xy_total, g5_deriv_x_xy_total, num_g5_x_xy, g5_x_xy_maxmin, &
  g5_x_xz_total, g5_deriv_x_xz_total, num_g5_x_xz, g5_x_xz_maxmin, &
  g5_x_yy_total, g5_deriv_x_yy_total, num_g5_x_yy, g5_x_yy_maxmin, &
  g5_x_yz_total, g5_deriv_x_yz_total, num_g5_x_yz, g5_x_yz_maxmin, &
  g5_x_zz_total, g5_deriv_x_zz_total, num_g5_x_zz, g5_x_zz_maxmin, &
  network_x, nlayer_x, num_g_x, &
  weight_x, nweight_x, hidden_x, nhidden_x, &
  g2_y_x_total,  g2_deriv_y_x_total,  num_g2_y_x,  g2_y_x_maxmin, &
  g2_y_y_total,  g2_deriv_y_y_total,  num_g2_y_y,  g2_y_y_maxmin, &
  g2_y_z_total,  g2_deriv_y_z_total,  num_g2_y_z,  g2_y_z_maxmin, &
  g5_y_xx_total, g5_deriv_y_xx_total, num_g5_y_xx, g5_y_xx_maxmin, &
  g5_y_xy_total, g5_deriv_y_xy_total, num_g5_y_xy, g5_y_xy_maxmin, &
  g5_y_xz_total, g5_deriv_y_xz_total, num_g5_y_xz, g5_y_xz_maxmin, &
  g5_y_yy_total, g5_deriv_y_yy_total, num_g5_y_yy, g5_y_yy_maxmin, &
  g5_y_yz_total, g5_deriv_y_yz_total, num_g5_y_yz, g5_y_yz_maxmin, &
  g5_y_zz_total, g5_deriv_y_zz_total, num_g5_y_zz, g5_y_zz_maxmin, &
  network_y, nlayer_y, num_g_y, &
  weight_y, nweight_y, hidden_y, nhidden_y, &
  g2_z_x_total,  g2_deriv_z_x_total,  num_g2_z_x,  g2_z_x_maxmin, &
  g2_z_y_total,  g2_deriv_z_y_total,  num_g2_z_y,  g2_z_y_maxmin, &
  g2_z_z_total,  g2_deriv_z_z_total,  num_g2_z_z,  g2_z_z_maxmin, &
  g5_z_xx_total, g5_deriv_z_xx_total, num_g5_z_xx, g5_z_xx_maxmin, &
  g5_z_xy_total, g5_deriv_z_xy_total, num_g5_z_xy, g5_z_xy_maxmin, &
  g5_z_xz_total, g5_deriv_z_xz_total, num_g5_z_xz, g5_z_xz_maxmin, &
  g5_z_yy_total, g5_deriv_z_yy_total, num_g5_z_yy, g5_z_yy_maxmin, &
  g5_z_yz_total, g5_deriv_z_yz_total, num_g5_z_yz, g5_z_yz_maxmin, &
  g5_z_zz_total, g5_deriv_z_zz_total, num_g5_z_zz, g5_z_zz_maxmin, &
  network_z, nlayer_z, num_g_z, &
  weight_z, nweight_z, hidden_z, nhidden_z, &
  energy, force )

implicit none
integer i, j, k, m
integer io_x, io_y, io_z
integer natom, natom_x, natom_y, natom_z
integer natom_x_total,   counter_natom_x,   natom_x_sum,   counter_deriv_x, &
        natom_xy_total,  counter_natom_xy,  natom_xy_sum,  counter_deriv_xy, &
        natom_xz_total,  counter_natom_xz,  natom_xz_sum,  counter_deriv_xz, &
        natom_xyz_total, counter_natom_xyz, natom_xyz_sum, counter_deriv_xyz, &
        natom_y_total,   counter_natom_y,   natom_y_sum,   counter_deriv_y, &
        natom_yx_total,  counter_natom_yx,  natom_yx_sum,  counter_deriv_yx, &
        natom_yz_total,  counter_natom_yz,  natom_yz_sum,  counter_deriv_yz, &
        natom_yxz_total, counter_natom_yxz, natom_yxz_sum, counter_deriv_yxz, &
        natom_z_total,   counter_natom_z,   natom_z_sum,   counter_deriv_z, &
        natom_zx_total,  counter_natom_zx,  natom_zx_sum,  counter_deriv_zx, &
        natom_zy_total,  counter_natom_zy,  natom_zy_sum,  counter_deriv_zy, &
        natom_zxy_total, counter_natom_zxy, natom_zxy_sum, counter_deriv_zxy
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
integer num_g_x(9), num_g_y(9), num_g_z(9)
double precision &
  g2_x_x_total(  natom_x_total,   num_g2_x_x ), &
  g2_x_y_total(  natom_xy_total,  num_g2_x_y ), &
  g2_x_z_total(  natom_xz_total,  num_g2_x_z ), &
  g5_x_xx_total( natom_x_total,   num_g5_x_xx ), &
  g5_x_xy_total( natom_xy_total,  num_g5_x_xy ), &
  g5_x_xz_total( natom_xz_total,  num_g5_x_xz ), &
  g5_x_yy_total( natom_xy_total,  num_g5_x_yy ), &
  g5_x_yz_total( natom_xyz_total, num_g5_x_yz ), &
  g5_x_zz_total( natom_xz_total,  num_g5_x_zz )
double precision &
  g2_y_x_total(  natom_yx_total,  num_g2_y_x ), &
  g2_y_y_total(  natom_y_total,   num_g2_y_y ), &
  g2_y_z_total(  natom_yz_total,  num_g2_y_z ), &
  g5_y_xx_total( natom_yx_total,  num_g5_y_xx ), &
  g5_y_xy_total( natom_yx_total,  num_g5_y_xy ), &
  g5_y_xz_total( natom_yxz_total, num_g5_y_xz ), &
  g5_y_yy_total( natom_y_total,   num_g5_y_yy ), &
  g5_y_yz_total( natom_yz_total,  num_g5_y_yz ), &
  g5_y_zz_total( natom_yz_total,  num_g5_y_zz )
double precision &
  g2_z_x_total(  natom_zx_total,  num_g2_z_x ), &
  g2_z_y_total(  natom_zy_total,  num_g2_z_y ), &
  g2_z_z_total(  natom_z_total,   num_g2_z_z ), &
  g5_z_xx_total( natom_zx_total,  num_g5_z_xx ), &
  g5_z_xy_total( natom_zxy_total, num_g5_z_xy ), &
  g5_z_xz_total( natom_zx_total,  num_g5_z_xz ), &
  g5_z_yy_total( natom_zy_total,  num_g5_z_yy ), &
  g5_z_yz_total( natom_zy_total,  num_g5_z_yz ), &
  g5_z_zz_total( natom_z_total,   num_g5_z_zz )
double precision &
  g2_deriv_x_x_total(  natom_x_sum,   3, num_g2_x_x ), &
  g2_deriv_x_y_total(  natom_xy_sum,  3, num_g2_x_y ), &
  g2_deriv_x_z_total(  natom_xz_sum,  3, num_g2_x_z ), &
  g5_deriv_x_xx_total( natom_x_sum,   3, num_g5_x_xx ), &
  g5_deriv_x_xy_total( natom_xy_sum,  3, num_g5_x_xy ), &
  g5_deriv_x_xz_total( natom_xz_sum,  3, num_g5_x_xz ), &
  g5_deriv_x_yy_total( natom_xy_sum,  3, num_g5_x_yy ), &
  g5_deriv_x_yz_total( natom_xyz_sum, 3, num_g5_x_yz ), &
  g5_deriv_x_zz_total( natom_xz_sum,  3, num_g5_x_zz )
double precision &
  g2_deriv_y_x_total(  natom_yx_sum,  3, num_g2_y_x ), &
  g2_deriv_y_y_total(  natom_y_sum,   3, num_g2_y_y ), &
  g2_deriv_y_z_total(  natom_yz_sum,  3, num_g2_y_z ), &
  g5_deriv_y_xx_total( natom_yx_sum,  3, num_g5_y_xx ), &
  g5_deriv_y_xy_total( natom_yx_sum,  3, num_g5_y_xy ), &
  g5_deriv_y_xz_total( natom_yxz_sum, 3, num_g5_y_xz ), &
  g5_deriv_y_yy_total( natom_y_sum,   3, num_g5_y_yy ), &
  g5_deriv_y_yz_total( natom_yz_sum,  3, num_g5_y_yz ), &
  g5_deriv_y_zz_total( natom_yz_sum,  3, num_g5_y_zz )
double precision &
  g2_deriv_z_x_total(  natom_zx_sum,  3, num_g2_z_x ), &
  g2_deriv_z_y_total(  natom_zy_sum,  3, num_g2_z_y ), &
  g2_deriv_z_z_total(  natom_z_sum,   3, num_g2_z_z ), &
  g5_deriv_z_xx_total( natom_zx_sum,  3, num_g5_z_xx ), &
  g5_deriv_z_xy_total( natom_zxy_sum, 3, num_g5_z_xy ), &
  g5_deriv_z_xz_total( natom_zx_sum,  3, num_g5_z_xz ), &
  g5_deriv_z_yy_total( natom_zy_sum,  3, num_g5_z_yy ), &
  g5_deriv_z_yz_total( natom_zy_sum,  3, num_g5_z_yz ), &
  g5_deriv_z_zz_total( natom_z_sum,   3, num_g5_z_zz )
double precision &
  g2_x_x(  natom_x, num_g2_x_x ), &
  g2_x_y(  natom_x, num_g2_x_y ), &
  g2_x_z(  natom_x, num_g2_x_z ), &
  g5_x_xx( natom_x, num_g5_x_xx ), &
  g5_x_xy( natom_x, num_g5_x_xy ), &
  g5_x_xz( natom_x, num_g5_x_xz ), &
  g5_x_yy( natom_x, num_g5_x_yy ), &
  g5_x_yz( natom_x, num_g5_x_yz ), &
  g5_x_zz( natom_x, num_g5_x_zz )
double precision &
  g2_y_x(  natom_y, num_g2_y_x ), &
  g2_y_y(  natom_y, num_g2_y_y ), &
  g2_y_z(  natom_y, num_g2_y_z ), &
  g5_y_xx( natom_y, num_g5_y_xx ), &
  g5_y_xy( natom_y, num_g5_y_xy ), &
  g5_y_xz( natom_y, num_g5_y_xz ), &
  g5_y_yy( natom_y, num_g5_y_yy ), &
  g5_y_yz( natom_y, num_g5_y_yz ), &
  g5_y_zz( natom_y, num_g5_y_zz )
double precision &
  g2_z_x(  natom_z, num_g2_z_x ), &
  g2_z_y(  natom_z, num_g2_z_y ), &
  g2_z_z(  natom_z, num_g2_z_z ), &
  g5_z_xx( natom_z, num_g5_z_xx ), &
  g5_z_xy( natom_z, num_g5_z_xy ), &
  g5_z_xz( natom_z, num_g5_z_xz ), &
  g5_z_yy( natom_z, num_g5_z_yy ), &
  g5_z_yz( natom_z, num_g5_z_yz ), &
  g5_z_zz( natom_z, num_g5_z_zz )
double precision &
  g2_deriv_x_x(  natom_x, natom, 3, num_g2_x_x ), &
  g2_deriv_x_y(  natom_x, natom, 3, num_g2_x_y ), &
  g2_deriv_x_z(  natom_x, natom, 3, num_g2_x_z ), &
  g5_deriv_x_xx( natom_x, natom, 3, num_g5_x_xx ), &
  g5_deriv_x_xy( natom_x, natom, 3, num_g5_x_xy ), &
  g5_deriv_x_xz( natom_x, natom, 3, num_g5_x_xz ), &
  g5_deriv_x_yy( natom_x, natom, 3, num_g5_x_yy ), &
  g5_deriv_x_yz( natom_x, natom, 3, num_g5_x_yz ), &
  g5_deriv_x_zz( natom_x, natom, 3, num_g5_x_zz )
double precision &
  g2_deriv_y_x(  natom_y, natom, 3, num_g2_y_x ), &
  g2_deriv_y_y(  natom_y, natom, 3, num_g2_y_y ), &
  g2_deriv_y_z(  natom_y, natom, 3, num_g2_y_z ), &
  g5_deriv_y_xx( natom_y, natom, 3, num_g5_y_xx ), &
  g5_deriv_y_xy( natom_y, natom, 3, num_g5_y_xy ), &
  g5_deriv_y_xz( natom_y, natom, 3, num_g5_y_xz ), &
  g5_deriv_y_yy( natom_y, natom, 3, num_g5_y_yy ), &
  g5_deriv_y_yz( natom_y, natom, 3, num_g5_y_yz ), &
  g5_deriv_y_zz( natom_y, natom, 3, num_g5_y_zz )
double precision &
  g2_deriv_z_x(  natom_z, natom, 3, num_g2_z_x ), &
  g2_deriv_z_y(  natom_z, natom, 3, num_g2_z_y ), &
  g2_deriv_z_z(  natom_z, natom, 3, num_g2_z_z ), &
  g5_deriv_z_xx( natom_z, natom, 3, num_g5_z_xx ), &
  g5_deriv_z_xy( natom_z, natom, 3, num_g5_z_xy ), &
  g5_deriv_z_xz( natom_z, natom, 3, num_g5_z_xz ), &
  g5_deriv_z_yy( natom_z, natom, 3, num_g5_z_yy ), &
  g5_deriv_z_yz( natom_z, natom, 3, num_g5_z_yz ), &
  g5_deriv_z_zz( natom_z, natom, 3, num_g5_z_zz )
double precision :: &
  g2_x_x_maxmin(  num_g2_x_x, 2 ), &
  g2_x_y_maxmin(  num_g2_x_y, 2 ), &
  g2_x_z_maxmin(  num_g2_x_z, 2 ), &
  g5_x_xx_maxmin( num_g5_x_xx, 2 ), &
  g5_x_xy_maxmin( num_g5_x_xy, 2 ), &
  g5_x_xz_maxmin( num_g5_x_xz, 2 ), &
  g5_x_yy_maxmin( num_g5_x_yy, 2 ), &
  g5_x_yz_maxmin( num_g5_x_yz, 2 ), &
  g5_x_zz_maxmin( num_g5_x_zz, 2 )
double precision :: &
  g2_y_x_maxmin(  num_g2_y_x, 2 ), &
  g2_y_y_maxmin(  num_g2_y_y, 2 ), &
  g2_y_z_maxmin(  num_g2_y_z, 2 ), &
  g5_y_xx_maxmin( num_g5_y_xx, 2 ), &
  g5_y_xy_maxmin( num_g5_y_xy, 2 ), &
  g5_y_xz_maxmin( num_g5_y_xz, 2 ), &
  g5_y_yy_maxmin( num_g5_y_yy, 2 ), &
  g5_y_yz_maxmin( num_g5_y_yz, 2 ), &
  g5_y_zz_maxmin( num_g5_y_zz, 2 )
double precision :: &
  g2_z_x_maxmin(  num_g2_z_x, 2 ), &
  g2_z_y_maxmin(  num_g2_z_y, 2 ), &
  g2_z_z_maxmin(  num_g2_z_z, 2 ), &
  g5_z_xx_maxmin( num_g5_z_xx, 2 ), &
  g5_z_xy_maxmin( num_g5_z_xy, 2 ), &
  g5_z_xz_maxmin( num_g5_z_xz, 2 ), &
  g5_z_yy_maxmin( num_g5_z_yy, 2 ), &
  g5_z_yz_maxmin( num_g5_z_yz, 2 ), &
  g5_z_zz_maxmin( num_g5_z_zz, 2 )
integer nlayer_x, nweight_x, nhidden_x
integer nlayer_y, nweight_y, nhidden_y
integer nlayer_z, nweight_z, nhidden_z
integer network_x( nlayer_x )
integer network_y( nlayer_y )
integer network_z( nlayer_z )
double precision weight_x( nweight_x )
double precision weight_y( nweight_y )
double precision weight_z( nweight_z )
double precision hidden_x( natom_x, nhidden_x )
double precision hidden_y( natom_y, nhidden_y )
double precision hidden_z( natom_z, nhidden_z )
double precision energy_x
double precision energy_y
double precision energy_z
double precision energy
double precision force( natom, 3 )


!------
! X_SF
!------

!X-X
 if( io_x == 1 .and. io_y == 0 .and. io_z == 0 )then

 do j = 1, num_g2_x_x
   do i = 1, natom_x
     g2_x_x(i,j) = g2_x_x_total( counter_natom_x + i, j )
   enddo
 enddo
 do j = 1, num_g2_x_y
   do i = 1, natom_x
     g2_x_y(i,j) = ( -2.0d0*g2_x_y_maxmin(j,2)/( g2_x_y_maxmin(j,1) - g2_x_y_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo
 do j = 1, num_g2_x_z
   do i = 1, natom_x
     g2_x_z(i,j) = ( -2.0d0*g2_x_z_maxmin(j,2)/( g2_x_z_maxmin(j,1) - g2_x_z_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo
 do j = 1, num_g5_x_xx
   do i = 1, natom_x
     g5_x_xx(i,j) = g5_x_xx_total( counter_natom_x + i, j )
   enddo
 enddo
 do j = 1, num_g5_x_xy
   do i = 1, natom_x
     g5_x_xy(i,j) = ( -2.0d0*g5_x_xy_maxmin(j,2)/( g5_x_xy_maxmin(j,1) - g5_x_xy_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo
 do j = 1, num_g5_x_xz
   do i = 1, natom_x
     g5_x_xz(i,j) = ( -2.0d0*g5_x_xz_maxmin(j,2)/( g5_x_xz_maxmin(j,1) - g5_x_xz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo
 do j = 1, num_g5_x_yy
   do i = 1, natom_x
     g5_x_yy(i,j) = ( -2.0d0*g5_x_yy_maxmin(j,2)/( g5_x_yy_maxmin(j,1) - g5_x_yy_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo
 do j = 1, num_g5_x_yz
   do i = 1, natom_x
     g5_x_yz(i,j) = ( -2.0d0*g5_x_yz_maxmin(j,2)/( g5_x_yz_maxmin(j,1) - g5_x_yz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo
 do j = 1, num_g5_x_zz
   do i = 1, natom_x
     g5_x_zz(i,j) = ( -2.0d0*g5_x_zz_maxmin(j,2)/( g5_x_zz_maxmin(j,1) - g5_x_zz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo

 counter_natom_x   = counter_natom_x   + natom_x

 endif

!X-XY
 if( io_x == 1 .and. io_y == 1 .and. io_z == 0 )then

 do j = 1, num_g2_x_x
   do i = 1, natom_x
     g2_x_x(i,j) = g2_x_x_total( counter_natom_x + i, j )
   enddo
 enddo
 do j = 1, num_g2_x_y
   do i = 1, natom_x
     g2_x_y(i,j) = g2_x_y_total( counter_natom_xy + i, j )
   enddo
 enddo
 do j = 1, num_g2_x_z
   do i = 1, natom_x
     g2_x_z(i,j) = ( -2.0d0*g2_x_z_maxmin(j,2)/( g2_x_z_maxmin(j,1) - g2_x_z_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo
 do j = 1, num_g5_x_xx
   do i = 1, natom_x
     g5_x_xx(i,j) = g5_x_xx_total( counter_natom_x + i, j )
   enddo
 enddo
 do j = 1, num_g5_x_xy
   do i = 1, natom_x
     g5_x_xy(i,j) = g5_x_xy_total( counter_natom_xy + i, j )
   enddo
 enddo
 do j = 1, num_g5_x_xz
   do i = 1, natom_x
     g5_x_xz(i,j) = ( -2.0d0*g5_x_xz_maxmin(j,2)/( g5_x_xz_maxmin(j,1) - g5_x_xz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo
 do j = 1, num_g5_x_yy
   do i = 1, natom_x
     g5_x_yy(i,j) = g5_x_yy_total( counter_natom_xy + i, j )
   enddo
 enddo
 do j = 1, num_g5_x_yz
   do i = 1, natom_x
     g5_x_yz(i,j) = ( -2.0d0*g5_x_yz_maxmin(j,2)/( g5_x_yz_maxmin(j,1) - g5_x_yz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo
 do j = 1, num_g5_x_zz
   do i = 1, natom_x
     g5_x_zz(i,j) = ( -2.0d0*g5_x_zz_maxmin(j,2)/( g5_x_zz_maxmin(j,1) - g5_x_zz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo

 counter_natom_x   = counter_natom_x   + natom_x
 counter_natom_xy  = counter_natom_xy  + natom_x

 endif

!X-XZ
 if( io_x == 1 .and. io_y == 0 .and. io_z == 1 )then


 do j = 1, num_g2_x_x
   do i = 1, natom_x
     g2_x_x(i,j) = g2_x_x_total( counter_natom_x + i, j )
   enddo
 enddo


 do j = 1, num_g2_x_y
   do i = 1, natom_x
     g2_x_y(i,j) = ( -2.0d0*g2_x_y_maxmin(j,2)/( g2_x_y_maxmin(j,1) - g2_x_y_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g2_x_z
   do i = 1, natom_x
     g2_x_z(i,j) = g2_x_z_total( counter_natom_xz + i, j )
   enddo
 enddo


 do j = 1, num_g5_x_xx
   do i = 1, natom_x
     g5_x_xx(i,j) = g5_x_xx_total( counter_natom_x + i, j )
   enddo
 enddo


 do j = 1, num_g5_x_xy
   do i = 1, natom_x
     g5_x_xy(i,j) = ( -2.0d0*g5_x_xy_maxmin(j,2)/( g5_x_xy_maxmin(j,1) - g5_x_xy_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_x_xz
   do i = 1, natom_x
     g5_x_xz(i,j) = g5_x_xz_total( counter_natom_xz + i, j )
   enddo
 enddo


 do j = 1, num_g5_x_yy
   do i = 1, natom_x
     g5_x_yy(i,j) = ( -2.0d0*g5_x_yy_maxmin(j,2)/( g5_x_yy_maxmin(j,1) - g5_x_yy_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_x_yz
   do i = 1, natom_x
     g5_x_yz(i,j) = ( -2.0d0*g5_x_yz_maxmin(j,2)/( g5_x_yz_maxmin(j,1) - g5_x_yz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_x_zz
   do i = 1, natom_x
     g5_x_zz(i,j) = g5_x_zz_total( counter_natom_xz + i, j )
   enddo
 enddo


 counter_natom_x   = counter_natom_x   + natom_x
 counter_natom_xz  = counter_natom_xz  + natom_x

 endif

!X-XYZ
 if( io_x == 1 .and. io_y == 1 .and. io_z == 1 )then


 do j = 1, num_g2_x_x
   do i = 1, natom_x
     g2_x_x(i,j) = g2_x_x_total( counter_natom_x + i, j )
   enddo
 enddo


 do j = 1, num_g2_x_y
   do i = 1, natom_x
     g2_x_y(i,j) = g2_x_y_total( counter_natom_xy + i, j )
   enddo
 enddo


 do j = 1, num_g2_x_z
   do i = 1, natom_x
     g2_x_z(i,j) = g2_x_z_total( counter_natom_xz + i, j )
   enddo
 enddo


 do j = 1, num_g5_x_xx
   do i = 1, natom_x
     g5_x_xx(i,j) = g5_x_xx_total( counter_natom_x + i, j )
   enddo
 enddo


 do j = 1, num_g5_x_xy
   do i = 1, natom_x
     g5_x_xy(i,j) = g5_x_xy_total( counter_natom_xy + i, j )
   enddo
 enddo


 do j = 1, num_g5_x_xz
   do i = 1, natom_x
     g5_x_xz(i,j) = g5_x_xz_total( counter_natom_xz + i, j )
   enddo
 enddo


 do j = 1, num_g5_x_yy
   do i = 1, natom_x
     g5_x_yy(i,j) = g5_x_yy_total( counter_natom_xy + i, j )
   enddo
 enddo


 do j = 1, num_g5_x_yz
   do i = 1, natom_x
     g5_x_yz(i,j) = g5_x_yz_total( counter_natom_xyz + i, j )
   enddo
 enddo


 do j = 1, num_g5_x_zz
   do i = 1, natom_x
     g5_x_zz(i,j) = g5_x_zz_total( counter_natom_xz + i, j )
   enddo
 enddo


 counter_natom_x   = counter_natom_x   + natom_x
 counter_natom_xy  = counter_natom_xy  + natom_x
 counter_natom_xz  = counter_natom_xz  + natom_x
 counter_natom_xyz = counter_natom_xyz + natom_x

 endif

!X-Y, X-Z, X-YZ
 if( io_x == 0 )then


 do j = 1, num_g2_x_x
   do i = 1, natom_x
     g2_x_x(i,j) = ( -2.0d0*g2_x_x_maxmin(j,2)/( g2_x_x_maxmin(j,1) - g2_x_x_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g2_x_y
   do i = 1, natom_x
     g2_x_y(i,j) = ( -2.0d0*g2_x_y_maxmin(j,2)/( g2_x_y_maxmin(j,1) - g2_x_y_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g2_x_z
   do i = 1, natom_x
     g2_x_z(i,j) = ( -2.0d0*g2_x_z_maxmin(j,2)/( g2_x_z_maxmin(j,1) - g2_x_z_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_x_xx
   do i = 1, natom_x
     g5_x_xx(i,j) = ( -2.0d0*g5_x_xx_maxmin(j,2)/( g5_x_xx_maxmin(j,1) - g5_x_xx_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_x_xy
   do i = 1, natom_x
     g5_x_xy(i,j) = ( -2.0d0*g5_x_xy_maxmin(j,2)/( g5_x_xy_maxmin(j,1) - g5_x_xy_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_x_xz
   do i = 1, natom_x
     g5_x_xz(i,j) = ( -2.0d0*g5_x_xz_maxmin(j,2)/( g5_x_xz_maxmin(j,1) - g5_x_xz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_x_yy
   do i = 1, natom_x
     g5_x_yy(i,j) = ( -2.0d0*g5_x_yy_maxmin(j,2)/( g5_x_yy_maxmin(j,1) - g5_x_yy_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_x_yz
   do i = 1, natom_x
     g5_x_yz(i,j) = ( -2.0d0*g5_x_yz_maxmin(j,2)/( g5_x_yz_maxmin(j,1) - g5_x_yz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_x_zz
   do i = 1, natom_x
     g5_x_zz(i,j) = ( -2.0d0*g5_x_zz_maxmin(j,2)/( g5_x_zz_maxmin(j,1) - g5_x_zz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 endif

!------------
! X_SF_deriv
!------------

!X-X
 if( io_x == 1 .and. io_y == 0 .and. io_z == 0 )then


 do k = 1, num_g2_x_x
   do j = 1, natom
     do i = 1, natom_x
       m = i + natom_x * ( j - 1 )
       g2_deriv_x_x(i,j,1,k) = g2_deriv_x_x_total( counter_deriv_x + m, 1, k )
       g2_deriv_x_x(i,j,2,k) = g2_deriv_x_x_total( counter_deriv_x + m, 2, k )
       g2_deriv_x_x(i,j,3,k) = g2_deriv_x_x_total( counter_deriv_x + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g2_x_y
   do j = 1, natom
     do i = 1, natom_x
       g2_deriv_x_y(i,j,1,k) = 0.0d0
       g2_deriv_x_y(i,j,2,k) = 0.0d0
       g2_deriv_x_y(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g2_x_z
   do j = 1, natom
     do i = 1, natom_x
       g2_deriv_x_z(i,j,1,k) = 0.0d0
       g2_deriv_x_z(i,j,2,k) = 0.0d0
       g2_deriv_x_z(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_xx
   do j = 1, natom
     do i = 1, natom_x
       m = i + natom_x * ( j - 1 )
       g5_deriv_x_xx(i,j,1,k) = g5_deriv_x_xx_total( counter_deriv_x + m, 1, k )
       g5_deriv_x_xx(i,j,2,k) = g5_deriv_x_xx_total( counter_deriv_x + m, 2, k )
       g5_deriv_x_xx(i,j,3,k) = g5_deriv_x_xx_total( counter_deriv_x + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_xy
   do j = 1, natom
     do i = 1, natom_x
       g5_deriv_x_xy(i,j,1,k) = 0.0d0
       g5_deriv_x_xy(i,j,2,k) = 0.0d0
       g5_deriv_x_xy(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_xz
   do j = 1, natom
     do i = 1, natom_x
       g5_deriv_x_xz(i,j,1,k) = 0.0d0
       g5_deriv_x_xz(i,j,2,k) = 0.0d0
       g5_deriv_x_xz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_yy
   do j = 1, natom
     do i = 1, natom_x
       g5_deriv_x_yy(i,j,1,k) = 0.0d0
       g5_deriv_x_yy(i,j,2,k) = 0.0d0
       g5_deriv_x_yy(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_yz
   do j = 1, natom
     do i = 1, natom_x
       g5_deriv_x_yz(i,j,1,k) = 0.0d0
       g5_deriv_x_yz(i,j,2,k) = 0.0d0
       g5_deriv_x_yz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_zz
   do j = 1, natom
     do i = 1, natom_x
       g5_deriv_x_zz(i,j,1,k) = 0.0d0
       g5_deriv_x_zz(i,j,2,k) = 0.0d0
       g5_deriv_x_zz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 counter_deriv_x   = counter_deriv_x   + natom_x * natom

 endif

!X-XY
 if( io_x == 1 .and. io_y == 1 .and. io_z == 0 )then


 do k = 1, num_g2_x_x
   do j = 1, natom
     do i = 1, natom_x
       m = i + natom_x * ( j - 1 )
       g2_deriv_x_x(i,j,1,k) = g2_deriv_x_x_total( counter_deriv_x + m, 1, k )
       g2_deriv_x_x(i,j,2,k) = g2_deriv_x_x_total( counter_deriv_x + m, 2, k )
       g2_deriv_x_x(i,j,3,k) = g2_deriv_x_x_total( counter_deriv_x + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g2_x_y
   do j = 1, natom
     do i = 1, natom_x
       m = i + natom_x * ( j - 1 )
       g2_deriv_x_y(i,j,1,k) = g2_deriv_x_y_total( counter_deriv_xy + m, 1, k )
       g2_deriv_x_y(i,j,2,k) = g2_deriv_x_y_total( counter_deriv_xy + m, 2, k )
       g2_deriv_x_y(i,j,3,k) = g2_deriv_x_y_total( counter_deriv_xy + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g2_x_z
   do j = 1, natom
     do i = 1, natom_x
       g2_deriv_x_z(i,j,1,k) = 0.0d0
       g2_deriv_x_z(i,j,2,k) = 0.0d0
       g2_deriv_x_z(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_xx
   do j = 1, natom
     do i = 1, natom_x
       m = i + natom_x * ( j - 1 )
       g5_deriv_x_xx(i,j,1,k) = g5_deriv_x_xx_total( counter_deriv_x + m, 1, k )
       g5_deriv_x_xx(i,j,2,k) = g5_deriv_x_xx_total( counter_deriv_x + m, 2, k )
       g5_deriv_x_xx(i,j,3,k) = g5_deriv_x_xx_total( counter_deriv_x + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_xy
   do j = 1, natom
     do i = 1, natom_x
       m = i + natom_x * ( j - 1 )
       g5_deriv_x_xy(i,j,1,k) = g5_deriv_x_xy_total( counter_deriv_xy + m, 1, k )
       g5_deriv_x_xy(i,j,2,k) = g5_deriv_x_xy_total( counter_deriv_xy + m, 2, k )
       g5_deriv_x_xy(i,j,3,k) = g5_deriv_x_xy_total( counter_deriv_xy + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_xz
   do j = 1, natom
     do i = 1, natom_x
       g5_deriv_x_xz(i,j,1,k) = 0.0d0
       g5_deriv_x_xz(i,j,2,k) = 0.0d0
       g5_deriv_x_xz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_yy
   do j = 1, natom
     do i = 1, natom_x
       m = i + natom_x * ( j - 1 )
       g5_deriv_x_yy(i,j,1,k) = g5_deriv_x_yy_total( counter_deriv_xy + m, 1, k )
       g5_deriv_x_yy(i,j,2,k) = g5_deriv_x_yy_total( counter_deriv_xy + m, 2, k )
       g5_deriv_x_yy(i,j,3,k) = g5_deriv_x_yy_total( counter_deriv_xy + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_yz
   do j = 1, natom
     do i = 1, natom_x
       g5_deriv_x_yz(i,j,1,k) = 0.0d0
       g5_deriv_x_yz(i,j,2,k) = 0.0d0
       g5_deriv_x_yz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_zz
   do j = 1, natom
     do i = 1, natom_x
       g5_deriv_x_zz(i,j,1,k) = 0.0d0
       g5_deriv_x_zz(i,j,2,k) = 0.0d0
       g5_deriv_x_zz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 counter_deriv_x   = counter_deriv_x   + natom_x * natom
 counter_deriv_xy  = counter_deriv_xy  + natom_x * natom

 endif

!X-XZ
 if( io_x == 1 .and. io_y == 0 .and. io_z == 1 )then


 do k = 1, num_g2_x_x
   do j = 1, natom
     do i = 1, natom_x
       m = i + natom_x * ( j - 1 )
       g2_deriv_x_x(i,j,1,k) = g2_deriv_x_x_total( counter_deriv_x + m, 1, k )
       g2_deriv_x_x(i,j,2,k) = g2_deriv_x_x_total( counter_deriv_x + m, 2, k )
       g2_deriv_x_x(i,j,3,k) = g2_deriv_x_x_total( counter_deriv_x + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g2_x_y
   do j = 1, natom
     do i = 1, natom_x
       g2_deriv_x_y(i,j,1,k) = 0.0d0
       g2_deriv_x_y(i,j,2,k) = 0.0d0
       g2_deriv_x_y(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g2_x_z
   do j = 1, natom
     do i = 1, natom_x
       m = i + natom_x * ( j - 1 )
       g2_deriv_x_z(i,j,1,k) = g2_deriv_x_z_total( counter_deriv_xz + m, 1, k )
       g2_deriv_x_z(i,j,2,k) = g2_deriv_x_z_total( counter_deriv_xz + m, 2, k )
       g2_deriv_x_z(i,j,3,k) = g2_deriv_x_z_total( counter_deriv_xz + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_xx
   do j = 1, natom
     do i = 1, natom_x
       m = i + natom_x * ( j - 1 )
       g5_deriv_x_xx(i,j,1,k) = g5_deriv_x_xx_total( counter_deriv_x + m, 1, k )
       g5_deriv_x_xx(i,j,2,k) = g5_deriv_x_xx_total( counter_deriv_x + m, 2, k )
       g5_deriv_x_xx(i,j,3,k) = g5_deriv_x_xx_total( counter_deriv_x + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_xy
   do j = 1, natom
     do i = 1, natom_x
       g5_deriv_x_xy(i,j,1,k) = 0.0d0
       g5_deriv_x_xy(i,j,2,k) = 0.0d0
       g5_deriv_x_xy(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_xz
   do j = 1, natom
     do i = 1, natom_x
       m = i + natom_x * ( j - 1 )
       g5_deriv_x_xz(i,j,1,k) = g5_deriv_x_xz_total( counter_deriv_xz + m, 1, k )
       g5_deriv_x_xz(i,j,2,k) = g5_deriv_x_xz_total( counter_deriv_xz + m, 2, k )
       g5_deriv_x_xz(i,j,3,k) = g5_deriv_x_xz_total( counter_deriv_xz + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_yy
   do j = 1, natom
     do i = 1, natom_x
       g5_deriv_x_yy(i,j,1,k) = 0.0d0
       g5_deriv_x_yy(i,j,2,k) = 0.0d0
       g5_deriv_x_yy(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_yz
   do j = 1, natom
     do i = 1, natom_x
       g5_deriv_x_yz(i,j,1,k) = 0.0d0
       g5_deriv_x_yz(i,j,2,k) = 0.0d0
       g5_deriv_x_yz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_zz
   do j = 1, natom
     do i = 1, natom_x
       m = i + natom_x * ( j - 1 )
       g5_deriv_x_zz(i,j,1,k) = g5_deriv_x_zz_total( counter_deriv_xz + m, 1, k )
       g5_deriv_x_zz(i,j,2,k) = g5_deriv_x_zz_total( counter_deriv_xz + m, 2, k )
       g5_deriv_x_zz(i,j,3,k) = g5_deriv_x_zz_total( counter_deriv_xz + m, 3, k )
     enddo
   enddo
 enddo


 counter_deriv_x   = counter_deriv_x   + natom_x * natom
 counter_deriv_xz  = counter_deriv_xz  + natom_x * natom

 endif

!X-XYZ
 if( io_x == 1 .and. io_y == 1 .and. io_z == 1 )then


 do k = 1, num_g2_x_x
   do j = 1, natom
     do i = 1, natom_x
       m = i + natom_x * ( j - 1 )
       g2_deriv_x_x(i,j,1,k) = g2_deriv_x_x_total( counter_deriv_x + m, 1, k )
       g2_deriv_x_x(i,j,2,k) = g2_deriv_x_x_total( counter_deriv_x + m, 2, k )
       g2_deriv_x_x(i,j,3,k) = g2_deriv_x_x_total( counter_deriv_x + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g2_x_y
   do j = 1, natom
     do i = 1, natom_x
       m = i + natom_x * ( j - 1 )
       g2_deriv_x_y(i,j,1,k) = g2_deriv_x_y_total( counter_deriv_xy + m, 1, k )
       g2_deriv_x_y(i,j,2,k) = g2_deriv_x_y_total( counter_deriv_xy + m, 2, k )
       g2_deriv_x_y(i,j,3,k) = g2_deriv_x_y_total( counter_deriv_xy + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g2_x_z
   do j = 1, natom
     do i = 1, natom_x
       m = i + natom_x * ( j - 1 )
       g2_deriv_x_z(i,j,1,k) = g2_deriv_x_z_total( counter_deriv_xz + m, 1, k )
       g2_deriv_x_z(i,j,2,k) = g2_deriv_x_z_total( counter_deriv_xz + m, 2, k )
       g2_deriv_x_z(i,j,3,k) = g2_deriv_x_z_total( counter_deriv_xz + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_xx
   do j = 1, natom
     do i = 1, natom_x
       m = i + natom_x * ( j - 1 )
       g5_deriv_x_xx(i,j,1,k) = g5_deriv_x_xx_total( counter_deriv_x + m, 1, k )
       g5_deriv_x_xx(i,j,2,k) = g5_deriv_x_xx_total( counter_deriv_x + m, 2, k )
       g5_deriv_x_xx(i,j,3,k) = g5_deriv_x_xx_total( counter_deriv_x + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_xy
   do j = 1, natom
     do i = 1, natom_x
       m = i + natom_x * ( j - 1 )
       g5_deriv_x_xy(i,j,1,k) = g5_deriv_x_xy_total( counter_deriv_xy + m, 1, k )
       g5_deriv_x_xy(i,j,2,k) = g5_deriv_x_xy_total( counter_deriv_xy + m, 2, k )
       g5_deriv_x_xy(i,j,3,k) = g5_deriv_x_xy_total( counter_deriv_xy + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_xz
   do j = 1, natom
     do i = 1, natom_x
       m = i + natom_x * ( j - 1 )
       g5_deriv_x_xz(i,j,1,k) = g5_deriv_x_xz_total( counter_deriv_xz + m, 1, k )
       g5_deriv_x_xz(i,j,2,k) = g5_deriv_x_xz_total( counter_deriv_xz + m, 2, k )
       g5_deriv_x_xz(i,j,3,k) = g5_deriv_x_xz_total( counter_deriv_xz + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_yy
   do j = 1, natom
     do i = 1, natom_x
       m = i + natom_x * ( j - 1 )
       g5_deriv_x_yy(i,j,1,k) = g5_deriv_x_yy_total( counter_deriv_xy + m, 1, k )
       g5_deriv_x_yy(i,j,2,k) = g5_deriv_x_yy_total( counter_deriv_xy + m, 2, k )
       g5_deriv_x_yy(i,j,3,k) = g5_deriv_x_yy_total( counter_deriv_xy + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_yz
   do j = 1, natom
     do i = 1, natom_x
       m = i + natom_x * ( j - 1 )
       g5_deriv_x_yz(i,j,1,k) = g5_deriv_x_yz_total( counter_deriv_xyz + m, 1, k )
       g5_deriv_x_yz(i,j,2,k) = g5_deriv_x_yz_total( counter_deriv_xyz + m, 2, k )
       g5_deriv_x_yz(i,j,3,k) = g5_deriv_x_yz_total( counter_deriv_xyz + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_zz
   do j = 1, natom
     do i = 1, natom_x
       m = i + natom_x * ( j - 1 )
       g5_deriv_x_zz(i,j,1,k) = g5_deriv_x_zz_total( counter_deriv_xz + m, 1, k )
       g5_deriv_x_zz(i,j,2,k) = g5_deriv_x_zz_total( counter_deriv_xz + m, 2, k )
       g5_deriv_x_zz(i,j,3,k) = g5_deriv_x_zz_total( counter_deriv_xz + m, 3, k )
     enddo
   enddo
 enddo


 counter_deriv_x   = counter_deriv_x   + natom_x * natom
 counter_deriv_xy  = counter_deriv_xy  + natom_x * natom
 counter_deriv_xz  = counter_deriv_xz  + natom_x * natom
 counter_deriv_xyz = counter_deriv_xyz + natom_x * natom 

 endif

!X-Y, X-Z, X-YZ
 if( io_x == 0 )then


 do k = 1, num_g2_x_x
   do j = 1, natom
     do i = 1, natom_x
       g2_deriv_x_x(i,j,1,k) = 0.0d0
       g2_deriv_x_x(i,j,2,k) = 0.0d0
       g2_deriv_x_x(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g2_x_y
   do j = 1, natom
     do i = 1, natom_x
       g2_deriv_x_y(i,j,1,k) = 0.0d0
       g2_deriv_x_y(i,j,2,k) = 0.0d0
       g2_deriv_x_y(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g2_x_z
   do j = 1, natom
     do i = 1, natom_x
       g2_deriv_x_z(i,j,1,k) = 0.0d0
       g2_deriv_x_z(i,j,2,k) = 0.0d0
       g2_deriv_x_z(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_xx
   do j = 1, natom
     do i = 1, natom_x
       g5_deriv_x_xx(i,j,1,k) = 0.0d0
       g5_deriv_x_xx(i,j,2,k) = 0.0d0
       g5_deriv_x_xx(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_xy
   do j = 1, natom
     do i = 1, natom_x
       g5_deriv_x_xy(i,j,1,k) = 0.0d0
       g5_deriv_x_xy(i,j,2,k) = 0.0d0
       g5_deriv_x_xy(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_xz
   do j = 1, natom
     do i = 1, natom_x
       g5_deriv_x_xz(i,j,1,k) = 0.0d0
       g5_deriv_x_xz(i,j,2,k) = 0.0d0
       g5_deriv_x_xz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_yy
   do j = 1, natom
     do i = 1, natom_x
       g5_deriv_x_yy(i,j,1,k) = 0.0d0
       g5_deriv_x_yy(i,j,2,k) = 0.0d0
       g5_deriv_x_yy(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_yz
   do j = 1, natom
     do i = 1, natom_x
       g5_deriv_x_yz(i,j,1,k) = 0.0d0
       g5_deriv_x_yz(i,j,2,k) = 0.0d0
       g5_deriv_x_yz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_x_zz
   do j = 1, natom
     do i = 1, natom_x
       g5_deriv_x_zz(i,j,1,k) = 0.0d0
       g5_deriv_x_zz(i,j,2,k) = 0.0d0
       g5_deriv_x_zz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 endif


!------
! Y_SF
!------

!Y-Y
 if( io_x == 0 .and. io_y == 1 .and. io_z == 0 )then


 do j = 1, num_g2_y_x
   do i = 1, natom_y
     g2_y_x(i,j) = ( -2.0d0*g2_y_x_maxmin(j,2)/( g2_y_x_maxmin(j,1) - g2_y_x_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g2_y_y
   do i = 1, natom_y
     g2_y_y(i,j) = g2_y_y_total( counter_natom_y + i, j )
   enddo
 enddo


 do j = 1, num_g2_y_z
   do i = 1, natom_y
     g2_y_z(i,j) = ( -2.0d0*g2_y_z_maxmin(j,2)/( g2_y_z_maxmin(j,1) - g2_y_z_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_y_xx
   do i = 1, natom_y
     g5_y_xx(i,j) = ( -2.0d0*g5_y_xx_maxmin(j,2)/( g5_y_xx_maxmin(j,1) - g5_y_xx_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_y_xy
   do i = 1, natom_y
     g5_y_xy(i,j) = ( -2.0d0*g5_y_xy_maxmin(j,2)/( g5_y_xy_maxmin(j,1) - g5_y_xy_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_y_xz
   do i = 1, natom_y
     g5_y_xz(i,j) = ( -2.0d0*g5_y_xz_maxmin(j,2)/( g5_y_xz_maxmin(j,1) - g5_y_xz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_y_yy
   do i = 1, natom_y
     g5_y_yy(i,j) = g5_y_yy_total( counter_natom_y + i, j )
   enddo
 enddo


 do j = 1, num_g5_y_yz
   do i = 1, natom_y
     g5_y_yz(i,j) = ( -2.0d0*g5_y_yz_maxmin(j,2)/( g5_y_yz_maxmin(j,1) - g5_y_yz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_y_zz
   do i = 1, natom_y
     g5_y_zz(i,j) = ( -2.0d0*g5_y_zz_maxmin(j,2)/( g5_y_zz_maxmin(j,1) - g5_y_zz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 counter_natom_y   = counter_natom_y   + natom_y

 endif

!Y-XY
 if( io_x == 1 .and. io_y == 1 .and. io_z == 0 )then


 do j = 1, num_g2_y_x
   do i = 1, natom_y
     g2_y_x(i,j) = g2_y_x_total( counter_natom_yx + i, j )
   enddo
 enddo


 do j = 1, num_g2_y_y
   do i = 1, natom_y
     g2_y_y(i,j) = g2_y_y_total( counter_natom_y + i, j )
   enddo
 enddo


 do j = 1, num_g2_y_z
   do i = 1, natom_y
     g2_y_z(i,j) = ( -2.0d0*g2_y_z_maxmin(j,2)/( g2_y_z_maxmin(j,1) - g2_y_z_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_y_xx
   do i = 1, natom_y
     g5_y_xx(i,j) = g5_y_xx_total( counter_natom_yx + i, j )
   enddo
 enddo


 do j = 1, num_g5_y_xy
   do i = 1, natom_y
     g5_y_xy(i,j) = g5_y_xy_total( counter_natom_yx + i, j )
   enddo
 enddo


 do j = 1, num_g5_y_xz
   do i = 1, natom_y
     g5_y_xz(i,j) = ( -2.0d0*g5_y_xz_maxmin(j,2)/( g5_y_xz_maxmin(j,1) - g5_y_xz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_y_yy
   do i = 1, natom_y
     g5_y_yy(i,j) = g5_y_yy_total( counter_natom_y + i, j )
   enddo
 enddo


 do j = 1, num_g5_y_yz
   do i = 1, natom_y
     g5_y_yz(i,j) = ( -2.0d0*g5_y_yz_maxmin(j,2)/( g5_y_yz_maxmin(j,1) - g5_y_yz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_y_zz
   do i = 1, natom_y
     g5_y_zz(i,j) = ( -2.0d0*g5_y_zz_maxmin(j,2)/( g5_y_zz_maxmin(j,1) - g5_y_zz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 counter_natom_y   = counter_natom_y   + natom_y
 counter_natom_yx  = counter_natom_yx  + natom_y

 endif

!Y-YZ
 if( io_x == 0 .and. io_y == 1 .and. io_z == 1 )then


 do j = 1, num_g2_y_x
   do i = 1, natom_y
     g2_y_x(i,j) = ( -2.0d0*g2_y_x_maxmin(j,2)/( g2_y_x_maxmin(j,1) - g2_y_x_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g2_y_y
   do i = 1, natom_y
     g2_y_y(i,j) = g2_y_y_total( counter_natom_y + i, j )
   enddo
 enddo


 do j = 1, num_g2_y_z
   do i = 1, natom_y
     g2_y_z(i,j) = g2_y_z_total( counter_natom_yz + i, j )
   enddo
 enddo


 do j = 1, num_g5_y_xx
   do i = 1, natom_y
     g5_y_xx(i,j) = ( -2.0d0*g5_y_xx_maxmin(j,2)/( g5_y_xx_maxmin(j,1) - g5_y_xx_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_y_xy
   do i = 1, natom_y
     g5_y_xy(i,j) = ( -2.0d0*g5_y_xy_maxmin(j,2)/( g5_y_xy_maxmin(j,1) - g5_y_xy_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_y_xz
   do i = 1, natom_y
     g5_y_xz(i,j) = ( -2.0d0*g5_y_xz_maxmin(j,2)/( g5_y_xz_maxmin(j,1) - g5_y_xz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_y_yy
   do i = 1, natom_y
     g5_y_yy(i,j) = g5_y_yy_total( counter_natom_y + i, j )
   enddo
 enddo


 do j = 1, num_g5_y_yz
   do i = 1, natom_y
     g5_y_yz(i,j) = g5_y_yz_total( counter_natom_yz + i, j )
   enddo
 enddo


 do j = 1, num_g5_y_zz
   do i = 1, natom_y
     g5_y_zz(i,j) = g5_y_zz_total( counter_natom_yz + i, j )
   enddo
 enddo


 counter_natom_y   = counter_natom_y   + natom_y
 counter_natom_yz  = counter_natom_yz  + natom_y

 endif

!Y-XYZ
 if( io_x == 1 .and. io_y == 1 .and. io_z == 1 )then


 do j = 1, num_g2_y_x
   do i = 1, natom_y
     g2_y_x(i,j) = g2_y_x_total( counter_natom_yx + i, j )
   enddo
 enddo


 do j = 1, num_g2_y_y
   do i = 1, natom_y
     g2_y_y(i,j) = g2_y_y_total( counter_natom_y + i, j )
   enddo
 enddo


 do j = 1, num_g2_y_z
   do i = 1, natom_y
     g2_y_z(i,j) = g2_y_z_total( counter_natom_yz + i, j )
   enddo
 enddo


 do j = 1, num_g5_y_xx
   do i = 1, natom_y
     g5_y_xx(i,j) = g5_y_xx_total( counter_natom_yx + i, j )
   enddo
 enddo


 do j = 1, num_g5_y_xy
   do i = 1, natom_y
     g5_y_xy(i,j) = g5_y_xy_total( counter_natom_yx + i, j )
   enddo
 enddo


 do j = 1, num_g5_y_xz
   do i = 1, natom_y
     g5_y_xz(i,j) = g5_y_xz_total( counter_natom_yxz + i, j )
   enddo
 enddo


 do j = 1, num_g5_y_yy
   do i = 1, natom_y
     g5_y_yy(i,j) = g5_y_yy_total( counter_natom_y + i, j )
   enddo
 enddo


 do j = 1, num_g5_y_yz
   do i = 1, natom_y
     g5_y_yz(i,j) = g5_y_yz_total( counter_natom_yz + i, j )
   enddo
 enddo


 do j = 1, num_g5_y_zz
   do i = 1, natom_y
     g5_y_zz(i,j) = g5_y_zz_total( counter_natom_yz + i, j )
   enddo
 enddo


 counter_natom_y   = counter_natom_y   + natom_y
 counter_natom_yx  = counter_natom_yx  + natom_y
 counter_natom_yz  = counter_natom_yz  + natom_y
 counter_natom_yxz = counter_natom_yxz + natom_y

 endif

!Y-X, Y-Z, Y-XZ
 if( io_y == 0 )then


 do j = 1, num_g2_y_x
   do i = 1, natom_y
     g2_y_x(i,j) = ( -2.0d0*g2_y_x_maxmin(j,2)/( g2_y_x_maxmin(j,1) - g2_y_x_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g2_y_y
   do i = 1, natom_y
     g2_y_y(i,j) = ( -2.0d0*g2_y_y_maxmin(j,2)/( g2_y_y_maxmin(j,1) - g2_y_y_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g2_y_z
   do i = 1, natom_y
     g2_y_z(i,j) = ( -2.0d0*g2_y_z_maxmin(j,2)/( g2_y_z_maxmin(j,1) - g2_y_z_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_y_xx
   do i = 1, natom_y
     g5_y_xx(i,j) = ( -2.0d0*g5_y_xx_maxmin(j,2)/( g5_y_xx_maxmin(j,1) - g5_y_xx_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_y_xy
   do i = 1, natom_y
     g5_y_xy(i,j) = ( -2.0d0*g5_y_xy_maxmin(j,2)/( g5_y_xy_maxmin(j,1) - g5_y_xy_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_y_xz
   do i = 1, natom_y
     g5_y_xz(i,j) = ( -2.0d0*g5_y_xz_maxmin(j,2)/( g5_y_xz_maxmin(j,1) - g5_y_xz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_y_yy
   do i = 1, natom_y
     g5_y_yy(i,j) = ( -2.0d0*g5_y_yy_maxmin(j,2)/( g5_y_yy_maxmin(j,1) - g5_y_yy_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_y_yz
   do i = 1, natom_y
     g5_y_yz(i,j) = ( -2.0d0*g5_y_yz_maxmin(j,2)/( g5_y_yz_maxmin(j,1) - g5_y_yz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_y_zz
   do i = 1, natom_y
     g5_y_zz(i,j) = ( -2.0d0*g5_y_zz_maxmin(j,2)/( g5_y_zz_maxmin(j,1) - g5_y_zz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 endif

!------------
! Y_SF_deriv
!------------

!Y-Y
 if( io_x == 0 .and. io_y == 1 .and. io_z == 0 )then


 do k = 1, num_g2_y_x
   do j = 1, natom
     do i = 1, natom_y
       g2_deriv_y_x(i,j,1,k) = 0.0d0
       g2_deriv_y_x(i,j,2,k) = 0.0d0
       g2_deriv_y_x(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g2_y_y
   do j = 1, natom
     do i = 1, natom_y
       m = i + natom_y * ( j - 1 )
       g2_deriv_y_y(i,j,1,k) = g2_deriv_y_y_total( counter_deriv_y + m, 1, k )
       g2_deriv_y_y(i,j,2,k) = g2_deriv_y_y_total( counter_deriv_y + m, 2, k )
       g2_deriv_y_y(i,j,3,k) = g2_deriv_y_y_total( counter_deriv_y + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g2_y_z
   do j = 1, natom
     do i = 1, natom_y
       g2_deriv_y_z(i,j,1,k) = 0.0d0
       g2_deriv_y_z(i,j,2,k) = 0.0d0
       g2_deriv_y_z(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_xx
   do j = 1, natom
     do i = 1, natom_y
       g5_deriv_y_xx(i,j,1,k) = 0.0d0
       g5_deriv_y_xx(i,j,2,k) = 0.0d0
       g5_deriv_y_xx(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_xy
   do j = 1, natom
     do i = 1, natom_y
       g5_deriv_y_xy(i,j,1,k) = 0.0d0
       g5_deriv_y_xy(i,j,2,k) = 0.0d0
       g5_deriv_y_xy(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_xz
   do j = 1, natom
     do i = 1, natom_y
       g5_deriv_y_xz(i,j,1,k) = 0.0d0
       g5_deriv_y_xz(i,j,2,k) = 0.0d0
       g5_deriv_y_xz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_yy
   do j = 1, natom
     do i = 1, natom_y
       m = i + natom_y * ( j - 1 )
       g5_deriv_y_yy(i,j,1,k) = g5_deriv_y_yy_total( counter_deriv_y + m, 1, k )
       g5_deriv_y_yy(i,j,2,k) = g5_deriv_y_yy_total( counter_deriv_y + m, 2, k )
       g5_deriv_y_yy(i,j,3,k) = g5_deriv_y_yy_total( counter_deriv_y + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_yz
   do j = 1, natom
     do i = 1, natom_y
       g5_deriv_y_yz(i,j,1,k) = 0.0d0
       g5_deriv_y_yz(i,j,2,k) = 0.0d0
       g5_deriv_y_yz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_zz
   do j = 1, natom
     do i = 1, natom_y
       g5_deriv_y_zz(i,j,1,k) = 0.0d0
       g5_deriv_y_zz(i,j,2,k) = 0.0d0
       g5_deriv_y_zz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 counter_deriv_y   = counter_deriv_y   + natom_y * natom

 endif

!Y-XY
 if( io_x == 1 .and. io_y == 1 .and. io_z == 0 )then


 do k = 1, num_g2_y_x
   do j = 1, natom
     do i = 1, natom_y
       m = i + natom_y * ( j - 1 )
       g2_deriv_y_x(i,j,1,k) = g2_deriv_y_x_total( counter_deriv_yx + m, 1, k )
       g2_deriv_y_x(i,j,2,k) = g2_deriv_y_x_total( counter_deriv_yx + m, 2, k )
       g2_deriv_y_x(i,j,3,k) = g2_deriv_y_x_total( counter_deriv_yx + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g2_y_y
   do j = 1, natom
     do i = 1, natom_y
       m = i + natom_y * ( j - 1 )
       g2_deriv_y_y(i,j,1,k) = g2_deriv_y_y_total( counter_deriv_y + m, 1, k )
       g2_deriv_y_y(i,j,2,k) = g2_deriv_y_y_total( counter_deriv_y + m, 2, k )
       g2_deriv_y_y(i,j,3,k) = g2_deriv_y_y_total( counter_deriv_y + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g2_y_z
   do j = 1, natom
     do i = 1, natom_y
       g2_deriv_y_z(i,j,1,k) = 0.0d0
       g2_deriv_y_z(i,j,2,k) = 0.0d0
       g2_deriv_y_z(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_xx
   do j = 1, natom
     do i = 1, natom_y
       m = i + natom_y * ( j - 1 )
       g5_deriv_y_xx(i,j,1,k) = g5_deriv_y_xx_total( counter_deriv_yx + m, 1, k )
       g5_deriv_y_xx(i,j,2,k) = g5_deriv_y_xx_total( counter_deriv_yx + m, 2, k )
       g5_deriv_y_xx(i,j,3,k) = g5_deriv_y_xx_total( counter_deriv_yx + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_xy
   do j = 1, natom
     do i = 1, natom_y
       m = i + natom_y * ( j - 1 )
       g5_deriv_y_xy(i,j,1,k) = g5_deriv_y_xy_total( counter_deriv_yx + m, 1, k )
       g5_deriv_y_xy(i,j,2,k) = g5_deriv_y_xy_total( counter_deriv_yx + m, 2, k )
       g5_deriv_y_xy(i,j,3,k) = g5_deriv_y_xy_total( counter_deriv_yx + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_xz
   do j = 1, natom
     do i = 1, natom_y
       g5_deriv_y_xz(i,j,1,k) = 0.0d0
       g5_deriv_y_xz(i,j,2,k) = 0.0d0
       g5_deriv_y_xz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_yy
   do j = 1, natom
     do i = 1, natom_y
       m = i + natom_y * ( j - 1 )
       g5_deriv_y_yy(i,j,1,k) = g5_deriv_y_yy_total( counter_deriv_y + m, 1, k )
       g5_deriv_y_yy(i,j,2,k) = g5_deriv_y_yy_total( counter_deriv_y + m, 2, k )
       g5_deriv_y_yy(i,j,3,k) = g5_deriv_y_yy_total( counter_deriv_y + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_yz
   do j = 1, natom
     do i = 1, natom_y
       g5_deriv_y_yz(i,j,1,k) = 0.0d0
       g5_deriv_y_yz(i,j,2,k) = 0.0d0
       g5_deriv_y_yz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_zz
   do j = 1, natom
     do i = 1, natom_y
       g5_deriv_y_zz(i,j,1,k) = 0.0d0
       g5_deriv_y_zz(i,j,2,k) = 0.0d0
       g5_deriv_y_zz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 counter_deriv_y   = counter_deriv_y   + natom_y * natom
 counter_deriv_yx  = counter_deriv_yx  + natom_y * natom

 endif

!Y-YZ
 if( io_x == 0 .and. io_y == 1 .and. io_z == 1 )then


 do k = 1, num_g2_y_x
   do j = 1, natom
     do i = 1, natom_y
       g2_deriv_y_x(i,j,1,k) = 0.0d0
       g2_deriv_y_x(i,j,2,k) = 0.0d0
       g2_deriv_y_x(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g2_y_y
   do j = 1, natom
     do i = 1, natom_y
       m = i + natom_y * ( j - 1 )
       g2_deriv_y_y(i,j,1,k) = g2_deriv_y_y_total( counter_deriv_y + m, 1, k )
       g2_deriv_y_y(i,j,2,k) = g2_deriv_y_y_total( counter_deriv_y + m, 2, k )
       g2_deriv_y_y(i,j,3,k) = g2_deriv_y_y_total( counter_deriv_y + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g2_y_z
   do j = 1, natom
     do i = 1, natom_y
       m = i + natom_y * ( j - 1 )
       g2_deriv_y_z(i,j,1,k) = g2_deriv_y_z_total( counter_deriv_yz + m, 1, k )
       g2_deriv_y_z(i,j,2,k) = g2_deriv_y_z_total( counter_deriv_yz + m, 2, k )
       g2_deriv_y_z(i,j,3,k) = g2_deriv_y_z_total( counter_deriv_yz + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_xx
   do j = 1, natom
     do i = 1, natom_y
       g5_deriv_y_xx(i,j,1,k) = 0.0d0
       g5_deriv_y_xx(i,j,2,k) = 0.0d0
       g5_deriv_y_xx(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_xy
   do j = 1, natom
     do i = 1, natom_y
       g5_deriv_y_xy(i,j,1,k) = 0.0d0
       g5_deriv_y_xy(i,j,2,k) = 0.0d0
       g5_deriv_y_xy(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_xz
   do j = 1, natom
     do i = 1, natom_y
       g5_deriv_y_xz(i,j,1,k) = 0.0d0
       g5_deriv_y_xz(i,j,2,k) = 0.0d0
       g5_deriv_y_xz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_yy
   do j = 1, natom
     do i = 1, natom_y
       m = i + natom_y * ( j - 1 )
       g5_deriv_y_yy(i,j,1,k) = g5_deriv_y_yy_total( counter_deriv_y + m, 1, k )
       g5_deriv_y_yy(i,j,2,k) = g5_deriv_y_yy_total( counter_deriv_y + m, 2, k )
       g5_deriv_y_yy(i,j,3,k) = g5_deriv_y_yy_total( counter_deriv_y + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_yz
   do j = 1, natom
     do i = 1, natom_y
       m = i + natom_y * ( j - 1 )
       g5_deriv_y_yz(i,j,1,k) = g5_deriv_y_yz_total( counter_deriv_yz + m, 1, k )
       g5_deriv_y_yz(i,j,2,k) = g5_deriv_y_yz_total( counter_deriv_yz + m, 2, k )
       g5_deriv_y_yz(i,j,3,k) = g5_deriv_y_yz_total( counter_deriv_yz + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_zz
   do j = 1, natom
     do i = 1, natom_y
       m = i + natom_y * ( j - 1 )
       g5_deriv_y_zz(i,j,1,k) = g5_deriv_y_zz_total( counter_deriv_yz + m, 1, k )
       g5_deriv_y_zz(i,j,2,k) = g5_deriv_y_zz_total( counter_deriv_yz + m, 2, k )
       g5_deriv_y_zz(i,j,3,k) = g5_deriv_y_zz_total( counter_deriv_yz + m, 3, k )
     enddo
   enddo
 enddo


 counter_deriv_y   = counter_deriv_y   + natom_y * natom
 counter_deriv_yz  = counter_deriv_yz  + natom_y * natom

 endif

!Y-XYZ
 if( io_x == 1 .and. io_y == 1 .and. io_z == 1 )then


 do k = 1, num_g2_y_x
   do j = 1, natom
     do i = 1, natom_y
       m = i + natom_y * ( j - 1 )
       g2_deriv_y_x(i,j,1,k) = g2_deriv_y_x_total( counter_deriv_yx + m, 1, k )
       g2_deriv_y_x(i,j,2,k) = g2_deriv_y_x_total( counter_deriv_yx + m, 2, k )
       g2_deriv_y_x(i,j,3,k) = g2_deriv_y_x_total( counter_deriv_yx + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g2_y_y
   do j = 1, natom
     do i = 1, natom_y
       m = i + natom_y * ( j - 1 )
       g2_deriv_y_y(i,j,1,k) = g2_deriv_y_y_total( counter_deriv_y + m, 1, k )
       g2_deriv_y_y(i,j,2,k) = g2_deriv_y_y_total( counter_deriv_y + m, 2, k )
       g2_deriv_y_y(i,j,3,k) = g2_deriv_y_y_total( counter_deriv_y + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g2_y_z
   do j = 1, natom
     do i = 1, natom_y
       m = i + natom_y * ( j - 1 )
       g2_deriv_y_z(i,j,1,k) = g2_deriv_y_z_total( counter_deriv_yz + m, 1, k )
       g2_deriv_y_z(i,j,2,k) = g2_deriv_y_z_total( counter_deriv_yz + m, 2, k )
       g2_deriv_y_z(i,j,3,k) = g2_deriv_y_z_total( counter_deriv_yz + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_xx
   do j = 1, natom
     do i = 1, natom_y
       m = i + natom_y * ( j - 1 )
       g5_deriv_y_xx(i,j,1,k) = g5_deriv_y_xx_total( counter_deriv_yx + m, 1, k )
       g5_deriv_y_xx(i,j,2,k) = g5_deriv_y_xx_total( counter_deriv_yx + m, 2, k )
       g5_deriv_y_xx(i,j,3,k) = g5_deriv_y_xx_total( counter_deriv_yx + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_xy
   do j = 1, natom
     do i = 1, natom_y
       m = i + natom_y * ( j - 1 )
       g5_deriv_y_xy(i,j,1,k) = g5_deriv_y_xy_total( counter_deriv_yx + m, 1, k )
       g5_deriv_y_xy(i,j,2,k) = g5_deriv_y_xy_total( counter_deriv_yx + m, 2, k )
       g5_deriv_y_xy(i,j,3,k) = g5_deriv_y_xy_total( counter_deriv_yx + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_xz
   do j = 1, natom
     do i = 1, natom_y
       m = i + natom_y * ( j - 1 )
       g5_deriv_y_xz(i,j,1,k) = g5_deriv_y_xz_total( counter_deriv_yxz + m, 1, k )
       g5_deriv_y_xz(i,j,2,k) = g5_deriv_y_xz_total( counter_deriv_yxz + m, 2, k )
       g5_deriv_y_xz(i,j,3,k) = g5_deriv_y_xz_total( counter_deriv_yxz + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_yy
   do j = 1, natom
     do i = 1, natom_y
       m = i + natom_y * ( j - 1 )
       g5_deriv_y_yy(i,j,1,k) = g5_deriv_y_yy_total( counter_deriv_y + m, 1, k )
       g5_deriv_y_yy(i,j,2,k) = g5_deriv_y_yy_total( counter_deriv_y + m, 2, k )
       g5_deriv_y_yy(i,j,3,k) = g5_deriv_y_yy_total( counter_deriv_y + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_yz
   do j = 1, natom
     do i = 1, natom_y
       m = i + natom_y * ( j - 1 )
       g5_deriv_y_yz(i,j,1,k) = g5_deriv_y_yz_total( counter_deriv_yz + m, 1, k )
       g5_deriv_y_yz(i,j,2,k) = g5_deriv_y_yz_total( counter_deriv_yz + m, 2, k )
       g5_deriv_y_yz(i,j,3,k) = g5_deriv_y_yz_total( counter_deriv_yz + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_zz
   do j = 1, natom
     do i = 1, natom_y
       m = i + natom_y * ( j - 1 )
       g5_deriv_y_zz(i,j,1,k) = g5_deriv_y_zz_total( counter_deriv_yz + m, 1, k )
       g5_deriv_y_zz(i,j,2,k) = g5_deriv_y_zz_total( counter_deriv_yz + m, 2, k )
       g5_deriv_y_zz(i,j,3,k) = g5_deriv_y_zz_total( counter_deriv_yz + m, 3, k )
     enddo
   enddo
 enddo


 counter_deriv_y   = counter_deriv_y   + natom_y * natom
 counter_deriv_yx  = counter_deriv_yx  + natom_y * natom
 counter_deriv_yz  = counter_deriv_yz  + natom_y * natom
 counter_deriv_yxz = counter_deriv_yxz + natom_y * natom 

 endif

!Y-X, Y-Z, Y-XZ
 if( io_y == 0 )then


 do k = 1, num_g2_y_x
   do j = 1, natom
     do i = 1, natom_y
       g2_deriv_y_x(i,j,1,k) = 0.0d0
       g2_deriv_y_x(i,j,2,k) = 0.0d0
       g2_deriv_y_x(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g2_y_y
   do j = 1, natom
     do i = 1, natom_y
       g2_deriv_y_y(i,j,1,k) = 0.0d0
       g2_deriv_y_y(i,j,2,k) = 0.0d0
       g2_deriv_y_y(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g2_y_z
   do j = 1, natom
     do i = 1, natom_y
       g2_deriv_y_z(i,j,1,k) = 0.0d0
       g2_deriv_y_z(i,j,2,k) = 0.0d0
       g2_deriv_y_z(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_xx
   do j = 1, natom
     do i = 1, natom_y
       g5_deriv_y_xx(i,j,1,k) = 0.0d0
       g5_deriv_y_xx(i,j,2,k) = 0.0d0
       g5_deriv_y_xx(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_xy
   do j = 1, natom
     do i = 1, natom_y
       g5_deriv_y_xy(i,j,1,k) = 0.0d0
       g5_deriv_y_xy(i,j,2,k) = 0.0d0
       g5_deriv_y_xy(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_xz
   do j = 1, natom
     do i = 1, natom_y
       g5_deriv_y_xz(i,j,1,k) = 0.0d0
       g5_deriv_y_xz(i,j,2,k) = 0.0d0
       g5_deriv_y_xz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_yy
   do j = 1, natom
     do i = 1, natom_y
       g5_deriv_y_yy(i,j,1,k) = 0.0d0
       g5_deriv_y_yy(i,j,2,k) = 0.0d0
       g5_deriv_y_yy(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_yz
   do j = 1, natom
     do i = 1, natom_y
       g5_deriv_y_yz(i,j,1,k) = 0.0d0
       g5_deriv_y_yz(i,j,2,k) = 0.0d0
       g5_deriv_y_yz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_y_zz
   do j = 1, natom
     do i = 1, natom_y
       g5_deriv_y_zz(i,j,1,k) = 0.0d0
       g5_deriv_y_zz(i,j,2,k) = 0.0d0
       g5_deriv_y_zz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 endif


!------
! Z_SF
!------

!Z-Z
 if( io_x == 0 .and. io_y == 0 .and. io_z == 1 )then


 do j = 1, num_g2_z_x
   do i = 1, natom_z
     g2_z_x(i,j) = ( -2.0d0*g2_z_x_maxmin(j,2)/( g2_z_x_maxmin(j,1) - g2_z_x_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g2_z_y
   do i = 1, natom_z
     g2_z_y(i,j) = ( -2.0d0*g2_z_y_maxmin(j,2)/( g2_z_y_maxmin(j,1) - g2_z_y_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g2_z_z
   do i = 1, natom_z
     g2_z_z(i,j) = g2_z_z_total( counter_natom_z + i, j )
   enddo
 enddo


 do j = 1, num_g5_z_xx
   do i = 1, natom_z
     g5_z_xx(i,j) = ( - 2.0d0*g5_z_xx_maxmin(j,2)/( g5_z_xx_maxmin(j,1) - g5_z_xx_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_z_xy
   do i = 1, natom_z
     g5_z_xy(i,j) = ( -2.0d0*g5_z_xy_maxmin(j,2)/( g5_z_xy_maxmin(j,1) - g5_z_xy_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_z_xz
   do i = 1, natom_z
     g5_z_xz(i,j) = ( -2.0d0*g5_z_xz_maxmin(j,2)/( g5_z_xz_maxmin(j,1) - g5_z_xz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_z_yy
   do i = 1, natom_z
     g5_z_yy(i,j) = ( -2.0d0*g5_z_yy_maxmin(j,2)/( g5_z_yy_maxmin(j,1) - g5_z_yy_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_z_yz
   do i = 1, natom_z
     g5_z_yz(i,j) = ( -2.0d0*g5_z_yz_maxmin(j,2)/( g5_z_yz_maxmin(j,1) - g5_z_yz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_z_zz
   do i = 1, natom_z
     g5_z_zz(i,j) = g5_z_zz_total( counter_natom_z + i, j )
   enddo
 enddo


 counter_natom_z   = counter_natom_z   + natom_z

 endif

!Z-XZ
 if( io_x == 1 .and. io_y == 0 .and. io_z == 1 )then


 do j = 1, num_g2_z_x
   do i = 1, natom_z
     g2_z_x(i,j) = g2_z_x_total( counter_natom_zx + i, j )
   enddo
 enddo


 do j = 1, num_g2_z_y
   do i = 1, natom_z
     g2_z_y(i,j) = ( -2.0d0*g2_z_y_maxmin(j,2)/( g2_z_y_maxmin(j,1) - g2_z_y_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g2_z_z
   do i = 1, natom_z
     g2_z_z(i,j) = g2_z_z_total( counter_natom_z + i, j )
   enddo
 enddo


 do j = 1, num_g5_z_xx
   do i = 1, natom_z
     g5_z_xx(i,j) = g5_z_xx_total( counter_natom_zx + i, j )
   enddo
 enddo


 do j = 1, num_g5_z_xy
   do i = 1, natom_z
     g5_z_xy(i,j) = ( -2.0d0*g5_z_xy_maxmin(j,2)/( g5_z_xy_maxmin(j,1) - g5_z_xy_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_z_xz
   do i = 1, natom_z
     g5_z_xz(i,j) = g5_z_xz_total( counter_natom_zx + i, j )
   enddo
 enddo


 do j = 1, num_g5_z_yy
   do i = 1, natom_z
     g5_z_yy(i,j) = ( -2.0d0*g5_z_yy_maxmin(j,2)/( g5_z_yy_maxmin(j,1) - g5_z_yy_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_z_yz
   do i = 1, natom_z
     g5_z_yz(i,j) = ( -2.0d0*g5_z_yz_maxmin(j,2)/( g5_z_yz_maxmin(j,1) - g5_z_yz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_z_zz
   do i = 1, natom_z
     g5_z_zz(i,j) = g5_z_zz_total( counter_natom_z + i, j )
   enddo
 enddo


 counter_natom_z   = counter_natom_z   + natom_z
 counter_natom_zx  = counter_natom_zx  + natom_z

 endif

!Z-YZ
 if( io_x == 0 .and. io_y == 1 .and. io_z == 1 )then


 do j = 1, num_g2_z_x
   do i = 1, natom_z
     g2_z_x(i,j) = ( -2.0d0*g2_z_x_maxmin(j,2)/( g2_z_x_maxmin(j,1) - g2_z_x_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g2_z_y
   do i = 1, natom_z
     g2_z_y(i,j) = g2_z_y_total( counter_natom_zy + i, j )
   enddo
 enddo


 do j = 1, num_g2_z_z
   do i = 1, natom_z
     g2_z_z(i,j) = g2_z_z_total( counter_natom_z + i, j )
   enddo
 enddo


 do j = 1, num_g5_z_xx
   do i = 1, natom_z
     g5_z_xx(i,j) = ( - 2.0d0*g5_z_xx_maxmin(j,2)/( g5_z_xx_maxmin(j,1) - g5_z_xx_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_z_xy
   do i = 1, natom_z
     g5_z_xy(i,j) = ( -2.0d0*g5_z_xy_maxmin(j,2)/( g5_z_xy_maxmin(j,1) - g5_z_xy_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_z_xz
   do i = 1, natom_z
     g5_z_xz(i,j) = ( -2.0d0*g5_z_xz_maxmin(j,2)/( g5_z_xz_maxmin(j,1) - g5_z_xz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_z_yy
   do i = 1, natom_z
     g5_z_yy(i,j) = g5_z_yy_total( counter_natom_zy + i, j )
   enddo
 enddo


 do j = 1, num_g5_z_yz
   do i = 1, natom_z
     g5_z_yz(i,j) = g5_z_yz_total( counter_natom_zy + i, j )
   enddo
 enddo


 do j = 1, num_g5_z_zz
   do i = 1, natom_z
     g5_z_zz(i,j) = g5_z_zz_total( counter_natom_z + i, j )
   enddo
 enddo


 counter_natom_z   = counter_natom_z   + natom_z
 counter_natom_zy  = counter_natom_zy  + natom_z

 endif

!Z-XYZ
 if( io_x == 1 .and. io_y == 1 .and. io_z == 1 )then


 do j = 1, num_g2_z_x
   do i = 1, natom_z
     g2_z_x(i,j) = g2_z_x_total( counter_natom_zx + i, j )
   enddo
 enddo


 do j = 1, num_g2_z_y
   do i = 1, natom_z
     g2_z_y(i,j) = g2_z_y_total( counter_natom_zy + i, j )
   enddo
 enddo


 do j = 1, num_g2_z_z
   do i = 1, natom_z
     g2_z_z(i,j) = g2_z_z_total( counter_natom_z + i, j )
   enddo
 enddo


 do j = 1, num_g5_z_xx
   do i = 1, natom_z
     g5_z_xx(i,j) = g5_z_xx_total( counter_natom_zx + i, j )
   enddo
 enddo


 do j = 1, num_g5_z_xy
   do i = 1, natom_z
     g5_z_xy(i,j) = g5_z_xy_total( counter_natom_zxy + i, j )
   enddo
 enddo


 do j = 1, num_g5_z_xz
   do i = 1, natom_z
     g5_z_xz(i,j) = g5_z_xz_total( counter_natom_zx + i, j )
   enddo
 enddo


 do j = 1, num_g5_z_yy
   do i = 1, natom_z
     g5_z_yy(i,j) = g5_z_yy_total( counter_natom_zy + i, j )
   enddo
 enddo


 do j = 1, num_g5_z_yz
   do i = 1, natom_z
     g5_z_yz(i,j) = g5_z_yz_total( counter_natom_zy + i, j )
   enddo
 enddo


 do j = 1, num_g5_z_zz
   do i = 1, natom_z
     g5_z_zz(i,j) = g5_z_zz_total( counter_natom_z + i, j )
   enddo
 enddo


 counter_natom_z   = counter_natom_z   + natom_z
 counter_natom_zx  = counter_natom_zx  + natom_z
 counter_natom_zy  = counter_natom_zy  + natom_z
 counter_natom_zxy = counter_natom_zxy + natom_z

 endif

!Z-X, Z-Y, Z-XY
 if( io_z == 0 )then


 do j = 1, num_g2_z_x
   do i = 1, natom_z
     g2_z_x(i,j) = ( -2.0d0*g2_z_x_maxmin(j,2)/( g2_z_x_maxmin(j,1) - g2_z_x_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g2_z_y
   do i = 1, natom_z
     g2_z_y(i,j) = ( -2.0d0*g2_z_y_maxmin(j,2)/( g2_z_y_maxmin(j,1) - g2_z_y_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g2_z_z
   do i = 1, natom_z
     g2_z_z(i,j) = ( -2.0d0*g2_z_z_maxmin(j,2)/( g2_z_z_maxmin(j,1) - g2_z_z_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_z_xx
   do i = 1, natom_z
     g5_z_xx(i,j) = ( - 2.0d0*g5_z_xx_maxmin(j,2)/( g5_z_xx_maxmin(j,1) - g5_z_xx_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_z_xy
   do i = 1, natom_z
     g5_z_xy(i,j) = ( -2.0d0*g5_z_xy_maxmin(j,2)/( g5_z_xy_maxmin(j,1) - g5_z_xy_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_z_xz
   do i = 1, natom_z
     g5_z_xz(i,j) = ( -2.0d0*g5_z_xz_maxmin(j,2)/( g5_z_xz_maxmin(j,1) - g5_z_xz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_z_yy
   do i = 1, natom_z
     g5_z_yy(i,j) = ( -2.0d0*g5_z_yy_maxmin(j,2)/( g5_z_yy_maxmin(j,1) - g5_z_yy_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_z_yz
   do i = 1, natom_z
     g5_z_yz(i,j) = ( -2.0d0*g5_z_yz_maxmin(j,2)/( g5_z_yz_maxmin(j,1) - g5_z_yz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 do j = 1, num_g5_z_zz
   do i = 1, natom_z
     g5_z_zz(i,j) = ( -2.0d0*g5_z_zz_maxmin(j,2)/( g5_z_zz_maxmin(j,1) - g5_z_zz_maxmin(j,2) ) ) - 1.0d0
   enddo
 enddo


 endif

!------------
! Z_SF_deriv
!------------

!Z-Z
 if( io_x == 0 .and. io_y == 0 .and. io_z == 1 )then


 do k = 1, num_g2_z_x
   do j = 1, natom
     do i = 1, natom_z
       g2_deriv_z_x(i,j,1,k) = 0.0d0
       g2_deriv_z_x(i,j,2,k) = 0.0d0
       g2_deriv_z_x(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g2_z_y
   do j = 1, natom
     do i = 1, natom_z
       g2_deriv_z_y(i,j,1,k) = 0.0d0
       g2_deriv_z_y(i,j,2,k) = 0.0d0
       g2_deriv_z_y(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g2_z_z
   do j = 1, natom
     do i = 1, natom_z
       m = i + natom_z * ( j - 1 )
       g2_deriv_z_z(i,j,1,k) = g2_deriv_z_z_total( counter_deriv_z + m, 1, k )
       g2_deriv_z_z(i,j,2,k) = g2_deriv_z_z_total( counter_deriv_z + m, 2, k )
       g2_deriv_z_z(i,j,3,k) = g2_deriv_z_z_total( counter_deriv_z + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_xx
   do j = 1, natom
     do i = 1, natom_z
       g5_deriv_z_xx(i,j,1,k) = 0.0d0
       g5_deriv_z_xx(i,j,2,k) = 0.0d0
       g5_deriv_z_xx(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_xy
   do j = 1, natom
     do i = 1, natom_z
       g5_deriv_z_xy(i,j,1,k) = 0.0d0
       g5_deriv_z_xy(i,j,2,k) = 0.0d0
       g5_deriv_z_xy(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_xz
   do j = 1, natom
     do i = 1, natom_z
       g5_deriv_z_xz(i,j,1,k) = 0.0d0
       g5_deriv_z_xz(i,j,2,k) = 0.0d0
       g5_deriv_z_xz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_yy
   do j = 1, natom
     do i = 1, natom_z
       g5_deriv_z_yy(i,j,1,k) = 0.0d0
       g5_deriv_z_yy(i,j,2,k) = 0.0d0
       g5_deriv_z_yy(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_yz
   do j = 1, natom
     do i = 1, natom_z
       g5_deriv_z_yz(i,j,1,k) = 0.0d0
       g5_deriv_z_yz(i,j,2,k) = 0.0d0
       g5_deriv_z_yz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_zz
   do j = 1, natom
     do i = 1, natom_z
       m = i + natom_z * ( j - 1 )
       g5_deriv_z_zz(i,j,1,k) = g5_deriv_z_zz_total( counter_deriv_z + m, 1, k )
       g5_deriv_z_zz(i,j,2,k) = g5_deriv_z_zz_total( counter_deriv_z + m, 2, k )
       g5_deriv_z_zz(i,j,3,k) = g5_deriv_z_zz_total( counter_deriv_z + m, 3, k )
     enddo
   enddo
 enddo


 counter_deriv_z   = counter_deriv_z   + natom_z * natom

 endif ! io_z

!Z-XZ
 if( io_x == 1 .and. io_y == 0 .and. io_z == 1 )then


 do k = 1, num_g2_z_x
   do j = 1, natom
     do i = 1, natom_z
       m = i + natom_z * ( j - 1 )
       g2_deriv_z_x(i,j,1,k) = g2_deriv_z_x_total( counter_deriv_zx + m, 1, k )
       g2_deriv_z_x(i,j,2,k) = g2_deriv_z_x_total( counter_deriv_zx + m, 2, k )
       g2_deriv_z_x(i,j,3,k) = g2_deriv_z_x_total( counter_deriv_zx + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g2_z_y
   do j = 1, natom
     do i = 1, natom_z
       g2_deriv_z_y(i,j,1,k) = 0.0d0
       g2_deriv_z_y(i,j,2,k) = 0.0d0
       g2_deriv_z_y(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g2_z_z
   do j = 1, natom
     do i = 1, natom_z
       m = i + natom_z * ( j - 1 )
       g2_deriv_z_z(i,j,1,k) = g2_deriv_z_z_total( counter_deriv_z + m, 1, k )
       g2_deriv_z_z(i,j,2,k) = g2_deriv_z_z_total( counter_deriv_z + m, 2, k )
       g2_deriv_z_z(i,j,3,k) = g2_deriv_z_z_total( counter_deriv_z + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_xx
   do j = 1, natom
     do i = 1, natom_z
       m = i + natom_z * ( j - 1 )
       g5_deriv_z_xx(i,j,1,k) = g5_deriv_z_xx_total( counter_deriv_zx + m, 1, k )
       g5_deriv_z_xx(i,j,2,k) = g5_deriv_z_xx_total( counter_deriv_zx + m, 2, k )
       g5_deriv_z_xx(i,j,3,k) = g5_deriv_z_xx_total( counter_deriv_zx + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_xy
   do j = 1, natom
     do i = 1, natom_z
       g5_deriv_z_xy(i,j,1,k) = 0.0d0
       g5_deriv_z_xy(i,j,2,k) = 0.0d0
       g5_deriv_z_xy(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_xz
   do j = 1, natom
     do i = 1, natom_z
       m = i + natom_z * ( j - 1 )
       g5_deriv_z_xz(i,j,1,k) = g5_deriv_z_xz_total( counter_deriv_zx + m, 1, k )
       g5_deriv_z_xz(i,j,2,k) = g5_deriv_z_xz_total( counter_deriv_zx + m, 2, k )
       g5_deriv_z_xz(i,j,3,k) = g5_deriv_z_xz_total( counter_deriv_zx + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_yy
   do j = 1, natom
     do i = 1, natom_z
       g5_deriv_z_yy(i,j,1,k) = 0.0d0
       g5_deriv_z_yy(i,j,2,k) = 0.0d0
       g5_deriv_z_yy(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_yz
   do j = 1, natom
     do i = 1, natom_z
       g5_deriv_z_yz(i,j,1,k) = 0.0d0
       g5_deriv_z_yz(i,j,2,k) = 0.0d0
       g5_deriv_z_yz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_zz
   do j = 1, natom
     do i = 1, natom_z
       m = i + natom_z * ( j - 1 )
       g5_deriv_z_zz(i,j,1,k) = g5_deriv_z_zz_total( counter_deriv_z + m, 1, k )
       g5_deriv_z_zz(i,j,2,k) = g5_deriv_z_zz_total( counter_deriv_z + m, 2, k )
       g5_deriv_z_zz(i,j,3,k) = g5_deriv_z_zz_total( counter_deriv_z + m, 3, k )
     enddo
   enddo
 enddo


 counter_deriv_z   = counter_deriv_z   + natom_z * natom
 counter_deriv_zx  = counter_deriv_zx  + natom_z * natom

 endif ! io_z

!Z-YZ
 if( io_x == 0 .and. io_y == 1 .and. io_z == 1 )then


 do k = 1, num_g2_z_x
   do j = 1, natom
     do i = 1, natom_z
       g2_deriv_z_x(i,j,1,k) = 0.0d0
       g2_deriv_z_x(i,j,2,k) = 0.0d0
       g2_deriv_z_x(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g2_z_y
   do j = 1, natom
     do i = 1, natom_z
       m = i + natom_z * ( j - 1 )
       g2_deriv_z_y(i,j,1,k) = g2_deriv_z_y_total( counter_deriv_zy + m, 1, k )
       g2_deriv_z_y(i,j,2,k) = g2_deriv_z_y_total( counter_deriv_zy + m, 2, k )
       g2_deriv_z_y(i,j,3,k) = g2_deriv_z_y_total( counter_deriv_zy + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g2_z_z
   do j = 1, natom
     do i = 1, natom_z
       m = i + natom_z * ( j - 1 )
       g2_deriv_z_z(i,j,1,k) = g2_deriv_z_z_total( counter_deriv_z + m, 1, k )
       g2_deriv_z_z(i,j,2,k) = g2_deriv_z_z_total( counter_deriv_z + m, 2, k )
       g2_deriv_z_z(i,j,3,k) = g2_deriv_z_z_total( counter_deriv_z + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_xx
   do j = 1, natom
     do i = 1, natom_z
       g5_deriv_z_xx(i,j,1,k) = 0.0d0
       g5_deriv_z_xx(i,j,2,k) = 0.0d0
       g5_deriv_z_xx(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_xy
   do j = 1, natom
     do i = 1, natom_z
       g5_deriv_z_xy(i,j,1,k) = 0.0d0
       g5_deriv_z_xy(i,j,2,k) = 0.0d0
       g5_deriv_z_xy(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_xz
   do j = 1, natom
     do i = 1, natom_z
       g5_deriv_z_xz(i,j,1,k) = 0.0d0
       g5_deriv_z_xz(i,j,2,k) = 0.0d0
       g5_deriv_z_xz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_yy
   do j = 1, natom
     do i = 1, natom_z
       m = i + natom_z * ( j - 1 )
       g5_deriv_z_yy(i,j,1,k) = g5_deriv_z_yy_total( counter_deriv_zy + m, 1, k )
       g5_deriv_z_yy(i,j,2,k) = g5_deriv_z_yy_total( counter_deriv_zy + m, 2, k )
       g5_deriv_z_yy(i,j,3,k) = g5_deriv_z_yy_total( counter_deriv_zy + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_yz
   do j = 1, natom
     do i = 1, natom_z
       m = i + natom_z * ( j - 1 )
       g5_deriv_z_yz(i,j,1,k) = g5_deriv_z_yz_total( counter_deriv_zy + m, 1, k )
       g5_deriv_z_yz(i,j,2,k) = g5_deriv_z_yz_total( counter_deriv_zy + m, 2, k )
       g5_deriv_z_yz(i,j,3,k) = g5_deriv_z_yz_total( counter_deriv_zy + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_zz
   do j = 1, natom
     do i = 1, natom_z
       m = i + natom_z * ( j - 1 )
       g5_deriv_z_zz(i,j,1,k) = g5_deriv_z_zz_total( counter_deriv_z + m, 1, k )
       g5_deriv_z_zz(i,j,2,k) = g5_deriv_z_zz_total( counter_deriv_z + m, 2, k )
       g5_deriv_z_zz(i,j,3,k) = g5_deriv_z_zz_total( counter_deriv_z + m, 3, k )
     enddo
   enddo
 enddo


 counter_deriv_z   = counter_deriv_z   + natom_z * natom
 counter_deriv_zy  = counter_deriv_zy  + natom_z * natom

 endif ! io_z

!Z-XYZ
 if( io_x == 1 .and. io_y == 1 .and. io_z == 1 )then


 do k = 1, num_g2_z_x
   do j = 1, natom
     do i = 1, natom_z
       m = i + natom_z * ( j - 1 )
       g2_deriv_z_x(i,j,1,k) = g2_deriv_z_x_total( counter_deriv_zx + m, 1, k )
       g2_deriv_z_x(i,j,2,k) = g2_deriv_z_x_total( counter_deriv_zx + m, 2, k )
       g2_deriv_z_x(i,j,3,k) = g2_deriv_z_x_total( counter_deriv_zx + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g2_z_y
   do j = 1, natom
     do i = 1, natom_z
       m = i + natom_z * ( j - 1 )
       g2_deriv_z_y(i,j,1,k) = g2_deriv_z_y_total( counter_deriv_zy + m, 1, k )
       g2_deriv_z_y(i,j,2,k) = g2_deriv_z_y_total( counter_deriv_zy + m, 2, k )
       g2_deriv_z_y(i,j,3,k) = g2_deriv_z_y_total( counter_deriv_zy + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g2_z_z
   do j = 1, natom
     do i = 1, natom_z
       m = i + natom_z * ( j - 1 )
       g2_deriv_z_z(i,j,1,k) = g2_deriv_z_z_total( counter_deriv_z + m, 1, k )
       g2_deriv_z_z(i,j,2,k) = g2_deriv_z_z_total( counter_deriv_z + m, 2, k )
       g2_deriv_z_z(i,j,3,k) = g2_deriv_z_z_total( counter_deriv_z + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_xx
   do j = 1, natom
     do i = 1, natom_z
       m = i + natom_z * ( j - 1 )
       g5_deriv_z_xx(i,j,1,k) = g5_deriv_z_xx_total( counter_deriv_zx + m, 1, k )
       g5_deriv_z_xx(i,j,2,k) = g5_deriv_z_xx_total( counter_deriv_zx + m, 2, k )
       g5_deriv_z_xx(i,j,3,k) = g5_deriv_z_xx_total( counter_deriv_zx + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_xy
   do j = 1, natom
     do i = 1, natom_z
       m = i + natom_z * ( j - 1 )
       g5_deriv_z_xy(i,j,1,k) = g5_deriv_z_xy_total( counter_deriv_zxy + m, 1, k )
       g5_deriv_z_xy(i,j,2,k) = g5_deriv_z_xy_total( counter_deriv_zxy + m, 2, k )
       g5_deriv_z_xy(i,j,3,k) = g5_deriv_z_xy_total( counter_deriv_zxy + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_xz
   do j = 1, natom
     do i = 1, natom_z
       m = i + natom_z * ( j - 1 )
       g5_deriv_z_xz(i,j,1,k) = g5_deriv_z_xz_total( counter_deriv_zx + m, 1, k )
       g5_deriv_z_xz(i,j,2,k) = g5_deriv_z_xz_total( counter_deriv_zx + m, 2, k )
       g5_deriv_z_xz(i,j,3,k) = g5_deriv_z_xz_total( counter_deriv_zx + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_yy
   do j = 1, natom
     do i = 1, natom_z
       m = i + natom_z * ( j - 1 )
       g5_deriv_z_yy(i,j,1,k) = g5_deriv_z_yy_total( counter_deriv_zy + m, 1, k )
       g5_deriv_z_yy(i,j,2,k) = g5_deriv_z_yy_total( counter_deriv_zy + m, 2, k )
       g5_deriv_z_yy(i,j,3,k) = g5_deriv_z_yy_total( counter_deriv_zy + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_yz
   do j = 1, natom
     do i = 1, natom_z
       m = i + natom_z * ( j - 1 )
       g5_deriv_z_yz(i,j,1,k) = g5_deriv_z_yz_total( counter_deriv_zy + m, 1, k )
       g5_deriv_z_yz(i,j,2,k) = g5_deriv_z_yz_total( counter_deriv_zy + m, 2, k )
       g5_deriv_z_yz(i,j,3,k) = g5_deriv_z_yz_total( counter_deriv_zy + m, 3, k )
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_zz
   do j = 1, natom
     do i = 1, natom_z
       m = i + natom_z * ( j - 1 )
       g5_deriv_z_zz(i,j,1,k) = g5_deriv_z_zz_total( counter_deriv_z + m, 1, k )
       g5_deriv_z_zz(i,j,2,k) = g5_deriv_z_zz_total( counter_deriv_z + m, 2, k )
       g5_deriv_z_zz(i,j,3,k) = g5_deriv_z_zz_total( counter_deriv_z + m, 3, k )
     enddo
   enddo
 enddo


 counter_deriv_z   = counter_deriv_z   + natom_z * natom
 counter_deriv_zx  = counter_deriv_zx  + natom_z * natom
 counter_deriv_zy  = counter_deriv_zy  + natom_z * natom
 counter_deriv_zxy = counter_deriv_zxy + natom_z * natom

 endif ! io_z

!Z-X, Z-Y, Z-XY
 if( io_z == 0 )then


 do k = 1, num_g2_z_x
   do j = 1, natom
     do i = 1, natom_z
       g2_deriv_z_x(i,j,1,k) = 0.0d0
       g2_deriv_z_x(i,j,2,k) = 0.0d0
       g2_deriv_z_x(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g2_z_y
   do j = 1, natom
     do i = 1, natom_z
       g2_deriv_z_y(i,j,1,k) = 0.0d0
       g2_deriv_z_y(i,j,2,k) = 0.0d0
       g2_deriv_z_y(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g2_z_z
   do j = 1, natom
     do i = 1, natom_z
       g2_deriv_z_z(i,j,1,k) = 0.0d0
       g2_deriv_z_z(i,j,2,k) = 0.0d0
       g2_deriv_z_z(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_xx
   do j = 1, natom
     do i = 1, natom_z
       g5_deriv_z_xx(i,j,1,k) = 0.0d0
       g5_deriv_z_xx(i,j,2,k) = 0.0d0
       g5_deriv_z_xx(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_xy
   do j = 1, natom
     do i = 1, natom_z
       g5_deriv_z_xy(i,j,1,k) = 0.0d0
       g5_deriv_z_xy(i,j,2,k) = 0.0d0
       g5_deriv_z_xy(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_xz
   do j = 1, natom
     do i = 1, natom_z
       g5_deriv_z_xz(i,j,1,k) = 0.0d0
       g5_deriv_z_xz(i,j,2,k) = 0.0d0
       g5_deriv_z_xz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_yy
   do j = 1, natom
     do i = 1, natom_z
       g5_deriv_z_yy(i,j,1,k) = 0.0d0
       g5_deriv_z_yy(i,j,2,k) = 0.0d0
       g5_deriv_z_yy(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_yz
   do j = 1, natom
     do i = 1, natom_z
       g5_deriv_z_yz(i,j,1,k) = 0.0d0
       g5_deriv_z_yz(i,j,2,k) = 0.0d0
       g5_deriv_z_yz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 do k = 1, num_g5_z_zz
   do j = 1, natom
     do i = 1, natom_z
       g5_deriv_z_zz(i,j,1,k) = 0.0d0
       g5_deriv_z_zz(i,j,2,k) = 0.0d0
       g5_deriv_z_zz(i,j,3,k) = 0.0d0
     enddo
   enddo
 enddo


 endif ! io_z


 ! Energy & force calculation through HDNNP
 !------------------------------------------
 energy_x = 0.0d0; energy_y = 0.0d0; energy_z = 0.0d0
 force = 0.0d0
 !X
 if( io_x == 1 )then
 call energy_force_nnp3abc( &
      g2_x_x,  g2_deriv_x_x,  num_g2_x_x, &
      g2_x_y,  g2_deriv_x_y,  num_g2_x_y, &
      g2_x_z,  g2_deriv_x_z,  num_g2_x_z, &
      g5_x_xx, g5_deriv_x_xx, num_g5_x_xx, &
      g5_x_xy, g5_deriv_x_xy, num_g5_x_xy, &
      g5_x_xz, g5_deriv_x_xz, num_g5_x_xz, &
      g5_x_yy, g5_deriv_x_yy, num_g5_x_yy, &
      g5_x_yz, g5_deriv_x_yz, num_g5_x_yz, &
      g5_x_zz, g5_deriv_x_zz, num_g5_x_zz, &
      num_g_x, &
      network_x, nlayer_x, weight_x, nweight_x, hidden_x, nhidden_x, &
      natom, natom_x, energy_x, force )
 endif
 !Y
 if( io_y == 1 )then
 call energy_force_nnp3abc( &
      g2_y_x,  g2_deriv_y_x,  num_g2_y_x, &
      g2_y_y,  g2_deriv_y_y,  num_g2_y_y, &
      g2_y_z,  g2_deriv_y_z,  num_g2_y_z, &
      g5_y_xx, g5_deriv_y_xx, num_g5_y_xx, &
      g5_y_xy, g5_deriv_y_xy, num_g5_y_xy, &
      g5_y_xz, g5_deriv_y_xz, num_g5_y_xz, &
      g5_y_yy, g5_deriv_y_yy, num_g5_y_yy, &
      g5_y_yz, g5_deriv_y_yz, num_g5_y_yz, &
      g5_y_zz, g5_deriv_y_zz, num_g5_y_zz, &
      num_g_y, &
      network_y, nlayer_y, weight_y, nweight_y, hidden_y, nhidden_y, &
      natom, natom_y, energy_y, force )
 endif
 !Z
 if( io_z == 1 )then
 call energy_force_nnp3abc( &
      g2_z_x,  g2_deriv_z_x,  num_g2_z_x, &
      g2_z_y,  g2_deriv_z_y,  num_g2_z_y, &
      g2_z_z,  g2_deriv_z_z,  num_g2_z_z, &
      g5_z_xx, g5_deriv_z_xx, num_g5_z_xx, &
      g5_z_xy, g5_deriv_z_xy, num_g5_z_xy, &
      g5_z_xz, g5_deriv_z_xz, num_g5_z_xz, &
      g5_z_yy, g5_deriv_z_yy, num_g5_z_yy, &
      g5_z_yz, g5_deriv_z_yz, num_g5_z_yz, &
      g5_z_zz, g5_deriv_z_zz, num_g5_z_zz, &
      num_g_z, &
      network_z, nlayer_z, weight_z, nweight_z, hidden_z, nhidden_z, &
      natom, natom_z, energy_z, force )
 endif

 energy = energy_x + energy_y + energy_z 


end subroutine valid3abc
