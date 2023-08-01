module allocarray
contains

!************************************************************
! ALLOCATE_3
!************************************************************
subroutine allocate_3( &
           tag2_x_x,  eta2_x_x,  rs2_x_x, &
           tag2_x_y,  eta2_x_y,  rs2_x_y, &
           tag2_x_z,  eta2_x_z,  rs2_x_z, &
           tag5_x_xx, eta5_x_xx, theta5_x_xx, zeta5_x_xx, R_s_x_xx, &
           tag5_x_xy, eta5_x_xy, theta5_x_xy, zeta5_x_xy, R_s_x_xy, &
           tag5_x_xz, eta5_x_xz, theta5_x_xz, zeta5_x_xz, R_s_x_xz, &
           tag5_x_yy, eta5_x_yy, theta5_x_yy, zeta5_x_yy, R_s_x_yy, &
           tag5_x_yz, eta5_x_yz, theta5_x_yz, zeta5_x_yz, R_s_x_yz, &
           tag5_x_zz, eta5_x_zz, theta5_x_zz, zeta5_x_zz, R_s_x_zz, &
           num_g2_x_x,  num_g2_x_y,  num_g2_x_z, &
           num_g5_x_xx, num_g5_x_xy, num_g5_x_xz, &
           num_g5_x_yy, num_g5_x_yz, &
           num_g5_x_zz, &
           g2_x_x_maxmin,  g2_x_y_maxmin,  g2_x_z_maxmin, &
           g5_x_xx_maxmin, g5_x_xy_maxmin, g5_x_xz_maxmin, &
           g5_x_yy_maxmin, g5_x_yz_maxmin, &
           g5_x_zz_maxmin )

implicit none
integer natom, natom_x
character(len=20) sf_name
character(len=6) tag2_x_x,  tag2_x_y,  tag2_x_z
character(len=7) tag5_x_xx, tag5_x_xy, tag5_x_xz, &
                 tag5_x_yy, tag5_x_yz, &
                 tag5_x_zz
integer num_g2_x_x,  num_g2_x_y,  num_g2_x_z
integer num_g5_x_xx, num_g5_x_xy, num_g5_x_xz, &
        num_g5_x_yy, num_g5_x_yz, &
        num_g5_x_zz
double precision,allocatable,dimension(:,:) :: &
  g2_x_x_maxmin,  g2_x_y_maxmin,  g2_x_z_maxmin, &
  g5_x_xx_maxmin, g5_x_xy_maxmin, g5_x_xz_maxmin, &
  g5_x_yy_maxmin, g5_x_yz_maxmin, &
  g5_x_zz_maxmin
double precision,allocatable,dimension(:) :: &
  eta2_x_x,  rs2_x_x, &
  eta2_x_y,  rs2_x_y, &
  eta2_x_z,  rs2_x_z, &
  eta5_x_xx, theta5_x_xx, zeta5_x_xx, R_s_x_xx, &
  eta5_x_xy, theta5_x_xy, zeta5_x_xy, R_s_x_xy, &
  eta5_x_xz, theta5_x_xz, zeta5_x_xz, R_s_x_xz, &
  eta5_x_yy, theta5_x_yy, zeta5_x_yy, R_s_x_yy, &
  eta5_x_yz, theta5_x_yz, zeta5_x_yz, R_s_x_yz, &
  eta5_x_zz, theta5_x_zz, zeta5_x_zz, R_s_x_zz


 allocate( g2_x_x_maxmin( num_g2_x_x, 2 ), &
           eta2_x_x( num_g2_x_x ), &
           rs2_x_x( num_g2_x_x ) )
 allocate( g2_x_y_maxmin( num_g2_x_y, 2 ), &
           eta2_x_y( num_g2_x_y ), &
           rs2_x_y( num_g2_x_y ) )
 allocate( g2_x_z_maxmin( num_g2_x_z, 2 ), &
           eta2_x_z( num_g2_x_z ), &
           rs2_x_z( num_g2_x_z ) )

 allocate( g5_x_xx_maxmin( num_g5_x_xx, 2 ), &
           eta5_x_xx( num_g5_x_xx ), &
           theta5_x_xx( num_g5_x_xx ), &
           zeta5_x_xx( num_g5_x_xx ), &
           R_s_x_xx( num_g5_x_xx ) )
 allocate( g5_x_xy_maxmin( num_g5_x_xy, 2 ), &
           eta5_x_xy( num_g5_x_xy ), &
           theta5_x_xy( num_g5_x_xy ), &
           zeta5_x_xy( num_g5_x_xy ), &
           R_s_x_xy( num_g5_x_xy ) )
 allocate( g5_x_xz_maxmin( num_g5_x_xz, 2 ), &
           eta5_x_xz( num_g5_x_xz ), &
           theta5_x_xz( num_g5_x_xz ), &
           zeta5_x_xz( num_g5_x_xz ), &
           R_s_x_xz( num_g5_x_xz ) )
 allocate( g5_x_yy_maxmin( num_g5_x_yy, 2 ), &
           eta5_x_yy( num_g5_x_yy ), &
           theta5_x_yy( num_g5_x_yy ), &
           zeta5_x_yy( num_g5_x_yy ), &
           R_s_x_yy( num_g5_x_yy ) )
 allocate( g5_x_yz_maxmin( num_g5_x_yz, 2 ), &
           eta5_x_yz( num_g5_x_yz ), &
           theta5_x_yz( num_g5_x_yz ), &
           zeta5_x_yz( num_g5_x_yz ), &
           R_s_x_yz( num_g5_x_yz ) )
 allocate( g5_x_zz_maxmin( num_g5_x_zz, 2 ), &
           eta5_x_zz( num_g5_x_zz ), &
           theta5_x_zz( num_g5_x_zz ), &
           zeta5_x_zz( num_g5_x_zz ), &
           R_s_x_zz( num_g5_x_zz ) )

 sf_name = trim(adjustl( tag2_x_x ))
 call param_read_g2( sf_name, g2_x_x_maxmin, eta2_x_x, rs2_x_x, num_g2_x_x )
 sf_name = trim(adjustl( tag2_x_y ))
 call param_read_g2( sf_name, g2_x_y_maxmin, eta2_x_y, rs2_x_y, num_g2_x_y )
 sf_name = trim(adjustl( tag2_x_z ))
 call param_read_g2( sf_name, g2_x_z_maxmin, eta2_x_z, rs2_x_z, num_g2_x_z )

 sf_name = trim(adjustl( tag5_x_xx ))
 call param_read_g5( sf_name, g5_x_xx_maxmin, eta5_x_xx, theta5_x_xx, zeta5_x_xx, R_s_x_xx, num_g5_x_xx )
 sf_name = trim(adjustl( tag5_x_xy ))
 call param_read_g5( sf_name, g5_x_xy_maxmin, eta5_x_xy, theta5_x_xy, zeta5_x_xy, R_s_x_xy, num_g5_x_xy )
 sf_name = trim(adjustl( tag5_x_xz ))
 call param_read_g5( sf_name, g5_x_xz_maxmin, eta5_x_xz, theta5_x_xz, zeta5_x_xz, R_s_x_xz, num_g5_x_xz )
 sf_name = trim(adjustl( tag5_x_yy ))
 call param_read_g5( sf_name, g5_x_yy_maxmin, eta5_x_yy, theta5_x_yy, zeta5_x_yy, R_s_x_yy, num_g5_x_yy )
 sf_name = trim(adjustl( tag5_x_yz ))
 call param_read_g5( sf_name, g5_x_yz_maxmin, eta5_x_yz, theta5_x_yz, zeta5_x_yz, R_s_x_yz, num_g5_x_yz )
 sf_name = trim(adjustl( tag5_x_zz ))
 call param_read_g5( sf_name, g5_x_zz_maxmin, eta5_x_zz, theta5_x_zz, zeta5_x_zz, R_s_x_zz, num_g5_x_zz )


end subroutine allocate_3


!************************************************************
! ALLOCATE_4
!************************************************************
subroutine allocate_4( &
           tag2_x_x,  eta2_x_x,  rs2_x_x, &
           tag2_x_y,  eta2_x_y,  rs2_x_y, &
           tag2_x_z,  eta2_x_z,  rs2_x_z, &
           tag2_x_v,  eta2_x_v,  rs2_x_v, &
           tag5_x_xx, eta5_x_xx, theta5_x_xx, zeta5_x_xx, R_s_x_xx, &
           tag5_x_xy, eta5_x_xy, theta5_x_xy, zeta5_x_xy, R_s_x_xy, &
           tag5_x_xz, eta5_x_xz, theta5_x_xz, zeta5_x_xz, R_s_x_xz, &
           tag5_x_xv, eta5_x_xv, theta5_x_xv, zeta5_x_xv, R_s_x_xv, &
           tag5_x_yy, eta5_x_yy, theta5_x_yy, zeta5_x_yy, R_s_x_yy, &
           tag5_x_yz, eta5_x_yz, theta5_x_yz, zeta5_x_yz, R_s_x_yz, &
           tag5_x_yv, eta5_x_yv, theta5_x_yv, zeta5_x_yv, R_s_x_yv, &
           tag5_x_zz, eta5_x_zz, theta5_x_zz, zeta5_x_zz, R_s_x_zz, &
           tag5_x_zv, eta5_x_zv, theta5_x_zv, zeta5_x_zv, R_s_x_zv, &
           tag5_x_vv, eta5_x_vv, theta5_x_vv, zeta5_x_vv, R_s_x_vv, &
           num_g2_x_x,  num_g2_x_y,  num_g2_x_z,  num_g2_x_v, &
           num_g5_x_xx, num_g5_x_xy, num_g5_x_xz, num_g5_x_xv, &
           num_g5_x_yy, num_g5_x_yz, num_g5_x_yv, &
           num_g5_x_zz, num_g5_x_zv, &
           num_g5_x_vv, &
           g2_x_x_maxmin,  g2_x_y_maxmin,  g2_x_z_maxmin,  g2_x_v_maxmin, &
           g5_x_xx_maxmin, g5_x_xy_maxmin, g5_x_xz_maxmin, g5_x_xv_maxmin, &
           g5_x_yy_maxmin, g5_x_yz_maxmin, g5_x_yv_maxmin, &
           g5_x_zz_maxmin, g5_x_zv_maxmin, &
           g5_x_vv_maxmin )

implicit none
integer natom, natom_x
character(len=20) sf_name
character(len=6) tag2_x_x,  tag2_x_y,  tag2_x_z,  tag2_x_v
character(len=7) tag5_x_xx, tag5_x_xy, tag5_x_xz, tag5_x_xv, &
                 tag5_x_yy, tag5_x_yz, tag5_x_yv, &
                 tag5_x_zz, tag5_x_zv, &
                 tag5_x_vv
integer num_g2_x_x,  num_g2_x_y,  num_g2_x_z,  num_g2_x_v
integer num_g5_x_xx, num_g5_x_xy, num_g5_x_xz, num_g5_x_xv, &
        num_g5_x_yy, num_g5_x_yz, num_g5_x_yv, &
        num_g5_x_zz, num_g5_x_zv, &
        num_g5_x_vv
double precision,allocatable,dimension(:,:) :: &
  g2_x_x_maxmin,  g2_x_y_maxmin,  g2_x_z_maxmin,  g2_x_v_maxmin, &
  g5_x_xx_maxmin, g5_x_xy_maxmin, g5_x_xz_maxmin, g5_x_xv_maxmin, &
  g5_x_yy_maxmin, g5_x_yz_maxmin, g5_x_yv_maxmin, &
  g5_x_zz_maxmin, g5_x_zv_maxmin, &
  g5_x_vv_maxmin
double precision,allocatable,dimension(:) :: &
  eta2_x_x,  rs2_x_x, &
  eta2_x_y,  rs2_x_y, &
  eta2_x_z,  rs2_x_z, &
  eta2_x_v,  rs2_x_v, &
  eta5_x_xx, theta5_x_xx, zeta5_x_xx, R_s_x_xx, &
  eta5_x_xy, theta5_x_xy, zeta5_x_xy, R_s_x_xy, &
  eta5_x_xz, theta5_x_xz, zeta5_x_xz, R_s_x_xz, &
  eta5_x_xv, theta5_x_xv, zeta5_x_xv, R_s_x_xv, &
  eta5_x_yy, theta5_x_yy, zeta5_x_yy, R_s_x_yy, &
  eta5_x_yz, theta5_x_yz, zeta5_x_yz, R_s_x_yz, &
  eta5_x_yv, theta5_x_yv, zeta5_x_yv, R_s_x_yv, &
  eta5_x_zz, theta5_x_zz, zeta5_x_zz, R_s_x_zz, &
  eta5_x_zv, theta5_x_zv, zeta5_x_zv, R_s_x_zv, &
  eta5_x_vv, theta5_x_vv, zeta5_x_vv, R_s_x_vv


 allocate( g2_x_x_maxmin( num_g2_x_x, 2 ), &
           eta2_x_x( num_g2_x_x ), &
           rs2_x_x( num_g2_x_x ) )
 allocate( g2_x_y_maxmin( num_g2_x_y, 2 ), &
           eta2_x_y( num_g2_x_y ), &
           rs2_x_y( num_g2_x_y ) )
 allocate( g2_x_z_maxmin( num_g2_x_z, 2 ), &
           eta2_x_z( num_g2_x_z ), &
           rs2_x_z( num_g2_x_z ) )
 allocate( g2_x_v_maxmin( num_g2_x_v, 2 ), &
           eta2_x_v( num_g2_x_v ), &
           rs2_x_v( num_g2_x_v ) )

 allocate( g5_x_xx_maxmin( num_g5_x_xx, 2 ), &
           eta5_x_xx( num_g5_x_xx ), &
           theta5_x_xx( num_g5_x_xx ), &
           zeta5_x_xx( num_g5_x_xx ), &
           R_s_x_xx( num_g5_x_xx ) )
 allocate( g5_x_xy_maxmin( num_g5_x_xy, 2 ), &
           eta5_x_xy( num_g5_x_xy ), &
           theta5_x_xy( num_g5_x_xy ), &
           zeta5_x_xy( num_g5_x_xy ), &
           R_s_x_xy( num_g5_x_xy ) )
 allocate( g5_x_xz_maxmin( num_g5_x_xz, 2 ), &
           eta5_x_xz( num_g5_x_xz ), &
           theta5_x_xz( num_g5_x_xz ), &
           zeta5_x_xz( num_g5_x_xz ), &
           R_s_x_xz( num_g5_x_xz ) )
 allocate( g5_x_xv_maxmin( num_g5_x_xv, 2 ), &
           eta5_x_xv( num_g5_x_xv ), &
           theta5_x_xv( num_g5_x_xv ), &
           zeta5_x_xv( num_g5_x_xv ), &
           R_s_x_xv( num_g5_x_xv ) )
 allocate( g5_x_yy_maxmin( num_g5_x_yy, 2 ), &
           eta5_x_yy( num_g5_x_yy ), &
           theta5_x_yy( num_g5_x_yy ), &
           zeta5_x_yy( num_g5_x_yy ), &
           R_s_x_yy( num_g5_x_yy ) )
 allocate( g5_x_yz_maxmin( num_g5_x_yz, 2 ), &
           eta5_x_yz( num_g5_x_yz ), &
           theta5_x_yz( num_g5_x_yz ), &
           zeta5_x_yz( num_g5_x_yz ), &
           R_s_x_yz( num_g5_x_yz ) )
 allocate( g5_x_yv_maxmin( num_g5_x_yv, 2 ), &
           eta5_x_yv( num_g5_x_yv ), &
           theta5_x_yv( num_g5_x_yv ), &
           zeta5_x_yv( num_g5_x_yv ), &
           R_s_x_yv( num_g5_x_yv ) )
 allocate( g5_x_zz_maxmin( num_g5_x_zz, 2 ), &
           eta5_x_zz( num_g5_x_zz ), &
           theta5_x_zz( num_g5_x_zz ), &
           zeta5_x_zz( num_g5_x_zz ), &
           R_s_x_zz( num_g5_x_zz ) )
 allocate( g5_x_zv_maxmin( num_g5_x_zv, 2 ), &
           eta5_x_zv( num_g5_x_zv ), &
           theta5_x_zv( num_g5_x_zv ), &
           zeta5_x_zv( num_g5_x_zv ), &
           R_s_x_zv( num_g5_x_zv ) )
 allocate( g5_x_vv_maxmin( num_g5_x_vv, 2 ), &
           eta5_x_vv( num_g5_x_vv ), &
           theta5_x_vv( num_g5_x_vv ), &
           zeta5_x_vv( num_g5_x_vv ), &
           R_s_x_vv( num_g5_x_vv ) )

 sf_name = trim(adjustl( tag2_x_x ))
 call param_read_g2( sf_name, g2_x_x_maxmin, eta2_x_x, rs2_x_x, num_g2_x_x )
 sf_name = trim(adjustl( tag2_x_y ))
 call param_read_g2( sf_name, g2_x_y_maxmin, eta2_x_y, rs2_x_y, num_g2_x_y )
 sf_name = trim(adjustl( tag2_x_z ))
 call param_read_g2( sf_name, g2_x_z_maxmin, eta2_x_z, rs2_x_z, num_g2_x_z )
 sf_name = trim(adjustl( tag2_x_v ))
 call param_read_g2( sf_name, g2_x_v_maxmin, eta2_x_v, rs2_x_v, num_g2_x_v )

 sf_name = trim(adjustl( tag5_x_xx ))
 call param_read_g5( sf_name, g5_x_xx_maxmin, eta5_x_xx, theta5_x_xx, zeta5_x_xx, R_s_x_xx, num_g5_x_xx )
 sf_name = trim(adjustl( tag5_x_xy ))
 call param_read_g5( sf_name, g5_x_xy_maxmin, eta5_x_xy, theta5_x_xy, zeta5_x_xy, R_s_x_xy, num_g5_x_xy )
 sf_name = trim(adjustl( tag5_x_xz ))
 call param_read_g5( sf_name, g5_x_xz_maxmin, eta5_x_xz, theta5_x_xz, zeta5_x_xz, R_s_x_xz, num_g5_x_xz )
 sf_name = trim(adjustl( tag5_x_xv ))
 call param_read_g5( sf_name, g5_x_xv_maxmin, eta5_x_xv, theta5_x_xv, zeta5_x_xv, R_s_x_xv, num_g5_x_xv )
 sf_name = trim(adjustl( tag5_x_yy ))
 call param_read_g5( sf_name, g5_x_yy_maxmin, eta5_x_yy, theta5_x_yy, zeta5_x_yy, R_s_x_yy, num_g5_x_yy )
 sf_name = trim(adjustl( tag5_x_yz ))
 call param_read_g5( sf_name, g5_x_yz_maxmin, eta5_x_yz, theta5_x_yz, zeta5_x_yz, R_s_x_yz, num_g5_x_yz )
 sf_name = trim(adjustl( tag5_x_yv ))
 call param_read_g5( sf_name, g5_x_yv_maxmin, eta5_x_yv, theta5_x_yv, zeta5_x_yv, R_s_x_yv, num_g5_x_yv )
 sf_name = trim(adjustl( tag5_x_zz ))
 call param_read_g5( sf_name, g5_x_zz_maxmin, eta5_x_zz, theta5_x_zz, zeta5_x_zz, R_s_x_zz, num_g5_x_zz )
 sf_name = trim(adjustl( tag5_x_zv ))
 call param_read_g5( sf_name, g5_x_zv_maxmin, eta5_x_zv, theta5_x_zv, zeta5_x_zv, R_s_x_zv, num_g5_x_zv )
 sf_name = trim(adjustl( tag5_x_vv ))
 call param_read_g5( sf_name, g5_x_vv_maxmin, eta5_x_vv, theta5_x_vv, zeta5_x_vv, R_s_x_vv, num_g5_x_vv )


end subroutine allocate_4


end module
