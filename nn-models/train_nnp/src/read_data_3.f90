subroutine read_data_3( &
  num_atom_x, num_atom_y, num_atom_z, &
  num_data_x, num_data_y, num_data_z, &
  natom_x_total, natom_xy_total, natom_xz_total, natom_xyz_total, &
  natom_y_total, natom_yx_total, natom_yz_total, natom_yxz_total, &
  natom_z_total, natom_zx_total, natom_zy_total, natom_zxy_total, &
  natom_x_sum, natom_xy_sum, natom_xz_sum, natom_xyz_sum, &
  natom_y_sum, natom_yx_sum, natom_yz_sum, natom_yxz_sum, &
  natom_z_sum, natom_zx_sum, natom_zy_sum, natom_zxy_sum, & 
  counter_x, counter_y, counter_z, &
  counter_natom_x, counter_natom_xy, counter_natom_xz, counter_natom_xyz, &
  counter_natom_y, counter_natom_yx, counter_natom_yz, counter_natom_yxz, &
  counter_natom_z, counter_natom_zx, counter_natom_zy, counter_natom_zxy, &
  counter_deriv_x, counter_deriv_xy, counter_deriv_xz, counter_deriv_xyz, &
  counter_deriv_y, counter_deriv_yx, counter_deriv_yz, counter_deriv_yxz, &
  counter_deriv_z, counter_deriv_zx, counter_deriv_zy, counter_deriv_zxy, &
  imd, natom_x, natom_y, natom_z, natom, &
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
integer i, j, k, m, imd
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
integer num_data_x, num_data_y, num_data_z 
integer natom, natom_x, natom_y, natom_z
integer natom_x_total, natom_xy_total, natom_xz_total, natom_xyz_total, &
        natom_y_total, natom_yx_total, natom_yz_total, natom_yxz_total, &
        natom_z_total, natom_zx_total, natom_zy_total, natom_zxy_total
integer natom_x_sum, natom_xy_sum, natom_xz_sum, natom_xyz_sum, &
        natom_y_sum, natom_yx_sum, natom_yz_sum, natom_yxz_sum, &
        natom_z_sum, natom_zx_sum, natom_zy_sum, natom_zxy_sum
integer num_atom_x( num_data_x ), &
        num_atom_y( num_data_y ), &
        num_atom_z( num_data_z )
integer counter_x, counter_y, counter_z
integer counter_natom_x, counter_natom_xy, counter_natom_xz, counter_natom_xyz, &
        counter_natom_y, counter_natom_yx, counter_natom_yz, counter_natom_yxz, &
        counter_natom_z, counter_natom_zx, counter_natom_zy, counter_natom_zxy
integer counter_deriv_x, counter_deriv_xy, counter_deriv_xz, counter_deriv_xyz, &
        counter_deriv_y, counter_deriv_yx, counter_deriv_yz, counter_deriv_yxz, &
        counter_deriv_z, counter_deriv_zx, counter_deriv_zy, counter_deriv_zxy
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
double precision,allocatable,dimension(:,:) :: g2, g5
double precision :: &
  g2_x_x(  natom_x_total,   num_g2_x_x ), &
  g2_x_y(  natom_xy_total,  num_g2_x_y ), &
  g2_x_z(  natom_xz_total,  num_g2_x_z ), &
  g5_x_xx( natom_x_total,   num_g5_x_xx ), &
  g5_x_xy( natom_xy_total,  num_g5_x_xy ), &
  g5_x_xz( natom_xz_total,  num_g5_x_xz ), &
  g5_x_yy( natom_xy_total,  num_g5_x_yy ), &
  g5_x_yz( natom_xyz_total, num_g5_x_yz ), &
  g5_x_zz( natom_xz_total,  num_g5_x_zz )
double precision :: &
  g2_y_x(  natom_yx_total,  num_g2_y_x ), &
  g2_y_y(  natom_y_total,   num_g2_y_y ), &
  g2_y_z(  natom_yz_total,  num_g2_y_z ), &
  g5_y_xx( natom_yx_total,  num_g5_y_xx ), &
  g5_y_xy( natom_yx_total,  num_g5_y_xy ), &
  g5_y_xz( natom_yxz_total, num_g5_y_xz ), &
  g5_y_yy( natom_y_total,   num_g5_y_yy ), &
  g5_y_yz( natom_yz_total,  num_g5_y_yz ), &
  g5_y_zz( natom_yz_total,  num_g5_y_zz )
double precision :: &
  g2_z_x(  natom_zx_total,  num_g2_z_x ), &
  g2_z_y(  natom_zy_total,  num_g2_z_y ), &
  g2_z_z(  natom_z_total,   num_g2_z_z ), &
  g5_z_xx( natom_zx_total,  num_g5_z_xx ), &
  g5_z_xy( natom_zxy_total, num_g5_z_xy ), &
  g5_z_xz( natom_zx_total,  num_g5_z_xz ), &
  g5_z_yy( natom_zy_total,  num_g5_z_yy ), &
  g5_z_yz( natom_zy_total,  num_g5_z_yz ), &
  g5_z_zz( natom_z_total,   num_g5_z_zz )
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
double precision,allocatable,dimension(:,:,:,:) :: g2_deriv, g5_deriv
double precision :: &
  g2_deriv_x_x(  natom_x_sum,   3, num_g2_x_x ), &
  g2_deriv_x_y(  natom_xy_sum,  3, num_g2_x_y ), &
  g2_deriv_x_z(  natom_xz_sum,  3, num_g2_x_z ), &
  g5_deriv_x_xx( natom_x_sum,   3, num_g5_x_xx ),&
  g5_deriv_x_xy( natom_xy_sum,  3, num_g5_x_xy ), &
  g5_deriv_x_xz( natom_xz_sum,  3, num_g5_x_xz ), &
  g5_deriv_x_yy( natom_xy_sum,  3, num_g5_x_yy ), &
  g5_deriv_x_yz( natom_xyz_sum, 3, num_g5_x_yz ), &
  g5_deriv_x_zz( natom_xz_sum,  3, num_g5_x_zz )
double precision :: &
  g2_deriv_y_x(  natom_yx_sum,  3, num_g2_y_x ), &
  g2_deriv_y_y(  natom_y_sum,   3, num_g2_y_y ), &
  g2_deriv_y_z(  natom_yz_sum,  3, num_g2_y_z ), &
  g5_deriv_y_xx( natom_yx_sum,  3, num_g5_y_xx ), &
  g5_deriv_y_xy( natom_yx_sum,  3, num_g5_y_xy ), &
  g5_deriv_y_xz( natom_yxz_sum, 3, num_g5_y_xz ), &
  g5_deriv_y_yy( natom_y_sum,   3, num_g5_y_yy ), &
  g5_deriv_y_yz( natom_yz_sum,  3, num_g5_y_yz ), &
  g5_deriv_y_zz( natom_yz_sum,  3, num_g5_y_zz )
double precision :: &
  g2_deriv_z_x(  natom_zx_sum,  3, num_g2_z_x ), &
  g2_deriv_z_y(  natom_zy_sum,  3, num_g2_z_y ), &
  g2_deriv_z_z(  natom_z_sum,   3, num_g2_z_z ), &
  g5_deriv_z_xx( natom_zx_sum,  3, num_g5_z_xx ), &
  g5_deriv_z_xy( natom_zxy_sum, 3, num_g5_z_xy ), &
  g5_deriv_z_xz( natom_zx_sum,  3, num_g5_z_xz ), &
  g5_deriv_z_yy( natom_zy_sum,  3, num_g5_z_yy ), &
  g5_deriv_z_yz( natom_zy_sum,  3, num_g5_z_yz ), &
  g5_deriv_z_zz( natom_z_sum,   3, num_g5_z_zz )


!X
  num_atom_x( imd + counter_x ) = natom_x

  allocate( g2( natom_x, num_g2_x_x ) )
  read(11) g2
  do j = 1, num_g2_x_x
    do i = 1, natom_x
      g2_x_x( i + counter_natom_x, j ) = ( 2.0d0*( g2(i,j) - g2_x_x_maxmin(j,2) )/ &
                                         ( g2_x_x_maxmin(j,1) - g2_x_x_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g2 )

  allocate( g2( natom_x, num_g2_x_y ) )
  read(11) g2
  do j = 1, num_g2_x_y
    do i = 1, natom_x
      g2_x_y( i + counter_natom_xy, j ) = ( 2.0d0*( g2(i,j) - g2_x_y_maxmin(j,2) )/ &
                                          ( g2_x_y_maxmin(j,1) - g2_x_y_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g2 )

  allocate( g2( natom_x, num_g2_x_z ) )
  read(11) g2
  do j = 1, num_g2_x_z
    do i = 1, natom_x
      g2_x_z( i + counter_natom_xz, j ) = ( 2.0d0*( g2(i,j) - g2_x_z_maxmin(j,2) )/ &
                                          ( g2_x_z_maxmin(j,1) - g2_x_z_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g2 )

  allocate( g5( natom_x, num_g5_x_xx ) )
  read(11) g5
  do j = 1, num_g5_x_xx
    do i = 1, natom_x
      g5_x_xx( i + counter_natom_x, j ) = ( 2.0d0*( g5(i,j) - g5_x_xx_maxmin(j,2) )/ &
                                          ( g5_x_xx_maxmin(j,1) - g5_x_xx_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  allocate( g5( natom_x, num_g5_x_xy ) )
  read(11) g5
  do j = 1, num_g5_x_xy
    do i = 1, natom_x
      g5_x_xy( i + counter_natom_xy, j ) = ( 2.0d0*( g5(i,j) - g5_x_xy_maxmin(j,2) )/ &
                                           ( g5_x_xy_maxmin(j,1) - g5_x_xy_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  allocate( g5( natom_x, num_g5_x_xz ) )
  read(11) g5
  do j = 1, num_g5_x_xz
    do i = 1, natom_x
      g5_x_xz( i + counter_natom_xz, j ) = ( 2.0d0*( g5(i,j) - g5_x_xz_maxmin(j,2) )/&
                                           ( g5_x_xz_maxmin(j,1) - g5_x_xz_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  allocate( g5( natom_x, num_g5_x_yy ) )
  read(11) g5
  do j = 1, num_g5_x_yy
    do i = 1, natom_x
      g5_x_yy( i + counter_natom_xy, j ) = ( 2.0d0*( g5(i,j) - g5_x_yy_maxmin(j,2) )/ &
                                           ( g5_x_yy_maxmin(j,1) - g5_x_yy_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  allocate( g5( natom_x, num_g5_x_yz ) )
  read(11) g5
  do j = 1, num_g5_x_yz
    do i = 1, natom_x
      g5_x_yz( i + counter_natom_xyz, j ) = ( 2.0d0*( g5(i,j) - g5_x_yz_maxmin(j,2) )/ &
                                            ( g5_x_yz_maxmin(j,1) - g5_x_yz_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  allocate( g5( natom_x, num_g5_x_zz ) )
  read(11) g5
  do j = 1, num_g5_x_zz
    do i = 1, natom_x
      g5_x_zz( i + counter_natom_xz, j ) = ( 2.0d0*( g5(i,j) - g5_x_zz_maxmin(j,2) )/ &
                                           ( g5_x_zz_maxmin(j,1) - g5_x_zz_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  counter_natom_x   = counter_natom_x   + natom_x
  counter_natom_xy  = counter_natom_xy  + natom_x
  counter_natom_xz  = counter_natom_xz  + natom_x
  counter_natom_xyz = counter_natom_xyz + natom_x


  ! SF_deriv
  allocate( g2_deriv( natom_x, natom, 3, num_g2_x_x ) )
  read(12) g2_deriv
  do k = 1, num_g2_x_x
    do j = 1, natom
      do i = 1, natom_x
        m = i + natom_x * ( j - 1 )
        g2_deriv_x_x( counter_deriv_x + m,1,k) = 2.0d0*g2_deriv(i,j,1,k)/ &
                                                 ( g2_x_x_maxmin(k,1) - g2_x_x_maxmin(k,2) )
        g2_deriv_x_x( counter_deriv_x + m,2,k) = 2.0d0*g2_deriv(i,j,2,k)/ &
                                                 ( g2_x_x_maxmin(k,1) - g2_x_x_maxmin(k,2) )
        g2_deriv_x_x( counter_deriv_x + m,3,k) = 2.0d0*g2_deriv(i,j,3,k)/ &
                                                 ( g2_x_x_maxmin(k,1) - g2_x_x_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g2_deriv )

  allocate( g2_deriv( natom_x, natom, 3, num_g2_x_y ) )
  read(12) g2_deriv
  do k = 1, num_g2_x_y
    do j = 1, natom
      do i = 1, natom_x
        m = i + natom_x * ( j - 1 )
        g2_deriv_x_y( counter_deriv_xy + m,1,k) = 2.0d0*g2_deriv(i,j,1,k)/ &
                                                 ( g2_x_y_maxmin(k,1) - g2_x_y_maxmin(k,2) )
        g2_deriv_x_y( counter_deriv_xy + m,2,k) = 2.0d0*g2_deriv(i,j,2,k)/ &
                                                 ( g2_x_y_maxmin(k,1) - g2_x_y_maxmin(k,2) )
        g2_deriv_x_y( counter_deriv_xy + m,3,k) = 2.0d0*g2_deriv(i,j,3,k)/ &
                                                 ( g2_x_y_maxmin(k,1) - g2_x_y_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g2_deriv )

  allocate( g2_deriv( natom_x, natom, 3, num_g2_x_z ) )
  read(12) g2_deriv
  do k = 1, num_g2_x_z
    do j = 1, natom
      do i = 1, natom_x
        m = i + natom_x * ( j - 1 )
        g2_deriv_x_z( counter_deriv_xz + m,1,k) = 2.0d0*g2_deriv(i,j,1,k)/ &
                                                 ( g2_x_z_maxmin(k,1) - g2_x_z_maxmin(k,2) )
        g2_deriv_x_z( counter_deriv_xz + m,2,k) = 2.0d0*g2_deriv(i,j,2,k)/ &
                                                 ( g2_x_z_maxmin(k,1) - g2_x_z_maxmin(k,2) )
        g2_deriv_x_z( counter_deriv_xz + m,3,k) = 2.0d0*g2_deriv(i,j,3,k)/ &
                                                 ( g2_x_z_maxmin(k,1) - g2_x_z_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g2_deriv )

  allocate( g5_deriv( natom_x, natom, 3, num_g5_x_xx ) )
  read(12) g5_deriv
  do k = 1, num_g5_x_xx
    do j = 1, natom
      do i = 1, natom_x
        m = i + natom_x * ( j - 1 )
        g5_deriv_x_xx( counter_deriv_x + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                 ( g5_x_xx_maxmin(k,1) - g5_x_xx_maxmin(k,2) )
        g5_deriv_x_xx( counter_deriv_x + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                 ( g5_x_xx_maxmin(k,1) - g5_x_xx_maxmin(k,2) )
        g5_deriv_x_xx( counter_deriv_x + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                 ( g5_x_xx_maxmin(k,1) - g5_x_xx_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  allocate( g5_deriv( natom_x, natom, 3, num_g5_x_xy ) )
  read(12) g5_deriv
  do k = 1, num_g5_x_xy
    do j = 1, natom
      do i = 1, natom_x
        m = i + natom_x * ( j - 1 )
        g5_deriv_x_xy( counter_deriv_xy + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                 ( g5_x_xy_maxmin(k,1) - g5_x_xy_maxmin(k,2) )
        g5_deriv_x_xy( counter_deriv_xy + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                 ( g5_x_xy_maxmin(k,1) - g5_x_xy_maxmin(k,2) )
        g5_deriv_x_xy( counter_deriv_xy + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                 ( g5_x_xy_maxmin(k,1) - g5_x_xy_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  allocate( g5_deriv( natom_x, natom, 3, num_g5_x_xz ) )
  read(12) g5_deriv
  do k = 1, num_g5_x_xz
    do j = 1, natom
      do i = 1, natom_x
        m = i + natom_x * ( j - 1 )
        g5_deriv_x_xz( counter_deriv_xz + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                 ( g5_x_xz_maxmin(k,1) - g5_x_xz_maxmin(k,2) )
        g5_deriv_x_xz( counter_deriv_xz + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                 ( g5_x_xz_maxmin(k,1) - g5_x_xz_maxmin(k,2) )
        g5_deriv_x_xz( counter_deriv_xz + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                 ( g5_x_xz_maxmin(k,1) - g5_x_xz_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  allocate( g5_deriv( natom_x, natom, 3, num_g5_x_yy ) )
  read(12) g5_deriv
  do k = 1, num_g5_x_yy
    do j = 1, natom
      do i = 1, natom_x
        m = i + natom_x * ( j - 1 )
        g5_deriv_x_yy( counter_deriv_xy + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                 ( g5_x_yy_maxmin(k,1) - g5_x_yy_maxmin(k,2) )
        g5_deriv_x_yy( counter_deriv_xy + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                 ( g5_x_yy_maxmin(k,1) - g5_x_yy_maxmin(k,2) )
        g5_deriv_x_yy( counter_deriv_xy + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                 ( g5_x_yy_maxmin(k,1) - g5_x_yy_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  allocate( g5_deriv( natom_x, natom, 3, num_g5_x_yz ) )
  read(12) g5_deriv
  do k = 1, num_g5_x_yz
    do j = 1, natom
      do i = 1, natom_x
        m = i + natom_x * ( j - 1 )
        g5_deriv_x_yz( counter_deriv_xyz + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                 ( g5_x_yz_maxmin(k,1) - g5_x_yz_maxmin(k,2) )
        g5_deriv_x_yz( counter_deriv_xyz + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                 ( g5_x_yz_maxmin(k,1) - g5_x_yz_maxmin(k,2) )
        g5_deriv_x_yz( counter_deriv_xyz + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                 ( g5_x_yz_maxmin(k,1) - g5_x_yz_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  allocate( g5_deriv( natom_x, natom, 3, num_g5_x_zz ) )
  read(12) g5_deriv
  do k = 1, num_g5_x_zz
    do j = 1, natom
      do i = 1, natom_x
        m = i + natom_x * ( j - 1 )
        g5_deriv_x_zz( counter_deriv_xz + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                 ( g5_x_zz_maxmin(k,1) - g5_x_zz_maxmin(k,2) )
        g5_deriv_x_zz( counter_deriv_xz + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                 ( g5_x_zz_maxmin(k,1) - g5_x_zz_maxmin(k,2) )
        g5_deriv_x_zz( counter_deriv_xz + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                 ( g5_x_zz_maxmin(k,1) - g5_x_zz_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  counter_deriv_x   = counter_deriv_x   + natom_x * natom
  counter_deriv_xy  = counter_deriv_xy  + natom_x * natom
  counter_deriv_xz  = counter_deriv_xz  + natom_x * natom
  counter_deriv_xyz = counter_deriv_xyz + natom_x * natom


!Y
  num_atom_y( imd + counter_y ) = natom_y

  allocate( g2( natom_y, num_g2_y_x ) )
  read(11) g2
  do j = 1, num_g2_y_x
    do i = 1, natom_y
      g2_y_x( i + counter_natom_yx, j ) = ( 2.0d0*( g2(i,j) - g2_y_x_maxmin(j,2) )/ &
                                          ( g2_y_x_maxmin(j,1) - g2_y_x_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g2 )

  allocate( g2( natom_y, num_g2_y_y ) )
  read(11) g2
  do j = 1, num_g2_y_y
    do i = 1, natom_y
      g2_y_y( i + counter_natom_y, j ) = ( 2.0d0*( g2(i,j) - g2_y_y_maxmin(j,2) )/ &
                                         ( g2_y_y_maxmin(j,1) - g2_y_y_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g2 )

  allocate( g2( natom_y, num_g2_y_z ) )
  read(11) g2
  do j = 1, num_g2_y_z
    do i = 1, natom_y
      g2_y_z( i + counter_natom_yz, j ) = ( 2.0d0*( g2(i,j) - g2_y_z_maxmin(j,2) )/ &
                                          ( g2_y_z_maxmin(j,1) - g2_y_z_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g2 )

  allocate( g5( natom_y, num_g5_y_xx ) )
  read(11) g5
  do j = 1, num_g5_y_xx
    do i = 1, natom_y
      g5_y_xx( i + counter_natom_yx, j ) = ( 2.0d0*( g5(i,j) - g5_y_xx_maxmin(j,2) )/ &
                                           ( g5_y_xx_maxmin(j,1) - g5_y_xx_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  allocate( g5( natom_y, num_g5_y_xy ) )
  read(11) g5
  do j = 1, num_g5_y_xy
    do i = 1, natom_y
      g5_y_xy( i + counter_natom_yx, j ) = ( 2.0d0*( g5(i,j) - g5_y_xy_maxmin(j,2) )/ &
                                           ( g5_y_xy_maxmin(j,1) - g5_y_xy_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  allocate( g5( natom_y, num_g5_y_xz ) )
  read(11) g5
  do j = 1, num_g5_y_xz
    do i = 1, natom_y
      g5_y_xz( i + counter_natom_yxz, j ) = ( 2.0d0*( g5(i,j) - g5_y_xz_maxmin(j,2) )/ &
                                            ( g5_y_xz_maxmin(j,1) - g5_y_xz_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  allocate( g5( natom_y, num_g5_y_yy ) )
  read(11) g5
  do j = 1, num_g5_y_yy
    do i = 1, natom_y
      g5_y_yy( i + counter_natom_y, j ) = ( 2.0d0*( g5(i,j) - g5_y_yy_maxmin(j,2) )/ &
                                          ( g5_y_yy_maxmin(j,1) - g5_y_yy_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  allocate( g5( natom_y, num_g5_y_yz ) )
  read(11) g5
  do j = 1, num_g5_y_yz
    do i = 1, natom_y
      g5_y_yz( i + counter_natom_yz, j ) = ( 2.0d0*( g5(i,j) - g5_y_yz_maxmin(j,2) )/ &
                                           ( g5_y_yz_maxmin(j,1) - g5_y_yz_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  allocate( g5( natom_y, num_g5_y_zz ) )
  read(11) g5
  do j = 1, num_g5_y_zz
    do i = 1, natom_y
      g5_y_zz( i + counter_natom_yz, j ) = ( 2.0d0*( g5(i,j) - g5_y_zz_maxmin(j,2) )/ &
                                           ( g5_y_zz_maxmin(j,1) - g5_y_zz_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  counter_natom_y   = counter_natom_y   + natom_y
  counter_natom_yx  = counter_natom_yx  + natom_y
  counter_natom_yz  = counter_natom_yz  + natom_y
  counter_natom_yxz = counter_natom_yxz + natom_y


  ! SF_deriv
  allocate( g2_deriv( natom_y, natom, 3, num_g2_y_x ) )
  read(12) g2_deriv
  do k = 1, num_g2_y_x
    do j = 1, natom
      do i = 1, natom_y
        m = i + natom_y * ( j - 1 )
        g2_deriv_y_x( counter_deriv_yx + m,1,k) = 2.0d0*g2_deriv(i,j,1,k)/ &
                                                  ( g2_y_x_maxmin(k,1) - g2_y_x_maxmin(k,2) )
        g2_deriv_y_x( counter_deriv_yx + m,2,k) = 2.0d0*g2_deriv(i,j,2,k)/ &
                                                  ( g2_y_x_maxmin(k,1) - g2_y_x_maxmin(k,2) )
        g2_deriv_y_x( counter_deriv_yx + m,3,k) = 2.0d0*g2_deriv(i,j,3,k)/ &
                                                  ( g2_y_x_maxmin(k,1) - g2_y_x_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g2_deriv )

  allocate( g2_deriv( natom_y, natom, 3, num_g2_y_y ) )
  read(12) g2_deriv
  do k = 1, num_g2_y_y
    do j = 1, natom
      do i = 1, natom_y
        m = i + natom_y * ( j - 1 )
        g2_deriv_y_y( counter_deriv_y + m,1,k) = 2.0d0*g2_deriv(i,j,1,k)/ &
                                                  ( g2_y_y_maxmin(k,1) - g2_y_y_maxmin(k,2) )
        g2_deriv_y_y( counter_deriv_y + m,2,k) = 2.0d0*g2_deriv(i,j,2,k)/ &
                                                  ( g2_y_y_maxmin(k,1) - g2_y_y_maxmin(k,2) )
        g2_deriv_y_y( counter_deriv_y + m,3,k) = 2.0d0*g2_deriv(i,j,3,k)/ &
                                                  ( g2_y_y_maxmin(k,1) - g2_y_y_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g2_deriv )

  allocate( g2_deriv( natom_y, natom, 3, num_g2_y_z ) )
  read(12) g2_deriv
  do k = 1, num_g2_y_z
    do j = 1, natom
      do i = 1, natom_y
        m = i + natom_y * ( j - 1 )
        g2_deriv_y_z( counter_deriv_yz + m,1,k) = 2.0d0*g2_deriv(i,j,1,k)/ &
                                                  ( g2_y_z_maxmin(k,1) - g2_y_z_maxmin(k,2) )
        g2_deriv_y_z( counter_deriv_yz + m,2,k) = 2.0d0*g2_deriv(i,j,2,k)/ &
                                                  ( g2_y_z_maxmin(k,1) - g2_y_z_maxmin(k,2) )
        g2_deriv_y_z( counter_deriv_yz + m,3,k) = 2.0d0*g2_deriv(i,j,3,k)/ &
                                                  ( g2_y_z_maxmin(k,1) - g2_y_z_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g2_deriv )

  allocate( g5_deriv( natom_y, natom, 3, num_g5_y_xx ) )
  read(12) g5_deriv
  do k = 1, num_g5_y_xx
    do j = 1, natom
      do i = 1, natom_y
        m = i + natom_y * ( j - 1 )
        g5_deriv_y_xx( counter_deriv_yx + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                  ( g5_y_xx_maxmin(k,1) - g5_y_xx_maxmin(k,2) )
        g5_deriv_y_xx( counter_deriv_yx + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                  ( g5_y_xx_maxmin(k,1) - g5_y_xx_maxmin(k,2) )
        g5_deriv_y_xx( counter_deriv_yx + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                  ( g5_y_xx_maxmin(k,1) - g5_y_xx_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  allocate( g5_deriv( natom_y, natom, 3, num_g5_y_xy ) )
  read(12) g5_deriv
  do k = 1, num_g5_y_xy
    do j = 1, natom
      do i = 1, natom_y
        m = i + natom_y * ( j - 1 )
        g5_deriv_y_xy( counter_deriv_yx + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                  ( g5_y_xy_maxmin(k,1) - g5_y_xy_maxmin(k,2) )
        g5_deriv_y_xy( counter_deriv_yx + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                  ( g5_y_xy_maxmin(k,1) - g5_y_xy_maxmin(k,2) )
        g5_deriv_y_xy( counter_deriv_yx + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                  ( g5_y_xy_maxmin(k,1) - g5_y_xy_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  allocate( g5_deriv( natom_y, natom, 3, num_g5_y_xz ) )
  read(12) g5_deriv
  do k = 1, num_g5_y_xz
    do j = 1, natom
      do i = 1, natom_y
        m = i + natom_y * ( j - 1 )
        g5_deriv_y_xz( counter_deriv_yxz + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                  ( g5_y_xz_maxmin(k,1) - g5_y_xz_maxmin(k,2) )
        g5_deriv_y_xz( counter_deriv_yxz + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                  ( g5_y_xz_maxmin(k,1) - g5_y_xz_maxmin(k,2) )
        g5_deriv_y_xz( counter_deriv_yxz + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                  ( g5_y_xz_maxmin(k,1) - g5_y_xz_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  allocate( g5_deriv( natom_y, natom, 3, num_g5_y_yy ) )
  read(12) g5_deriv
  do k = 1, num_g5_y_yy
    do j = 1, natom
      do i = 1, natom_y
        m = i + natom_y * ( j - 1 )
        g5_deriv_y_yy( counter_deriv_y + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                  ( g5_y_yy_maxmin(k,1) - g5_y_yy_maxmin(k,2) )
        g5_deriv_y_yy( counter_deriv_y + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                  ( g5_y_yy_maxmin(k,1) - g5_y_yy_maxmin(k,2) )
        g5_deriv_y_yy( counter_deriv_y + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                  ( g5_y_yy_maxmin(k,1) - g5_y_yy_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  allocate( g5_deriv( natom_y, natom, 3, num_g5_y_yz ) )
  read(12) g5_deriv
  do k = 1, num_g5_y_yz
    do j = 1, natom
      do i = 1, natom_y
        m = i + natom_y * ( j - 1 )
        g5_deriv_y_yz( counter_deriv_yz + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                  ( g5_y_yz_maxmin(k,1) - g5_y_yz_maxmin(k,2) )
        g5_deriv_y_yz( counter_deriv_yz + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                  ( g5_y_yz_maxmin(k,1) - g5_y_yz_maxmin(k,2) )
        g5_deriv_y_yz( counter_deriv_yz + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                  ( g5_y_yz_maxmin(k,1) - g5_y_yz_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  allocate( g5_deriv( natom_y, natom, 3, num_g5_y_zz ) )
  read(12) g5_deriv
  do k = 1, num_g5_y_zz
    do j = 1, natom
      do i = 1, natom_y
        m = i + natom_y * ( j - 1 )
        g5_deriv_y_zz( counter_deriv_yz + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                  ( g5_y_zz_maxmin(k,1) - g5_y_zz_maxmin(k,2) )
        g5_deriv_y_zz( counter_deriv_yz + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                  ( g5_y_zz_maxmin(k,1) - g5_y_zz_maxmin(k,2) )
        g5_deriv_y_zz( counter_deriv_yz + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                  ( g5_y_zz_maxmin(k,1) - g5_y_zz_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  counter_deriv_y   = counter_deriv_y   + natom_y * natom
  counter_deriv_yx  = counter_deriv_yx  + natom_y * natom
  counter_deriv_yz  = counter_deriv_yz  + natom_y * natom
  counter_deriv_yxz = counter_deriv_yxz + natom_y * natom


!Z
  num_atom_z( imd + counter_z ) = natom_z

  allocate( g2( natom_z, num_g2_z_x ) )
  read(11) g2
  do j = 1, num_g2_z_x
    do i = 1, natom_z
      g2_z_x( i + counter_natom_zx, j ) = ( 2.0d0*( g2(i,j) - g2_z_x_maxmin(j,2) )/ &
                                          ( g2_z_x_maxmin(j,1) - g2_z_x_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g2 )

  allocate( g2( natom_z, num_g2_z_y ) )
  read(11) g2
  do j = 1, num_g2_z_y
    do i = 1, natom_z
      g2_z_y( i + counter_natom_zy, j ) = ( 2.0d0*( g2(i,j) - g2_z_y_maxmin(j,2) )/ &
                                          ( g2_z_y_maxmin(j,1) - g2_z_y_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g2 )

  allocate( g2( natom_z, num_g2_z_z ) )
  read(11) g2
  do j = 1, num_g2_z_z
    do i = 1, natom_z
      g2_z_z( i + counter_natom_z, j ) = ( 2.0d0*( g2(i,j) - g2_z_z_maxmin(j,2) )/ &
                                         ( g2_z_z_maxmin(j,1) - g2_z_z_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g2 )

  allocate( g5( natom_z, num_g5_z_xx ) )
  read(11) g5
  do j = 1, num_g5_z_xx
    do i = 1, natom_z
      g5_z_xx( i + counter_natom_zx, j ) = ( 2.0d0*( g5(i,j) - g5_z_xx_maxmin(j,2) )/ &
                                           ( g5_z_xx_maxmin(j,1) - g5_z_xx_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  allocate( g5( natom_z, num_g5_z_xy ) )
  read(11) g5
  do j = 1, num_g5_z_xy
    do i = 1, natom_z
      g5_z_xy( i + counter_natom_zxy, j ) = ( 2.0d0*( g5(i,j) - g5_z_xy_maxmin(j,2) )/ &
                                            ( g5_z_xy_maxmin(j,1) - g5_z_xy_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  allocate( g5( natom_z, num_g5_z_xz ) )
  read(11) g5
  do j = 1, num_g5_z_xz
    do i = 1, natom_z
      g5_z_xz( i + counter_natom_zx, j ) = ( 2.0d0*( g5(i,j) - g5_z_xz_maxmin(j,2) )/ &
                                           ( g5_z_xz_maxmin(j,1) - g5_z_xz_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  allocate( g5( natom_z, num_g5_z_yy ) )
  read(11) g5
  do j = 1, num_g5_z_yy
    do i = 1, natom_z
      g5_z_yy( i + counter_natom_zy, j ) = ( 2.0d0*( g5(i,j) - g5_z_yy_maxmin(j,2) )/ &
                                           ( g5_z_yy_maxmin(j,1) - g5_z_yy_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  allocate( g5( natom_z, num_g5_z_yz ) )
  read(11) g5
  do j = 1, num_g5_z_yz
    do i = 1, natom_z
      g5_z_yz( i + counter_natom_zy, j ) = ( 2.0d0*( g5(i,j) - g5_z_yz_maxmin(j,2) )/ &
                                           ( g5_z_yz_maxmin(j,1) - g5_z_yz_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  allocate( g5( natom_z, num_g5_z_zz ) )
  read(11) g5
  do j = 1, num_g5_z_zz
    do i = 1, natom_z
      g5_z_zz( i + counter_natom_z, j ) = ( 2.0d0*( g5(i,j) - g5_z_zz_maxmin(j,2) )/ &
                                          ( g5_z_zz_maxmin(j,1) - g5_z_zz_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  counter_natom_z   = counter_natom_z   + natom_z
  counter_natom_zx  = counter_natom_zx  + natom_z
  counter_natom_zy  = counter_natom_zy  + natom_z
  counter_natom_zxy = counter_natom_zxy + natom_z


  ! SF_deriv
  allocate( g2_deriv( natom_z, natom, 3, num_g2_z_x ) )
  read(12) g2_deriv
  do k = 1, num_g2_z_x
    do j = 1, natom
      do i = 1, natom_z
        m = i + natom_z * ( j - 1 )
        g2_deriv_z_x( counter_deriv_zx + m,1,k) = 2.0d0*g2_deriv(i,j,1,k)/ &
                                                  ( g2_z_x_maxmin(k,1) - g2_z_x_maxmin(k,2) ) 
        g2_deriv_z_x( counter_deriv_zx + m,2,k) = 2.0d0*g2_deriv(i,j,2,k)/ &
                                                  ( g2_z_x_maxmin(k,1) - g2_z_x_maxmin(k,2) ) 
        g2_deriv_z_x( counter_deriv_zx + m,3,k) = 2.0d0*g2_deriv(i,j,3,k)/ & 
                                                  ( g2_z_x_maxmin(k,1) - g2_z_x_maxmin(k,2) ) 
      enddo
    enddo
  enddo
  deallocate( g2_deriv )

  allocate( g2_deriv( natom_z, natom, 3, num_g2_z_y ) )
  read(12) g2_deriv
  do k = 1, num_g2_z_y
    do j = 1, natom
      do i = 1, natom_z
        m = i + natom_z * ( j - 1 )
        g2_deriv_z_y( counter_deriv_zy + m,1,k) = 2.0d0*g2_deriv(i,j,1,k)/ &
                                                  ( g2_z_y_maxmin(k,1) - g2_z_y_maxmin(k,2) ) 
        g2_deriv_z_y( counter_deriv_zy + m,2,k) = 2.0d0*g2_deriv(i,j,2,k)/ &
                                                  ( g2_z_y_maxmin(k,1) - g2_z_y_maxmin(k,2) ) 
        g2_deriv_z_y( counter_deriv_zy + m,3,k) = 2.0d0*g2_deriv(i,j,3,k)/ &
                                                  ( g2_z_y_maxmin(k,1) - g2_z_y_maxmin(k,2) ) 
      enddo
    enddo
  enddo
  deallocate( g2_deriv )

  allocate( g2_deriv( natom_z, natom, 3, num_g2_z_z ) )
  read(12) g2_deriv
  do k = 1, num_g2_z_z
    do j = 1, natom
      do i = 1, natom_z
        m = i + natom_z * ( j - 1 )
        g2_deriv_z_z( counter_deriv_z + m,1,k) = 2.0d0*g2_deriv(i,j,1,k)/ &
                                                  ( g2_z_z_maxmin(k,1) - g2_z_z_maxmin(k,2) ) 
        g2_deriv_z_z( counter_deriv_z + m,2,k) = 2.0d0*g2_deriv(i,j,2,k)/ &
                                                  ( g2_z_z_maxmin(k,1) - g2_z_z_maxmin(k,2) ) 
        g2_deriv_z_z( counter_deriv_z + m,3,k) = 2.0d0*g2_deriv(i,j,3,k)/ &
                                                  ( g2_z_z_maxmin(k,1) - g2_z_z_maxmin(k,2) ) 
      enddo
    enddo
  enddo
  deallocate( g2_deriv )

  allocate( g5_deriv( natom_z, natom, 3, num_g5_z_xx ) )
  read(12) g5_deriv
  do k = 1, num_g5_z_xx
    do j = 1, natom
      do i = 1, natom_z
        m = i + natom_z * ( j - 1 )
        g5_deriv_z_xx( counter_deriv_zx + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                  ( g5_z_xx_maxmin(k,1) - g5_z_xx_maxmin(k,2) ) 
        g5_deriv_z_xx( counter_deriv_zx + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                  ( g5_z_xx_maxmin(k,1) - g5_z_xx_maxmin(k,2) ) 
        g5_deriv_z_xx( counter_deriv_zx + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                  ( g5_z_xx_maxmin(k,1) - g5_z_xx_maxmin(k,2) ) 
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  allocate( g5_deriv( natom_z, natom, 3, num_g5_z_xy ) )
  read(12) g5_deriv
  do k = 1, num_g5_z_xy
    do j = 1, natom
      do i = 1, natom_z
        m = i + natom_z * ( j - 1 )
        g5_deriv_z_xy( counter_deriv_zxy + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                  ( g5_z_xy_maxmin(k,1) - g5_z_xy_maxmin(k,2) ) 
        g5_deriv_z_xy( counter_deriv_zxy + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                  ( g5_z_xy_maxmin(k,1) - g5_z_xy_maxmin(k,2) ) 
        g5_deriv_z_xy( counter_deriv_zxy + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                  ( g5_z_xy_maxmin(k,1) - g5_z_xy_maxmin(k,2) ) 
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  allocate( g5_deriv( natom_z, natom, 3, num_g5_z_xz ) )
  read(12) g5_deriv
  do k = 1, num_g5_z_xz
    do j = 1, natom
      do i = 1, natom_z
        m = i + natom_z * ( j - 1 )
        g5_deriv_z_xz( counter_deriv_zx + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                  ( g5_z_xz_maxmin(k,1) - g5_z_xz_maxmin(k,2) ) 
        g5_deriv_z_xz( counter_deriv_zx + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                  ( g5_z_xz_maxmin(k,1) - g5_z_xz_maxmin(k,2) ) 
        g5_deriv_z_xz( counter_deriv_zx + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                  ( g5_z_xz_maxmin(k,1) - g5_z_xz_maxmin(k,2) ) 
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  allocate( g5_deriv( natom_z, natom, 3, num_g5_z_yy ) )
  read(12) g5_deriv
  do k = 1, num_g5_z_yy
    do j = 1, natom
      do i = 1, natom_z
        m = i + natom_z * ( j - 1 )
        g5_deriv_z_yy( counter_deriv_zy + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                  ( g5_z_yy_maxmin(k,1) - g5_z_yy_maxmin(k,2) ) 
        g5_deriv_z_yy( counter_deriv_zy + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                  ( g5_z_yy_maxmin(k,1) - g5_z_yy_maxmin(k,2) ) 
        g5_deriv_z_yy( counter_deriv_zy + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                  ( g5_z_yy_maxmin(k,1) - g5_z_yy_maxmin(k,2) ) 
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  allocate( g5_deriv( natom_z, natom, 3, num_g5_z_yz ) )
  read(12) g5_deriv
  do k = 1, num_g5_z_yz
    do j = 1, natom
      do i = 1, natom_z
        m = i + natom_z * ( j - 1 )
        g5_deriv_z_yz( counter_deriv_zy + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                  ( g5_z_yz_maxmin(k,1) - g5_z_yz_maxmin(k,2) ) 
        g5_deriv_z_yz( counter_deriv_zy + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                  ( g5_z_yz_maxmin(k,1) - g5_z_yz_maxmin(k,2) ) 
        g5_deriv_z_yz( counter_deriv_zy + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                  ( g5_z_yz_maxmin(k,1) - g5_z_yz_maxmin(k,2) ) 
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  allocate( g5_deriv( natom_z, natom, 3, num_g5_z_zz ) )
  read(12) g5_deriv
  do k = 1, num_g5_z_zz
    do j = 1, natom
      do i = 1, natom_z
        m = i + natom_z * ( j - 1 )
        g5_deriv_z_zz( counter_deriv_z + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                  ( g5_z_zz_maxmin(k,1) - g5_z_zz_maxmin(k,2) ) 
        g5_deriv_z_zz( counter_deriv_z + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                  ( g5_z_zz_maxmin(k,1) - g5_z_zz_maxmin(k,2) ) 
        g5_deriv_z_zz( counter_deriv_z + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                  ( g5_z_zz_maxmin(k,1) - g5_z_zz_maxmin(k,2) ) 
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  counter_deriv_z   = counter_deriv_z   + natom_z * natom
  counter_deriv_zx  = counter_deriv_zx  + natom_z * natom
  counter_deriv_zy  = counter_deriv_zy  + natom_z * natom
  counter_deriv_zxy = counter_deriv_zxy + natom_z * natom


end subroutine read_data_3
