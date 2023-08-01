module allocarray
contains


!************************************************************
! ALLOCATE_1
!************************************************************
subroutine allocate_1( &
           tag2_x_x,  g2_x_x,  eta2_x_x,  rs2_x_x, &
           tag5_x_xx, g5_x_xx, eta5_x_xx, theta5_x_xx, zeta5_x_xx, lambda5_x_xx, &
           num_g2_x_x, num_g5_x_xx, &
           natom, natom_x, &
           g2_deriv_x_x, g5_deriv_x_xx )

implicit none
integer natom, natom_x
character(len=20) sf_name
character(len=6) tag2_x_x
character(len=7) tag5_x_xx
integer num_g2_x_x, num_g5_x_xx
double precision,allocatable,dimension(:,:) :: &
  g2_x_x, g5_x_xx
double precision,allocatable,dimension(:) :: &
  eta2_x_x,  rs2_x_x, &
  eta5_x_xx, theta5_x_xx, zeta5_x_xx, lambda5_x_xx
double precision,allocatable,dimension(:,:,:,:) :: &
  g2_deriv_x_x, g5_deriv_x_xx


 allocate(   g2_x_x( natom_x, num_g2_x_x ), &
           eta2_x_x( num_g2_x_x ), &
            rs2_x_x( num_g2_x_x ) )

 allocate(      g5_x_xx( natom_x, num_g5_x_xx ), &
              eta5_x_xx( num_g5_x_xx ), &
            theta5_x_xx( num_g5_x_xx ), &
             zeta5_x_xx( num_g5_x_xx ), &
           lambda5_x_xx( num_g5_x_xx ) )

 sf_name = trim( adjustl( tag2_x_x ) )
 call param_read_g2( sf_name, eta2_x_x, rs2_x_x, num_g2_x_x )

 sf_name = trim( adjustl( tag5_x_xx ) )
 call param_read_g5( sf_name, eta5_x_xx, theta5_x_xx, zeta5_x_xx, lambda5_x_xx, num_g5_x_xx )

 allocate( g2_deriv_x_x(  natom_x, natom, 3, num_g2_x_x ) )
 allocate( g5_deriv_x_xx( natom_x, natom, 3, num_g5_x_xx ) )


end subroutine allocate_1


!************************************************************
! ALLOCATE_2
!************************************************************
subroutine allocate_2( &
           tag2_x_x, g2_x_x, eta2_x_x, rs2_x_x, &
           tag2_x_y, g2_x_y, eta2_x_y, rs2_x_y, &
           tag5_x_xx, g5_x_xx, eta5_x_xx, theta5_x_xx, zeta5_x_xx, lambda5_x_xx, &
           tag5_x_xy, g5_x_xy, eta5_x_xy, theta5_x_xy, zeta5_x_xy, lambda5_x_xy, &
           tag5_x_yy, g5_x_yy, eta5_x_yy, theta5_x_yy, zeta5_x_yy, lambda5_x_yy, &
           num_g2_x_x, num_g2_x_y, num_g5_x_xx, num_g5_x_xy, num_g5_x_yy, &
           natom, natom_x, &
           g2_deriv_x_x, g2_deriv_x_y, g5_deriv_x_xx, g5_deriv_x_xy, g5_deriv_x_yy )

implicit none
integer natom, natom_x
character(len=20) sf_name
character(len=6) tag2_x_x, tag2_x_y
character(len=7) tag5_x_xx, tag5_x_xy, tag5_x_yy
integer num_g2_x_x, num_g2_x_y
integer num_g5_x_xx, num_g5_x_xy, num_g5_x_yy
double precision,allocatable,dimension(:,:) :: &
  g2_x_x, g2_x_y, &
  g5_x_xx, g5_x_xy, g5_x_yy
double precision,allocatable,dimension(:) :: &
  eta2_x_x, rs2_x_x, &
  eta2_x_y, rs2_x_y, &
  eta5_x_xx, theta5_x_xx, zeta5_x_xx, lambda5_x_xx, &
  eta5_x_xy, theta5_x_xy, zeta5_x_xy, lambda5_x_xy, &
  eta5_x_yy, theta5_x_yy, zeta5_x_yy, lambda5_x_yy
double precision,allocatable,dimension(:,:,:,:) :: &
  g2_deriv_x_x, g2_deriv_x_y, &
  g5_deriv_x_xx, g5_deriv_x_xy, g5_deriv_x_yy


 allocate(   g2_x_x( natom_x, num_g2_x_x ), &
           eta2_x_x( num_g2_x_x ), &
            rs2_x_x( num_g2_x_x ) )
 allocate(   g2_x_y( natom_x, num_g2_x_y ), &
           eta2_x_y( num_g2_x_y ), &
            rs2_x_y( num_g2_x_y ) )

 allocate(      g5_x_xx( natom_x, num_g5_x_xx ), &
              eta5_x_xx( num_g5_x_xx ), &
            theta5_x_xx( num_g5_x_xx ), &
             zeta5_x_xx( num_g5_x_xx ), &
           lambda5_x_xx( num_g5_x_xx ) )
 allocate(      g5_x_xy( natom_x, num_g5_x_xy ), &
              eta5_x_xy( num_g5_x_xy ), &
            theta5_x_xy( num_g5_x_xy ), &
             zeta5_x_xy( num_g5_x_xy ), &
           lambda5_x_xy( num_g5_x_xy ) )
 allocate(      g5_x_yy( natom_x, num_g5_x_yy ), &
              eta5_x_yy( num_g5_x_yy ), &
            theta5_x_yy( num_g5_x_yy ), &
             zeta5_x_yy( num_g5_x_yy ), &
           lambda5_x_yy( num_g5_x_yy ) )

 sf_name = trim(adjustl( tag2_x_x ))
 call param_read_g2( sf_name, eta2_x_x, rs2_x_x, num_g2_x_x )
 sf_name = trim(adjustl( tag2_x_y ))
 call param_read_g2( sf_name, eta2_x_y, rs2_x_y, num_g2_x_y )

 sf_name = trim(adjustl( tag5_x_xx ))
 call param_read_g5( sf_name, eta5_x_xx, theta5_x_xx, zeta5_x_xx, lambda5_x_xx, num_g5_x_xx )
 sf_name = trim(adjustl( tag5_x_xy ))
 call param_read_g5( sf_name, eta5_x_xy, theta5_x_xy, zeta5_x_xy, lambda5_x_xy, num_g5_x_xy )
 sf_name = trim(adjustl( tag5_x_yy ))
 call param_read_g5( sf_name, eta5_x_yy, theta5_x_yy, zeta5_x_yy, lambda5_x_yy, num_g5_x_yy )

 allocate( g2_deriv_x_x( natom_x, natom, 3, num_g2_x_x ) )
 allocate( g2_deriv_x_y( natom_x, natom, 3, num_g2_x_y ) )

 allocate( g5_deriv_x_xx( natom_x, natom, 3, num_g5_x_xx ) )
 allocate( g5_deriv_x_xy( natom_x, natom, 3, num_g5_x_xy ) )
 allocate( g5_deriv_x_yy( natom_x, natom, 3, num_g5_x_yy ) )


end subroutine allocate_2


!************************************************************
! ALLOCATE_3
!************************************************************
subroutine allocate_3( &
           tag2_x_x, g2_x_x, eta2_x_x, rs2_x_x, &
           tag2_x_y, g2_x_y, eta2_x_y, rs2_x_y, &
           tag2_x_z, g2_x_z, eta2_x_z, rs2_x_z, &
           tag5_x_xx, g5_x_xx, eta5_x_xx, theta5_x_xx, zeta5_x_xx, lambda5_x_xx, &
           tag5_x_xy, g5_x_xy, eta5_x_xy, theta5_x_xy, zeta5_x_xy, lambda5_x_xy, &
           tag5_x_xz, g5_x_xz, eta5_x_xz, theta5_x_xz, zeta5_x_xz, lambda5_x_xz, &
           tag5_x_yy, g5_x_yy, eta5_x_yy, theta5_x_yy, zeta5_x_yy, lambda5_x_yy, &
           tag5_x_yz, g5_x_yz, eta5_x_yz, theta5_x_yz, zeta5_x_yz, lambda5_x_yz, &
           tag5_x_zz, g5_x_zz, eta5_x_zz, theta5_x_zz, zeta5_x_zz, lambda5_x_zz, &
           num_g2_x_x, num_g2_x_y, num_g2_x_z, &
           num_g5_x_xx, num_g5_x_xy, num_g5_x_xz, &
           num_g5_x_yy, num_g5_x_yz, &
           num_g5_x_zz, &
           natom, natom_x, &
           g2_deriv_x_x, g2_deriv_x_y, g2_deriv_x_z, &
           g5_deriv_x_xx, g5_deriv_x_xy, g5_deriv_x_xz, &
           g5_deriv_x_yy, g5_deriv_x_yz, &
           g5_deriv_x_zz )

implicit none
integer natom, natom_x
character(len=20) sf_name
character(len=6) tag2_x_x, tag2_x_y, tag2_x_z
character(len=7) tag5_x_xx, tag5_x_xy, tag5_x_xz, &
                 tag5_x_yy, tag5_x_yz, &
                 tag5_x_zz
integer num_g2_x_x, num_g2_x_y, num_g2_x_z
integer num_g5_x_xx, num_g5_x_xy, num_g5_x_xz, &
        num_g5_x_yy, num_g5_x_yz, &
        num_g5_x_zz
double precision,allocatable,dimension(:,:) :: &
  g2_x_x, g2_x_y, g2_x_z, &
  g5_x_xx, g5_x_xy, g5_x_xz, g5_x_yy, g5_x_yz, g5_x_zz
double precision,allocatable,dimension(:) :: &
  eta2_x_x, rs2_x_x, &
  eta2_x_y, rs2_x_y, &
  eta2_x_z, rs2_x_z, &
  eta5_x_xx, theta5_x_xx, zeta5_x_xx, lambda5_x_xx, &
  eta5_x_xy, theta5_x_xy, zeta5_x_xy, lambda5_x_xy, &
  eta5_x_xz, theta5_x_xz, zeta5_x_xz, lambda5_x_xz, &
  eta5_x_yy, theta5_x_yy, zeta5_x_yy, lambda5_x_yy, &
  eta5_x_yz, theta5_x_yz, zeta5_x_yz, lambda5_x_yz, &
  eta5_x_zz, theta5_x_zz, zeta5_x_zz, lambda5_x_zz
double precision,allocatable,dimension(:,:,:,:) :: &
  g2_deriv_x_x, g2_deriv_x_y, g2_deriv_x_z, &
  g5_deriv_x_xx, g5_deriv_x_xy, g5_deriv_x_xz, &
  g5_deriv_x_yy, g5_deriv_x_yz, &
  g5_deriv_x_zz


 allocate(   g2_x_x( natom_x, num_g2_x_x ), &
           eta2_x_x( num_g2_x_x ), &
            rs2_x_x( num_g2_x_x ) )
 allocate(   g2_x_y( natom_x, num_g2_x_y ), &
           eta2_x_y( num_g2_x_y ), &
            rs2_x_y( num_g2_x_y ) )
 allocate(   g2_x_z( natom_x, num_g2_x_z ), &
           eta2_x_z( num_g2_x_z ), &
            rs2_x_z( num_g2_x_z ) )

 allocate(      g5_x_xx( natom_x, num_g5_x_xx ), &
              eta5_x_xx( num_g5_x_xx ), &
            theta5_x_xx( num_g5_x_xx ), &
             zeta5_x_xx( num_g5_x_xx ), &
           lambda5_x_xx( num_g5_x_xx ) )
 allocate(      g5_x_xy( natom_x, num_g5_x_xy ), &
              eta5_x_xy( num_g5_x_xy ), &
            theta5_x_xy( num_g5_x_xy ), &
             zeta5_x_xy( num_g5_x_xy ), &
           lambda5_x_xy( num_g5_x_xy ) )
 allocate(      g5_x_xz( natom_x, num_g5_x_xz ), &
              eta5_x_xz( num_g5_x_xz ), &
            theta5_x_xz( num_g5_x_xz ), &
             zeta5_x_xz( num_g5_x_xz ), &
           lambda5_x_xz( num_g5_x_xz ) )
 allocate(      g5_x_yy( natom_x, num_g5_x_yy ), &
              eta5_x_yy( num_g5_x_yy ), &
            theta5_x_yy( num_g5_x_yy ), &
             zeta5_x_yy( num_g5_x_yy ), &
           lambda5_x_yy( num_g5_x_yy ) )
 allocate(      g5_x_yz( natom_x, num_g5_x_yz ), &
              eta5_x_yz( num_g5_x_yz ), &
            theta5_x_yz( num_g5_x_yz ), &
             zeta5_x_yz( num_g5_x_yz ), &
           lambda5_x_yz( num_g5_x_yz ) )
 allocate(      g5_x_zz( natom_x, num_g5_x_zz ), &
              eta5_x_zz( num_g5_x_zz ), &
            theta5_x_zz( num_g5_x_zz ), &
             zeta5_x_zz( num_g5_x_zz ), &
           lambda5_x_zz( num_g5_x_zz ) )

 sf_name = trim(adjustl( tag2_x_x ))
 call param_read_g2( sf_name, eta2_x_x, rs2_x_x, num_g2_x_x )
 sf_name = trim(adjustl( tag2_x_y ))
 call param_read_g2( sf_name, eta2_x_y, rs2_x_y, num_g2_x_y )
 sf_name = trim(adjustl( tag2_x_z ))
 call param_read_g2( sf_name, eta2_x_z, rs2_x_z, num_g2_x_z )

 sf_name = trim(adjustl( tag5_x_xx ))
 call param_read_g5( sf_name, eta5_x_xx, theta5_x_xx, zeta5_x_xx, lambda5_x_xx, num_g5_x_xx )
 sf_name = trim(adjustl( tag5_x_xy ))
 call param_read_g5( sf_name, eta5_x_xy, theta5_x_xy, zeta5_x_xy, lambda5_x_xy, num_g5_x_xy )
 sf_name = trim(adjustl( tag5_x_xz ))
 call param_read_g5( sf_name, eta5_x_xz, theta5_x_xz, zeta5_x_xz, lambda5_x_xz, num_g5_x_xz )
 sf_name = trim(adjustl( tag5_x_yy ))
 call param_read_g5( sf_name, eta5_x_yy, theta5_x_yy, zeta5_x_yy, lambda5_x_yy, num_g5_x_yy )
 sf_name = trim(adjustl( tag5_x_yz ))
 call param_read_g5( sf_name, eta5_x_yz, theta5_x_yz, zeta5_x_yz, lambda5_x_yz, num_g5_x_yz )
 sf_name = trim(adjustl( tag5_x_zz ))
 call param_read_g5( sf_name, eta5_x_zz, theta5_x_zz, zeta5_x_zz, lambda5_x_zz, num_g5_x_zz )

 allocate( g2_deriv_x_x( natom_x, natom, 3, num_g2_x_x ) )
 allocate( g2_deriv_x_y( natom_x, natom, 3, num_g2_x_y ) )
 allocate( g2_deriv_x_z( natom_x, natom, 3, num_g2_x_z ) )

 allocate( g5_deriv_x_xx( natom_x, natom, 3, num_g5_x_xx ) )
 allocate( g5_deriv_x_xy( natom_x, natom, 3, num_g5_x_xy ) )
 allocate( g5_deriv_x_xz( natom_x, natom, 3, num_g5_x_xz ) )
 allocate( g5_deriv_x_yy( natom_x, natom, 3, num_g5_x_yy ) )
 allocate( g5_deriv_x_yz( natom_x, natom, 3, num_g5_x_yz ) )
 allocate( g5_deriv_x_zz( natom_x, natom, 3, num_g5_x_zz ) )


end subroutine allocate_3


!************************************************************
! ALLOCATE_4
!************************************************************
subroutine allocate_4( &
           tag2_x_x,  g2_x_x,  eta2_x_x,  rs2_x_x, &
           tag2_x_y,  g2_x_y,  eta2_x_y,  rs2_x_y, &
           tag2_x_z,  g2_x_z,  eta2_x_z,  rs2_x_z, &
           tag2_x_v,  g2_x_v,  eta2_x_v,  rs2_x_v, &
           tag5_x_xx, g5_x_xx, eta5_x_xx, theta5_x_xx, zeta5_x_xx, lambda5_x_xx, &
           tag5_x_xy, g5_x_xy, eta5_x_xy, theta5_x_xy, zeta5_x_xy, lambda5_x_xy, &
           tag5_x_xz, g5_x_xz, eta5_x_xz, theta5_x_xz, zeta5_x_xz, lambda5_x_xz, &
           tag5_x_xv, g5_x_xv, eta5_x_xv, theta5_x_xv, zeta5_x_xv, lambda5_x_xv, &
           tag5_x_yy, g5_x_yy, eta5_x_yy, theta5_x_yy, zeta5_x_yy, lambda5_x_yy, &
           tag5_x_yz, g5_x_yz, eta5_x_yz, theta5_x_yz, zeta5_x_yz, lambda5_x_yz, &
           tag5_x_yv, g5_x_yv, eta5_x_yv, theta5_x_yv, zeta5_x_yv, lambda5_x_yv, &
           tag5_x_zz, g5_x_zz, eta5_x_zz, theta5_x_zz, zeta5_x_zz, lambda5_x_zz, &
           tag5_x_zv, g5_x_zv, eta5_x_zv, theta5_x_zv, zeta5_x_zv, lambda5_x_zv, &
           tag5_x_vv, g5_x_vv, eta5_x_vv, theta5_x_vv, zeta5_x_vv, lambda5_x_vv, &
           num_g2_x_x,  num_g2_x_y,  num_g2_x_z,  num_g2_x_v, &
           num_g5_x_xx, num_g5_x_xy, num_g5_x_xz, num_g5_x_xv, &
           num_g5_x_yy, num_g5_x_yz, num_g5_x_yv, &
           num_g5_x_zz, num_g5_x_zv, &
           num_g5_x_vv, &
           natom, natom_x, &
           g2_deriv_x_x,  g2_deriv_x_y,  g2_deriv_x_z,  g2_deriv_x_v, &
           g5_deriv_x_xx, g5_deriv_x_xy, g5_deriv_x_xz, g5_deriv_x_xv, &
           g5_deriv_x_yy, g5_deriv_x_yz, g5_deriv_x_yv, &
           g5_deriv_x_zz, g5_deriv_x_zv, &
           g5_deriv_x_vv )

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
  g2_x_x,  g2_x_y,  g2_x_z,  g2_x_v, &
  g5_x_xx, g5_x_xy, g5_x_xz, g5_x_xv, &
  g5_x_yy, g5_x_yz, g5_x_yv, &
  g5_x_zz, g5_x_zv, &
  g5_x_vv
double precision,allocatable,dimension(:) :: &
  eta2_x_x,  rs2_x_x, &
  eta2_x_y,  rs2_x_y, &
  eta2_x_z,  rs2_x_z, &
  eta2_x_v,  rs2_x_v, &
  eta5_x_xx, theta5_x_xx, zeta5_x_xx, lambda5_x_xx, &
  eta5_x_xy, theta5_x_xy, zeta5_x_xy, lambda5_x_xy, &
  eta5_x_xz, theta5_x_xz, zeta5_x_xz, lambda5_x_xz, &
  eta5_x_xv, theta5_x_xv, zeta5_x_xv, lambda5_x_xv, &
  eta5_x_yy, theta5_x_yy, zeta5_x_yy, lambda5_x_yy, &
  eta5_x_yz, theta5_x_yz, zeta5_x_yz, lambda5_x_yz, &
  eta5_x_yv, theta5_x_yv, zeta5_x_yv, lambda5_x_yv, &
  eta5_x_zz, theta5_x_zz, zeta5_x_zz, lambda5_x_zz, &
  eta5_x_zv, theta5_x_zv, zeta5_x_zv, lambda5_x_zv, &
  eta5_x_vv, theta5_x_vv, zeta5_x_vv, lambda5_x_vv
double precision,allocatable,dimension(:,:,:,:) :: &
  g2_deriv_x_x,  g2_deriv_x_y,  g2_deriv_x_z,  g2_deriv_x_v, &
  g5_deriv_x_xx, g5_deriv_x_xy, g5_deriv_x_xz, g5_deriv_x_xv, &
  g5_deriv_x_yy, g5_deriv_x_yz, g5_deriv_x_yv, &
  g5_deriv_x_zz, g5_deriv_x_zv, &
  g5_deriv_x_vv


 allocate(   g2_x_x( natom_x, num_g2_x_x ), &
           eta2_x_x( num_g2_x_x ), &
            rs2_x_x( num_g2_x_x ) )
 allocate(   g2_x_y( natom_x, num_g2_x_y ), &
           eta2_x_y( num_g2_x_y ), &
            rs2_x_y( num_g2_x_y ) )
 allocate(   g2_x_z( natom_x, num_g2_x_z ), &
           eta2_x_z( num_g2_x_z ), &
            rs2_x_z( num_g2_x_z ) )
 allocate(   g2_x_v( natom_x, num_g2_x_v ), &
           eta2_x_v( num_g2_x_v ), &
            rs2_x_v( num_g2_x_v ) )

 allocate(      g5_x_xx( natom_x, num_g5_x_xx ), &
              eta5_x_xx( num_g5_x_xx ), &
            theta5_x_xx( num_g5_x_xx ), &
             zeta5_x_xx( num_g5_x_xx ), &
           lambda5_x_xx( num_g5_x_xx ) )
 allocate(      g5_x_xy( natom_x, num_g5_x_xy ), &
              eta5_x_xy( num_g5_x_xy ), &
            theta5_x_xy( num_g5_x_xy ), &
             zeta5_x_xy( num_g5_x_xy ), &
           lambda5_x_xy( num_g5_x_xy ) )
 allocate(      g5_x_xz( natom_x, num_g5_x_xz ), &
              eta5_x_xz( num_g5_x_xz ), &
            theta5_x_xz( num_g5_x_xz ), &
             zeta5_x_xz( num_g5_x_xz ), &
           lambda5_x_xz( num_g5_x_xz ) )
 allocate(      g5_x_xv( natom_x, num_g5_x_xv ), &
              eta5_x_xv( num_g5_x_xv ), &
            theta5_x_xv( num_g5_x_xv ), &
             zeta5_x_xv( num_g5_x_xv ), &
           lambda5_x_xv( num_g5_x_xv ) )
 allocate(      g5_x_yy( natom_x, num_g5_x_yy ), &
              eta5_x_yy( num_g5_x_yy ), &
            theta5_x_yy( num_g5_x_yy ), &
             zeta5_x_yy( num_g5_x_yy ), &
           lambda5_x_yy( num_g5_x_yy ) )
 allocate(      g5_x_yz( natom_x, num_g5_x_yz ), &
              eta5_x_yz( num_g5_x_yz ), &
            theta5_x_yz( num_g5_x_yz ), &
             zeta5_x_yz( num_g5_x_yz ), &
           lambda5_x_yz( num_g5_x_yz ) )
 allocate(      g5_x_yv( natom_x, num_g5_x_yv ), &
              eta5_x_yv( num_g5_x_yv ), &
            theta5_x_yv( num_g5_x_yv ), &
             zeta5_x_yv( num_g5_x_yv ), &
           lambda5_x_yv( num_g5_x_yv ) )
 allocate(      g5_x_zz( natom_x, num_g5_x_zz ), &
              eta5_x_zz( num_g5_x_zz ), &
            theta5_x_zz( num_g5_x_zz ), &
             zeta5_x_zz( num_g5_x_zz ), &
           lambda5_x_zz( num_g5_x_zz ) )
 allocate(      g5_x_zv( natom_x, num_g5_x_zv ), &
              eta5_x_zv( num_g5_x_zv ), &
            theta5_x_zv( num_g5_x_zv ), &
             zeta5_x_zv( num_g5_x_zv ), &
           lambda5_x_zv( num_g5_x_zv ) )
 allocate(      g5_x_vv( natom_x, num_g5_x_vv ), &
              eta5_x_vv( num_g5_x_vv ), &
            theta5_x_vv( num_g5_x_vv ), &
             zeta5_x_vv( num_g5_x_vv ), &
           lambda5_x_vv( num_g5_x_vv ) )

 sf_name = trim(adjustl( tag2_x_x ))
 call param_read_g2( sf_name, eta2_x_x, rs2_x_x, num_g2_x_x )
 sf_name = trim(adjustl( tag2_x_y ))
 call param_read_g2( sf_name, eta2_x_y, rs2_x_y, num_g2_x_y )
 sf_name = trim(adjustl( tag2_x_z ))
 call param_read_g2( sf_name, eta2_x_z, rs2_x_z, num_g2_x_z )
 sf_name = trim(adjustl( tag2_x_v ))
 call param_read_g2( sf_name, eta2_x_v, rs2_x_v, num_g2_x_v )

 sf_name = trim(adjustl( tag5_x_xx ))
 call param_read_g5( sf_name, eta5_x_xx, theta5_x_xx, zeta5_x_xx, lambda5_x_xx, num_g5_x_xx )
 sf_name = trim(adjustl( tag5_x_xy ))
 call param_read_g5( sf_name, eta5_x_xy, theta5_x_xy, zeta5_x_xy, lambda5_x_xy, num_g5_x_xy )
 sf_name = trim(adjustl( tag5_x_xz ))
 call param_read_g5( sf_name, eta5_x_xz, theta5_x_xz, zeta5_x_xz, lambda5_x_xz, num_g5_x_xz )
 sf_name = trim(adjustl( tag5_x_xv ))
 call param_read_g5( sf_name, eta5_x_xv, theta5_x_xv, zeta5_x_xv, lambda5_x_xv, num_g5_x_xv )
 sf_name = trim(adjustl( tag5_x_yy ))
 call param_read_g5( sf_name, eta5_x_yy, theta5_x_yy, zeta5_x_yy, lambda5_x_yy, num_g5_x_yy )
 sf_name = trim(adjustl( tag5_x_yz ))
 call param_read_g5( sf_name, eta5_x_yz, theta5_x_yz, zeta5_x_yz, lambda5_x_yz, num_g5_x_yz )
 sf_name = trim(adjustl( tag5_x_yv ))
 call param_read_g5( sf_name, eta5_x_yv, theta5_x_yv, zeta5_x_yv, lambda5_x_yv, num_g5_x_yv )
 sf_name = trim(adjustl( tag5_x_zz ))
 call param_read_g5( sf_name, eta5_x_zz, theta5_x_zz, zeta5_x_zz, lambda5_x_zz, num_g5_x_zz )
 sf_name = trim(adjustl( tag5_x_zv ))
 call param_read_g5( sf_name, eta5_x_zv, theta5_x_zv, zeta5_x_zv, lambda5_x_zv, num_g5_x_zv )
 sf_name = trim(adjustl( tag5_x_vv ))
 call param_read_g5( sf_name, eta5_x_vv, theta5_x_vv, zeta5_x_vv, lambda5_x_vv, num_g5_x_vv )

 allocate( g2_deriv_x_x( natom_x, natom, 3, num_g2_x_x ) )
 allocate( g2_deriv_x_y( natom_x, natom, 3, num_g2_x_y ) )
 allocate( g2_deriv_x_z( natom_x, natom, 3, num_g2_x_z ) )
 allocate( g2_deriv_x_v( natom_x, natom, 3, num_g2_x_v ) )

 allocate( g5_deriv_x_xx( natom_x, natom, 3, num_g5_x_xx ) )
 allocate( g5_deriv_x_xy( natom_x, natom, 3, num_g5_x_xy ) )
 allocate( g5_deriv_x_xz( natom_x, natom, 3, num_g5_x_xz ) )
 allocate( g5_deriv_x_xv( natom_x, natom, 3, num_g5_x_xv ) )
 allocate( g5_deriv_x_yy( natom_x, natom, 3, num_g5_x_yy ) )
 allocate( g5_deriv_x_yz( natom_x, natom, 3, num_g5_x_yz ) )
 allocate( g5_deriv_x_yv( natom_x, natom, 3, num_g5_x_yv ) )
 allocate( g5_deriv_x_zz( natom_x, natom, 3, num_g5_x_zz ) )
 allocate( g5_deriv_x_zv( natom_x, natom, 3, num_g5_x_zv ) )
 allocate( g5_deriv_x_vv( natom_x, natom, 3, num_g5_x_vv ) )


end subroutine allocate_4


end module
