subroutine sf_4( posi, xlattice, ylattice, zlattice, &
                 g2_x_x,  eta2_x_x,  rs2_x_x, &
                 g2_x_y,  eta2_x_y,  rs2_x_y, &
                 g2_x_z,  eta2_x_z,  rs2_x_z, &
                 g2_x_v,  eta2_x_v,  rs2_x_v, &
                 g5_x_xx, eta5_x_xx, theta5_x_xx, zeta5_x_xx, lambda5_x_xx, &
                 g5_x_xy, eta5_x_xy, theta5_x_xy, zeta5_x_xy, lambda5_x_xy, &
                 g5_x_xz, eta5_x_xz, theta5_x_xz, zeta5_x_xz, lambda5_x_xz, &
                 g5_x_xv, eta5_x_xv, theta5_x_xv, zeta5_x_xv, lambda5_x_xv, &
                 g5_x_yy, eta5_x_yy, theta5_x_yy, zeta5_x_yy, lambda5_x_yy, &
                 g5_x_yz, eta5_x_yz, theta5_x_yz, zeta5_x_yz, lambda5_x_yz, &
                 g5_x_yv, eta5_x_yv, theta5_x_yv, zeta5_x_yv, lambda5_x_yv, &
                 g5_x_zz, eta5_x_zz, theta5_x_zz, zeta5_x_zz, lambda5_x_zz, &
                 g5_x_zv, eta5_x_zv, theta5_x_zv, zeta5_x_zv, lambda5_x_zv, &
                 g5_x_vv, eta5_x_vv, theta5_x_vv, zeta5_x_vv, lambda5_x_vv, &
                 g2_y_x,  eta2_y_x,  rs2_y_x, &
                 g2_y_y,  eta2_y_y,  rs2_y_y, &
                 g2_y_z,  eta2_y_z,  rs2_y_z, &
                 g2_y_v,  eta2_y_v,  rs2_y_v, &
                 g5_y_xx, eta5_y_xx, theta5_y_xx, zeta5_y_xx, lambda5_y_xx, &
                 g5_y_xy, eta5_y_xy, theta5_y_xy, zeta5_y_xy, lambda5_y_xy, &
                 g5_y_xz, eta5_y_xz, theta5_y_xz, zeta5_y_xz, lambda5_y_xz, &
                 g5_y_xv, eta5_y_xv, theta5_y_xv, zeta5_y_xv, lambda5_y_xv, &
                 g5_y_yy, eta5_y_yy, theta5_y_yy, zeta5_y_yy, lambda5_y_yy, &
                 g5_y_yz, eta5_y_yz, theta5_y_yz, zeta5_y_yz, lambda5_y_yz, &
                 g5_y_yv, eta5_y_yv, theta5_y_yv, zeta5_y_yv, lambda5_y_yv, &
                 g5_y_zz, eta5_y_zz, theta5_y_zz, zeta5_y_zz, lambda5_y_zz, &
                 g5_y_zv, eta5_y_zv, theta5_y_zv, zeta5_y_zv, lambda5_y_zv, &
                 g5_y_vv, eta5_y_vv, theta5_y_vv, zeta5_y_vv, lambda5_y_vv, &
                 g2_z_x,  eta2_z_x,  rs2_z_x, &
                 g2_z_y,  eta2_z_y,  rs2_z_y, &
                 g2_z_z,  eta2_z_z,  rs2_z_z, &
                 g2_z_v,  eta2_z_v,  rs2_z_v, &
                 g5_z_xx, eta5_z_xx, theta5_z_xx, zeta5_z_xx, lambda5_z_xx, &
                 g5_z_xy, eta5_z_xy, theta5_z_xy, zeta5_z_xy, lambda5_z_xy, &
                 g5_z_xz, eta5_z_xz, theta5_z_xz, zeta5_z_xz, lambda5_z_xz, &
                 g5_z_xv, eta5_z_xv, theta5_z_xv, zeta5_z_xv, lambda5_z_xv, &
                 g5_z_yy, eta5_z_yy, theta5_z_yy, zeta5_z_yy, lambda5_z_yy, &
                 g5_z_yz, eta5_z_yz, theta5_z_yz, zeta5_z_yz, lambda5_z_yz, &
                 g5_z_yv, eta5_z_yv, theta5_z_yv, zeta5_z_yv, lambda5_z_yv, &
                 g5_z_zz, eta5_z_zz, theta5_z_zz, zeta5_z_zz, lambda5_z_zz, &
                 g5_z_zv, eta5_z_zv, theta5_z_zv, zeta5_z_zv, lambda5_z_zv, &
                 g5_z_vv, eta5_z_vv, theta5_z_vv, zeta5_z_vv, lambda5_z_vv, &
                 g2_v_x,  eta2_v_x,  rs2_v_x, &
                 g2_v_y,  eta2_v_y,  rs2_v_y, &
                 g2_v_z,  eta2_v_z,  rs2_v_z, &
                 g2_v_v,  eta2_v_v,  rs2_v_v, &
                 g5_v_xx, eta5_v_xx, theta5_v_xx, zeta5_v_xx, lambda5_v_xx, &
                 g5_v_xy, eta5_v_xy, theta5_v_xy, zeta5_v_xy, lambda5_v_xy, &
                 g5_v_xz, eta5_v_xz, theta5_v_xz, zeta5_v_xz, lambda5_v_xz, &
                 g5_v_xv, eta5_v_xv, theta5_v_xv, zeta5_v_xv, lambda5_v_xv, &
                 g5_v_yy, eta5_v_yy, theta5_v_yy, zeta5_v_yy, lambda5_v_yy, &
                 g5_v_yz, eta5_v_yz, theta5_v_yz, zeta5_v_yz, lambda5_v_yz, &
                 g5_v_yv, eta5_v_yv, theta5_v_yv, zeta5_v_yv, lambda5_v_yv, &
                 g5_v_zz, eta5_v_zz, theta5_v_zz, zeta5_v_zz, lambda5_v_zz, &
                 g5_v_zv, eta5_v_zv, theta5_v_zv, zeta5_v_zv, lambda5_v_zv, &
                 g5_v_vv, eta5_v_vv, theta5_v_vv, zeta5_v_vv, lambda5_v_vv, &
                 g2_deriv_x_x,  g2_deriv_x_y,  g2_deriv_x_z,  g2_deriv_x_v, &
                 g5_deriv_x_xx, g5_deriv_x_xy, g5_deriv_x_xz, g5_deriv_x_xv, &
                 g5_deriv_x_yy, g5_deriv_x_yz, g5_deriv_x_yv, &
                 g5_deriv_x_zz, g5_deriv_x_zv, &
                 g5_deriv_x_vv, &
                 g2_deriv_y_x,  g2_deriv_y_y,  g2_deriv_y_z,  g2_deriv_y_v, &
                 g5_deriv_y_xx, g5_deriv_y_xy, g5_deriv_y_xz, g5_deriv_y_xv, &
                 g5_deriv_y_yy, g5_deriv_y_yz, g5_deriv_y_yv, &
                 g5_deriv_y_zz, g5_deriv_y_zv, &
                 g5_deriv_y_vv, &
                 g2_deriv_z_x,  g2_deriv_z_y,  g2_deriv_z_z,  g2_deriv_z_v, &
                 g5_deriv_z_xx, g5_deriv_z_xy, g5_deriv_z_xz, g5_deriv_z_xv, &
                 g5_deriv_z_yy, g5_deriv_z_yz, g5_deriv_z_yv, &
                 g5_deriv_z_zz, g5_deriv_z_zv, &
                 g5_deriv_z_vv, &
                 g2_deriv_v_x,  g2_deriv_v_y,  g2_deriv_v_z,  g2_deriv_v_v, &
                 g5_deriv_v_xx, g5_deriv_v_xy, g5_deriv_v_xz, g5_deriv_v_xv, &
                 g5_deriv_v_yy, g5_deriv_v_yz, g5_deriv_v_yv, &
                 g5_deriv_v_zz, g5_deriv_v_zv, &
                 g5_deriv_v_vv, &
                 num_g2_x_x,  num_g2_x_y,  num_g2_x_z,  num_g2_x_v, &
                 num_g5_x_xx, num_g5_x_xy, num_g5_x_xz, num_g5_x_xv, &
                 num_g5_x_yy, num_g5_x_yz, num_g5_x_yv, &
                 num_g5_x_zz, num_g5_x_zv, &
                 num_g5_x_vv, &
                 num_g2_y_x,  num_g2_y_y,  num_g2_y_z,  num_g2_y_v, &
                 num_g5_y_xx, num_g5_y_xy, num_g5_y_xz, num_g5_y_xv, &
                 num_g5_y_yy, num_g5_y_yz, num_g5_y_yv, &
                 num_g5_y_zz, num_g5_y_zv, &
                 num_g5_y_vv, &
                 num_g2_z_x,  num_g2_z_y,  num_g2_z_z,  num_g2_z_v, &
                 num_g5_z_xx, num_g5_z_xy, num_g5_z_xz, num_g5_z_xv, &
                 num_g5_z_yy, num_g5_z_yz, num_g5_z_yv, &
                 num_g5_z_zz, num_g5_z_zv, &
                 num_g5_z_vv, &
                 num_g2_v_x,  num_g2_v_y,  num_g2_v_z,  num_g2_v_v, &
                 num_g5_v_xx, num_g5_v_xy, num_g5_v_xz, num_g5_v_xv, &
                 num_g5_v_yy, num_g5_v_yz, num_g5_v_yv, &
                 num_g5_v_zz, num_g5_v_zv, &
                 num_g5_v_vv, &
                 natom_x, natom_y, natom_z, natom_v, &
                 natom, rc, neighbor, io_sf, num_mpi, myrank )
implicit none
include 'mpif.h'
integer ierr, num_mpi, myrank
integer i, j, k, l, m, n, ii, ig5, i0
integer ix1, iy1, iz1, ix2, lattice_box, boxsize, boxcenter, neighbor, center
integer natom, atom1, atom2 !, time, jamp
integer natom_x, natom_y, natom_z, natom_v
! natom : total number of atoms
double precision xlattice, ylattice, zlattice
double precision posi(  natom, 3 ), posi1( natom, 3 )
double precision,allocatable :: posi2( : , : ) ! posi1(natom,3)
integer :: list( natom, neighbor )
!X
integer num_g2_x_x,  num_g2_x_y,  num_g2_x_z,  num_g2_x_v, &
        num_g5_x_xx, num_g5_x_xy, num_g5_x_xz, num_g5_x_xv, &
        num_g5_x_yy, num_g5_x_yz, num_g5_x_yv, &
        num_g5_x_zz, num_g5_x_zv, &
        num_g5_x_vv
!Y
integer num_g2_y_x,  num_g2_y_y,  num_g2_y_z,  num_g2_y_v, &
        num_g5_y_xx, num_g5_y_xy, num_g5_y_xz, num_g5_y_xv, &
        num_g5_y_yy, num_g5_y_yz, num_g5_y_yv, &
        num_g5_y_zz, num_g5_y_zv, &
        num_g5_y_vv
!Z
integer num_g2_z_x,  num_g2_z_y,  num_g2_z_z,  num_g2_z_v, &
        num_g5_z_xx, num_g5_z_xy, num_g5_z_xz, num_g5_z_xv, &
        num_g5_z_yy, num_g5_z_yz, num_g5_z_yv, &
        num_g5_z_zz, num_g5_z_zv, &
        num_g5_z_vv
!V
integer num_g2_v_x,  num_g2_v_y,  num_g2_v_z,  num_g2_v_v, &
        num_g5_v_xx, num_g5_v_xy, num_g5_v_xz, num_g5_v_xv, &
        num_g5_v_yy, num_g5_v_yz, num_g5_v_yv, &
        num_g5_v_zz, num_g5_v_zv, &
        num_g5_v_vv
!X
double precision &
  g2_x_x( natom_x, num_g2_x_x ), eta2_x_x( num_g2_x_x ), rs2_x_x( num_g2_x_x ), &
  g2_x_y( natom_x, num_g2_x_y ), eta2_x_y( num_g2_x_y ), rs2_x_y( num_g2_x_y ), &
  g2_x_z( natom_x, num_g2_x_z ), eta2_x_z( num_g2_x_z ), rs2_x_z( num_g2_x_z ), &
  g2_x_v( natom_x, num_g2_x_v ), eta2_x_v( num_g2_x_v ), rs2_x_v( num_g2_x_v )
double precision &
  g5_x_xx( natom_x, num_g5_x_xx ), &
   eta5_x_xx( num_g5_x_xx ),  theta5_x_xx( num_g5_x_xx ), &
  zeta5_x_xx( num_g5_x_xx ), lambda5_x_xx( num_g5_x_xx ), &
  g5_x_xy( natom_x, num_g5_x_xy ), &
   eta5_x_xy( num_g5_x_xy ),  theta5_x_xy( num_g5_x_xy ), &
  zeta5_x_xy( num_g5_x_xy ), lambda5_x_xy( num_g5_x_xy ), &
  g5_x_xz( natom_x, num_g5_x_xz ), &
   eta5_x_xz( num_g5_x_xz ),  theta5_x_xz( num_g5_x_xz ), &
  zeta5_x_xz( num_g5_x_xz ), lambda5_x_xz( num_g5_x_xz ), &
  g5_x_xv( natom_x, num_g5_x_xv ), &
   eta5_x_xv( num_g5_x_xv ),  theta5_x_xv( num_g5_x_xv ), &
  zeta5_x_xv( num_g5_x_xv ), lambda5_x_xv( num_g5_x_xv ), &
  g5_x_yy( natom_x, num_g5_x_yy ), &
   eta5_x_yy( num_g5_x_yy ),  theta5_x_yy( num_g5_x_yy ), &
  zeta5_x_yy( num_g5_x_yy ), lambda5_x_yy( num_g5_x_yy ), &
  g5_x_yz( natom_x, num_g5_x_yz ), &
   eta5_x_yz( num_g5_x_yz ),  theta5_x_yz( num_g5_x_yz ), &
  zeta5_x_yz( num_g5_x_yz ), lambda5_x_yz( num_g5_x_yz ), &
  g5_x_yv( natom_x, num_g5_x_yv ), &
   eta5_x_yv( num_g5_x_yv ),  theta5_x_yv( num_g5_x_yv ), &
  zeta5_x_yv( num_g5_x_yv ), lambda5_x_yv( num_g5_x_yv ), &
  g5_x_zz( natom_x, num_g5_x_zz ), &
   eta5_x_zz( num_g5_x_zz ),  theta5_x_zz( num_g5_x_zz ), &
  zeta5_x_zz( num_g5_x_zz ), lambda5_x_zz( num_g5_x_zz ), &
  g5_x_zv( natom_x, num_g5_x_zv ), &
   eta5_x_zv( num_g5_x_zv ),  theta5_x_zv( num_g5_x_zv ), &
  zeta5_x_zv( num_g5_x_zv ), lambda5_x_zv( num_g5_x_zv ), &
  g5_x_vv( natom_x, num_g5_x_vv ), &
   eta5_x_vv( num_g5_x_vv ),  theta5_x_vv( num_g5_x_vv ), &
  zeta5_x_vv( num_g5_x_vv ), lambda5_x_vv( num_g5_x_vv )
double precision &
  g2_x_x_mpi( natom_x, num_g2_x_x ), &
  g2_x_y_mpi( natom_x, num_g2_x_x ), &
  g2_x_z_mpi( natom_x, num_g2_x_x ), &
  g2_x_v_mpi( natom_x, num_g2_x_x ), &
  g5_x_xx_mpi( natom_x, num_g5_x_xx ), &
  g5_x_xy_mpi( natom_x, num_g5_x_xy ), &
  g5_x_xz_mpi( natom_x, num_g5_x_xz ), &
  g5_x_xv_mpi( natom_x, num_g5_x_xv ), &
  g5_x_yy_mpi( natom_x, num_g5_x_yy ), &
  g5_x_yz_mpi( natom_x, num_g5_x_yz ), &
  g5_x_yv_mpi( natom_x, num_g5_x_yv ), &
  g5_x_zz_mpi( natom_x, num_g5_x_zz ), &
  g5_x_zv_mpi( natom_x, num_g5_x_zv ), &
  g5_x_vv_mpi( natom_x, num_g5_x_vv )
!Y
double precision &
  g2_y_x( natom_y, num_g2_y_x ), eta2_y_x( num_g2_y_x ), rs2_y_x( num_g2_y_x ), &
  g2_y_y( natom_y, num_g2_y_y ), eta2_y_y( num_g2_y_y ), rs2_y_y( num_g2_y_y ), &
  g2_y_z( natom_y, num_g2_y_z ), eta2_y_z( num_g2_y_z ), rs2_y_z( num_g2_y_z ), &
  g2_y_v( natom_y, num_g2_y_v ), eta2_y_v( num_g2_y_v ), rs2_y_v( num_g2_y_v )
double precision &
  g5_y_xx( natom_y, num_g5_y_xx ), &
   eta5_y_xx( num_g5_y_xx ),  theta5_y_xx( num_g5_y_xx ), &
  zeta5_y_xx( num_g5_y_xx ), lambda5_y_xx( num_g5_y_xx ), &
  g5_y_xy( natom_y, num_g5_y_xy ), &
   eta5_y_xy( num_g5_y_xy ),  theta5_y_xy( num_g5_y_xy ), &
  zeta5_y_xy( num_g5_y_xy ), lambda5_y_xy( num_g5_y_xy ), &
  g5_y_xz( natom_y, num_g5_y_xz ), &
   eta5_y_xz( num_g5_y_xz ),  theta5_y_xz( num_g5_y_xz ),&
  zeta5_y_xz( num_g5_y_xz ), lambda5_y_xz( num_g5_y_xz ), &
  g5_y_xv( natom_y, num_g5_y_xv ), &
   eta5_y_xv( num_g5_y_xv ),  theta5_y_xv( num_g5_y_xv ),&
  zeta5_y_xv( num_g5_y_xv ), lambda5_y_xv( num_g5_y_xv ), &
  g5_y_yy( natom_y, num_g5_y_yy ), &
   eta5_y_yy( num_g5_y_yy ),  theta5_y_yy( num_g5_y_yy ), &
  zeta5_y_yy( num_g5_y_yy ), lambda5_y_yy( num_g5_y_yy ), &
  g5_y_yz( natom_y, num_g5_y_yz ), &
   eta5_y_yz( num_g5_y_yz ),  theta5_y_yz( num_g5_y_yz ), &
  zeta5_y_yz( num_g5_y_yz ), lambda5_y_yz( num_g5_y_yz ), &
  g5_y_yv( natom_y, num_g5_y_yv ), &
   eta5_y_yv( num_g5_y_yv ),  theta5_y_yv( num_g5_y_yv ), &
  zeta5_y_yv( num_g5_y_yv ), lambda5_y_yv( num_g5_y_yv ), &
  g5_y_zz( natom_y, num_g5_y_zz ), &
   eta5_y_zz( num_g5_y_zz ),  theta5_y_zz( num_g5_y_zz ), &
  zeta5_y_zz( num_g5_y_zz ), lambda5_y_zz( num_g5_y_zz ), &
  g5_y_zv( natom_y, num_g5_y_zv ), &
   eta5_y_zv( num_g5_y_zv ),  theta5_y_zv( num_g5_y_zv ), &
  zeta5_y_zv( num_g5_y_zv ), lambda5_y_zv( num_g5_y_zv ), &
  g5_y_vv( natom_y, num_g5_y_vv ), &
   eta5_y_vv( num_g5_y_vv ),  theta5_y_vv( num_g5_y_vv ), &
  zeta5_y_vv( num_g5_y_vv ), lambda5_y_vv( num_g5_y_vv )
double precision &
  g2_y_x_mpi( natom_y, num_g2_y_x ), &
  g2_y_y_mpi( natom_y, num_g2_y_x ), &
  g2_y_z_mpi( natom_y, num_g2_y_x ), &
  g2_y_v_mpi( natom_y, num_g2_y_x ), &
  g5_y_xx_mpi( natom_y, num_g5_y_xx ), &
  g5_y_xy_mpi( natom_y, num_g5_y_xy ), &
  g5_y_xz_mpi( natom_y, num_g5_y_xz ), &
  g5_y_xv_mpi( natom_y, num_g5_y_xv ), &
  g5_y_yy_mpi( natom_y, num_g5_y_yy ), &
  g5_y_yz_mpi( natom_y, num_g5_y_yz ), &
  g5_y_yv_mpi( natom_y, num_g5_y_yv ), &
  g5_y_zz_mpi( natom_y, num_g5_y_zz ), &
  g5_y_zv_mpi( natom_y, num_g5_y_zv ), &
  g5_y_vv_mpi( natom_y, num_g5_y_vv )
!Z
double precision &
  g2_z_x( natom_z, num_g2_z_x ), eta2_z_x( num_g2_z_x ), rs2_z_x( num_g2_z_x ), &
  g2_z_y( natom_z, num_g2_z_y ), eta2_z_y( num_g2_z_y ), rs2_z_y( num_g2_z_y ), &
  g2_z_z( natom_z, num_g2_z_z ), eta2_z_z( num_g2_z_z ), rs2_z_z( num_g2_z_z ), &
  g2_z_v( natom_z, num_g2_z_v ), eta2_z_v( num_g2_z_v ), rs2_z_v( num_g2_z_v )
double precision &
  g5_z_xx( natom_z, num_g5_z_xx ), &
   eta5_z_xx( num_g5_z_xx ),  theta5_z_xx( num_g5_z_xx ), &
  zeta5_z_xx( num_g5_z_xx ), lambda5_z_xx( num_g5_z_xx ), &
  g5_z_xy( natom_z, num_g5_z_xy ), &
   eta5_z_xy( num_g5_z_xy ),  theta5_z_xy( num_g5_z_xy ), &
  zeta5_z_xy( num_g5_z_xy ), lambda5_z_xy( num_g5_z_xy ), &
  g5_z_xz( natom_z, num_g5_z_xz ), &
   eta5_z_xz( num_g5_z_xz ),  theta5_z_xz( num_g5_z_xz ), &
  zeta5_z_xz( num_g5_z_xz ), lambda5_z_xz( num_g5_z_xz ), &
  g5_z_xv( natom_z, num_g5_z_xv ), &
   eta5_z_xv( num_g5_z_xv ),  theta5_z_xv( num_g5_z_xv ), &
  zeta5_z_xv( num_g5_z_xv ), lambda5_z_xv( num_g5_z_xv ), &
  g5_z_yy( natom_z, num_g5_z_yy ), &
   eta5_z_yy( num_g5_z_yy ),  theta5_z_yy( num_g5_z_yy ), &
  zeta5_z_yy( num_g5_z_yy ), lambda5_z_yy( num_g5_z_yy ), &
  g5_z_yz( natom_z, num_g5_z_yz ), &
   eta5_z_yz( num_g5_z_yz ),  theta5_z_yz( num_g5_z_yz ), &
  zeta5_z_yz( num_g5_z_yz ), lambda5_z_yz( num_g5_z_yz ), &
  g5_z_yv( natom_z, num_g5_z_yv ), &
   eta5_z_yv( num_g5_z_yv ),  theta5_z_yv( num_g5_z_yv ), &
  zeta5_z_yv( num_g5_z_yv ), lambda5_z_yv( num_g5_z_yv ), &
  g5_z_zz( natom_z, num_g5_z_zz ), &
   eta5_z_zz( num_g5_z_zz ),  theta5_z_zz( num_g5_z_zz ), &
  zeta5_z_zz( num_g5_z_zz ), lambda5_z_zz( num_g5_z_zz ), &
  g5_z_zv( natom_z, num_g5_z_zv ), &
   eta5_z_zv( num_g5_z_zv ),  theta5_z_zv( num_g5_z_zv ), &
  zeta5_z_zv( num_g5_z_zv ), lambda5_z_zv( num_g5_z_zv ), &
  g5_z_vv( natom_z, num_g5_z_vv ), &
   eta5_z_vv( num_g5_z_vv ),  theta5_z_vv( num_g5_z_vv ), &
  zeta5_z_vv( num_g5_z_vv ), lambda5_z_vv( num_g5_z_vv )
double precision &
  g2_z_x_mpi( natom_z, num_g2_z_x ), &
  g2_z_y_mpi( natom_z, num_g2_z_x ), &
  g2_z_z_mpi( natom_z, num_g2_z_x ), &
  g2_z_v_mpi( natom_z, num_g2_z_x ), &
  g5_z_xx_mpi( natom_z, num_g5_z_xx ), &
  g5_z_xy_mpi( natom_z, num_g5_z_xy ), &
  g5_z_xz_mpi( natom_z, num_g5_z_xz ), &
  g5_z_xv_mpi( natom_z, num_g5_z_xv ), &
  g5_z_yy_mpi( natom_z, num_g5_z_yy ), &
  g5_z_yz_mpi( natom_z, num_g5_z_yz ), &
  g5_z_yv_mpi( natom_z, num_g5_z_yv ), &
  g5_z_zz_mpi( natom_z, num_g5_z_zz ), &
  g5_z_zv_mpi( natom_z, num_g5_z_zv ), &
  g5_z_vv_mpi( natom_z, num_g5_z_vv )
!V
double precision &
  g2_v_x( natom_v, num_g2_v_x ), eta2_v_x( num_g2_v_x ), rs2_v_x( num_g2_v_x ), &
  g2_v_y( natom_v, num_g2_v_y ), eta2_v_y( num_g2_v_y ), rs2_v_y( num_g2_v_y ), &
  g2_v_z( natom_v, num_g2_v_z ), eta2_v_z( num_g2_v_z ), rs2_v_z( num_g2_v_z ), &
  g2_v_v( natom_v, num_g2_v_v ), eta2_v_v( num_g2_v_v ), rs2_v_v( num_g2_v_v )
double precision &
  g5_v_xx( natom_v, num_g5_v_xx ), &
   eta5_v_xx( num_g5_v_xx ),  theta5_v_xx( num_g5_v_xx ), &
  zeta5_v_xx( num_g5_v_xx ), lambda5_v_xx( num_g5_v_xx ), &
  g5_v_xy( natom_v, num_g5_v_xy ), &
   eta5_v_xy( num_g5_v_xy ),  theta5_v_xy( num_g5_v_xy ), &
  zeta5_v_xy( num_g5_v_xy ), lambda5_v_xy( num_g5_v_xy ), &
  g5_v_xz( natom_v, num_g5_v_xz ), &
   eta5_v_xz( num_g5_v_xz ),  theta5_v_xz( num_g5_v_xz ), &
  zeta5_v_xz( num_g5_v_xz ), lambda5_v_xz( num_g5_v_xz ), &
  g5_v_xv( natom_v, num_g5_v_xv ), &
   eta5_v_xv( num_g5_v_xv ),  theta5_v_xv( num_g5_v_xv ), &
  zeta5_v_xv( num_g5_v_xv ), lambda5_v_xv( num_g5_v_xv ), &
  g5_v_yy( natom_v, num_g5_v_yy ), &
   eta5_v_yy( num_g5_v_yy ),  theta5_v_yy( num_g5_v_yy ), &
  zeta5_v_yy( num_g5_v_yy ), lambda5_v_yy( num_g5_v_yy ), &
  g5_v_yz( natom_v, num_g5_v_yz ), &
   eta5_v_yz( num_g5_v_yz ),  theta5_v_yz( num_g5_v_yz ), &
  zeta5_v_yz( num_g5_v_yz ), lambda5_v_yz( num_g5_v_yz ), &
  g5_v_yv( natom_v, num_g5_v_yv ), &
   eta5_v_yv( num_g5_v_yv ),  theta5_v_yv( num_g5_v_yv ), &
  zeta5_v_yv( num_g5_v_yv ), lambda5_v_yv( num_g5_v_yv ), &
  g5_v_zz( natom_v, num_g5_v_zz ), &
   eta5_v_zz( num_g5_v_zz ),  theta5_v_zz( num_g5_v_zz ), &
  zeta5_v_zz( num_g5_v_zz ), lambda5_v_zz( num_g5_v_zz ), &
  g5_v_zv( natom_v, num_g5_v_zv ), &
   eta5_v_zv( num_g5_v_zv ),  theta5_v_zv( num_g5_v_zv ), &
  zeta5_v_zv( num_g5_v_zv ), lambda5_v_zv( num_g5_v_zv ), &
  g5_v_vv( natom_v, num_g5_v_vv ), &
   eta5_v_vv( num_g5_v_vv ),  theta5_v_vv( num_g5_v_vv ), &
  zeta5_v_vv( num_g5_v_vv ), lambda5_v_vv( num_g5_v_vv )
double precision &
  g2_v_x_mpi( natom_v, num_g2_v_x ), &
  g2_v_y_mpi( natom_v, num_g2_v_x ), &
  g2_v_z_mpi( natom_v, num_g2_v_x ), &
  g2_v_v_mpi( natom_v, num_g2_v_x ), &
  g5_v_xx_mpi( natom_v, num_g5_v_xx ), &
  g5_v_xy_mpi( natom_v, num_g5_v_xy ), &
  g5_v_xz_mpi( natom_v, num_g5_v_xz ), &
  g5_v_xv_mpi( natom_v, num_g5_v_xv ), &
  g5_v_yy_mpi( natom_v, num_g5_v_yy ), &
  g5_v_yz_mpi( natom_v, num_g5_v_yz ), &
  g5_v_yv_mpi( natom_v, num_g5_v_yv ), &
  g5_v_zz_mpi( natom_v, num_g5_v_zz ), &
  g5_v_zv_mpi( natom_v, num_g5_v_zv ), &
  g5_v_vv_mpi( natom_v, num_g5_v_vv )

!X
double precision &
 g2_deriv_x_x(  natom_x, natom, 3, num_g2_x_x  ), &
 g2_deriv_x_y(  natom_x, natom, 3, num_g2_x_y  ), &
 g2_deriv_x_z(  natom_x, natom, 3, num_g2_x_z  ), &
 g2_deriv_x_v(  natom_x, natom, 3, num_g2_x_v  ), &
 g5_deriv_x_xx( natom_x, natom, 3, num_g5_x_xx ), &
 g5_deriv_x_xy( natom_x, natom, 3, num_g5_x_xy ), &
 g5_deriv_x_xz( natom_x, natom, 3, num_g5_x_xz ), &
 g5_deriv_x_xv( natom_x, natom, 3, num_g5_x_xv ), &
 g5_deriv_x_yy( natom_x, natom, 3, num_g5_x_yy ), &
 g5_deriv_x_yz( natom_x, natom, 3, num_g5_x_yz ), &
 g5_deriv_x_yv( natom_x, natom, 3, num_g5_x_yv ), &
 g5_deriv_x_zz( natom_x, natom, 3, num_g5_x_zz ), &
 g5_deriv_x_zv( natom_x, natom, 3, num_g5_x_zv ), &
 g5_deriv_x_vv( natom_x, natom, 3, num_g5_x_vv )
!Y
double precision &
 g2_deriv_y_x(  natom_y, natom, 3, num_g2_y_x  ), &
 g2_deriv_y_y(  natom_y, natom, 3, num_g2_y_y  ), &
 g2_deriv_y_z(  natom_y, natom, 3, num_g2_y_z  ), &
 g2_deriv_y_v(  natom_y, natom, 3, num_g2_y_v  ), &
 g5_deriv_y_xx( natom_y, natom, 3, num_g5_y_xx ), &
 g5_deriv_y_xy( natom_y, natom, 3, num_g5_y_xy ), &
 g5_deriv_y_xz( natom_y, natom, 3, num_g5_y_xz ), &
 g5_deriv_y_xv( natom_y, natom, 3, num_g5_y_xv ), &
 g5_deriv_y_yy( natom_y, natom, 3, num_g5_y_yy ), &
 g5_deriv_y_yz( natom_y, natom, 3, num_g5_y_yz ), &
 g5_deriv_y_yv( natom_y, natom, 3, num_g5_y_yv ), &
 g5_deriv_y_zz( natom_y, natom, 3, num_g5_y_zz ), &
 g5_deriv_y_zv( natom_y, natom, 3, num_g5_y_zv ), &
 g5_deriv_y_vv( natom_y, natom, 3, num_g5_y_vv )
!Z
double precision &
 g2_deriv_z_x(  natom_z, natom, 3, num_g2_z_x  ), &
 g2_deriv_z_y(  natom_z, natom, 3, num_g2_z_y  ), &
 g2_deriv_z_z(  natom_z, natom, 3, num_g2_z_z  ), &
 g2_deriv_z_v(  natom_z, natom, 3, num_g2_z_v  ), &
 g5_deriv_z_xx( natom_z, natom, 3, num_g5_z_xx ), &
 g5_deriv_z_xy( natom_z, natom, 3, num_g5_z_xy ), &
 g5_deriv_z_xz( natom_z, natom, 3, num_g5_z_xz ), &
 g5_deriv_z_xv( natom_z, natom, 3, num_g5_z_xv ), &
 g5_deriv_z_yy( natom_z, natom, 3, num_g5_z_yy ), &
 g5_deriv_z_yz( natom_z, natom, 3, num_g5_z_yz ), &
 g5_deriv_z_yv( natom_z, natom, 3, num_g5_z_yv ), &
 g5_deriv_z_zz( natom_z, natom, 3, num_g5_z_zz ), &
 g5_deriv_z_zv( natom_z, natom, 3, num_g5_z_zv ), &
 g5_deriv_z_vv( natom_z, natom, 3, num_g5_z_vv )
!V
double precision &
 g2_deriv_v_x(  natom_v, natom, 3, num_g2_v_x  ), &
 g2_deriv_v_y(  natom_v, natom, 3, num_g2_v_y  ), &
 g2_deriv_v_z(  natom_v, natom, 3, num_g2_v_z  ), &
 g2_deriv_v_v(  natom_v, natom, 3, num_g2_v_v  ), &
 g5_deriv_v_xx( natom_v, natom, 3, num_g5_v_xx ), &
 g5_deriv_v_xy( natom_v, natom, 3, num_g5_v_xy ), &
 g5_deriv_v_xz( natom_v, natom, 3, num_g5_v_xz ), &
 g5_deriv_v_xv( natom_v, natom, 3, num_g5_v_xv ), &
 g5_deriv_v_yy( natom_v, natom, 3, num_g5_v_yy ), &
 g5_deriv_v_yz( natom_v, natom, 3, num_g5_v_yz ), &
 g5_deriv_v_yv( natom_v, natom, 3, num_g5_v_yv ), &
 g5_deriv_v_zz( natom_v, natom, 3, num_g5_v_zz ), &
 g5_deriv_v_zv( natom_v, natom, 3, num_g5_v_zv ), &
 g5_deriv_v_vv( natom_v, natom, 3, num_g5_v_vv )

double precision &
 g2_deriv_x_x_mpi(  natom_x, natom, 3, num_g2_x_x  ), &
 g2_deriv_x_y_mpi(  natom_x, natom, 3, num_g2_x_y  ), &
 g2_deriv_x_z_mpi(  natom_x, natom, 3, num_g2_x_z  ), &
 g2_deriv_x_v_mpi(  natom_x, natom, 3, num_g2_x_v  ), &
 g5_deriv_x_xx_mpi( natom_x, natom, 3, num_g5_x_xx ), &
 g5_deriv_x_xy_mpi( natom_x, natom, 3, num_g5_x_xy ), &
 g5_deriv_x_xz_mpi( natom_x, natom, 3, num_g5_x_xz ), &
 g5_deriv_x_xv_mpi( natom_x, natom, 3, num_g5_x_xv ), &
 g5_deriv_x_yy_mpi( natom_x, natom, 3, num_g5_x_yy ), &
 g5_deriv_x_yz_mpi( natom_x, natom, 3, num_g5_x_yz ), &
 g5_deriv_x_yv_mpi( natom_x, natom, 3, num_g5_x_yv ), &
 g5_deriv_x_zz_mpi( natom_x, natom, 3, num_g5_x_zz ), &
 g5_deriv_x_zv_mpi( natom_x, natom, 3, num_g5_x_zv ), &
 g5_deriv_x_vv_mpi( natom_x, natom, 3, num_g5_x_vv )
!Y
double precision &
 g2_deriv_y_x_mpi(  natom_y, natom, 3, num_g2_y_x  ), &
 g2_deriv_y_y_mpi(  natom_y, natom, 3, num_g2_y_y  ), &
 g2_deriv_y_z_mpi(  natom_y, natom, 3, num_g2_y_z  ), &
 g2_deriv_y_v_mpi(  natom_y, natom, 3, num_g2_y_v  ), &
 g5_deriv_y_xx_mpi( natom_y, natom, 3, num_g5_y_xx ), &
 g5_deriv_y_xy_mpi( natom_y, natom, 3, num_g5_y_xy ), &
 g5_deriv_y_xz_mpi( natom_y, natom, 3, num_g5_y_xz ), &
 g5_deriv_y_xv_mpi( natom_y, natom, 3, num_g5_y_xv ), &
 g5_deriv_y_yy_mpi( natom_y, natom, 3, num_g5_y_yy ), &
 g5_deriv_y_yz_mpi( natom_y, natom, 3, num_g5_y_yz ), &
 g5_deriv_y_yv_mpi( natom_y, natom, 3, num_g5_y_yv ), &
 g5_deriv_y_zz_mpi( natom_y, natom, 3, num_g5_y_zz ), &
 g5_deriv_y_zv_mpi( natom_y, natom, 3, num_g5_y_zv ), &
 g5_deriv_y_vv_mpi( natom_y, natom, 3, num_g5_y_vv )
!Z
double precision &
 g2_deriv_z_x_mpi(  natom_z, natom, 3, num_g2_z_x  ), &
 g2_deriv_z_y_mpi(  natom_z, natom, 3, num_g2_z_y  ), &
 g2_deriv_z_z_mpi(  natom_z, natom, 3, num_g2_z_z  ), &
 g2_deriv_z_v_mpi(  natom_z, natom, 3, num_g2_z_v  ), &
 g5_deriv_z_xx_mpi( natom_z, natom, 3, num_g5_z_xx ), &
 g5_deriv_z_xy_mpi( natom_z, natom, 3, num_g5_z_xy ), &
 g5_deriv_z_xz_mpi( natom_z, natom, 3, num_g5_z_xz ), &
 g5_deriv_z_xv_mpi( natom_z, natom, 3, num_g5_z_xv ), &
 g5_deriv_z_yy_mpi( natom_z, natom, 3, num_g5_z_yy ), &
 g5_deriv_z_yz_mpi( natom_z, natom, 3, num_g5_z_yz ), &
 g5_deriv_z_yv_mpi( natom_z, natom, 3, num_g5_z_yv ), &
 g5_deriv_z_zz_mpi( natom_z, natom, 3, num_g5_z_zz ), &
 g5_deriv_z_zv_mpi( natom_z, natom, 3, num_g5_z_zv ), &
 g5_deriv_z_vv_mpi( natom_z, natom, 3, num_g5_z_vv )
!V
double precision &
 g2_deriv_v_x_mpi(  natom_v, natom, 3, num_g2_v_x  ), &
 g2_deriv_v_y_mpi(  natom_v, natom, 3, num_g2_v_y  ), &
 g2_deriv_v_z_mpi(  natom_v, natom, 3, num_g2_v_z  ), &
 g2_deriv_v_v_mpi(  natom_v, natom, 3, num_g2_v_v  ), &
 g5_deriv_v_xx_mpi( natom_v, natom, 3, num_g5_v_xx ), &
 g5_deriv_v_xy_mpi( natom_v, natom, 3, num_g5_v_xy ), &
 g5_deriv_v_xz_mpi( natom_v, natom, 3, num_g5_v_xz ), &
 g5_deriv_v_xv_mpi( natom_v, natom, 3, num_g5_v_xv ), &
 g5_deriv_v_yy_mpi( natom_v, natom, 3, num_g5_v_yy ), &
 g5_deriv_v_yz_mpi( natom_v, natom, 3, num_g5_v_yz ), &
 g5_deriv_v_yv_mpi( natom_v, natom, 3, num_g5_v_yv ), &
 g5_deriv_v_zz_mpi( natom_v, natom, 3, num_g5_v_zz ), &
 g5_deriv_v_zv_mpi( natom_v, natom, 3, num_g5_v_zv ), &
 g5_deriv_v_vv_mpi( natom_v, natom, 3, num_g5_v_vv )

double precision Rij, Rik, g5_const_3
double precision rc
integer io_sf


 do i = 1, natom
   posi1(i,1) = posi(i,1) * xlattice
   posi1(i,2) = posi(i,2) * ylattice
   posi1(i,3) = posi(i,3) * zlattice
 enddo


!Repeating box
 if( min( xlattice, ylattice, zlattice ) >= rc )then
   lattice_box = 27
   boxsize     = 3
   boxcenter   = 2
   center      = natom * 13
 elseif( rc > min( xlattice, ylattice, zlattice ) &
         .and. 2.0d0 * min( xlattice, ylattice, zlattice ) >= rc )then
   lattice_box = 125
   boxsize     = 5
   boxcenter   = 3
   center      = natom * 62
 else
   write(*,*) "Error in sf_4-1 : revise program !!"
   stop
 endif

 if( myrank == 0 .and. io_sf == 0 )then 
   write(*,*) " "
   write(*,*) " Repeating box size: ", lattice_box
   write(*,*) " "
   io_sf = 1
 endif


 allocate( posi2( lattice_box * natom, 3 ) )

 m = 0
 do ix1 = 1, boxsize
   do iy1 = 1, boxsize
     do iz1 = 1, boxsize
       do i = 1, natom
         m = m + 1
         posi2(m,1) = posi1(i,1) + xlattice * ( ix1 - boxcenter )
         posi2(m,2) = posi1(i,2) + ylattice * ( iy1 - boxcenter )
         posi2(m,3) = posi1(i,3) + zlattice * ( iz1 - boxcenter )
       enddo
     enddo
   enddo
 enddo


!Neighboring list
!if( mod(time,neighbor_check)==1 .or. jamp == 1 )then
     list = 0
  do i = 1, natom
     m = 0
     n = 0
    do ix1 = 1, boxsize
      do iy1 = 1, boxsize
        do iz1 = 1, boxsize
          do j = 1, natom

            m = m + 1

            Rij = dsqrt(   ( posi1(i,1) - posi2(m,1) )**2 &
                         + ( posi1(i,2) - posi2(m,2) )**2 &
                         + ( posi1(i,3) - posi2(m,3) )**2 )

            if( Rij <= rc )then !+ Rn )then

              n = n + 1

              list(i,n) = m

              if( n == neighbor )then
                write(*,*) "Insufficient neighbor list !!"
                stop
              endif

            endif
          enddo
        enddo
      enddo
    enddo
  enddo
!endif


!  g2_x & g5_x
   g2_x_x  = 0.0d0 ; g2_x_y  = 0.0d0 ; g2_x_z  = 0.0d0 ; g2_x_v  = 0.0d0
   g5_x_xx = 0.0d0 ; g5_x_xy = 0.0d0 ; g5_x_xz = 0.0d0 ; g5_x_xv = 0.0d0
   g5_x_yy = 0.0d0 ; g5_x_yz = 0.0d0 ; g5_x_yv = 0.0d0
   g5_x_zz = 0.0d0 ; g5_x_zv = 0.0d0
   g5_x_vv = 0.0d0
!  g2_y & g5_y
   g2_y_x  = 0.0d0 ; g2_y_y  = 0.0d0 ; g2_y_z  = 0.0d0 ; g2_y_v  = 0.0d0
   g5_y_xx = 0.0d0 ; g5_y_xy = 0.0d0 ; g5_y_xz = 0.0d0 ; g5_y_xv = 0.0d0
   g5_y_yy = 0.0d0 ; g5_y_yz = 0.0d0 ; g5_y_yv = 0.0d0
   g5_y_zz = 0.0d0 ; g5_y_zv = 0.0d0
   g5_y_vv = 0.0d0
!  g2_z & g5_z
   g2_z_x  = 0.0d0 ; g2_z_y  = 0.0d0 ; g2_z_z  = 0.0d0 ; g2_z_v  = 0.0d0
   g5_z_xx = 0.0d0 ; g5_z_xy = 0.0d0 ; g5_z_xz = 0.0d0 ; g5_z_xv = 0.0d0
   g5_z_yy = 0.0d0 ; g5_z_yz = 0.0d0 ; g5_z_yv = 0.0d0
   g5_z_zz = 0.0d0 ; g5_z_zv = 0.0d0
   g5_z_vv = 0.0d0
!  g2_v & g5_v
   g2_v_x  = 0.0d0 ; g2_v_y  = 0.0d0 ; g2_v_z  = 0.0d0 ; g2_v_v  = 0.0d0
   g5_v_xx = 0.0d0 ; g5_v_xy = 0.0d0 ; g5_v_xz = 0.0d0 ; g5_v_xv = 0.0d0
   g5_v_yy = 0.0d0 ; g5_v_yz = 0.0d0 ; g5_v_yv = 0.0d0
   g5_v_zz = 0.0d0 ; g5_v_zv = 0.0d0
   g5_v_vv = 0.0d0

!  g2_deriv_x & g5_deriv_x
   g2_deriv_x_x  = 0.0d0 ; g2_deriv_x_y  = 0.0d0 ; g2_deriv_x_z  = 0.0d0 ; g2_deriv_x_v  = 0.0d0 
   g5_deriv_x_xx = 0.0d0 ; g5_deriv_x_xy = 0.0d0 ; g5_deriv_x_xz = 0.0d0 ; g5_deriv_x_xv = 0.0d0
   g5_deriv_x_yy = 0.0d0 ; g5_deriv_x_yz = 0.0d0 ; g5_deriv_x_yv = 0.0d0
   g5_deriv_x_zz = 0.0d0 ; g5_deriv_x_zv = 0.0d0
   g5_deriv_x_vv = 0.0d0
!  g2_deriv_y & g5_deriv_y
   g2_deriv_y_x  = 0.0d0 ; g2_deriv_y_y  = 0.0d0 ; g2_deriv_y_z  = 0.0d0 ; g2_deriv_y_v  = 0.0d0
   g5_deriv_y_xx = 0.0d0 ; g5_deriv_y_xy = 0.0d0 ; g5_deriv_y_xz = 0.0d0 ; g5_deriv_y_xv = 0.0d0
   g5_deriv_y_yy = 0.0d0 ; g5_deriv_y_yz = 0.0d0 ; g5_deriv_y_yv = 0.0d0
   g5_deriv_y_zz = 0.0d0 ; g5_deriv_y_zv = 0.0d0
   g5_deriv_y_vv = 0.0d0
!  g2_deriv_z & g5_deriv_z
   g2_deriv_z_x  = 0.0d0 ; g2_deriv_z_y  = 0.0d0 ; g2_deriv_z_z  = 0.0d0 ; g2_deriv_z_v  = 0.0d0
   g5_deriv_z_xx = 0.0d0 ; g5_deriv_z_xy = 0.0d0 ; g5_deriv_z_xz = 0.0d0 ; g5_deriv_z_xv = 0.0d0
   g5_deriv_z_yy = 0.0d0 ; g5_deriv_z_yz = 0.0d0 ; g5_deriv_z_yv = 0.0d0
   g5_deriv_z_zz = 0.0d0 ; g5_deriv_z_zv = 0.0d0
   g5_deriv_z_vv = 0.0d0
!  g2_deriv_v & g5_deriv_v
   g2_deriv_v_x  = 0.0d0 ; g2_deriv_v_y  = 0.0d0 ; g2_deriv_v_z  = 0.0d0 ; g2_deriv_v_v  = 0.0d0
   g5_deriv_v_xx = 0.0d0 ; g5_deriv_v_xy = 0.0d0 ; g5_deriv_v_xz = 0.0d0 ; g5_deriv_v_xv = 0.0d0
   g5_deriv_v_yy = 0.0d0 ; g5_deriv_v_yz = 0.0d0 ; g5_deriv_v_yv = 0.0d0
   g5_deriv_v_zz = 0.0d0 ; g5_deriv_v_zv = 0.0d0
   g5_deriv_v_vv = 0.0d0

! i : taget atom
! j : 1st neighbor atom
! k : 2nd neighbot atom
 do i0 = 1, natom, num_mpi
   i = i0 + myrank
   if( i <= natom )then
 !do i = 1, natom ! Target atom

   ! G2_X_... & G5_X_...
   if( i <= natom_x )then

     do ix1 = 1, neighbor

       atom1 = list(i,ix1)

       if( atom1 == 0 )cycle
       if( atom1 - center == i )cycle

       j = atom1 - natom * ( ( atom1 - 1 ) / natom )
       Rij = dsqrt(   ( posi1(i,1) - posi2(atom1,1) )**2 &
                    + ( posi1(i,2) - posi2(atom1,2) )**2 &
                    + ( posi1(i,3) - posi2(atom1,3) )**2 )

       if( Rij < rc )then

         ! G2_X_X
         if( j <= natom_x )then
           call sf_4_g2( g2_x_x, g2_deriv_x_x, eta2_x_x, rs2_x_x, num_g2_x_x, &
                         posi1(i,1),     posi1(i,2),     posi1(i,3), &
                         posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                         natom, natom_x, natom_y, natom_z, natom_v, &
                         Rij, i, j, rc, natom_x )
         ! G2_X_Y
         elseif( natom_x < j .and. j <= natom_x + natom_y )then
           call sf_4_g2( g2_x_y, g2_deriv_x_y, eta2_x_y, rs2_x_y, num_g2_x_y, &
                         posi1(i,1),     posi1(i,2),     posi1(i,3), &
                         posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                         natom, natom_x, natom_y, natom_z, natom_v, &
                         Rij, i, j, rc, natom_x )
         ! G2_X_Z
         elseif( natom_x + natom_y < j .and. j <= natom_x + natom_y + natom_z )then
           call sf_4_g2( g2_x_z, g2_deriv_x_z, eta2_x_z, rs2_x_z, num_g2_x_z, &
                         posi1(i,1),     posi1(i,2),     posi1(i,3), &
                         posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                         natom, natom_x, natom_y, natom_z, natom_v, &
                         Rij, i, j, rc, natom_x )
         ! G2_X_V
         elseif( natom_x + natom_y + natom_z < j &
                 .and. j <= natom_x + natom_y + natom_z + natom_v )then
           call sf_4_g2( g2_x_v, g2_deriv_x_v, eta2_x_v, rs2_x_v, num_g2_x_v, &
                         posi1(i,1),     posi1(i,2),     posi1(i,3), &
                         posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                         natom, natom_x, natom_y, natom_z, natom_v, &
                         Rij, i, j, rc, natom_x )
         else
           write(*,*) "Error(sf_4): G2_X_*"
           stop
         endif ! j

         do ix2 = 1, neighbor

           atom2 = list(i,ix2)

           if( atom2 == 0 )cycle
           if( atom2 - center == i )cycle
           if( atom1 == atom2 )cycle

           k = atom2 - natom * ( ( atom2 - 1 ) / natom )
           Rik = dsqrt(   ( posi1(i,1) - posi2(atom2,1) )**2 &
                        + ( posi1(i,2) - posi2(atom2,2) )**2 &
                        + ( posi1(i,3) - posi2(atom2,3) )**2 )

           if( Rik <= rc )then
             
             ! G5_X_XX
             if( j <= natom_x .and. &
                 k <= natom_x )then
               call sf_4_g5( g5_x_xx, g5_deriv_x_xx, &
                             eta5_x_xx, lambda5_x_xx, zeta5_x_xx, num_g5_x_xx, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_x )
             ! G5_X_XY
             elseif( j <= natom_x .and. &
                     natom_x < k .and. k <= natom_x + natom_y )then
               call sf_4_g5( g5_x_xy, g5_deriv_x_xy, &
                             eta5_x_xy, lambda5_x_xy, zeta5_x_xy, num_g5_x_xy, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_x )
             ! G5_X_XZ
             elseif( j <= natom_x .and. &
                     natom_x + natom_y < k .and. k <= natom_x + natom_y + natom_z )then
               call sf_4_g5( g5_x_xz, g5_deriv_x_xz, &
                             eta5_x_xz, lambda5_x_xz, zeta5_x_xz, num_g5_x_xz, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_x )
             ! G5_X_XV
             elseif( j <= natom_x .and. &
                     natom_x + natom_y + natom_z < k &
                     .and. k <= natom_x + natom_y + natom_z + natom_v )then
               call sf_4_g5( g5_x_xv, g5_deriv_x_xv, &
                             eta5_x_xv, lambda5_x_xv, zeta5_x_xv, num_g5_x_xv, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_x )
             ! G5_X_YY
             elseif( natom_x < j .and. j <= natom_x + natom_y .and. &
                     natom_x < k .and. k <= natom_x + natom_y )then
               call sf_4_g5( g5_x_yy, g5_deriv_x_yy, &
                             eta5_x_yy, lambda5_x_yy, zeta5_x_yy, num_g5_x_yy, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_x )
             ! G5_X_YZ
             elseif( natom_x < j .and. j <= natom_x + natom_y .and. &
                     natom_x + natom_y < k .and. k <= natom_x + natom_y + natom_z )then
               call sf_4_g5( g5_x_yz, g5_deriv_x_yz, &
                             eta5_x_yz, lambda5_x_yz, zeta5_x_yz, num_g5_x_yz, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_x )
             ! G5_X_YV
             elseif( natom_x < j .and. j <= natom_x + natom_y .and. &
                     natom_x + natom_y + natom_z < k &
                     .and. k <= natom_x + natom_y + natom_z + natom_v )then
               call sf_4_g5( g5_x_yv, g5_deriv_x_yv, &
                             eta5_x_yv, lambda5_x_yv, zeta5_x_yv, num_g5_x_yv, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_x )
             ! G5_X_ZZ
             elseif( natom_x + natom_y < j .and. j <= natom_x + natom_y + natom_z .and. &
                     natom_x + natom_y < k .and. k <= natom_x + natom_y + natom_z )then
               call sf_4_g5( g5_x_zz, g5_deriv_x_zz, &
                             eta5_x_zz, lambda5_x_zz, zeta5_x_zz, num_g5_x_zz, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_x )
             ! G5_X_ZV
             elseif( natom_x + natom_y < j .and. j <= natom_x + natom_y + natom_z .and. &
                     natom_x + natom_y + natom_z < k &
                     .and. k <= natom_x + natom_y + natom_z + natom_v )then
               call sf_4_g5( g5_x_zv, g5_deriv_x_zv, &
                             eta5_x_zv, lambda5_x_zv, zeta5_x_zv, num_g5_x_zv, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_x )
             ! G5_X_VV
             elseif( natom_x + natom_y + natom_z < j &
                     .and. j <= natom_x + natom_y + natom_z + natom_v .and. &
                     natom_x + natom_y + natom_z < k &
                     .and. k <= natom_x + natom_y + natom_z + natom_v )then
               call sf_4_g5( g5_x_vv, g5_deriv_x_vv, &
                             eta5_x_vv, lambda5_x_vv, zeta5_x_vv, num_g5_x_vv, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_x )
             ! G5_X_YX
             elseif( natom_x < j .and. j <= natom_x + natom_y .and. &
                     k <= natom_x )then
               cycle
             ! G5_X_ZX
             elseif( natom_x + natom_y < j .and. j <= natom_x + natom_y + natom_z .and. &
                     k <= natom_x )then
               cycle
             ! G5_X_ZY
             elseif( natom_x + natom_y < j .and. j <= natom_x + natom_y + natom_z .and. &
                     natom_x < k .and. k <= natom_x + natom_y )then
               cycle
             ! G5_X_VX
             elseif( natom_x + natom_y + natom_z < j &
                     .and. j <= natom_x + natom_y + natom_z + natom_v .and. &
                     k <= natom_x )then
               cycle
             ! G5_X_VY
             elseif( natom_x + natom_y + natom_z < j &
                     .and. j <= natom_x + natom_y + natom_z + natom_v .and. &
                     natom_x < k .and. k <= natom_x + natom_y )then
               cycle
             ! G5_X_VZ
             elseif( natom_x + natom_y + natom_z < j &
                     .and. j <= natom_x + natom_y + natom_z + natom_v .and. &
                     natom_x + natom_y < k .and. k <= natom_x + natom_y + natom_z )then
               cycle
             else
               write(*,*) "Error(sf_4): G5_X_*"
               stop
             endif ! j,k

           endif ! if( Rik <= rc )
         enddo ! do ix2
       endif ! if( Rij < rc )
     enddo ! do ix1 (g2_x... & g5_x...)

   ! G2_Y_... & G5_Y_...
   elseif( natom_x < i .and. i <= natom_x + natom_y )then

     do ix1 = 1, neighbor

       atom1 = list(i,ix1)

       if( atom1 == 0 )cycle
       if( atom1 - center == i )cycle

       j = atom1 - natom * ( ( atom1 - 1 ) / natom )
       Rij = dsqrt(   ( posi1(i,1) - posi2(atom1,1) )**2 &
                    + ( posi1(i,2) - posi2(atom1,2) )**2 &
                    + ( posi1(i,3) - posi2(atom1,3) )**2 )

       if( Rij < rc )then

         ! G2_Y_X
         if( j <= natom_x )then
           call sf_4_g2( g2_y_x, g2_deriv_y_x, eta2_y_x, rs2_y_x, num_g2_y_x, &
                         posi1(i,1),     posi1(i,2),     posi1(i,3), &
                         posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                         natom, natom_x, natom_y, natom_z, natom_v, &
                         Rij, i, j, rc, natom_y )
         ! G2_Y_Y
         elseif( natom_x < j .and. j <= natom_x + natom_y )then
           call sf_4_g2( g2_y_y, g2_deriv_y_y, eta2_y_y, rs2_y_y, num_g2_y_y, &
                         posi1(i,1),     posi1(i,2),     posi1(i,3), &
                         posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                         natom, natom_x, natom_y, natom_z, natom_v, &
                         Rij, i, j, rc, natom_y )
         ! G2_Y_Z
         elseif( natom_x + natom_y < j .and. j <= natom_x + natom_y + natom_z )then
           call sf_4_g2( g2_y_z, g2_deriv_y_z, eta2_y_z, rs2_y_z, num_g2_y_z, &
                         posi1(i,1),     posi1(i,2),     posi1(i,3), &
                         posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                         natom, natom_x, natom_y, natom_z, natom_v, &
                         Rij, i, j, rc, natom_y )
         ! G2_Y_V
         elseif( natom_x + natom_y + natom_z < j &
                 .and. j <= natom_x + natom_y + natom_z + natom_v )then
           call sf_4_g2( g2_y_v, g2_deriv_y_v, eta2_y_v, rs2_y_v, num_g2_y_v, &
                         posi1(i,1),     posi1(i,2),     posi1(i,3), &
                         posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                         natom, natom_x, natom_y, natom_z, natom_v, &
                         Rij, i, j, rc, natom_y )
         else
           write(*,*) "Error(sf_4): G2_Y_*"
           stop
         endif

         do ix2 = 1, neighbor

           atom2 = list(i,ix2)

           if( atom2 == 0 )cycle
           if( atom2 - center == i )cycle
           if( atom1 == atom2 )cycle

           k = atom2 - natom * ( ( atom2 - 1 ) / natom )
           Rik = dsqrt(   ( posi1(i,1) - posi2(atom2,1) )**2 &
                        + ( posi1(i,2) - posi2(atom2,2) )**2 &
                        + ( posi1(i,3) - posi2(atom2,3) )**2 )

           if( Rik <= rc )then
             
             ! G5_Y_XX
             if( j <= natom_x .and. &
                 k <= natom_x )then
               call sf_4_g5( g5_y_xx, g5_deriv_y_xx, &
                             eta5_y_xx, lambda5_y_xx, zeta5_y_xx, num_g5_y_xx, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_y )
             ! G5_Y_XY
             elseif( j <= natom_x .and. &
                     natom_x < k .and. k <= natom_x + natom_y )then
               call sf_4_g5( g5_y_xy, g5_deriv_y_xy, &
                             eta5_y_xy, lambda5_y_xy, zeta5_y_xy, num_g5_y_xy, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_y )
             ! G5_Y_XZ
             elseif( j <= natom_x .and. &
                     natom_x + natom_y < k .and. k <= natom_x + natom_y + natom_z )then
               call sf_4_g5( g5_y_xz, g5_deriv_y_xz, &
                             eta5_y_xz, lambda5_y_xz, zeta5_y_xz, num_g5_y_xz, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_y )
             ! G5_Y_XV
             elseif( j <= natom_x .and. &
                     natom_x + natom_y + natom_z < k &
                     .and. k <= natom_x + natom_y + natom_z + natom_v )then
               call sf_4_g5( g5_y_xv, g5_deriv_y_xv, &
                             eta5_y_xv, lambda5_y_xv, zeta5_y_xv, num_g5_y_xv, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_y )
             ! G5_Y_YY
             elseif( natom_x < j .and. j <= natom_x + natom_y .and. &
                     natom_x < k .and. k <= natom_x + natom_y )then
               call sf_4_g5( g5_y_yy, g5_deriv_y_yy, &
                             eta5_y_yy, lambda5_y_yy, zeta5_y_yy, num_g5_y_yy, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_y )
             ! G5_Y_YZ
             elseif( natom_x < j .and. j <= natom_x + natom_y .and. &
                     natom_x + natom_y < k .and. k <= natom_x + natom_y + natom_z )then
               call sf_4_g5( g5_y_yz, g5_deriv_y_yz, &
                             eta5_y_yz, lambda5_y_yz, zeta5_y_yz, num_g5_y_yz, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_y )
             ! G5_Y_YV
             elseif( natom_x < j .and. j <= natom_x + natom_y .and. &
                     natom_x + natom_y + natom_z < k &
                     .and. k <= natom_x + natom_y + natom_z + natom_v )then
               call sf_4_g5( g5_y_yv, g5_deriv_y_yv, &
                             eta5_y_yv, lambda5_y_yv, zeta5_y_yv, num_g5_y_yv, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_y )
             ! G5_Y_ZZ
             elseif( natom_x + natom_y < j .and. j <= natom_x + natom_y + natom_z .and. &
                     natom_x + natom_y < k .and. k <= natom_x + natom_y + natom_z )then
               call sf_4_g5( g5_y_zz, g5_deriv_y_zz, &
                             eta5_y_zz, lambda5_y_zz, zeta5_y_zz, num_g5_y_zz, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_y )
             ! G5_Y_ZV
             elseif( natom_x + natom_y < j .and. j <= natom_x + natom_y + natom_z .and. &
                     natom_x + natom_y + natom_z < k &
                     .and. k <= natom_x + natom_y + natom_z + natom_v )then
               call sf_4_g5( g5_y_zv, g5_deriv_y_zv, &
                             eta5_y_zv, lambda5_y_zv, zeta5_y_zv, num_g5_y_zv, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v,  &
                             Rij, Rik, i, j, k, rc, natom_y )
             ! G5_Y_VV
             elseif( natom_x + natom_y + natom_z < j &
                     .and. j <= natom_x + natom_y + natom_z + natom_v .and. &
                     natom_x + natom_y + natom_z < k &
                     .and. k <= natom_x + natom_y + natom_z + natom_v )then
               call sf_4_g5( g5_y_vv, g5_deriv_y_vv, &
                             eta5_y_vv, lambda5_y_vv, zeta5_y_vv, num_g5_y_vv, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v,  &
                             Rij, Rik, i, j, k, rc, natom_y )
             ! G5_Y_YX
             elseif( natom_x < j .and. j <= natom_x + natom_y .and. &
                     k <= natom_x )then
               cycle
             ! G5_Y_ZX
             elseif( natom_x + natom_y < j .and. j <= natom_x + natom_y + natom_z .and. &
                     k <= natom_x )then
               cycle
             ! G5_Y_ZY
             elseif( natom_x + natom_y < j .and. j <= natom_x + natom_y + natom_z .and. &
                     natom_x < k .and. k <= natom_x + natom_y )then
               cycle
             ! G5_Y_VX
             elseif( natom_x + natom_y + natom_z < j &
                     .and. j <= natom_x + natom_y + natom_z + natom_v .and. &
                     k <= natom_x )then
               cycle
             ! G5_Y_VY
             elseif( natom_x + natom_y + natom_z < j &
                     .and. j <= natom_x + natom_y + natom_z + natom_v .and. &
                     natom_x < k .and. k <= natom_x + natom_y )then
               cycle
             ! G5_Y_VZ
             elseif( natom_x + natom_y + natom_z < j &
                     .and. j <= natom_x + natom_y + natom_z + natom_v .and. &
                     natom_x + natom_y < k .and. k <= natom_x + natom_y + natom_z )then
               cycle
             else
               write(*,*) "Error(sf_4): G5_Y_*"
               stop
             endif ! j,k

           endif ! if( Rik <= rc )
         enddo ! do ix2
       endif ! if( Rij < rc )
     enddo ! do ix1 (g2_y... & g5_y...)

   ! G2_Z_... & G5_Z_...
   elseif( natom_x + natom_y < i .and. i <= natom_x + natom_y + natom_z )then

     do ix1 = 1, neighbor

       atom1 = list(i,ix1)

       if( atom1 == 0 )cycle
       if( atom1 - center == i )cycle

       j = atom1 - natom * ( ( atom1 - 1 ) / natom )
       Rij = dsqrt(   ( posi1(i,1) - posi2(atom1,1) )**2 &
                    + ( posi1(i,2) - posi2(atom1,2) )**2 &
                    + ( posi1(i,3) - posi2(atom1,3) )**2 )

       if( Rij < rc )then

         ! G2_Z_X
         if( j <= natom_x )then
           call sf_4_g2( g2_z_x, g2_deriv_z_x, eta2_z_x, rs2_z_x, num_g2_z_x, &
                         posi1(i,1),     posi1(i,2),     posi1(i,3), &
                         posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                         natom, natom_x, natom_y, natom_z, natom_v, &
                         Rij, i, j, rc, natom_z )
         ! G2_Z_Y
         elseif( natom_x < j .and. j <= natom_x + natom_y )then
           call sf_4_g2( g2_z_y, g2_deriv_z_y, eta2_z_y, rs2_z_y, num_g2_z_y, &
                         posi1(i,1),     posi1(i,2),     posi1(i,3), &
                         posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                         natom, natom_x, natom_y, natom_z, natom_v, &
                         Rij, i, j, rc, natom_z )
         ! G2_Z_Z
         elseif( natom_x + natom_y < j .and. j <= natom_x + natom_y + natom_z )then
           call sf_4_g2( g2_z_z, g2_deriv_z_z, eta2_z_z, rs2_z_z, num_g2_z_z, &
                         posi1(i,1),     posi1(i,2),     posi1(i,3), &
                         posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                         natom, natom_x, natom_y, natom_z, natom_v, &
                         Rij, i, j, rc, natom_z )
         ! G2_Z_V
         elseif( natom_x + natom_y + natom_z < j &
                 .and. j <= natom_x + natom_y + natom_z + natom_v )then
           call sf_4_g2( g2_z_v, g2_deriv_z_v, eta2_z_v, rs2_z_v, num_g2_z_v, &
                         posi1(i,1),     posi1(i,2),     posi1(i,3), &
                         posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                         natom, natom_x, natom_y, natom_z, natom_v, &
                         Rij, i, j, rc, natom_z )
         else
           write(*,*) "Error(sf_4): G2_Z_*"
           stop
         endif

         do ix2 = 1, neighbor

           atom2 = list(i,ix2)

           if( atom2 == 0 )cycle
           if( atom2 - center == i )cycle
           if( atom1 == atom2 )cycle

           k = atom2 - natom * ( ( atom2 - 1 ) / natom )
           Rik = dsqrt(   ( posi1(i,1) - posi2(atom2,1) )**2 &
                        + ( posi1(i,2) - posi2(atom2,2) )**2 &
                        + ( posi1(i,3) - posi2(atom2,3) )**2 )

           if( Rik <= rc )then
             
             ! G5_Z_XX
             if( j <= natom_x .and. &
                 k <= natom_x )then
               call sf_4_g5( g5_z_xx, g5_deriv_z_xx, &
                             eta5_z_xx, lambda5_z_xx, zeta5_z_xx, num_g5_z_xx, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_z )
             ! G5_Z_XY
             elseif( j <= natom_x .and. &
                     natom_x < k .and. k <= natom_x + natom_y )then
               call sf_4_g5( g5_z_xy, g5_deriv_z_xy, &
                             eta5_z_xy, lambda5_z_xy, zeta5_z_xy, num_g5_z_xy, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_z )
             ! G5_Z_XZ
             elseif( j <= natom_x .and. &
                     natom_x + natom_y < k .and. k <= natom_x + natom_y + natom_z )then
               call sf_4_g5( g5_z_xz, g5_deriv_z_xz, &
                             eta5_z_xz, lambda5_z_xz, zeta5_z_xz, num_g5_z_xz, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_z )
             ! G5_Z_XV
             elseif( j <= natom_x .and. &
                     natom_x + natom_y + natom_z < k &
                     .and. k <= natom_x + natom_y + natom_z + natom_v )then
               call sf_4_g5( g5_z_xv, g5_deriv_z_xv, &
                             eta5_z_xv, lambda5_z_xv, zeta5_z_xv, num_g5_z_xv, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_z )
             ! G5_Z_YY
             elseif( natom_x < j .and. j <= natom_x + natom_y .and. &
                     natom_x < k .and. k <= natom_x + natom_y )then
               call sf_4_g5( g5_z_yy, g5_deriv_z_yy, &
                             eta5_z_yy, lambda5_z_yy, zeta5_z_yy, num_g5_z_yy, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_z )
             ! G5_Z_YZ
             elseif( natom_x < j .and. j <= natom_x + natom_y .and. &
                     natom_x + natom_y < k .and. k <= natom_x + natom_y + natom_z )then
               call sf_4_g5( g5_z_yz, g5_deriv_z_yz, &
                             eta5_z_yz, lambda5_z_yz, zeta5_z_yz, num_g5_z_yz, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_z )
             ! G5_Z_YV
             elseif( natom_x < j .and. j <= natom_x + natom_y .and. &
                     natom_x + natom_y + natom_z < k &
                     .and. k <= natom_x + natom_y + natom_z + natom_v )then
               call sf_4_g5( g5_z_yv, g5_deriv_z_yv, &
                             eta5_z_yv, lambda5_z_yv, zeta5_z_yv, num_g5_z_yv, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_z )
             ! G5_Z_ZZ
             elseif( natom_x + natom_y < j .and. j <= natom_x + natom_y + natom_z .and. &
                     natom_x + natom_y < k .and. k <= natom_x + natom_y + natom_z )then
               call sf_4_g5( g5_z_zz, g5_deriv_z_zz, &
                             eta5_z_zz, lambda5_z_zz, zeta5_z_zz, num_g5_z_zz, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_z )
             ! G5_Z_ZV
             elseif( natom_x + natom_y < j .and. j <= natom_x + natom_y + natom_z .and. &
                     natom_x + natom_y + natom_z < k &
                     .and. k <= natom_x + natom_y + natom_z + natom_v )then
               call sf_4_g5( g5_z_zv, g5_deriv_z_zv, &
                             eta5_z_zv, lambda5_z_zv, zeta5_z_zv, num_g5_z_zv, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_z )
             ! G5_Z_VV
             elseif( natom_x + natom_y + natom_z < j &
                     .and. j <= natom_x + natom_y + natom_z + natom_v .and. &
                     natom_x + natom_y + natom_z < k &
                     .and. k <= natom_x + natom_y + natom_z + natom_v )then
               call sf_4_g5( g5_z_vv, g5_deriv_z_vv, &
                             eta5_z_vv, lambda5_z_vv, zeta5_z_vv, num_g5_z_vv, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_z )
             ! G5_Z_YX
             elseif( natom_x < j .and. j <= natom_x + natom_y .and. &
                     k <= natom_x )then
               cycle
             ! G5_Z_ZX
             elseif( natom_x + natom_y < j .and. j <= natom_x + natom_y + natom_z .and. &
                     k <= natom_x )then
               cycle
             ! G5_Z_ZY
             elseif( natom_x + natom_y < j .and. j <= natom_x + natom_y + natom_z .and. &
                     natom_x < k .and. k <= natom_x + natom_y )then
               cycle
             ! G5_Z_VX
             elseif( natom_x + natom_y + natom_z < j &
                     .and. j <= natom_x + natom_y + natom_z + natom_v .and. &
                     k <= natom_x )then
               cycle
             ! G5_Z_VY
             elseif( natom_x + natom_y + natom_z < j &
                     .and. j <= natom_x + natom_y + natom_z + natom_v .and. &
                     natom_x < k .and. k <= natom_x + natom_y )then
               cycle
             ! G5_Z_VZ
             elseif( natom_x + natom_y + natom_z < j &
                     .and. j <= natom_x + natom_y + natom_z + natom_v .and. &
                     natom_x + natom_y < k .and. k <= natom_x + natom_y + natom_z )then
               cycle
             else
               write(*,*) "Error(sf_4): G5_Z_*"
               stop
             endif

           endif ! if( Rik <= rc )
         enddo ! do ix2
       endif ! if( Rij < rc )
     enddo ! do ix1 (g2_z... & g5_z...)

   ! G2_V_... & G5_V_...
   elseif( natom_x + natom_y + natom_z < i &
           .and. i <= natom_x + natom_y + natom_z + natom_v )then

     do ix1 = 1, neighbor

       atom1 = list(i,ix1)

       if( atom1 == 0 )cycle
       if( atom1 - center == i )cycle

       j = atom1 - natom * ( ( atom1 - 1) / natom )
       Rij = dsqrt(   ( posi1(i,1) - posi2(atom1,1) )**2 &
                    + ( posi1(i,2) - posi2(atom1,2) )**2 &
                    + ( posi1(i,3) - posi2(atom1,3) )**2 )

       if( Rij < rc )then

         ! G2_V_X
         if( j <= natom_x )then
           call sf_4_g2( g2_v_x, g2_deriv_v_x, eta2_v_x, rs2_v_x, num_g2_v_x, &
                         posi1(i,1),     posi1(i,2),     posi1(i,3), &
                         posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                         natom, natom_x, natom_y, natom_z, natom_v, &
                         Rij, i, j, rc, natom_v )
         ! G2_V_Y
         elseif( natom_x < j .and. j <= natom_x + natom_y )then
           call sf_4_g2( g2_v_y, g2_deriv_v_y, eta2_v_y, rs2_v_y, num_g2_v_y, &
                         posi1(i,1),     posi1(i,2),     posi1(i,3), &
                         posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                         natom, natom_x, natom_y, natom_z, natom_v, &
                         Rij, i, j, rc, natom_v )
         ! G2_V_Z
         elseif( natom_x + natom_y < j .and. j <= natom_x + natom_y + natom_z )then
           call sf_4_g2( g2_v_z, g2_deriv_v_z, eta2_v_z, rs2_v_z, num_g2_v_z, &
                         posi1(i,1),     posi1(i,2),     posi1(i,3), &
                         posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                         natom, natom_x, natom_y, natom_z, natom_v, &
                         Rij, i, j, rc, natom_v )
         ! G2_V_V
         elseif( natom_x + natom_y + natom_z < j &
                 .and. j <= natom_x + natom_y + natom_z + natom_v )then
           call sf_4_g2( g2_v_v, g2_deriv_v_v, eta2_v_v, rs2_v_v, num_g2_v_v, &
                         posi1(i,1),     posi1(i,2),     posi1(i,3), &
                         posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                         natom, natom_x, natom_y, natom_z, natom_v, &
                         Rij, i, j, rc, natom_v )
         else
           write(*,*) "Error(sf_4): G2_V_*"
           stop
         endif

         do ix2 = 1, neighbor

           atom2 = list(i,ix2)

           if( atom2 == 0 )cycle
           if( atom2 - center == i )cycle
           if( atom1 == atom2 )cycle

           k = atom2 - natom * ( ( atom2 - 1 ) / natom )
           Rik = dsqrt(   ( posi1(i,1) - posi2(atom2,1) )**2 &
                        + ( posi1(i,2) - posi2(atom2,2) )**2 &
                        + ( posi1(i,3) - posi2(atom2,3) )**2 )

           if( Rik <= rc )then
             
             ! G5_V_XX
             if( j <= natom_x .and. &
                 k <= natom_x )then
               call sf_4_g5( g5_v_xx, g5_deriv_v_xx, &
                             eta5_v_xx, lambda5_v_xx, zeta5_v_xx, num_g5_v_xx, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_v )
             ! G5_V_XY
             elseif( j <= natom_x .and. &
                     natom_x < k .and. k <= natom_x + natom_y )then
               call sf_4_g5( g5_v_xy, g5_deriv_v_xy, &
                             eta5_v_xy, lambda5_v_xy, zeta5_v_xy, num_g5_v_xy, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_v )
             ! G5_V_XZ
             elseif( j <= natom_x .and. &
                     natom_x + natom_y < k .and. k <= natom_x + natom_y + natom_z )then
               call sf_4_g5( g5_v_xz, g5_deriv_v_xz, &
                             eta5_v_xz, lambda5_v_xz, zeta5_v_xz, num_g5_v_xz, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_v )
             ! G5_V_XV
             elseif( j <= natom_x .and. &
                     natom_x + natom_y + natom_z < k &
                     .and. k <= natom_x + natom_y + natom_z + natom_v )then
               call sf_4_g5( g5_v_xv, g5_deriv_v_xv, &
                             eta5_v_xv, lambda5_v_xv, zeta5_v_xv, num_g5_v_xv, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_v )
             ! G5_V_YY
             elseif( natom_x < j .and. j <= natom_x + natom_y .and. &
                     natom_x < k .and. k <= natom_x + natom_y )then
               call sf_4_g5( g5_v_yy, g5_deriv_v_yy, &
                             eta5_v_yy, lambda5_v_yy, zeta5_v_yy, num_g5_v_yy, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_v )
             ! G5_V_YZ
             elseif( natom_x < j .and. j <= natom_x + natom_y .and. &
                     natom_x + natom_y < k .and. k <= natom_x + natom_y + natom_z )then
               call sf_4_g5( g5_v_yz, g5_deriv_v_yz, &
                             eta5_v_yz, lambda5_v_yz, zeta5_v_yz, num_g5_v_yz, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_v )
             ! G5_V_YV
             elseif( natom_x < j .and. j <= natom_x + natom_y .and. &
                     natom_x + natom_y + natom_z < k &
                     .and. k <= natom_x + natom_y + natom_z + natom_v )then
               call sf_4_g5( g5_v_yv, g5_deriv_v_yv, &
                             eta5_v_yv, lambda5_v_yv, zeta5_v_yv, num_g5_v_yv, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_v )
             ! G5_V_ZZ
             elseif( natom_x + natom_y < j .and. j <= natom_x + natom_y + natom_z .and. &
                     natom_x + natom_y < k .and. k <= natom_x + natom_y + natom_z )then
               call sf_4_g5( g5_v_zz, g5_deriv_v_zz, &
                             eta5_v_zz, lambda5_v_zz, zeta5_v_zz, num_g5_v_zz, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_v )
             ! G5_V_ZV
             elseif( natom_x + natom_y < j .and. j <= natom_x + natom_y + natom_z .and. &
                     natom_x + natom_y + natom_z < k &
                     .and. k <= natom_x + natom_y + natom_z + natom_v )then
               call sf_4_g5( g5_v_zv, g5_deriv_v_zv, &
                             eta5_v_zv, lambda5_v_zv, zeta5_v_zv, num_g5_v_zv, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_v )
             ! G5_V_VV
             elseif( natom_x + natom_y + natom_z < j &
                     .and. j <= natom_x + natom_y + natom_z + natom_v .and. &
                     natom_x + natom_y + natom_z < k &
                     .and. k <= natom_x + natom_y + natom_z + natom_v )then
               call sf_4_g5( g5_v_vv, g5_deriv_v_vv, &
                             eta5_v_vv, lambda5_v_vv, zeta5_v_vv, num_g5_v_vv, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, natom_z, natom_v, &
                             Rij, Rik, i, j, k, rc, natom_v )
             ! G5_V_YX
             elseif( natom_x < j .and. j <= natom_x + natom_y .and. &
                     k <= natom_x )then
               cycle
             ! G5_V_ZX
             elseif( natom_x + natom_y < j .and. j <= natom_x + natom_y + natom_z .and. &
                     k <= natom_x )then
               cycle
             ! G5_V_ZY
             elseif( natom_x + natom_y < j .and. j <= natom_x + natom_y + natom_z .and. &
                     natom_x < k .and. k <= natom_x + natom_y )then
               cycle
             ! G5_V_VX
             elseif( natom_x + natom_y + natom_z < j &
                     .and. j <= natom_x + natom_y + natom_z + natom_v .and. &
                     k <= natom_x )then
               cycle
             ! G5_V_VY
             elseif( natom_x + natom_y + natom_z < j &
                     .and. j <= natom_x + natom_y + natom_z + natom_v .and. &
                     natom_x < k .and. k <= natom_x + natom_y )then
               cycle
             ! G5_V_VZ
             elseif( natom_x + natom_y + natom_z < j &
                     .and. j <= natom_x + natom_y + natom_z + natom_v .and. &
                     natom_x + natom_y < k .and. k <= natom_x + natom_y + natom_z )then
               cycle
             else
               write(*,*) "Error(sf_4): G5_V_*"
               stop
             endif

           endif ! if( Rik <= rc )
         enddo ! do ix2
       endif ! if( Rij < rc )
     enddo ! do ix1 (g2_v... & g5_v...)

   endif ! if( i <= natom_x )

   endif ! i <= natom
 enddo ! i

 call mpi_barrier( mpi_comm_world, ierr )


!*****************
!Coefficient of g5
!*****************
 do i0 = 1, natom, num_mpi
   i = i0 + myrank
   if( i <= natom )then

!  A
   if( i <= natom_x )then

     do ig5 = 1, num_g5_x_xx
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_x_xx(ig5) )
         g5_x_xx(i,ig5) = g5_x_xx(i,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_x_xx(i,j,1,ig5) = g5_deriv_x_xx(i,j,1,ig5) * g5_const_3
         g5_deriv_x_xx(i,j,2,ig5) = g5_deriv_x_xx(i,j,2,ig5) * g5_const_3
         g5_deriv_x_xx(i,j,3,ig5) = g5_deriv_x_xx(i,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_x_xy
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_x_xy(ig5) )
         g5_x_xy(i,ig5) = g5_x_xy(i,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_x_xy(i,j,1,ig5) = g5_deriv_x_xy(i,j,1,ig5) * g5_const_3
         g5_deriv_x_xy(i,j,2,ig5) = g5_deriv_x_xy(i,j,2,ig5) * g5_const_3
         g5_deriv_x_xy(i,j,3,ig5) = g5_deriv_x_xy(i,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_x_xz
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_x_xz(ig5) )
         g5_x_xz(i,ig5) = g5_x_xz(i,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_x_xz(i,j,1,ig5) = g5_deriv_x_xz(i,j,1,ig5) * g5_const_3
         g5_deriv_x_xz(i,j,2,ig5) = g5_deriv_x_xz(i,j,2,ig5) * g5_const_3
         g5_deriv_x_xz(i,j,3,ig5) = g5_deriv_x_xz(i,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_x_xv
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_x_xv(ig5) )
         g5_x_xv(i,ig5) = g5_x_xv(i,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_x_xv(i,j,1,ig5) = g5_deriv_x_xv(i,j,1,ig5) * g5_const_3
         g5_deriv_x_xv(i,j,2,ig5) = g5_deriv_x_xv(i,j,2,ig5) * g5_const_3
         g5_deriv_x_xv(i,j,3,ig5) = g5_deriv_x_xv(i,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_x_yy
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_x_yy(ig5) )
         g5_x_yy(i,ig5) = g5_x_yy(i,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_x_yy(i,j,1,ig5) = g5_deriv_x_yy(i,j,1,ig5) * g5_const_3
         g5_deriv_x_yy(i,j,2,ig5) = g5_deriv_x_yy(i,j,2,ig5) * g5_const_3
         g5_deriv_x_yy(i,j,3,ig5) = g5_deriv_x_yy(i,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_x_yz
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_x_yz(ig5) )
         g5_x_yz(i,ig5) = g5_x_yz(i,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_x_yz(i,j,1,ig5) = g5_deriv_x_yz(i,j,1,ig5) * g5_const_3
         g5_deriv_x_yz(i,j,2,ig5) = g5_deriv_x_yz(i,j,2,ig5) * g5_const_3
         g5_deriv_x_yz(i,j,3,ig5) = g5_deriv_x_yz(i,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_x_yv
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_x_yv(ig5) )
         g5_x_yv(i,ig5) = g5_x_yv(i,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_x_yv(i,j,1,ig5) = g5_deriv_x_yv(i,j,1,ig5) * g5_const_3
         g5_deriv_x_yv(i,j,2,ig5) = g5_deriv_x_yv(i,j,2,ig5) * g5_const_3
         g5_deriv_x_yv(i,j,3,ig5) = g5_deriv_x_yv(i,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_x_zz
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_x_zz(ig5) )
         g5_x_zz(i,ig5) = g5_x_zz(i,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_x_zz(i,j,1,ig5) = g5_deriv_x_zz(i,j,1,ig5) * g5_const_3
         g5_deriv_x_zz(i,j,2,ig5) = g5_deriv_x_zz(i,j,2,ig5) * g5_const_3
         g5_deriv_x_zz(i,j,3,ig5) = g5_deriv_x_zz(i,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_x_zv
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_x_zv(ig5) )
         g5_x_zv(i,ig5) = g5_x_zv(i,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_x_zv(i,j,1,ig5) = g5_deriv_x_zv(i,j,1,ig5) * g5_const_3
         g5_deriv_x_zv(i,j,2,ig5) = g5_deriv_x_zv(i,j,2,ig5) * g5_const_3
         g5_deriv_x_zv(i,j,3,ig5) = g5_deriv_x_zv(i,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_x_vv
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_x_vv(ig5) )
         g5_x_vv(i,ig5) = g5_x_vv(i,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_x_vv(i,j,1,ig5) = g5_deriv_x_vv(i,j,1,ig5) * g5_const_3
         g5_deriv_x_vv(i,j,2,ig5) = g5_deriv_x_vv(i,j,2,ig5) * g5_const_3
         g5_deriv_x_vv(i,j,3,ig5) = g5_deriv_x_vv(i,j,3,ig5) * g5_const_3
       enddo
     enddo

!  B
   elseif( natom_x < i .and. i <= natom_x + natom_y )then

     ii = i - natom_x
  
     do ig5 = 1, num_g5_y_xx
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_y_xx(ig5) )
         g5_y_xx(ii,ig5) = g5_y_xx(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_y_xx(ii,j,1,ig5) = g5_deriv_y_xx(ii,j,1,ig5) * g5_const_3
         g5_deriv_y_xx(ii,j,2,ig5) = g5_deriv_y_xx(ii,j,2,ig5) * g5_const_3
         g5_deriv_y_xx(ii,j,3,ig5) = g5_deriv_y_xx(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_y_xy
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_y_xy(ig5) )
         g5_y_xy(ii,ig5) = g5_y_xy(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_y_xy(ii,j,1,ig5) = g5_deriv_y_xy(ii,j,1,ig5) * g5_const_3
         g5_deriv_y_xy(ii,j,2,ig5) = g5_deriv_y_xy(ii,j,2,ig5) * g5_const_3
         g5_deriv_y_xy(ii,j,3,ig5) = g5_deriv_y_xy(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_y_xz
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_y_xz(ig5) )
         g5_y_xz(ii,ig5) = g5_y_xz(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_y_xz(ii,j,1,ig5) = g5_deriv_y_xz(ii,j,1,ig5) * g5_const_3
         g5_deriv_y_xz(ii,j,2,ig5) = g5_deriv_y_xz(ii,j,2,ig5) * g5_const_3
         g5_deriv_y_xz(ii,j,3,ig5) = g5_deriv_y_xz(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_y_xv
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_y_xv(ig5) )
         g5_y_xv(ii,ig5) = g5_y_xv(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_y_xv(ii,j,1,ig5) = g5_deriv_y_xv(ii,j,1,ig5) * g5_const_3
         g5_deriv_y_xv(ii,j,2,ig5) = g5_deriv_y_xv(ii,j,2,ig5) * g5_const_3
         g5_deriv_y_xv(ii,j,3,ig5) = g5_deriv_y_xv(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_y_yy
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_y_yy(ig5) )
         g5_y_yy(ii,ig5) = g5_y_yy(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_y_yy(ii,j,1,ig5) = g5_deriv_y_yy(ii,j,1,ig5) * g5_const_3
         g5_deriv_y_yy(ii,j,2,ig5) = g5_deriv_y_yy(ii,j,2,ig5) * g5_const_3
         g5_deriv_y_yy(ii,j,3,ig5) = g5_deriv_y_yy(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_y_yz
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_y_yz(ig5) )
         g5_y_yz(ii,ig5) = g5_y_yz(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_y_yz(ii,j,1,ig5) = g5_deriv_y_yz(ii,j,1,ig5) * g5_const_3
         g5_deriv_y_yz(ii,j,2,ig5) = g5_deriv_y_yz(ii,j,2,ig5) * g5_const_3
         g5_deriv_y_yz(ii,j,3,ig5) = g5_deriv_y_yz(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_y_yv
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_y_yv(ig5) )
         g5_y_yv(ii,ig5) = g5_y_yv(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_y_yv(ii,j,1,ig5) = g5_deriv_y_yv(ii,j,1,ig5) * g5_const_3
         g5_deriv_y_yv(ii,j,2,ig5) = g5_deriv_y_yv(ii,j,2,ig5) * g5_const_3
         g5_deriv_y_yv(ii,j,3,ig5) = g5_deriv_y_yv(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_y_zz
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_y_zz(ig5) )
         g5_y_zz(ii,ig5) = g5_y_zz(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_y_zz(ii,j,1,ig5) = g5_deriv_y_zz(ii,j,1,ig5) * g5_const_3
         g5_deriv_y_zz(ii,j,2,ig5) = g5_deriv_y_zz(ii,j,2,ig5) * g5_const_3
         g5_deriv_y_zz(ii,j,3,ig5) = g5_deriv_y_zz(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_y_zv
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_y_zv(ig5) )
         g5_y_zv(ii,ig5) = g5_y_zv(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_y_zv(ii,j,1,ig5) = g5_deriv_y_zv(ii,j,1,ig5) * g5_const_3
         g5_deriv_y_zv(ii,j,2,ig5) = g5_deriv_y_zv(ii,j,2,ig5) * g5_const_3
         g5_deriv_y_zv(ii,j,3,ig5) = g5_deriv_y_zv(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_y_vv
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_y_vv(ig5) )
         g5_y_vv(ii,ig5) = g5_y_vv(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_y_vv(ii,j,1,ig5) = g5_deriv_y_vv(ii,j,1,ig5) * g5_const_3
         g5_deriv_y_vv(ii,j,2,ig5) = g5_deriv_y_vv(ii,j,2,ig5) * g5_const_3
         g5_deriv_y_vv(ii,j,3,ig5) = g5_deriv_y_vv(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

!  C
   elseif( natom_x + natom_y < i .and. i <= natom_x + natom_y + natom_z )then

     ii = i - natom_x - natom_y

     do ig5 = 1, num_g5_z_xx
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_z_xx(ig5) )
         g5_z_xx(ii,ig5) = g5_z_xx(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_z_xx(ii,j,1,ig5) = g5_deriv_z_xx(ii,j,1,ig5) * g5_const_3
         g5_deriv_z_xx(ii,j,2,ig5) = g5_deriv_z_xx(ii,j,2,ig5) * g5_const_3
         g5_deriv_z_xx(ii,j,3,ig5) = g5_deriv_z_xx(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_z_xy
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_z_xy(ig5) )
         g5_z_xy(ii,ig5) = g5_z_xy(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_z_xy(ii,j,1,ig5) = g5_deriv_z_xy(ii,j,1,ig5) * g5_const_3
         g5_deriv_z_xy(ii,j,2,ig5) = g5_deriv_z_xy(ii,j,2,ig5) * g5_const_3
         g5_deriv_z_xy(ii,j,3,ig5) = g5_deriv_z_xy(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_z_xz
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_z_xz(ig5) )
         g5_z_xz(ii,ig5) = g5_z_xz(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_z_xz(ii,j,1,ig5) = g5_deriv_z_xz(ii,j,1,ig5) * g5_const_3
         g5_deriv_z_xz(ii,j,2,ig5) = g5_deriv_z_xz(ii,j,2,ig5) * g5_const_3
         g5_deriv_z_xz(ii,j,3,ig5) = g5_deriv_z_xz(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_z_xv
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_z_xv(ig5) )
         g5_z_xv(ii,ig5) = g5_z_xv(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_z_xv(ii,j,1,ig5) = g5_deriv_z_xv(ii,j,1,ig5) * g5_const_3
         g5_deriv_z_xv(ii,j,2,ig5) = g5_deriv_z_xv(ii,j,2,ig5) * g5_const_3
         g5_deriv_z_xv(ii,j,3,ig5) = g5_deriv_z_xv(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_z_yy
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_z_yy(ig5) )
         g5_z_yy(ii,ig5) = g5_z_yy(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_z_yy(ii,j,1,ig5) = g5_deriv_z_yy(ii,j,1,ig5) * g5_const_3
         g5_deriv_z_yy(ii,j,2,ig5) = g5_deriv_z_yy(ii,j,2,ig5) * g5_const_3
         g5_deriv_z_yy(ii,j,3,ig5) = g5_deriv_z_yy(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_z_yz
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_z_yz(ig5) )
         g5_z_yz(ii,ig5) = g5_z_yz(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_z_yz(ii,j,1,ig5) = g5_deriv_z_yz(ii,j,1,ig5) * g5_const_3
         g5_deriv_z_yz(ii,j,2,ig5) = g5_deriv_z_yz(ii,j,2,ig5) * g5_const_3
         g5_deriv_z_yz(ii,j,3,ig5) = g5_deriv_z_yz(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_z_yv
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_z_yv(ig5) )
         g5_z_yv(ii,ig5) = g5_z_yv(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_z_yv(ii,j,1,ig5) = g5_deriv_z_yv(ii,j,1,ig5) * g5_const_3
         g5_deriv_z_yv(ii,j,2,ig5) = g5_deriv_z_yv(ii,j,2,ig5) * g5_const_3
         g5_deriv_z_yv(ii,j,3,ig5) = g5_deriv_z_yv(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_z_zz
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_z_zz(ig5) )
         g5_z_zz(ii,ig5) = g5_z_zz(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_z_zz(ii,j,1,ig5) = g5_deriv_z_zz(ii,j,1,ig5) * g5_const_3
         g5_deriv_z_zz(ii,j,2,ig5) = g5_deriv_z_zz(ii,j,2,ig5) * g5_const_3
         g5_deriv_z_zz(ii,j,3,ig5) = g5_deriv_z_zz(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_z_zv
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_z_zv(ig5) )
         g5_z_zv(ii,ig5) = g5_z_zv(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_z_zv(ii,j,1,ig5) = g5_deriv_z_zv(ii,j,1,ig5) * g5_const_3
         g5_deriv_z_zv(ii,j,2,ig5) = g5_deriv_z_zv(ii,j,2,ig5) * g5_const_3
         g5_deriv_z_zv(ii,j,3,ig5) = g5_deriv_z_zv(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_z_vv
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_z_vv(ig5) )
         g5_z_vv(ii,ig5) = g5_z_vv(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_z_vv(ii,j,1,ig5) = g5_deriv_z_vv(ii,j,1,ig5) * g5_const_3
         g5_deriv_z_vv(ii,j,2,ig5) = g5_deriv_z_vv(ii,j,2,ig5) * g5_const_3
         g5_deriv_z_vv(ii,j,3,ig5) = g5_deriv_z_vv(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

!  D
   elseif( natom_x + natom_y + natom_z < i &
           .and. i <= natom_x + natom_y + natom_z + natom_v )then

     ii = i - natom_x - natom_y - natom_z

     do ig5 = 1, num_g5_v_xx
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_v_xx(ig5) )
         g5_v_xx(ii,ig5) = g5_v_xx(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_v_xx(ii,j,1,ig5) = g5_deriv_v_xx(ii,j,1,ig5) * g5_const_3
         g5_deriv_v_xx(ii,j,2,ig5) = g5_deriv_v_xx(ii,j,2,ig5) * g5_const_3
         g5_deriv_v_xx(ii,j,3,ig5) = g5_deriv_v_xx(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_v_xy
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_v_xy(ig5) )
         g5_v_xy(ii,ig5) = g5_v_xy(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_v_xy(ii,j,1,ig5) = g5_deriv_v_xy(ii,j,1,ig5) * g5_const_3
         g5_deriv_v_xy(ii,j,2,ig5) = g5_deriv_v_xy(ii,j,2,ig5) * g5_const_3
         g5_deriv_v_xy(ii,j,3,ig5) = g5_deriv_v_xy(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_v_xz
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_v_xz(ig5) )
         g5_v_xz(ii,ig5) = g5_v_xz(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_v_xz(ii,j,1,ig5) = g5_deriv_v_xz(ii,j,1,ig5) * g5_const_3
         g5_deriv_v_xz(ii,j,2,ig5) = g5_deriv_v_xz(ii,j,2,ig5) * g5_const_3
         g5_deriv_v_xz(ii,j,3,ig5) = g5_deriv_v_xz(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_v_xv
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_v_xv(ig5) )
         g5_v_xv(ii,ig5) = g5_v_xv(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_v_xv(ii,j,1,ig5) = g5_deriv_v_xv(ii,j,1,ig5) * g5_const_3
         g5_deriv_v_xv(ii,j,2,ig5) = g5_deriv_v_xv(ii,j,2,ig5) * g5_const_3
         g5_deriv_v_xv(ii,j,3,ig5) = g5_deriv_v_xv(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_v_yy
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_v_yy(ig5) )
         g5_v_yy(ii,ig5) = g5_v_yy(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_v_yy(ii,j,1,ig5) = g5_deriv_v_yy(ii,j,1,ig5) * g5_const_3
         g5_deriv_v_yy(ii,j,2,ig5) = g5_deriv_v_yy(ii,j,2,ig5) * g5_const_3
         g5_deriv_v_yy(ii,j,3,ig5) = g5_deriv_v_yy(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_v_yz
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_v_yz(ig5) )
         g5_v_yz(ii,ig5) = g5_v_yz(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_v_yz(ii,j,1,ig5) = g5_deriv_v_yz(ii,j,1,ig5) * g5_const_3
         g5_deriv_v_yz(ii,j,2,ig5) = g5_deriv_v_yz(ii,j,2,ig5) * g5_const_3
         g5_deriv_v_yz(ii,j,3,ig5) = g5_deriv_v_yz(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_v_yv
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_v_yv(ig5) )
         g5_v_yv(ii,ig5) = g5_v_yv(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_v_yv(ii,j,1,ig5) = g5_deriv_v_yv(ii,j,1,ig5) * g5_const_3
         g5_deriv_v_yv(ii,j,2,ig5) = g5_deriv_v_yv(ii,j,2,ig5) * g5_const_3
         g5_deriv_v_yv(ii,j,3,ig5) = g5_deriv_v_yv(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_v_zz
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_v_zz(ig5) )
         g5_v_zz(ii,ig5) = g5_v_zz(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_v_zz(ii,j,1,ig5) = g5_deriv_v_zz(ii,j,1,ig5) * g5_const_3
         g5_deriv_v_zz(ii,j,2,ig5) = g5_deriv_v_zz(ii,j,2,ig5) * g5_const_3
         g5_deriv_v_zz(ii,j,3,ig5) = g5_deriv_v_zz(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_v_zv
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_v_zv(ig5) )
         g5_v_zv(ii,ig5) = g5_v_zv(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_v_zv(ii,j,1,ig5) = g5_deriv_v_zv(ii,j,1,ig5) * g5_const_3
         g5_deriv_v_zv(ii,j,2,ig5) = g5_deriv_v_zv(ii,j,2,ig5) * g5_const_3
         g5_deriv_v_zv(ii,j,3,ig5) = g5_deriv_v_zv(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_v_vv
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_v_vv(ig5) )
         g5_v_vv(ii,ig5) = g5_v_vv(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_v_vv(ii,j,1,ig5) = g5_deriv_v_vv(ii,j,1,ig5) * g5_const_3
         g5_deriv_v_vv(ii,j,2,ig5) = g5_deriv_v_vv(ii,j,2,ig5) * g5_const_3
         g5_deriv_v_vv(ii,j,3,ig5) = g5_deriv_v_vv(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

   endif ! i

   endif ! i < natom

 enddo ! i0

 call mpi_barrier( mpi_comm_world, ierr )

 call mpi_reduce( g2_x_x, g2_x_x_mpi, natom_x*num_g2_x_x, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_x_y, g2_x_y_mpi, natom_x*num_g2_x_y, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_x_z, g2_x_z_mpi, natom_x*num_g2_x_z, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_x_v, g2_x_v_mpi, natom_x*num_g2_x_v, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_x_xx, g5_x_xx_mpi, natom_x*num_g5_x_xx, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_x_xy, g5_x_xy_mpi, natom_x*num_g5_x_xy, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_x_xz, g5_x_xz_mpi, natom_x*num_g5_x_xz, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_x_xv, g5_x_xv_mpi, natom_x*num_g5_x_xv, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_x_yy, g5_x_yy_mpi, natom_x*num_g5_x_yy, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_x_yz, g5_x_yz_mpi, natom_x*num_g5_x_yz, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_x_yv, g5_x_yv_mpi, natom_x*num_g5_x_yv, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_x_zz, g5_x_zz_mpi, natom_x*num_g5_x_zz, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_x_zv, g5_x_zv_mpi, natom_x*num_g5_x_zv, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_x_vv, g5_x_vv_mpi, natom_x*num_g5_x_vv, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )

 call mpi_reduce( g2_y_x, g2_y_x_mpi, natom_y*num_g2_y_x, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_y_y, g2_y_y_mpi, natom_y*num_g2_y_y, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_y_z, g2_y_z_mpi, natom_y*num_g2_y_z, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_y_v, g2_y_v_mpi, natom_y*num_g2_y_v, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_y_xx, g5_y_xx_mpi, natom_y*num_g5_y_xx, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_y_xy, g5_y_xy_mpi, natom_y*num_g5_y_xy, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_y_xz, g5_y_xz_mpi, natom_y*num_g5_y_xz, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_y_xv, g5_y_xv_mpi, natom_y*num_g5_y_xv, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_y_yy, g5_y_yy_mpi, natom_y*num_g5_y_yy, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_y_yz, g5_y_yz_mpi, natom_y*num_g5_y_yz, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_y_yv, g5_y_yv_mpi, natom_y*num_g5_y_yv, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_y_zz, g5_y_zz_mpi, natom_y*num_g5_y_zz, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_y_zv, g5_y_zv_mpi, natom_y*num_g5_y_zv, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_y_vv, g5_y_vv_mpi, natom_y*num_g5_y_vv, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )

 call mpi_reduce( g2_z_x, g2_z_x_mpi, natom_z*num_g2_z_x, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_z_y, g2_z_y_mpi, natom_z*num_g2_z_y, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_z_z, g2_z_z_mpi, natom_z*num_g2_z_z, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_z_v, g2_z_v_mpi, natom_z*num_g2_z_v, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_z_xx, g5_z_xx_mpi, natom_z*num_g5_z_xx, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_z_xy, g5_z_xy_mpi, natom_z*num_g5_z_xy, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_z_xz, g5_z_xz_mpi, natom_z*num_g5_z_xz, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_z_xv, g5_z_xv_mpi, natom_z*num_g5_z_xv, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_z_yy, g5_z_yy_mpi, natom_z*num_g5_z_yy, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_z_yz, g5_z_yz_mpi, natom_z*num_g5_z_yz, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_z_yv, g5_z_yv_mpi, natom_z*num_g5_z_yv, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_z_zz, g5_z_zz_mpi, natom_z*num_g5_z_zz, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_z_zv, g5_z_zv_mpi, natom_z*num_g5_z_zv, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_z_vv, g5_z_vv_mpi, natom_z*num_g5_z_vv, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )

 call mpi_reduce( g2_v_x, g2_v_x_mpi, natom_v*num_g2_v_x, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_v_y, g2_v_y_mpi, natom_v*num_g2_v_y, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_v_z, g2_v_z_mpi, natom_v*num_g2_v_z, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_v_v, g2_v_v_mpi, natom_v*num_g2_v_v, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_v_xx, g5_v_xx_mpi, natom_v*num_g5_v_xx, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_v_xy, g5_v_xy_mpi, natom_v*num_g5_v_xy, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_v_xz, g5_v_xz_mpi, natom_v*num_g5_v_xz, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_v_xv, g5_v_xv_mpi, natom_v*num_g5_v_xv, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_v_yy, g5_v_yy_mpi, natom_v*num_g5_v_yy, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_v_yz, g5_v_yz_mpi, natom_v*num_g5_v_yz, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_v_yv, g5_v_yv_mpi, natom_v*num_g5_v_yv, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_v_zz, g5_v_zz_mpi, natom_v*num_g5_v_zz, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_v_zv, g5_v_zv_mpi, natom_v*num_g5_v_zv, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_v_vv, g5_v_vv_mpi, natom_v*num_g5_v_vv, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )


 call mpi_reduce( g2_deriv_x_x, g2_deriv_x_x_mpi, natom_x*natom*3*num_g2_x_x, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_deriv_x_y, g2_deriv_x_y_mpi, natom_x*natom*3*num_g2_x_y, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_deriv_x_z, g2_deriv_x_z_mpi, natom_x*natom*3*num_g2_x_z, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_deriv_x_v, g2_deriv_x_v_mpi, natom_x*natom*3*num_g2_x_v, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_x_xx, g5_deriv_x_xx_mpi, natom_x*natom*3*num_g5_x_xx, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_x_xy, g5_deriv_x_xy_mpi, natom_x*natom*3*num_g5_x_xy, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_x_xz, g5_deriv_x_xz_mpi, natom_x*natom*3*num_g5_x_xz, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_x_xv, g5_deriv_x_xv_mpi, natom_x*natom*3*num_g5_x_xv, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_x_yy, g5_deriv_x_yy_mpi, natom_x*natom*3*num_g5_x_yy, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_x_yz, g5_deriv_x_yz_mpi, natom_x*natom*3*num_g5_x_yz, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_x_yv, g5_deriv_x_yv_mpi, natom_x*natom*3*num_g5_x_yv, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_x_zz, g5_deriv_x_zz_mpi, natom_x*natom*3*num_g5_x_zz, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_x_zv, g5_deriv_x_zv_mpi, natom_x*natom*3*num_g5_x_zv, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_x_vv, g5_deriv_x_vv_mpi, natom_x*natom*3*num_g5_x_vv, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )

 call mpi_reduce( g2_deriv_y_x, g2_deriv_y_x_mpi, natom_y*natom*3*num_g2_y_x, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_deriv_y_y, g2_deriv_y_y_mpi, natom_y*natom*3*num_g2_y_y, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_deriv_y_z, g2_deriv_y_z_mpi, natom_y*natom*3*num_g2_y_z, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_deriv_y_v, g2_deriv_y_v_mpi, natom_y*natom*3*num_g2_y_v, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_y_xx, g5_deriv_y_xx_mpi, natom_y*natom*3*num_g5_y_xx, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_y_xy, g5_deriv_y_xy_mpi, natom_y*natom*3*num_g5_y_xy, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_y_xz, g5_deriv_y_xz_mpi, natom_y*natom*3*num_g5_y_xz, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_y_xv, g5_deriv_y_xv_mpi, natom_y*natom*3*num_g5_y_xv, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_y_yy, g5_deriv_y_yy_mpi, natom_y*natom*3*num_g5_y_yy, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_y_yz, g5_deriv_y_yz_mpi, natom_y*natom*3*num_g5_y_yz, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_y_yv, g5_deriv_y_yv_mpi, natom_y*natom*3*num_g5_y_yv, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_y_zz, g5_deriv_y_zz_mpi, natom_y*natom*3*num_g5_y_zz, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_y_zv, g5_deriv_y_zv_mpi, natom_y*natom*3*num_g5_y_zv, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_y_vv, g5_deriv_y_vv_mpi, natom_y*natom*3*num_g5_y_vv, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )

 call mpi_reduce( g2_deriv_z_x, g2_deriv_z_x_mpi, natom_z*natom*3*num_g2_z_x, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_deriv_z_y, g2_deriv_z_y_mpi, natom_z*natom*3*num_g2_z_y, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_deriv_z_z, g2_deriv_z_z_mpi, natom_z*natom*3*num_g2_z_z, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_deriv_z_v, g2_deriv_z_v_mpi, natom_z*natom*3*num_g2_z_v, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_z_xx, g5_deriv_z_xx_mpi, natom_z*natom*3*num_g5_z_xx, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_z_xy, g5_deriv_z_xy_mpi, natom_z*natom*3*num_g5_z_xy, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_z_xz, g5_deriv_z_xz_mpi, natom_z*natom*3*num_g5_z_xz, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_z_xv, g5_deriv_z_xv_mpi, natom_z*natom*3*num_g5_z_xv, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_z_yy, g5_deriv_z_yy_mpi, natom_z*natom*3*num_g5_z_yy, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_z_yz, g5_deriv_z_yz_mpi, natom_z*natom*3*num_g5_z_yz, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_z_yv, g5_deriv_z_yv_mpi, natom_z*natom*3*num_g5_z_yv, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_z_zz, g5_deriv_z_zz_mpi, natom_z*natom*3*num_g5_z_zz, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_z_zv, g5_deriv_z_zv_mpi, natom_z*natom*3*num_g5_z_zv, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_z_vv, g5_deriv_z_vv_mpi, natom_z*natom*3*num_g5_z_vv, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )

 call mpi_reduce( g2_deriv_v_x, g2_deriv_v_x_mpi, natom_v*natom*3*num_g2_v_x, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_deriv_v_y, g2_deriv_v_y_mpi, natom_v*natom*3*num_g2_v_y, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_deriv_v_z, g2_deriv_v_z_mpi, natom_v*natom*3*num_g2_v_z, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_deriv_v_v, g2_deriv_v_v_mpi, natom_v*natom*3*num_g2_v_v, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_v_xx, g5_deriv_v_xx_mpi, natom_v*natom*3*num_g5_v_xx, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_v_xy, g5_deriv_v_xy_mpi, natom_v*natom*3*num_g5_v_xy, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_v_xz, g5_deriv_v_xz_mpi, natom_v*natom*3*num_g5_v_xz, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_v_xv, g5_deriv_v_xv_mpi, natom_v*natom*3*num_g5_v_xv, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_v_yy, g5_deriv_v_yy_mpi, natom_v*natom*3*num_g5_v_yy, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_v_yz, g5_deriv_v_yz_mpi, natom_v*natom*3*num_g5_v_yz, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_v_yv, g5_deriv_v_yv_mpi, natom_v*natom*3*num_g5_v_yv, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_v_zz, g5_deriv_v_zz_mpi, natom_v*natom*3*num_g5_v_zz, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_v_zv, g5_deriv_v_zv_mpi, natom_v*natom*3*num_g5_v_zv, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_v_vv, g5_deriv_v_vv_mpi, natom_v*natom*3*num_g5_v_vv, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )

 call mpi_barrier( mpi_comm_world, ierr )


 deallocate( posi2 )


!--------
!Check sf
!--------

!X
do j = 1, num_g2_x_x
  do i = 1, natom_x
    if( isnan( g2_x_x(i,j) ) ) stop 'NaN g2_x_x'
  enddo
enddo

do j = 1, num_g2_x_y
  do i = 1, natom_x
    if( isnan( g2_x_y(i,j) ) ) stop 'NaN g2_x_y'
  enddo
enddo

do j = 1, num_g2_x_z
  do i = 1, natom_x
    if( isnan( g2_x_z(i,j) ) ) stop 'NaN g2_x_z'
  enddo
enddo

do j = 1, num_g2_x_v
  do i = 1, natom_x
    if( isnan( g2_x_v(i,j) ) ) stop 'NaN g2_x_v'
  enddo
enddo

do j = 1, num_g5_x_xx
  do i = 1, natom_x
    if( isnan( g5_x_xx(i,j) ) ) stop 'NaN g5_x_xx'
  enddo
enddo

do j = 1, num_g5_x_xy
  do i = 1, natom_x
    if( isnan( g5_x_xy(i,j) ) ) stop 'NaN g5_x_xy'
  enddo
enddo

do j = 1, num_g5_x_xz
  do i = 1, natom_x
    if( isnan( g5_x_xz(i,j) ) ) stop 'NaN g5_x_xz'
  enddo
enddo

do j = 1, num_g5_x_xv
  do i = 1, natom_x
    if( isnan( g5_x_xv(i,j) ) ) stop 'NaN g5_x_xv'
  enddo
enddo

do j = 1, num_g5_x_yy
  do i = 1, natom_x
    if( isnan( g5_x_yy(i,j) ) ) stop 'NaN g5_x_yy'
  enddo
enddo

do j = 1, num_g5_x_yz
  do i = 1, natom_x
    if( isnan( g5_x_yz(i,j) ) ) stop 'NaN g5_x_yz'
  enddo
enddo

do j = 1, num_g5_x_yv
  do i = 1, natom_x
    if( isnan( g5_x_yv(i,j) ) ) stop 'NaN g5_x_yv'
  enddo
enddo

do j = 1, num_g5_x_zz
  do i = 1, natom_x
    if( isnan( g5_x_zz(i,j) ) ) stop 'NaN g5_x_zz'
  enddo
enddo

do j = 1, num_g5_x_zv
  do i = 1, natom_x
    if( isnan( g5_x_zv(i,j) ) ) stop 'NaN g5_x_zv'
  enddo
enddo

do j = 1, num_g5_x_vv
  do i = 1, natom_x
    if( isnan( g5_x_vv(i,j) ) ) stop 'NaN g5_x_vv'
  enddo
enddo

!Y
do j = 1, num_g2_y_x
  do i = 1, natom_y
    if( isnan( g2_y_x(i,j) ) ) stop 'NaN g2_y_x'
  enddo
enddo

do j = 1, num_g2_y_y
  do i = 1, natom_y
    if( isnan( g2_y_y(i,j) ) ) stop 'NaN g2_y_y'
  enddo
enddo

do j = 1, num_g2_y_z
  do i = 1, natom_y
    if( isnan( g2_y_z(i,j) ) ) stop 'NaN g2_y_z'
  enddo
enddo

do j = 1, num_g2_y_v
  do i = 1, natom_y
    if( isnan( g2_y_v(i,j) ) ) stop 'NaN g2_y_v'
  enddo
enddo

do j = 1, num_g5_y_xx
  do i = 1, natom_y
    if( isnan( g5_y_xx(i,j) ) ) stop 'NaN g5_y_xx'
  enddo
enddo

do j = 1, num_g5_y_xy
  do i = 1, natom_y
    if( isnan( g5_y_xy(i,j) ) ) stop 'NaN g5_y_xy'
  enddo
enddo

do j = 1, num_g5_y_xz
  do i = 1, natom_y
    if( isnan( g5_y_xz(i,j) ) ) stop 'NaN g5_y_xz'
  enddo
enddo

do j = 1, num_g5_y_xv
  do i = 1, natom_y
    if( isnan( g5_y_xv(i,j) ) ) stop 'NaN g5_y_xv'
  enddo
enddo

do j = 1, num_g5_y_yy
  do i = 1, natom_y
    if( isnan( g5_y_yy(i,j) ) ) stop 'NaN g5_y_yy'
  enddo
enddo

do j = 1, num_g5_y_yz
  do i = 1, natom_y
    if( isnan( g5_y_yz(i,j) ) ) stop 'NaN g5_y_yz'
  enddo
enddo

do j = 1, num_g5_y_yv
  do i = 1, natom_y
    if( isnan( g5_y_yv(i,j) ) ) stop 'NaN g5_y_yv'
  enddo
enddo

do j = 1, num_g5_y_zz
  do i = 1, natom_y
    if( isnan( g5_y_zz(i,j) ) ) stop 'NaN g5_y_zz'
  enddo
enddo

do j = 1, num_g5_y_zv
  do i = 1, natom_y
    if( isnan( g5_y_zv(i,j) ) ) stop 'NaN g5_y_zv'
  enddo
enddo

do j = 1, num_g5_y_vv
  do i = 1, natom_y
    if( isnan( g5_y_vv(i,j) ) ) stop 'NaN g5_y_vv'
  enddo
enddo

!Z
do j = 1, num_g2_z_x
  do i = 1, natom_z
    if( isnan( g2_z_x(i,j) ) ) stop 'NaN g2_z_x'
  enddo
enddo

do j = 1, num_g2_z_y
  do i = 1, natom_z
    if( isnan( g2_z_y(i,j) ) ) stop 'NaN g2_z_y'
  enddo
enddo

do j = 1, num_g2_z_z
  do i = 1, natom_z
    if( isnan( g2_z_z(i,j) ) ) stop 'NaN g2_z_z'
  enddo
enddo

do j = 1, num_g2_z_v
  do i = 1, natom_z
    if( isnan( g2_z_v(i,j) ) ) stop 'NaN g2_z_v'
  enddo
enddo

do j = 1, num_g5_z_xx
  do i = 1, natom_z
    if( isnan( g5_z_xx(i,j) ) ) stop 'NaN g5_z_xx'
  enddo
enddo

do j = 1, num_g5_z_xy
  do i = 1, natom_z
    if( isnan( g5_z_xy(i,j) ) ) stop 'NaN g5_z_xy'
  enddo
enddo

do j = 1, num_g5_z_xz
  do i = 1, natom_z
    if( isnan( g5_z_xz(i,j) ) ) stop 'NaN g5_z_xz'
  enddo
enddo

do j = 1, num_g5_z_xv
  do i = 1, natom_z
    if( isnan( g5_z_xv(i,j) ) ) stop 'NaN g5_z_xv'
  enddo
enddo

do j = 1, num_g5_z_yy
  do i = 1, natom_z
    if( isnan( g5_z_yy(i,j) ) ) stop 'NaN g5_z_yy'
  enddo
enddo

do j = 1, num_g5_z_yz
  do i = 1, natom_z
    if( isnan( g5_z_yz(i,j) ) ) stop 'NaN g5_z_yz'
  enddo
enddo

do j = 1, num_g5_z_yv
  do i = 1, natom_z
    if( isnan( g5_z_yv(i,j) ) ) stop 'NaN g5_z_yv'
  enddo
enddo

do j = 1, num_g5_z_zz
  do i = 1, natom_z
    if( isnan( g5_z_zz(i,j) ) ) stop 'NaN g5_z_zz'
  enddo
enddo

do j = 1, num_g5_z_zv
  do i = 1, natom_z
    if( isnan( g5_z_zv(i,j) ) ) stop 'NaN g5_z_zv'
  enddo
enddo

do j = 1, num_g5_z_vv
  do i = 1, natom_z
    if( isnan( g5_z_vv(i,j) ) ) stop 'NaN g5_z_vv'
  enddo
enddo

!V
do j = 1, num_g2_v_x
  do i = 1, natom_v
    if( isnan( g2_v_x(i,j) ) ) stop 'NaN g2_v_x'
  enddo
enddo

do j = 1, num_g2_v_y
  do i = 1, natom_v
    if( isnan( g2_v_y(i,j) ) ) stop 'NaN g2_v_y'
  enddo
enddo

do j = 1, num_g2_v_z
  do i = 1, natom_v
    if( isnan( g2_v_z(i,j) ) ) stop 'NaN g2_v_z'
  enddo
enddo

do j = 1, num_g2_v_v
  do i = 1, natom_v
    if( isnan( g2_v_v(i,j) ) ) stop 'NaN g2_v_v'
  enddo
enddo

do j = 1, num_g5_v_xx
  do i = 1, natom_v
    if( isnan( g5_v_xx(i,j) ) ) stop 'NaN g5_v_xx'
  enddo
enddo

do j = 1, num_g5_v_xy
  do i = 1, natom_v
    if( isnan( g5_v_xy(i,j) ) ) stop 'NaN g5_v_xy'
  enddo
enddo

do j = 1, num_g5_v_xz
  do i = 1, natom_v
    if( isnan( g5_v_xz(i,j) ) ) stop 'NaN g5_v_xz'
  enddo
enddo

do j = 1, num_g5_v_xv
  do i = 1, natom_v
    if( isnan( g5_v_xv(i,j) ) ) stop 'NaN g5_v_xv'
  enddo
enddo

do j = 1, num_g5_v_yy
  do i = 1, natom_v
    if( isnan( g5_v_yy(i,j) ) ) stop 'NaN g5_v_yy'
  enddo
enddo

do j = 1, num_g5_v_yz
  do i = 1, natom_v
    if( isnan( g5_v_yz(i,j) ) ) stop 'NaN g5_v_yz'
  enddo
enddo

do j = 1, num_g5_v_yv
  do i = 1, natom_v
    if( isnan( g5_v_yv(i,j) ) ) stop 'NaN g5_v_yv'
  enddo
enddo

do j = 1, num_g5_v_zz
  do i = 1, natom_v
    if( isnan( g5_v_zz(i,j) ) ) stop 'NaN g5_v_zz'
  enddo
enddo

do j = 1, num_g5_v_zv
  do i = 1, natom_v
    if( isnan( g5_v_zv(i,j) ) ) stop 'NaN g5_v_zv'
  enddo
enddo

do j = 1, num_g5_v_vv
  do i = 1, natom_v
    if( isnan( g5_v_vv(i,j) ) ) stop 'NaN g5_v_vv'
  enddo
enddo


!Deriv X
do l = 1, num_g2_x_x
do k = 1, 3
do j = 1, natom
do i = 1, natom_x
  if( isnan( g2_deriv_x_x(i,j,k,l) ) ) stop 'NaN g2_deriv_x_x'
enddo
enddo
enddo
enddo

do l = 1, num_g2_x_y
do k = 1, 3
do j = 1, natom
do i = 1, natom_x
  if( isnan( g2_deriv_x_y(i,j,k,l) ) ) stop 'NaN g2_deriv_x_y'
enddo
enddo
enddo
enddo

do l = 1, num_g2_x_z
do k = 1, 3
do j = 1, natom
do i = 1, natom_x
  if( isnan( g2_deriv_x_z(i,j,k,l) ) ) stop 'NaN g2_deriv_x_z'
enddo
enddo
enddo
enddo

do l = 1, num_g2_x_v
do k = 1, 3
do j = 1, natom
do i = 1, natom_x
  if( isnan( g2_deriv_x_v(i,j,k,l) ) ) stop 'NaN g2_deriv_x_v'
enddo
enddo
enddo
enddo

do l = 1, num_g5_x_xx
do k = 1, 3
do j = 1, natom
do i = 1, natom_x
  if( isnan( g5_deriv_x_xx(i,j,k,l) ) ) stop 'NaN g5_deriv_x_xx'
enddo
enddo
enddo
enddo

do l = 1, num_g5_x_xy
do k = 1, 3
do j = 1, natom
do i = 1, natom_x
  if( isnan( g5_deriv_x_xy(i,j,k,l) ) ) stop 'NaN g5_deriv_x_xy'
enddo
enddo
enddo
enddo

do l = 1, num_g5_x_xz
do k = 1, 3
do j = 1, natom
do i = 1, natom_x
  if( isnan( g5_deriv_x_xz(i,j,k,l) ) ) stop 'NaN g5_deriv_x_xz'
enddo
enddo
enddo
enddo

do l = 1, num_g5_x_xv
do k = 1, 3
do j = 1, natom
do i = 1, natom_x
  if( isnan( g5_deriv_x_xv(i,j,k,l) ) ) stop 'NaN g5_deriv_x_xv'
enddo
enddo
enddo
enddo

do l = 1, num_g5_x_yy
do k = 1, 3
do j = 1, natom
do i = 1, natom_x
  if( isnan( g5_deriv_x_yy(i,j,k,l) ) ) stop 'NaN g5_deriv_x_yy'
enddo
enddo
enddo
enddo

do l = 1, num_g5_x_yz
do k = 1, 3
do j = 1, natom
do i = 1, natom_x
  if( isnan( g5_deriv_x_yz(i,j,k,l) ) ) stop 'NaN g5_deriv_x_yz'
enddo
enddo
enddo
enddo

do l = 1, num_g5_x_yv
do k = 1, 3
do j = 1, natom
do i = 1, natom_x
  if( isnan( g5_deriv_x_yv(i,j,k,l) ) ) stop 'NaN g5_deriv_x_yv'
enddo
enddo
enddo
enddo

do l = 1, num_g5_x_zz
do k = 1, 3
do j = 1, natom
do i = 1, natom_x
  if( isnan( g5_deriv_x_zz(i,j,k,l) ) ) stop 'NaN g5_deriv_x_zz'
enddo
enddo
enddo
enddo

do l = 1, num_g5_x_zv
do k = 1, 3
do j = 1, natom
do i = 1, natom_x
  if( isnan( g5_deriv_x_zv(i,j,k,l) ) ) stop 'NaN g5_deriv_x_zv'
enddo
enddo
enddo
enddo

do l = 1, num_g5_x_vv
do k = 1, 3
do j = 1, natom
do i = 1, natom_x
  if( isnan( g5_deriv_x_vv(i,j,k,l) ) ) stop 'NaN g5_deriv_x_vv'
enddo
enddo
enddo
enddo

!Deriv Y
do l = 1, num_g2_y_x
do k = 1, 3
do j = 1, natom
do i = 1, natom_y
  if( isnan( g2_deriv_y_x(i,j,k,l) ) ) stop 'NaN g2_deriv_y_x'
enddo
enddo
enddo
enddo

do l = 1, num_g2_y_y
do k = 1, 3
do j = 1, natom
do i = 1, natom_y
  if( isnan( g2_deriv_y_y(i,j,k,l) ) ) stop 'NaN g2_deriv_y_y'
enddo
enddo
enddo
enddo

do l = 1, num_g2_y_z
do k = 1, 3
do j = 1, natom
do i = 1, natom_y
  if( isnan( g2_deriv_y_z(i,j,k,l) ) ) stop 'NaN g2_deriv_y_z'
enddo
enddo
enddo
enddo

do l = 1, num_g2_y_v
do k = 1, 3
do j = 1, natom
do i = 1, natom_y
  if( isnan( g2_deriv_y_v(i,j,k,l) ) ) stop 'NaN g2_deriv_y_v'
enddo
enddo
enddo
enddo

do l = 1, num_g5_y_xx
do k = 1, 3
do j = 1, natom
do i = 1, natom_y
  if( isnan( g5_deriv_y_xx(i,j,k,l) ) ) stop 'NaN g5_deriv_y_xx'
enddo
enddo
enddo
enddo

do l = 1, num_g5_y_xy
do k = 1, 3
do j = 1, natom
do i = 1, natom_y
  if( isnan( g5_deriv_y_xy(i,j,k,l) ) ) stop 'NaN g5_deriv_y_xy'
enddo
enddo
enddo
enddo

do l = 1, num_g5_y_xz
do k = 1, 3
do j = 1, natom
do i = 1, natom_y
  if( isnan( g5_deriv_y_xz(i,j,k,l) ) ) stop 'NaN g5_deriv_y_xz'
enddo
enddo
enddo
enddo

do l = 1, num_g5_y_xv
do k = 1, 3
do j = 1, natom
do i = 1, natom_y
  if( isnan( g5_deriv_y_xv(i,j,k,l) ) ) stop 'NaN g5_deriv_y_xv'
enddo
enddo
enddo
enddo

do l = 1, num_g5_y_yy
do k = 1, 3
do j = 1, natom
do i = 1, natom_y
  if( isnan( g5_deriv_y_yy(i,j,k,l) ) ) stop 'NaN g5_deriv_y_yy'
enddo
enddo
enddo
enddo

do l = 1, num_g5_y_yz
do k = 1, 3
do j = 1, natom
do i = 1, natom_y
  if( isnan( g5_deriv_y_yz(i,j,k,l) ) ) stop 'NaN g5_deriv_y_yz'
enddo
enddo
enddo
enddo

do l = 1, num_g5_y_yv
do k = 1, 3
do j = 1, natom
do i = 1, natom_y
  if( isnan( g5_deriv_y_yv(i,j,k,l) ) ) stop 'NaN g5_deriv_y_yv'
enddo
enddo
enddo
enddo

do l = 1, num_g5_y_zz
do k = 1, 3
do j = 1, natom
do i = 1, natom_y
  if( isnan( g5_deriv_y_zz(i,j,k,l) ) ) stop 'NaN g5_deriv_y_zz'
enddo
enddo
enddo
enddo

do l = 1, num_g5_y_zv
do k = 1, 3
do j = 1, natom
do i = 1, natom_y
  if( isnan( g5_deriv_y_zv(i,j,k,l) ) ) stop 'NaN g5_deriv_y_zv'
enddo
enddo
enddo
enddo

do l = 1, num_g5_y_vv
do k = 1, 3
do j = 1, natom
do i = 1, natom_y
  if( isnan( g5_deriv_y_vv(i,j,k,l) ) ) stop 'NaN g5_deriv_y_vv'
enddo
enddo
enddo
enddo

!Deriv Z
do l = 1, num_g2_z_x
do k = 1, 3
do j = 1, natom
do i = 1, natom_z
  if( isnan( g2_deriv_z_x(i,j,k,l) ) ) stop 'NaN g2_deriv_z_x'
enddo
enddo
enddo
enddo

do l = 1, num_g2_z_y
do k = 1, 3
do j = 1, natom
do i = 1, natom_z
  if( isnan( g2_deriv_z_y(i,j,k,l) ) ) stop 'NaN g2_deriv_z_y'
enddo
enddo
enddo
enddo

do l = 1, num_g2_z_z
do k = 1, 3
do j = 1, natom
do i = 1, natom_z
  if( isnan( g2_deriv_z_z(i,j,k,l) ) ) stop 'NaN g2_deriv_z_z'
enddo
enddo
enddo
enddo

do l = 1, num_g2_z_v
do k = 1, 3
do j = 1, natom
do i = 1, natom_z
  if( isnan( g2_deriv_z_v(i,j,k,l) ) ) stop 'NaN g2_deriv_z_v'
enddo
enddo
enddo
enddo

do l = 1, num_g5_z_xx
do k = 1, 3
do j = 1, natom
do i = 1, natom_z
  if( isnan( g5_deriv_z_xx(i,j,k,l) ) ) stop 'NaN g5_deriv_z_xx'
enddo
enddo
enddo
enddo

do l = 1, num_g5_z_xy
do k = 1, 3
do j = 1, natom
do i = 1, natom_z
  if( isnan( g5_deriv_z_xy(i,j,k,l) ) ) stop 'NaN g5_deriv_z_xy'
enddo
enddo
enddo
enddo

do l = 1, num_g5_z_xz
do k = 1, 3
do j = 1, natom
do i = 1, natom_z
  if( isnan( g5_deriv_z_xz(i,j,k,l) ) ) stop 'NaN g5_deriv_z_xz'
enddo
enddo
enddo
enddo

do l = 1, num_g5_z_xv
do k = 1, 3
do j = 1, natom
do i = 1, natom_z
  if( isnan( g5_deriv_z_xv(i,j,k,l) ) ) stop 'NaN g5_deriv_z_xv'
enddo
enddo
enddo
enddo

do l = 1, num_g5_z_yy
do k = 1, 3
do j = 1, natom
do i = 1, natom_z
  if( isnan( g5_deriv_z_yy(i,j,k,l) ) ) stop 'NaN g5_deriv_z_yy'
enddo
enddo
enddo
enddo

do l = 1, num_g5_z_yz
do k = 1, 3
do j = 1, natom
do i = 1, natom_z
  if( isnan( g5_deriv_z_yz(i,j,k,l) ) ) stop 'NaN g5_deriv_z_yz'
enddo
enddo
enddo
enddo

do l = 1, num_g5_z_yv
do k = 1, 3
do j = 1, natom
do i = 1, natom_z
  if( isnan( g5_deriv_z_yv(i,j,k,l) ) ) stop 'NaN g5_deriv_z_yv'
enddo
enddo
enddo
enddo

do l = 1, num_g5_z_zz
do k = 1, 3
do j = 1, natom
do i = 1, natom_z
  if( isnan( g5_deriv_z_zz(i,j,k,l) ) ) stop 'NaN g5_deriv_z_zz'
enddo
enddo
enddo
enddo

do l = 1, num_g5_z_zv
do k = 1, 3
do j = 1, natom
do i = 1, natom_z
  if( isnan( g5_deriv_z_zv(i,j,k,l) ) ) stop 'NaN g5_deriv_z_zv'
enddo
enddo
enddo
enddo

do l = 1, num_g5_z_vv
do k = 1, 3
do j = 1, natom
do i = 1, natom_z
  if( isnan( g5_deriv_z_vv(i,j,k,l) ) ) stop 'NaN g5_deriv_z_vv'
enddo
enddo
enddo
enddo

!Deriv V
do l = 1, num_g2_v_x
do k = 1, 3
do j = 1, natom
do i = 1, natom_v
  if( isnan( g2_deriv_v_x(i,j,k,l) ) ) stop 'NaN g2_deriv_v_x'
enddo
enddo
enddo
enddo

do l = 1, num_g2_v_y
do k = 1, 3
do j = 1, natom
do i = 1, natom_v
  if( isnan( g2_deriv_v_y(i,j,k,l) ) ) stop 'NaN g2_deriv_v_y'
enddo
enddo
enddo
enddo

do l = 1, num_g2_v_z
do k = 1, 3
do j = 1, natom
do i = 1, natom_v
  if( isnan( g2_deriv_v_z(i,j,k,l) ) ) stop 'NaN g2_deriv_v_z'
enddo
enddo
enddo
enddo

do l = 1, num_g2_v_v
do k = 1, 3
do j = 1, natom
do i = 1, natom_v
  if( isnan( g2_deriv_v_v(i,j,k,l) ) ) stop 'NaN g2_deriv_v_v'
enddo
enddo
enddo
enddo

do l = 1, num_g5_v_xx
do k = 1, 3
do j = 1, natom
do i = 1, natom_v
  if( isnan( g5_deriv_v_xx(i,j,k,l) ) ) stop 'NaN g5_deriv_v_xx'
enddo
enddo
enddo
enddo

do l = 1, num_g5_v_xy
do k = 1, 3
do j = 1, natom
do i = 1, natom_v
  if( isnan( g5_deriv_v_xy(i,j,k,l) ) ) stop 'NaN g5_deriv_v_xy'
enddo
enddo
enddo
enddo

do l = 1, num_g5_v_xz
do k = 1, 3
do j = 1, natom
do i = 1, natom_v
  if( isnan( g5_deriv_v_xz(i,j,k,l) ) ) stop 'NaN g5_deriv_v_xz'
enddo
enddo
enddo
enddo

do l = 1, num_g5_v_xv
do k = 1, 3
do j = 1, natom
do i = 1, natom_v
  if( isnan( g5_deriv_v_xv(i,j,k,l) ) ) stop 'NaN g5_deriv_v_xv'
enddo
enddo
enddo
enddo

do l = 1, num_g5_v_yy
do k = 1, 3
do j = 1, natom
do i = 1, natom_v
  if( isnan( g5_deriv_v_yy(i,j,k,l) ) ) stop 'NaN g5_deriv_v_yy'
enddo
enddo
enddo
enddo

do l = 1, num_g5_v_yz
do k = 1, 3
do j = 1, natom
do i = 1, natom_v
  if( isnan( g5_deriv_v_yz(i,j,k,l) ) ) stop 'NaN g5_deriv_v_yz'
enddo
enddo
enddo
enddo

do l = 1, num_g5_v_yv
do k = 1, 3
do j = 1, natom
do i = 1, natom_v
  if( isnan( g5_deriv_v_yv(i,j,k,l) ) ) stop 'NaN g5_deriv_v_yv'
enddo
enddo
enddo
enddo

do l = 1, num_g5_v_zz
do k = 1, 3
do j = 1, natom
do i = 1, natom_v
  if( isnan( g5_deriv_v_zz(i,j,k,l) ) ) stop 'NaN g5_deriv_v_zz'
enddo
enddo
enddo
enddo

do l = 1, num_g5_v_zv
do k = 1, 3
do j = 1, natom
do i = 1, natom_v
  if( isnan( g5_deriv_v_zv(i,j,k,l) ) ) stop 'NaN g5_deriv_v_zv'
enddo
enddo
enddo
enddo

do l = 1, num_g5_v_vv
do k = 1, 3
do j = 1, natom
do i = 1, natom_v
  if( isnan( g5_deriv_v_vv(i,j,k,l) ) ) stop 'NaN g5_deriv_v_vv'
enddo
enddo
enddo
enddo


if( myrank == 0 )then

!--------
!Write sf
!--------

!g2_x-... & g5_x_...
 write(11)g2_x_x_mpi
 write(11)g2_x_y_mpi
 write(11)g2_x_z_mpi
 write(11)g2_x_v_mpi
 write(11)g5_x_xx_mpi
 write(11)g5_x_xy_mpi
 write(11)g5_x_xz_mpi
 write(11)g5_x_xv_mpi
 write(11)g5_x_yy_mpi
 write(11)g5_x_yz_mpi
 write(11)g5_x_yv_mpi
 write(11)g5_x_zz_mpi
 write(11)g5_x_zv_mpi
 write(11)g5_x_vv_mpi

!g2_y-... & g5_y_...
 write(11)g2_y_x_mpi
 write(11)g2_y_y_mpi
 write(11)g2_y_z_mpi
 write(11)g2_y_v_mpi
 write(11)g5_y_xx_mpi
 write(11)g5_y_xy_mpi
 write(11)g5_y_xz_mpi
 write(11)g5_y_xv_mpi
 write(11)g5_y_yy_mpi
 write(11)g5_y_yz_mpi
 write(11)g5_y_yv_mpi
 write(11)g5_y_zz_mpi
 write(11)g5_y_zv_mpi
 write(11)g5_y_vv_mpi

!g2_z-... & g5_z_...
 write(11)g2_z_x_mpi
 write(11)g2_z_y_mpi
 write(11)g2_z_z_mpi
 write(11)g2_z_v_mpi
 write(11)g5_z_xx_mpi
 write(11)g5_z_xy_mpi
 write(11)g5_z_xz_mpi
 write(11)g5_z_xv_mpi
 write(11)g5_z_yy_mpi
 write(11)g5_z_yz_mpi
 write(11)g5_z_yv_mpi
 write(11)g5_z_zz_mpi
 write(11)g5_z_zv_mpi
 write(11)g5_z_vv_mpi

!g2_v-... & g5_v_...
 write(11)g2_v_x_mpi
 write(11)g2_v_y_mpi
 write(11)g2_v_z_mpi
 write(11)g2_v_v_mpi
 write(11)g5_v_xx_mpi
 write(11)g5_v_xy_mpi
 write(11)g5_v_xz_mpi
 write(11)g5_v_xv_mpi
 write(11)g5_v_yy_mpi
 write(11)g5_v_yz_mpi
 write(11)g5_v_yv_mpi
 write(11)g5_v_zz_mpi
 write(11)g5_v_zv_mpi
 write(11)g5_v_vv_mpi

!--------------
!Write sf_deriv
!--------------

!g2_deriv_x-... & g5_deriv_x-...
 write(12)g2_deriv_x_x_mpi
 write(12)g2_deriv_x_y_mpi
 write(12)g2_deriv_x_z_mpi
 write(12)g2_deriv_x_v_mpi
 write(12)g5_deriv_x_xx_mpi
 write(12)g5_deriv_x_xy_mpi
 write(12)g5_deriv_x_xz_mpi
 write(12)g5_deriv_x_xv_mpi
 write(12)g5_deriv_x_yy_mpi
 write(12)g5_deriv_x_yz_mpi
 write(12)g5_deriv_x_yv_mpi
 write(12)g5_deriv_x_zz_mpi
 write(12)g5_deriv_x_zv_mpi
 write(12)g5_deriv_x_vv_mpi

!g2_deriv_y-... & g5_deriv_y-...
 write(12)g2_deriv_y_x_mpi
 write(12)g2_deriv_y_y_mpi
 write(12)g2_deriv_y_z_mpi
 write(12)g2_deriv_y_v_mpi
 write(12)g5_deriv_y_xx_mpi
 write(12)g5_deriv_y_xy_mpi
 write(12)g5_deriv_y_xz_mpi
 write(12)g5_deriv_y_xv_mpi
 write(12)g5_deriv_y_yy_mpi
 write(12)g5_deriv_y_yz_mpi
 write(12)g5_deriv_y_yv_mpi
 write(12)g5_deriv_y_zz_mpi
 write(12)g5_deriv_y_zv_mpi
 write(12)g5_deriv_y_vv_mpi

!g2_deriv_z-... & g5_deriv_z-...
 write(12)g2_deriv_z_x_mpi
 write(12)g2_deriv_z_y_mpi
 write(12)g2_deriv_z_z_mpi
 write(12)g2_deriv_z_v_mpi
 write(12)g5_deriv_z_xx_mpi
 write(12)g5_deriv_z_xy_mpi
 write(12)g5_deriv_z_xz_mpi
 write(12)g5_deriv_z_xv_mpi
 write(12)g5_deriv_z_yy_mpi
 write(12)g5_deriv_z_yz_mpi
 write(12)g5_deriv_z_yv_mpi
 write(12)g5_deriv_z_zz_mpi
 write(12)g5_deriv_z_zv_mpi
 write(12)g5_deriv_z_vv_mpi

!g2_deriv_v-... & g5_deriv_v-...
 write(12)g2_deriv_v_x_mpi
 write(12)g2_deriv_v_y_mpi
 write(12)g2_deriv_v_z_mpi
 write(12)g2_deriv_v_v_mpi
 write(12)g5_deriv_v_xx_mpi
 write(12)g5_deriv_v_xy_mpi
 write(12)g5_deriv_v_xz_mpi
 write(12)g5_deriv_v_xv_mpi
 write(12)g5_deriv_v_yy_mpi
 write(12)g5_deriv_v_yz_mpi
 write(12)g5_deriv_v_yv_mpi
 write(12)g5_deriv_v_zz_mpi
 write(12)g5_deriv_v_zv_mpi
 write(12)g5_deriv_v_vv_mpi

 endif ! myrank == 0


!do i = 1,12
!print'(A5,3F20.15)',"posi ",posi1(i+natom_x+natom_y+natom_z,1),posi1(i+natom_x+natom_y+natom_z,2),posi1(i+natom_x+natom_y+natom_z,3)
!do j = 1,132
!do k = 1,22
!print'(A4,2F20.15,F10.5,2I3)',"g5d ",g2_x_x_mpi(i,k),g2_deriv_x_x_mpi(i,j,1,k),eta2_x_x(k),j,k
!print'(A4,2F20.15,3F10.5,2I3)',"g5d ",g5_x_xx_mpi(i,k),g5_deriv_x_xx_mpi(i,j,1,k),eta5_x_xx(k),zeta5_x_xx(k),lambda5_x_xx(k),j,k
!print'(A4,2F20.15,3F10.5,2I3)',"g5d ",g5_x_yy_mpi(i,k),g5_deriv_x_yy_mpi(i,j,1,k),eta5_x_yy(k),zeta5_x_yy(k),lambda5_x_yy(k),j,k
!print'(A4,2F20.15,3F10.5,2I3)',"g5d ",g5_z_xx_mpi(i,k),g5_deriv_z_xx_mpi(i,j,1,k),eta5_z_xx(k),zeta5_z_xx(k),lambda5_z_xx(k),j,k
!print'(A4,2F20.15,3F10.5,2I3)',"g5d ",g5_v_xx_mpi(i,k),g5_deriv_v_xx_mpi(i,j,1,k),eta5_v_xx(k),zeta5_v_xx(k),lambda5_v_xx(k),j,k
!print'(A4,2F20.15,3F10.5,2I3)',"g5d ",g5_x_zv_mpi(i,k),g5_deriv_x_zv_mpi(i,j,1,k),eta5_x_zv(k),zeta5_x_zv(k),lambda5_x_zv(k),j,k
!print'(A4,2F20.15,3F10.5,2I3)',"g5d ",g5_x_vv_mpi(i,k),g5_deriv_x_vv_mpi(i,j,1,k),eta5_x_vv(k),zeta5_x_vv(k),lambda5_x_vv(k),j,k
!enddo
!enddo
!enddo


end subroutine sf_4
