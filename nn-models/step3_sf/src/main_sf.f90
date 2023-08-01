program main_sf
use allocarray
implicit none
include 'mpif.h'
integer ierr, num_mpi, myrank
integer i, j, k
integer stat, counter, nfile, imd, md_steps, neighbor

integer nelem, natom
integer natom_a, natom_b, natom_c, natom_d
integer natom_1, natom_2, natom_3, natom_4
integer io_a, io_b, io_c, io_d
integer io_sf
character(len=4) elem_a, elem_b, elem_c, elem_d
character(len=4) elem_1, elem_2, elem_3, elem_4

double precision lattice, x(3), y(3), z(3), xlattice, ylattice, zlattice
double precision rc, rn
double precision,allocatable :: posi(:,:)
character dir*6, system*10, phi*10
character(len=120) tag, filename, dir_data, dir_sf

double precision  t0, t1
integer           dummyi
double precision  dummy
character(len=20) dummyc

!Variables for symmetry function
!A
integer &
  num_g2_a_a,  num_g2_a_b,  num_g2_a_c,  num_g2_a_d, &
  num_g5_a_aa, num_g5_a_ab, num_g5_a_ac, num_g5_a_ad, &
  num_g5_a_bb, num_g5_a_bc, num_g5_a_bd, &
  num_g5_a_cc, num_g5_a_cd, &
  num_g5_a_dd
!B
integer &
  num_g2_b_a,  num_g2_b_b,  num_g2_b_c,  num_g2_b_d, &
  num_g5_b_aa, num_g5_b_ab, num_g5_b_ac, num_g5_b_ad, &
  num_g5_b_bb, num_g5_b_bc, num_g5_b_bd, &
  num_g5_b_cc, num_g5_b_cd, &
  num_g5_b_dd
!C
integer &
  num_g2_c_a,  num_g2_c_b,  num_g2_c_c,  num_g2_c_d, &
  num_g5_c_aa, num_g5_c_ab, num_g5_c_ac, num_g5_c_ad, &
  num_g5_c_bb, num_g5_c_bc, num_g5_c_bd, &
  num_g5_c_cc, num_g5_c_cd, &
  num_g5_c_dd
!D
integer &
  num_g2_d_a,  num_g2_d_b,  num_g2_d_c,  num_g2_d_d, &
  num_g5_d_aa, num_g5_d_ab, num_g5_d_ac, num_g5_d_ad, &
  num_g5_d_bb, num_g5_d_bc, num_g5_d_bd, &
  num_g5_d_cc, num_g5_d_cd, &
  num_g5_d_dd
!A
double precision,allocatable,dimension(:,:) :: &
  g2_a_a,  g2_a_b,  g2_a_c,  g2_a_d, &
  g5_a_aa, g5_a_ab, g5_a_ac, g5_a_ad, &
  g5_a_bb, g5_a_bc, g5_a_bd, &
  g5_a_cc, g5_a_cd, &
  g5_a_dd
double precision,allocatable,dimension(:) :: &
  eta2_a_a,  rs2_a_a, &
  eta2_a_b,  rs2_a_b, &
  eta2_a_c,  rs2_a_c, &
  eta2_a_d,  rs2_a_d, &
  eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
  eta5_a_ab, theta5_a_ab, zeta5_a_ab, lambda5_a_ab, &
  eta5_a_ac, theta5_a_ac, zeta5_a_ac, lambda5_a_ac, &
  eta5_a_ad, theta5_a_ad, zeta5_a_ad, lambda5_a_ad, &
  eta5_a_bb, theta5_a_bb, zeta5_a_bb, lambda5_a_bb, &
  eta5_a_bc, theta5_a_bc, zeta5_a_bc, lambda5_a_bc, &
  eta5_a_bd, theta5_a_bd, zeta5_a_bd, lambda5_a_bd, &
  eta5_a_cc, theta5_a_cc, zeta5_a_cc, lambda5_a_cc, &
  eta5_a_cd, theta5_a_cd, zeta5_a_cd, lambda5_a_cd, &
  eta5_a_dd, theta5_a_dd, zeta5_a_dd, lambda5_a_dd
!B
double precision,allocatable,dimension(:,:) :: &
  g2_b_a,  g2_b_b,  g2_b_c,  g2_b_d, &
  g5_b_aa, g5_b_ab, g5_b_ac, g5_b_ad, &
  g5_b_bb, g5_b_bc, g5_b_bd, &
  g5_b_cc, g5_b_cd, &
  g5_b_dd
double precision,allocatable,dimension(:) :: &
  eta2_b_a,  rs2_b_a, &
  eta2_b_b,  rs2_b_b, &
  eta2_b_c,  rs2_b_c, &
  eta2_b_d,  rs2_b_d, &
  eta5_b_aa, theta5_b_aa, zeta5_b_aa, lambda5_b_aa, &
  eta5_b_ab, theta5_b_ab, zeta5_b_ab, lambda5_b_ab, &
  eta5_b_ac, theta5_b_ac, zeta5_b_ac, lambda5_b_ac, &
  eta5_b_ad, theta5_b_ad, zeta5_b_ad, lambda5_b_ad, &
  eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
  eta5_b_bc, theta5_b_bc, zeta5_b_bc, lambda5_b_bc, &
  eta5_b_bd, theta5_b_bd, zeta5_b_bd, lambda5_b_bd, &
  eta5_b_cc, theta5_b_cc, zeta5_b_cc, lambda5_b_cc, &
  eta5_b_cd, theta5_b_cd, zeta5_b_cd, lambda5_b_cd, &
  eta5_b_dd, theta5_b_dd, zeta5_b_dd, lambda5_b_dd
!C
double precision,allocatable,dimension(:,:) :: &
  g2_c_a,  g2_c_b,  g2_c_c,  g2_c_d, &
  g5_c_aa, g5_c_ab, g5_c_ac, g5_c_ad, &
  g5_c_bb, g5_c_bc, g5_c_bd, &
  g5_c_cc, g5_c_cd, &
  g5_c_dd
double precision,allocatable,dimension(:) :: &
  eta2_c_a,  rs2_c_a, &
  eta2_c_b,  rs2_c_b, &
  eta2_c_c,  rs2_c_c, &
  eta2_c_d,  rs2_c_d, &
  eta5_c_aa, theta5_c_aa, zeta5_c_aa, lambda5_c_aa, &
  eta5_c_ab, theta5_c_ab, zeta5_c_ab, lambda5_c_ab, &
  eta5_c_ac, theta5_c_ac, zeta5_c_ac, lambda5_c_ac, &
  eta5_c_ad, theta5_c_ad, zeta5_c_ad, lambda5_c_ad, &
  eta5_c_bb, theta5_c_bb, zeta5_c_bb, lambda5_c_bb, &
  eta5_c_bc, theta5_c_bc, zeta5_c_bc, lambda5_c_bc, &
  eta5_c_bd, theta5_c_bd, zeta5_c_bd, lambda5_c_bd, &
  eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc, &
  eta5_c_cd, theta5_c_cd, zeta5_c_cd, lambda5_c_cd, &
  eta5_c_dd, theta5_c_dd, zeta5_c_dd, lambda5_c_dd
!D
double precision,allocatable,dimension(:,:) :: &
  g2_d_a,  g2_d_b,  g2_d_c,  g2_d_d, &
  g5_d_aa, g5_d_ab, g5_d_ac, g5_d_ad, &
  g5_d_bb, g5_d_bc, g5_d_bd, &
  g5_d_cc, g5_d_cd, &
  g5_d_dd
double precision,allocatable,dimension(:) :: &
  eta2_d_a, rs2_d_a, &
  eta2_d_b, rs2_d_b, &
  eta2_d_c, rs2_d_c, &
  eta2_d_d, rs2_d_d, &
  eta5_d_aa, theta5_d_aa, zeta5_d_aa, lambda5_d_aa, &
  eta5_d_ab, theta5_d_ab, zeta5_d_ab, lambda5_d_ab, &
  eta5_d_ac, theta5_d_ac, zeta5_d_ac, lambda5_d_ac, &
  eta5_d_ad, theta5_d_ad, zeta5_d_ad, lambda5_d_ad, &
  eta5_d_bb, theta5_d_bb, zeta5_d_bb, lambda5_d_bb, &
  eta5_d_bc, theta5_d_bc, zeta5_d_bc, lambda5_d_bc, &
  eta5_d_bd, theta5_d_bd, zeta5_d_bd, lambda5_d_bd, &
  eta5_d_cc, theta5_d_cc, zeta5_d_cc, lambda5_d_cc, &
  eta5_d_cd, theta5_d_cd, zeta5_d_cd, lambda5_d_cd, &
  eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd
!A
double precision,allocatable,dimension(:,:,:,:) :: &
  g2_deriv_a_a,  g2_deriv_a_b,  g2_deriv_a_c,  g2_deriv_a_d, &
  g5_deriv_a_aa, g5_deriv_a_ab, g5_deriv_a_ac, g5_deriv_a_ad, &
  g5_deriv_a_bb, g5_deriv_a_bc, g5_deriv_a_bd, &
  g5_deriv_a_cc, g5_deriv_a_cd, &
  g5_deriv_a_dd
!B
double precision,allocatable,dimension(:,:,:,:) :: &
  g2_deriv_b_a,  g2_deriv_b_b,  g2_deriv_b_c,  g2_deriv_b_d, &
  g5_deriv_b_aa, g5_deriv_b_ab, g5_deriv_b_ac, g5_deriv_b_ad, &
  g5_deriv_b_bb, g5_deriv_b_bc, g5_deriv_b_bd, &
  g5_deriv_b_cc, g5_deriv_b_cd, &
  g5_deriv_b_dd
!C
double precision,allocatable,dimension(:,:,:,:) :: &
  g2_deriv_c_a,  g2_deriv_c_b,  g2_deriv_c_c,  g2_deriv_c_d, &
  g5_deriv_c_aa, g5_deriv_c_ab, g5_deriv_c_ac, g5_deriv_c_ad, &
  g5_deriv_c_bb, g5_deriv_c_bc, g5_deriv_c_bd, &
  g5_deriv_c_cc, g5_deriv_c_cd, &
  g5_deriv_c_dd
!D
double precision,allocatable,dimension(:,:,:,:) :: &
  g2_deriv_d_a,  g2_deriv_d_b,  g2_deriv_d_c,  g2_deriv_d_d, &
  g5_deriv_d_aa, g5_deriv_d_ab, g5_deriv_d_ac, g5_deriv_d_ad, &
  g5_deriv_d_bb, g5_deriv_d_bc, g5_deriv_d_bd, &
  g5_deriv_d_cc, g5_deriv_d_cd, &
  g5_deriv_d_dd


 !MPI INIT
 call mpi_init( ierr )
 call mpi_comm_size( mpi_comm_world, num_mpi, ierr )
 call mpi_comm_rank( mpi_comm_world, myrank, ierr )
 write(*,*) " NUM_MPI: ", myrank, num_mpi
 call mpi_barrier( mpi_comm_world, ierr )

 if( myrank == 0 ) call cpu_time( t0 )


!*************************************************************
! Read input_nnp.dat file
!*************************************************************

 if( myrank == 0 )then
   call read_input_nnp4( &
        elem_a, &
        num_g2_a_a,  num_g2_a_b,  num_g2_a_c,  num_g2_a_d, &
        num_g5_a_aa, num_g5_a_ab, num_g5_a_ac, num_g5_a_ad, &
        num_g5_a_bb, num_g5_a_bc, num_g5_a_bd, &
        num_g5_a_cc, num_g5_a_cd, &
        num_g5_a_dd, &
        elem_b, &
        num_g2_b_a,  num_g2_b_b,  num_g2_b_c,  num_g2_b_d, &
        num_g5_b_aa, num_g5_b_ab, num_g5_b_ac, num_g5_b_ad, &
        num_g5_b_bb, num_g5_b_bc, num_g5_b_bd, &
        num_g5_b_cc, num_g5_b_cd, &
        num_g5_b_dd, &
        elem_c, &
        num_g2_c_a,  num_g2_c_b,  num_g2_c_c,  num_g2_c_d, &
        num_g5_c_aa, num_g5_c_ab, num_g5_c_ac, num_g5_c_ad, &
        num_g5_c_bb, num_g5_c_bc, num_g5_c_bd, &
        num_g5_c_cc, num_g5_c_cd, &
        num_g5_c_dd, &
        elem_d, &
        num_g2_d_a,  num_g2_d_b,  num_g2_d_c,  num_g2_d_d, &
        num_g5_d_aa, num_g5_d_ab, num_g5_d_ac, num_g5_d_ad, &
        num_g5_d_bb, num_g5_d_bc, num_g5_d_bd, &
        num_g5_d_cc, num_g5_d_cd, &
        num_g5_d_dd, &
        rc, rn, phi, neighbor, &
        dir_data, dir_sf )
 endif ! myrank == 0

 call mpi_barrier( mpi_comm_world, ierr )
 call mpi_bcast( elem_a, 4, mpi_character, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_a_a, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_a_b, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_a_c, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_a_d, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_a_aa, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_a_ab, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_a_ac, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_a_ad, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_a_bb, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_a_bc, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_a_bd, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_a_cc, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_a_cd, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_a_dd, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( elem_b, 4, mpi_character, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_b_a, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_b_b, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_b_c, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_b_d, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_b_aa, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_b_ab, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_b_ac, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_b_ad, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_b_bb, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_b_bc, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_b_bd, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_b_cc, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_b_cd, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_b_dd, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( elem_c, 4, mpi_character, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_c_a, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_c_b, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_c_c, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_c_d, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_c_aa, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_c_ab, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_c_ac, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_c_ad, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_c_bb, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_c_bc, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_c_bd, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_c_cc, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_c_cd, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_c_dd, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( elem_d, 4, mpi_character, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_d_a, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_d_b, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_d_c, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_d_d, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_d_aa, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_d_ab, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_d_ac, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_d_ad, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_d_bb, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_d_bc, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_d_bd, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_d_cc, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_d_cd, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_d_dd, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( rc, 1, mpi_double_precision, 0, mpi_comm_world, ierr )
 call mpi_bcast( rn, 1, mpi_double_precision, 0, mpi_comm_world, ierr )
 call mpi_bcast( neighbor, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( dir_data, 120, mpi_character, 0, mpi_comm_world, ierr )
 call mpi_bcast( dir_sf, 120, mpi_character, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )


!*************************************************************
! Read input_tag.dat file
!*************************************************************

 if( myrank == 0 )then
   open(unit=99,file="input_tag.dat",action="read")
   ! input_tag.dat includes
   !----------------------
   ! bulk_au_n32_300k
   ! au111_n32_300k
   ! .
   ! .
   ! .
   !---

   !Number of lines in input_tag.dat
   counter = 0
   rewind(99)
   do
     read( 99, '(A5)', iostat=stat ) dummyc
     if( stat /= 0 ) exit
     if( trim( dummyc ) /= '' ) counter = counter + 1
   enddo
   write(*,*) " "
   write(*,*)  counter, "tags are recognized"
   write(*,*) " "
 endif ! myrank == 0

 call mpi_barrier( mpi_comm_world, ierr )
 call mpi_bcast( counter, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )

 rewind(99)
 !*************************************************************
 ! DO LOOP : input_tag
 !*************************************************************
 
 do nfile = 1, counter

   if( myrank == 0 )then
     write(*,*)  " "
     write(*,*)  " # ", nfile
     write(*,*) " ------------------------------------------------------------"

     read(99,*) tag
     write(*,*)  " "
     write(*,*)  " Tag_name: ", trim( tag )
     write(*,*)  " "
   endif ! myrank == 0

   call mpi_barrier( mpi_comm_world, ierr )
   call mpi_bcast( tag, 120, mpi_character, 0, mpi_comm_world, ierr )
   call mpi_barrier( mpi_comm_world, ierr )

1111 format( A120 )

   !Info READ
   write( filename, 1111 ) tag
   !filename = '../step2_data_binary/data/info_'//trim(adjustl(filename))//'.dat'
   filename = trim(adjustl(dir_data))//'/info_'//trim(adjustl(filename))//'.dat'
   open(unit=1,file=filename,action="read",form="unformatted")
   !positions READ
   write( filename, 1111 ) tag
   filename = trim(adjustl(dir_data))//'/positions_'//trim(adjustl(filename))//'.dat'
   open(unit=2,file=filename,action="read",form="unformatted")

   !sf WRITE
   write( filename, 1111 ) tag
   !write( dummyc, * ) int(rc)
   !filename = './data_sf/rc_'//trim(adjustl(dummyc))//'/sf_'//trim(adjustl(filename))//'.dat'
   filename = trim(adjustl(dir_sf))//'/sf_'//trim(adjustl(filename))//'.dat'
   open(unit=11,file=filename,action="write",form="unformatted")
   !sf_deriv WRITE
   write( filename, 1111 ) tag
   !write( dummyc, * ) int(rc)
   !filename = './data_sf/rc_'//trim(adjustl(dummyc))//'/sf_deriv_'//trim(adjustl(filename))//'.dat'
   filename = trim(adjustl(dir_sf))//'/sf_deriv_'//trim(adjustl(filename))//'.dat'
   open(unit=12,file=filename,action="write",form="unformatted")


   !*************************************************************
   ! Read info_tag.dat file (binary)
   !*************************************************************
   ! info_tag.dat includes
   !------------------------   -----------------------
   ! Total_number_of_ions 32   Total_number_of_ions 32
   ! Ions_per_type 1           Ions_per_type 3
   ! Au 32                     Li 12
   !                           P 4
   !                           O 16
   ! MD_steps 400              MD_steps 2000
   !------------------------   -----------------------

   ! 4 elements
   if( myrank == 0 )then
     call read_io4( natom, nelem, &
                    elem_a,  elem_b,  elem_c,  elem_d, &
                    natom_a, natom_b, natom_c, natom_d, &
                    io_a,    io_b,    io_c,    io_d )
   endif ! myrank == 0
 
   call mpi_barrier( mpi_comm_world, ierr )
   call mpi_bcast( natom, 1, mpi_integer, 0, mpi_comm_world, ierr )
   call mpi_bcast( nelem, 1, mpi_integer, 0, mpi_comm_world, ierr )
   call mpi_bcast( natom_a, 1, mpi_integer, 0, mpi_comm_world, ierr )
   call mpi_bcast( natom_b, 1, mpi_integer, 0, mpi_comm_world, ierr )
   call mpi_bcast( natom_c, 1, mpi_integer, 0, mpi_comm_world, ierr )
   call mpi_bcast( natom_d, 1, mpi_integer, 0, mpi_comm_world, ierr )
   call mpi_bcast( io_a, 1, mpi_integer, 0, mpi_comm_world, ierr )
   call mpi_bcast( io_b, 1, mpi_integer, 0, mpi_comm_world, ierr )
   call mpi_bcast( io_c, 1, mpi_integer, 0, mpi_comm_world, ierr )
   call mpi_bcast( io_d, 1, mpi_integer, 0, mpi_comm_world, ierr )
   call mpi_barrier( mpi_comm_world, ierr )

   !MD steps
   if( myrank == 0 )then
     read(1) md_steps
     write(*,*) " "
     write(*,*) " MD steps: ", md_steps
     write(*,*) " "
   endif ! myrank == 0
   call mpi_barrier( mpi_comm_world, ierr )
   call mpi_bcast( md_steps, 1, mpi_integer, 0, mpi_comm_world, ierr )
   call mpi_barrier( mpi_comm_world, ierr )

   !*************************************************************
   ! Allocate G2, G5 & read paramters
   !*************************************************************
   ! 1 element
   if( nelem == 1 )then
     !A
     if( io_a == 1 .and. io_b == 0 .and. io_c == 0 .and. io_d == 0 )then
       if( myrank == 0 ) write(*,*) " Allocate G2(a-a),G5(a-aa): ", elem_a
       call allocate_1( "g2_a_a", g2_a_a, eta2_a_a, rs2_a_a, &
                        "g5_a_aa", g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                        num_g2_a_a, num_g5_a_aa, &
                        natom, natom_a, &
                        g2_deriv_a_a, g5_deriv_a_aa )
     !B
     elseif( io_a == 0 .and. io_b == 1 .and. io_c == 0 .and. io_d == 0 )then
       if( myrank == 0 ) write(*,*) " Allocate G2(b-b),G5(b-bb): ", elem_b
       call allocate_1( "g2_b_b", g2_b_b, eta2_b_b, rs2_b_b, &
                        "g5_b_bb", g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                        num_g2_b_b, num_g5_b_bb, &
                        natom, natom_b, &
                        g2_deriv_b_b, g5_deriv_b_bb )
     !C
     elseif( io_a == 0 .and. io_b == 0 .and. io_c == 1 .and. io_d == 0 )then
       if( myrank == 0 ) write(*,*) " Allocate G2(c-c),G5(c-cc): ", elem_c
       call allocate_1( "g2_c_c", g2_c_c, eta2_c_c, rs2_c_c, &
                        "g5_c_cc", g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc, &
                        num_g2_c_c, num_g5_c_cc, &
                        natom, natom_c, &
                        g2_deriv_c_c, g5_deriv_c_cc )
     !D
     elseif( io_a == 0 .and. io_b == 0 .and. io_c == 0 .and. io_d == 1 )then
       if( myrank == 0 ) write(*,*) " Allocate G2(d-d),G5(d-dd): ", elem_d
       call allocate_1( "g2_d_d", g2_d_d, eta2_d_d, rs2_d_d, &
                        "g5_d_dd", g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd, &
                        num_g2_d_d, num_g5_d_dd, &
                        natom, natom_d, &
                        g2_deriv_d_d, g5_deriv_d_dd )
     else
       write(*,*) "Error(main): allocate nelem = 1"
       stop
     endif
   endif ! nelem == 1
   ! 2 elements
   if( nelem == 2 )then
     !A,B
     if( io_a == 1 .and. io_b == 1 .and. io_c == 0 .and. io_d == 0 )then
       if( myrank == 0 ) write(*,*) " Allocate G2(a-a,a-b),G5(a-aa,a-ab,a-bb): ", elem_a
       call allocate_2( "g2_a_a",  g2_a_a,  eta2_a_a,  rs2_a_a, &
                        "g2_a_b",  g2_a_b,  eta2_a_b,  rs2_a_b, &
                        "g5_a_aa", g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                        "g5_a_ab", g5_a_ab, eta5_a_ab, theta5_a_ab, zeta5_a_ab, lambda5_a_ab, &
                        "g5_a_bb", g5_a_bb, eta5_a_bb, theta5_a_bb, zeta5_a_bb, lambda5_a_bb, &
                        num_g2_a_a, num_g2_a_b, num_g5_a_aa, num_g5_a_ab, num_g5_a_bb, &
                        natom, natom_a, &
                        g2_deriv_a_a, g2_deriv_a_b, g5_deriv_a_aa, g5_deriv_a_ab, g5_deriv_a_bb )
       if( myrank == 0 ) write(*,*) " Allocate G2(b-a,b-b),G5(b-aa,b-ab,b-bb): ", elem_b
       call allocate_2( "g2_b_a",  g2_b_a,  eta2_b_a,  rs2_b_a, &
                        "g2_b_b",  g2_b_b,  eta2_b_b,  rs2_b_b, &
                        "g5_b_aa", g5_b_aa, eta5_b_aa, theta5_b_aa, zeta5_b_aa, lambda5_b_aa, &
                        "g5_b_ab", g5_b_ab, eta5_b_ab, theta5_b_ab, zeta5_b_ab, lambda5_b_ab, &
                        "g5_b_bb", g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                        num_g2_b_a, num_g2_b_b, num_g5_b_aa, num_g5_b_ab, num_g5_b_bb, &
                        natom, natom_b, &
                        g2_deriv_b_a, g2_deriv_b_b, g5_deriv_b_aa, g5_deriv_b_ab, g5_deriv_b_bb )
     !A,C
     elseif( io_a == 1 .and. io_b == 0 .and. io_c == 1 .and. io_d == 0 )then
       if( myrank == 0 ) write(*,*) " Allocate G2(a-a,a-c),G5(a-aa,a-ac,a-cc): ", elem_a
       call allocate_2( "g2_a_a",  g2_a_a,  eta2_a_a,  rs2_a_a, &
                        "g2_a_c",  g2_a_c,  eta2_a_c,  rs2_a_c, &
                        "g5_a_aa", g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                        "g5_a_ac", g5_a_ac, eta5_a_ac, theta5_a_ac, zeta5_a_ac, lambda5_a_ac, &
                        "g5_a_cc", g5_a_cc, eta5_a_cc, theta5_a_cc, zeta5_a_cc, lambda5_a_cc, &
                        num_g2_a_a, num_g2_a_c, num_g5_a_aa, num_g5_a_ac, num_g5_a_cc, &
                        natom, natom_a, &
                        g2_deriv_a_a, g2_deriv_a_c, g5_deriv_a_aa, g5_deriv_a_ac, g5_deriv_a_cc )
       if( myrank == 0 ) write(*,*) " Allocate G2(c-a,c-c),G5(c-aa,c-ac,c-cc): ", elem_c
       call allocate_2( "g2_c_a",  g2_c_a,  eta2_c_a,  rs2_c_a, &
                        "g2_c_c",  g2_c_c,  eta2_c_c,  rs2_c_c, &
                        "g5_c_aa", g5_c_aa, eta5_c_aa, theta5_c_aa, zeta5_c_aa, lambda5_c_aa, &
                        "g5_c_ac", g5_c_ac, eta5_c_ac, theta5_c_ac, zeta5_c_ac, lambda5_c_ac, &
                        "g5_c_cc", g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc, &
                        num_g2_c_a, num_g2_c_c, num_g5_c_aa, num_g5_c_ac, num_g5_c_cc, &
                        natom, natom_c, &
                        g2_deriv_c_a, g2_deriv_c_c, g5_deriv_c_aa, g5_deriv_c_ac, g5_deriv_c_cc )
     !A,D
     elseif( io_a == 1 .and. io_b == 0 .and. io_c == 0 .and. io_d == 1 )then
       if( myrank == 0 ) write(*,*) " Allocate G2(a-a,a-d),G5(a-aa,a-ad,a-dd): ", elem_a
       call allocate_2( "g2_a_a",  g2_a_a,  eta2_a_a,  rs2_a_a, &
                        "g2_a_d",  g2_a_d,  eta2_a_d,  rs2_a_d, &
                        "g5_a_aa", g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                        "g5_a_ad", g5_a_ad, eta5_a_ad, theta5_a_ad, zeta5_a_ad, lambda5_a_ad, &
                        "g5_a_dd", g5_a_dd, eta5_a_dd, theta5_a_dd, zeta5_a_dd, lambda5_a_dd, &
                        num_g2_a_a, num_g2_a_d, num_g5_a_aa, num_g5_a_ad, num_g5_a_dd, &
                        natom, natom_a, &
                        g2_deriv_a_a, g2_deriv_a_d, g5_deriv_a_aa, g5_deriv_a_ad, g5_deriv_a_dd )
       if( myrank == 0 ) write(*,*) " Allocate G2(d-a,d-d),G5(d-aa,d-ad,d-dd): ", elem_d
       call allocate_2( "g2_d_a",  g2_d_a,  eta2_d_a,  rs2_d_a, &
                        "g2_d_d",  g2_d_d,  eta2_d_d,  rs2_d_d, &
                        "g5_d_aa", g5_d_aa, eta5_d_aa, theta5_d_aa, zeta5_d_aa, lambda5_d_aa, &
                        "g5_d_ad", g5_d_ad, eta5_d_ad, theta5_d_ad, zeta5_d_ad, lambda5_d_ad, &
                        "g5_d_dd", g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd, &
                        num_g2_d_a, num_g2_d_d, num_g5_d_aa, num_g5_d_ad, num_g5_d_dd, &
                        natom, natom_d, &
                        g2_deriv_d_a, g2_deriv_d_d, g5_deriv_d_aa, g5_deriv_d_ad, g5_deriv_d_dd )
     !B,C
     elseif( io_a == 0 .and. io_b == 1 .and. io_c == 1 .and. io_d == 0 )then
       if( myrank == 0 ) write(*,*) " Allocate G2(b-b,b-c),G5(b-bb,b-bc,b-cc): ", elem_b
       call allocate_2( "g2_b_b",  g2_b_b,  eta2_b_b,  rs2_b_b, &
                        "g2_b_c",  g2_b_c,  eta2_b_c,  rs2_b_c, &
                        "g5_b_bb", g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                        "g5_b_bc", g5_b_bc, eta5_b_bc, theta5_b_bc, zeta5_b_bc, lambda5_b_bc, &
                        "g5_b_cc", g5_b_cc, eta5_b_cc, theta5_b_cc, zeta5_b_cc, lambda5_b_cc, &
                        num_g2_b_b, num_g2_b_c, num_g5_b_bb, num_g5_b_bc, num_g5_b_cc, &
                        natom, natom_b, &
                        g2_deriv_b_b, g2_deriv_b_c, g5_deriv_b_bb, g5_deriv_b_bc, g5_deriv_b_cc )
       if( myrank == 0 ) write(*,*) " Allocate G2(c-b,c-c),G5(c-bb,c-bc,c-cc): ", elem_c
       call allocate_2( "g2_c_b",  g2_c_b,  eta2_c_b,  rs2_c_b, &
                        "g2_c_c",  g2_c_c,  eta2_c_c,  rs2_c_c, &
                        "g5_c_bb", g5_c_bb, eta5_c_bb, theta5_c_bb, zeta5_c_bb, lambda5_c_bb, &
                        "g5_c_bc", g5_c_bc, eta5_c_bc, theta5_c_bc, zeta5_c_bc, lambda5_c_bc, &
                        "g5_c_cc", g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc, &
                        num_g2_c_b, num_g2_c_c, num_g5_c_bb, num_g5_c_bc, num_g5_c_cc, &
                        natom, natom_c, &
                        g2_deriv_c_b, g2_deriv_c_c, g5_deriv_c_bb, g5_deriv_c_bc, g5_deriv_c_cc )
     !B,D
     elseif( io_a == 0 .and. io_b == 1 .and. io_c == 0 .and. io_d == 1 )then
       if( myrank == 0 ) write(*,*) " Allocate G2(b-b,b-d),G5(b-bb,b-bd,b-dd): ", elem_b
       call allocate_2( "g2_b_b",  g2_b_b,  eta2_b_b,  rs2_b_b, &
                        "g2_b_d",  g2_b_d,  eta2_b_d,  rs2_b_d, &
                        "g5_b_bb", g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                        "g5_b_bd", g5_b_bd, eta5_b_bd, theta5_b_bd, zeta5_b_bd, lambda5_b_bd, &
                        "g5_b_dd", g5_b_dd, eta5_b_dd, theta5_b_dd, zeta5_b_dd, lambda5_b_dd, &
                        num_g2_b_b, num_g2_b_d, num_g5_b_bb, num_g5_b_bd, num_g5_b_dd, &
                        natom, natom_b, &
                        g2_deriv_b_b, g2_deriv_b_d, g5_deriv_b_bb, g5_deriv_b_bd, g5_deriv_b_dd )
       if( myrank == 0 ) write(*,*) " Allocate G2(d-b,d-d),G5(d-bb,d-bd,d-dd): ", elem_d
       call allocate_2( "g2_d_b",  g2_d_b,  eta2_d_b,  rs2_d_b, &
                        "g2_d_d",  g2_d_d,  eta2_d_d,  rs2_d_d, &
                        "g5_d_bb", g5_d_bb, eta5_d_bb, theta5_d_bb, zeta5_d_bb, lambda5_d_bb, &
                        "g5_d_bd", g5_d_bd, eta5_d_bd, theta5_d_bd, zeta5_d_bd, lambda5_d_bd, &
                        "g5_d_dd", g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd, &
                        num_g2_d_b, num_g2_d_d, num_g5_d_bb, num_g5_d_bd, num_g5_d_dd, &
                        natom, natom_d, &
                        g2_deriv_d_b, g2_deriv_d_d, g5_deriv_d_bb, g5_deriv_d_bd, g5_deriv_d_dd )
     !C,D
     elseif( io_a == 0 .and. io_b == 0 .and. io_c == 1 .and. io_d == 1 )then
       if( myrank == 0 ) write(*,*) " Allocate G2(c-c,c-d),G5(c-cc,c-cd,c-dd): ", elem_c
       call allocate_2( "g2_c_c",  g2_c_c,  eta2_c_c,  rs2_c_c, &
                        "g2_c_d",  g2_c_d,  eta2_c_d,  rs2_c_d, &
                        "g5_c_cc", g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc, &
                        "g5_c_cd", g5_c_cd, eta5_c_cd, theta5_c_cd, zeta5_c_cd, lambda5_c_cd, &
                        "g5_c_dd", g5_c_dd, eta5_c_dd, theta5_c_dd, zeta5_c_dd, lambda5_c_dd, &
                        num_g2_c_c, num_g2_c_d, num_g5_c_cc, num_g5_c_cd, num_g5_c_dd, &
                        natom, natom_c, &
                        g2_deriv_c_c, g2_deriv_c_d, g5_deriv_c_cc, g5_deriv_c_cd, g5_deriv_c_dd )
       if( myrank == 0 ) write(*,*) " Allocate G2(d-c,d-d),G5(d-cc,d-cd,d-dd): ", elem_d
       call allocate_2( "g2_d_c",  g2_d_c,  eta2_d_c,  rs2_d_c, &
                        "g2_d_d",  g2_d_d,  eta2_d_d,  rs2_d_d, &
                        "g5_d_cc", g5_d_cc, eta5_d_cc, theta5_d_cc, zeta5_d_cc, lambda5_d_cc, &
                        "g5_d_cd", g5_d_cd, eta5_d_cd, theta5_d_cd, zeta5_d_cd, lambda5_d_cd, &
                        "g5_d_dd", g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd, &
                        num_g2_d_c, num_g2_d_d, num_g5_d_cc, num_g5_d_cd, num_g5_d_dd, &
                        natom, natom_d, &
                        g2_deriv_d_c, g2_deriv_d_d, g5_deriv_d_cc, g5_deriv_d_cd, g5_deriv_d_dd )
     else
       if( myrank == 0 ) write(*,*) "Error(main): allocate_2"
       stop
     endif
   endif ! nelem == 2
   ! 3 elements
   if( nelem == 3 )then
     !A,B,C
     if( io_a == 1 .and. io_b == 1 .and. io_c == 1 .and. io_d == 0 )then
       if( myrank == 0 ) write(*,*) " Allocate G2(a-a,a-b,a-c),G5(a-aa,a-ab,a-ac,a-bb,a-bc,a-cc): ", elem_a
       call allocate_3( "g2_a_a",  g2_a_a,  eta2_a_a,  rs2_a_a, &
                        "g2_a_b",  g2_a_b,  eta2_a_b,  rs2_a_b, &
                        "g2_a_c",  g2_a_c,  eta2_a_c,  rs2_a_c, &
                        "g5_a_aa", g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                        "g5_a_ab", g5_a_ab, eta5_a_ab, theta5_a_ab, zeta5_a_ab, lambda5_a_ab, &
                        "g5_a_ac", g5_a_ac, eta5_a_ac, theta5_a_ac, zeta5_a_ac, lambda5_a_ac, &
                        "g5_a_bb", g5_a_bb, eta5_a_bb, theta5_a_bb, zeta5_a_bb, lambda5_a_bb, &
                        "g5_a_bc", g5_a_bc, eta5_a_bc, theta5_a_bc, zeta5_a_bc, lambda5_a_bc, &
                        "g5_a_cc", g5_a_cc, eta5_a_cc, theta5_a_cc, zeta5_a_cc, lambda5_a_cc, &
                        num_g2_a_a,  num_g2_a_b,  num_g2_a_c, &
                        num_g5_a_aa, num_g5_a_ab, num_g5_a_ac, &
                        num_g5_a_bb, num_g5_a_bc, &
                        num_g5_a_cc, &
                        natom, natom_a, &
                        g2_deriv_a_a,  g2_deriv_a_b,  g2_deriv_a_c, &
                        g5_deriv_a_aa, g5_deriv_a_ab, g5_deriv_a_ac, &
                        g5_deriv_a_bb, g5_deriv_a_bc, &
                        g5_deriv_a_cc )
       if( myrank == 0 ) write(*,*) " Allocate G2(b-a,b-b,b-c),G5(b-aa,b-ab,b-ac,b-bb,b-bc,b-cc): ", elem_b
       call allocate_3( "g2_b_a",  g2_b_a,  eta2_b_a,  rs2_b_a, &
                        "g2_b_b",  g2_b_b,  eta2_b_b,  rs2_b_b, &
                        "g2_b_c",  g2_b_c,  eta2_b_c,  rs2_b_c, &
                        "g5_b_aa", g5_b_aa, eta5_b_aa, theta5_b_aa, zeta5_b_aa, lambda5_b_aa, &
                        "g5_b_ab", g5_b_ab, eta5_b_ab, theta5_b_ab, zeta5_b_ab, lambda5_b_ab, &
                        "g5_b_ac", g5_b_ac, eta5_b_ac, theta5_b_ac, zeta5_b_ac, lambda5_b_ac, &
                        "g5_b_bb", g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                        "g5_b_bc", g5_b_bc, eta5_b_bc, theta5_b_bc, zeta5_b_bc, lambda5_b_bc, &
                        "g5_b_cc", g5_b_cc, eta5_b_cc, theta5_b_cc, zeta5_b_cc, lambda5_b_cc, &
                        num_g2_b_a,  num_g2_b_b,  num_g2_b_c, &
                        num_g5_b_aa, num_g5_b_ab, num_g5_b_ac, &
                        num_g5_b_bb, num_g5_b_bc, &
                        num_g5_b_cc, &
                        natom, natom_b, &
                        g2_deriv_b_a,  g2_deriv_b_b,  g2_deriv_b_c, &
                        g5_deriv_b_aa, g5_deriv_b_ab, g5_deriv_b_ac, &
                        g5_deriv_b_bb, g5_deriv_b_bc, &
                        g5_deriv_b_cc )
       if( myrank == 0 ) write(*,*) " Allocate G2(c-a,c-b,c-c),G5(c-aa,c-ab,c-ac,c-bb,c-bc,c-cc): ", elem_c
       call allocate_3( "g2_c_a",  g2_c_a,  eta2_c_a,  rs2_c_a, &
                        "g2_c_b",  g2_c_b,  eta2_c_b,  rs2_c_b, &
                        "g2_c_c",  g2_c_c,  eta2_c_c,  rs2_c_c, &
                        "g5_c_aa", g5_c_aa, eta5_c_aa, theta5_c_aa, zeta5_c_aa, lambda5_c_aa, &
                        "g5_c_ab", g5_c_ab, eta5_c_ab, theta5_c_ab, zeta5_c_ab, lambda5_c_ab, &
                        "g5_c_ac", g5_c_ac, eta5_c_ac, theta5_c_ac, zeta5_c_ac, lambda5_c_ac, &
                        "g5_c_bb", g5_c_bb, eta5_c_bb, theta5_c_bb, zeta5_c_bb, lambda5_c_bb, &
                        "g5_c_bc", g5_c_bc, eta5_c_bc, theta5_c_bc, zeta5_c_bc, lambda5_c_bc, &
                        "g5_c_cc", g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc, &
                        num_g2_c_a,  num_g2_c_b,  num_g2_c_c, &
                        num_g5_c_aa, num_g5_c_ab, num_g5_c_ac, &
                        num_g5_c_bb, num_g5_c_bc, &
                        num_g5_c_cc, &
                        natom, natom_c, &
                        g2_deriv_c_a,  g2_deriv_c_b,  g2_deriv_c_c, &
                        g5_deriv_c_aa, g5_deriv_c_ab, g5_deriv_c_ac, &
                        g5_deriv_c_bb, g5_deriv_c_bc, &
                        g5_deriv_c_cc )
     !A,B,D
     elseif( io_a == 1 .and. io_b == 1 .and. io_c == 0 .and. io_d == 1 )then
       if( myrank == 0 ) write(*,*) " Allocate G2(a-a,a-b,a-d),G5(a-aa,a-ab,a-ad,a-bb,a-bd,a-dd): ", elem_a
       call allocate_3( "g2_a_a",  g2_a_a,  eta2_a_a,  rs2_a_a, &
                        "g2_a_b",  g2_a_b,  eta2_a_b,  rs2_a_b, &
                        "g2_a_d",  g2_a_d,  eta2_a_d,  rs2_a_d, &
                        "g5_a_aa", g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                        "g5_a_ab", g5_a_ab, eta5_a_ab, theta5_a_ab, zeta5_a_ab, lambda5_a_ab, &
                        "g5_a_ad", g5_a_ad, eta5_a_ad, theta5_a_ad, zeta5_a_ad, lambda5_a_ad, &
                        "g5_a_bb", g5_a_bb, eta5_a_bb, theta5_a_bb, zeta5_a_bb, lambda5_a_bb, &
                        "g5_a_bd", g5_a_bd, eta5_a_bd, theta5_a_bd, zeta5_a_bd, lambda5_a_bd, &
                        "g5_a_dd", g5_a_dd, eta5_a_dd, theta5_a_dd, zeta5_a_dd, lambda5_a_dd, &
                        num_g2_a_a,  num_g2_a_b,  num_g2_a_d, &
                        num_g5_a_aa, num_g5_a_ab, num_g5_a_ad, &
                        num_g5_a_bb, num_g5_a_bd, &
                        num_g5_a_dd, &
                        natom, natom_a, &
                        g2_deriv_a_a,  g2_deriv_a_b,  g2_deriv_a_d, &
                        g5_deriv_a_aa, g5_deriv_a_ab, g5_deriv_a_ad, &
                        g5_deriv_a_bb, g5_deriv_a_bd, &
                        g5_deriv_a_dd )
       if( myrank == 0 ) write(*,*) " Allocate G2(b-a,b-b,b-d),G5(b-aa,b-ab,b-ad,b-bb,b-bd,b-dd): ", elem_b
       call allocate_3( "g2_b_a",  g2_b_a,  eta2_b_a,  rs2_b_a, &
                        "g2_b_b",  g2_b_b,  eta2_b_b,  rs2_b_b, &
                        "g2_b_d",  g2_b_d,  eta2_b_d,  rs2_b_d, &
                        "g5_b_aa", g5_b_aa, eta5_b_aa, theta5_b_aa, zeta5_b_aa, lambda5_b_aa, &
                        "g5_b_ab", g5_b_ab, eta5_b_ab, theta5_b_ab, zeta5_b_ab, lambda5_b_ab, &
                        "g5_b_ad", g5_b_ad, eta5_b_ad, theta5_b_ad, zeta5_b_ad, lambda5_b_ad, &
                        "g5_b_bb", g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                        "g5_b_bd", g5_b_bd, eta5_b_bd, theta5_b_bd, zeta5_b_bd, lambda5_b_bd, &
                        "g5_b_dd", g5_b_dd, eta5_b_dd, theta5_b_dd, zeta5_b_dd, lambda5_b_dd, &
                        num_g2_b_a,  num_g2_b_b,  num_g2_b_d, &
                        num_g5_b_aa, num_g5_b_ab, num_g5_b_ad, &
                        num_g5_b_bb, num_g5_b_bd, &
                        num_g5_b_dd, &
                        natom, natom_b, &
                        g2_deriv_b_a,  g2_deriv_b_b,  g2_deriv_b_d, &
                        g5_deriv_b_aa, g5_deriv_b_ab, g5_deriv_b_ad, &
                        g5_deriv_b_bb, g5_deriv_b_bd, &
                        g5_deriv_b_dd )
       if( myrank == 0 ) write(*,*) " Allocate G2(d-a,d-b,d-d),G5(d-aa,d-ab,d-ad,d-bb,d-bd,d-dd): ", elem_d
       call allocate_3( "g2_d_a",  g2_d_a,  eta2_d_a,  rs2_d_a, &
                        "g2_d_b",  g2_d_b,  eta2_d_b,  rs2_d_b, &
                        "g2_d_d",  g2_d_d,  eta2_d_d,  rs2_d_d, &
                        "g5_d_aa", g5_d_aa, eta5_d_aa, theta5_d_aa, zeta5_d_aa, lambda5_d_aa, &
                        "g5_d_ab", g5_d_ab, eta5_d_ab, theta5_d_ab, zeta5_d_ab, lambda5_d_ab, &
                        "g5_d_ad", g5_d_ad, eta5_d_ad, theta5_d_ad, zeta5_d_ad, lambda5_d_ad, &
                        "g5_d_bb", g5_d_bb, eta5_d_bb, theta5_d_bb, zeta5_d_bb, lambda5_d_bb, &
                        "g5_d_bd", g5_d_bd, eta5_d_bd, theta5_d_bd, zeta5_d_bd, lambda5_d_bd, &
                        "g5_d_dd", g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd, &
                        num_g2_d_a,  num_g2_d_b,  num_g2_d_d, &
                        num_g5_d_aa, num_g5_d_ab, num_g5_d_ad, &
                        num_g5_d_bb, num_g5_d_bd, &
                        num_g5_d_dd, &
                        natom, natom_d, &
                        g2_deriv_d_a,  g2_deriv_d_b,  g2_deriv_d_d, &
                        g5_deriv_d_aa, g5_deriv_d_ab, g5_deriv_d_ad, &
                        g5_deriv_d_bb, g5_deriv_d_bd, &
                        g5_deriv_d_dd )
     !A,C,D
     elseif( io_a == 1 .and. io_b == 0 .and. io_c == 1 .and. io_d == 1 )then
       if( myrank == 0 ) write(*,*) " Allocate G2(a-a,a-c,a-d),G5(a-aa,a-ac,a-ad,a-cc,a-cd,a-dd): ", elem_a
       call allocate_3( "g2_a_a",  g2_a_a,  eta2_a_a,  rs2_a_a, &
                        "g2_a_c",  g2_a_c,  eta2_a_c,  rs2_a_c, &
                        "g2_a_d",  g2_a_d,  eta2_a_d,  rs2_a_d, &
                        "g5_a_aa", g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                        "g5_a_ac", g5_a_ac, eta5_a_ac, theta5_a_ac, zeta5_a_ac, lambda5_a_ac, &
                        "g5_a_ad", g5_a_ad, eta5_a_ad, theta5_a_ad, zeta5_a_ad, lambda5_a_ad, &
                        "g5_a_cc", g5_a_cc, eta5_a_cc, theta5_a_cc, zeta5_a_cc, lambda5_a_cc, &
                        "g5_a_cd", g5_a_cd, eta5_a_cd, theta5_a_cd, zeta5_a_cd, lambda5_a_cd, &
                        "g5_a_dd", g5_a_dd, eta5_a_dd, theta5_a_dd, zeta5_a_dd, lambda5_a_dd, &
                        num_g2_a_a,  num_g2_a_c,  num_g2_a_d, &
                        num_g5_a_aa, num_g5_a_ac, num_g5_a_ad, &
                        num_g5_a_cc, num_g5_a_cd, &
                        num_g5_a_dd, &
                        natom, natom_a, &
                        g2_deriv_a_a,  g2_deriv_a_c,  g2_deriv_a_d, &
                        g5_deriv_a_aa, g5_deriv_a_ac, g5_deriv_a_ad, &
                        g5_deriv_a_cc, g5_deriv_a_cd, &
                        g5_deriv_a_dd )
       if( myrank == 0 ) write(*,*) " Allocate G2(c-a,c-c,c-d),G5(c-aa,c-ac,c-ad,c-cc,c-cd,c-dd): ", elem_c
       call allocate_3( "g2_c_a",  g2_c_a,  eta2_c_a,  rs2_c_a, &
                        "g2_c_c",  g2_c_c,  eta2_c_c,  rs2_c_c, &
                        "g2_c_d",  g2_c_d,  eta2_c_d,  rs2_c_d, &
                        "g5_c_aa", g5_c_aa, eta5_c_aa, theta5_c_aa, zeta5_c_aa, lambda5_c_aa, &
                        "g5_c_ac", g5_c_ac, eta5_c_ac, theta5_c_ac, zeta5_c_ac, lambda5_c_ac, &
                        "g5_c_ad", g5_c_ad, eta5_c_ad, theta5_c_ad, zeta5_c_ad, lambda5_c_ad, &
                        "g5_c_cc", g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc, &
                        "g5_c_cd", g5_c_cd, eta5_c_cd, theta5_c_cd, zeta5_c_cd, lambda5_c_cd, &
                        "g5_c_dd", g5_c_dd, eta5_c_dd, theta5_c_dd, zeta5_c_dd, lambda5_c_dd, &
                        num_g2_c_a,  num_g2_c_c,  num_g2_c_d, &
                        num_g5_c_aa, num_g5_c_ac, num_g5_c_ad, &
                        num_g5_c_cc, num_g5_c_cd, &
                        num_g5_c_dd, &
                        natom, natom_c, &
                        g2_deriv_c_a,  g2_deriv_c_c,  g2_deriv_c_d, &
                        g5_deriv_c_aa, g5_deriv_c_ac, g5_deriv_c_ad, &
                        g5_deriv_c_cc, g5_deriv_c_cd, &
                        g5_deriv_c_dd )
       if( myrank == 0 ) write(*,*) " Allocate G2(d-a,d-c,d-d),G5(d-aa,d-ac,d-ad,d-cc,d-cd,d-dd): ", elem_d
       call allocate_3( "g2_d_a",  g2_d_a,  eta2_d_a,  rs2_d_a, &
                        "g2_d_c",  g2_d_c,  eta2_d_c,  rs2_d_c, &
                        "g2_d_d",  g2_d_d,  eta2_d_d,  rs2_d_d, &
                        "g5_d_aa", g5_d_aa, eta5_d_aa, theta5_d_aa, zeta5_d_aa, lambda5_d_aa, &
                        "g5_d_ac", g5_d_ac, eta5_d_ac, theta5_d_ac, zeta5_d_ac, lambda5_d_ac, &
                        "g5_d_ad", g5_d_ad, eta5_d_ad, theta5_d_ad, zeta5_d_ad, lambda5_d_ad, &
                        "g5_d_cc", g5_d_cc, eta5_d_cc, theta5_d_cc, zeta5_d_cc, lambda5_d_cc, &
                        "g5_d_cd", g5_d_cd, eta5_d_cd, theta5_d_cd, zeta5_d_cd, lambda5_d_cd, &
                        "g5_d_dd", g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd, &
                        num_g2_d_a,  num_g2_d_c,  num_g2_d_d, &
                        num_g5_d_aa, num_g5_d_ac, num_g5_d_ad, &
                        num_g5_d_cc, num_g5_d_cd, &
                        num_g5_d_dd, &
                        natom, natom_d, &
                        g2_deriv_d_a,  g2_deriv_d_c,  g2_deriv_d_d, &
                        g5_deriv_d_aa, g5_deriv_d_ac, g5_deriv_d_ad, &
                        g5_deriv_d_cc, g5_deriv_d_cd, &
                        g5_deriv_d_dd )
     !B,C,D
     elseif( io_a == 0 .and. io_b == 1 .and. io_c == 1 .and. io_d == 1 )then
       if( myrank == 0 ) write(*,*) " Allocate G2(b-b,b-c,b-d),G5(b-bb,b-bc,b-bd,b-cc,b-cd,b-dd): ", elem_b
       call allocate_3( "g2_b_b",  g2_b_b,  eta2_b_b,  rs2_b_b, &
                        "g2_b_c",  g2_b_c,  eta2_b_c,  rs2_b_c, &
                        "g2_b_d",  g2_b_d,  eta2_b_d,  rs2_b_d, &
                        "g5_b_bb", g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                        "g5_b_bc", g5_b_bc, eta5_b_bc, theta5_b_bc, zeta5_b_bc, lambda5_b_bc, &
                        "g5_b_bd", g5_b_bd, eta5_b_bd, theta5_b_bd, zeta5_b_bd, lambda5_b_bd, &
                        "g5_b_cc", g5_b_cc, eta5_b_cc, theta5_b_cc, zeta5_b_cc, lambda5_b_cc, &
                        "g5_b_cd", g5_b_cd, eta5_b_cd, theta5_b_cd, zeta5_b_cd, lambda5_b_cd, &
                        "g5_b_dd", g5_b_dd, eta5_b_dd, theta5_b_dd, zeta5_b_dd, lambda5_b_dd, &
                        num_g2_b_b,  num_g2_b_c,  num_g2_b_d, &
                        num_g5_b_bb, num_g5_b_bc, num_g5_b_bd, &
                        num_g5_b_cc, num_g5_b_cd, &
                        num_g5_b_dd, &
                        natom, natom_b, &
                        g2_deriv_b_b,  g2_deriv_b_c,  g2_deriv_b_d, &
                        g5_deriv_b_bb, g5_deriv_b_bc, g5_deriv_b_bd, &
                        g5_deriv_b_cc, g5_deriv_b_cd, &
                        g5_deriv_b_dd )
       if( myrank == 0 ) write(*,*) " Allocate G2(c-b,c-c,c-d),G5(c-bb,c-bc,c-bd,c-cc,c-cd,c-dd): ", elem_c
       call allocate_3( "g2_c_b",  g2_c_b,  eta2_c_b,  rs2_c_b, &
                        "g2_c_c",  g2_c_c,  eta2_c_c,  rs2_c_c, &
                        "g2_c_d",  g2_c_d,  eta2_c_d,  rs2_c_d, &
                        "g5_c_bb", g5_c_bb, eta5_c_bb, theta5_c_bb, zeta5_c_bb, lambda5_c_bb, &
                        "g5_c_bc", g5_c_bc, eta5_c_bc, theta5_c_bc, zeta5_c_bc, lambda5_c_bc, &
                        "g5_c_bd", g5_c_bd, eta5_c_bd, theta5_c_bd, zeta5_c_bd, lambda5_c_bd, &
                        "g5_c_cc", g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc, &
                        "g5_c_cd", g5_c_cd, eta5_c_cd, theta5_c_cd, zeta5_c_cd, lambda5_c_cd, &
                        "g5_c_dd", g5_c_dd, eta5_c_dd, theta5_c_dd, zeta5_c_dd, lambda5_c_dd, &
                        num_g2_c_b,  num_g2_c_c,  num_g2_c_d, &
                        num_g5_c_bb, num_g5_c_bc, num_g5_c_bd, &
                        num_g5_c_cc, num_g5_c_cd, &
                        num_g5_c_dd, &
                        natom, natom_c, &
                        g2_deriv_c_b,  g2_deriv_c_c,  g2_deriv_c_d, &
                        g5_deriv_c_bb, g5_deriv_c_bc, g5_deriv_c_bd, &
                        g5_deriv_c_cc, g5_deriv_c_cd, &
                        g5_deriv_c_dd )
       if( myrank == 0 ) write(*,*) " Allocate G2(d-b,d-c,d-d),G5(d-bb,d-bc,d-bd,d-cc,d-cd,d-dd): ", elem_d
       call allocate_3( "g2_d_b",  g2_d_b,  eta2_d_b,  rs2_d_b, &
                        "g2_d_c",  g2_d_c,  eta2_d_c,  rs2_d_c, &
                        "g2_d_d",  g2_d_d,  eta2_d_d,  rs2_d_d, &
                        "g5_d_bb", g5_d_bb, eta5_d_bb, theta5_d_bb, zeta5_d_bb, lambda5_d_bb, &
                        "g5_d_bc", g5_d_bc, eta5_d_bc, theta5_d_bc, zeta5_d_bc, lambda5_d_bc, &
                        "g5_d_bd", g5_d_bd, eta5_d_bd, theta5_d_bd, zeta5_d_bd, lambda5_d_bd, &
                        "g5_d_cc", g5_d_cc, eta5_d_cc, theta5_d_cc, zeta5_d_cc, lambda5_d_cc, &
                        "g5_d_cd", g5_d_cd, eta5_d_cd, theta5_d_cd, zeta5_d_cd, lambda5_d_cd, &
                        "g5_d_dd", g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd, &
                        num_g2_d_b,  num_g2_d_c,  num_g2_d_d, &
                        num_g5_d_bb, num_g5_d_bc, num_g5_d_bd, &
                        num_g5_d_cc, num_g5_d_cd, &
                        num_g5_d_dd, &
                        natom, natom_d, &
                        g2_deriv_d_b,  g2_deriv_d_c,  g2_deriv_d_d, &
                        g5_deriv_d_bb, g5_deriv_d_bc, g5_deriv_d_bd, &
                        g5_deriv_d_cc, g5_deriv_d_cd, &
                        g5_deriv_d_dd )
     else
       write(*,*) "Error(main): allocate_3"
       stop
     endif 
   endif! nelem == 3
   ! 4 elements
   if( nelem == 4 )then
     !A,B,C,D
     if( io_a == 1 .and. io_b == 1 .and. io_c == 1 .and. io_d == 1 )then
       if( myrank == 0 ) write(*,*) " Allocate &
                  &G2(a-a,a-b,a-c,a-d),&
                  &G5(a-aa,a-ab,a-ac,a-ad,a-bb,a-bc,a-bd,a-cc,a-cd,a-dd): ", elem_a
       call allocate_4( "g2_a_a",  g2_a_a,  eta2_a_a,  rs2_a_a, &
                        "g2_a_b",  g2_a_b,  eta2_a_b,  rs2_a_b, &
                        "g2_a_c",  g2_a_c,  eta2_a_c,  rs2_a_c, &
                        "g2_a_d",  g2_a_d,  eta2_a_d,  rs2_a_d, &
                        "g5_a_aa", g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                        "g5_a_ab", g5_a_ab, eta5_a_ab, theta5_a_ab, zeta5_a_ab, lambda5_a_ab, &
                        "g5_a_ac", g5_a_ac, eta5_a_ac, theta5_a_ac, zeta5_a_ac, lambda5_a_ac, &
                        "g5_a_ad", g5_a_ad, eta5_a_ad, theta5_a_ad, zeta5_a_ad, lambda5_a_ad, &
                        "g5_a_bb", g5_a_bb, eta5_a_bb, theta5_a_bb, zeta5_a_bb, lambda5_a_bb, &
                        "g5_a_bc", g5_a_bc, eta5_a_bc, theta5_a_bc, zeta5_a_bc, lambda5_a_bc, &
                        "g5_a_bd", g5_a_bd, eta5_a_bd, theta5_a_bd, zeta5_a_bd, lambda5_a_bd, &
                        "g5_a_cc", g5_a_cc, eta5_a_cc, theta5_a_cc, zeta5_a_cc, lambda5_a_cc, &
                        "g5_a_cd", g5_a_cd, eta5_a_cd, theta5_a_cd, zeta5_a_cd, lambda5_a_cd, &
                        "g5_a_dd", g5_a_dd, eta5_a_dd, theta5_a_dd, zeta5_a_dd, lambda5_a_dd, &
                        num_g2_a_a,  num_g2_a_b,  num_g2_a_c,  num_g2_a_d, &
                        num_g5_a_aa, num_g5_a_ab, num_g5_a_ac, num_g5_a_ad, &
                        num_g5_a_bb, num_g5_a_bc, num_g5_a_bd, &
                        num_g5_a_cc, num_g5_a_cd, &
                        num_g5_a_dd, &
                        natom, natom_a, &
                        g2_deriv_a_a,  g2_deriv_a_b,  g2_deriv_a_c,  g2_deriv_a_d, &
                        g5_deriv_a_aa, g5_deriv_a_ab, g5_deriv_a_ac, g5_deriv_a_ad, &
                        g5_deriv_a_bb, g5_deriv_a_bc, g5_deriv_a_bd, &
                        g5_deriv_a_cc, g5_deriv_a_cd, &
                        g5_deriv_a_dd )
       if( myrank == 0 ) write(*,*) " Allocate &
                   &G2(b-a,b-b,b-c,b-d),&
                   &G5(b-aa,b-ab,b-ac,b-ad,b-bb,b-bc,b-bd,b-cc,b-cd,b-dd): ", elem_b
       call allocate_4( "g2_b_a",  g2_b_a,  eta2_b_a,  rs2_b_a, &
                        "g2_b_b",  g2_b_b,  eta2_b_b,  rs2_b_b, &
                        "g2_b_c",  g2_b_c,  eta2_b_c,  rs2_b_c, &
                        "g2_b_d",  g2_b_d,  eta2_b_d,  rs2_b_d, &
                        "g5_b_aa", g5_b_aa, eta5_b_aa, theta5_b_aa, zeta5_b_aa, lambda5_b_aa, &
                        "g5_b_ab", g5_b_ab, eta5_b_ab, theta5_b_ab, zeta5_b_ab, lambda5_b_ab, &
                        "g5_b_ac", g5_b_ac, eta5_b_ac, theta5_b_ac, zeta5_b_ac, lambda5_b_ac, &
                        "g5_b_ad", g5_b_ad, eta5_b_ad, theta5_b_ad, zeta5_b_ad, lambda5_b_ad, &
                        "g5_b_bb", g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                        "g5_b_bc", g5_b_bc, eta5_b_bc, theta5_b_bc, zeta5_b_bc, lambda5_b_bc, &
                        "g5_b_bd", g5_b_bd, eta5_b_bd, theta5_b_bd, zeta5_b_bd, lambda5_b_bd, &
                        "g5_b_cc", g5_b_cc, eta5_b_cc, theta5_b_cc, zeta5_b_cc, lambda5_b_cc, &
                        "g5_b_cd", g5_b_cd, eta5_b_cd, theta5_b_cd, zeta5_b_cd, lambda5_b_cd, &
                        "g5_b_dd", g5_b_dd, eta5_b_dd, theta5_b_dd, zeta5_b_dd, lambda5_b_dd, &
                        num_g2_b_a,  num_g2_b_b,  num_g2_b_c,  num_g2_b_d, &
                        num_g5_b_aa, num_g5_b_ab, num_g5_b_ac, num_g5_b_ad, &
                        num_g5_b_bb, num_g5_b_bc, num_g5_b_bd, &
                        num_g5_b_cc, num_g5_b_cd, &
                        num_g5_b_dd, &
                        natom, natom_b, &
                        g2_deriv_b_a,  g2_deriv_b_b,  g2_deriv_b_c,  g2_deriv_b_d, &
                        g5_deriv_b_aa, g5_deriv_b_ab, g5_deriv_b_ac, g5_deriv_b_ad, &
                        g5_deriv_b_bb, g5_deriv_b_bc, g5_deriv_b_bd, &
                        g5_deriv_b_cc, g5_deriv_b_cd, &
                        g5_deriv_b_dd )
       if( myrank == 0 ) write(*,*) " Allocate &
                   &G2(c-a,c-b,c-c,c-d),&
                   &G5(c-aa,c-ab,c-ac,c-ad,c-bb,c-bc,c-bd,c-cc,c-cd,c-dd): ", elem_c
       call allocate_4( "g2_c_a",  g2_c_a,  eta2_c_a,  rs2_c_a, &
                        "g2_c_b",  g2_c_b,  eta2_c_b,  rs2_c_b, &
                        "g2_c_c",  g2_c_c,  eta2_c_c,  rs2_c_c, &
                        "g2_c_d",  g2_c_d,  eta2_c_d,  rs2_c_d, &
                        "g5_c_aa", g5_c_aa, eta5_c_aa, theta5_c_aa, zeta5_c_aa, lambda5_c_aa, &
                        "g5_c_ab", g5_c_ab, eta5_c_ab, theta5_c_ab, zeta5_c_ab, lambda5_c_ab, &
                        "g5_c_ac", g5_c_ac, eta5_c_ac, theta5_c_ac, zeta5_c_ac, lambda5_c_ac, &
                        "g5_c_ad", g5_c_ad, eta5_c_ad, theta5_c_ad, zeta5_c_ad, lambda5_c_ad, &
                        "g5_c_bb", g5_c_bb, eta5_c_bb, theta5_c_bb, zeta5_c_bb, lambda5_c_bb, &
                        "g5_c_bc", g5_c_bc, eta5_c_bc, theta5_c_bc, zeta5_c_bc, lambda5_c_bc, &
                        "g5_c_bd", g5_c_bd, eta5_c_bd, theta5_c_bd, zeta5_c_bd, lambda5_c_bd, &
                        "g5_c_cc", g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc, &
                        "g5_c_cd", g5_c_cd, eta5_c_cd, theta5_c_cd, zeta5_c_cd, lambda5_c_cd, &
                        "g5_c_dd", g5_c_dd, eta5_c_dd, theta5_c_dd, zeta5_c_dd, lambda5_c_dd, &
                        num_g2_c_a,  num_g2_c_b,  num_g2_c_c,  num_g2_c_d, &
                        num_g5_c_aa, num_g5_c_ab, num_g5_c_ac, num_g5_c_ad, &
                        num_g5_c_bb, num_g5_c_bc, num_g5_c_bd, &
                        num_g5_c_cc, num_g5_c_cd, &
                        num_g5_c_dd, &
                        natom, natom_c, &
                        g2_deriv_c_a,  g2_deriv_c_b,  g2_deriv_c_c,  g2_deriv_c_d, &
                        g5_deriv_c_aa, g5_deriv_c_ab, g5_deriv_c_ac, g5_deriv_c_ad, &
                        g5_deriv_c_bb, g5_deriv_c_bc, g5_deriv_c_bd, &
                        g5_deriv_c_cc, g5_deriv_c_cd, &
                        g5_deriv_c_dd )
       if( myrank == 0 ) write(*,*) " Allocate &
                   &G2(d-a,d-b,d-c,d-d),&
                   &G5(d-aa,d-ab,d-ac,d-ad,d-bb,d-bc,d-bd,d-cc,d-cd,d-dd): ", elem_d
       call allocate_4( "g2_d_a",  g2_d_a,  eta2_d_a,  rs2_d_a, &
                        "g2_d_b",  g2_d_b,  eta2_d_b,  rs2_d_b, &
                        "g2_d_c",  g2_d_c,  eta2_d_c,  rs2_d_c, &
                        "g2_d_d",  g2_d_d,  eta2_d_d,  rs2_d_d, &
                        "g5_d_aa", g5_d_aa, eta5_d_aa, theta5_d_aa, zeta5_d_aa, lambda5_d_aa, &
                        "g5_d_ab", g5_d_ab, eta5_d_ab, theta5_d_ab, zeta5_d_ab, lambda5_d_ab, &
                        "g5_d_ac", g5_d_ac, eta5_d_ac, theta5_d_ac, zeta5_d_ac, lambda5_d_ac, &
                        "g5_d_ad", g5_d_ad, eta5_d_ad, theta5_d_ad, zeta5_d_ad, lambda5_d_ad, &
                        "g5_d_bb", g5_d_bb, eta5_d_bb, theta5_d_bb, zeta5_d_bb, lambda5_d_bb, &
                        "g5_d_bc", g5_d_bc, eta5_d_bc, theta5_d_bc, zeta5_d_bc, lambda5_d_bc, &
                        "g5_d_bd", g5_d_bd, eta5_d_bd, theta5_d_bd, zeta5_d_bd, lambda5_d_bd, &
                        "g5_d_cc", g5_d_cc, eta5_d_cc, theta5_d_cc, zeta5_d_cc, lambda5_d_cc, &
                        "g5_d_cd", g5_d_cd, eta5_d_cd, theta5_d_cd, zeta5_d_cd, lambda5_d_cd, &
                        "g5_d_dd", g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd, &
                        num_g2_d_a,  num_g2_d_b,  num_g2_d_c,  num_g2_d_d, &
                        num_g5_d_aa, num_g5_d_ab, num_g5_d_ac, num_g5_d_ad, &
                        num_g5_d_bb, num_g5_d_bc, num_g5_d_bd, &
                        num_g5_d_cc, num_g5_d_cd, &
                        num_g5_d_dd, &
                        natom, natom_d, &
                        g2_deriv_d_a,  g2_deriv_d_b,  g2_deriv_d_c,  g2_deriv_d_d, &
                        g5_deriv_d_aa, g5_deriv_d_ab, g5_deriv_d_ac, g5_deriv_d_ad, &
                        g5_deriv_d_bb, g5_deriv_d_bc, g5_deriv_d_bd, &
                        g5_deriv_d_cc, g5_deriv_d_cd, &
                        g5_deriv_d_dd )
     else
       write(*,*) "Error(main): allocate_4"
       stop
     endif
   endif ! nelem == 4
   ! more than 5 elements
   if( nelem > 4 )then
     write(*,*) "Error(main): allocate"
     stop
   endif ! nelem > 4


   !*************************************************************
   ! Read positions file
   !*************************************************************

   !Basis
   if( myrank == 0 )then
     read(2) lattice
     read(2) x
     read(2) y
     read(2) z
   endif ! myrank == 0
   call mpi_barrier( mpi_comm_world, ierr )
   call mpi_bcast( lattice, 1, mpi_double_precision, 0, mpi_comm_world, ierr )
   call mpi_bcast( x, 3, mpi_double_precision, 0, mpi_comm_world, ierr )
   call mpi_bcast( y, 3, mpi_double_precision, 0, mpi_comm_world, ierr )
   call mpi_bcast( z, 3, mpi_double_precision, 0, mpi_comm_world, ierr )
   call mpi_barrier( mpi_comm_world, ierr )

   !Number of atoms and elements
   if( myrank == 0 )then
     ! 1 element
     if( nelem == 1 )then
       read(2) elem_1
       read(2) natom_1
     ! 2 elements
     elseif( nelem == 2 )then
       read(2) elem_1,  elem_2
       read(2) natom_1, natom_2
     ! 3 elements
     elseif( nelem == 3 )then
       read(2) elem_1,  elem_2,  elem_3
       read(2) natom_1, natom_2, natom_3
     ! 4 elements
     elseif( nelem == 4 )then
       read(2) elem_1,  elem_2,  elem_3,  elem_4
       read(2) natom_1, natom_2, natom_3, natom_4
     ! more than 5 elements
     else
       write(*,*) "Error(main): read POSCAR"
       stop
     endif
   endif ! myrank == 0

   call mpi_barrier( mpi_comm_world, ierr )
   call mpi_bcast( elem_1, 4, mpi_character, 0, mpi_comm_world, ierr )
   call mpi_bcast( natom_1, 1, mpi_integer, 0, mpi_comm_world, ierr )
   call mpi_bcast( elem_2, 4, mpi_character, 0, mpi_comm_world, ierr )
   call mpi_bcast( natom_2, 1, mpi_integer, 0, mpi_comm_world, ierr )
   call mpi_bcast( elem_3, 4, mpi_character, 0, mpi_comm_world, ierr )
   call mpi_bcast( natom_3, 1, mpi_integer, 0, mpi_comm_world, ierr )
   call mpi_bcast( elem_4, 4, mpi_character, 0, mpi_comm_world, ierr )
   call mpi_bcast( natom_4, 1, mpi_integer, 0, mpi_comm_world, ierr )
   call mpi_barrier( mpi_comm_world, ierr )

   xlattice = dsqrt( ( lattice**2 )*( x(1)*x(1) + x(2)*x(2) + x(3)*x(3) ) )
   ylattice = dsqrt( ( lattice**2 )*( y(1)*y(1) + y(2)*y(2) + y(3)*y(3) ) )
   zlattice = dsqrt( ( lattice**2 )*( z(1)*z(1) + z(2)*z(2) + z(3)*z(3) ) )

   allocate( posi( natom, 3 ) )


   !*************************************************************
   ! Structure(MD) loop
   !*************************************************************

   if( myrank == 0 )then
     write(*,*) " "
     write(*,*) " * SF calculations start *"
     write(*,*) " "
   endif ! myrank == 0 

   io_sf = 0

   do imd = 1, md_steps

     if( myrank == 0 )then
       if( dble(imd)/dble(md_steps) == 0.10d0 ) write(*,*) " progress: 10 %",  imd, "/", md_steps
       if( dble(imd)/dble(md_steps) == 0.20d0 ) write(*,*) " progress: 20 %",  imd, "/", md_steps
       if( dble(imd)/dble(md_steps) == 0.30d0 ) write(*,*) " progress: 30 %",  imd, "/", md_steps
       if( dble(imd)/dble(md_steps) == 0.40d0 ) write(*,*) " progress: 40 %",  imd, "/", md_steps
       if( dble(imd)/dble(md_steps) == 0.50d0 ) write(*,*) " progress: 50 %",  imd, "/", md_steps
       if( dble(imd)/dble(md_steps) == 0.60d0 ) write(*,*) " progress: 60 %",  imd, "/", md_steps
       if( dble(imd)/dble(md_steps) == 0.70d0 ) write(*,*) " progress: 70 %",  imd, "/", md_steps
       if( dble(imd)/dble(md_steps) == 0.80d0 ) write(*,*) " progress: 80 %",  imd, "/", md_steps
       if( dble(imd)/dble(md_steps) == 0.90d0 ) write(*,*) " progress: 90 %",  imd, "/", md_steps
       if( dble(imd)/dble(md_steps) == 1.0d0  ) write(*,*) " progress: 100 %", imd, "/", md_steps
     endif ! myrank == 0 

     !Read positions
     if( myrank == 0 )then
       read(2) dir  ! Direct configuration
       read(2) posi ! posi(i,1), posi(i,2), posi(i,3)
     endif ! myrank == 0 

     call mpi_barrier( mpi_comm_world, ierr )
     call mpi_bcast( dir, 6, mpi_character, 0, mpi_comm_world, ierr )
     call mpi_bcast( posi, natom*3, mpi_double_precision, 0, mpi_comm_world, ierr )
     call mpi_barrier( mpi_comm_world, ierr )

     !*************************************************************
     ! Symmetry function calculation, 
     ! G2 & G5, G2_deriv & G5_deriv
     !*************************************************************

     ! 1 element
     if( nelem == 1 )then
       !A
       if( io_a == 1 .and. io_b == 0 .and. io_c == 0 .and. io_d == 0 )then
         if( myrank == 0 .and. imd == 1 ) write(*,*) " SF_1", nelem, elem_a
         call sf_1( lattice, x, y, z, posi, xlattice, ylattice, zlattice, &
                    g2_a_a, eta2_a_a, rs2_a_a, &
                    g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                    g2_deriv_a_a, g5_deriv_a_aa, &
                    num_g2_a_a, num_g5_a_aa, &
                    natom_a, &
                    natom, rc, neighbor, io_sf, num_mpi, myrank )
       !B
       elseif( io_a == 0 .and. io_b == 1 .and. io_c == 0 .and. io_d == 0 )then
         if( myrank == 0 .and. imd == 1 ) write(*,*) " SF_1", nelem, elem_b
         call sf_1( lattice, x, y, z, posi, xlattice, ylattice, zlattice, &
                    g2_b_b, eta2_b_b, rs2_b_b, &
                    g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                    g2_deriv_b_b, g5_deriv_b_bb, &
                    num_g2_b_b, num_g5_b_bb, &
                    natom_b, &
                    natom, rc, neighbor, io_sf, num_mpi, myrank )
       !C
       elseif( io_a == 0 .and. io_b == 0 .and. io_c == 1 .and. io_d == 0 )then
         if( myrank == 0 .and. imd == 1 ) write(*,*) " SF_1", nelem, elem_c
         call sf_1( lattice, x, y, z, posi, xlattice, ylattice, zlattice, &
                    g2_c_c, eta2_c_c, rs2_c_c, &
                    g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc, &
                    g2_deriv_c_c, g5_deriv_c_cc, &
                    num_g2_c_c, num_g5_c_cc, &
                    natom_c, &
                    natom, rc, neighbor, io_sf, num_mpi, myrank )
       !D
       elseif( io_a == 0 .and. io_b == 0 .and. io_c == 0 .and. io_d == 1 )then
         if( myrank == 0 .and. imd == 1 ) write(*,*) " SF_1", nelem, elem_d
         call sf_1( lattice, x, y, z, posi, xlattice, ylattice, zlattice, &
                    g2_d_d, eta2_d_d, rs2_d_d, &
                    g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd, &
                    g2_deriv_d_d, g5_deriv_d_dd, &
                    num_g2_d_d, num_g5_d_dd, &
                    natom_d, &
                    natom, rc, neighbor, io_sf, num_mpi, myrank )
       else
         write(*,*) "Error(main): sf_1"
         stop
       endif
     endif ! nelem == 1
     ! 2 elements
     if( nelem == 2 )then
       !A,B
       if( io_a == 1 .and. io_b == 1 .and. io_c == 0 .and. io_d == 0 )then
         if( myrank == 0 .and. imd == 1 ) write(*,*) " SF_2", nelem, elem_a, elem_b
         call sf_2( lattice, x, y, z, posi, xlattice, ylattice, zlattice, &
                    g2_a_a,  eta2_a_a,  rs2_a_a, &
                    g2_a_b,  eta2_a_b,  rs2_a_b, &
                    g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                    g5_a_ab, eta5_a_ab, theta5_a_ab, zeta5_a_ab, lambda5_a_ab, &
                    g5_a_bb, eta5_a_bb, theta5_a_bb, zeta5_a_bb, lambda5_a_bb, &
                    g2_b_a,  eta2_b_a,  rs2_b_a, &
                    g2_b_b,  eta2_b_b,  rs2_b_b, &
                    g5_b_aa, eta5_b_aa, theta5_b_aa, zeta5_b_aa, lambda5_b_aa, &
                    g5_b_ab, eta5_b_ab, theta5_b_ab, zeta5_b_ab, lambda5_b_ab, &
                    g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                    g2_deriv_a_a,  g2_deriv_a_b, &
                    g5_deriv_a_aa, g5_deriv_a_ab, g5_deriv_a_bb, &
                    g2_deriv_b_a,  g2_deriv_b_b, &
                    g5_deriv_b_aa, g5_deriv_b_ab, g5_deriv_b_bb, &
                    num_g2_a_a,  num_g2_a_b, &
                    num_g5_a_aa, num_g5_a_ab, num_g5_a_bb, &
                    num_g2_b_a,  num_g2_b_b, &
                    num_g5_b_aa, num_g5_b_ab, num_g5_b_bb, &
                    natom_a, natom_b, &
                    natom, rc, neighbor, io_sf, num_mpi, myrank )
       !A,C
       elseif( io_a == 1 .and. io_b == 0 .and. io_c == 1 .and. io_d == 0 )then
         if( myrank == 0 .and. imd == 1 ) write(*,*) " SF_2", nelem, elem_a, elem_c
         call sf_2( lattice, x, y, z, posi, xlattice, ylattice, zlattice, &
                    g2_a_a,  eta2_a_a,  rs2_a_a, &
                    g2_a_c,  eta2_a_c,  rs2_a_c, &
                    g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                    g5_a_ac, eta5_a_ac, theta5_a_ac, zeta5_a_ac, lambda5_a_ac, &
                    g5_a_cc, eta5_a_cc, theta5_a_cc, zeta5_a_cc, lambda5_a_cc, &
                    g2_c_a,  eta2_c_a,  rs2_c_a, &
                    g2_c_c,  eta2_c_c,  rs2_c_c, &
                    g5_c_aa, eta5_c_aa, theta5_c_aa, zeta5_c_aa, lambda5_c_aa, &
                    g5_c_ac, eta5_c_ac, theta5_c_ac, zeta5_c_ac, lambda5_c_ac, &
                    g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc, &
                    g2_deriv_a_a,  g2_deriv_a_c, &
                    g5_deriv_a_aa, g5_deriv_a_ac, g5_deriv_a_cc, &
                    g2_deriv_c_a,  g2_deriv_c_c, &
                    g5_deriv_c_aa, g5_deriv_c_ac, g5_deriv_c_cc, &
                    num_g2_a_a,  num_g2_a_c, &
                    num_g5_a_aa, num_g5_a_ac, num_g5_a_cc, &
                    num_g2_c_a,  num_g2_c_c, &
                    num_g5_c_aa, num_g5_c_ac, num_g5_c_cc, &
                    natom_a, natom_c, &
                    natom, rc, neighbor, io_sf, num_mpi, myrank )
       !A,D
       elseif( io_a == 1 .and. io_b == 0 .and. io_c == 0 .and. io_d == 1 )then
         if( myrank == 0 .and. imd == 1 ) write(*,*) " SF_2", nelem, elem_a, elem_d
         call sf_2( lattice, x, y, z, posi, xlattice, ylattice, zlattice, &
                    g2_a_a,  eta2_a_a,  rs2_a_a, &
                    g2_a_d,  eta2_a_d,  rs2_a_d, &
                    g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                    g5_a_ad, eta5_a_ad, theta5_a_ad, zeta5_a_ad, lambda5_a_ad, &
                    g5_a_dd, eta5_a_dd, theta5_a_dd, zeta5_a_dd, lambda5_a_dd, &
                    g2_d_a,  eta2_d_a,  rs2_d_a, &
                    g2_d_d,  eta2_d_d,  rs2_d_d, &
                    g5_d_aa, eta5_d_aa, theta5_d_aa, zeta5_d_aa, lambda5_d_aa, &
                    g5_d_ad, eta5_d_ad, theta5_d_ad, zeta5_d_ad, lambda5_d_ad, &
                    g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd, &
                    g2_deriv_a_a,  g2_deriv_a_d, &
                    g5_deriv_a_aa, g5_deriv_a_ad, g5_deriv_a_dd, &
                    g2_deriv_d_a,  g2_deriv_d_d, &
                    g5_deriv_d_aa, g5_deriv_d_ad, g5_deriv_d_dd, &
                    num_g2_a_a,  num_g2_a_d, &
                    num_g5_a_aa, num_g5_a_ad, num_g5_a_dd, &
                    num_g2_d_a,  num_g2_d_d, &
                    num_g5_d_aa, num_g5_d_ad, num_g5_d_dd, &
                    natom_a, natom_d, &
                    natom, rc, neighbor, io_sf, num_mpi, myrank )
       !B,C
       elseif( io_a == 0 .and. io_b == 1 .and. io_c == 1 .and. io_d == 0 )then
         if( myrank == 0 .and. imd == 1 ) write(*,*) " SF_2", nelem, elem_b, elem_c
         call sf_2( lattice, x, y, z, posi, xlattice, ylattice, zlattice, &
                    g2_b_b,  eta2_b_b,  rs2_b_b, &
                    g2_b_c,  eta2_b_c,  rs2_b_c, &
                    g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                    g5_b_bc, eta5_b_bc, theta5_b_bc, zeta5_b_bc, lambda5_b_bc, &
                    g5_b_cc, eta5_b_cc, theta5_b_cc, zeta5_b_cc, lambda5_b_cc, &
                    g2_c_b,  eta2_c_b,  rs2_c_b, &
                    g2_c_c,  eta2_c_c,  rs2_c_c, &
                    g5_c_bb, eta5_c_bb, theta5_c_bb, zeta5_c_bb, lambda5_c_bb, &
                    g5_c_bc, eta5_c_bc, theta5_c_bc, zeta5_c_bc, lambda5_c_bc, &
                    g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc, &
                    g2_deriv_b_b,  g2_deriv_b_c, &
                    g5_deriv_b_bb, g5_deriv_b_bc, g5_deriv_b_cc, &
                    g2_deriv_c_b,  g2_deriv_c_c, &
                    g5_deriv_c_bb, g5_deriv_c_bc, g5_deriv_c_cc, &
                    num_g2_b_b,  num_g2_b_c, &
                    num_g5_b_bb, num_g5_b_bc, num_g5_b_cc, &
                    num_g2_c_b,  num_g2_c_c, &
                    num_g5_c_bb, num_g5_c_bc, num_g5_c_cc, &
                    natom_b, natom_c, &
                    natom, rc, neighbor, io_sf, num_mpi, myrank )
       !B,D
       elseif( io_a == 0 .and. io_b == 1 .and. io_c == 0 .and. io_d == 1 )then
         if( myrank == 0 .and. imd == 1 ) write(*,*) " SF_2", nelem, elem_b, elem_d
         call sf_2( lattice, x, y, z, posi, xlattice, ylattice, zlattice, &
                    g2_b_b,  eta2_b_b,  rs2_b_b, &
                    g2_b_d,  eta2_b_d,  rs2_b_d, &
                    g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                    g5_b_bd, eta5_b_bd, theta5_b_bd, zeta5_b_bd, lambda5_b_bd, &
                    g5_b_dd, eta5_b_dd, theta5_b_dd, zeta5_b_dd, lambda5_b_dd, &
                    g2_d_b,  eta2_d_b,  rs2_d_b, &
                    g2_d_d,  eta2_d_d,  rs2_d_d, &
                    g5_d_bb, eta5_d_bb, theta5_d_bb, zeta5_d_bb, lambda5_d_bb, &
                    g5_d_bd, eta5_d_bd, theta5_d_bd, zeta5_d_bd, lambda5_d_bd, &
                    g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd, &
                    g2_deriv_b_b,  g2_deriv_b_d, &
                    g5_deriv_b_bb, g5_deriv_b_bd, g5_deriv_b_dd, &
                    g2_deriv_d_b,  g2_deriv_d_d, &
                    g5_deriv_d_bb, g5_deriv_d_bd, g5_deriv_d_dd, &
                    num_g2_b_b,  num_g2_b_d, &
                    num_g5_b_bb, num_g5_b_bd, num_g5_b_dd, &
                    num_g2_d_b,  num_g2_d_d, &
                    num_g5_d_bb, num_g5_d_bd, num_g5_d_dd, &
                    natom_b, natom_d, &
                    natom, rc, neighbor, io_sf, num_mpi, myrank )
       !C,D
       elseif( io_a == 0 .and. io_b == 0 .and. io_c == 1 .and. io_d == 1 )then
         if( myrank == 0 .and. imd == 1 ) write(*,*) " SF_2", nelem, elem_c, elem_d
         call sf_2( lattice, x, y, z, posi, xlattice, ylattice, zlattice, &
                    g2_c_c,  eta2_c_c,  rs2_c_c, &
                    g2_c_d,  eta2_c_d,  rs2_c_d, &
                    g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc, &
                    g5_c_cd, eta5_c_cd, theta5_c_cd, zeta5_c_cd, lambda5_c_cd, &
                    g5_c_dd, eta5_c_dd, theta5_c_dd, zeta5_c_dd, lambda5_c_dd, &
                    g2_d_c, eta2_d_c, rs2_d_c, &
                    g2_d_d, eta2_d_d, rs2_d_d, &
                    g5_d_cc, eta5_d_cc, theta5_d_cc, zeta5_d_cc, lambda5_d_cc, &
                    g5_d_cd, eta5_d_cd, theta5_d_cd, zeta5_d_cd, lambda5_d_cd, &
                    g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd, &
                    g2_deriv_c_c,  g2_deriv_c_d, &
                    g5_deriv_c_cc, g5_deriv_c_cd, g5_deriv_c_dd, &
                    g2_deriv_d_c,  g2_deriv_d_d, &
                    g5_deriv_d_cc, g5_deriv_d_cd, g5_deriv_d_dd, &
                    num_g2_c_c,  num_g2_c_d, &
                    num_g5_c_cc, num_g5_c_cd, num_g5_c_dd, &
                    num_g2_d_c,  num_g2_d_d, &
                    num_g5_d_cc, num_g5_d_cd, num_g5_d_dd, &
                    natom_c, natom_d, &
                    natom, rc, neighbor, io_sf, num_mpi, myrank )
       else
         write(*,*) "Error(main): sf_2"
         stop
       endif
     endif ! nelem == 2
     ! 3 elements
     if( nelem == 3 )then
       !A,B,C
       if( io_a == 1 .and. io_b == 1 .and. io_c == 1 .and. io_d == 0 )then
         if( myrank == 0 .and. imd == 1 ) write(*,*) " SF_3", nelem, elem_a, elem_b, elem_c
         call sf_3( lattice, x, y, z, posi, xlattice, ylattice, zlattice, &
                    g2_a_a,  eta2_a_a,  rs2_a_a, &
                    g2_a_b,  eta2_a_b,  rs2_a_b, &
                    g2_a_c,  eta2_a_c,  rs2_a_c, &
                    g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                    g5_a_ab, eta5_a_ab, theta5_a_ab, zeta5_a_ab, lambda5_a_ab, &
                    g5_a_ac, eta5_a_ac, theta5_a_ac, zeta5_a_ac, lambda5_a_ac, &
                    g5_a_bb, eta5_a_bb, theta5_a_bb, zeta5_a_bb, lambda5_a_bb, &
                    g5_a_bc, eta5_a_bc, theta5_a_bc, zeta5_a_bc, lambda5_a_bc, &
                    g5_a_cc, eta5_a_cc, theta5_a_cc, zeta5_a_cc, lambda5_a_cc, &
                    g2_b_a,  eta2_b_a,  rs2_b_a, &
                    g2_b_b,  eta2_b_b,  rs2_b_b, &
                    g2_b_c,  eta2_b_c,  rs2_b_c, &
                    g5_b_aa, eta5_b_aa, theta5_b_aa, zeta5_b_aa, lambda5_b_aa, &
                    g5_b_ab, eta5_b_ab, theta5_b_ab, zeta5_b_ab, lambda5_b_ab, &
                    g5_b_ac, eta5_b_ac, theta5_b_ac, zeta5_b_ac, lambda5_b_ac, &
                    g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                    g5_b_bc, eta5_b_bc, theta5_b_bc, zeta5_b_bc, lambda5_b_bc, &
                    g5_b_cc, eta5_b_cc, theta5_b_cc, zeta5_b_cc, lambda5_b_cc, &
                    g2_c_a,  eta2_c_a,  rs2_c_a, &
                    g2_c_b,  eta2_c_b,  rs2_c_b, &
                    g2_c_c,  eta2_c_c,  rs2_c_c, &
                    g5_c_aa, eta5_c_aa, theta5_c_aa, zeta5_c_aa, lambda5_c_aa, &
                    g5_c_ab, eta5_c_ab, theta5_c_ab, zeta5_c_ab, lambda5_c_ab, &
                    g5_c_ac, eta5_c_ac, theta5_c_ac, zeta5_c_ac, lambda5_c_ac, &
                    g5_c_bb, eta5_c_bb, theta5_c_bb, zeta5_c_bb, lambda5_c_bb, &
                    g5_c_bc, eta5_c_bc, theta5_c_bc, zeta5_c_bc, lambda5_c_bc, &
                    g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc, &
                    g2_deriv_a_a,  g2_deriv_a_b,  g2_deriv_a_c, &
                    g5_deriv_a_aa, g5_deriv_a_ab, g5_deriv_a_ac, &
                    g5_deriv_a_bb, g5_deriv_a_bc, &
                    g5_deriv_a_cc, &
                    g2_deriv_b_a,  g2_deriv_b_b,  g2_deriv_b_c, &
                    g5_deriv_b_aa, g5_deriv_b_ab, g5_deriv_b_ac, &
                    g5_deriv_b_bb, g5_deriv_b_bc, &
                    g5_deriv_b_cc, &
                    g2_deriv_c_a,  g2_deriv_c_b,  g2_deriv_c_c, &
                    g5_deriv_c_aa, g5_deriv_c_ab, g5_deriv_c_ac, &
                    g5_deriv_c_bb, g5_deriv_c_bc, &
                    g5_deriv_c_cc, &
                    num_g2_a_a,  num_g2_a_b,  num_g2_a_c, &
                    num_g5_a_aa, num_g5_a_ab, num_g5_a_ac, &
                    num_g5_a_bb, num_g5_a_bc, &
                    num_g5_a_cc, &
                    num_g2_b_a,  num_g2_b_b,  num_g2_b_c, &
                    num_g5_b_aa, num_g5_b_ab, num_g5_b_ac, &
                    num_g5_b_bb, num_g5_b_bc, &
                    num_g5_b_cc, &
                    num_g2_c_a,  num_g2_c_b,  num_g2_c_c, &
                    num_g5_c_aa, num_g5_c_ab, num_g5_c_ac, &
                    num_g5_c_bb, num_g5_c_bc, &
                    num_g5_c_cc, &
                    natom_a, natom_b, natom_c, &
                    natom, rc, neighbor, io_sf, num_mpi, myrank )
       !A,B,D
       elseif( io_a == 1 .and. io_b == 1 .and. io_c == 0 .and. io_d == 1 )then
         if( myrank == 0 .and. imd == 1 ) write(*,*) " SF_3", nelem, elem_a, elem_b, elem_d
         call sf_3( lattice, x, y, z, posi, xlattice, ylattice, zlattice, &
                    g2_a_a,  eta2_a_a,  rs2_a_a, &
                    g2_a_b,  eta2_a_b,  rs2_a_b, &
                    g2_a_d,  eta2_a_d,  rs2_a_d, &
                    g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                    g5_a_ab, eta5_a_ab, theta5_a_ab, zeta5_a_ab, lambda5_a_ab, &
                    g5_a_ad, eta5_a_ad, theta5_a_ad, zeta5_a_ad, lambda5_a_ad, &
                    g5_a_bb, eta5_a_bb, theta5_a_bb, zeta5_a_bb, lambda5_a_bb, &
                    g5_a_bd, eta5_a_bd, theta5_a_bd, zeta5_a_bd, lambda5_a_bd, &
                    g5_a_dd, eta5_a_dd, theta5_a_dd, zeta5_a_dd, lambda5_a_dd, &
                    g2_b_a,  eta2_b_a,  rs2_b_a, &
                    g2_b_b,  eta2_b_b,  rs2_b_b, &
                    g2_b_d,  eta2_b_d,  rs2_b_d, &
                    g5_b_aa, eta5_b_aa, theta5_b_aa, zeta5_b_aa, lambda5_b_aa, &
                    g5_b_ab, eta5_b_ab, theta5_b_ab, zeta5_b_ab, lambda5_b_ab, &
                    g5_b_ad, eta5_b_ad, theta5_b_ad, zeta5_b_ad, lambda5_b_ad, &
                    g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                    g5_b_bd, eta5_b_bd, theta5_b_bd, zeta5_b_bd, lambda5_b_bd, &
                    g5_b_dd, eta5_b_dd, theta5_b_dd, zeta5_b_dd, lambda5_b_dd, &
                    g2_d_a,  eta2_d_a,  rs2_d_a, &
                    g2_d_b,  eta2_d_b,  rs2_d_b, &
                    g2_d_d,  eta2_d_d,  rs2_d_d, &
                    g5_d_aa, eta5_d_aa, theta5_d_aa, zeta5_d_aa, lambda5_d_aa, &
                    g5_d_ab, eta5_d_ab, theta5_d_ab, zeta5_d_ab, lambda5_d_ab, &
                    g5_d_ad, eta5_d_ad, theta5_d_ad, zeta5_d_ad, lambda5_d_ad, &
                    g5_d_bb, eta5_d_bb, theta5_d_bb, zeta5_d_bb, lambda5_d_bb, &
                    g5_d_bd, eta5_d_bd, theta5_d_bd, zeta5_d_bd, lambda5_d_bd, &
                    g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd, &
                    g2_deriv_a_a,  g2_deriv_a_b,  g2_deriv_a_d, &
                    g5_deriv_a_aa, g5_deriv_a_ab, g5_deriv_a_ad, &
                    g5_deriv_a_bb, g5_deriv_a_bd, &
                    g5_deriv_a_dd, &
                    g2_deriv_b_a,  g2_deriv_b_b,  g2_deriv_b_d, &
                    g5_deriv_b_aa, g5_deriv_b_ab, g5_deriv_b_ad, &
                    g5_deriv_b_bb, g5_deriv_b_bd, &
                    g5_deriv_b_dd, &
                    g2_deriv_d_a,  g2_deriv_d_b,  g2_deriv_d_d, &
                    g5_deriv_d_aa, g5_deriv_d_ab, g5_deriv_d_ad, &
                    g5_deriv_d_bb, g5_deriv_d_bd, &
                    g5_deriv_d_dd, &
                    num_g2_a_a,  num_g2_a_b,  num_g2_a_d, &
                    num_g5_a_aa, num_g5_a_ab, num_g5_a_ad, &
                    num_g5_a_bb, num_g5_a_bd, &
                    num_g5_a_dd, &
                    num_g2_b_a,  num_g2_b_b,  num_g2_b_d, &
                    num_g5_b_aa, num_g5_b_ab, num_g5_b_ad, &
                    num_g5_b_bb, num_g5_b_bd, &
                    num_g5_b_dd, &
                    num_g2_d_a,  num_g2_d_b,  num_g2_d_d, &
                    num_g5_d_aa, num_g5_d_ab, num_g5_d_ad, &
                    num_g5_d_bb, num_g5_d_bd, &
                    num_g5_d_dd, &
                    natom_a, natom_b, natom_d, &
                    natom, rc, neighbor, io_sf, num_mpi, myrank )
       !A,C,D
       elseif( io_a == 1 .and. io_b == 0 .and. io_c == 1 .and. io_d == 1 )then
         if( myrank == 0 .and. imd == 1 ) write(*,*) " SF_3", nelem, elem_a, elem_c, elem_d
         call sf_3( lattice, x, y, z, posi, xlattice, ylattice, zlattice, &
                    g2_a_a,  eta2_a_a,  rs2_a_a, &
                    g2_a_c,  eta2_a_c,  rs2_a_c, &
                    g2_a_d,  eta2_a_d,  rs2_a_d, &
                    g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                    g5_a_ac, eta5_a_ac, theta5_a_ac, zeta5_a_ac, lambda5_a_ac, &
                    g5_a_ad, eta5_a_ad, theta5_a_ad, zeta5_a_ad, lambda5_a_ad, &
                    g5_a_cc, eta5_a_cc, theta5_a_cc, zeta5_a_cc, lambda5_a_cc, &
                    g5_a_cd, eta5_a_cd, theta5_a_cd, zeta5_a_cd, lambda5_a_cd, &
                    g5_a_dd, eta5_a_dd, theta5_a_dd, zeta5_a_dd, lambda5_a_dd, &
                    g2_c_a,  eta2_c_a,  rs2_b_a, &
                    g2_c_c,  eta2_c_c,  rs2_b_c, &
                    g2_c_d,  eta2_c_d,  rs2_b_d, &
                    g5_c_aa, eta5_c_aa, theta5_c_aa, zeta5_c_aa, lambda5_c_aa, &
                    g5_c_ac, eta5_c_ac, theta5_c_ac, zeta5_c_ac, lambda5_c_ac, &
                    g5_c_ad, eta5_c_ad, theta5_c_ad, zeta5_c_ad, lambda5_c_ad, &
                    g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc, &
                    g5_c_cd, eta5_c_cd, theta5_c_cd, zeta5_c_cd, lambda5_c_cd, &
                    g5_c_dd, eta5_c_dd, theta5_c_dd, zeta5_c_dd, lambda5_c_dd, &
                    g2_d_a,  eta2_d_a,  rs2_d_a, &
                    g2_d_c,  eta2_d_c,  rs2_d_c, &
                    g2_d_d,  eta2_d_d,  rs2_d_d, &
                    g5_d_aa, eta5_d_aa, theta5_d_aa, zeta5_d_aa, lambda5_d_aa, &
                    g5_d_ac, eta5_d_ac, theta5_d_ac, zeta5_d_ac, lambda5_d_ac, &
                    g5_d_ad, eta5_d_ad, theta5_d_ad, zeta5_d_ad, lambda5_d_ad, &
                    g5_d_cc, eta5_d_cc, theta5_d_cc, zeta5_d_cc, lambda5_d_cc, &
                    g5_d_cd, eta5_d_cd, theta5_d_cd, zeta5_d_cd, lambda5_d_cd, &
                    g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd, &
                    g2_deriv_a_a,  g2_deriv_a_c,  g2_deriv_a_d, &
                    g5_deriv_a_aa, g5_deriv_a_ac, g5_deriv_a_ad, &
                    g5_deriv_a_cc, g5_deriv_a_cd, &
                    g5_deriv_a_dd, &
                    g2_deriv_c_a,  g2_deriv_c_c,  g2_deriv_c_d, &
                    g5_deriv_c_aa, g5_deriv_c_ac, g5_deriv_c_ad, &
                    g5_deriv_c_cc, g5_deriv_c_cd, &
                    g5_deriv_c_dd, &
                    g2_deriv_d_a,  g2_deriv_d_c,  g2_deriv_d_d, &
                    g5_deriv_d_aa, g5_deriv_d_ac, g5_deriv_d_ad, &
                    g5_deriv_d_cc, g5_deriv_d_cd, &
                    g5_deriv_d_dd, &
                    num_g2_a_a,  num_g2_a_c,  num_g2_a_d, &
                    num_g5_a_aa, num_g5_a_ac, num_g5_a_ad, &
                    num_g5_a_cc, num_g5_a_cd, &
                    num_g5_a_dd, &
                    num_g2_c_a,  num_g2_c_c,  num_g2_c_d, &
                    num_g5_c_aa, num_g5_c_ac, num_g5_c_ad, &
                    num_g5_c_cc, num_g5_c_cd, &
                    num_g5_c_dd, &
                    num_g2_d_a,  num_g2_d_c,  num_g2_d_d, &
                    num_g5_d_aa, num_g5_d_ac, num_g5_d_ad, &
                    num_g5_d_cc, num_g5_d_cd, &
                    num_g5_d_dd, &
                    natom_a, natom_c, natom_d, &
                    natom, rc, neighbor, io_sf, num_mpi, myrank )
       !B,C,D
       elseif( io_a == 0 .and. io_b == 1 .and. io_c == 1 .and. io_d == 1 )then
         if( myrank == 0 .and. imd == 1 ) write(*,*) " SF_3", nelem, elem_b, elem_c, elem_d
         call sf_3( lattice, x, y, z, posi, xlattice, ylattice, zlattice, &
                    g2_b_b,  eta2_b_b,  rs2_b_b, &
                    g2_b_c,  eta2_b_c,  rs2_b_c, &
                    g2_b_d,  eta2_b_d,  rs2_b_d, &
                    g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                    g5_b_bc, eta5_b_bc, theta5_b_bc, zeta5_b_bc, lambda5_b_bc, &
                    g5_b_bd, eta5_b_bd, theta5_b_bd, zeta5_b_bd, lambda5_b_bd, &
                    g5_b_cc, eta5_b_cc, theta5_b_cc, zeta5_b_cc, lambda5_b_cc, &
                    g5_b_cd, eta5_b_cd, theta5_b_cd, zeta5_b_cd, lambda5_b_cd, &
                    g5_b_dd, eta5_b_dd, theta5_b_dd, zeta5_b_dd, lambda5_b_dd, &
                    g2_c_b,  eta2_c_b,  rs2_c_b, &
                    g2_c_c,  eta2_c_c,  rs2_c_c, &
                    g2_c_d,  eta2_c_d,  rs2_c_d, &
                    g5_c_bb, eta5_c_bb, theta5_c_bb, zeta5_c_bb, lambda5_c_bb, &
                    g5_c_bc, eta5_c_bc, theta5_c_bc, zeta5_c_bc, lambda5_c_bc, &
                    g5_c_bd, eta5_c_bd, theta5_c_bd, zeta5_c_bd, lambda5_c_bd, &
                    g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc, &
                    g5_c_cd, eta5_c_cd, theta5_c_cd, zeta5_c_cd, lambda5_c_cd, &
                    g5_c_dd, eta5_c_dd, theta5_c_dd, zeta5_c_dd, lambda5_c_dd, &
                    g2_d_b,  eta2_d_b,  rs2_d_b, &
                    g2_d_c,  eta2_d_c,  rs2_d_c, &
                    g2_d_d,  eta2_d_d,  rs2_d_d, &
                    g5_d_bb, eta5_d_bb, theta5_d_bb, zeta5_d_bb, lambda5_d_bb, &
                    g5_d_bc, eta5_d_bc, theta5_d_bc, zeta5_d_bc, lambda5_d_bc, &
                    g5_d_bd, eta5_d_bd, theta5_d_bd, zeta5_d_bd, lambda5_d_bd, &
                    g5_d_cc, eta5_d_cc, theta5_d_cc, zeta5_d_cc, lambda5_d_cc, &
                    g5_d_cd, eta5_d_cd, theta5_d_cd, zeta5_d_cd, lambda5_d_cd, &
                    g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd, &
                    g2_deriv_b_b,  g2_deriv_b_c,  g2_deriv_b_d, &
                    g5_deriv_b_bb, g5_deriv_b_bc, g5_deriv_b_bd, &
                    g5_deriv_b_cc, g5_deriv_b_cd, &
                    g5_deriv_b_dd, &
                    g2_deriv_c_b,  g2_deriv_c_c,  g2_deriv_c_d, &
                    g5_deriv_c_bb, g5_deriv_c_bc, g5_deriv_c_bd, &
                    g5_deriv_c_cc, g5_deriv_c_cd, &
                    g5_deriv_c_dd, &
                    g2_deriv_d_b,  g2_deriv_d_c,  g2_deriv_d_d, &
                    g5_deriv_d_bb, g5_deriv_d_bc, g5_deriv_d_bd, &
                    g5_deriv_d_cc, g5_deriv_d_cd, &
                    g5_deriv_d_dd, &
                    num_g2_b_b,  num_g2_b_c,  num_g2_b_d, &
                    num_g5_b_bb, num_g5_b_bc, num_g5_b_bd, &
                    num_g5_b_cc, num_g5_b_cd, &
                    num_g5_b_dd, &
                    num_g2_c_b,  num_g2_c_c,  num_g2_c_d, &
                    num_g5_c_bb, num_g5_c_bc, num_g5_c_bd, &
                    num_g5_c_cc, num_g5_c_cd, &
                    num_g5_c_dd, &
                    num_g2_d_b,  num_g2_d_c,  num_g2_d_d, &
                    num_g5_d_bb, num_g5_d_bc, num_g5_d_bd, &
                    num_g5_d_cc, num_g5_d_cd, &
                    num_g5_d_dd, &
                    natom_b, natom_c, natom_d, &
                    natom, rc, neighbor, io_sf, num_mpi, myrank )
       else
         write(*,*) "Error(main): sf_3"
         stop
       endif
     endif ! nelem == 3
     ! 4 elements
     if( nelem == 4 )then
       !A,B,C,D
       if( io_a == 1 .and. io_b == 1 .and. io_c == 1 .and. io_d == 1 )then
         if( myrank == 0 .and. imd == 1 ) write(*,*) " SF_4", nelem, elem_a, elem_b, elem_c, elem_d
         call sf_4( lattice, x, y, z, posi, xlattice, ylattice, zlattice, &
                    g2_a_a,  eta2_a_a,  rs2_a_a, &
                    g2_a_b,  eta2_a_b,  rs2_a_b, &
                    g2_a_c,  eta2_a_c,  rs2_a_c, &
                    g2_a_d,  eta2_a_d,  rs2_a_d, &
                    g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                    g5_a_ab, eta5_a_ab, theta5_a_ab, zeta5_a_ab, lambda5_a_ab, &
                    g5_a_ac, eta5_a_ac, theta5_a_ac, zeta5_a_ac, lambda5_a_ac, &
                    g5_a_ad, eta5_a_ad, theta5_a_ad, zeta5_a_ad, lambda5_a_ad, &
                    g5_a_bb, eta5_a_bb, theta5_a_bb, zeta5_a_bb, lambda5_a_bb, &
                    g5_a_bc, eta5_a_bc, theta5_a_bc, zeta5_a_bc, lambda5_a_bc, &
                    g5_a_bd, eta5_a_bd, theta5_a_bd, zeta5_a_bd, lambda5_a_bd, &
                    g5_a_cc, eta5_a_cc, theta5_a_cc, zeta5_a_cc, lambda5_a_cc, &
                    g5_a_cd, eta5_a_cd, theta5_a_cd, zeta5_a_cd, lambda5_a_cd, &
                    g5_a_dd, eta5_a_dd, theta5_a_dd, zeta5_a_dd, lambda5_a_dd, &
                    g2_b_a,  eta2_b_a,  rs2_b_a, &
                    g2_b_b,  eta2_b_b,  rs2_b_b, &
                    g2_b_c,  eta2_b_c,  rs2_b_c, &
                    g2_b_d,  eta2_b_d,  rs2_b_d, &
                    g5_b_aa, eta5_b_aa, theta5_b_aa, zeta5_b_aa, lambda5_b_aa, &
                    g5_b_ab, eta5_b_ab, theta5_b_ab, zeta5_b_ab, lambda5_b_ab, &
                    g5_b_ac, eta5_b_ac, theta5_b_ac, zeta5_b_ac, lambda5_b_ac, &
                    g5_b_ad, eta5_b_ad, theta5_b_ad, zeta5_b_ad, lambda5_b_ad, &
                    g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                    g5_b_bc, eta5_b_bc, theta5_b_bc, zeta5_b_bc, lambda5_b_bc, &
                    g5_b_bd, eta5_b_bd, theta5_b_bd, zeta5_b_bd, lambda5_b_bd, &
                    g5_b_cc, eta5_b_cc, theta5_b_cc, zeta5_b_cc, lambda5_b_cc, &
                    g5_b_cd, eta5_b_cd, theta5_b_cd, zeta5_b_cd, lambda5_b_cd, &
                    g5_b_dd, eta5_b_dd, theta5_b_dd, zeta5_b_dd, lambda5_b_dd, &
                    g2_c_a,  eta2_c_a,  rs2_c_a, &
                    g2_c_b,  eta2_c_b,  rs2_c_b, &
                    g2_c_c,  eta2_c_c,  rs2_c_c, &
                    g2_c_d,  eta2_c_d,  rs2_c_d, &
                    g5_c_aa, eta5_c_aa, theta5_c_aa, zeta5_c_aa, lambda5_c_aa, &
                    g5_c_ab, eta5_c_ab, theta5_c_ab, zeta5_c_ab, lambda5_c_ab, &
                    g5_c_ac, eta5_c_ac, theta5_c_ac, zeta5_c_ac, lambda5_c_ac, &
                    g5_c_ad, eta5_c_ad, theta5_c_ad, zeta5_c_ad, lambda5_c_ad, &
                    g5_c_bb, eta5_c_bb, theta5_c_bb, zeta5_c_bb, lambda5_c_bb, &
                    g5_c_bc, eta5_c_bc, theta5_c_bc, zeta5_c_bc, lambda5_c_bc, &
                    g5_c_bd, eta5_c_bd, theta5_c_bd, zeta5_c_bd, lambda5_c_bd, &
                    g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc, &
                    g5_c_cd, eta5_c_cd, theta5_c_cd, zeta5_c_cd, lambda5_c_cd, &
                    g5_c_dd, eta5_c_dd, theta5_c_dd, zeta5_c_dd, lambda5_c_dd, &
                    g2_d_a,  eta2_d_a,  rs2_d_a, &
                    g2_d_b,  eta2_d_b,  rs2_d_b, &
                    g2_d_c,  eta2_d_c,  rs2_d_c, &
                    g2_d_d,  eta2_d_d,  rs2_d_d, &
                    g5_d_aa, eta5_d_aa, theta5_d_aa, zeta5_d_aa, lambda5_d_aa, &
                    g5_d_ab, eta5_d_ab, theta5_d_ab, zeta5_d_ab, lambda5_d_ab, &
                    g5_d_ac, eta5_d_ac, theta5_d_ac, zeta5_d_ac, lambda5_d_ac, &
                    g5_d_ad, eta5_d_ad, theta5_d_ad, zeta5_d_ad, lambda5_d_ad, &
                    g5_d_bb, eta5_d_bb, theta5_d_bb, zeta5_d_bb, lambda5_d_bb, &
                    g5_d_bc, eta5_d_bc, theta5_d_bc, zeta5_d_bc, lambda5_d_bc, &
                    g5_d_bd, eta5_d_bd, theta5_d_bd, zeta5_d_bd, lambda5_d_bd, &
                    g5_d_cc, eta5_d_cc, theta5_d_cc, zeta5_d_cc, lambda5_d_cc, &
                    g5_d_cd, eta5_d_cd, theta5_d_cd, zeta5_d_cd, lambda5_d_cd, &
                    g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd, &
                    g2_deriv_a_a,  g2_deriv_a_b,  g2_deriv_a_c,  g2_deriv_a_d, &
                    g5_deriv_a_aa, g5_deriv_a_ab, g5_deriv_a_ac, g5_deriv_a_ad, &
                    g5_deriv_a_bb, g5_deriv_a_bc, g5_deriv_a_bd, &
                    g5_deriv_a_cc, g5_deriv_a_cd, &
                    g5_deriv_a_dd, &
                    g2_deriv_b_a,  g2_deriv_b_b,  g2_deriv_b_c,  g2_deriv_b_d, &
                    g5_deriv_b_aa, g5_deriv_b_ab, g5_deriv_b_ac, g5_deriv_b_ad, &
                    g5_deriv_b_bb, g5_deriv_b_bc, g5_deriv_b_bd, &
                    g5_deriv_b_cc, g5_deriv_b_cd, &
                    g5_deriv_b_dd, &
                    g2_deriv_c_a,  g2_deriv_c_b,  g2_deriv_c_c,  g2_deriv_c_d, &
                    g5_deriv_c_aa, g5_deriv_c_ab, g5_deriv_c_ac, g5_deriv_c_ad, &
                    g5_deriv_c_bb, g5_deriv_c_bc, g5_deriv_c_bd, &
                    g5_deriv_c_cc, g5_deriv_c_cd, &
                    g5_deriv_c_dd, &
                    g2_deriv_d_a,  g2_deriv_d_b,  g2_deriv_d_c,  g2_deriv_d_d, &
                    g5_deriv_d_aa, g5_deriv_d_ab, g5_deriv_d_ac, g5_deriv_d_ad, &
                    g5_deriv_d_bb, g5_deriv_d_bc, g5_deriv_d_bd, &
                    g5_deriv_d_cc, g5_deriv_d_cd, &
                    g5_deriv_d_dd, &
                    num_g2_a_a,  num_g2_a_b,  num_g2_a_c,  num_g2_a_d, &
                    num_g5_a_aa, num_g5_a_ab, num_g5_a_ac, num_g5_a_ad, &
                    num_g5_a_bb, num_g5_a_bc, num_g5_a_bd, &
                    num_g5_a_cc, num_g5_a_cd, &
                    num_g5_a_dd, &
                    num_g2_b_a,  num_g2_b_b,  num_g2_b_c,  num_g2_b_d, &
                    num_g5_b_aa, num_g5_b_ab, num_g5_b_ac, num_g5_b_ad, &
                    num_g5_b_bb, num_g5_b_bc, num_g5_b_bd, &
                    num_g5_b_cc, num_g5_b_cd, &
                    num_g5_b_dd, &
                    num_g2_c_a,  num_g2_c_b,  num_g2_c_c,  num_g2_c_d, &
                    num_g5_c_aa, num_g5_c_ab, num_g5_c_ac, num_g5_c_ad, &
                    num_g5_c_bb, num_g5_c_bc, num_g5_c_bd, &
                    num_g5_c_cc, num_g5_c_cd, &
                    num_g5_c_dd, &
                    num_g2_d_a,  num_g2_d_b,  num_g2_d_c,  num_g2_d_d, &
                    num_g5_d_aa, num_g5_d_ab, num_g5_d_ac, num_g5_d_ad, &
                    num_g5_d_bb, num_g5_d_bc, num_g5_d_bd, &
                    num_g5_d_cc, num_g5_d_cd, &
                    num_g5_d_dd, &
                    natom_a, natom_b, natom_c, natom_d, &
                    natom, rc, neighbor, io_sf, num_mpi, myrank )
       else
         write(*,*) "Error(main): sf_4"
         stop
       endif
     endif ! nelem == 4
     ! more than 5 elements
     if( nelem > 4 )then
       write(*,*) "Error(main): sf"
       stop
     endif ! nelem > 4


   enddo ! imd


   if( myrank == 0 )then
     write(*,*) " "
     write(*,*) " * Finish SF calculations *"
     write(*,*) " "
   endif ! myrank == 0


   close(1) ! info
   close(2) ! positions
   close(11) ! sf_
   close(12) ! sf_deriv_


   !*************************************************************
   ! Deallocate variables
   !*************************************************************
   ! 1 element
   if( nelem == 1 )then
     !A
     if( io_a == 1 .and. io_b == 0 .and. io_c == 0 .and. io_d == 0 )then
       deallocate( g2_a_a,  eta2_a_a,  rs2_a_a, &
                   g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa )
       deallocate( g2_deriv_a_a, g5_deriv_a_aa )
       if( myrank == 0 ) write(*,*) " Deallocate G2(a-a),G5(a-aa): ", elem_a
     !B
     elseif( io_a == 0 .and. io_b == 1 .and. io_c == 0 .and. io_d == 0 )then
       deallocate( g2_b_b,  eta2_b_b,  rs2_b_b, &
                   g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb )
       deallocate( g2_deriv_b_b, g5_deriv_b_bb )
       if( myrank == 0 ) write(*,*) " Deallocate G2(b-b),G5(b-bb): ", elem_b
     !C
     elseif( io_a == 0 .and. io_b == 0 .and. io_c == 1 .and. io_d == 0 )then
       deallocate( g2_c_c,  eta2_c_c,  rs2_c_c, &
                   g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc )
       deallocate( g2_deriv_c_c, g5_deriv_c_cc )
       if( myrank == 0 ) write(*,*) " Deallocate G2(c-c),G5(c-cc): ", elem_c
     !D
     elseif( io_a == 0 .and. io_b == 0 .and. io_c == 0 .and. io_d == 1 )then
       deallocate( g2_d_d,  eta2_d_d,  rs2_d_d, &
                   g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd )
       deallocate( g2_deriv_d_d, g5_deriv_d_dd )
       if( myrank == 0 ) write(*,*) " Deallocate G2(d-d),G5(d-dd): ", elem_d
     else
       write(*,*) "Error(main): deallocate_1"
       stop
     endif
   endif ! neleme == 1
   ! 2 elements
   if( nelem == 2 )then
     !A,B
     if( io_a == 1 .and. io_b == 1 .and. io_c == 0 .and. io_d == 0 )then
       deallocate( g2_a_a,  eta2_a_a,  rs2_a_a, &
                   g2_a_b,  eta2_a_b,  rs2_a_b, &
                   g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                   g5_a_ab, eta5_a_ab, theta5_a_ab, zeta5_a_ab, lambda5_a_ab, &
                   g5_a_bb, eta5_a_bb, theta5_a_bb, zeta5_a_bb, lambda5_a_bb )
       deallocate( g2_b_a,  eta2_b_a,  rs2_b_a, &
                   g2_b_b,  eta2_b_b,  rs2_b_b, &
                   g5_b_aa, eta5_b_aa, theta5_b_aa, zeta5_b_aa, lambda5_b_aa, &
                   g5_b_ab, eta5_b_ab, theta5_b_ab, zeta5_b_ab, lambda5_b_ab, &
                   g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb )
       deallocate( g2_deriv_a_a,  g2_deriv_a_b, &
                   g5_deriv_a_aa, g5_deriv_a_ab, g5_deriv_a_bb )
       deallocate( g2_deriv_b_a,  g2_deriv_b_b, &
                   g5_deriv_b_aa, g5_deriv_b_ab, g5_deriv_b_bb )
       if( myrank == 0 ) write(*,*) " Deallocate G2(a-a,a-b),G5(a-aa,a-ab,a-bb): ", elem_a
       if( myrank == 0 ) write(*,*) " Deallocate G2(b-a,b-b),G5(b-aa,b-ab,b-bb): ", elem_b
     !A,C
     elseif( io_a == 1 .and. io_b == 0 .and. io_c == 1 .and. io_d == 0 )then
       deallocate( g2_a_a,  eta2_a_a,  rs2_a_a, &
                   g2_a_c,  eta2_a_c,  rs2_a_c, &
                   g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                   g5_a_ac, eta5_a_ac, theta5_a_ac, zeta5_a_ac, lambda5_a_ac, &
                   g5_a_cc, eta5_a_cc, theta5_a_cc, zeta5_a_cc, lambda5_a_cc )
       deallocate( g2_c_a,  eta2_c_a,  rs2_c_a, &
                   g2_c_c,  eta2_c_c,  rs2_c_c, &
                   g5_c_aa, eta5_c_aa, theta5_c_aa, zeta5_c_aa, lambda5_c_aa, &
                   g5_c_ac, eta5_c_ac, theta5_c_ac, zeta5_c_ac, lambda5_c_ac, &
                   g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc )
       deallocate( g2_deriv_a_a,  g2_deriv_a_c, &
                   g5_deriv_a_aa, g5_deriv_a_ac, g5_deriv_a_cc )
       deallocate( g2_deriv_c_a,  g2_deriv_c_c, &
                   g5_deriv_c_aa, g5_deriv_c_ac, g5_deriv_c_cc )
       if( myrank == 0 ) write(*,*) " Deallocate G2(a-a,a-c),G5(a-aa,a-ac,a-cc): ", elem_a
       if( myrank == 0 ) write(*,*) " Deallocate G2(c-a,c-c),G5(c-aa,c-ac,c-cc): ", elem_c
     !A,D
     elseif( io_a == 1 .and. io_b == 0 .and. io_c == 0 .and. io_d == 1 )then
       deallocate( g2_a_a,  eta2_a_a,  rs2_a_a, &
                   g2_a_d,  eta2_a_d,  rs2_a_d, &
                   g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                   g5_a_ad, eta5_a_ad, theta5_a_ad, zeta5_a_ad, lambda5_a_ad, &
                   g5_a_dd, eta5_a_dd, theta5_a_dd, zeta5_a_dd, lambda5_a_dd )
       deallocate( g2_d_a,  eta2_d_a,  rs2_d_a, &
                   g2_d_d,  eta2_d_d,  rs2_d_d, &
                   g5_d_aa, eta5_d_aa, theta5_d_aa, zeta5_d_aa, lambda5_d_aa, &
                   g5_d_ad, eta5_d_ad, theta5_d_ad, zeta5_d_ad, lambda5_d_ad, &
                   g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd )
       deallocate( g2_deriv_a_a,  g2_deriv_a_d, &
                   g5_deriv_a_aa, g5_deriv_a_ad, g5_deriv_a_dd )
       deallocate( g2_deriv_d_a,  g2_deriv_d_d, &
                   g5_deriv_d_aa, g5_deriv_d_ad, g5_deriv_d_dd )
       if( myrank == 0 ) write(*,*) " Deallocate G2(a-a,a-d),G5(a-aa,a-ad,a-dd): ", elem_a
       if( myrank == 0 ) write(*,*) " Deallocate G2(d-a,d-d),G5(d-aa,d-ad,d-dd): ", elem_d
     !B,C
     elseif( io_a == 0 .and. io_b == 1 .and. io_c == 1 .and. io_d == 0 )then
       deallocate( g2_b_b,  eta2_b_b,  rs2_b_b, &
                   g2_b_c,  eta2_b_c,  rs2_b_c, &
                   g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                   g5_b_bc, eta5_b_bc, theta5_b_bc, zeta5_b_bc, lambda5_b_bc, &
                   g5_b_cc, eta5_b_cc, theta5_b_cc, zeta5_b_cc, lambda5_b_cc )
       deallocate( g2_c_b,  eta2_c_b,  rs2_c_b, &
                   g2_c_c,  eta2_c_c,  rs2_c_c, &
                   g5_c_bb, eta5_c_bb, theta5_c_bb, zeta5_c_bb, lambda5_c_bb, &
                   g5_c_bc, eta5_c_bc, theta5_c_bc, zeta5_c_bc, lambda5_c_bc, &
                   g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc )
       deallocate( g2_deriv_b_b,  g2_deriv_b_c, &
                   g5_deriv_b_bb, g5_deriv_b_bc, g5_deriv_b_cc )
       deallocate( g2_deriv_c_b,  g2_deriv_c_c, &
                   g5_deriv_c_bb, g5_deriv_c_bc, g5_deriv_c_cc )
       if( myrank == 0 ) write(*,*) " Deallocate G2(b-b,b-c),G5(b-bb,b-bc,b-cc): ", elem_b
       if( myrank == 0 ) write(*,*) " Deallocate G2(c-b,c-c),G5(c-bb,c-bc,c-cc): ", elem_c
     !B,D
     elseif( io_a == 0 .and. io_b == 1 .and. io_c == 0 .and. io_d == 1 )then
       deallocate( g2_b_b,  eta2_b_b,  rs2_b_b, &
                   g2_b_d,  eta2_b_d,  rs2_b_d, &
                   g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                   g5_b_bd, eta5_b_bd, theta5_b_bd, zeta5_b_bd, lambda5_b_bd, &
                   g5_b_dd, eta5_b_dd, theta5_b_dd, zeta5_b_dd, lambda5_b_dd )
       deallocate( g2_d_b,  eta2_d_b,  rs2_d_b, &
                   g2_d_d,  eta2_d_d,  rs2_d_d, &
                   g5_d_bb, eta5_d_bb, theta5_d_bb, zeta5_d_bb, lambda5_d_bb, &
                   g5_d_bd, eta5_d_bd, theta5_d_bd, zeta5_d_bd, lambda5_d_bd, &
                   g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd )
       deallocate( g2_deriv_b_b,  g2_deriv_b_d, &
                   g5_deriv_b_bb, g5_deriv_b_bd, g5_deriv_b_dd )
       deallocate( g2_deriv_d_b,  g2_deriv_d_d, &
                   g5_deriv_d_bb, g5_deriv_d_bd, g5_deriv_d_dd )
       if( myrank == 0 ) write(*,*) " Deallocate G2(b-b,b-d),G5(b-bb,b-bd,b-dd): ", elem_b
       if( myrank == 0 ) write(*,*) " Deallocate G2(d-b,d-d),G5(d-bb,d-bd,d-dd): ", elem_d
     !C,D
     elseif( io_a == 0 .and. io_b == 0 .and. io_c == 1 .and. io_d == 1 )then

       deallocate( g2_c_c,  eta2_c_c,  rs2_c_c, &
                   g2_c_d,  eta2_c_d,  rs2_c_d, &
                   g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc, &
                   g5_c_cd, eta5_c_cd, theta5_c_cd, zeta5_c_cd, lambda5_c_cd, &
                   g5_c_dd, eta5_c_dd, theta5_c_dd, zeta5_c_dd, lambda5_b_dd )
       deallocate( g2_d_c,  eta2_d_c,  rs2_d_c, &
                   g2_d_d,  eta2_d_d,  rs2_d_d, &
                   g5_d_cc, eta5_d_cc, theta5_d_cc, zeta5_d_cc, lambda5_d_cc, &
                   g5_d_cd, eta5_d_cd, theta5_d_cd, zeta5_d_cd, lambda5_d_cd, &
                   g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd )
       deallocate( g2_deriv_c_c,  g2_deriv_c_d, &
                   g5_deriv_c_cc, g5_deriv_c_cd, g5_deriv_c_dd )
       deallocate( g2_deriv_d_c,  g2_deriv_d_d, &
                   g5_deriv_d_cc, g5_deriv_d_cd, g5_deriv_d_dd )
       if( myrank == 0 ) write(*,*) " Deallocate G2(c-c,c-d),G5(c-cc,c-cd,c-dd): ", elem_c
       if( myrank == 0 ) write(*,*) " Deallocate G2(d-c,d-d),G5(d-cc,d-cd,d-dd): ", elem_d
     else
       write(*,*) "Error(main): deallocate_2"
       stop
     endif
   endif ! nelem = 2
   ! 3 elements
   if( nelem == 3 )then
     !A,B,C
     if( io_a == 1 .and. io_b == 1 .and. io_c == 1 .and. io_d == 0 )then
       deallocate( g2_a_a,  eta2_a_a,  rs2_a_a, &
                   g2_a_b,  eta2_a_b,  rs2_a_b, &
                   g2_a_c,  eta2_a_c,  rs2_a_c, &
                   g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                   g5_a_ab, eta5_a_ab, theta5_a_ab, zeta5_a_ab, lambda5_a_ab, &
                   g5_a_ac, eta5_a_ac, theta5_a_ac, zeta5_a_ac, lambda5_a_ac, &
                   g5_a_bb, eta5_a_bb, theta5_a_bb, zeta5_a_bb, lambda5_a_bb, &
                   g5_a_bc, eta5_a_bc, theta5_a_bc, zeta5_a_bc, lambda5_a_bc, &
                   g5_a_cc, eta5_a_cc, theta5_a_cc, zeta5_a_cc, lambda5_a_cc )
       deallocate( g2_b_a,  eta2_b_a,  rs2_b_a, &
                   g2_b_b,  eta2_b_b,  rs2_b_b, &
                   g2_b_c,  eta2_b_c,  rs2_b_c, &
                   g5_b_aa, eta5_b_aa, theta5_b_aa, zeta5_b_aa, lambda5_b_aa, &
                   g5_b_ab, eta5_b_ab, theta5_b_ab, zeta5_b_ab, lambda5_b_ab, &
                   g5_b_ac, eta5_b_ac, theta5_b_ac, zeta5_b_ac, lambda5_b_ac, &
                   g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                   g5_b_bc, eta5_b_bc, theta5_b_bc, zeta5_b_bc, lambda5_b_bc, &
                   g5_b_cc, eta5_b_cc, theta5_b_cc, zeta5_b_cc, lambda5_b_cc )
       deallocate( g2_c_a,  eta2_c_a,  rs2_c_a, &
                   g2_c_b,  eta2_c_b,  rs2_c_b, &
                   g2_c_c,  eta2_c_c,  rs2_c_c, &
                   g5_c_aa, eta5_c_aa, theta5_c_aa, zeta5_c_aa, lambda5_c_aa, &
                   g5_c_ab, eta5_c_ab, theta5_c_ab, zeta5_c_ab, lambda5_c_ab, &
                   g5_c_ac, eta5_c_ac, theta5_c_ac, zeta5_c_ac, lambda5_c_ac, &
                   g5_c_bb, eta5_c_bb, theta5_c_bb, zeta5_c_bb, lambda5_c_bb, &
                   g5_c_bc, eta5_c_bc, theta5_c_bc, zeta5_c_bc, lambda5_c_bc, &
                   g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc )
       deallocate( g2_deriv_a_a,  g2_deriv_a_b,  g2_deriv_a_c, &
                   g5_deriv_a_aa, g5_deriv_a_ab, g5_deriv_a_ac, &
                   g5_deriv_a_bb, g5_deriv_a_bc, g5_deriv_a_cc )
       deallocate( g2_deriv_b_a,  g2_deriv_b_b,  g2_deriv_b_c, &
                   g5_deriv_b_aa, g5_deriv_b_ab, g5_deriv_b_ac, &
                   g5_deriv_b_bb, g5_deriv_b_bc, g5_deriv_b_cc )
       deallocate( g2_deriv_c_a,  g2_deriv_c_b,  g2_deriv_c_c, &
                   g5_deriv_c_aa, g5_deriv_c_ab, g5_deriv_c_ac, &
                   g5_deriv_c_bb, g5_deriv_c_bc, g5_deriv_c_cc )
       if( myrank == 0 ) write(*,*) " Deallocate G2(a-a,a-b,a-c),G5(a-aa,a-ab,a-ac,a-bb,a-bc,a-cc): ", elem_a
       if( myrank == 0 ) write(*,*) " Deallocate G2(b-a,b-b,b-c),G5(b-aa,b-ab,b-ac,b-bb,b-bc,b-cc): ", elem_b
       if( myrank == 0 ) write(*,*) " Deallocate G2(c-a,c-b,c-c),G5(c-aa,c-ab,c-ac,c-bb,c-bc,c-cc): ", elem_c
     !A,B,D
     elseif( io_a == 1 .and. io_b == 1 .and. io_c == 0 .and. io_d == 1 )then
       deallocate( g2_a_a,  eta2_a_a,  rs2_a_a, &
                   g2_a_b,  eta2_a_b,  rs2_a_b, &
                   g2_a_d,  eta2_a_d,  rs2_a_d, &
                   g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                   g5_a_ab, eta5_a_ab, theta5_a_ab, zeta5_a_ab, lambda5_a_ab, &
                   g5_a_ad, eta5_a_ad, theta5_a_ad, zeta5_a_ad, lambda5_a_ad, &
                   g5_a_bb, eta5_a_bb, theta5_a_bb, zeta5_a_bb, lambda5_a_bb, &
                   g5_a_bd, eta5_a_bd, theta5_a_bd, zeta5_a_bd, lambda5_a_bd, &
                   g5_a_dd, eta5_a_dd, theta5_a_dd, zeta5_a_dd, lambda5_a_dd )
       deallocate( g2_b_a,  eta2_b_a,  rs2_b_a, &
                   g2_b_b,  eta2_b_b,  rs2_b_b, &
                   g2_b_d,  eta2_b_d,  rs2_b_d, &
                   g5_b_aa, eta5_b_aa, theta5_b_aa, zeta5_b_aa, lambda5_b_aa, &
                   g5_b_ab, eta5_b_ab, theta5_b_ab, zeta5_b_ab, lambda5_b_ab, &
                   g5_b_ad, eta5_b_ad, theta5_b_ad, zeta5_b_ad, lambda5_b_ad, &
                   g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                   g5_b_bd, eta5_b_bd, theta5_b_bd, zeta5_b_bd, lambda5_b_bd, &
                   g5_b_dd, eta5_b_dd, theta5_b_dd, zeta5_b_dd, lambda5_b_dd )
       deallocate( g2_d_a,  eta2_d_a,  rs2_d_a, &
                   g2_d_b,  eta2_d_b,  rs2_d_b, &
                   g2_d_d,  eta2_d_d,  rs2_d_d, &
                   g5_d_aa, eta5_d_aa, theta5_d_aa, zeta5_d_aa, lambda5_d_aa, &
                   g5_d_ab, eta5_d_ab, theta5_d_ab, zeta5_d_ab, lambda5_d_ab, &
                   g5_d_ad, eta5_d_ad, theta5_d_ad, zeta5_d_ad, lambda5_d_ad, &
                   g5_d_bb, eta5_d_bb, theta5_d_bb, zeta5_d_bb, lambda5_d_bb, &
                   g5_d_bd, eta5_d_bd, theta5_d_bd, zeta5_d_bd, lambda5_d_bd, &
                   g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd )
       deallocate( g2_deriv_a_a,  g2_deriv_a_b,  g2_deriv_a_d, &
                   g5_deriv_a_aa, g5_deriv_a_ab, g5_deriv_a_ad, &
                   g5_deriv_a_bb, g5_deriv_a_bd, g5_deriv_a_dd )
       deallocate( g2_deriv_b_a,  g2_deriv_b_b,  g2_deriv_b_d, &
                   g5_deriv_b_aa, g5_deriv_b_ab, g5_deriv_b_ad, &
                   g5_deriv_b_bb, g5_deriv_b_bd, g5_deriv_b_dd )
       deallocate( g2_deriv_d_a,  g2_deriv_d_b,  g2_deriv_d_d, &
                   g5_deriv_d_aa, g5_deriv_d_ab, g5_deriv_d_ad, &
                   g5_deriv_d_bb, g5_deriv_d_bd, g5_deriv_d_dd )
       if( myrank == 0 ) write(*,*) " Deallocate G2(a-a,a-b,a-d),G5(a-aa,a-ab,a-ad,a-bb,a-bd,a-dd): ", elem_a
       if( myrank == 0 ) write(*,*) " Deallocate G2(b-a,b-b,b-d),G5(b-aa,b-ab,b-ad,b-bb,b-bd,b-dd): ", elem_b
       if( myrank == 0 ) write(*,*) " Deallocate G2(d-a,d-b,d-d),G5(d-aa,d-ab,d-ad,d-bb,d-bd,d-dd): ", elem_d
     !A,C,D
     elseif( io_a == 1 .and. io_b == 0 .and. io_c == 1 .and. io_d == 1 )then
       deallocate( g2_a_a,  eta2_a_a,  rs2_a_a, &
                   g2_a_c,  eta2_a_c,  rs2_a_c, &
                   g2_a_d,  eta2_a_d,  rs2_a_d, &
                   g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                   g5_a_ac, eta5_a_ac, theta5_a_ac, zeta5_a_ac, lambda5_a_ac, &
                   g5_a_ad, eta5_a_ad, theta5_a_ad, zeta5_a_ad, lambda5_a_ad, &
                   g5_a_cc, eta5_a_cc, theta5_a_cc, zeta5_a_cc, lambda5_a_cc, &
                   g5_a_cd, eta5_a_cd, theta5_a_cd, zeta5_a_cd, lambda5_a_cd, &
                   g5_a_dd, eta5_a_dd, theta5_a_dd, zeta5_a_dd, lambda5_a_dd )
       deallocate( g2_c_a,  eta2_c_a,  rs2_c_a, &
                   g2_c_c,  eta2_c_c,  rs2_c_c, &
                   g2_c_d,  eta2_c_d,  rs2_c_d, &
                   g5_c_aa, eta5_c_aa, theta5_c_aa, zeta5_c_aa, lambda5_c_aa, &
                   g5_c_ac, eta5_c_ac, theta5_c_ac, zeta5_c_ac, lambda5_c_ac, &
                   g5_c_ad, eta5_c_ad, theta5_c_ad, zeta5_c_ad, lambda5_c_ad, &
                   g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc, &
                   g5_c_cd, eta5_c_cd, theta5_c_cd, zeta5_c_cd, lambda5_c_cd, &
                   g5_c_dd, eta5_c_dd, theta5_c_dd, zeta5_c_dd, lambda5_c_dd )
       deallocate( g2_d_a,  eta2_d_a,  rs2_d_a, &
                   g2_d_c,  eta2_d_c,  rs2_d_c, &
                   g2_d_d,  eta2_d_d,  rs2_d_d, &
                   g5_d_aa, eta5_d_aa, theta5_d_aa, zeta5_d_aa, lambda5_d_aa, &
                   g5_d_ac, eta5_d_ac, theta5_d_ac, zeta5_d_ac, lambda5_d_ac, &
                   g5_d_ad, eta5_d_ad, theta5_d_ad, zeta5_d_ad, lambda5_d_ad, &
                   g5_d_cc, eta5_d_cc, theta5_d_cc, zeta5_d_cc, lambda5_d_cc, &
                   g5_d_cd, eta5_d_cd, theta5_d_cd, zeta5_d_cd, lambda5_d_cd, &
                   g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd )
       deallocate( g2_deriv_a_a,  g2_deriv_a_c,  g2_deriv_a_d, &
                   g5_deriv_a_aa, g5_deriv_a_ac, g5_deriv_a_ad, &
                   g5_deriv_a_cc, g5_deriv_a_cd, g5_deriv_a_dd )
       deallocate( g2_deriv_c_a,  g2_deriv_c_c,  g2_deriv_c_d, &
                   g5_deriv_c_aa, g5_deriv_c_ac, g5_deriv_c_ad, &
                   g5_deriv_c_cc, g5_deriv_c_cd, g5_deriv_c_dd )
       deallocate( g2_deriv_d_a,  g2_deriv_d_c,  g2_deriv_d_d, &
                   g5_deriv_d_aa, g5_deriv_d_ac, g5_deriv_d_ad, &
                   g5_deriv_d_cc, g5_deriv_d_cd, g5_deriv_d_dd )
       if( myrank == 0 ) write(*,*) " Deallocate G2(a-a,a-c,a-d),G5(a-aa,a-ac,a-ad,a-cc,a-cd,a-dd): ", elem_a
       if( myrank == 0 ) write(*,*) " Deallocate G2(c-a,c-c,c-d),G5(c-aa,c-ac,c-ad,c-cc,c-cd,c-dd): ", elem_c
       if( myrank == 0 ) write(*,*) " Deallocate G2(d-a,d-c,d-d),G5(d-aa,d-ac,d-ad,d-cc,d-cd,d-dd): ", elem_d
     !B,C,D
     elseif( io_a == 0 .and. io_b == 1 .and. io_c == 1 .and. io_d == 1 )then
       deallocate( g2_b_b,  eta2_b_b,  rs2_b_b, &
                   g2_b_c,  eta2_b_c,  rs2_b_c, &
                   g2_b_d,  eta2_b_d,  rs2_b_d, &
                   g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                   g5_b_bc, eta5_b_bc, theta5_b_bc, zeta5_b_bc, lambda5_b_bc, &
                   g5_b_bd, eta5_b_bd, theta5_b_bd, zeta5_b_bd, lambda5_b_bd, &
                   g5_b_cc, eta5_b_cc, theta5_b_cc, zeta5_b_cc, lambda5_b_cc, &
                   g5_b_cd, eta5_b_cd, theta5_b_cd, zeta5_b_cd, lambda5_b_cd, &
                   g5_b_dd, eta5_b_dd, theta5_b_dd, zeta5_b_dd, lambda5_b_dd )
       deallocate( g2_c_b,  eta2_c_b,  rs2_c_b, &
                   g2_c_c,  eta2_c_c,  rs2_c_c, &
                   g2_c_d,  eta2_c_d,  rs2_c_d, &
                   g5_c_bb, eta5_c_bb, theta5_c_bb, zeta5_c_bb, lambda5_c_bb, &
                   g5_c_bc, eta5_c_bc, theta5_c_bc, zeta5_c_bc, lambda5_c_bc, &
                   g5_c_bd, eta5_c_bd, theta5_c_bd, zeta5_c_bd, lambda5_c_bd, &
                   g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc, &
                   g5_c_cd, eta5_c_cd, theta5_c_cd, zeta5_c_cd, lambda5_c_cd, &
                   g5_c_dd, eta5_c_dd, theta5_c_dd, zeta5_c_dd, lambda5_c_dd )
       deallocate( g2_d_b,  eta2_d_b,  rs2_d_b, &
                   g2_d_c,  eta2_d_c,  rs2_d_c, &
                   g2_d_d,  eta2_d_d,  rs2_d_d, &
                   g5_d_bb, eta5_d_bb, theta5_d_bb, zeta5_d_bb, lambda5_d_bb, &
                   g5_d_bc, eta5_d_bc, theta5_d_bc, zeta5_d_bc, lambda5_d_bc, &
                   g5_d_bd, eta5_d_bd, theta5_d_bd, zeta5_d_bd, lambda5_d_bd, &
                   g5_d_cc, eta5_d_cc, theta5_d_cc, zeta5_d_cc, lambda5_d_cc, &
                   g5_d_cd, eta5_d_cd, theta5_d_cd, zeta5_d_cd, lambda5_d_cd, &
                   g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd )
       deallocate( g2_deriv_b_b,  g2_deriv_b_c,  g2_deriv_b_d, &
                   g5_deriv_b_bb, g5_deriv_b_bc, g5_deriv_b_bd, &
                   g5_deriv_b_cc, g5_deriv_b_cd, g5_deriv_b_dd )
       deallocate( g2_deriv_c_b,  g2_deriv_c_c,  g2_deriv_c_d, &
                   g5_deriv_c_bb, g5_deriv_c_bc, g5_deriv_c_bd, &
                   g5_deriv_c_cc, g5_deriv_c_cd, g5_deriv_c_dd )
       deallocate( g2_deriv_d_b,  g2_deriv_d_c,  g2_deriv_d_d, &
                   g5_deriv_d_bb, g5_deriv_d_bc, g5_deriv_d_bd, &
                   g5_deriv_d_cc, g5_deriv_d_cd, g5_deriv_d_dd )
       if( myrank == 0 ) write(*,*) " Deallocate G2(b-b,b-c,b-d),G5(b-bb,b-bc,b-bd,b-cc,b-cd,b-dd): ", elem_b
       if( myrank == 0 ) write(*,*) " Deallocate G2(c-b,c-c,c-d),G5(c-bb,c-bc,c-bd,c-cc,c-cd,c-dd): ", elem_c
       if( myrank == 0 ) write(*,*) " Deallocate G2(d-b,d-c,d-d),G5(d-bb,d-bc,d-bd,d-cc,d-cd,d-dd): ", elem_d
     else
       write(*,*) "Error(main): deallocate_3"
       stop
     endif
   endif ! nelem == 3
   ! 4 elements
   if( nelem == 4 )then
     !A,B,C,D
     if( io_a == 1 .and. io_b == 1 .and. io_c == 1 .and. io_d == 1 )then
       deallocate( g2_a_a,  eta2_a_a,  rs2_a_a, &
                   g2_a_b,  eta2_a_b,  rs2_a_b, &
                   g2_a_c,  eta2_a_c,  rs2_a_c, &
                   g2_a_d,  eta2_a_d,  rs2_a_d, &
                   g5_a_aa, eta5_a_aa, theta5_a_aa, zeta5_a_aa, lambda5_a_aa, &
                   g5_a_ab, eta5_a_ab, theta5_a_ab, zeta5_a_ab, lambda5_a_ab, &
                   g5_a_ac, eta5_a_ac, theta5_a_ac, zeta5_a_ac, lambda5_a_ac, &
                   g5_a_ad, eta5_a_ad, theta5_a_ad, zeta5_a_ad, lambda5_a_ad, &
                   g5_a_bb, eta5_a_bb, theta5_a_bb, zeta5_a_bb, lambda5_a_bb, &
                   g5_a_bc, eta5_a_bc, theta5_a_bc, zeta5_a_bc, lambda5_a_bc, &
                   g5_a_bd, eta5_a_bd, theta5_a_bd, zeta5_a_bd, lambda5_a_bd, &
                   g5_a_cc, eta5_a_cc, theta5_a_cc, zeta5_a_cc, lambda5_a_cc, &
                   g5_a_cd, eta5_a_cd, theta5_a_cd, zeta5_a_cd, lambda5_a_cd, &
                   g5_a_dd, eta5_a_dd, theta5_a_dd, zeta5_a_dd, lambda5_a_dd )
       deallocate( g2_b_a,  eta2_b_a,  rs2_b_a, &
                   g2_b_b,  eta2_b_b,  rs2_b_b, &
                   g2_b_c,  eta2_b_c,  rs2_b_c, &
                   g2_b_d,  eta2_b_d,  rs2_b_d, &
                   g5_b_aa, eta5_b_aa, theta5_b_aa, zeta5_b_aa, lambda5_b_aa, &
                   g5_b_ab, eta5_b_ab, theta5_b_ab, zeta5_b_ab, lambda5_b_ab, &
                   g5_b_ac, eta5_b_ac, theta5_b_ac, zeta5_b_ac, lambda5_b_ac, &
                   g5_b_ad, eta5_b_ad, theta5_b_ad, zeta5_b_ad, lambda5_b_ad, &
                   g5_b_bb, eta5_b_bb, theta5_b_bb, zeta5_b_bb, lambda5_b_bb, &
                   g5_b_bc, eta5_b_bc, theta5_b_bc, zeta5_b_bc, lambda5_b_bc, &
                   g5_b_bd, eta5_b_bd, theta5_b_bd, zeta5_b_bd, lambda5_b_bd, &
                   g5_b_cc, eta5_b_cc, theta5_b_cc, zeta5_b_cc, lambda5_b_cc, &
                   g5_b_cd, eta5_b_cd, theta5_b_cd, zeta5_b_cd, lambda5_b_cd, &
                   g5_b_dd, eta5_b_dd, theta5_b_dd, zeta5_b_dd, lambda5_b_dd )
       deallocate( g2_c_a,  eta2_c_a,  rs2_c_a, &
                   g2_c_b,  eta2_c_b,  rs2_c_b, &
                   g2_c_c,  eta2_c_c,  rs2_c_c, &
                   g2_c_d,  eta2_c_d,  rs2_c_d, &
                   g5_c_aa, eta5_c_aa, theta5_c_aa, zeta5_c_aa, lambda5_c_aa, &
                   g5_c_ab, eta5_c_ab, theta5_c_ab, zeta5_c_ab, lambda5_c_ab, &
                   g5_c_ac, eta5_c_ac, theta5_c_ac, zeta5_c_ac, lambda5_c_ac, &
                   g5_c_ad, eta5_c_ad, theta5_c_ad, zeta5_c_ad, lambda5_c_ad, &
                   g5_c_bb, eta5_c_bb, theta5_c_bb, zeta5_c_bb, lambda5_c_bb, &
                   g5_c_bc, eta5_c_bc, theta5_c_bc, zeta5_c_bc, lambda5_c_bc, &
                   g5_c_bd, eta5_c_bd, theta5_c_bd, zeta5_c_bd, lambda5_c_bd, &
                   g5_c_cc, eta5_c_cc, theta5_c_cc, zeta5_c_cc, lambda5_c_cc, &
                   g5_c_cd, eta5_c_cd, theta5_c_cd, zeta5_c_cd, lambda5_c_cd, &
                   g5_c_dd, eta5_c_dd, theta5_c_dd, zeta5_c_dd, lambda5_c_dd )
       deallocate( g2_d_a,  eta2_d_a,  rs2_d_a, &
                   g2_d_b,  eta2_d_b,  rs2_d_b, &
                   g2_d_c,  eta2_d_c,  rs2_d_c, &
                   g2_d_d,  eta2_d_d,  rs2_d_d, &
                   g5_d_aa, eta5_d_aa, theta5_d_aa, zeta5_d_aa, lambda5_d_aa, &
                   g5_d_ab, eta5_d_ab, theta5_d_ab, zeta5_d_ab, lambda5_d_ab, &
                   g5_d_ac, eta5_d_ac, theta5_d_ac, zeta5_d_ac, lambda5_d_ac, &
                   g5_d_ad, eta5_d_ad, theta5_d_ad, zeta5_d_ad, lambda5_d_ad, &
                   g5_d_bb, eta5_d_bb, theta5_d_bb, zeta5_d_bb, lambda5_d_bb, &
                   g5_d_bc, eta5_d_bc, theta5_d_bc, zeta5_d_bc, lambda5_d_bc, &
                   g5_d_bd, eta5_d_bd, theta5_d_bd, zeta5_d_bd, lambda5_d_bd, &
                   g5_d_cc, eta5_d_cc, theta5_d_cc, zeta5_d_cc, lambda5_d_cc, &
                   g5_d_cd, eta5_d_cd, theta5_d_cd, zeta5_d_cd, lambda5_d_cd, &
                   g5_d_dd, eta5_d_dd, theta5_d_dd, zeta5_d_dd, lambda5_d_dd )
       deallocate( g2_deriv_a_a,  g2_deriv_a_b,  g2_deriv_a_c,  g2_deriv_a_d, &
                   g5_deriv_a_aa, g5_deriv_a_ab, g5_deriv_a_ac, g5_deriv_a_ad, &
                   g5_deriv_a_bb, g5_deriv_a_bc, g5_deriv_a_bd, &
                   g5_deriv_a_cc, g5_deriv_a_cd, &
                   g5_deriv_a_dd )
       deallocate( g2_deriv_b_a,  g2_deriv_b_b,  g2_deriv_b_c,  g2_deriv_b_d, &
                   g5_deriv_b_aa, g5_deriv_b_ab, g5_deriv_b_ac, g5_deriv_b_ad, &
                   g5_deriv_b_bb, g5_deriv_b_bc, g5_deriv_b_bd, &
                   g5_deriv_b_cc, g5_deriv_b_cd, &
                   g5_deriv_b_dd )
       deallocate( g2_deriv_c_a,  g2_deriv_c_b,  g2_deriv_c_c,  g2_deriv_c_d, &
                   g5_deriv_c_aa, g5_deriv_c_ab, g5_deriv_c_ac, g5_deriv_c_ad, &
                   g5_deriv_c_bb, g5_deriv_c_bc, g5_deriv_c_bd, &
                   g5_deriv_c_cc, g5_deriv_c_cd, &
                   g5_deriv_c_dd )
       deallocate( g2_deriv_d_a,  g2_deriv_d_b,  g2_deriv_d_c,  g2_deriv_d_d, &
                   g5_deriv_d_aa, g5_deriv_d_ab, g5_deriv_d_ac, g5_deriv_d_ad, &
                   g5_deriv_d_bb, g5_deriv_d_bc, g5_deriv_d_bd, &
                   g5_deriv_d_cc, g5_deriv_d_cd, &
                   g5_deriv_d_dd )
       if( myrank == 0 ) write(*,*) " Deallocate &
                   &G2(a-a,a-b,a-c,a-d),&
                   &G5(a-aa,a-ab,a-ac,a-ad,a-bb,a-bc,a-bd,a-cc,a-cd,a-dd): ", elem_a
       if( myrank == 0 ) write(*,*) " Deallocate &
                   &G2(b-a,b-b,b-c,b-d),&
                   &G5(b-aa,b-ab,b-ac,b-ad,b-bb,b-bc,b-bd,b-cc,b-cd,b-dd): ", elem_b
       if( myrank == 0 ) write(*,*) " Deallocate &
                   &G2(c-a,c-b,c-c,c-d),&
                   &G5(c-aa,c-ab,c-ac,c-ad,c-bb,c-bc,c-bd,c-cc,c-cd,c-dd): ", elem_c
       if( myrank == 0 ) write(*,*) " Deallocate &
                   &G2(d-a,d-b,d-c,d-d),&
                   &G5(d-aa,d-ab,d-ac,d-ad,d-bb,d-bc,d-bd,d-cc,d-cd,d-dd): ", elem_d
     else
       write(*,*) "Error(main): deallocate_4"
       stop
     endif
   endif ! nelem == 4
   ! more than 5 elements
   if( nelem > 4 )then
     write(*,*) "Error(main): deallocate"
     stop
   endif ! more than 5 elements


   deallocate( posi )


   if( myrank == 0 ) write(*,*) " "


 enddo ! nfile, number of tag: bulk_au_n32_300k ...


 close(99) ! input_tag.dat

 if( myrank == 0 )then
   call cpu_time( t1 )
   write(*,*) "Total CPU time ",t1-t0,"[s]"
 endif ! myrank == 0

 call mpi_finalize( ierr )


end
