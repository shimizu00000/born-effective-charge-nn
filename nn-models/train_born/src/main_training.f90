program main_training
use allocarray
implicit none
include 'mpif.h'
integer ierr, num_mpi, myrank, dummyi
integer,allocatable,dimension(:) :: displs, recvcount
integer i, j, k, m
integer stat, counter, nfile, infile
integer imd, md_steps
double precision t0, t1
character(len=20) dummyc
character(len=120) tag, filename, dir_data, dir_sf
character(len=120),allocatable,dimension(:) :: tag_total
! num_data_trn   : number of training data
! num_data_tst   : number of test data
! num_data_total : num_data_trn + num_data_tst
! num_data_x     : number of data for each element
! list           : list of trn and tst data
integer num_data_trn, num_data_tst, num_data_total
integer num_data_trn_mpi, num_data_tst_mpi, num_data_total_mpi
integer num_data_a, num_data_b, num_data_c
integer num_data_a_trn, num_data_b_trn, num_data_c_trn
integer num_data_a_trn_mpi, num_data_b_trn_mpi, num_data_c_trn_mpi
integer,allocatable,dimension(:) :: &
  num_atom, num_atom_a, num_atom_b, num_atom_c, num_atom_mpi
integer natom_total, natom_total_mpi
integer natom_a_total, natom_ab_total, natom_ac_total, natom_abc_total
integer natom_b_total, natom_ba_total, natom_bc_total, natom_bac_total
integer natom_c_total, natom_ca_total, natom_cb_total, natom_cab_total
integer natom_a_sum, natom_ab_sum, natom_ac_sum, natom_abc_sum
integer natom_b_sum, natom_ba_sum, natom_bc_sum, natom_bac_sum
integer natom_c_sum, natom_ca_sum, natom_cb_sum, natom_cab_sum
integer,allocatable,dimension(:)   :: list, list_mpi
integer,allocatable,dimension(:) :: info_io, info_io_mpi
double precision seed(2)
double precision,allocatable,dimension(:) :: ransu_x, ransu_y
double precision,allocatable,dimension(:) :: ransu_r, ransu_t
double precision,allocatable,dimension(:) :: nransu
double precision,parameter :: pi = 4.0d0 * datan(1.0d0)
! nelem   : total number of elements
! natom   : total number of atoms
! natom_x : distribute natom to each element
integer nelem, natom
integer natom_a, natom_b, natom_c
integer io_a, io_b, io_c
integer counter_ene, counter_io
integer counter_a, counter_b, counter_c
integer counter_natom
integer counter_natom_a,   counter_natom_ab,  counter_natom_ac,  counter_natom_abc
integer counter_natom_b,   counter_natom_ba,  counter_natom_bc,  counter_natom_bac
integer counter_natom_c,   counter_natom_ca,  counter_natom_cb,  counter_natom_cab
integer counter_deriv_a,   counter_deriv_ab,  counter_deriv_ac,  counter_deriv_abc
integer counter_deriv_b,   counter_deriv_ba,  counter_deriv_bc,  counter_deriv_bac
integer counter_deriv_c,   counter_deriv_ca,  counter_deriv_cb,  counter_deriv_cab
character(len=4) elem_a, elem_b, elem_c

! Variables for L-BFGS-B
!------------------------
! nparam : number of parameters
! factr, pgtol : convergence criteria
integer nparam
!integer,parameter :: memory = 5
!integer,parameter :: iprint = 1
integer,parameter :: dp = kind( 1.0d0 )
!double precision,parameter :: factr = 1.0d+7
!double precision,parameter :: pgtol = 1.0d-5
character(len=60) :: task, csave
logical           :: lsave(4)
integer           :: isave(44), memory, iprint, weight_counter, itest
double precision  :: dsave(29), factr, pgtol, f_store
integer,allocatable,dimension(:) :: nbd, iwa
double precision,allocatable,dimension(:) :: weight, l, u, wa, g, weight_store
! f : error function
double precision  :: f, f_mpi
! g  : diffential of error function with respect to weight parameters
! g_x: distribute g to each element
double precision,allocatable,dimension(:) :: &
  g_a, g_e_a, g_f_a, g_a_mpi, &
  g_b, g_e_b, g_f_b, g_b_mpi, &
  g_c, g_e_c, g_f_c, g_c_mpi

! Variables for NNP
!-------------------
integer io_nnp_a,  io_nnp_b,  io_nnp_c
integer nlayer_a,  nlayer_b,  nlayer_c
integer nweight_a, nweight_b, nweight_c
integer nhidden_a, nhidden_b, nhidden_c
integer,allocatable,dimension(:) :: &
  network_a, network_b, network_c
double precision,allocatable,dimension(:) :: &
  weight_a, weight_b, weight_c
double precision,allocatable,dimension(:) :: &
  weight_inv_a, weight_inv_b, weight_inv_c
double precision,allocatable,dimension(:) :: &
  atomic_ene_a, atomic_ene_b, atomic_ene_c
double precision,allocatable,dimension(:,:) :: &
  hidden_a, hidden_b, hidden_c
double precision energy, energy0
double precision,allocatable,dimension(:) :: &
  energy0_all, energy_all
double precision,allocatable,dimension(:,:) :: force, force0
double precision,allocatable,dimension(:) :: &
  force_all, force0_all
double precision alpha, beta
! mse_energy_trn : mean square error of energy in trn
! mse_energy_tst : mean square error of energy in tst
double precision mse_energy_trn, mse_energy_tst
double precision mse_energy_trn_mpi, mse_energy_tst_mpi
double precision mse_force_trn, mse_force_tst
double precision mse_force_trn_mpi, mse_force_tst_mpi
! io_weight   : weight from scratch or file
! percent_tst : fraction of tst data
integer io_weight, percent_tst
! rc : cutoff radius
! rn : cutoff radius + rn
! neighbor : number of neighboring atoms
! phi : activation function
double precision rc, rn
integer neighbor
character(len=10) phi

! Variables for symmetry function
!---------------------------------
!A
integer &
  num_g2_a_a,  num_g2_a_b,  num_g2_a_c, &
  num_g5_a_aa, num_g5_a_ab, num_g5_a_ac, &
  num_g5_a_bb, num_g5_a_bc, &
  num_g5_a_cc, &
  num_g_a(9)
!B
integer &
  num_g2_b_a,  num_g2_b_b,  num_g2_b_c, &
  num_g5_b_aa, num_g5_b_ab, num_g5_b_ac, &
  num_g5_b_bb, num_g5_b_bc, &
  num_g5_b_cc, &
  num_g_b(9)
!C
integer &
  num_g2_c_a,  num_g2_c_b,  num_g2_c_c, &
  num_g5_c_aa, num_g5_c_ab, num_g5_c_ac, &
  num_g5_c_bb, num_g5_c_bc, &
  num_g5_c_cc, &
  num_g_c(9)
!A
double precision,allocatable,dimension(:,:) :: &
  g2_a_a,  g2_a_b,  g2_a_c, &
  g5_a_aa, g5_a_ab, g5_a_ac, &
  g5_a_bb, g5_a_bc, &
  g5_a_cc
!B
double precision,allocatable,dimension(:,:) :: &
  g2_b_a,  g2_b_b,  g2_b_c, &
  g5_b_aa, g5_b_ab, g5_b_ac, &
  g5_b_bb, g5_b_bc, &
  g5_b_cc
!C
double precision,allocatable,dimension(:,:) :: &
  g2_c_a,  g2_c_b,  g2_c_c, &
  g5_c_aa, g5_c_ab, g5_c_ac, &
  g5_c_bb, g5_c_bc, &
  g5_c_cc
!A
double precision,allocatable,dimension(:,:,:) :: &
  g2_deriv_a_a,  g2_deriv_a_b,  g2_deriv_a_c, &
  g5_deriv_a_aa, g5_deriv_a_ab, g5_deriv_a_ac, &
  g5_deriv_a_bb, g5_deriv_a_bc, &
  g5_deriv_a_cc
!B
double precision,allocatable,dimension(:,:,:) :: &
  g2_deriv_b_a,  g2_deriv_b_b,  g2_deriv_b_c, &
  g5_deriv_b_aa, g5_deriv_b_ab, g5_deriv_b_ac, &
  g5_deriv_b_bb, g5_deriv_b_bc, &
  g5_deriv_b_cc
!C
double precision,allocatable,dimension(:,:,:) :: &
  g2_deriv_c_a,  g2_deriv_c_b,  g2_deriv_c_c, &
  g5_deriv_c_aa, g5_deriv_c_ab, g5_deriv_c_ac, &
  g5_deriv_c_bb, g5_deriv_c_bc, &
  g5_deriv_c_cc
!A
double precision,allocatable,dimension(:,:) :: &
  g2_a_a_maxmin,  g2_a_b_maxmin,  g2_a_c_maxmin, &
  g5_a_aa_maxmin, g5_a_ab_maxmin, g5_a_ac_maxmin, &
  g5_a_bb_maxmin, g5_a_bc_maxmin, &
  g5_a_cc_maxmin
!B
double precision,allocatable,dimension(:,:) :: &
  g2_b_a_maxmin,  g2_b_b_maxmin,  g2_b_c_maxmin, &
  g5_b_aa_maxmin, g5_b_ab_maxmin, g5_b_ac_maxmin, &
  g5_b_bb_maxmin, g5_b_bc_maxmin, &
  g5_b_cc_maxmin
!C
double precision,allocatable,dimension(:,:) :: &
  g2_c_a_maxmin,  g2_c_b_maxmin,  g2_c_c_maxmin, &
  g5_c_aa_maxmin, g5_c_ab_maxmin, g5_c_ac_maxmin, &
  g5_c_bb_maxmin, g5_c_bc_maxmin, &
  g5_c_cc_maxmin


 !MPI INIT
 call mpi_init( ierr )
 call mpi_comm_size( mpi_comm_world, num_mpi, ierr )
 call mpi_comm_rank( mpi_comm_world, myrank, ierr )
 write(*,*) " NUM_MPI: ", myrank, num_mpi
 call mpi_barrier( mpi_comm_world, ierr )

 if( myrank == 0 ) call cpu_time( t0 )

 !*************************************************************
 ! READ input_nnp.dat FILE
 !*************************************************************
 if( myrank == 0 )then

   ! 3 elements
   call read_input_nnp3( &
        elem_a, &
        num_g2_a_a,  num_g2_a_b,  num_g2_a_c, &
        num_g5_a_aa, num_g5_a_ab, num_g5_a_ac, &
        num_g5_a_bb, num_g5_a_bc, &
        num_g5_a_cc, &
        elem_b, &
        num_g2_b_a,  num_g2_b_b,  num_g2_b_c, &
        num_g5_b_aa, num_g5_b_ab, num_g5_b_ac, &
        num_g5_b_bb, num_g5_b_bc, &
        num_g5_b_cc, &
        elem_c, &
        num_g2_c_a,  num_g2_c_b,  num_g2_c_c, &
        num_g5_c_aa, num_g5_c_ab, num_g5_c_ac, &
        num_g5_c_bb, num_g5_c_bc, &
        num_g5_c_cc, &
        rc, rn, phi, neighbor, &
        dir_data, dir_sf, itest, &
        io_weight, percent_tst, alpha, beta, &
        factr, pgtol, memory, iprint )

   num_g_a(1) = num_g2_a_a  ; num_g_a(2) = num_g2_a_b  ; num_g_a(3) = num_g2_a_c
   num_g_a(4) = num_g5_a_aa ; num_g_a(5) = num_g5_a_ab ; num_g_a(6) = num_g5_a_ac
   num_g_a(7) = num_g5_a_bb ; num_g_a(8) = num_g5_a_bc
   num_g_a(9) = num_g5_a_cc

   num_g_b(1) = num_g2_b_a  ; num_g_b(2) = num_g2_b_b  ; num_g_b(3) = num_g2_b_c
   num_g_b(4) = num_g5_b_aa ; num_g_b(5) = num_g5_b_ab ; num_g_b(6) = num_g5_b_ac
   num_g_b(7) = num_g5_b_bb ; num_g_b(8) = num_g5_b_bc
   num_g_b(9) = num_g5_b_cc

   num_g_c(1) = num_g2_c_a  ; num_g_c(2) = num_g2_c_b  ; num_g_c(3) = num_g2_c_c
   num_g_c(4) = num_g5_c_aa ; num_g_c(5) = num_g5_c_ab ; num_g_c(6) = num_g5_c_ac
   num_g_c(7) = num_g5_c_bb ; num_g_c(8) = num_g5_c_bc
   num_g_c(9) = num_g5_c_cc

   !------------------------------
   ! Read number of hidden layers
   !------------------------------
   call read_layer_node3( nlayer_a, nlayer_b, nlayer_c )

 endif ! myrank == 0

 call mpi_barrier( mpi_comm_world, ierr )
 call mpi_bcast( elem_a, 4, mpi_character, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_a_a, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_a_b, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_a_c, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_a_aa, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_a_ab, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_a_ac, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_a_bb, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_a_bc, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_a_cc, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( elem_b, 4, mpi_character, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_b_a, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_b_b, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_b_c, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_b_aa, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_b_ab, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_b_ac, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_b_bb, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_b_bc, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_b_cc, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( elem_c, 4, mpi_character, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_c_a, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_c_b, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g2_c_c, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_c_aa, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_c_ab, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_c_ac, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_c_bb, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_c_bc, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g5_c_cc, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( rc, 1, mpi_double_precision, 0, mpi_comm_world, ierr )
 call mpi_bcast( dir_data, 120, mpi_character, 0, mpi_comm_world, ierr )
 call mpi_bcast( dir_sf, 120, mpi_character, 0, mpi_comm_world, ierr )
 call mpi_bcast( percent_tst, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( alpha, 1, mpi_double_precision, 0, mpi_comm_world, ierr )
 call mpi_bcast( beta, 1, mpi_double_precision, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g_a, 9, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g_b, 9, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_g_c, 9, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( nlayer_a, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( nlayer_b, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( nlayer_c, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )

 allocate( network_a( nlayer_a ) )
 allocate( network_b( nlayer_b ) )
 allocate( network_c( nlayer_c ) )

 !--------------------------------
 ! Read number of nodes per layer
 !--------------------------------
 if( myrank == 0 )then
   call read_layer_node_frex3( network_a, network_b, network_c, &
                               nlayer_a,  nlayer_b,  nlayer_c, &
                               io_nnp_a,  io_nnp_b,  io_nnp_c )
 endif ! myrank == 0

 call mpi_bcast( network_a, nlayer_a, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( network_b, nlayer_b, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( network_c, nlayer_c, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( io_nnp_a, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( io_nnp_b, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( io_nnp_c, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )

 !if( io_nnp_a == 0 ) deallocate( network_a )
 !if( io_nnp_b == 0 ) deallocate( network_b )
 !if( io_nnp_c == 0 ) deaalocate( network_c )

! Network structure
!-------------------------------------------------------------------------
! network(1)    network(2)      network(3)      network(4)      network(5)
!-------------------------------------------------------------------------
! Bias          Bias            Bias            Bias
! G1            H1-1            H2-1            H3-1            Ea
! G2            H1-2            H2-2            H3-2
! ...           ...             ...             ...
!
! #node = network(layer) + 1(bias)

!------------------------------------------------------------
! Total number of weight parameters
! Total number of nodes (excluding input nodes and output(atomic_ene) node )
!------------------------------------------------------------
!A
 if( io_nnp_a == 1 )then
   nweight_a = 0
   do i = 1, nlayer_a - 1
     nweight_a = nweight_a + ( network_a(i) + 1 ) * network_a(i+1)
   enddo
   allocate( weight_a( nweight_a ), weight_inv_a( nweight_a ) )
   allocate(      g_a( nweight_a ), g_e_a( nweight_a ), g_f_a( nweight_a ), &
              g_a_mpi( nweight_a ) )
   nhidden_a = 0
   do i = 2, nlayer_a - 1
     nhidden_a = nhidden_a + ( network_a(i) + 1 )
   enddo
 endif
!B
 if( io_nnp_b == 1 )then
   nweight_b = 0
   do i = 1, nlayer_b - 1
     nweight_b = nweight_b + ( network_b(i) + 1 ) * network_b(i+1)
   enddo
   allocate( weight_b( nweight_b ), weight_inv_b( nweight_b ) )
   allocate(      g_b( nweight_b ), g_e_b( nweight_b ), g_f_b( nweight_b ), &
              g_b_mpi( nweight_b ) )
   nhidden_b = 0
   do i = 2, nlayer_b - 1
     nhidden_b = nhidden_b + ( network_b(i) + 1 )
   enddo
 endif
!C
 if( io_nnp_c == 1 )then
   nweight_c = 0
   do i = 1, nlayer_c - 1
     nweight_c = nweight_c + ( network_c(i) + 1 ) * network_c(i+1)
   enddo
   allocate( weight_c( nweight_c ), weight_inv_c( nweight_c ) )
   allocate(      g_c( nweight_c ), g_e_c( nweight_c ), g_f_c( nweight_c ), &
              g_c_mpi( nweight_c ) )
   nhidden_c = 0
   do i = 2, nlayer_c - 1
     nhidden_c = nhidden_c + ( network_c(i) + 1 )
   enddo
 endif
 call mpi_barrier( mpi_comm_world, ierr )

 !********************
 ! PARAMETER FOR BFGS
 !********************
 !----------------------
 ! Number of parameters
 !----------------------
 if( myrank == 0 )then
   nparam = 0
   if( io_nnp_a == 1 ) nparam = nparam + nweight_a
   if( io_nnp_b == 1 ) nparam = nparam + nweight_b
   if( io_nnp_c == 1 ) nparam = nparam + nweight_c
   write(*,*) " Number of parameters : ", nparam
   write(*,*) " ===================================="
   if( io_nnp_a == 1 ) write(*,*) " nweight_a: ", nweight_a
   if( io_nnp_b == 1 ) write(*,*) " nweight_b: ", nweight_b
   if( io_nnp_c == 1 ) write(*,*) " nweight_c: ", nweight_c
   write(*,*) " ===================================="
   write(*,*) " "
 endif ! myrank == 0
 call mpi_bcast( nparam, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )

 if( myrank == 0 )then
   allocate( nbd( nparam ), l( nparam ), u( nparam ) )
   allocate( iwa( 3*nparam ) )
   allocate( wa( 2*memory*nparam + 5*nparam + 11*memory*memory + 8*memory ) )
 endif ! myrank == 0
 allocate( g( nparam ), weight( nparam ) )
 call mpi_barrier( mpi_comm_world, ierr )

 !------------------------
 ! Upper and lower bounds
 !------------------------
 if( myrank == 0 )then
   do i = 1, nparam
     nbd(i) =  0 ! 0 nobounded, 2 both bounded
     l(i)   = -1.0d2
     u(i)   =  1.0d2
   enddo
 endif ! myrank == 0
!************************************************************
! DEFINE INITIAL WEIGHT PARAMETERS
!************************************************************
if( myrank == 0 )then
  ! From scrach
  if( io_weight == 0 )then
    allocate( weight_store( nparam ) )
    weight_counter = 0; f_store = 1.0d10; f_mpi = 1.0d10
    allocate ( ransu_x( (nparam + 2) / 2) )
    allocate ( ransu_y( (nparam + 2) / 2) )
    allocate ( ransu_r( (nparam + 2) / 2) )
    allocate ( ransu_t( (nparam + 2) / 2) )
    allocate ( nransu( nparam) )
    call pre_random
    call random_number( ransu_x )
    call pre_random
    call random_number( ransu_y )

    ransu_r = dsqrt(-2.0 * dlog(ransu_x))
    ransu_t = 2.0 * pi * ransu_y

    nransu(1:((nparam + 2) / 2)) = ransu_r * dcos(ransu_t)

    do i = ((nparam + 2) / 2) + 1, nparam
      nransu(i) = ransu_r(i - (nparam + 2) / 2) * dsin(ransu_t(i - (nparam + 2) / 2))
    enddo

    if (nlayer_a == 4 .and. nlayer_b == 4 .and. nlayer_c == 4)then
      do i = 1, (network_a(1) + 1) * network_a(2)
        weight(i) = nransu(i) / dsqrt(dble(network_a(1) + 1))
      enddo
      do i = 1, (network_a(2) + 1) * network_a(3)
        weight((network_a(1) + 1) * network_a(2) + 1 + i) = nransu((network_a(1) + 1) * network_a(2) + i) / dsqrt(dble(network_a(2) + 1))
      enddo
      do i = 1, (network_a(3) + 1) * network_a(4)
        weight((network_a(1) + 1) * network_a(2) + (network_a(2) + 1) * network_a(3) + i) = nransu((network_a(1) + 1) * network_a(2) + (network_a(2) + 1) * network_a(3) + i) / dsqrt(dble(network_a(3) + 1))
      enddo
      do i = 1, (network_b(1) + 1) * network_b(2)
        weight(nweight_a + i) = nransu(nweight_a + i) / dsqrt(dble(network_a(1)))
      enddo
      do i = 1, (network_b(2) + 1) * network_b(3)
        weight(nweight_a + (network_b(1) + 1) * network_b(2) + i) = nransu(nweight_a + (network_b(1) + 1)* network_b(2) + i) / dsqrt(dble(network_b(2) + 1))
      enddo
      do i = 1, (network_b(3) + 1) * network_b(4)
        weight(nweight_a + (network_b(1) + 1) * network_b(2) + (network_b(2) + 1) * network_b(3) + i) = nransu(nweight_a + (network_b(1) + 1) * network_b(2) + (network_b(2) + 1) * network_b(3) + i) / dsqrt(dble(network_b(3) + 1))
      enddo
      do i = 1, (network_c(1) + 1) * network_c(2)
        weight(nweight_a + nweight_b + i) = nransu(nweight_a + nweight_b + i) / dsqrt(dble(network_c(1) + 1))
      enddo
      do i = 1, (network_c(2) + 1) * network_c(3)
        weight(nweight_a + nweight_b + (network_c(1) + 1) * network_c(2) + i) = nransu(nweight_a + nweight_b + (network_c(1) + 1) * network_c(2) + i) / dsqrt(dble(network_c(2) + 1))
      enddo
      do i = 1, (network_a(3) + 1) * network_a(4)
        weight(nweight_a + nweight_b + (network_c(1) + 1) * network_c(2) + (network_c(2) + 1)* network_c(3) + i) = nransu(nweight_a + nweight_b + (network_c(1) + 1) * network_c(2) + (network_c(2) + 1)* network_c(3) + i) / dsqrt(dble(network_c(3) + 1))
      enddo

    endif

    deallocate( ransu_x, ransu_y, ransu_r, ransu_t, nransu )
  ! From file
  elseif( io_weight == 1 )then
    open(unit=7,file="./data_weight/weight.dat",action="read",form="unformatted")
    read(7) weight
    close(7)
  ! Bug
  else
    write(*,*) "Error(main): weight parameter setting."
    stop
  endif
endif ! myrank == 0

call mpi_bcast( weight, nparam, mpi_double_precision, 0, mpi_comm_world, ierr )
call mpi_barrier( mpi_comm_world, ierr )


 !*************************************************************
 ! COUNTING NUMBER OF DATA from input_tag_training.dat FILE
 !*************************************************************
 if( myrank == 0 )then
   open(unit=99,file="input_tag_training.dat",action="read")
   counter = 0
   rewind(99)
   do
     read( 99, '(A5)', iostat=stat ) dummyc
     if( stat /= 0 ) exit
     if( trim( dummyc ) /= '' ) counter = counter + 1
   enddo
 endif ! myrank == 0
 call mpi_bcast( counter, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )


 !--------------------------------------
 ! DO LOOP : nfile (input_tag_training)
 !--------------------------------------
 allocate( tag_total( counter ) )
 if( myrank == 0 )then
   rewind(99)
   do nfile = 1, counter
     read(99,*) tag
     tag_total( nfile ) = trim(adjustl(tag))
   enddo
 endif ! myrank == 0
 close(99) ! input_tag_training.dat

 call mpi_bcast( tag_total, counter*120, mpi_character, 0, mpi_comm_world, ierr )
 !bcast character: #data*len
 call mpi_barrier( mpi_comm_world, ierr )

 num_data_total = 0
 if( io_nnp_a == 1 ) num_data_a = 0
 if( io_nnp_b == 1 ) num_data_b = 0
 if( io_nnp_c == 1 ) num_data_c = 0

 natom_total = 0
 !A,B,C
 if( io_nnp_a == 1 .and. io_nnp_b == 1 .and. io_nnp_c == 1 )then
   natom_a_total   = 0 ; natom_a_sum   = 0
   natom_ab_total  = 0 ; natom_ab_sum  = 0
   natom_ac_total  = 0 ; natom_ac_sum  = 0
   natom_abc_total = 0 ; natom_abc_sum = 0
   natom_b_total   = 0 ; natom_b_sum   = 0
   natom_ba_total  = 0 ; natom_ba_sum  = 0
   natom_bc_total  = 0 ; natom_bc_sum  = 0
   natom_bac_total = 0 ; natom_bac_sum = 0
   natom_c_total   = 0 ; natom_c_sum   = 0
   natom_ca_total  = 0 ; natom_ca_sum  = 0
   natom_cb_total  = 0 ; natom_cb_sum  = 0
   natom_cab_total = 0 ; natom_cab_sum = 0
   !total: natom_x (g2, g5)
   !sum: natom_x * natom (g2_deriv, g5_deriv)
 endif

 do nfile = 1, ceiling( dble(counter)/dble(num_mpi) )
   !infile = nfile + myrank*ceiling( dble(counter)/dble(num_mpi) )
   infile = ( myrank + 1 ) + ( ( nfile - 1 ) * num_mpi )
   if( infile <= counter )then

     tag = tag_total( infile )
     ! info (binary), READ
     !filename = '../step2_data_binary/data/info_'//trim(adjustl(tag))//'.dat'
     filename = trim(adjustl(dir_data))//'/info_'//trim(adjustl(tag))//'.dat'
     open(unit=1,file=filename,action="read",form="unformatted")
     !---------------------------------
     ! READ info_tag.dat FILE (BINARY)
     !---------------------------------
     ! 3 elements
     call read_io3( natom,   nelem, &
                    elem_a,  elem_b,  elem_c, &
                    natom_a, natom_b, natom_c, &
                    io_a,    io_b,    io_c )
     ! MD steps
     read(1) md_steps

     num_data_total = num_data_total + md_steps
     natom_total = natom_total + natom * md_steps

     !A
     if(     io_a == 1 .and. io_b == 0 .and. io_c == 0 )then
       num_data_a    = num_data_a    + md_steps
       natom_a_total = natom_a_total + natom_a * md_steps
       natom_a_sum   = natom_a_sum   + natom_a * natom * md_steps
     !B
     elseif( io_a == 0 .and. io_b == 1 .and. io_c == 0 )then
       num_data_b    = num_data_b    + md_steps
       natom_b_total = natom_b_total + natom_b * md_steps
       natom_b_sum   = natom_b_sum   + natom_b * natom * md_steps
     !C
     elseif( io_a == 0 .and. io_b == 0 .and. io_c == 1 )then
       num_data_c    = num_data_c    + md_steps
       natom_c_total = natom_c_total + natom_c * md_steps
       natom_c_sum   = natom_c_sum   + natom_c * natom * md_steps
     !A,B
     elseif( io_a == 1 .and. io_b == 1 .and. io_c == 0 )then
       num_data_a     = num_data_a     + md_steps
       natom_a_total  = natom_a_total  + natom_a * md_steps
       natom_a_sum    = natom_a_sum    + natom_a * natom * md_steps
       natom_ab_total = natom_ab_total + natom_a * md_steps
       natom_ab_sum   = natom_ab_sum   + natom_a * natom * md_steps
       num_data_b     = num_data_b     + md_steps
       natom_b_total  = natom_b_total  + natom_b * md_steps
       natom_b_sum    = natom_b_sum    + natom_b * natom * md_steps
       natom_ba_total = natom_ba_total + natom_b * md_steps
       natom_ba_sum   = natom_ba_sum   + natom_b * natom * md_steps
     !A,C
     elseif( io_a == 1 .and. io_b == 0 .and. io_c == 1 )then
       num_data_a     = num_data_a     + md_steps
       natom_a_total  = natom_a_total  + natom_a * md_steps
       natom_a_sum    = natom_a_sum    + natom_a * natom * md_steps
       natom_ac_total = natom_ac_total + natom_a * md_steps
       natom_ac_sum   = natom_ac_sum   + natom_a * natom * md_steps
       num_data_c     = num_data_c     + md_steps
       natom_c_total  = natom_c_total  + natom_c * md_steps
       natom_c_sum    = natom_c_sum    + natom_c * natom * md_steps
       natom_ca_total = natom_ca_total + natom_c * md_steps
       natom_ca_sum   = natom_ca_sum   + natom_c * natom * md_steps
     !B,C
     elseif( io_a == 0 .and. io_b == 1 .and. io_c == 1 )then
       num_data_b     = num_data_b     + md_steps
       natom_b_total  = natom_b_total  + natom_b * md_steps
       natom_b_sum    = natom_b_sum    + natom_b * natom * md_steps
       natom_bc_total = natom_bc_total + natom_b * md_steps
       natom_bc_sum   = natom_bc_sum   + natom_b * natom * md_steps
       num_data_c     = num_data_c     + md_steps
       natom_c_total  = natom_c_total  + natom_c * md_steps
       natom_c_sum    = natom_c_sum    + natom_c * natom * md_steps
       natom_cb_total = natom_cb_total + natom_c * md_steps
       natom_cb_sum   = natom_cb_sum   + natom_c * natom * md_steps
     !A,B,C
     elseif( io_a == 1 .and. io_b == 1 .and. io_c == 1 )then
       num_data_a      = num_data_a      + md_steps
       natom_a_total   = natom_a_total   + natom_a * md_steps
       natom_a_sum     = natom_a_sum     + natom_a * natom * md_steps
       natom_ab_total  = natom_ab_total  + natom_a * md_steps
       natom_ab_sum    = natom_ab_sum    + natom_a * natom * md_steps
       natom_ac_total  = natom_ac_total  + natom_a * md_steps
       natom_ac_sum    = natom_ac_sum    + natom_a * natom * md_steps
       natom_abc_total = natom_abc_total + natom_a * md_steps
       natom_abc_sum   = natom_abc_sum   + natom_a * natom * md_steps
       num_data_b      = num_data_b      + md_steps
       natom_b_total   = natom_b_total   + natom_b * md_steps
       natom_b_sum     = natom_b_sum     + natom_b * natom * md_steps
       natom_ba_total  = natom_ba_total  + natom_b * md_steps
       natom_ba_sum    = natom_ba_sum    + natom_b * natom * md_steps
       natom_bc_total  = natom_bc_total  + natom_b * md_steps
       natom_bc_sum    = natom_bc_sum    + natom_b * natom * md_steps
       natom_bac_total = natom_bac_total + natom_b * md_steps
       natom_bac_sum   = natom_bac_sum   + natom_b * natom * md_steps
       num_data_c      = num_data_c      + md_steps
       natom_c_total   = natom_c_total   + natom_c * md_steps
       natom_c_sum     = natom_c_sum     + natom_c * natom * md_steps
       natom_ca_total  = natom_ca_total  + natom_c * md_steps
       natom_ca_sum    = natom_ca_sum    + natom_c * natom * md_steps
       natom_cb_total  = natom_cb_total  + natom_c * md_steps
       natom_cb_sum    = natom_cb_sum    + natom_c * natom * md_steps
       natom_cab_total = natom_cab_total + natom_c * md_steps
       natom_cab_sum   = natom_cab_sum   + natom_c * natom * md_steps
     else
       write(*,*)"Error(main): counting data io."
       stop
     endif

     close(1) ! info(binary)

   endif ! infile <= counter
 enddo ! nfile
 call mpi_barrier( mpi_comm_world, ierr )
 !--------------
 ! END COUNTING
 !--------------

 !***********
 ! READ DATA
 !***********
 !A,B,C
 if( io_nnp_a == 1 .and. io_nnp_b == 1 .and. io_nnp_c == 1 )then
   call alloc_3( num_atom, num_data_total, &
                 num_atom_a, num_atom_b, num_atom_c, &
                 num_data_a, num_data_b, num_data_c, &
                 natom_a_total, natom_ab_total, natom_ac_total, natom_abc_total, &
                 natom_b_total, natom_ba_total, natom_bc_total, natom_bac_total, &
                 natom_c_total, natom_ca_total, natom_cb_total, natom_cab_total, &
                 natom_a_sum,   natom_ab_sum,  natom_ac_sum,  natom_abc_sum, &
                 natom_b_sum,   natom_ba_sum,  natom_bc_sum,  natom_bac_sum, &
                 natom_c_sum,   natom_ca_sum,  natom_cb_sum,  natom_cab_sum, &
                 "g2_a_a",  g2_a_a,  g2_deriv_a_a,  g2_a_a_maxmin,  num_g2_a_a, &
                 "g2_a_b",  g2_a_b,  g2_deriv_a_b,  g2_a_b_maxmin,  num_g2_a_b, &
                 "g2_a_c",  g2_a_c,  g2_deriv_a_c,  g2_a_c_maxmin,  num_g2_a_c, &
                 "g5_a_aa", g5_a_aa, g5_deriv_a_aa, g5_a_aa_maxmin, num_g5_a_aa, &
                 "g5_a_ab", g5_a_ab, g5_deriv_a_ab, g5_a_ab_maxmin, num_g5_a_ab, &
                 "g5_a_ac", g5_a_ac, g5_deriv_a_ac, g5_a_ac_maxmin, num_g5_a_ac, &
                 "g5_a_bb", g5_a_bb, g5_deriv_a_bb, g5_a_bb_maxmin, num_g5_a_bb, &
                 "g5_a_bc", g5_a_bc, g5_deriv_a_bc, g5_a_bc_maxmin, num_g5_a_bc, &
                 "g5_a_cc", g5_a_cc, g5_deriv_a_cc, g5_a_cc_maxmin, num_g5_a_cc, &
                 "g2_b_a",  g2_b_a,  g2_deriv_b_a,  g2_b_a_maxmin,  num_g2_b_a, &
                 "g2_b_b",  g2_b_b,  g2_deriv_b_b,  g2_b_b_maxmin,  num_g2_b_b, &
                 "g2_b_c",  g2_b_c,  g2_deriv_b_c,  g2_b_c_maxmin,  num_g2_b_c, &
                 "g5_b_aa", g5_b_aa, g5_deriv_b_aa, g5_b_aa_maxmin, num_g5_b_aa, &
                 "g5_b_ab", g5_b_ab, g5_deriv_b_ab, g5_b_ab_maxmin, num_g5_b_ab, &
                 "g5_b_ac", g5_b_ac, g5_deriv_b_ac, g5_b_ac_maxmin, num_g5_b_ac, &
                 "g5_b_bb", g5_b_bb, g5_deriv_b_bb, g5_b_bb_maxmin, num_g5_b_bb, &
                 "g5_b_bc", g5_b_bc, g5_deriv_b_bc, g5_b_bc_maxmin, num_g5_b_bc, &
                 "g5_b_cc", g5_b_cc, g5_deriv_b_cc, g5_b_cc_maxmin, num_g5_b_cc, &
                 "g2_c_a",  g2_c_a,  g2_deriv_c_a,  g2_c_a_maxmin,  num_g2_c_a, &
                 "g2_c_b",  g2_c_b,  g2_deriv_c_b,  g2_c_b_maxmin,  num_g2_c_b, &
                 "g2_c_c",  g2_c_c,  g2_deriv_c_c,  g2_c_c_maxmin,  num_g2_c_c, &
                 "g5_c_aa", g5_c_aa, g5_deriv_c_aa, g5_c_aa_maxmin, num_g5_c_aa, &
                 "g5_c_ab", g5_c_ab, g5_deriv_c_ab, g5_c_ab_maxmin, num_g5_c_ab, &
                 "g5_c_ac", g5_c_ac, g5_deriv_c_ac, g5_c_ac_maxmin, num_g5_c_ac, &
                 "g5_c_bb", g5_c_bb, g5_deriv_c_bb, g5_c_bb_maxmin, num_g5_c_bb, &
                 "g5_c_bc", g5_c_bc, g5_deriv_c_bc, g5_c_bc_maxmin, num_g5_c_bc, &
                 "g5_c_cc", g5_c_cc, g5_deriv_c_cc, g5_c_cc_maxmin, num_g5_c_cc )
 endif

 ! 3 elements system
 allocate( info_io( num_data_total*3 ) )
 allocate( energy0_all( num_data_total ) )
 allocate( force0_all( natom_total*3 ) )

 !--------------------------------------
 ! DO LOOP : nfile (input_tag_training)
 !--------------------------------------
 if( myrank == 0 )then
   write(*,*)" Reading data"
   write(*,*)" "
 endif ! myrank == 0

 counter_io = 0
 counter_ene = 0
 if( io_nnp_a == 1 ) counter_a = 0
 if( io_nnp_b == 1 ) counter_b = 0
 if( io_nnp_c == 1 ) counter_c = 0

 counter_natom = 0
 !A,B,C
 if( io_nnp_a == 1 .and. io_nnp_b == 1 .and. io_nnp_c == 1 )then
   counter_natom_a   = 0 ; counter_deriv_a   = 0
   counter_natom_ab  = 0 ; counter_deriv_ab  = 0
   counter_natom_ac  = 0 ; counter_deriv_ac  = 0
   counter_natom_abc = 0 ; counter_deriv_abc = 0
   counter_natom_b   = 0 ; counter_deriv_b   = 0
   counter_natom_ba  = 0 ; counter_deriv_ba  = 0
   counter_natom_bc  = 0 ; counter_deriv_bc  = 0
   counter_natom_bac = 0 ; counter_deriv_bac = 0
   counter_natom_c   = 0 ; counter_deriv_c   = 0
   counter_natom_ca  = 0 ; counter_deriv_ca  = 0
   counter_natom_cb  = 0 ; counter_deriv_cb  = 0
   counter_natom_cab = 0 ; counter_deriv_cab = 0
   ! counter_deriv_x: natom_x * natom (g2_deriv & g5_deriv)
 endif

 do nfile = 1, ceiling( dble(counter)/dble(num_mpi) )
   !infile = nfile + myrank*ceiling( dble(counter)/dble(num_mpi) )
   infile = ( myrank + 1 ) + ( ( nfile - 1 ) * num_mpi )
   if( infile <= counter )then

     tag = tag_total( infile )
     ! info (binary), READ
     !filename = '../step2_data_binary/data/info_'//trim(adjustl(tag))//'.dat'
     filename = trim(adjustl(dir_data))//'/info_'//trim(adjustl(tag))//'.dat'
     open(unit=1,file=filename,action="read",form="unformatted")
     ! energies (binary), READ
     !filename = '../step2_data_binary/data/energies_'//trim(adjustl(tag))//'.dat'
     filename = trim(adjustl(dir_data))//'/energies_'//trim(adjustl(tag))//'.dat'
     open(unit=3,file=filename,action="read",form="unformatted")
     ! forces (binary), READ
     !filename = '../step2_data_binary/data/forces_'//trim(adjustl(tag))//'.dat'
     filename = trim(adjustl(dir_data))//'/forces_'//trim(adjustl(tag))//'.dat'
     open(unit=4,file=filename,action="read",form="unformatted")
     ! sf (binary), READ
     !write( dummyc, * ) int(rc)
     !filename = '../step3_sf_binary/data_sf/rc_'//trim(adjustl(dummyc))//&
     !           '/sf_'//trim(adjustl(tag))//'.dat'
     filename = trim(adjustl(dir_sf))//'/sf_'//trim(adjustl(tag))//'.dat'
     open(unit=11,file=filename,action="read",form="unformatted")
     ! sf_deriv (binary), READ
     !filename = '../step3_sf_binary/data_sf/rc_'//trim(adjustl(dummyc))//&
     !           '/sf_deriv_'//trim(adjustl(tag))//'.dat'
     filename = trim(adjustl(dir_sf))//'/sf_deriv_'//trim(adjustl(tag))//'.dat'
     open(unit=12,file=filename,action="read",form="unformatted")

     !---------------------------------
     ! READ info_tag.dat FILE (BINARY)
     !---------------------------------
     ! 3 elements system
     call read_io3( natom,   nelem, &
                    elem_a,  elem_b,  elem_c, &
                    natom_a, natom_b, natom_c, &
                    io_a,    io_b,    io_c )
     ! Force
     allocate( force0( natom, 3 ) )

     ! MD steps
     read(1) md_steps
     do imd = 1, md_steps

       read(3) energy0
       counter_ene = counter_ene + 1
       energy0_all( counter_ene ) = energy0
       num_atom( counter_ene ) = natom

       read(4) force0
       do i = 1, natom
         force0_all( counter_natom + 3*(i-1) + 1 ) = force0(i,1)
         force0_all( counter_natom + 3*(i-1) + 2 ) = force0(i,2)
         force0_all( counter_natom + 3*(i-1) + 3 ) = force0(i,3)
       enddo
       counter_natom = counter_natom + 3*natom

       info_io( counter_io + 1 ) = io_a
       info_io( counter_io + 2 ) = io_b
       info_io( counter_io + 3 ) = io_c
       counter_io = counter_io + 3

       ! 1 element
       !A
       if(     io_a == 1 .and. io_b == 0 .and. io_c == 0 )then
         call read_data_1( num_atom_a, num_data_a, natom_a_total, natom_a_sum, &
                           counter_a, counter_natom_a, counter_deriv_a, &
                           imd, natom_a, natom, &
                           "g2_a_a",  g2_a_a,  g2_deriv_a_a,  g2_a_a_maxmin,  num_g2_a_a, &
                           "g5_a_aa", g5_a_aa, g5_deriv_a_aa, g5_a_aa_maxmin, num_g5_a_aa )
       !B
       elseif( io_a == 0 .and. io_b == 1 .and. io_c == 0 )then
         call read_data_1( num_atom_b, num_data_b, natom_b_total, natom_b_sum, &
                           counter_b, counter_natom_b, counter_deriv_b, &
                           imd, natom_b, natom, &
                           "g2_b_b",  g2_b_b,  g2_deriv_b_b,  g2_b_b_maxmin,  num_g2_b_b, &
                           "g5_b_bb", g5_b_bb, g5_deriv_b_bb, g5_b_bb_maxmin, num_g5_b_bb )
       !C
       elseif( io_a == 0 .and. io_b == 0 .and. io_c == 1 )then
         call read_data_1( num_atom_c, num_data_c, natom_c_total, natom_c_sum, &
                           counter_c, counter_natom_c, counter_deriv_c, &
                           imd, natom_c, natom, &
                           "g2_c_c",  g2_c_c,  g2_deriv_c_c,  g2_c_c_maxmin,  num_g2_c_c, &
                           "g5_c_cc", g5_c_cc, g5_deriv_c_cc, g5_c_cc_maxmin, num_g5_c_cc )
       !A,B
       elseif( io_a == 1 .and. io_b == 1 .and. io_c == 0 )then
         call read_data_2( num_atom_a, num_atom_b, &
                           num_data_a, num_data_b, &
                           natom_a_total, natom_ab_total, &
                           natom_b_total, natom_ba_total, &
                           natom_a_sum, natom_ab_sum, &
                           natom_b_sum, natom_ba_sum, &
                           counter_a, counter_b, &
                           counter_natom_a, counter_natom_ab, &
                           counter_natom_b, counter_natom_ba, &
                           counter_deriv_a, counter_deriv_ab, &
                           counter_deriv_b, counter_deriv_ba, &
                           imd, natom_a, natom_b, natom, &
                           "g2_a_a",  g2_a_a,  g2_deriv_a_a,  g2_a_a_maxmin,  num_g2_a_a, &
                           "g2_a_b",  g2_a_b,  g2_deriv_a_b,  g2_a_b_maxmin,  num_g2_a_b, &
                           "g5_a_aa", g5_a_aa, g5_deriv_a_aa, g5_a_aa_maxmin, num_g5_a_aa, &
                           "g5_a_ab", g5_a_ab, g5_deriv_a_ab, g5_a_ab_maxmin, num_g5_a_ab, &
                           "g5_a_bb", g5_a_bb, g5_deriv_a_bb, g5_a_bb_maxmin, num_g5_a_bb, &
                           "g2_b_a",  g2_b_a,  g2_deriv_b_a,  g2_b_a_maxmin,  num_g2_b_a, &
                           "g2_b_b",  g2_b_b,  g2_deriv_b_b,  g2_b_b_maxmin,  num_g2_b_b, &
                           "g5_b_aa", g5_b_aa, g5_deriv_b_aa, g5_b_aa_maxmin, num_g5_b_aa, &
                           "g5_b_ab", g5_b_ab, g5_deriv_b_ab, g5_b_ab_maxmin, num_g5_b_ab, &
                           "g5_b_bb", g5_b_bb, g5_deriv_b_bb, g5_b_bb_maxmin, num_g5_b_bb )
       !A,C
       elseif( io_a == 1 .and. io_b == 0 .and. io_c == 1 )then
         call read_data_2( num_atom_a, num_atom_c, &
                           num_data_a, num_data_c, &
                           natom_a_total, natom_ac_total, &
                           natom_c_total, natom_ca_total, &
                           natom_a_sum, natom_ac_sum, &
                           natom_c_sum, natom_ca_sum, &
                           counter_a, counter_c, &
                           counter_natom_a, counter_natom_ac, &
                           counter_natom_c, counter_natom_ca, &
                           counter_deriv_a, counter_deriv_ac, &
                           counter_deriv_c, counter_deriv_ca, &
                           imd, natom_a, natom_c, natom, &
                           "g2_a_a",  g2_a_a,  g2_deriv_a_a,  g2_a_a_maxmin,  num_g2_a_a, &
                           "g2_a_c",  g2_a_c,  g2_deriv_a_c,  g2_a_c_maxmin,  num_g2_a_c, &
                           "g5_a_aa", g5_a_aa, g5_deriv_a_aa, g5_a_aa_maxmin, num_g5_a_aa, &
                           "g5_a_ac", g5_a_ac, g5_deriv_a_ac, g5_a_ac_maxmin, num_g5_a_ac, &
                           "g5_a_cc", g5_a_cc, g5_deriv_a_cc, g5_a_cc_maxmin, num_g5_a_cc, &
                           "g2_c_a",  g2_c_a,  g2_deriv_c_a,  g2_c_a_maxmin,  num_g2_c_a, &
                           "g2_c_c",  g2_c_c,  g2_deriv_c_c,  g2_c_c_maxmin,  num_g2_c_c, &
                           "g5_c_aa", g5_c_aa, g5_deriv_c_aa, g5_c_aa_maxmin, num_g5_c_aa, &
                           "g5_c_ac", g5_c_ac, g5_deriv_c_ac, g5_c_ac_maxmin, num_g5_c_ac, &
                           "g5_c_cc", g5_c_cc, g5_deriv_c_cc, g5_c_cc_maxmin, num_g5_c_cc )
       !B,C
       elseif( io_a == 0 .and. io_b == 1 .and. io_c == 1 )then
         call read_data_2( num_atom_b, num_atom_c, &
                           num_data_b, num_data_c, &
                           natom_b_total, natom_bc_total, &
                           natom_c_total, natom_cb_total, &
                           natom_b_sum, natom_bc_sum, &
                           natom_c_sum, natom_cb_sum, &
                           counter_b, counter_c, &
                           counter_natom_b, counter_natom_bc, &
                           counter_natom_c, counter_natom_cb, &
                           counter_deriv_b, counter_deriv_bc, &
                           counter_deriv_c, counter_deriv_cb, &
                           imd, natom_b, natom_c, natom, &
                           "g2_b_b",  g2_b_b,  g2_deriv_b_b,  g2_b_b_maxmin,  num_g2_b_b, &
                           "g2_b_c",  g2_b_c,  g2_deriv_b_c,  g2_b_c_maxmin,  num_g2_b_c, &
                           "g5_b_bb", g5_b_bb, g5_deriv_b_bb, g5_b_bb_maxmin, num_g5_b_bb, &
                           "g5_b_bc", g5_b_bc, g5_deriv_b_bc, g5_b_bc_maxmin, num_g5_b_bc, &
                           "g5_b_cc", g5_b_cc, g5_deriv_b_cc, g5_b_cc_maxmin, num_g5_b_cc, &
                           "g2_c_b",  g2_c_b,  g2_deriv_c_b,  g2_c_b_maxmin,  num_g2_c_b, &
                           "g2_c_c",  g2_c_c,  g2_deriv_c_c,  g2_c_c_maxmin,  num_g2_c_c, &
                           "g5_c_bb", g5_c_bb, g5_deriv_c_bb, g5_c_bb_maxmin, num_g5_c_bb, &
                           "g5_c_bc", g5_c_bc, g5_deriv_c_bc, g5_c_bc_maxmin, num_g5_c_bc, &
                           "g5_c_cc", g5_c_cc, g5_deriv_c_cc, g5_c_cc_maxmin, num_g5_c_cc )
       !A,B,C
       elseif( io_a == 1 .and. io_b == 1 .and. io_c == 1 )then
         call read_data_3( num_atom_a, num_atom_b, num_atom_c, &
                           num_data_a, num_data_b, num_data_c, &
                           natom_a_total, natom_ab_total, natom_ac_total, natom_abc_total, &
                           natom_b_total, natom_ba_total, natom_bc_total, natom_bac_total, &
                           natom_c_total, natom_ca_total, natom_cb_total, natom_cab_total, &
                           natom_a_sum, natom_ab_sum, natom_ac_sum, natom_abc_sum, &
                           natom_b_sum, natom_ba_sum, natom_bc_sum, natom_bac_sum, &
                           natom_c_sum, natom_ca_sum, natom_cb_sum, natom_cab_sum, &
                           counter_a, counter_b, counter_c, &
                           counter_natom_a, counter_natom_ab, counter_natom_ac, counter_natom_abc, &
                           counter_natom_b, counter_natom_ba, counter_natom_bc, counter_natom_bac, &
                           counter_natom_c, counter_natom_ca, counter_natom_cb, counter_natom_cab, &
                           counter_deriv_a, counter_deriv_ab, counter_deriv_ac, counter_deriv_abc, &
                           counter_deriv_b, counter_deriv_ba, counter_deriv_bc, counter_deriv_bac, &
                           counter_deriv_c, counter_deriv_ca, counter_deriv_cb, counter_deriv_cab, &
                           imd, natom_a, natom_b, natom_c, natom, &
                           "g2_a_a",  g2_a_a,  g2_deriv_a_a,  g2_a_a_maxmin,  num_g2_a_a, &
                           "g2_a_b",  g2_a_b,  g2_deriv_a_b,  g2_a_b_maxmin,  num_g2_a_b, &
                           "g2_a_c",  g2_a_c,  g2_deriv_a_c,  g2_a_c_maxmin,  num_g2_a_c, &
                           "g5_a_aa", g5_a_aa, g5_deriv_a_aa, g5_a_aa_maxmin, num_g5_a_aa, &
                           "g5_a_ab", g5_a_ab, g5_deriv_a_ab, g5_a_ab_maxmin, num_g5_a_ab, &
                           "g5_a_ac", g5_a_ac, g5_deriv_a_ac, g5_a_ac_maxmin, num_g5_a_ac, &
                           "g5_a_bb", g5_a_bb, g5_deriv_a_bb, g5_a_bb_maxmin, num_g5_a_bb, &
                           "g5_a_bc", g5_a_bc, g5_deriv_a_bc, g5_a_bc_maxmin, num_g5_a_bc, &
                           "g5_a_cc", g5_a_cc, g5_deriv_a_cc, g5_a_cc_maxmin, num_g5_a_cc, &
                           "g2_b_a",  g2_b_a,  g2_deriv_b_a,  g2_b_a_maxmin,  num_g2_b_a, &
                           "g2_b_b",  g2_b_b,  g2_deriv_b_b,  g2_b_b_maxmin,  num_g2_b_b, &
                           "g2_b_c",  g2_b_c,  g2_deriv_b_c,  g2_b_c_maxmin,  num_g2_b_c, &
                           "g5_b_aa", g5_b_aa, g5_deriv_b_aa, g5_b_aa_maxmin, num_g5_b_aa, &
                           "g5_b_ab", g5_b_ab, g5_deriv_b_ab, g5_b_ab_maxmin, num_g5_b_ab, &
                           "g5_b_ac", g5_b_ac, g5_deriv_b_ac, g5_b_ac_maxmin, num_g5_b_ac, &
                           "g5_b_bb", g5_b_bb, g5_deriv_b_bb, g5_b_bb_maxmin, num_g5_b_bb, &
                           "g5_b_bc", g5_b_bc, g5_deriv_b_bc, g5_b_bc_maxmin, num_g5_b_bc, &
                           "g5_b_cc", g5_b_cc, g5_deriv_b_cc, g5_b_cc_maxmin, num_g5_b_cc, &
                           "g2_c_a",  g2_c_a,  g2_deriv_c_a,  g2_c_a_maxmin,  num_g2_c_a, &
                           "g2_c_b",  g2_c_b,  g2_deriv_c_b,  g2_c_b_maxmin,  num_g2_c_b, &
                           "g2_c_c",  g2_c_c,  g2_deriv_c_c,  g2_c_c_maxmin,  num_g2_c_c, &
                           "g5_c_aa", g5_c_aa, g5_deriv_c_aa, g5_c_aa_maxmin, num_g5_c_aa, &
                           "g5_c_ab", g5_c_ab, g5_deriv_c_ab, g5_c_ab_maxmin, num_g5_c_ab, &
                           "g5_c_ac", g5_c_ac, g5_deriv_c_ac, g5_c_ac_maxmin, num_g5_c_ac, &
                           "g5_c_bb", g5_c_bb, g5_deriv_c_bb, g5_c_bb_maxmin, num_g5_c_bb, &
                           "g5_c_bc", g5_c_bc, g5_deriv_c_bc, g5_c_bc_maxmin, num_g5_c_bc, &
                           "g5_c_cc", g5_c_cc, g5_deriv_c_cc, g5_c_cc_maxmin, num_g5_c_cc )
       else
         write(*,*)"Error(main): read data."
         stop
       endif

     enddo ! imd

     deallocate( force0 )

     if( io_a == 1 ) counter_a = counter_a + md_steps
     if( io_b == 1 ) counter_b = counter_b + md_steps
     if( io_c == 1 ) counter_c = counter_c + md_steps

     close(1)  ! info
     close(3)  ! energies
     close(4)  ! forces
     close(11) ! sf
     close(12) ! sf_deriv

   endif ! infile <= counter
 enddo ! nfile
 call mpi_barrier( mpi_comm_world, ierr )

 if( myrank == 0 )then
   write(*,*)" Finish reading data"
   write(*,*)" "
 endif ! myrank == 0
 !---------------
 ! END READ DATA
 !---------------

 !---------------------
 ! List of trn and tst
 !---------------------
 allocate( list( num_data_total ) )
 list = 0
 num_data_trn = 0
 if( io_nnp_a == 1 ) num_data_a_trn = 0
 if( io_nnp_b == 1 ) num_data_b_trn = 0
 if( io_nnp_c == 1 ) num_data_c_trn = 0
 call pre_random

 counter_io = 0
 do i = 1, num_data_total
   call random_number( seed )
   seed(2) = seed(2) * 100.0d0
   ! list = 1: training, list = 0: test
   if( seed(2) > percent_tst )then
     list(i) = 1
     num_data_trn = num_data_trn + 1
     if( info_io( counter_io + 1 ) == 1 ) num_data_a_trn = num_data_a_trn + 1
     if( info_io( counter_io + 2 ) == 1 ) num_data_b_trn = num_data_b_trn + 1
     if( info_io( counter_io + 3 ) == 1 ) num_data_c_trn = num_data_c_trn + 1
   endif
   counter_io = counter_io + 3
 enddo ! i

 call mpi_reduce( natom_total, natom_total_mpi, 1, mpi_integer, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( num_data_total, num_data_total_mpi, 1, mpi_integer, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( num_data_trn, num_data_trn_mpi, 1, mpi_integer, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( num_data_a_trn, num_data_a_trn_mpi, 1, mpi_integer, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( num_data_b_trn, num_data_b_trn_mpi, 1, mpi_integer, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( num_data_c_trn, num_data_c_trn_mpi, 1, mpi_integer, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_bcast( num_data_total_mpi, 1, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )


 ! DISPLS for MPI_GATHERV
 allocate( displs( num_mpi ), recvcount( num_mpi ) )
 displs = 0 ; recvcount = 0
 if( myrank == 0 )then
   m = 0
   do j = 1, num_mpi
     do i = 1, ceiling( dble(counter)/dble(num_mpi) )
       m = j + ( ( i - 1 ) * num_mpi )
       !m = m + 1
       if( m <= counter ) recvcount(j) = recvcount(j) + 1
     enddo
   enddo
   write(*,*) " Data devided into (MPI): "
   write(*,*) recvcount
   write(*,*) " "
   do i = 1, num_mpi
     j = num_mpi + 1 - i
     do k = 1, j - 1
       displs(j) = displs(j) + recvcount(k)
     enddo
   enddo
   displs(1) = 0
 endif ! myrank == 0
 call mpi_bcast( displs, num_mpi, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( recvcount, num_mpi, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )

 allocate( list_mpi( num_data_total_mpi ) )
 call mpi_barrier( mpi_comm_world, ierr )
 call mpi_gatherv( list, num_data_total, mpi_integer, &
                   list_mpi, recvcount, displs, &
                   mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )

 allocate( num_atom_mpi( num_data_total_mpi ) )
 call mpi_barrier( mpi_comm_world, ierr )
 call mpi_gatherv( num_atom, num_data_total, mpi_integer, &
                   num_atom_mpi, recvcount, displs, &
                   mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )

 deallocate( displs, recvcount )


 ! DISPLS for MPI_GATHERV
 allocate( displs( num_mpi ), recvcount( num_mpi ) )
 displs = 0 ; recvcount = 0
 if( myrank == 0 )then
   m = 0
   do j = 1, num_mpi
     do i = 1, ceiling( dble(counter)/dble(num_mpi) )
       m = j + ( ( i - 1 ) * num_mpi )
       !m = m + 1
       if( m <= counter ) recvcount(j) = recvcount(j) + 3
     enddo
   enddo
   do i = 1, num_mpi
     j = num_mpi + 1 - i
     do k = 1, j - 1
       displs(j) = displs(j) + recvcount(k)
     enddo
   enddo
   displs(1) = 0
 endif ! myrank == 0
 call mpi_bcast( displs, num_mpi, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_bcast( recvcount, num_mpi, mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )

 allocate( info_io_mpi( num_data_total_mpi*3 ) )
 call mpi_gatherv( info_io, 3*num_data_total, mpi_integer, &
                   info_io_mpi, recvcount, displs, &
                   mpi_integer, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )

 deallocate( displs, recvcount )


 if( myrank == 0 )then
   num_data_tst_mpi = num_data_total_mpi - num_data_trn_mpi
   write(*,*) " Total number of data: ",    num_data_total_mpi
   write(*,*) " ======================================"
   write(*,*) " Number of training data: ", num_data_trn_mpi
   write(*,*) " Number of test data: ",     num_data_tst_mpi
   write(*,*) " ======================================"
   write(*,*) percent_tst, "% is allocated for test set."
   write(*,*) " "
1112 format(I10,I3,A3,A4,I2,A4,I2,A4,I2,I4)
   open(unit=98,file="list.dat",action="write")
   write(98,*) "# total number of data"
   write(98,*) num_data_total_mpi
   write(98,*) "# 1: training, 0: test, A, B, C 1: on, 0:off"
   counter_io = 0
   do i = 1, num_data_total_mpi
     write(98,1112) i, list_mpi(i), " : ", &
                    "  A ", info_io_mpi(counter_io+1), "  B ", info_io_mpi(counter_io+2),&
                    "  C ", info_io_mpi(counter_io+3), &
                    num_atom_mpi(i)
     counter_io = counter_io + 3
   enddo
   close(98)

   open(unit=98,file="tag_out.dat",action="write")
   do i = 1, num_mpi
     do j = 1, ceiling(dble(counter)/dble(num_mpi))
       infile = i + ( ( j - 1 ) * num_mpi )
       if( infile <= counter )then
         write(98,*) trim(adjustl( tag_total(infile)) )
       endif
     enddo
   enddo
   close(98)

 endif ! myrank == 0


 !***************
 ! L-BFGS-B LOOP
 !***************
 if( myrank == 0 )then
   open(unit=13,file="./mse_training.dat",action="write")
   open(unit=14,file="./mse_test.dat",action="write")
   write(13,*) "# MSE energy (eV^2/atom^2) & force in training (eV^2/ang^2)"
   write(14,*) "# MSE energy (eV^2/atom^2) & force in test (eV^2/ang^2)"
 endif ! myrank == 0

 task = 'START'
 do while( task(1:2) .eq. 'FG' .or. task .eq. 'NEW_X' .or. task .eq. 'START' )
   !---------------
   ! L-BFGS-B code
   !---------------
   if( myrank == 0 )then

     ! Initial weight parameters trial, From scrach
     if( io_weight == 0 .and. weight_counter <= itest )then
       if( weight_counter /= 0 ) write(*,*) " Initial weight test:", weight_counter, f_mpi
       task(1:2) = 'FG'
       weight_counter = weight_counter + 1
       if( weight_counter == 1 ) go to 91
       !if( f_store > f_mpi )then
       !  f_store = f_mpi
       if( f_store > mse_energy_trn_mpi + mse_force_trn_mpi )then
         f_store = mse_energy_trn_mpi + mse_force_trn_mpi
         do i = 1, nparam
           weight_store(i) = weight(i)
         enddo
       endif

       allocate ( ransu_x( (nparam + 2) / 2) )
       allocate ( ransu_y( (nparam + 2) / 2) )
       allocate ( ransu_r( (nparam + 2) / 2) )
       allocate ( ransu_t( (nparam + 2) / 2) )
       allocate ( nransu( nparam) )
       call pre_random
       call random_number( ransu_x )
       call pre_random
       call random_number( ransu_y )

       ransu_r = dsqrt(-2.0 * log(ransu_x))
       ransu_t = 2.0 * pi * ransu_y

       nransu(1:((nparam + 2) / 2)) = ransu_r * dcos(ransu_t)

       do i = ((nparam + 2) / 2) + 1, nparam
         nransu(i) = ransu_r(i - (nparam + 2) / 2) * dsin(ransu_t(i - (nparam + 2) / 2))
       enddo

       if (nlayer_a == 4 .and. nlayer_b == 4 .and. nlayer_c == 4)then
         do i = 1, (network_a(1) + 1) * network_a(2)
           weight(i) = nransu(i) / dsqrt(dble(network_a(1) + 1))
         enddo
         do i = 1, (network_a(2) + 1) * network_a(3)
           weight((network_a(1) + 1) * network_a(2) + 1 + i) = nransu((network_a(1) + 1) * network_a(2) + i) / dsqrt(dble(network_a(2) + 1))
         enddo
         do i = 1, (network_a(3) + 1) * network_a(4)
           weight((network_a(1) + 1) * network_a(2) + (network_a(2) + 1) * network_a(3) + i) = nransu((network_a(1) + 1) * network_a(2) + (network_a(2) + 1) * network_a(3) + i) / dsqrt(dble(network_a(3) + 1))
         enddo
         do i = 1, (network_b(1) + 1) * network_b(2)
           weight(nweight_a + i) = nransu(nweight_a + i) / dsqrt(dble(network_a(1)))
         enddo
         do i = 1, (network_b(2) + 1) * network_b(3)
           weight(nweight_a + (network_b(1) + 1) * network_b(2) + i) = nransu(nweight_a + (network_b(1) + 1)* network_b(2) + i) / dsqrt(dble(network_b(2) + 1))
         enddo
         do i = 1, (network_b(3) + 1) * network_b(4)
           weight(nweight_a + (network_b(1) + 1) * network_b(2) + (network_b(2) + 1) * network_b(3) + i) = nransu(nweight_a + (network_b(1) + 1) * network_b(2) + (network_b(2) + 1) * network_b(3) + i) / dsqrt(dble(network_b(3) + 1))
         enddo
         do i = 1, (network_c(1) + 1) * network_c(2)
           weight(nweight_a + nweight_b + i) = nransu(nweight_a + nweight_b + i) / dsqrt(dble(network_c(1) + 1))
         enddo
         do i = 1, (network_c(2) + 1) * network_c(3)
           weight(nweight_a + nweight_b + (network_c(1) + 1) * network_c(2) + i) = nransu(nweight_a + nweight_b + (network_c(1) + 1) * network_c(2) + i) / dsqrt(dble(network_c(2) + 1))
         enddo
         do i = 1, (network_a(3) + 1) * network_a(4)
           weight(nweight_a + nweight_b + (network_c(1) + 1) * network_c(2) + (network_c(2) + 1)* network_c(3) + i) = nransu(nweight_a + nweight_b + (network_c(1) + 1) * network_c(2) + (network_c(2) + 1)* network_c(3) + i) / dsqrt(dble(network_c(3) + 1))
         enddo
         if (nweight_a + nweight_b + (network_c(1) + 1) * network_c(2) + (network_c(2) + 1)* network_c(3) + (network_a(3) + 1) * network_a(4) .ne. nparam)then
           write(*,*) "nparam error"
         endif
       endif

       deallocate( ransu_x, ransu_y, ransu_r, ransu_t, nransu )

       if( weight_counter <= itest ) go to 91
     elseif( io_weight == 0 .and. weight_counter == itest + 1 )then
       do i = 1, nparam
         weight(i) = weight_store(i)
       enddo
       deallocate( weight_store )
       weight_counter = weight_counter + 1
       task(1:2) = 'START'
     elseif( io_weight == 0 .and. weight_counter == itest + 2 )then
       call setulb( nparam, memory, weight, l, u, nbd, f_mpi, g, factr, pgtol, &
                    wa, iwa, task, iprint, csave, lsave, isave, dsave )
     elseif( io_weight == 1 )then
       call setulb( nparam, memory, weight, l, u, nbd, f_mpi, g, factr, pgtol, &
                    wa, iwa, task, iprint, csave, lsave, isave, dsave )
     endif

   endif ! myrank == 0
91 call mpi_barrier( mpi_comm_world, ierr )
   call mpi_bcast( weight, nparam, mpi_double_precision, 0, mpi_comm_world, ierr )
   call mpi_bcast( task, 60, mpi_character, 0, mpi_comm_world, ierr )
   call mpi_barrier( mpi_comm_world, ierr )

   if( task(1:2) .eq. 'FG' )then

     ! Evaluate f and g
     f = 0.0d0 ; g = 0.0d0
     if( io_nnp_a == 1 ) g_a = 0.0d0
     if( io_nnp_b == 1 ) g_b = 0.0d0
     if( io_nnp_c == 1 ) g_c = 0.0d0

     ! Evaluate RMSE of energy
     mse_energy_trn = 0.0d0 ; mse_energy_tst = 0.0d0
     mse_force_trn  = 0.0d0 ; mse_force_tst  = 0.0d0

     !------------------------------
     ! Distribute weight parameters
     !------------------------------
     !A,B,C
     if( io_nnp_a == 1 .and. io_nnp_b == 1 .and. io_nnp_c == 1 )then
       do i = 1, nweight_a
         weight_a(i) = weight(i)
         weight_inv_a(i) = 0.50d0 / weight_a(i)
       enddo
       dummyi = nweight_a
       do i = 1, nweight_b
         weight_b(i) = weight( i + dummyi )
         weight_inv_b(i) = 0.50d0 / weight_b(i)
       enddo
       dummyi = nweight_a + nweight_b
       do i = 1, nweight_c
         weight_c(i) = weight( i + dummyi )
         weight_inv_c(i) = 0.50d0 / weight_c(i)
       enddo
     endif

     !*************************************************************
     ! DO LOOP MD_STEPS
     !*************************************************************
     counter_io = 0
     if( io_nnp_a == 1 ) counter_a = 0
     if( io_nnp_b == 1 ) counter_b = 0
     if( io_nnp_c == 1 ) counter_c = 0
     counter_natom = 0
     !A,B,C
     if( io_nnp_a == 1 .and. io_nnp_b == 1 .and. io_nnp_c == 1 )then
       counter_natom_a   = 0 ; counter_deriv_a   = 0
       counter_natom_ab  = 0 ; counter_deriv_ab  = 0
       counter_natom_ac  = 0 ; counter_deriv_ac  = 0
       counter_natom_abc = 0 ; counter_deriv_abc = 0
       counter_natom_b   = 0 ; counter_deriv_b   = 0
       counter_natom_ba  = 0 ; counter_deriv_ba  = 0
       counter_natom_bc  = 0 ; counter_deriv_bc  = 0
       counter_natom_bac = 0 ; counter_deriv_bac = 0
       counter_natom_c   = 0 ; counter_deriv_c   = 0
       counter_natom_ca  = 0 ; counter_deriv_ca  = 0
       counter_natom_cb  = 0 ; counter_deriv_cb  = 0
       counter_natom_cab = 0 ; counter_deriv_cab = 0
       ! g2_deriv & g5_deriv
     endif

     do imd = 1, num_data_total

       energy0 = energy0_all( imd )
       io_a = info_io( counter_io + 1 )
       io_b = info_io( counter_io + 2 )
       io_c = info_io( counter_io + 3 )
       counter_io = counter_io + 3

       !---------------
       ! TRAINING DATA
       !---------------
       if( list( imd ) == 1 )then
         !------------------------------------------------------------
         ! Energy(f) &
         ! Diffential of energy with respsct to weight parameters (g)
         !------------------------------------------------------------
         !*** A,B,C NNP only ***
         !A,B,C
         !elseif( io_a == 1 .and. io_b == 1 .and. io_c == 1 )then
         natom_a = 0; natom_b = 0; natom_c = 0
         if( io_a == 1 )then
           counter_a = counter_a + 1
           natom_a   = num_atom_a( counter_a )
           allocate( hidden_a( natom_a, nhidden_a ), atomic_ene_a( natom_a ) )
         endif
         if( io_b == 1 )then
           counter_b = counter_b + 1
           natom_b   = num_atom_b( counter_b )
           allocate( hidden_b( natom_b, nhidden_b ), atomic_ene_b( natom_b ) )
         endif
         if( io_c == 1 )then
           counter_c = counter_c + 1
           natom_c   = num_atom_c( counter_c )
           allocate( hidden_c( natom_c, nhidden_c ), atomic_ene_c( natom_c ) )
         endif

         natom = num_atom( imd )

         allocate( force( natom, 3 ), force0( natom, 3 ) )
         do i = 1, natom
           force0(i,1) = force0_all( counter_natom + 3*(i-1) + 1 )
           force0(i,2) = force0_all( counter_natom + 3*(i-1) + 2 )
           force0(i,3) = force0_all( counter_natom + 3*(i-1) + 3 )
         enddo

         call train3abc( io_a, io_b, io_c, &
                         natom, natom_a, natom_b, natom_c, &
                         natom_a_total,   counter_natom_a,   natom_a_sum,   counter_deriv_a, &
                         natom_ab_total,  counter_natom_ab,  natom_ab_sum,  counter_deriv_ab, &
                         natom_ac_total,  counter_natom_ac,  natom_ac_sum,  counter_deriv_ac, &
                         natom_abc_total, counter_natom_abc, natom_abc_sum, counter_deriv_abc, &
                         natom_b_total,   counter_natom_b,   natom_b_sum,   counter_deriv_b, &
                         natom_ba_total,  counter_natom_ba,  natom_ba_sum,  counter_deriv_ba, &
                         natom_bc_total,  counter_natom_bc,  natom_bc_sum,  counter_deriv_bc, &
                         natom_bac_total, counter_natom_bac, natom_bac_sum, counter_deriv_bac, &
                         natom_c_total,   counter_natom_c,   natom_c_sum,   counter_deriv_c, &
                         natom_ca_total,  counter_natom_ca,  natom_ca_sum,  counter_deriv_ca, &
                         natom_cb_total,  counter_natom_cb,  natom_cb_sum,  counter_deriv_cb, &
                         natom_cab_total, counter_natom_cab, natom_cab_sum, counter_deriv_cab, &
                         g2_a_a,  g2_deriv_a_a,  num_g2_a_a,  g2_a_a_maxmin, &
                         g2_a_b,  g2_deriv_a_b,  num_g2_a_b,  g2_a_b_maxmin, &
                         g2_a_c,  g2_deriv_a_c,  num_g2_a_c,  g2_a_c_maxmin, &
                         g5_a_aa, g5_deriv_a_aa, num_g5_a_aa, g5_a_aa_maxmin, &
                         g5_a_ab, g5_deriv_a_ab, num_g5_a_ab, g5_a_ab_maxmin, &
                         g5_a_ac, g5_deriv_a_ac, num_g5_a_ac, g5_a_ac_maxmin, &
                         g5_a_bb, g5_deriv_a_bb, num_g5_a_bb, g5_a_bb_maxmin, &
                         g5_a_bc, g5_deriv_a_bc, num_g5_a_bc, g5_a_bc_maxmin, &
                         g5_a_cc, g5_deriv_a_cc, num_g5_a_cc, g5_a_cc_maxmin, &
                         network_a, nlayer_a, num_g_a, &
                         weight_a, weight_inv_a, nweight_a, hidden_a, nhidden_a, &
                         g2_b_a,  g2_deriv_b_a,  num_g2_b_a,  g2_b_a_maxmin, &
                         g2_b_b,  g2_deriv_b_b,  num_g2_b_b,  g2_b_b_maxmin, &
                         g2_b_c,  g2_deriv_b_c,  num_g2_b_c,  g2_b_c_maxmin, &
                         g5_b_aa, g5_deriv_b_aa, num_g5_b_aa, g5_b_aa_maxmin, &
                         g5_b_ab, g5_deriv_b_ab, num_g5_b_ab, g5_b_ab_maxmin, &
                         g5_b_ac, g5_deriv_b_ac, num_g5_b_ac, g5_b_ac_maxmin, &
                         g5_b_bb, g5_deriv_b_bb, num_g5_b_bb, g5_b_bb_maxmin, &
                         g5_b_bc, g5_deriv_b_bc, num_g5_b_bc, g5_b_bc_maxmin, &
                         g5_b_cc, g5_deriv_b_cc, num_g5_b_cc, g5_b_cc_maxmin, &
                         network_b, nlayer_b, num_g_b, &
                         weight_b, weight_inv_b, nweight_b, hidden_b, nhidden_b, &
                         g2_c_a,  g2_deriv_c_a,  num_g2_c_a,  g2_c_a_maxmin, &
                         g2_c_b,  g2_deriv_c_b,  num_g2_c_b,  g2_c_b_maxmin, &
                         g2_c_c,  g2_deriv_c_c,  num_g2_c_c,  g2_c_c_maxmin, &
                         g5_c_aa, g5_deriv_c_aa, num_g5_c_aa, g5_c_aa_maxmin, &
                         g5_c_ab, g5_deriv_c_ab, num_g5_c_ab, g5_c_ab_maxmin, &
                         g5_c_ac, g5_deriv_c_ac, num_g5_c_ac, g5_c_ac_maxmin, &
                         g5_c_bb, g5_deriv_c_bb, num_g5_c_bb, g5_c_bb_maxmin, &
                         g5_c_bc, g5_deriv_c_bc, num_g5_c_bc, g5_c_bc_maxmin, &
                         g5_c_cc, g5_deriv_c_cc, num_g5_c_cc, g5_c_cc_maxmin, &
                         network_c, nlayer_c, num_g_c, &
                         weight_c, weight_inv_c, nweight_c, hidden_c, nhidden_c, &
                         beta, g_e_a, g_f_a, g_e_b, g_f_b, g_e_c, g_f_c, &
                         energy, energy0, force, force0 )
         if( io_a == 1 )then
           deallocate( hidden_a, atomic_ene_a )
           g_a =   g_a &
                 + alpha*( g_e_a ) + beta*( g_f_a )
                 !+ alpha*( g_e_a/dble(natom_a) ) + beta*( g_f_a/dble(3*natom_a) )
         endif
         if( io_b == 1 )then
           deallocate( hidden_b, atomic_ene_b )
           g_b =   g_b &
                 + alpha*( g_e_b ) + beta*( g_f_b )
                 !+ alpha*( g_e_b/dble(natom_b) ) + beta*( g_f_b/dble(3*natom_b) )
         endif
         if( io_c == 1 )then
           deallocate( hidden_c, atomic_ene_c )
           g_c =   g_c + &
                 + alpha*( g_e_c ) + beta*( g_f_c )
                 !+ alpha*( g_e_c/dble(natom_c) ) + beta*( g_f_c/dble(3*natom_c) )
         endif
         !-----------
         ! END f & g
         !-----------

         f = f + alpha*( ( energy - energy0 ) )**2
         !f = f + alpha*( ( energy - energy0 )/dble(natom) )**2

         do i = 1, natom
           f = f + beta*( force(i,1)-force0(i,1) )**2
           f = f + beta*( force(i,2)-force0(i,2) )**2
           f = f + beta*( force(i,3)-force0(i,3) )**2
           !f = f + beta*( force(i,1)-force0(i,1) )**2/dble(3*natom)
           !f = f + beta*( force(i,2)-force0(i,2) )**2/dble(3*natom)
           !f = f + beta*( force(i,3)-force0(i,3) )**2/dble(3*natom)
         enddo

         mse_energy_trn =   mse_energy_trn &
                          + ( ( energy - energy0 )/dble(natom) )**2

         do i = 1, natom
           mse_force_trn = mse_force_trn + ( ( force(i,1)-force0(i,1) )**2/dble(3*natom) )
           mse_force_trn = mse_force_trn + ( ( force(i,2)-force0(i,2) )**2/dble(3*natom) )
           mse_force_trn = mse_force_trn + ( ( force(i,3)-force0(i,3) )**2/dble(3*natom) )
         enddo

       !-----------
       ! TEST DATA
       !-----------
       elseif( list( imd ) == 0 )then

         !--------
         ! Energy
         !--------
         !A,B,C
         natom_a = 0; natom_b = 0; natom_c = 0
         if( io_a == 1 )then
           counter_a = counter_a + 1
           natom_a   = num_atom_a( counter_a )
           allocate( hidden_a( natom_a, nhidden_a ), atomic_ene_a( natom_a ) )
         endif
         if( io_b == 1 )then
           counter_b = counter_b + 1
           natom_b   = num_atom_b( counter_b )
           allocate( hidden_b( natom_b, nhidden_b ), atomic_ene_b( natom_b ) )
         endif
         if( io_c == 1 )then
           counter_c = counter_c + 1
           natom_c   = num_atom_c( counter_c )
           allocate( hidden_c( natom_c, nhidden_c ), atomic_ene_c( natom_c ) )
         endif

         natom = num_atom( imd )

         allocate( force( natom, 3 ), force0( natom, 3 ) )
         do i = 1, natom
           force0(i,1) = force0_all( counter_natom + 3*(i-1) + 1 )
           force0(i,2) = force0_all( counter_natom + 3*(i-1) + 2 )
           force0(i,3) = force0_all( counter_natom + 3*(i-1) + 3 )
         enddo

         call valid3abc( io_a, io_b, io_c, &
                         natom, natom_a, natom_b, natom_c, &
                         natom_a_total,   counter_natom_a,   natom_a_sum,   counter_deriv_a, &
                         natom_ab_total,  counter_natom_ab,  natom_ab_sum,  counter_deriv_ab, &
                         natom_ac_total,  counter_natom_ac,  natom_ac_sum,  counter_deriv_ac, &
                         natom_abc_total, counter_natom_abc, natom_abc_sum, counter_deriv_abc, &
                         natom_b_total,   counter_natom_b,   natom_b_sum,   counter_deriv_b, &
                         natom_ba_total,  counter_natom_ba,  natom_ba_sum,  counter_deriv_ba, &
                         natom_bc_total,  counter_natom_bc,  natom_bc_sum,  counter_deriv_bc, &
                         natom_bac_total, counter_natom_bac, natom_bac_sum, counter_deriv_bac, &
                         natom_c_total,   counter_natom_c,   natom_c_sum,   counter_deriv_c, &
                         natom_ca_total,  counter_natom_ca,  natom_ca_sum,  counter_deriv_ca, &
                         natom_cb_total,  counter_natom_cb,  natom_cb_sum,  counter_deriv_cb, &
                         natom_cab_total, counter_natom_cab, natom_cab_sum, counter_deriv_cab, &
                         g2_a_a,  g2_deriv_a_a,  num_g2_a_a,  g2_a_a_maxmin, &
                         g2_a_b,  g2_deriv_a_b,  num_g2_a_b,  g2_a_b_maxmin, &
                         g2_a_c,  g2_deriv_a_c,  num_g2_a_c,  g2_a_c_maxmin, &
                         g5_a_aa, g5_deriv_a_aa, num_g5_a_aa, g5_a_aa_maxmin, &
                         g5_a_ab, g5_deriv_a_ab, num_g5_a_ab, g5_a_ab_maxmin, &
                         g5_a_ac, g5_deriv_a_ac, num_g5_a_ac, g5_a_ac_maxmin, &
                         g5_a_bb, g5_deriv_a_bb, num_g5_a_bb, g5_a_bb_maxmin, &
                         g5_a_bc, g5_deriv_a_bc, num_g5_a_bc, g5_a_bc_maxmin, &
                         g5_a_cc, g5_deriv_a_cc, num_g5_a_cc, g5_a_cc_maxmin, &
                         network_a, nlayer_a, num_g_a, &
                         weight_a, nweight_a, hidden_a, nhidden_a, &
                         g2_b_a,  g2_deriv_b_a,  num_g2_b_a,  g2_b_a_maxmin, &
                         g2_b_b,  g2_deriv_b_b,  num_g2_b_b,  g2_b_b_maxmin, &
                         g2_b_c,  g2_deriv_b_c,  num_g2_b_c,  g2_b_c_maxmin, &
                         g5_b_aa, g5_deriv_b_aa, num_g5_b_aa, g5_b_aa_maxmin, &
                         g5_b_ab, g5_deriv_b_ab, num_g5_b_ab, g5_b_ab_maxmin, &
                         g5_b_ac, g5_deriv_b_ac, num_g5_b_ac, g5_b_ac_maxmin, &
                         g5_b_bb, g5_deriv_b_bb, num_g5_b_bb, g5_b_bb_maxmin, &
                         g5_b_bc, g5_deriv_b_bc, num_g5_b_bc, g5_b_bc_maxmin, &
                         g5_b_cc, g5_deriv_b_cc, num_g5_b_cc, g5_b_cc_maxmin, &
                         network_b, nlayer_b, num_g_b, &
                         weight_b, nweight_b, hidden_b, nhidden_b, &
                         g2_c_a,  g2_deriv_c_a,  num_g2_c_a,  g2_c_a_maxmin, &
                         g2_c_b,  g2_deriv_c_b,  num_g2_c_b,  g2_c_b_maxmin, &
                         g2_c_c,  g2_deriv_c_c,  num_g2_c_c,  g2_c_c_maxmin, &
                         g5_c_aa, g5_deriv_c_aa, num_g5_c_aa, g5_c_aa_maxmin, &
                         g5_c_ab, g5_deriv_c_ab, num_g5_c_ab, g5_c_ab_maxmin, &
                         g5_c_ac, g5_deriv_c_ac, num_g5_c_ac, g5_c_ac_maxmin, &
                         g5_c_bb, g5_deriv_c_bb, num_g5_c_bb, g5_c_bb_maxmin, &
                         g5_c_bc, g5_deriv_c_bc, num_g5_c_bc, g5_c_bc_maxmin, &
                         g5_c_cc, g5_deriv_c_cc, num_g5_c_cc, g5_c_cc_maxmin, &
                         network_c, nlayer_c, num_g_c, &
                         weight_c, nweight_c, hidden_c, nhidden_c, &
                         energy, force )
         if( io_a == 1 ) deallocate( hidden_a, atomic_ene_a )
         if( io_b == 1 ) deallocate( hidden_b, atomic_ene_b )
         if( io_c == 1 ) deallocate( hidden_c, atomic_ene_c )

         mse_energy_tst =   mse_energy_tst &
                          + ( ( energy - energy0 )/dble(natom) )**2

         do i = 1, natom
           mse_force_tst = mse_force_tst + ( ( force(i,1)-force0(i,1) )**2/dble(3*natom) )
           mse_force_tst = mse_force_tst + ( ( force(i,2)-force0(i,2) )**2/dble(3*natom) )
           mse_force_tst = mse_force_tst + ( ( force(i,3)-force0(i,3) )**2/dble(3*natom) )
         enddo

       endif ! list, training or test

       counter_natom = counter_natom + 3*natom
       deallocate( force, force0 )

     enddo ! do imd = 1, num_data_total ! md_steps
     call mpi_barrier( mpi_comm_world, ierr )
     !------------
     ! END do imd
     !------------

     call mpi_reduce( f, f_mpi, 1, mpi_double_precision, mpi_sum, 0, &
                      mpi_comm_world, ierr )
     call mpi_reduce( mse_energy_trn, mse_energy_trn_mpi, 1, mpi_double_precision, &
                      mpi_sum, 0, mpi_comm_world, ierr )
     call mpi_reduce( mse_energy_tst, mse_energy_tst_mpi, 1, mpi_double_precision, &
                      mpi_sum, 0, mpi_comm_world, ierr )
     call mpi_reduce( mse_force_trn, mse_force_trn_mpi, 1, mpi_double_precision, &
                      mpi_sum, 0, mpi_comm_world, ierr )
     call mpi_reduce( mse_force_tst, mse_force_tst_mpi, 1, mpi_double_precision, &
                      mpi_sum, 0, mpi_comm_world, ierr )
     call mpi_barrier( mpi_comm_world, ierr )


     if( myrank == 0 )then

       f_mpi = f_mpi
       !f_mpi = f_mpi / dble( num_data_trn_mpi )

       mse_energy_trn_mpi = mse_energy_trn_mpi / dble( num_data_trn_mpi )
       mse_force_trn_mpi  = mse_force_trn_mpi  / dble( num_data_trn_mpi )
       write(13,fmt='(2F25.12)') mse_energy_trn_mpi, mse_force_trn_mpi ! mse_training

       mse_energy_tst_mpi = mse_energy_tst_mpi / dble( num_data_tst_mpi )
       mse_force_tst_mpi  = mse_force_tst_mpi  / dble( num_data_tst_mpi )
       write(14,fmt='(2F25.12)') mse_energy_tst_mpi, mse_force_tst_mpi ! mse_test

     endif ! myrank == 0

     !----------------------------
     ! Concatenate g_x parameters
     !----------------------------
     call mpi_reduce( g_a, g_a_mpi, nweight_a, mpi_double_precision, &
                      mpi_sum, 0, mpi_comm_world, ierr )
     call mpi_reduce( g_b, g_b_mpi, nweight_b, mpi_double_precision, &
                      mpi_sum, 0, mpi_comm_world, ierr )
     call mpi_reduce( g_c, g_c_mpi, nweight_c, mpi_double_precision, &
                      mpi_sum, 0, mpi_comm_world, ierr )
     call mpi_barrier( mpi_comm_world, ierr )

     if( myrank == 0 )then
       !A,B,C
       if( io_nnp_a == 1 .and. io_nnp_b == 1 .and. io_nnp_c == 1 )then
         if( num_data_a_trn_mpi /= 0 )then
           !g_a_mpi = g_a_mpi / dble( num_data_a_trn_mpi )
           do i = 1, nweight_a
             g(i) = g_a_mpi(i)
           enddo
         endif
         if( num_data_b_trn_mpi /= 0 )then
           !g_b_mpi = g_b_mpi / dble( num_data_b_trn_mpi )
           do i = 1, nweight_b
             g( nweight_a + i ) = g_b_mpi(i)
           enddo
         endif
         if( num_data_c_trn_mpi /= 0 )then
           !g_c_mpi = g_c_mpi / dble( num_data_c_trn_mpi )
           do i = 1, nweight_c
             g( nweight_a + nweight_b + i ) = g_c_mpi(i)
           enddo
         endif
       endif
     endif ! myrank == 0

   endif ! if( task(1:2) .eq. 'FG' )then

   ! Temporary write weight parameters
   if( myrank == 0 )then
     open(unit=7,file="./data_weight/weight.dat",action="write",form="unformatted")
     write(7) weight
     close(7)
   endif ! myrank == 0

 enddo !end of loop do while
 !---------------
 ! END BFGS LOOP
 !---------------


 if( myrank == 0 )then

   mse_energy_trn_mpi = dsqrt( mse_energy_trn_mpi )
   write(*,*)  "# RMSE of training: ", mse_energy_trn_mpi, "eV/atom"
   write(13,*) "# RMSE of training: ", mse_energy_trn_mpi, "eV/atom"
   mse_force_trn_mpi = dsqrt( mse_force_trn_mpi )
   write(*,*)  "# RMSE of training: ", mse_force_trn_mpi, "eV/ang"
   write(13,*) "# RMSE of training: ", mse_force_trn_mpi, "eV/ang"

   mse_energy_tst_mpi = dsqrt( mse_energy_tst_mpi )
   write(*,*)  "# RMSE of test: ", mse_energy_tst_mpi, "eV/atom"
   write(14,*) "# RMSE of test: ", mse_energy_tst_mpi, "eV/atom"
   mse_force_tst_mpi = dsqrt( mse_force_tst_mpi )
   write(*,*)  "# RMSE of test: ", mse_force_tst_mpi, "eV/ang"
   write(14,*) "# RMSE of test: ", mse_force_tst_mpi, "eV/ang"

   close(13) ! mse_trn.dat
   close(14) ! mse_tst.dat

   !-------------------------
   ! Write weight parameters
   !-------------------------
   open(unit=7,file="./data_weight/weight.dat",action="write",form="unformatted")
   write(7) weight
   close(7)

 endif ! myrank == 0


 !***************
 ! FINALIZE DATA
 !***************
 !------------------------------
 ! Distribute weight parameters
 !------------------------------
 !A,B,C
 if( io_nnp_a == 1 .and. io_nnp_b == 1 .and. io_nnp_c == 1 )then
   do i = 1, nweight_a
     weight_a(i) = weight(i)
   enddo
   dummyi = nweight_a
   do i = 1, nweight_b
     weight_b(i) = weight( dummyi + i )
   enddo
   dummyi = nweight_a + nweight_b
   do i = 1, nweight_c
     weight_c(i) = weight( dummyi + i )
   enddo
 endif

 !******************
 ! DO LOOP MD_STEPS
 !******************
 counter_io = 0
 if( io_nnp_a == 1 ) counter_a = 0
 if( io_nnp_b == 1 ) counter_b = 0
 if( io_nnp_c == 1 ) counter_c = 0
 counter_natom = 0
 !A,B,C
 if( io_nnp_a == 1 .and. io_nnp_b == 1 .and. io_nnp_c == 1 )then
   counter_natom_a   = 0 ; counter_deriv_a   = 0
   counter_natom_ab  = 0 ; counter_deriv_ab  = 0
   counter_natom_ac  = 0 ; counter_deriv_ac  = 0
   counter_natom_abc = 0 ; counter_deriv_abc = 0
   counter_natom_b   = 0 ; counter_deriv_b   = 0
   counter_natom_ba  = 0 ; counter_deriv_ba  = 0
   counter_natom_bc  = 0 ; counter_deriv_bc  = 0
   counter_natom_bac = 0 ; counter_deriv_bac = 0
   counter_natom_c   = 0 ; counter_deriv_c   = 0
   counter_natom_ca  = 0 ; counter_deriv_ca  = 0
   counter_natom_cb  = 0 ; counter_deriv_cb  = 0
   counter_natom_cab = 0 ; counter_deriv_cab = 0
   ! g2_deriv & g5_deriv
 endif

 allocate( energy_all( num_data_total ) )
 allocate( force_all( natom_total*3 ) )

 do imd = 1, num_data_total

   energy0 = energy0_all( imd )

   io_a = info_io( counter_io + 1 )
   io_b = info_io( counter_io + 2 )
   io_c = info_io( counter_io + 3 )
   counter_io = counter_io + 3

   !A,B,C
   natom_a = 0; natom_b = 0; natom_c = 0
   if( io_a == 1 )then
     counter_a = counter_a + 1
     natom_a   = num_atom_a( counter_a )
     allocate( hidden_a( natom_a, nhidden_a ), atomic_ene_a( natom_a ) )
   endif
   if( io_b == 1 )then
     counter_b = counter_b + 1
     natom_b   = num_atom_b( counter_b )
     allocate( hidden_b( natom_b, nhidden_b ), atomic_ene_b( natom_b ) )
   endif
   if( io_c == 1 )then
     counter_c = counter_c + 1
     natom_c   = num_atom_c( counter_c )
     allocate( hidden_c( natom_c, nhidden_c ), atomic_ene_c( natom_c ) )
   endif

   natom = num_atom( imd )
   if( natom /= natom_a + natom_b + natom_c ) write(*,*)" natom does not match"
   if( natom /= natom_a + natom_b + natom_c ) stop

   allocate( force( natom, 3 ), force0( natom, 3 ) )
   do i = 1, natom
     force0(i,1) = force0_all( counter_natom + 3*(i-1) + 1 )
     force0(i,2) = force0_all( counter_natom + 3*(i-1) + 2 )
     force0(i,3) = force0_all( counter_natom + 3*(i-1) + 3 )
   enddo

   call valid3abc( io_a, io_b, io_c, &
                   natom, natom_a, natom_b, natom_c, &
                   natom_a_total,   counter_natom_a,   natom_a_sum,   counter_deriv_a, &
                   natom_ab_total,  counter_natom_ab,  natom_ab_sum,  counter_deriv_ab, &
                   natom_ac_total,  counter_natom_ac,  natom_ac_sum,  counter_deriv_ac, &
                   natom_abc_total, counter_natom_abc, natom_abc_sum, counter_deriv_abc, &
                   natom_b_total,   counter_natom_b,   natom_b_sum,   counter_deriv_b, &
                   natom_ba_total,  counter_natom_ba,  natom_ba_sum,  counter_deriv_ba, &
                   natom_bc_total,  counter_natom_bc,  natom_bc_sum,  counter_deriv_bc, &
                   natom_bac_total, counter_natom_bac, natom_bac_sum, counter_deriv_bac, &
                   natom_c_total,   counter_natom_c,   natom_c_sum,   counter_deriv_c, &
                   natom_ca_total,  counter_natom_ca,  natom_ca_sum,  counter_deriv_ca, &
                   natom_cb_total,  counter_natom_cb,  natom_cb_sum,  counter_deriv_cb, &
                   natom_cab_total, counter_natom_cab, natom_cab_sum, counter_deriv_cab, &
                   g2_a_a,  g2_deriv_a_a,  num_g2_a_a,  g2_a_a_maxmin, &
                   g2_a_b,  g2_deriv_a_b,  num_g2_a_b,  g2_a_b_maxmin, &
                   g2_a_c,  g2_deriv_a_c,  num_g2_a_c,  g2_a_c_maxmin, &
                   g5_a_aa, g5_deriv_a_aa, num_g5_a_aa, g5_a_aa_maxmin, &
                   g5_a_ab, g5_deriv_a_ab, num_g5_a_ab, g5_a_ab_maxmin, &
                   g5_a_ac, g5_deriv_a_ac, num_g5_a_ac, g5_a_ac_maxmin, &
                   g5_a_bb, g5_deriv_a_bb, num_g5_a_bb, g5_a_bb_maxmin, &
                   g5_a_bc, g5_deriv_a_bc, num_g5_a_bc, g5_a_bc_maxmin, &
                   g5_a_cc, g5_deriv_a_cc, num_g5_a_cc, g5_a_cc_maxmin, &
                   network_a, nlayer_a, num_g_a, &
                   weight_a, nweight_a, hidden_a, nhidden_a, &
                   g2_b_a,  g2_deriv_b_a,  num_g2_b_a,  g2_b_a_maxmin, &
                   g2_b_b,  g2_deriv_b_b,  num_g2_b_b,  g2_b_b_maxmin, &
                   g2_b_c,  g2_deriv_b_c,  num_g2_b_c,  g2_b_c_maxmin, &
                   g5_b_aa, g5_deriv_b_aa, num_g5_b_aa, g5_b_aa_maxmin, &
                   g5_b_ab, g5_deriv_b_ab, num_g5_b_ab, g5_b_ab_maxmin, &
                   g5_b_ac, g5_deriv_b_ac, num_g5_b_ac, g5_b_ac_maxmin, &
                   g5_b_bb, g5_deriv_b_bb, num_g5_b_bb, g5_b_bb_maxmin, &
                   g5_b_bc, g5_deriv_b_bc, num_g5_b_bc, g5_b_bc_maxmin, &
                   g5_b_cc, g5_deriv_b_cc, num_g5_b_cc, g5_b_cc_maxmin, &
                   network_b, nlayer_b, num_g_b, &
                   weight_b, nweight_b, hidden_b, nhidden_b, &
                   g2_c_a,  g2_deriv_c_a,  num_g2_c_a,  g2_c_a_maxmin, &
                   g2_c_b,  g2_deriv_c_b,  num_g2_c_b,  g2_c_b_maxmin, &
                   g2_c_c,  g2_deriv_c_c,  num_g2_c_c,  g2_c_c_maxmin, &
                   g5_c_aa, g5_deriv_c_aa, num_g5_c_aa, g5_c_aa_maxmin, &
                   g5_c_ab, g5_deriv_c_ab, num_g5_c_ab, g5_c_ab_maxmin, &
                   g5_c_ac, g5_deriv_c_ac, num_g5_c_ac, g5_c_ac_maxmin, &
                   g5_c_bb, g5_deriv_c_bb, num_g5_c_bb, g5_c_bb_maxmin, &
                   g5_c_bc, g5_deriv_c_bc, num_g5_c_bc, g5_c_bc_maxmin, &
                   g5_c_cc, g5_deriv_c_cc, num_g5_c_cc, g5_c_cc_maxmin, &
                   network_c, nlayer_c, num_g_c, &
                   weight_c, nweight_c, hidden_c, nhidden_c, &
                   energy, force )

   if( io_a == 1 ) deallocate( hidden_a, atomic_ene_a )
   if( io_b == 1 ) deallocate( hidden_b, atomic_ene_b )
   if( io_c == 1 ) deallocate( hidden_c, atomic_ene_c )

   energy_all(  imd ) = energy  / dble( natom )
   energy0_all( imd ) = energy0 / dble( natom )

   do i = 1, natom
     force_all( counter_natom + 3*(i-1) + 1 ) = force(i,1)
     force_all( counter_natom + 3*(i-1) + 2 ) = force(i,2)
     force_all( counter_natom + 3*(i-1) + 3 ) = force(i,3)
   enddo

   counter_natom = counter_natom + 3*natom

   deallocate( force, force0 )

 enddo ! do imd = 1, num_data_total
 !-----------------
 ! END DO LOOP imd
 !-----------------
 call mpi_barrier( mpi_comm_world, ierr )


 write(filename,'("energy_training",i4.4,".dat")') myrank
 open(unit=13,file=filename,action="write")
 write(filename,'("energy_test",i4.4,".dat")') myrank
 open(unit=14,file=filename,action="write")
 write(13,*) "# DFT(eV/atom) NNP(eV/atom) Natom"
 write(14,*) "# DFT(eV/atom) NNP(eV/atom) Natom"

 write(filename,'("force_training",i4.4,".dat")') myrank
 open(unit=15,file=filename,action="write")
 write(filename,'("force_test",i4.4,".dat")') myrank
 open(unit=16,file=filename,action="write")
 write(15,*) "# DFT(eV/ang) NNP(eV/ang) Natom"
 write(16,*) "# DFT(eV/ang) NNP(eV/ang) Natom"


 counter_natom = 0

 do imd = 1, num_data_total

   natom = num_atom( imd )

   !---------------
   ! TRAINING DATA
   !---------------
   if( list( imd ) == 1 )then

     write(13,*) energy0_all( imd ), energy_all( imd )
     do i = 1, natom
       write(15,fmt='(6F16.9,I5)') &
       force0_all( counter_natom + 3*(i-1) + 1 ), &
       force0_all( counter_natom + 3*(i-1) + 2 ), &
       force0_all( counter_natom + 3*(i-1) + 3 ), &
       force_all(  counter_natom + 3*(i-1) + 1 ), &
       force_all(  counter_natom + 3*(i-1) + 2 ), &
       force_all(  counter_natom + 3*(i-1) + 3 ), &
       imd
     enddo

   !-----------
   ! TEST DATA
   !-----------
   elseif( list( imd ) == 0 )then

     write(14,*) energy0_all( imd ), energy_all( imd )
     do i = 1, natom
       write(16,fmt='(6F16.9,I5)') &
       force0_all( counter_natom + 3*(i-1) + 1 ), &
       force0_all( counter_natom + 3*(i-1) + 2 ), &
       force0_all( counter_natom + 3*(i-1) + 3 ), &
       force_all(  counter_natom + 3*(i-1) + 1 ), &
       force_all(  counter_natom + 3*(i-1) + 2 ), &
       force_all(  counter_natom + 3*(i-1) + 3 ), &
       imd
     enddo

   else
     write(*,*)" Error(main): Finalize"
     stop
   endif ! list

   counter_natom = counter_natom + 3*natom

 enddo ! do imd = 1, num_data_total
 call mpi_barrier( mpi_comm_world, ierr )

 if( myrank == 0 ) write(13,*) "# RMSE of training: ", mse_energy_trn_mpi, "eV/atom"
 if( myrank == 0 ) write(14,*) "# RMSE of test: ", mse_energy_tst_mpi, "eV/atom"
 if( myrank == 0 ) write(15,*) "# RMSE of training: ", mse_force_trn_mpi, "eV/ang"
 if( myrank == 0 ) write(16,*) "# RMSE of test: ", mse_force_tst_mpi, "eV/ang"
 close(13) ! energy_training.dat
 close(14) ! energy_test.dat
 close(15) ! force_training.dat
 close(16) ! force_test.dat


 !***************************
 ! END FINALIZE TRANING DATA
 !***************************

 if( myrank == 0 )then
   write(*,*) " Finish calculations."
   call cpu_time( t1 )
   write(*,*) " Total time: ", t1-t0, "seconds"
 endif ! myrank == 0

 call mpi_finalize( ierr )





end
