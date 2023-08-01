program main_training
use allocarray
implicit none
integer i, j, k, m, dummyi
character(len=20) dummyc
character work*8,work1*100,work2*100
! Variables for NNP
!-------------------
integer nparam
double precision,allocatable,dimension(:) :: weight
character(len=4) elem_a, elem_b, elem_c
integer io_nnp_a,  io_nnp_b,  io_nnp_c
integer nlayer_a,  nlayer_b,  nlayer_c
integer nweight_a, nweight_b, nweight_c
integer nhidden_a, nhidden_b, nhidden_c
integer,allocatable,dimension(:) :: &
  network_a, network_b, network_c
double precision,allocatable,dimension(:) :: &
  weight_a, weight_b, weight_c
double precision,allocatable,dimension(:) :: &
  atomic_ene_a, atomic_ene_b, atomic_ene_c
double precision,allocatable,dimension(:,:) :: &
  hidden_a, hidden_b, hidden_c
double precision  :: factr, pgtol
integer io_weight, percent_tst
! rc : cutoff radius
double precision rc, rn, eforce_z
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
!A
double precision,allocatable,dimension(:) :: &
  eta2_a_a,  rs2_a_a, &
  eta2_a_b,  rs2_a_b, &
  eta2_a_c,  rs2_a_c, &
  eta5_a_aa, theta5_a_aa, zeta5_a_aa, R_s_a_aa, &
  eta5_a_ab, theta5_a_ab, zeta5_a_ab, R_s_a_ab, &
  eta5_a_ac, theta5_a_ac, zeta5_a_ac, R_s_a_ac, &
  eta5_a_bb, theta5_a_bb, zeta5_a_bb, R_s_a_bb, &
  eta5_a_bc, theta5_a_bc, zeta5_a_bc, R_s_a_bc, &
  eta5_a_cc, theta5_a_cc, zeta5_a_cc, R_s_a_cc
!B
double precision,allocatable,dimension(:) :: &
  eta2_b_a,  rs2_b_a, &
  eta2_b_b,  rs2_b_b, &
  eta2_b_c,  rs2_b_c, &
  eta5_b_aa, theta5_b_aa, zeta5_b_aa, R_s_b_aa, &
  eta5_b_ab, theta5_b_ab, zeta5_b_ab, R_s_b_ab, &
  eta5_b_ac, theta5_b_ac, zeta5_b_ac, R_s_b_ac, &
  eta5_b_bb, theta5_b_bb, zeta5_b_bb, R_s_b_bb, &
  eta5_b_bc, theta5_b_bc, zeta5_b_bc, R_s_b_bc, &
  eta5_b_cc, theta5_b_cc, zeta5_b_cc, R_s_b_cc
!C
double precision,allocatable,dimension(:) :: &
  eta2_c_a,  rs2_c_a, &
  eta2_c_b,  rs2_c_b, &
  eta2_c_c,  rs2_c_c, &
  eta5_c_aa, theta5_c_aa, zeta5_c_aa, R_s_c_aa, &
  eta5_c_ab, theta5_c_ab, zeta5_c_ab, R_s_c_ab, &
  eta5_c_ac, theta5_c_ac, zeta5_c_ac, R_s_c_ac, &
  eta5_c_bb, theta5_c_bb, zeta5_c_bb, R_s_c_bb, &
  eta5_c_bc, theta5_c_bc, zeta5_c_bc, R_s_c_bc, &
  eta5_c_cc, theta5_c_cc, zeta5_c_cc, R_s_c_cc


 !*************************
 ! READ input_nnp.dat FILE
 !*************************
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
        rc, rn, eforce_z, phi, neighbor, &
        io_weight, percent_tst, &
        factr, pgtol )

   num_g_a(1) = num_g2_a_a  ; num_g_a(2) = num_g2_a_b  ; num_g_a(3) = num_g2_a_c
   num_g_a(4) = num_g5_a_aa ; num_g_a(5) = num_g5_a_ab ; num_g_a(6) = num_g5_a_ac
   num_g_a(7) = num_g5_a_bb ; num_g_a(8) = num_g5_a_bc
   num_g_a(9) = num_g5_a_cc

   num_g_b(1) = num_g2_b_a  ; num_g_b(2) = num_g2_b_b  ; num_g_b(3) = num_g2_a_c
   num_g_b(4) = num_g5_b_aa ; num_g_b(5) = num_g5_b_ab ; num_g_b(6) = num_g5_a_ac
   num_g_b(7) = num_g5_b_bb ; num_g_b(8) = num_g5_b_bc
   num_g_b(9) = num_g5_b_cc

   num_g_c(1) = num_g2_c_a  ; num_g_c(2) = num_g2_c_b  ; num_g_c(3) = num_g2_a_c
   num_g_c(4) = num_g5_c_aa ; num_g_c(5) = num_g5_c_ab ; num_g_c(6) = num_g5_a_ac
   num_g_c(7) = num_g5_c_bb ; num_g_c(8) = num_g5_c_bc
   num_g_c(9) = num_g5_c_cc


   !------------------------------
   ! Read number of hidden layers
   !------------------------------
   call read_layer_node3( nlayer_a, nlayer_b, nlayer_c )

   allocate( network_a( nlayer_a ) )
   allocate( network_b( nlayer_b ) )
   allocate( network_c( nlayer_c ) )

   !--------------------------------
   ! Read number of nodes per layer
   !--------------------------------
   call read_layer_node_frex3( network_a, network_b, network_c, &
                               nlayer_a,  nlayer_b,  nlayer_c, &
                               io_nnp_a,  io_nnp_b,  io_nnp_c )

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
   allocate( weight_a( nweight_a ) )
 endif
!B
 if( io_nnp_b == 1 )then
   nweight_b = 0
   do i = 1, nlayer_b - 1
     nweight_b = nweight_b + ( network_b(i) + 1 ) * network_b(i+1)
   enddo
   allocate( weight_b( nweight_b ) )
 endif
!C
 if( io_nnp_c == 1 )then
   nweight_c = 0
   do i = 1, nlayer_c - 1
     nweight_c = nweight_c + ( network_c(i) + 1 ) * network_c(i+1)
   enddo
   allocate( weight_c( nweight_c ) )
 endif


 !----------------------
 ! Number of parameters
 !----------------------
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

 allocate( weight( nparam ) )

 !************************************************************
 ! DEFINE INITIAL WEIGHT PARAMETERS
 !************************************************************
   ! From file
   open(unit=7,file="./data_weight/weight.dat",action="read",form="unformatted")
   read(7) weight
   close(7)

   !A,B,C
   if( io_nnp_a == 1 .and. io_nnp_b == 1 .and. io_nnp_c == 1 )then
     do i = 1, nweight_a
       weight_a(i) = weight( i )
     enddo
     dummyi = nweight_a
     do i = 1, nweight_b
       weight_b(i) = weight( i + dummyi )
     enddo
     dummyi = nweight_a + nweight_b
     do i = 1, nweight_c
       weight_c(i) = weight( i + dummyi )
     enddo
   endif


 !***********
 ! READ DATA
 !***********
 !A,B,C
 if( io_nnp_a == 1 .and. io_nnp_b == 1 .and. io_nnp_c == 1 )then
       call allocate_3( "g2_a_a",  eta2_a_a,  rs2_a_a, &
                        "g2_a_b",  eta2_a_b,  rs2_a_b, &
                        "g2_a_c",  eta2_a_c,  rs2_a_c, &
                        "g5_a_aa", eta5_a_aa, theta5_a_aa, zeta5_a_aa, R_s_a_aa, &
                        "g5_a_ab", eta5_a_ab, theta5_a_ab, zeta5_a_ab, R_s_a_ab, &
                        "g5_a_ac", eta5_a_ac, theta5_a_ac, zeta5_a_ac, R_s_a_ac, &
                        "g5_a_bb", eta5_a_bb, theta5_a_bb, zeta5_a_bb, R_s_a_bb, &
                        "g5_a_bc", eta5_a_bc, theta5_a_bc, zeta5_a_bc, R_s_a_bc, &
                        "g5_a_cc", eta5_a_cc, theta5_a_cc, zeta5_a_cc, R_s_a_cc, &
                        num_g2_a_a,  num_g2_a_b,  num_g2_a_c, &
                        num_g5_a_aa, num_g5_a_ab, num_g5_a_ac, &
                        num_g5_a_bb, num_g5_a_bc, &
                        num_g5_a_cc, &
                        g2_a_a_maxmin,  g2_a_b_maxmin,  g2_a_c_maxmin, &
                        g5_a_aa_maxmin, g5_a_ab_maxmin, g5_a_ac_maxmin, &
                        g5_a_bb_maxmin, g5_a_bc_maxmin, &
                        g5_a_cc_maxmin )
       call allocate_3( "g2_b_a",  eta2_b_a,  rs2_b_a, &
                        "g2_b_b",  eta2_b_b,  rs2_b_b, &
                        "g2_b_c",  eta2_b_c,  rs2_b_c, &
                        "g5_b_aa", eta5_b_aa, theta5_b_aa, zeta5_b_aa, R_s_b_aa, &
                        "g5_b_ab", eta5_b_ab, theta5_b_ab, zeta5_b_ab, R_s_b_ab, &
                        "g5_b_ac", eta5_b_ac, theta5_b_ac, zeta5_b_ac, R_s_b_ac, &
                        "g5_b_bb", eta5_b_bb, theta5_b_bb, zeta5_b_bb, R_s_b_bb, &
                        "g5_b_bc", eta5_b_bc, theta5_b_bc, zeta5_b_bc, R_s_b_bc, &
                        "g5_b_cc", eta5_b_cc, theta5_b_cc, zeta5_b_cc, R_s_b_cc, &
                        num_g2_b_a,  num_g2_b_b,  num_g2_b_c, &
                        num_g5_b_aa, num_g5_b_ab, num_g5_b_ac, &
                        num_g5_b_bb, num_g5_b_bc, &
                        num_g5_b_cc, &
                        g2_b_a_maxmin,  g2_b_b_maxmin,  g2_b_c_maxmin, &
                        g5_b_aa_maxmin, g5_b_ab_maxmin, g5_b_ac_maxmin, &
                        g5_b_bb_maxmin, g5_b_bc_maxmin, &
                        g5_b_cc_maxmin )
       call allocate_3( "g2_c_a",  eta2_c_a,  rs2_c_a, &
                        "g2_c_b",  eta2_c_b,  rs2_c_b, &
                        "g2_c_c",  eta2_c_c,  rs2_c_c, &
                        "g5_c_aa", eta5_c_aa, theta5_c_aa, zeta5_c_aa, R_s_c_aa, &
                        "g5_c_ab", eta5_c_ab, theta5_c_ab, zeta5_c_ab, R_s_c_ab, &
                        "g5_c_ac", eta5_c_ac, theta5_c_ac, zeta5_c_ac, R_s_c_ac, &
                        "g5_c_bb", eta5_c_bb, theta5_c_bb, zeta5_c_bb, R_s_c_bb, &
                        "g5_c_bc", eta5_c_bc, theta5_c_bc, zeta5_c_bc, R_s_c_bc, &
                        "g5_c_cc", eta5_c_cc, theta5_c_cc, zeta5_c_cc, R_s_c_cc, &
                        num_g2_c_a,  num_g2_c_b,  num_g2_c_c, &
                        num_g5_c_aa, num_g5_c_ab, num_g5_c_ac, &
                        num_g5_c_bb, num_g5_c_bc, &
                        num_g5_c_cc, &
                        g2_c_a_maxmin,  g2_c_b_maxmin,  g2_c_c_maxmin, &
                        g5_c_aa_maxmin, g5_c_ab_maxmin, g5_c_ac_maxmin, &
                        g5_c_bb_maxmin, g5_c_bc_maxmin, &
                        g5_c_cc_maxmin )
 endif


 !*************
 ! LAMMPS file
 !*************
     open(unit=99,file="./lammps.nnp",action="write")

     call date_and_time(dummyc)
     write(99,'(2A8)')"# DATE: ",trim(adjustl(dummyc))
     write(99,*)" "
     write(99,'(A15)')"#num_elements 3"
     write(99,'(A16,F10.5)')"cutoff_distance ",rc
     write(99,'(A9,F10.5)')"eforce_z ",eforce_z
     write(99,*)" "
     write(99,'(A8,I5)')"num_g_a ",9
     write(99,*)" "
     write(99,'(A11,I5)')"num_g2_a_a ",  num_g_a(1)
     write(99,'(A11,I5)')"num_g2_a_b ",  num_g_a(2)
     write(99,'(A11,I5)')"num_g2_a_c ",  num_g_a(3)
     write(99,'(A12,I5)')"num_g5_a_aa ", num_g_a(4)
     write(99,'(A12,I5)')"num_g5_a_ab ", num_g_a(5)
     write(99,'(A12,I5)')"num_g5_a_ac ", num_g_a(6)
     write(99,'(A12,I5)')"num_g5_a_bb ", num_g_a(7)
     write(99,'(A12,I5)')"num_g5_a_bc ", num_g_a(8)
     write(99,'(A12,I5)')"num_g5_a_cc ", num_g_a(9)
     write(99,*)" "
     write(99,'(A8,I5)')"num_g_b ",9
     write(99,*)" "
     write(99,'(A11,I5)')"num_g2_b_a ",  num_g_b(1)
     write(99,'(A11,I5)')"num_g2_b_b ",  num_g_b(2)
     write(99,'(A11,I5)')"num_g2_b_c ",  num_g_b(3)
     write(99,'(A12,I5)')"num_g5_b_aa ", num_g_b(4)
     write(99,'(A12,I5)')"num_g5_b_ab ", num_g_b(5)
     write(99,'(A12,I5)')"num_g5_b_ac ", num_g_b(6)
     write(99,'(A12,I5)')"num_g5_b_bb ", num_g_b(7)
     write(99,'(A12,I5)')"num_g5_b_bc ", num_g_b(8)
     write(99,'(A12,I5)')"num_g5_b_cc ", num_g_b(9)
     write(99,*)" "
     write(99,'(A8,I5)')"num_g_c ",9
     write(99,*)" "
     write(99,'(A11,I5)')"num_g2_c_a ",  num_g_c(1)
     write(99,'(A11,I5)')"num_g2_c_b ",  num_g_c(2)
     write(99,'(A11,I5)')"num_g2_c_c ",  num_g_c(3)
     write(99,'(A12,I5)')"num_g5_c_aa ", num_g_c(4)
     write(99,'(A12,I5)')"num_g5_c_ab ", num_g_c(5)
     write(99,'(A12,I5)')"num_g5_c_ac ", num_g_c(6)
     write(99,'(A12,I5)')"num_g5_c_bb ", num_g_c(7)
     write(99,'(A12,I5)')"num_g5_c_bc ", num_g_c(8)
     write(99,'(A12,I5)')"num_g5_c_cc ", num_g_c(9)
     write(99,*)" "
     write(99,'(A9,I5)')"nlayer_a ",   nlayer_a
     write(99,'(A10,I5)')"node_h1_a ",  network_a(2)
     write(99,'(A10,I5)')"node_h2_a ",  network_a(3)
     write(99,'(A11,I5)')"node_out_a ", network_a(4)
     write(99,*)" "
     write(99,'(A9,I5)')"nlayer_b ",   nlayer_b
     write(99,'(A10,I5)')"node_h1_b ",  network_b(2)
     write(99,'(A10,I5)')"node_h2_b ",  network_b(3)
     write(99,'(A11,I5)')"node_out_b ", network_b(4)
     write(99,*)" "
     write(99,'(A9,I5)')"nlayer_c ",   nlayer_c
     write(99,'(A10,I5)')"node_h1_c ",  network_c(2)
     write(99,'(A10,I5)')"node_h2_c ",  network_c(3)
     write(99,'(A11,I5)')"node_out_c ", network_c(4)
     write(99,*)" "


write(*,*)" Writing input parameters..."

write(99,*)" "
write(99,'(A12)')"param_g2_a_a"
do i = 1, num_g2_a_a
  write(99,'(2F15.10)') eta2_a_a(i),  rs2_a_a(i)
enddo
write(99,*)" "
write(99,'(A12)')"param_g2_a_b"
do i = 1, num_g2_a_b
  write(99,'(2F15.10)') eta2_a_b(i),  rs2_a_b(i)
enddo
write(99,*)" "
write(99,'(A12)')"param_g2_a_c"
do i = 1, num_g2_a_c
  write(99,'(2F15.10)') eta2_a_c(i),  rs2_a_c(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_a_aa"
do i = 1, num_g5_a_aa
  write(99,'(4F15.10)') eta5_a_aa(i), theta5_a_aa(i), zeta5_a_aa(i), R_s_a_aa(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_a_ab"
do i = 1, num_g5_a_ab
  write(99,'(4F15.10)') eta5_a_ab(i), theta5_a_ab(i), zeta5_a_ab(i), R_s_a_ab(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_a_ac"
do i = 1, num_g5_a_ac
  write(99,'(4F15.10)') eta5_a_ac(i), theta5_a_ac(i), zeta5_a_ac(i), R_s_a_ac(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_a_bb"
do i = 1, num_g5_a_bb
  write(99,'(4F15.10)') eta5_a_bb(i), theta5_a_bb(i), zeta5_a_bb(i), R_s_a_bb(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_a_bc"
do i = 1, num_g5_a_bc
  write(99,'(4F15.10)') eta5_a_bc(i), theta5_a_bc(i), zeta5_a_bc(i), R_s_a_bc(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_a_cc"
do i = 1, num_g5_a_cc
  write(99,'(4F15.10)') eta5_a_cc(i), theta5_a_cc(i), zeta5_a_cc(i), R_s_a_cc(i)
enddo
write(99,*)" "


write(99,*)" "
write(99,'(A12)')"param_g2_b_a"
do i = 1, num_g2_b_a
  write(99,'(2F15.10)') eta2_b_a(i),  rs2_b_a(i)
enddo
write(99,*)" "
write(99,'(A12)')"param_g2_b_b"
do i = 1, num_g2_b_b
  write(99,'(2F15.10)') eta2_b_b(i),  rs2_b_b(i)
enddo
write(99,*)" "
write(99,'(A12)')"param_g2_b_c"
do i = 1, num_g2_b_c
  write(99,'(2F15.10)') eta2_b_c(i),  rs2_b_c(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_b_aa"
do i = 1, num_g5_b_aa
  write(99,'(4F15.10)') eta5_b_aa(i), theta5_b_aa(i), zeta5_b_aa(i), R_s_b_aa(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_b_ab"
do i = 1, num_g5_b_ab
  write(99,'(4F15.10)') eta5_b_ab(i), theta5_b_ab(i), zeta5_b_ab(i), R_s_b_ab(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_b_ac"
do i = 1, num_g5_b_ac
  write(99,'(4F15.10)') eta5_b_ac(i), theta5_b_ac(i), zeta5_b_ac(i), R_s_b_ac(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_b_bb"
do i = 1, num_g5_b_bb
  write(99,'(4F15.10)') eta5_b_bb(i), theta5_b_bb(i), zeta5_b_bb(i), R_s_b_bb(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_b_bc"
do i = 1, num_g5_b_bc
  write(99,'(4F15.10)') eta5_b_bc(i), theta5_b_bc(i), zeta5_b_bc(i), R_s_b_bc(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_b_cc"
do i = 1, num_g5_b_cc
  write(99,'(4F15.10)') eta5_b_cc(i), theta5_b_cc(i), zeta5_b_cc(i), R_s_b_cc(i)
enddo
write(99,*)" "


write(99,*)" "
write(99,'(A12)')"param_g2_c_a"
do i = 1, num_g2_c_a
  write(99,'(2F15.10)') eta2_c_a(i),  rs2_c_a(i)
enddo
write(99,*)" "
write(99,'(A12)')"param_g2_c_b"
do i = 1, num_g2_c_b
  write(99,'(2F15.10)') eta2_c_b(i),  rs2_c_b(i)
enddo
write(99,*)" "
write(99,'(A12)')"param_g2_c_c"
do i = 1, num_g2_c_c
  write(99,'(2F15.10)') eta2_c_c(i),  rs2_c_c(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_c_aa"
do i = 1, num_g5_c_aa
  write(99,'(4F15.10)') eta5_c_aa(i), theta5_c_aa(i), zeta5_c_aa(i), R_s_c_aa(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_c_ab"
do i = 1, num_g5_c_ab
  write(99,'(4F15.10)') eta5_c_ab(i), theta5_c_ab(i), zeta5_c_ab(i), R_s_c_ab(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_c_ac"
do i = 1, num_g5_c_ac
  write(99,'(4F15.10)') eta5_c_ac(i), theta5_c_ac(i), zeta5_c_ac(i), R_s_c_ac(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_c_bb"
do i = 1, num_g5_c_bb
  write(99,'(4F15.10)') eta5_c_bb(i), theta5_c_bb(i), zeta5_c_bb(i), R_s_c_bb(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_c_bc"
do i = 1, num_g5_c_bc
  write(99,'(4F15.10)') eta5_c_bc(i), theta5_c_bc(i), zeta5_c_bc(i), R_s_c_bc(i)
enddo
write(99,*)" "
write(99,'(A13)')"param_g5_c_cc"
do i = 1, num_g5_c_cc
  write(99,'(4F15.10)') eta5_c_cc(i), theta5_c_cc(i), zeta5_c_cc(i), R_s_c_cc(i)
enddo
write(99,*)" "





write(*,*)" Writing max_min values..."

write(99,*)" "
write(99,'(A10)')"g2_a_a_max"
do i = 1, num_g2_a_a
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_a_a_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo

write(99,*)" "
write(99,'(A10)')"g2_a_a_min"
do i = 1, num_g2_a_a
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_a_a_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo

write(99,*)" "
write(99,'(A10)')"g2_a_b_max"
do i = 1, num_g2_a_b
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_a_b_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo

write(99,*)" "
write(99,'(A10)')"g2_a_b_min"
do i = 1, num_g2_a_b
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_a_b_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo

write(99,*)" "
write(99,'(A10)')"g2_a_c_max"
do i = 1, num_g2_a_c
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_a_c_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo

write(99,*)" "
write(99,'(A10)')"g2_a_c_min"
do i = 1, num_g2_a_c
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_a_c_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo


write(99,*)" "
write(99,'(A11)')"g5_a_aa_max"
do i = 1, num_g5_a_aa
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_a_aa_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_a_aa_min"
do i = 1, num_g5_a_aa
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_a_aa_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo

write(99,*)" "
write(99,'(A11)')"g5_a_ab_max"
do i = 1, num_g5_a_ab
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_a_ab_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_a_ab_min"
do i = 1, num_g5_a_ab
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_a_ab_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo

write(99,*)" "
write(99,'(A11)')"g5_a_ac_max"
do i = 1, num_g5_a_ac
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_a_ac_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_a_ac_min"
do i = 1, num_g5_a_ac
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_a_ac_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo


write(99,*)" "
write(99,'(A11)')"g5_a_bb_max"
do i = 1, num_g5_a_bb
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_a_bb_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_a_bb_min"
do i = 1, num_g5_a_bb
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_a_bb_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo

write(99,*)" "
write(99,'(A11)')"g5_a_bc_max"
do i = 1, num_g5_a_bc
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_a_bc_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_a_bc_min"
do i = 1, num_g5_a_bc
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_a_bc_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo


write(99,*)" "
write(99,'(A11)')"g5_a_cc_max"
do i = 1, num_g5_a_cc
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_a_cc_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_a_cc_min"
do i = 1, num_g5_a_cc
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_a_cc_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo




write(99,*)" "
write(99,'(A10)')"g2_b_a_max"
do i = 1, num_g2_b_a
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_b_a_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo
write(99,*)" "
write(99,'(A10)')"g2_b_a_min"
do i = 1, num_g2_b_a
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_b_a_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo

write(99,*)" "
write(99,'(A10)')"g2_b_b_max"
do i = 1, num_g2_b_b
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_b_b_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo
write(99,*)" "
write(99,'(A10)')"g2_b_b_min"
do i = 1, num_g2_b_b
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_b_b_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo

write(99,*)" "
write(99,'(A10)')"g2_b_c_max"
do i = 1, num_g2_b_c
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_b_c_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo
write(99,*)" "
write(99,'(A10)')"g2_b_c_min"
do i = 1, num_g2_b_c
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_b_c_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo


write(99,*)" "
write(99,'(A11)')"g5_b_aa_max"
do i = 1, num_g5_b_aa
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_b_aa_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_b_aa_min"
do i = 1, num_g5_b_aa
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_b_aa_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo

write(99,*)" "
write(99,'(A11)')"g5_b_ab_max"
do i = 1, num_g5_b_ab
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_b_ab_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_b_ab_min"
do i = 1, num_g5_b_ab
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_b_ab_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo

write(99,*)" "
write(99,'(A11)')"g5_b_ac_max"
do i = 1, num_g5_b_ac
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_b_ac_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_b_ac_min"
do i = 1, num_g5_b_ac
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_b_ac_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo


write(99,*)" "
write(99,'(A11)')"g5_b_bb_max"
do i = 1, num_g5_b_bb
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_b_bb_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_b_bb_min"
do i = 1, num_g5_b_bb
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_b_bb_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo

write(99,*)" "
write(99,'(A11)')"g5_b_bc_max"
do i = 1, num_g5_b_bc
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_b_bc_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_b_bc_min"
do i = 1, num_g5_b_bc
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_b_bc_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo


write(99,*)" "
write(99,'(A11)')"g5_b_cc_max"
do i = 1, num_g5_b_cc
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_b_cc_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_b_cc_min"
do i = 1, num_g5_b_cc
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_b_cc_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo


write(99,*)" "
write(99,'(A10)')"g2_c_a_max"
do i = 1, num_g2_c_a
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_c_a_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo
write(99,*)" "
write(99,'(A10)')"g2_c_a_min"
do i = 1, num_g2_c_a
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_c_a_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo

write(99,*)" "
write(99,'(A10)')"g2_c_b_max"
do i = 1, num_g2_c_b
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_c_b_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo
write(99,*)" "
write(99,'(A10)')"g2_c_b_min"
do i = 1, num_g2_c_b
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_c_b_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo

write(99,*)" "
write(99,'(A10)')"g2_c_c_max"
do i = 1, num_g2_c_c
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_c_c_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo
write(99,*)" "
write(99,'(A10)')"g2_c_c_min"
do i = 1, num_g2_c_c
  write(work,'(i8)')i-1
  write(work1,'(g60.45e3)')g2_c_c_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo


write(99,*)" "
write(99,'(A11)')"g5_c_aa_max"
do i = 1, num_g5_c_aa
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_c_aa_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_c_aa_min"
do i = 1, num_g5_c_aa
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_c_aa_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo

write(99,*)" "
write(99,'(A11)')"g5_c_ab_max"
do i = 1, num_g5_c_ab
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_c_ab_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_c_ab_min"
do i = 1, num_g5_c_ab
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_c_ab_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo

write(99,*)" "
write(99,'(A11)')"g5_c_ac_max"
do i = 1, num_g5_c_ac
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_c_ac_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_c_ac_min"
do i = 1, num_g5_c_ac
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_c_ac_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo


write(99,*)" "
write(99,'(A11)')"g5_c_bb_max"
do i = 1, num_g5_c_bb
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_c_bb_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_c_bb_min"
do i = 1, num_g5_c_bb
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_c_bb_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo

write(99,*)" "
write(99,'(A11)')"g5_c_bc_max"
do i = 1, num_g5_c_bc
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_c_bc_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_c_bc_min"
do i = 1, num_g5_c_bc
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_c_bc_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo


write(99,*)" "
write(99,'(A11)')"g5_c_cc_max"
do i = 1, num_g5_c_cc
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_c_cc_maxmin(i,1)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo
write(99,*)" "
write(99,'(A11)')"g5_c_cc_min"
do i = 1, num_g5_c_cc
  write(work,'(i8)')i-1
  write(work2,'(g60.45e3)')g5_c_cc_maxmin(i,2)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work2))
enddo




write(*,*)" Writing weight parameters..."

write(99,*)" "
write(99,'(A15,I8)')"weight_a_params ", nweight_a
write(99,*)" "
write(99,'(A8)')"weight_a"
do i = 1, nweight_a
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')weight_a(i)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo

write(99,*)" "
write(99,'(A15,I8)')"weight_b_params ", nweight_b
write(99,*)" "
write(99,'(A8)')"weight_b"
do i = 1, nweight_b
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')weight_b(i)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo

write(99,*)" "
write(99,'(A15,I8)')"weight_c_params ", nweight_c
write(99,*)" "
write(99,'(A8)')"weight_c"
do i = 1, nweight_c
  write(work,'(i8)')i-1
  write(work1,'(f50.45)')weight_c(i)
  write(99,'(a," ",$)')trim(adjustl(work)) !dummy(i)
  write(99,'(a)')trim(adjustl(work1))
enddo



write(*,*)" "
write(*,*)" Finish"



end
