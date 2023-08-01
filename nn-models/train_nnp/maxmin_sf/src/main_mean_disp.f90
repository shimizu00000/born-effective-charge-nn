program main_maxmin

implicit none
integer i, j, k
integer stat, counter, time, iter
integer nfile, imd, iparam

integer myid
character(len=120) tag, filename, dir_data, dir_sf
character(len=120) cread
integer            iread
double precision   dread

integer   nelem, natom
integer   natom_a, natom_b, natom_c
integer   natom_1, natom_2, natom_3
integer   io_a, io_b, io_c
integer   io_nnp_a, io_nnp_b, io_nnp_c
character(len=4) elem_a, elem_b, elem_c
character(len=4) elem_1, elem_2, elem_3
integer   md_steps

!Variables for symmetry function
integer node_in_a, node_in_b, node_in_c
!A
integer &
  num_g2_a_a,  num_g2_a_b,  num_g2_a_c, &
  num_g5_a_aa, num_g5_a_ab, num_g5_a_ac, &
  num_g5_a_bb, num_g5_a_bc, &
  num_g5_a_cc
!B
integer &
  num_g2_b_a,  num_g2_b_b,  num_g2_b_c, &
  num_g5_b_aa, num_g5_b_ab, num_g5_b_ac, &
  num_g5_b_bb, num_g5_b_bc, &
  num_g5_b_cc
!C
integer &
  num_g2_c_a,  num_g2_c_b,  num_g2_c_c, &
  num_g5_c_aa, num_g5_c_ab, num_g5_c_ac, &
  num_g5_c_bb, num_g5_c_bc, &
  num_g5_c_cc
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
double precision :: g2_max, g2_min, g5_max, g5_min
!A
integer :: &
  io_a_a,  io_a_b,  io_a_c, &
  io_a_aa, io_a_ab, io_a_ac, &
  io_a_bb, io_a_bc, &
  io_a_cc
!B
integer :: &
  io_b_a,  io_b_b,  io_b_c, &
  io_b_aa, io_b_ab, io_b_ac, &
  io_b_bb, io_b_bc, &
  io_b_cc
!C
integer :: &
  io_c_a,  io_c_b,  io_c_c, &
  io_c_aa, io_c_ab, io_c_ac, &
  io_c_bb, io_c_bc, &
  io_c_cc
integer io_mix

double precision rc, rn
integer neighbor
character(len=10) phi

double precision  t0, t1, t2, t_req
integer           dummyi
double precision  dummy
character(len=20) dummyc


!*************************************************************
! Read input_nnp.dat file
!*************************************************************

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
      dir_data, dir_sf, io_mix )


 call read_io_nnp3( io_nnp_a, io_nnp_b, io_nnp_c )

 write(*,*) " "
 if( io_nnp_a == 1 ) write(*,*) " NNP(a): ON"
 if( io_nnp_b == 1 ) write(*,*) " NNP(b): ON"
 if( io_nnp_c == 1 ) write(*,*) " NNP(c): ON"
 write(*,*) " "


!*************************************************************
! MEAN
!*************************************************************

!--------------------------------
! Allocate G2_maxmin & G5_maxmin
!--------------------------------

 ! 3 elements
 !A,B,C
 if( io_nnp_a == 1 .and. io_nnp_b == 1 .and. io_nnp_c == 1 )then
   !A
   allocate( g2_a_a_maxmin(  num_g2_a_a, 2 ), &
             g2_a_b_maxmin(  num_g2_a_b, 2 ), &
             g2_a_c_maxmin(  num_g2_a_c, 2 ), &
             g5_a_aa_maxmin( num_g5_a_aa, 2 ), &
             g5_a_ab_maxmin( num_g5_a_ab, 2 ), &
             g5_a_ac_maxmin( num_g5_a_ac, 2 ), &
             g5_a_bb_maxmin( num_g5_a_bb, 2 ), &
             g5_a_bc_maxmin( num_g5_a_bc, 2 ), &
             g5_a_cc_maxmin( num_g5_a_cc, 2 ) )
   g2_a_a_maxmin = 0.0d0 
   g2_a_b_maxmin = 0.0d0 
   g2_a_c_maxmin = 0.0d0 
   g5_a_aa_maxmin = 0.0d0 
   g5_a_ab_maxmin = 0.0d0 
   g5_a_ac_maxmin = 0.0d0 
   g5_a_bb_maxmin = 0.0d0 
   g5_a_bc_maxmin = 0.0d0 
   g5_a_cc_maxmin = 0.0d0 
   write(*,*) " Allocate G2_maxmin(a-a,a-b,a-c),&
              &G5_maxmin(a-aa,a-ab,a-ac,a-bb,a-bc,a-cc)"
   !B
   allocate( g2_b_a_maxmin(  num_g2_b_a, 2 ), &
             g2_b_b_maxmin(  num_g2_b_b, 2 ), &
             g2_b_c_maxmin(  num_g2_b_c, 2 ), &
             g5_b_aa_maxmin( num_g5_b_aa, 2 ), &
             g5_b_ab_maxmin( num_g5_b_ab, 2 ), &
             g5_b_ac_maxmin( num_g5_b_ac, 2 ), &
             g5_b_bb_maxmin( num_g5_b_bb, 2 ), &
             g5_b_bc_maxmin( num_g5_b_bc, 2 ), &
             g5_b_cc_maxmin( num_g5_b_cc, 2 ) )
   g2_b_a_maxmin = 0.0d0  
   g2_b_b_maxmin = 0.0d0  
   g2_b_c_maxmin = 0.0d0  
   g5_b_aa_maxmin = 0.0d0 
   g5_b_ab_maxmin = 0.0d0 
   g5_b_ac_maxmin = 0.0d0 
   g5_b_bb_maxmin = 0.0d0 
   g5_b_bc_maxmin = 0.0d0 
   g5_b_cc_maxmin = 0.0d0 
   write(*,*) " Allocate G2_maxmin(b-a,b-b,b-c),&
              &G5_maxmin(b-aa,b-ab,b-ac,b-bb,b-bc,b-cc)"
   !C
   allocate( g2_c_a_maxmin(  num_g2_c_a, 2 ), &
             g2_c_b_maxmin(  num_g2_c_b, 2 ), &
             g2_c_c_maxmin(  num_g2_c_c, 2 ), &
             g5_c_aa_maxmin( num_g5_c_aa, 2 ), &
             g5_c_ab_maxmin( num_g5_c_ab, 2 ), &
             g5_c_ac_maxmin( num_g5_c_ac, 2 ), &
             g5_c_bb_maxmin( num_g5_c_bb, 2 ), &
             g5_c_bc_maxmin( num_g5_c_bc, 2 ), &
             g5_c_cc_maxmin( num_g5_c_cc, 2 ) )
   g2_c_a_maxmin = 0.0d0 
   g2_c_b_maxmin = 0.0d0 
   g2_c_c_maxmin = 0.0d0 
   g5_c_aa_maxmin = 0.0d0 
   g5_c_ab_maxmin = 0.0d0 
   g5_c_ac_maxmin = 0.0d0 
   g5_c_bb_maxmin = 0.0d0 
   g5_c_bc_maxmin = 0.0d0 
   g5_c_cc_maxmin = 0.0d0 
   write(*,*) " Allocate G2_maxmin(c-a,c-b,c-c),&
              &G5_maxmin(c-aa,c-ab,c-ac,c-bb,c-bc,c-cc)"
 endif

 write(*,*) " "


!A
 io_a_a  = 0; io_a_b  = 0; io_a_c  = 0
 io_a_aa = 0; io_a_ab = 0; io_a_ac = 0
 io_a_bb = 0; io_a_bc = 0
 io_a_cc = 0
!B
 io_b_a  = 0; io_b_b  = 0; io_b_c  = 0
 io_b_aa = 0; io_b_ab = 0; io_b_ac = 0
 io_b_bb = 0; io_b_bc = 0
 io_b_cc = 0
!C
 io_c_a  = 0; io_c_b  = 0; io_c_c  = 0
 io_c_aa = 0; io_c_ab = 0; io_c_ac = 0
 io_c_bb = 0; io_c_bc = 0
 io_c_cc = 0

!*************************************************************
! Read input_tag.dat file
!*************************************************************
 open(unit=99,file="input_tag.dat",action="read")
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
 rewind(99)
 !*************************************************************
 ! DO LOOP : nfile (input_tag)
 !*************************************************************
 do nfile = 1, counter
   write(*,*)  " "
   write(*,*)  " # ", nfile
   write(*,*) " ------------------------------------------------------------"
   read(99,*) tag
   write(*,*)  " "
   write(*,*)  " Tag_name: ", trim( tag )
   write(*,*)  " "
1111 format( A120 )
   !Info READ
   write( filename, 1111 ) tag
   !filename = '../step2_data_binary/data/info_'//trim(adjustl(filename))//'.dat'
   filename = trim(adjustl(dir_data))//'/info_'//trim(adjustl(filename))//'.dat'
   open(unit=1,file=filename,action="read",form="unformatted")
   !sf READ
   write( filename, 1111 ) tag
   !write( dummyc, * ) int(rc)
   !filename = '../step3_sf_binary/data_sf/rc_'//trim(adjustl(dummyc))//&
   !            '/sf_'//trim(adjustl(filename))//'.dat'
   filename = trim(adjustl(dir_sf))//'/sf_'//trim(adjustl(filename))//'.dat'
   open(unit=11,file=filename,action="read",form="unformatted")
   !*************************************************************
   ! Read info_tag.dat file (binary)
   !*************************************************************

   call read_io3( natom, nelem, &
                  elem_a,  elem_b,  elem_c, &
                  natom_a, natom_b, natom_c, &
                  io_a,    io_b,    io_c )

   !*************************************************************
   ! Allocate G2 & G5
   !*************************************************************

   ! 1 element
   if( nelem == 1 )then
     !A
     if( io_a == 1 .and. io_b == 0 .and. io_c == 0 )then
       allocate( g2_a_a(  natom_a, num_g2_a_a ), &
                 g5_a_aa( natom_a, num_g5_a_aa ) )
       g2_a_a  = 0.0d0
       g5_a_aa = 0.0d0
       write(*,*) " Allocate G2(a-a),G5(a-aa)"
     !B
     elseif( io_a == 0 .and. io_b == 1 .and. io_c == 0 )then
       allocate( g2_b_b(  natom_b, num_g2_b_b ), &
                 g5_b_bb( natom_b, num_g5_b_bb ) )
       g2_b_b  = 0.0d0
       g5_b_bb = 0.0d0
       write(*,*) " Allocate G2(b-b),G5(b-bb)"
     !C
     elseif( io_a == 0 .and. io_b == 0 .and. io_c == 1 )then
       allocate( g2_c_c(  natom_c, num_g2_c_c ), &
                 g5_c_cc( natom_c, num_g5_c_cc ) )
       g2_c_c  = 0.0d0
       g5_c_cc = 0.0d0
       write(*,*) " Allocate G2(c-c),G5(c-cc)"
     else
       write(*,*) "Error in main : revise program!!"
       stop
     endif
   endif ! 1 element
   ! 2 elements
   if( nelem == 2 )then
     !A,B
     if( io_a == 1 .and. io_b == 1 .and. io_c == 0 )then
       !A
       allocate( g2_a_a(  natom_a, num_g2_a_a ), &
                 g2_a_b(  natom_a, num_g2_a_b ), &
                 g5_a_aa( natom_a, num_g5_a_aa ), &
                 g5_a_ab( natom_a, num_g5_a_ab ), &
                 g5_a_bb( natom_a, num_g5_a_bb ) )
       g2_a_a  = 0.0d0 ; g2_a_b  = 0.0d0
       g5_a_aa = 0.0d0 ; g5_a_ab = 0.0d0
       g5_a_bb = 0.0d0
       write(*,*) " Allocate G2(a-a,a-b),G5(a-aa,a-ab,a-bb)"
       !B
       allocate( g2_b_a(  natom_b, num_g2_b_a ), &
                 g2_b_b(  natom_b, num_g2_b_b ), &
                 g5_b_aa( natom_b, num_g5_b_aa ), &
                 g5_b_ab( natom_b, num_g5_b_ab ), &
                 g5_b_bb( natom_b, num_g5_b_bb ) )
       g2_b_a  = 0.0d0 ; g2_b_b  = 0.0d0
       g5_b_aa = 0.0d0 ; g5_b_ab = 0.0d0
       g5_b_bb = 0.0d0
       write(*,*) " Allocate G2(b-a,b-b),G5(b-aa,b-ab,b-bb)"
     !A,C
     elseif( io_a == 1 .and. io_b == 0 .and. io_c == 1 )then
       !A
       allocate( g2_a_a(  natom_a, num_g2_a_a ), &
                 g2_a_c(  natom_a, num_g2_a_c ), &
                 g5_a_aa( natom_a, num_g5_a_aa ), &
                 g5_a_ac( natom_a, num_g5_a_ac ), &
                 g5_a_cc( natom_a, num_g5_a_cc ) )
       g2_a_a  = 0.0d0 ; g2_a_c  = 0.0d0
       g5_a_aa = 0.0d0 ; g5_a_ac = 0.0d0
       g5_a_cc = 0.0d0
       write(*,*) " Allocate G2(a-a,a-c),G5(a-aa,a-ac,a-cc)"
       !C
       allocate( g2_c_a(  natom_c, num_g2_c_a ), &
                 g2_c_c(  natom_c, num_g2_c_c ), &
                 g5_c_aa( natom_c, num_g5_c_aa ), &
                 g5_c_ac( natom_c, num_g5_c_ac ), &
                 g5_c_cc( natom_c, num_g5_c_cc ) )
       g2_c_a  = 0.0d0 ; g2_c_c  = 0.0d0
       g5_c_aa = 0.0d0 ; g5_c_ac = 0.0d0
       g5_c_cc = 0.0d0
       write(*,*) " Allocate G2(c-a,c-c),G5(c-aa,c-ac,c-cc)"
     !B,C
     elseif( io_a == 0 .and. io_b == 1 .and. io_c == 1 )then
       !B
       allocate( g2_b_b(  natom_b, num_g2_b_b ), &
                 g2_b_c(  natom_b, num_g2_b_c ), &
                 g5_b_bb( natom_b, num_g5_b_bb ), &
                 g5_b_bc( natom_b, num_g5_b_bc ), &
                 g5_b_cc( natom_b, num_g5_b_cc ) )
       g2_b_b  = 0.0d0 ; g2_b_c  = 0.0d0
       g5_b_bb = 0.0d0 ; g5_b_bc = 0.0d0
       g5_b_cc = 0.0d0
       write(*,*) " Allocate G2(b-b,b-c),G5(b-bb,b-bc,b-cc)"
       !C
       allocate( g2_c_b(  natom_c, num_g2_c_b ), &
                 g2_c_c(  natom_c, num_g2_c_c ), &
                 g5_c_bb( natom_c, num_g5_c_bb ), &
                 g5_c_bc( natom_c, num_g5_c_bc ), &
                 g5_c_cc( natom_c, num_g5_c_cc ) )
       g2_c_b  = 0.0d0 ; g2_c_c  = 0.0d0
       g5_c_bb = 0.0d0 ; g5_c_bc = 0.0d0
       g5_c_cc = 0.0d0
       write(*,*) " Allocate G2(c-b,c-c),G5(c-bb,c-bc,c-cc)"
     else
       write(*,*) "Error in main : revise program!!"
       stop
     endif
   endif ! 2 elements
   ! 3 elements
   if( nelem == 3 )then
     !A,B,C
     if( io_a == 1 .and. io_b == 1 .and. io_c == 1 )then
       !A
       allocate( g2_a_a(  natom_a, num_g2_a_a ), &
                 g2_a_b(  natom_a, num_g2_a_b ), &
                 g2_a_c(  natom_a, num_g2_a_c ), &
                 g5_a_aa( natom_a, num_g5_a_aa ), &
                 g5_a_ab( natom_a, num_g5_a_ab ), &
                 g5_a_ac( natom_a, num_g5_a_ac ), &
                 g5_a_bb( natom_a, num_g5_a_bb ), &
                 g5_a_bc( natom_a, num_g5_a_bc ), &
                 g5_a_cc( natom_a, num_g5_a_cc ) )
       g2_a_a  = 0.0d0 ; g2_a_b  = 0.0d0 ; g2_a_c  = 0.0d0
       g5_a_aa = 0.0d0 ; g5_a_ab = 0.0d0 ; g5_a_ac = 0.0d0
       g5_a_bb = 0.0d0 ; g5_a_bc = 0.0d0
       g5_a_cc = 0.0d0
       write(*,*) " Allocate G2(a-a,a-b,a-c),G5(a-aa,a-ab,a-ac,a-bb,a-bc,a-cc)"
       !B
       allocate( g2_b_a(  natom_b, num_g2_b_a ), &
                 g2_b_b(  natom_b, num_g2_b_b ), &
                 g2_b_c(  natom_b, num_g2_b_c ), &
                 g5_b_aa( natom_b, num_g5_b_aa ), &
                 g5_b_ab( natom_b, num_g5_b_ab ), &
                 g5_b_ac( natom_b, num_g5_b_ac ), &
                 g5_b_bb( natom_b, num_g5_b_bb ), &
                 g5_b_bc( natom_b, num_g5_b_bc ), &
                 g5_b_cc( natom_b, num_g5_b_cc ) )
       g2_b_a  = 0.0d0 ; g2_b_b  = 0.0d0 ; g2_b_c  = 0.0d0
       g5_b_aa = 0.0d0 ; g5_b_ab = 0.0d0 ; g5_b_ac = 0.0d0
       g5_b_bb = 0.0d0 ; g5_b_bc = 0.0d0
       g5_b_cc = 0.0d0
       write(*,*) " Allocate G2(b-a,b-b,b-c),G5(b-aa,b-ab,b-ac,b-bb,b-bc,b-cc)"
       !C
       allocate( g2_c_a(  natom_c, num_g2_c_a ), &
                 g2_c_b(  natom_c, num_g2_c_b ), &
                 g2_c_c(  natom_c, num_g2_c_c ), &
                 g5_c_aa( natom_c, num_g5_c_aa ), &
                 g5_c_ab( natom_c, num_g5_c_ab ), &
                 g5_c_ac( natom_c, num_g5_c_ac ), &
                 g5_c_bb( natom_c, num_g5_c_bb ), &
                 g5_c_bc( natom_c, num_g5_c_bc ), &
                 g5_c_cc( natom_c, num_g5_c_cc ) )
       g2_c_a  = 0.0d0 ; g2_c_b  = 0.0d0 ; g2_c_c  = 0.0d0
       g5_c_aa = 0.0d0 ; g5_c_ab = 0.0d0 ; g5_c_ac = 0.0d0
       g5_c_bb = 0.0d0 ; g5_c_bc = 0.0d0
       g5_c_cc = 0.0d0
       write(*,*) " Allocate G2(c-a,c-b,c-c),G5(c-aa,c-ab,c-ac,c-bb,c-bc,c-cc)"
     else
       write(*,*) "Error in main : revise program!!"
       stop
     endif
   endif ! 3 elements
   ! more than 4 elements
   if( nelem > 3 )then
     write(*,*) "Error in main : revise program!!"
     stop
   endif ! more than 4 elements

   !*************************************************************
   ! End Allocate
   !*************************************************************


   write(*,*) " "
   write(*,*) " Ions per type:"
   if( io_a == 1 ) write(*,*) "   a: ", elem_a, natom_a
   if( io_b == 1 ) write(*,*) "   b: ", elem_b, natom_b
   if( io_c == 1 ) write(*,*) "   c: ", elem_c, natom_c

   !MD steps
   read(1) md_steps
   write(*,*) " "
   write(*,*) " MD steps: ", md_steps
   write(*,*) " "


   !*************************************************************
   ! DO LOOP imd (structures)
   !*************************************************************
   do imd = 1, md_steps


     !************************************************************
     ! MEAN
     !************************************************************

     ! 1 element
     if( nelem == 1 )then

       !A
       if( io_a == 1 .and. io_b == 0 .and. io_c == 0 )then

         read(11)g2_a_a
         read(11)g5_a_aa

         if( io_a_a == 0 )then
           do j = 1, num_g2_a_a
             g2_a_a_maxmin(j,1)  = maxval( g2_a_a(:,j) )
             g2_a_a_maxmin(j,2)  = minval( g2_a_a(:,j) )
           enddo
           io_a_a = 1
         else
           do j = 1, num_g2_a_a
             g2_max = maxval( g2_a_a(:,j) )
             g2_min = minval( g2_a_a(:,j) )
             if( g2_max > g2_a_a_maxmin(j,1) ) g2_a_a_maxmin(j,1) = g2_max
             if( g2_min < g2_a_a_maxmin(j,2) ) g2_a_a_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_a_aa == 0 )then
           do j = 1, num_g5_a_aa
             g5_a_aa_maxmin(j,1) = maxval( g5_a_aa(:,j) )
             g5_a_aa_maxmin(j,2) = minval( g5_a_aa(:,j) )
           enddo
           io_a_aa = 1
         else
           do j = 1, num_g5_a_aa
             g5_max = maxval( g5_a_aa(:,j) )
             g5_min = minval( g5_a_aa(:,j) )
             if( g5_max > g5_a_aa_maxmin(j,1) ) g5_a_aa_maxmin(j,1) = g5_max
             if( g5_min < g5_a_aa_maxmin(j,2) ) g5_a_aa_maxmin(j,2) = g5_min
           enddo
         endif

       !B
       elseif( io_a == 0 .and. io_b == 1 .and. io_c == 0 )then

         read(11)g2_b_b
         read(11)g5_b_bb

         if( io_b_b == 0 )then
           do j = 1, num_g2_b_b
             g2_b_b_maxmin(j,1)  = maxval( g2_b_b(:,j) )
             g2_b_b_maxmin(j,2)  = minval( g2_b_b(:,j) )
           enddo
           io_b_b = 1
         else
           do j = 1, num_g2_b_b
             g2_max = maxval( g2_b_b(:,j) )
             g2_min = minval( g2_b_b(:,j) )
             if( g2_max > g2_b_b_maxmin(j,1) ) g2_b_b_maxmin(j,1) = g2_max
             if( g2_min < g2_b_b_maxmin(j,2) ) g2_b_b_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_b_bb == 0 )then
           do j = 1, num_g5_b_bb
             g5_b_bb_maxmin(j,1) = maxval( g5_b_bb(:,j) )
             g5_b_bb_maxmin(j,2) = minval( g5_b_bb(:,j) )
           enddo
           io_b_bb = 1
         else
           do j = 1, num_g5_b_bb
             g5_max = maxval( g5_b_bb(:,j) )
             g5_min = minval( g5_b_bb(:,j) )
             if( g5_max > g5_b_bb_maxmin(j,1) ) g5_b_bb_maxmin(j,1) = g5_max
             if( g5_min < g5_b_bb_maxmin(j,2) ) g5_b_bb_maxmin(j,2) = g5_min
           enddo
         endif

       !C
       elseif( io_a == 0 .and. io_b == 0 .and. io_c == 1 )then

         read(11)g2_c_c
         read(11)g5_c_cc

         if( io_c_c == 0 )then
           do j = 1, num_g2_c_c
             g2_c_c_maxmin(j,1)  = maxval( g2_c_c(:,j) )
             g2_c_c_maxmin(j,2)  = minval( g2_c_c(:,j) )
           enddo
           io_c_c = 1
         else
           do j = 1, num_g2_c_c
             g2_max = maxval( g2_c_c(:,j) )
             g2_min = minval( g2_c_c(:,j) )
             if( g2_max > g2_c_c_maxmin(j,1) ) g2_c_c_maxmin(j,1) = g2_max
             if( g2_min < g2_c_c_maxmin(j,2) ) g2_c_c_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_c_cc == 0 )then
           do j = 1, num_g5_c_cc
             g5_c_cc_maxmin(j,1) = maxval( g5_c_cc(:,j) )
             g5_c_cc_maxmin(j,2) = minval( g5_c_cc(:,j) )
           enddo
           io_c_cc = 1
           do j = 1, num_g5_c_cc
             g5_max = maxval( g5_c_cc(:,j) )
             g5_min = minval( g5_c_cc(:,j) )
             if( g5_max > g5_c_cc_maxmin(j,1) ) g5_c_cc_maxmin(j,1) = g5_max
             if( g5_min < g5_c_cc_maxmin(j,2) ) g5_c_cc_maxmin(j,2) = g5_min
           enddo
         endif

       else
         write(*,*) "Error in main-6 : revise program!!"
         stop
       endif

     endif ! 1 element

     ! 2 elements
     if( nelem == 2 )then   

       !A,B
       if( io_a == 1 .and. io_b == 1 .and. io_c == 0 )then

         !A
         read(11)g2_a_a
         read(11)g2_a_b
         read(11)g5_a_aa
         read(11)g5_a_ab
         read(11)g5_a_bb

         if( io_a_a == 0 )then
           do j = 1, num_g2_a_a
             g2_a_a_maxmin(j,1)  = maxval( g2_a_a(:,j) )
             g2_a_a_maxmin(j,2)  = minval( g2_a_a(:,j) )
           enddo
           io_a_a = 1
         else
           do j = 1, num_g2_a_a
             g2_max = maxval( g2_a_a(:,j) )
             g2_min = minval( g2_a_a(:,j) )
             if( g2_max > g2_a_a_maxmin(j,1) ) g2_a_a_maxmin(j,1) = g2_max
             if( g2_min < g2_a_a_maxmin(j,2) ) g2_a_a_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_a_b == 0 )then
           do j = 1, num_g2_a_b
             g2_a_b_maxmin(j,1)  = maxval( g2_a_b(:,j) )
             g2_a_b_maxmin(j,2)  = minval( g2_a_b(:,j) )
           enddo
           io_a_b = 1
         else
           do j = 1, num_g2_a_b
             g2_max = maxval( g2_a_b(:,j) )
             g2_min = minval( g2_a_b(:,j) )
             if( g2_max > g2_a_b_maxmin(j,1) ) g2_a_b_maxmin(j,1) = g2_max
             if( g2_min < g2_a_b_maxmin(j,2) ) g2_a_b_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_a_aa == 0 )then
           do j = 1, num_g5_a_aa
             g5_a_aa_maxmin(j,1) = maxval( g5_a_aa(:,j) )
             g5_a_aa_maxmin(j,2) = minval( g5_a_aa(:,j) )
           enddo
           io_a_aa = 1
         else
           do j = 1, num_g5_a_aa
             g5_max = maxval( g5_a_aa(:,j) )
             g5_min = minval( g5_a_aa(:,j) )
             if( g5_max > g5_a_aa_maxmin(j,1) ) g5_a_aa_maxmin(j,1) = g5_max
             if( g5_min < g5_a_aa_maxmin(j,2) ) g5_a_aa_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_a_ab == 0 )then
           do j = 1, num_g5_a_ab
             g5_a_ab_maxmin(j,1) = maxval( g5_a_ab(:,j) )
             g5_a_ab_maxmin(j,2) = minval( g5_a_ab(:,j) )
           enddo
           io_a_ab = 1
         else
           do j = 1, num_g5_a_ab
             g5_max = maxval( g5_a_ab(:,j) )
             g5_min = minval( g5_a_ab(:,j) )
             if( g5_max > g5_a_ab_maxmin(j,1) ) g5_a_ab_maxmin(j,1) = g5_max
             if( g5_min < g5_a_ab_maxmin(j,2) ) g5_a_ab_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_a_bb == 0 )then
           do j = 1, num_g5_a_bb
             g5_a_bb_maxmin(j,1) = maxval( g5_a_bb(:,j) )
             g5_a_bb_maxmin(j,2) = minval( g5_a_bb(:,j) )
           enddo
           io_a_bb = 1
         else
           do j = 1, num_g5_a_bb
             g5_max = maxval( g5_a_bb(:,j) )
             g5_min = minval( g5_a_bb(:,j) )
             if( g5_max > g5_a_bb_maxmin(j,1) ) g5_a_bb_maxmin(j,1) = g5_max
             if( g5_min < g5_a_bb_maxmin(j,2) ) g5_a_bb_maxmin(j,2) = g5_min
           enddo
         endif

         !B
         read(11)g2_b_a 
         read(11)g2_b_b
         read(11)g5_b_aa 
         read(11)g5_b_ab
         read(11)g5_b_bb

         if( io_b_a == 0 )then
           do j = 1, num_g2_b_a
             g2_b_a_maxmin(j,1)  = maxval( g2_b_a(:,j) )
             g2_b_a_maxmin(j,2)  = minval( g2_b_a(:,j) )
           enddo
           io_b_a = 1
         else
           do j = 1, num_g2_b_a
             g2_max = maxval( g2_b_a(:,j) )
             g2_min = minval( g2_b_a(:,j) )
             if( g2_max > g2_b_a_maxmin(j,1) ) g2_b_a_maxmin(j,1) = g2_max
             if( g2_min < g2_b_a_maxmin(j,2) ) g2_b_a_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_b_b == 0 )then
           do j = 1, num_g2_b_b
             g2_b_b_maxmin(j,1)  = maxval( g2_b_b(:,j) )
             g2_b_b_maxmin(j,2)  = minval( g2_b_b(:,j) )
           enddo
           io_b_b = 1
         else
           do j = 1, num_g2_b_b
             g2_max = maxval( g2_b_b(:,j) )
             g2_min = minval( g2_b_b(:,j) )
             if( g2_max > g2_b_b_maxmin(j,1) ) g2_b_b_maxmin(j,1) = g2_max
             if( g2_min < g2_b_b_maxmin(j,2) ) g2_b_b_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_b_aa == 0 )then
           do j = 1, num_g5_b_aa
             g5_b_aa_maxmin(j,1) = maxval( g5_b_aa(:,j) )
             g5_b_aa_maxmin(j,2) = minval( g5_b_aa(:,j) )
           enddo
           io_b_aa = 1
         else
           do j = 1, num_g5_b_aa
             g5_max = maxval( g5_b_aa(:,j) )
             g5_min = minval( g5_b_aa(:,j) )
             if( g5_max > g5_b_aa_maxmin(j,1) ) g5_b_aa_maxmin(j,1) = g5_max
             if( g5_min < g5_b_aa_maxmin(j,2) ) g5_b_aa_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_b_ab == 0 )then
           do j = 1, num_g5_b_ab
             g5_b_ab_maxmin(j,1) = maxval( g5_b_ab(:,j) )
             g5_b_ab_maxmin(j,2) = minval( g5_b_ab(:,j) )
           enddo
           io_b_ab = 1
         else
           do j = 1, num_g5_b_ab
             g5_max = maxval( g5_b_ab(:,j) )
             g5_min = minval( g5_b_ab(:,j) )
             if( g5_max > g5_b_ab_maxmin(j,1) ) g5_b_ab_maxmin(j,1) = g5_max
             if( g5_min < g5_b_ab_maxmin(j,2) ) g5_b_ab_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_b_bb == 0 )then
           do j = 1, num_g5_b_bb
             g5_b_bb_maxmin(j,1) = maxval( g5_b_bb(:,j) )
             g5_b_bb_maxmin(j,2) = minval( g5_b_bb(:,j) )
           enddo
           io_b_bb = 1
         else
           do j = 1, num_g5_b_bb
             g5_max = maxval( g5_b_bb(:,j) )
             g5_min = minval( g5_b_bb(:,j) )
             if( g5_max > g5_b_bb_maxmin(j,1) ) g5_b_bb_maxmin(j,1) = g5_max
             if( g5_min < g5_b_bb_maxmin(j,2) ) g5_b_bb_maxmin(j,2) = g5_min
           enddo
         endif

       !A,C
       elseif( io_a == 1 .and. io_b == 0 .and. io_c == 1 )then

         !A
         read(11)g2_a_a
         read(11)g2_a_c
         read(11)g5_a_aa
         read(11)g5_a_ac
         read(11)g5_a_cc

         if( io_a_a == 0 )then
           do j = 1, num_g2_a_a
             g2_a_a_maxmin(j,1)  = maxval( g2_a_a(:,j) )
             g2_a_a_maxmin(j,2)  = minval( g2_a_a(:,j) )
           enddo
           io_a_a = 1
         else
           do j = 1, num_g2_a_a
             g2_max = maxval( g2_a_a(:,j) )
             g2_min = minval( g2_a_a(:,j) )
             if( g2_max > g2_a_a_maxmin(j,1) ) g2_a_a_maxmin(j,1) = g2_max
             if( g2_min < g2_a_a_maxmin(j,2) ) g2_a_a_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_a_c == 0 )then
           do j = 1, num_g2_a_c
             g2_a_c_maxmin(j,1)  = maxval( g2_a_c(:,j) )
             g2_a_c_maxmin(j,2)  = minval( g2_a_c(:,j) )
           enddo
           io_a_c = 1
         else
           do j = 1, num_g2_a_c
             g2_max = maxval( g2_a_c(:,j) )
             g2_min = minval( g2_a_c(:,j) )
             if( g2_max > g2_a_c_maxmin(j,1) ) g2_a_c_maxmin(j,1) = g2_max
             if( g2_min < g2_a_c_maxmin(j,2) ) g2_a_c_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_a_aa == 0 )then
           do j = 1, num_g5_a_aa
             g5_a_aa_maxmin(j,1) = maxval( g5_a_aa(:,j) )
             g5_a_aa_maxmin(j,2) = minval( g5_a_aa(:,j) )
           enddo
           io_a_aa = 1
         else
           do j = 1, num_g5_a_aa
             g5_max = maxval( g5_a_aa(:,j) )
             g5_min = minval( g5_a_aa(:,j) )
             if( g5_max > g5_a_aa_maxmin(j,1) ) g5_a_aa_maxmin(j,1) = g5_max
             if( g5_min < g5_a_aa_maxmin(j,2) ) g5_a_aa_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_a_ac == 0 )then
           do j = 1, num_g5_a_ac
             g5_a_ac_maxmin(j,1) = maxval( g5_a_ac(:,j) )
             g5_a_ac_maxmin(j,2) = minval( g5_a_ac(:,j) )
           enddo
           io_a_ac = 1
         else
           do j = 1, num_g5_a_ac
             g5_max = maxval( g5_a_ac(:,j) )
             g5_min = minval( g5_a_ac(:,j) )
             if( g5_max > g5_a_ac_maxmin(j,1) ) g5_a_ac_maxmin(j,1) = g5_max
             if( g5_min < g5_a_ac_maxmin(j,2) ) g5_a_ac_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_a_cc == 0 )then
           do j = 1, num_g5_a_cc
             g5_a_cc_maxmin(j,1) = maxval( g5_a_cc(:,j) )
             g5_a_cc_maxmin(j,2) = minval( g5_a_cc(:,j) )
           enddo
           io_a_cc = 1
         else
           do j = 1, num_g5_a_cc
             g5_max = maxval( g5_a_cc(:,j) )
             g5_min = minval( g5_a_cc(:,j) )
             if( g5_max > g5_a_cc_maxmin(j,1) ) g5_a_cc_maxmin(j,1) = g5_max
             if( g5_min < g5_a_cc_maxmin(j,2) ) g5_a_cc_maxmin(j,2) = g5_min
           enddo
         endif

         !C
         read(11)g2_c_a
         read(11)g2_c_c
         read(11)g5_c_aa
         read(11)g5_c_ac
         read(11)g5_c_cc

         if( io_c_a == 0 )then
           do j = 1, num_g2_c_a
             g2_c_a_maxmin(j,1)  = maxval( g2_c_a(:,j) )
             g2_c_a_maxmin(j,2)  = minval( g2_c_a(:,j) )
           enddo
           io_c_a = 1
         else
           do j = 1, num_g2_c_a
             g2_max = maxval( g2_c_a(:,j) )
             g2_min = minval( g2_c_a(:,j) )
             if( g2_max > g2_c_a_maxmin(j,1) ) g2_c_a_maxmin(j,1) = g2_max
             if( g2_min < g2_c_a_maxmin(j,2) ) g2_c_a_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_c_c == 0 )then
           do j = 1, num_g2_c_c
             g2_c_c_maxmin(j,1)  = maxval( g2_c_c(:,j) )
             g2_c_c_maxmin(j,2)  = minval( g2_c_c(:,j) )
           enddo
           io_c_c = 1
         else
           do j = 1, num_g2_c_c
             g2_max = maxval( g2_c_c(:,j) )
             g2_min = minval( g2_c_c(:,j) )
             if( g2_max > g2_c_c_maxmin(j,1) ) g2_c_c_maxmin(j,1) = g2_max
             if( g2_min < g2_c_c_maxmin(j,2) ) g2_c_c_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_c_aa == 0 )then
           do j = 1, num_g5_c_aa
             g5_c_aa_maxmin(j,1) = maxval( g5_c_aa(:,j) )
             g5_c_aa_maxmin(j,2) = minval( g5_c_aa(:,j) )
           enddo
           io_c_aa = 1
         else
           do j = 1, num_g5_c_aa
             g5_max = maxval( g5_c_aa(:,j) )
             g5_min = minval( g5_c_aa(:,j) )
             if( g5_max > g5_c_aa_maxmin(j,1) ) g5_c_aa_maxmin(j,1) = g5_max
             if( g5_min < g5_c_aa_maxmin(j,2) ) g5_c_aa_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_c_ac == 0 )then
           do j = 1, num_g5_c_ac
             g5_c_ac_maxmin(j,1) = maxval( g5_c_ac(:,j) )
             g5_c_ac_maxmin(j,2) = minval( g5_c_ac(:,j) )
           enddo
           io_c_ac = 1
         else
           do j = 1, num_g5_c_ac
             g5_max = maxval( g5_c_ac(:,j) )
             g5_min = minval( g5_c_ac(:,j) )
             if( g5_max > g5_c_ac_maxmin(j,1) ) g5_c_ac_maxmin(j,1) = g5_max
             if( g5_min < g5_c_ac_maxmin(j,2) ) g5_c_ac_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_c_cc == 0 )then
           do j = 1, num_g5_c_cc
             g5_c_cc_maxmin(j,1) = maxval( g5_c_cc(:,j) )
             g5_c_cc_maxmin(j,2) = minval( g5_c_cc(:,j) )
           enddo
           io_c_cc = 1
         else
           do j = 1, num_g5_c_cc
             g5_max = maxval( g5_c_cc(:,j) )
             g5_min = minval( g5_c_cc(:,j) )
             if( g5_max > g5_c_cc_maxmin(j,1) ) g5_c_cc_maxmin(j,1) = g5_max
             if( g5_min < g5_c_cc_maxmin(j,2) ) g5_c_cc_maxmin(j,2) = g5_min
           enddo
         endif

       !B,C
       elseif( io_a == 0 .and. io_b == 1 .and. io_c == 1 )then

         !B
         read(11)g2_b_b
         read(11)g2_b_c
         read(11)g5_b_bb
         read(11)g5_b_bc
         read(11)g5_b_cc

         if( io_b_b == 0 )then
           do j = 1, num_g2_b_b
             g2_b_b_maxmin(j,1)  = maxval( g2_b_b(:,j) )
             g2_b_b_maxmin(j,2)  = minval( g2_b_b(:,j) )
           enddo
           io_b_b = 1
         else
           do j = 1, num_g2_b_b
             g2_max = maxval( g2_b_b(:,j) )
             g2_min = minval( g2_b_b(:,j) )
             if( g2_max > g2_b_b_maxmin(j,1) ) g2_b_b_maxmin(j,1) = g2_max
             if( g2_min < g2_b_b_maxmin(j,2) ) g2_b_b_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_b_c == 0 )then
           do j = 1, num_g2_b_c
             g2_b_c_maxmin(j,1)  = maxval( g2_b_c(:,j) )
             g2_b_c_maxmin(j,2)  = minval( g2_b_c(:,j) )
           enddo
           io_b_c = 1
         else
           do j = 1, num_g2_b_c
             g2_max = maxval( g2_b_c(:,j) )
             g2_min = minval( g2_b_c(:,j) )
             if( g2_max > g2_b_c_maxmin(j,1) ) g2_b_c_maxmin(j,1) = g2_max
             if( g2_min < g2_b_c_maxmin(j,2) ) g2_b_c_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_b_bb == 0 )then
           do j = 1, num_g5_b_bb
             g5_b_bb_maxmin(j,1) = maxval( g5_b_bb(:,j) )
             g5_b_bb_maxmin(j,2) = minval( g5_b_bb(:,j) )
           enddo
           io_b_bb = 1
         else
           do j = 1, num_g5_b_bb
             g5_max = maxval( g5_b_bb(:,j) )
             g5_min = minval( g5_b_bb(:,j) )
             if( g5_max > g5_b_bb_maxmin(j,1) ) g5_b_bb_maxmin(j,1) = g5_max
             if( g5_min < g5_b_bb_maxmin(j,2) ) g5_b_bb_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_b_bc == 0 )then
           do j = 1, num_g5_b_bc
             g5_b_bc_maxmin(j,1) = maxval( g5_b_bc(:,j) )
             g5_b_bc_maxmin(j,2) = minval( g5_b_bc(:,j) )
           enddo
           io_b_bc = 1
         else
           do j = 1, num_g5_b_bc
             g5_max = maxval( g5_b_bc(:,j) )
             g5_min = minval( g5_b_bc(:,j) )
             if( g5_max > g5_b_bc_maxmin(j,1) ) g5_b_bc_maxmin(j,1) = g5_max
             if( g5_min < g5_b_bc_maxmin(j,2) ) g5_b_bc_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_b_cc == 0 )then
           do j = 1, num_g5_b_cc
             g5_b_cc_maxmin(j,1) = maxval( g5_b_cc(:,j) )
             g5_b_cc_maxmin(j,2) = minval( g5_b_cc(:,j) )
           enddo
           io_b_cc = 1
         else
           do j = 1, num_g5_b_cc
             g5_max = maxval( g5_b_cc(:,j) )
             g5_min = minval( g5_b_cc(:,j) )
             if( g5_max > g5_b_cc_maxmin(j,1) ) g5_b_cc_maxmin(j,1) = g5_max
             if( g5_min < g5_b_cc_maxmin(j,2) ) g5_b_cc_maxmin(j,2) = g5_min
           enddo
         endif

         !C
         read(11)g2_c_b  
         read(11)g2_c_c
         read(11)g5_c_bb
         read(11)g5_c_bc
         read(11)g5_c_cc

         if( io_c_b == 0 )then
           do j = 1, num_g2_c_b
             g2_c_b_maxmin(j,1)  = maxval( g2_c_b(:,j) )
             g2_c_b_maxmin(j,2)  = minval( g2_c_b(:,j) )
           enddo
           io_c_b = 1
         else
           do j = 1, num_g2_c_b
             g2_max = maxval( g2_c_b(:,j) )
             g2_min = minval( g2_c_b(:,j) )
             if( g2_max > g2_c_b_maxmin(j,1) ) g2_c_b_maxmin(j,1) = g2_max
             if( g2_min < g2_c_b_maxmin(j,2) ) g2_c_b_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_c_c == 0 )then
           do j = 1, num_g2_c_c
             g2_c_c_maxmin(j,1)  = maxval( g2_c_c(:,j) )
             g2_c_c_maxmin(j,2)  = minval( g2_c_c(:,j) )
           enddo
           io_c_c = 1
         else
           do j = 1, num_g2_c_c
             g2_max = maxval( g2_c_c(:,j) )
             g2_min = minval( g2_c_c(:,j) )
             if( g2_max > g2_c_c_maxmin(j,1) ) g2_c_c_maxmin(j,1) = g2_max
             if( g2_min < g2_c_c_maxmin(j,2) ) g2_c_c_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_c_bb == 0 )then
           do j = 1, num_g5_c_bb
             g5_c_bb_maxmin(j,1) = maxval( g5_c_bb(:,j) )
             g5_c_bb_maxmin(j,2) = minval( g5_c_bb(:,j) )
           enddo
           io_c_bb = 1
         else
           do j = 1, num_g5_c_bb
             g5_max = maxval( g5_c_bb(:,j) )
             g5_min = minval( g5_c_bb(:,j) )
             if( g5_max > g5_c_bb_maxmin(j,1) ) g5_c_bb_maxmin(j,1) = g5_max
             if( g5_min < g5_c_bb_maxmin(j,2) ) g5_c_bb_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_c_bc == 0 )then
           do j = 1, num_g5_c_bc
             g5_c_bc_maxmin(j,1) = maxval( g5_c_bc(:,j) )
             g5_c_bc_maxmin(j,2) = minval( g5_c_bc(:,j) )
           enddo
           io_c_bc = 1
         else
           do j = 1, num_g5_c_bc
             g5_max = maxval( g5_c_bc(:,j) )
             g5_min = minval( g5_c_bc(:,j) )
             if( g5_max > g5_c_bc_maxmin(j,1) ) g5_c_bc_maxmin(j,1) = g5_max
             if( g5_min < g5_c_bc_maxmin(j,2) ) g5_c_bc_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_c_cc == 0 )then
           do j = 1, num_g5_c_cc
             g5_c_cc_maxmin(j,1) = maxval( g5_c_cc(:,j) )
             g5_c_cc_maxmin(j,2) = minval( g5_c_cc(:,j) )
           enddo
           io_c_cc = 1
         else
           do j = 1, num_g5_c_cc
             g5_max = maxval( g5_c_cc(:,j) )
             g5_min = minval( g5_c_cc(:,j) )
             if( g5_max > g5_c_cc_maxmin(j,1) ) g5_c_cc_maxmin(j,1) = g5_max
             if( g5_min < g5_c_cc_maxmin(j,2) ) g5_c_cc_maxmin(j,2) = g5_min
           enddo
         endif

       else
         write(*,*) "Error in main-7 : revise program!!"
         stop
       endif

     endif ! 2 elements

     ! 3 elements
     if( nelem == 3 )then

       !A,B,C
       if( io_a == 1 .and. io_b == 1 .and. io_c == 1 )then

         !A
         read(11)g2_a_a
         read(11)g2_a_b  
         read(11)g2_a_c
         read(11)g5_a_aa
         read(11)g5_a_ab 
         read(11)g5_a_ac
         read(11)g5_a_bb 
         read(11)g5_a_bc
         read(11)g5_a_cc

         if( io_a_a == 0 )then
           do j = 1, num_g2_a_a
             g2_a_a_maxmin(j,1)  = maxval( g2_a_a(:,j) )
             g2_a_a_maxmin(j,2)  = minval( g2_a_a(:,j) )
           enddo
           io_a_a = 1
         else
           do j = 1, num_g2_a_a
             g2_max = maxval( g2_a_a(:,j) )
             g2_min = minval( g2_a_a(:,j) )
             if( g2_max > g2_a_a_maxmin(j,1) ) g2_a_a_maxmin(j,1) = g2_max
             if( g2_min < g2_a_a_maxmin(j,2) ) g2_a_a_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_a_b == 0 )then
           do j = 1, num_g2_a_b
             g2_a_b_maxmin(j,1)  = maxval( g2_a_b(:,j) )
             g2_a_b_maxmin(j,2)  = minval( g2_a_b(:,j) )
           enddo
           io_a_b = 1
         else
           do j = 1, num_g2_a_b
             g2_max = maxval( g2_a_b(:,j) )
             g2_min = minval( g2_a_b(:,j) )
             if( g2_max > g2_a_b_maxmin(j,1) ) g2_a_b_maxmin(j,1) = g2_max
             if( g2_min < g2_a_b_maxmin(j,2) ) g2_a_b_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_a_c == 0 )then
           do j = 1, num_g2_a_c
             g2_a_c_maxmin(j,1)  = maxval( g2_a_c(:,j) )
             g2_a_c_maxmin(j,2)  = minval( g2_a_c(:,j) )
           enddo
           io_a_c = 1
         else
           do j = 1, num_g2_a_c
             g2_max = maxval( g2_a_c(:,j) )
             g2_min = minval( g2_a_c(:,j) )
             if( g2_max > g2_a_c_maxmin(j,1) ) g2_a_c_maxmin(j,1) = g2_max
             if( g2_min < g2_a_c_maxmin(j,2) ) g2_a_c_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_a_aa == 0 )then
           do j = 1, num_g5_a_aa
             g5_a_aa_maxmin(j,1) = maxval( g5_a_aa(:,j) )
             g5_a_aa_maxmin(j,2) = minval( g5_a_aa(:,j) )
           enddo
           io_a_aa = 1
         else
           do j = 1, num_g5_a_aa
             g5_max = maxval( g5_a_aa(:,j) )
             g5_min = minval( g5_a_aa(:,j) )
             if( g5_max > g5_a_aa_maxmin(j,1) ) g5_a_aa_maxmin(j,1) = g5_max
             if( g5_min < g5_a_aa_maxmin(j,2) ) g5_a_aa_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_a_ab == 0 )then
           do j = 1, num_g5_a_ab
             g5_a_ab_maxmin(j,1) = maxval( g5_a_ab(:,j) )
             g5_a_ab_maxmin(j,2) = minval( g5_a_ab(:,j) )
           enddo
           io_a_ab = 1
         else
           do j = 1, num_g5_a_ab
             g5_max = maxval( g5_a_ab(:,j) )
             g5_min = minval( g5_a_ab(:,j) )
             if( g5_max > g5_a_ab_maxmin(j,1) ) g5_a_ab_maxmin(j,1) = g5_max
             if( g5_min < g5_a_ab_maxmin(j,2) ) g5_a_ab_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_a_ac == 0 )then
           do j = 1, num_g5_a_ac
             g5_a_ac_maxmin(j,1) = maxval( g5_a_ac(:,j) )
             g5_a_ac_maxmin(j,2) = minval( g5_a_ac(:,j) )
           enddo
           io_a_ac = 1
         else
           do j = 1, num_g5_a_ac
             g5_max = maxval( g5_a_ac(:,j) )
             g5_min = minval( g5_a_ac(:,j) )
             if( g5_max > g5_a_ac_maxmin(j,1) ) g5_a_ac_maxmin(j,1) = g5_max
             if( g5_min < g5_a_ac_maxmin(j,2) ) g5_a_ac_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_a_bb == 0 )then
           do j = 1, num_g5_a_bb
             g5_a_bb_maxmin(j,1) = maxval( g5_a_bb(:,j) )
             g5_a_bb_maxmin(j,2) = minval( g5_a_bb(:,j) )
           enddo
           io_a_bb = 1
         else
           do j = 1, num_g5_a_bb
             g5_max = maxval( g5_a_bb(:,j) )
             g5_min = minval( g5_a_bb(:,j) )
             if( g5_max > g5_a_bb_maxmin(j,1) ) g5_a_bb_maxmin(j,1) = g5_max
             if( g5_min < g5_a_bb_maxmin(j,2) ) g5_a_bb_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_a_bc == 0 )then
           do j = 1, num_g5_a_bc
             g5_a_bc_maxmin(j,1) = maxval( g5_a_bc(:,j) )
             g5_a_bc_maxmin(j,2) = minval( g5_a_bc(:,j) )
           enddo
           io_a_bc = 1
         else
           do j = 1, num_g5_a_bc
             g5_max = maxval( g5_a_bc(:,j) )
             g5_min = minval( g5_a_bc(:,j) )
             if( g5_max > g5_a_bc_maxmin(j,1) ) g5_a_bc_maxmin(j,1) = g5_max
             if( g5_min < g5_a_bc_maxmin(j,2) ) g5_a_bc_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_a_cc == 0 )then
           do j = 1, num_g5_a_cc
             g5_a_cc_maxmin(j,1) = maxval( g5_a_cc(:,j) )
             g5_a_cc_maxmin(j,2) = minval( g5_a_cc(:,j) )
           enddo
           io_a_cc = 1
         else
           do j = 1, num_g5_a_cc
             g5_max = maxval( g5_a_cc(:,j) )
             g5_min = minval( g5_a_cc(:,j) )
             if( g5_max > g5_a_cc_maxmin(j,1) ) g5_a_cc_maxmin(j,1) = g5_max
             if( g5_min < g5_a_cc_maxmin(j,2) ) g5_a_cc_maxmin(j,2) = g5_min
           enddo
         endif

         !B
         read(11)g2_b_a
         read(11)g2_b_b
         read(11)g2_b_c
         read(11)g5_b_aa
         read(11)g5_b_ab
         read(11)g5_b_ac
         read(11)g5_b_bb
         read(11)g5_b_bc
         read(11)g5_b_cc

         if( io_b_a == 0 )then
           do j = 1, num_g2_b_a
             g2_b_a_maxmin(j,1)  = maxval( g2_b_a(:,j) )
             g2_b_a_maxmin(j,2)  = minval( g2_b_a(:,j) )
           enddo
           io_b_a = 1
         else
           do j = 1, num_g2_b_a
             g2_max = maxval( g2_b_a(:,j) )
             g2_min = minval( g2_b_a(:,j) )
             if( g2_max > g2_b_a_maxmin(j,1) ) g2_b_a_maxmin(j,1) = g2_max
             if( g2_min < g2_b_a_maxmin(j,2) ) g2_b_a_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_b_b == 0 )then
           do j = 1, num_g2_b_b
             g2_b_b_maxmin(j,1)  = maxval( g2_b_b(:,j) )
             g2_b_b_maxmin(j,2)  = minval( g2_b_b(:,j) )
           enddo
           io_b_b = 1
         else
           do j = 1, num_g2_b_b
             g2_max = maxval( g2_b_b(:,j) )
             g2_min = minval( g2_b_b(:,j) )
             if( g2_max > g2_b_b_maxmin(j,1) ) g2_b_b_maxmin(j,1) = g2_max
             if( g2_min < g2_b_b_maxmin(j,2) ) g2_b_b_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_b_c == 0 )then
           do j = 1, num_g2_b_c
             g2_b_c_maxmin(j,1)  = maxval( g2_b_c(:,j) )
             g2_b_c_maxmin(j,2)  = minval( g2_b_c(:,j) )
           enddo
           io_b_c = 1
         else
           do j = 1, num_g2_b_c
             g2_max = maxval( g2_b_c(:,j) )
             g2_min = minval( g2_b_c(:,j) )
             if( g2_max > g2_b_c_maxmin(j,1) ) g2_b_c_maxmin(j,1) = g2_max
             if( g2_min < g2_b_c_maxmin(j,2) ) g2_b_c_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_b_aa == 0 )then
           do j = 1, num_g5_b_aa
             g5_b_aa_maxmin(j,1) = maxval( g5_b_aa(:,j) )
             g5_b_aa_maxmin(j,2) = minval( g5_b_aa(:,j) )
           enddo
           io_b_aa = 1
         else
           do j = 1, num_g5_b_aa
             g5_max = maxval( g5_b_aa(:,j) )
             g5_min = minval( g5_b_aa(:,j) )
             if( g5_max > g5_b_aa_maxmin(j,1) ) g5_b_aa_maxmin(j,1) = g5_max
             if( g5_min < g5_b_aa_maxmin(j,2) ) g5_b_aa_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_b_ab == 0 )then
           do j = 1, num_g5_b_ab
             g5_b_ab_maxmin(j,1) = maxval( g5_b_ab(:,j) )
             g5_b_ab_maxmin(j,2) = minval( g5_b_ab(:,j) )
           enddo
           io_b_ab = 1
         else
           do j = 1, num_g5_b_ab
             g5_max = maxval( g5_b_ab(:,j) )
             g5_min = minval( g5_b_ab(:,j) )
             if( g5_max > g5_b_ab_maxmin(j,1) ) g5_b_ab_maxmin(j,1) = g5_max
             if( g5_min < g5_b_ab_maxmin(j,2) ) g5_b_ab_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_b_ac == 0 )then
           do j = 1, num_g5_b_ac
             g5_b_ac_maxmin(j,1) = maxval( g5_b_ac(:,j) )
             g5_b_ac_maxmin(j,2) = minval( g5_b_ac(:,j) )
           enddo
           io_b_ac = 1
         else
           do j = 1, num_g5_b_ac
             g5_max = maxval( g5_b_ac(:,j) )
             g5_min = minval( g5_b_ac(:,j) )
             if( g5_max > g5_b_ac_maxmin(j,1) ) g5_b_ac_maxmin(j,1) = g5_max
             if( g5_min < g5_b_ac_maxmin(j,2) ) g5_b_ac_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_b_bb == 0 )then
           do j = 1, num_g5_b_bb
             g5_b_bb_maxmin(j,1) = maxval( g5_b_bb(:,j) )
             g5_b_bb_maxmin(j,2) = minval( g5_b_bb(:,j) )
           enddo
           io_b_bb = 1
         else
           do j = 1, num_g5_b_bb
             g5_max = maxval( g5_b_bb(:,j) )
             g5_min = minval( g5_b_bb(:,j) )
             if( g5_max > g5_b_bb_maxmin(j,1) ) g5_b_bb_maxmin(j,1) = g5_max
             if( g5_min < g5_b_bb_maxmin(j,2) ) g5_b_bb_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_b_bc == 0 )then
           do j = 1, num_g5_b_bc
             g5_b_bc_maxmin(j,1) = maxval( g5_b_bc(:,j) )
             g5_b_bc_maxmin(j,2) = minval( g5_b_bc(:,j) )
           enddo
           io_b_bc = 1
         else
           do j = 1, num_g5_b_bc
             g5_max = maxval( g5_b_bc(:,j) )
             g5_min = minval( g5_b_bc(:,j) )
             if( g5_max > g5_b_bc_maxmin(j,1) ) g5_b_bc_maxmin(j,1) = g5_max
             if( g5_min < g5_b_bc_maxmin(j,2) ) g5_b_bc_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_b_cc == 0 )then
           do j = 1, num_g5_b_cc
             g5_b_cc_maxmin(j,1) = maxval( g5_b_cc(:,j) )
             g5_b_cc_maxmin(j,2) = minval( g5_b_cc(:,j) )
           enddo
           io_b_cc = 1
         else
           do j = 1, num_g5_b_cc
             g5_max = maxval( g5_b_cc(:,j) )
             g5_min = minval( g5_b_cc(:,j) )
             if( g5_max > g5_b_cc_maxmin(j,1) ) g5_b_cc_maxmin(j,1) = g5_max
             if( g5_min < g5_b_cc_maxmin(j,2) ) g5_b_cc_maxmin(j,2) = g5_min
           enddo
         endif

         !C
         read(11)g2_c_a
         read(11)g2_c_b
         read(11)g2_c_c
         read(11)g5_c_aa
         read(11)g5_c_ab
         read(11)g5_c_ac
         read(11)g5_c_bb
         read(11)g5_c_bc
         read(11)g5_c_cc

         if( io_c_a == 0 )then
           do j = 1, num_g2_c_a
             g2_c_a_maxmin(j,1)  = maxval( g2_c_a(:,j) )
             g2_c_a_maxmin(j,2)  = minval( g2_c_a(:,j) )
           enddo
           io_c_a = 1
         else
           do j = 1, num_g2_c_a
             g2_max = maxval( g2_c_a(:,j) )
             g2_min = minval( g2_c_a(:,j) )
             if( g2_max > g2_c_a_maxmin(j,1) ) g2_c_a_maxmin(j,1) = g2_max
             if( g2_min < g2_c_a_maxmin(j,2) ) g2_c_a_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_c_b == 0 )then
           do j = 1, num_g2_c_b
             g2_c_b_maxmin(j,1)  = maxval( g2_c_b(:,j) )
             g2_c_b_maxmin(j,2)  = minval( g2_c_b(:,j) )
           enddo
           io_c_b = 1
         else
           do j = 1, num_g2_c_b
             g2_max = maxval( g2_c_b(:,j) )
             g2_min = minval( g2_c_b(:,j) )
             if( g2_max > g2_c_b_maxmin(j,1) ) g2_c_b_maxmin(j,1) = g2_max
             if( g2_min < g2_c_b_maxmin(j,2) ) g2_c_b_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_c_c == 0 )then
           do j = 1, num_g2_c_c
             g2_c_c_maxmin(j,1)  = maxval( g2_c_c(:,j) )
             g2_c_c_maxmin(j,2)  = minval( g2_c_c(:,j) )
           enddo
           io_c_c = 1
         else
           do j = 1, num_g2_c_c
             g2_max = maxval( g2_c_c(:,j) )
             g2_min = minval( g2_c_c(:,j) )
             if( g2_max > g2_c_c_maxmin(j,1) ) g2_c_c_maxmin(j,1) = g2_max
             if( g2_min < g2_c_c_maxmin(j,2) ) g2_c_c_maxmin(j,2) = g2_min
           enddo
         endif

         if( io_c_aa == 0 )then
           do j = 1, num_g5_c_aa
             g5_c_aa_maxmin(j,1) = maxval( g5_c_aa(:,j) )
             g5_c_aa_maxmin(j,2) = minval( g5_c_aa(:,j) )
           enddo
           io_c_aa = 1
         else
           do j = 1, num_g5_c_aa
             g5_max = maxval( g5_c_aa(:,j) )
             g5_min = minval( g5_c_aa(:,j) )
             if( g5_max > g5_c_aa_maxmin(j,1) ) g5_c_aa_maxmin(j,1) = g5_max
             if( g5_min < g5_c_aa_maxmin(j,2) ) g5_c_aa_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_c_ab == 0 )then
           do j = 1, num_g5_c_ab
             g5_c_ab_maxmin(j,1) = maxval( g5_c_ab(:,j) )
             g5_c_ab_maxmin(j,2) = minval( g5_c_ab(:,j) )
           enddo
           io_c_ab = 1
         else
           do j = 1, num_g5_c_ab
             g5_max = maxval( g5_c_ab(:,j) )
             g5_min = minval( g5_c_ab(:,j) )
             if( g5_max > g5_c_ab_maxmin(j,1) ) g5_c_ab_maxmin(j,1) = g5_max
             if( g5_min < g5_c_ab_maxmin(j,2) ) g5_c_ab_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_c_ac == 0 )then
           do j = 1, num_g5_c_ac
             g5_c_ac_maxmin(j,1) = maxval( g5_c_ac(:,j) )
             g5_c_ac_maxmin(j,2) = minval( g5_c_ac(:,j) )
           enddo
           io_c_ac = 1
         else
           do j = 1, num_g5_c_ac
             g5_max = maxval( g5_c_ac(:,j) )
             g5_min = minval( g5_c_ac(:,j) )
             if( g5_max > g5_c_ac_maxmin(j,1) ) g5_c_ac_maxmin(j,1) = g5_max
             if( g5_min < g5_c_ac_maxmin(j,2) ) g5_c_ac_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_c_bb == 0 )then
           do j = 1, num_g5_c_bb
             g5_c_bb_maxmin(j,1) = maxval( g5_c_bb(:,j) )
             g5_c_bb_maxmin(j,2) = minval( g5_c_bb(:,j) )
           enddo
           io_c_bb = 1
         else
           do j = 1, num_g5_c_bb
             g5_max = maxval( g5_c_bb(:,j) )
             g5_min = minval( g5_c_bb(:,j) )
             if( g5_max > g5_c_bb_maxmin(j,1) ) g5_c_bb_maxmin(j,1) = g5_max
             if( g5_min < g5_c_bb_maxmin(j,2) ) g5_c_bb_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_c_bc == 0 )then
           do j = 1, num_g5_c_bc
             g5_c_bc_maxmin(j,1) = maxval( g5_c_bc(:,j) )
             g5_c_bc_maxmin(j,2) = minval( g5_c_bc(:,j) )
           enddo
           io_c_bc = 1
         else
           do j = 1, num_g5_c_bc
             g5_max = maxval( g5_c_bc(:,j) )
             g5_min = minval( g5_c_bc(:,j) )
             if( g5_max > g5_c_bc_maxmin(j,1) ) g5_c_bc_maxmin(j,1) = g5_max
             if( g5_min < g5_c_bc_maxmin(j,2) ) g5_c_bc_maxmin(j,2) = g5_min
           enddo
         endif

         if( io_c_cc == 0 )then
           do j = 1, num_g5_c_cc
             g5_c_cc_maxmin(j,1) = maxval( g5_c_cc(:,j) )
             g5_c_cc_maxmin(j,2) = minval( g5_c_cc(:,j) )
           enddo
           io_c_cc = 1
         else
           do j = 1, num_g5_c_cc
             g5_max = maxval( g5_c_cc(:,j) )
             g5_min = minval( g5_c_cc(:,j) )
             if( g5_max > g5_c_cc_maxmin(j,1) ) g5_c_cc_maxmin(j,1) = g5_max
             if( g5_min < g5_c_cc_maxmin(j,2) ) g5_c_cc_maxmin(j,2) = g5_min
           enddo
         endif

       else
         write(*,*) "Error in main-8 : revise program!!"
         stop
       endif

     endif

     ! more than 4 elements
     if( nelem > 3 )then
       write(*,*) "Error in main-9 : revise program !!"
       stop
     endif


   enddo ! imd


   close(1) ! info
   close(11) ! sf_...


   !************************************************************
   ! Deallocate
   !************************************************************

   ! 1 element
   if( nelem == 1 )then
     !A
     if( io_a == 1 .and. io_b == 0 .and. io_c == 0 )then
       deallocate( g2_a_a, g5_a_aa )
       write(*,*) " Deallocate G2(a-a),G5(a-aa)"
     !B
     elseif( io_a == 0 .and. io_b == 1 .and. io_c == 0 )then
       deallocate( g2_b_b, g5_b_bb )
     !C
     elseif( io_a == 0 .and. io_b == 0 .and. io_c == 1 )then
       deallocate( g2_c_c, g5_c_cc )
     else
       write(*,*) "Error in main : revise program!!"
       stop
     endif
   endif
   ! 2 elements
   if( nelem == 2 )then
     !A,B
     if( io_a == 1 .and. io_b == 1 .and. io_c == 0 )then
       deallocate( g2_a_a,  g2_a_b, &
                   g5_a_aa, g5_a_ab, g5_a_bb )
       deallocate( g2_b_a,  g2_b_b, &
                   g5_b_aa, g5_b_ab, g5_b_bb )
     !A,C
     elseif( io_a == 1 .and. io_b == 0 .and. io_c == 1 )then
       deallocate( g2_a_a,  g2_a_c, &
                   g5_a_aa, g5_a_ac, g5_a_cc )
       deallocate( g2_c_a,  g2_c_c, &
                   g5_c_aa, g5_c_ac, g5_c_cc )
     !B,C
     elseif( io_a == 0 .and. io_b == 1 .and. io_c == 1 )then
       deallocate( g2_b_b,  g2_b_c, &
                   g5_b_bb, g5_b_bc, g5_b_cc )
       deallocate( g2_c_b,  g2_c_c, &
                   g5_c_bb, g5_c_bc, g5_c_cc )
     else
       write(*,*) "Error in main : revise program!!"
       stop
     endif
   endif
   ! 3 elements
   if( nelem == 3 )then
     !A,B,C
     if( io_a == 1 .and. io_b == 1 .and. io_c == 1 )then
       deallocate( g2_a_a,  g2_a_b,  g2_a_c, &
                   g5_a_aa, g5_a_ab, g5_a_ac, &
                   g5_a_bb, g5_a_bc, &
                   g5_a_cc )
       deallocate( g2_b_a,  g2_b_b,  g2_b_c, &
                   g5_b_aa, g5_b_ab, g5_b_ac, &
                   g5_b_bb, g5_b_bc, &
                   g5_b_cc )
       deallocate( g2_c_a,  g2_c_b,  g2_c_c, &
                   g5_c_aa, g5_c_ab, g5_c_ac, &
                   g5_c_bb, g5_c_bc, &
                   g5_c_cc )
     else
       write(*,*) "Error in main : revise program!!"
       stop
     endif
   endif
   ! more than 4 elements
   if( nelem > 3 )then
     write(*,*) "Error in main : revise program!!"
     stop
   endif


 enddo ! nfile


 close(99) ! input_tag.dat


 !if( io_mix == 1 )then

   !g2_a_b_maxmin(:,2) = 0.0d0
   !g2_a_c_maxmin(:,2) = 0.0d0
   !g5_a_ab_maxmin(:,2) = 0.0d0
   !g5_a_ac_maxmin(:,2) = 0.0d0
   !g5_a_bb_maxmin(:,2) = 0.0d0
   !g5_a_bc_maxmin(:,2) = 0.0d0
   !g5_a_cc_maxmin(:,2) = 0.0d0

   !g2_b_a_maxmin(:,2) = 0.0d0
   !g2_b_c_maxmin(:,2) = 0.0d0
   !g5_b_aa_maxmin(:,2) = 0.0d0
   !g5_b_ab_maxmin(:,2) = 0.0d0
   !g5_b_ac_maxmin(:,2) = 0.0d0
   !g5_b_bc_maxmin(:,2) = 0.0d0
   !g5_b_cc_maxmin(:,2) = 0.0d0

   !g2_c_a_maxmin(:,2) = 0.0d0
   !g2_c_b_maxmin(:,2) = 0.0d0
   !g5_c_aa_maxmin(:,2) = 0.0d0
   !g5_c_ab_maxmin(:,2) = 0.0d0
   !g5_c_ac_maxmin(:,2) = 0.0d0
   !g5_c_bb_maxmin(:,2) = 0.0d0
   !g5_c_bc_maxmin(:,2) = 0.0d0
   !g5_c_cc_maxmin(:,2) = 0.0d0

 !endif


 !************************************************************
 ! WRITE MEAN
 !************************************************************

 write(*,*) " "
 write(*,*) " Writing maxmin values"
 write(*,*) " "

 !A
 if( io_nnp_a == 1 )then
   open(unit=20,file="./data_maxmin/g2_a_a_maxmin.dat",action="write",form="unformatted")
   write(20)g2_a_a_maxmin
   close(20)
 endif
 if( io_nnp_a == 1 .and. io_nnp_b == 1 )then
   open(unit=20,file="./data_maxmin/g2_a_b_maxmin.dat",action="write",form="unformatted")
   write(20)g2_a_b_maxmin
   close(20)
 endif
 if( io_nnp_a == 1 .and. io_nnp_c == 1 )then
   open(unit=20,file="./data_maxmin/g2_a_c_maxmin.dat",action="write",form="unformatted")
   write(20)g2_a_c_maxmin
   close(20)
 endif
 if( io_nnp_a == 1 )then
   open(unit=20,file="./data_maxmin/g5_a_aa_maxmin.dat",action="write",form="unformatted")
   write(20)g5_a_aa_maxmin
   close(20)
 endif
 if( io_nnp_a == 1 .and. io_nnp_b == 1 )then
   open(unit=20,file="./data_maxmin/g5_a_ab_maxmin.dat",action="write",form="unformatted")
   write(20)g5_a_ab_maxmin
   close(20)
 endif
 if( io_nnp_a == 1 .and. io_nnp_c == 1 )then
   open(unit=20,file="./data_maxmin/g5_a_ac_maxmin.dat",action="write",form="unformatted")
   write(20)g5_a_ac_maxmin
   close(20)
 endif
 if( io_nnp_a == 1 .and. io_nnp_b == 1 )then
   open(unit=20,file="./data_maxmin/g5_a_bb_maxmin.dat",action="write",form="unformatted")
   write(20)g5_a_bb_maxmin
   close(20)
 endif
 if( io_nnp_a == 1 .and. io_nnp_b == 1 .and. io_nnp_c == 1 )then
   open(unit=20,file="./data_maxmin/g5_a_bc_maxmin.dat",action="write",form="unformatted")
   write(20)g5_a_bc_maxmin
   close(20)
 endif
 if( io_nnp_a == 1 .and. io_nnp_c == 1 )then
   open(unit=20,file="./data_maxmin/g5_a_cc_maxmin.dat",action="write",form="unformatted")
   write(20)g5_a_cc_maxmin
   close(20)
 endif
 
 !B
 if( io_nnp_b == 1 .and. io_nnp_a == 1 )then
   open(unit=20,file="./data_maxmin/g2_b_a_maxmin.dat",action="write",form="unformatted")
   write(20)g2_b_a_maxmin
   close(20)
 endif
 if( io_nnp_b == 1 )then
   open(unit=20,file="./data_maxmin/g2_b_b_maxmin.dat",action="write",form="unformatted")
   write(20)g2_b_b_maxmin
   close(20)
 endif
 if( io_nnp_b == 1 .and. io_nnp_c == 1 )then
   open(unit=20,file="./data_maxmin/g2_b_c_maxmin.dat",action="write",form="unformatted")
   write(20)g2_b_c_maxmin
   close(20)
 endif
 if( io_nnp_b == 1 .and. io_nnp_a == 1 )then
   open(unit=20,file="./data_maxmin/g5_b_aa_maxmin.dat",action="write",form="unformatted")
   write(20)g5_b_aa_maxmin
   close(20)
 endif
 if( io_nnp_b == 1 .and. io_nnp_a == 1 )then
   open(unit=20,file="./data_maxmin/g5_b_ab_maxmin.dat",action="write",form="unformatted")
   write(20)g5_b_ab_maxmin
   close(20)
 endif
 if( io_nnp_b == 1 .and. io_nnp_a == 1 .and. io_nnp_c == 1 )then
   open(unit=20,file="./data_maxmin/g5_b_ac_maxmin.dat",action="write",form="unformatted")
   write(20)g5_b_ac_maxmin
   close(20)
 endif
 if( io_nnp_b == 1 )then
   open(unit=20,file="./data_maxmin/g5_b_bb_maxmin.dat",action="write",form="unformatted")
   write(20)g5_b_bb_maxmin
   close(20)
 endif
 if( io_nnp_b == 1 .and. io_nnp_c == 1 )then
   open(unit=20,file="./data_maxmin/g5_b_bc_maxmin.dat",action="write",form="unformatted")
   write(20)g5_b_bc_maxmin
   close(20)
 endif
 if( io_nnp_b == 1 .and. io_nnp_c == 1 )then
   open(unit=20,file="./data_maxmin/g5_b_cc_maxmin.dat",action="write",form="unformatted")
   write(20)g5_b_cc_maxmin
   close(20)
 endif

 !C
 if( io_nnp_c == 1 .and. io_nnp_a == 1 )then
   open(unit=20,file="./data_maxmin/g2_c_a_maxmin.dat",action="write",form="unformatted")
   write(20)g2_c_a_maxmin
   close(20)
 endif
 if( io_nnp_c == 1 .and. io_nnp_b == 1 )then
   open(unit=20,file="./data_maxmin/g2_c_b_maxmin.dat",action="write",form="unformatted")
   write(20)g2_c_b_maxmin
   close(20)
 endif
 if( io_nnp_c == 1 )then
   open(unit=20,file="./data_maxmin/g2_c_c_maxmin.dat",action="write",form="unformatted")
   write(20)g2_c_c_maxmin
   close(20)
 endif
 if( io_nnp_c == 1 .and. io_nnp_a == 1 )then
   open(unit=20,file="./data_maxmin/g5_c_aa_maxmin.dat",action="write",form="unformatted")
   write(20)g5_c_aa_maxmin
   close(20)
 endif
 if( io_nnp_c == 1 .and. io_nnp_a == 1 .and. io_nnp_b == 1 )then
   open(unit=20,file="./data_maxmin/g5_c_ab_maxmin.dat",action="write",form="unformatted")
   write(20)g5_c_ab_maxmin
   close(20)
 endif
 if( io_nnp_c == 1 .and. io_nnp_a == 1 )then
   open(unit=20,file="./data_maxmin/g5_c_ac_maxmin.dat",action="write",form="unformatted")
   write(20)g5_c_ac_maxmin
   close(20)
 endif
 if( io_nnp_c == 1 .and. io_nnp_b == 1 )then
   open(unit=20,file="./data_maxmin/g5_c_bb_maxmin.dat",action="write",form="unformatted")
   write(20)g5_c_bb_maxmin
   close(20)
 endif
 if( io_nnp_c == 1 .and. io_nnp_b == 1 )then
   open(unit=20,file="./data_maxmin/g5_c_bc_maxmin.dat",action="write",form="unformatted")
   write(20)g5_c_bc_maxmin
   close(20)
 endif
 if( io_nnp_c == 1 )then
   open(unit=20,file="./data_maxmin/g5_c_cc_maxmin.dat",action="write",form="unformatted")
   write(20)g5_c_cc_maxmin
   close(20)
 endif
  

 write(*,*) " "
 write(*,*) " Finish !!"
 write(*,*) " "


   deallocate( g2_a_a_maxmin,  g2_a_b_maxmin,  g2_a_c_maxmin, &
               g5_a_aa_maxmin, g5_a_ab_maxmin, g5_a_ac_maxmin, &
               g5_a_bb_maxmin, g5_a_bc_maxmin, &
               g5_a_cc_maxmin )
   deallocate( g2_b_a_maxmin,  g2_b_b_maxmin,  g2_b_c_maxmin, &
               g5_b_aa_maxmin, g5_b_ab_maxmin, g5_b_ac_maxmin, &
               g5_b_bb_maxmin, g5_b_bc_maxmin, &
               g5_b_cc_maxmin )
   deallocate( g2_c_a_maxmin,  g2_c_b_maxmin,  g2_c_c_maxmin, &
               g5_c_aa_maxmin, g5_c_ab_maxmin, g5_c_ac_maxmin, &
               g5_c_bb_maxmin, g5_c_bc_maxmin, &
               g5_c_cc_maxmin )


end
