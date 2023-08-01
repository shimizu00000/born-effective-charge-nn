subroutine read_input_nnp4( &
           elem_a, &
           num_g2_a_a,  num_g2_a_b,  num_g2_a_c, num_g2_a_d, &
           num_g5_a_aa, num_g5_a_ab, num_g5_a_ac, num_g5_a_ad, &
           num_g5_a_bb, num_g5_a_bc, num_g5_a_bd, &
           num_g5_a_cc, num_g5_a_cd, &
           num_g5_a_dd, &
           elem_b, &
           num_g2_b_a,  num_g2_b_b,  num_g2_b_c, num_g2_b_d, &
           num_g5_b_aa, num_g5_b_ab, num_g5_b_ac, num_g5_b_ad, &
           num_g5_b_bb, num_g5_b_bc, num_g5_b_bd, &
           num_g5_b_cc, num_g5_b_cd, &
           num_g5_b_dd, &
           elem_c, &
           num_g2_c_a,  num_g2_c_b,  num_g2_c_c, num_g2_c_d, &
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

implicit none

character(len=120) tag
character(len=120) cread
integer            iread
double precision   dread

integer i, stat, counter

integer node_in_a, node_in_b, node_in_c, node_in_d

!A-*
integer num_g2_a_a,  num_g2_a_b,  num_g2_a_c,  num_g2_a_d, &
        num_g5_a_aa, num_g5_a_ab, num_g5_a_ac, num_g5_a_ad, &
        num_g5_a_bb, num_g5_a_bc, num_g5_a_bd, &
        num_g5_a_cc, num_g5_a_cd, &
        num_g5_a_dd
!B-*
integer num_g2_b_a,  num_g2_b_b,  num_g2_b_c,  num_g2_b_d, &
        num_g5_b_aa, num_g5_b_ab, num_g5_b_ac, num_g5_b_ad, &
        num_g5_b_bb, num_g5_b_bc, num_g5_b_bd, &
        num_g5_b_cc, num_g5_b_cd, &
        num_g5_b_dd
!C-*
integer num_g2_c_a,  num_g2_c_b,  num_g2_c_c,  num_g2_c_d, &
        num_g5_c_aa, num_g5_c_ab, num_g5_c_ac, num_g5_c_ad, &
        num_g5_c_bb, num_g5_c_bc, num_g5_c_bd, &
        num_g5_c_cc, num_g5_c_cd, &
        num_g5_c_dd
!D-*
integer num_g2_d_a,  num_g2_d_b,  num_g2_d_c,  num_g2_d_d, &
        num_g5_d_aa, num_g5_d_ab, num_g5_d_ac, num_g5_d_ad, &
        num_g5_d_bb, num_g5_d_bc, num_g5_d_bd, &
        num_g5_d_cc, num_g5_d_cd, &
        num_g5_d_dd

integer nelem, natom
integer natom_a, natom_b, natom_c, natom_d
character(len=4) elem_a, elem_b, elem_c, elem_d

integer md_steps
integer neighbor
double precision rc, rn
character(len=10) phi
character(len=120) dir_data, dir_sf
character dummyc


!***********************
!Read input_nnp.dat file
!***********************
!input_nnp.dat includes
!----------------------
!rc = 7.0d0
!rn = 1.0d0
!neighbor = 300
!phi = tanh
!
!elemen_a = Au
!node_in_a = 44
!num_g2_a_a = 8
!num_g2_a_b = 8
!num_g2_a_c = 8
!num_g5_a_aa = 22,
!num_g5_a_ab = 22
!num_g5_a_ac = 22
!num_g5_a_bb = 22
!num_g5_a_bc = 22,
!num_g5_a_cc = 22

!also b,c,d
!---
 open(unit=98,file="input_nnp.dat",action="read")
!Number of lines in input_nnp.dat
 counter = 0
 rewind(98)
 do
   read( 98, '(A1)', iostat=stat ) dummyc
   if( stat /= 0 ) exit
   if( trim( dummyc ) /= '' ) counter = counter + 1
 enddo

 rc = 0.0d0
 rn = 0.0d0
 neighbor = 0

!A
 elem_a = "none"
 node_in_a   = 0 
 num_g2_a_a  = 0 ; num_g2_a_b  = 0 ; num_g2_a_c  = 0 ; num_g2_a_d  = 0
 num_g5_a_aa = 0 ; num_g5_a_ab = 0 ; num_g5_a_ac = 0 ; num_g5_a_ad = 0
 num_g5_a_bb = 0 ; num_g5_a_bc = 0 ; num_g5_a_bd = 0
 num_g5_a_cc = 0 ; num_g5_a_cd = 0
 num_g5_a_dd = 0
!B
 elem_b = "none"
 node_in_b   = 0 
 num_g2_b_a  = 0 ; num_g2_b_b  = 0 ; num_g2_b_c  = 0 ; num_g2_b_d  = 0
 num_g5_b_aa = 0 ; num_g5_b_ab = 0 ; num_g5_b_ac = 0 ; num_g5_b_ad = 0
 num_g5_b_bb = 0 ; num_g5_b_bc = 0 ; num_g5_b_bd = 0
 num_g5_b_cc = 0 ; num_g5_b_cd = 0 
 num_g5_b_dd = 0
!C
 elem_c = "none"
 node_in_c   = 0 
 num_g2_c_a  = 0 ; num_g2_c_b  = 0 ; num_g2_c_c  = 0 ; num_g2_c_d  = 0
 num_g5_c_aa = 0 ; num_g5_c_ab = 0 ; num_g5_c_ac = 0 ; num_g5_c_ad = 0
 num_g5_c_bb = 0 ; num_g5_c_bc = 0 ; num_g5_c_bd = 0 
 num_g5_c_cc = 0 ; num_g5_c_cd = 0 
 num_g5_c_dd = 0
!D
 elem_d = "none"
 node_in_d   = 0
 num_g2_d_a  = 0 ; num_g2_d_b  = 0 ; num_g2_d_c  = 0 ; num_g2_d_d  = 0
 num_g5_d_aa = 0 ; num_g5_d_ab = 0 ; num_g5_d_ac = 0 ; num_g5_d_ad = 0
 num_g5_d_bb = 0 ; num_g5_d_bc = 0 ; num_g5_d_bd = 0
 num_g5_d_cc = 0 ; num_g5_d_cd = 0
 num_g5_d_dd = 0
 
 rewind(98)
 do i = 1, counter

   read(98,*) tag, dummyc, cread

   !parameters for A
   if( tag == "element_a" )then
     elem_a = cread
   endif
   if( tag == "node_in_a" )then
     read(cread,*) iread
     node_in_a   = iread
   endif
   if( tag == "num_g2_a_a" )then
     read(cread,*) iread
     num_g2_a_a  = iread
   endif
   if( tag == "num_g2_a_b" )then
     read(cread,*) iread
     num_g2_a_b  = iread
   endif
   if( tag == "num_g2_a_c" )then
     read(cread,*) iread
     num_g2_a_c  = iread
   endif
   if( tag == "num_g2_a_d" )then
     read(cread,*) iread
     num_g2_a_d  = iread
   endif
   if( tag == "num_g5_a_aa" )then
     read(cread,*) iread
     num_g5_a_aa = iread
   endif
   if( tag == "num_g5_a_ab" )then
     read(cread,*) iread
     num_g5_a_ab = iread
   endif
   if( tag == "num_g5_a_ac" )then
     read(cread,*) iread
     num_g5_a_ac = iread
   endif
   if( tag == "num_g5_a_ad" )then
     read(cread,*) iread
     num_g5_a_ad = iread
   endif
   if( tag == "num_g5_a_bb" )then
     read(cread,*) iread
     num_g5_a_bb = iread
   endif
   if( tag == "num_g5_a_bc" )then
     read(cread,*) iread
     num_g5_a_bc = iread
   endif
   if( tag == "num_g5_a_bd" )then
     read(cread,*) iread
     num_g5_a_bd = iread
   endif
   if( tag == "num_g5_a_cc" )then
     read(cread,*) iread
     num_g5_a_cc = iread
   endif
   if( tag == "num_g5_a_cd" )then
     read(cread,*) iread
     num_g5_a_cd = iread
   endif
   if( tag == "num_g5_a_dd" )then
     read(cread,*) iread
     num_g5_a_dd = iread
   endif

   !parameters for B
   if( tag == "element_b" )then
     elem_b = cread
   endif
   if( tag == "node_in_b" )then
     read(cread,*) iread
     node_in_b   = iread
   endif
   if( tag == "num_g2_b_a" )then
     read(cread,*) iread
     num_g2_b_a  = iread
   endif
   if( tag == "num_g2_b_b" )then
     read(cread,*) iread
     num_g2_b_b  = iread
   endif
   if( tag == "num_g2_b_c" )then
     read(cread,*) iread
     num_g2_b_c  = iread
   endif
   if( tag == "num_g2_b_d" )then
     read(cread,*) iread
     num_g2_b_d  = iread
   endif
   if( tag == "num_g5_b_aa" )then
     read(cread,*) iread
     num_g5_b_aa = iread
   endif
   if( tag == "num_g5_b_ab" )then
     read(cread,*) iread
     num_g5_b_ab = iread
   endif
   if( tag == "num_g5_b_ac" )then
     read(cread,*) iread
     num_g5_b_ac = iread
   endif
   if( tag == "num_g5_b_ad" )then
     read(cread,*) iread
     num_g5_b_ad = iread
   endif
   if( tag == "num_g5_b_bb" )then
     read(cread,*) iread
     num_g5_b_bb = iread
   endif
   if( tag == "num_g5_b_bc" )then
     read(cread,*) iread
     num_g5_b_bc = iread
   endif
   if( tag == "num_g5_b_bd" )then
     read(cread,*) iread
     num_g5_b_bd = iread
   endif
   if( tag == "num_g5_b_cc" )then
     read(cread,*) iread
     num_g5_b_cc = iread
   endif
   if( tag == "num_g5_b_cd" )then
     read(cread,*) iread
     num_g5_b_cd = iread
   endif
   if( tag == "num_g5_b_dd" )then
     read(cread,*) iread
     num_g5_b_dd = iread
   endif

   !parameters for C
   if( tag == "element_c" )then
     elem_c = cread
   endif
   if( tag == "node_in_c" )then
     read(cread,*) iread
     node_in_c   = iread
   endif
   if( tag == "num_g2_c_a" )then
     read(cread,*) iread
     num_g2_c_a  = iread
   endif
   if( tag == "num_g2_c_b" )then
     read(cread,*) iread
     num_g2_c_b  = iread
   endif
   if( tag == "num_g2_c_c" )then
     read(cread,*) iread
     num_g2_c_c  = iread
   endif
   if( tag == "num_g2_c_d" )then
     read(cread,*) iread
     num_g2_c_d  = iread
   endif
   if( tag == "num_g5_c_aa" )then
     read(cread,*) iread
     num_g5_c_aa = iread
   endif
   if( tag == "num_g5_c_ab" )then
     read(cread,*) iread
     num_g5_c_ab = iread
   endif
   if( tag == "num_g5_c_ac" )then
     read(cread,*) iread
     num_g5_c_ac = iread
   endif
   if( tag == "num_g5_c_ad" )then
     read(cread,*) iread
     num_g5_c_ad = iread
   endif
   if( tag == "num_g5_c_bb" )then
     read(cread,*) iread
     num_g5_c_bb = iread
   endif
   if( tag == "num_g5_c_bc" )then
     read(cread,*) iread
     num_g5_c_bc = iread
   endif
   if( tag == "num_g5_c_bd" )then
     read(cread,*) iread
     num_g5_c_bd = iread
   endif
   if( tag == "num_g5_c_cc" )then
     read(cread,*) iread
     num_g5_c_cc = iread
   endif
   if( tag == "num_g5_c_cd" )then
     read(cread,*) iread
     num_g5_c_cd = iread
   endif
   if( tag == "num_g5_c_dd" )then
     read(cread,*) iread
     num_g5_c_dd = iread
   endif

   !parameters for D
   if( tag == "element_d" )then
     elem_d = cread
   endif
   if( tag == "node_in_d" )then
     read(cread,*) iread
     node_in_d   = iread
   endif
   if( tag == "num_g2_d_a" )then
     read(cread,*) iread
     num_g2_d_a  = iread
   endif
   if( tag == "num_g2_d_b" )then
     read(cread,*) iread
     num_g2_d_b  = iread
   endif
   if( tag == "num_g2_d_c" )then
     read(cread,*) iread
     num_g2_d_c  = iread
   endif
   if( tag == "num_g2_d_d" )then
     read(cread,*) iread
     num_g2_d_d  = iread
   endif
   if( tag == "num_g5_d_aa" )then
     read(cread,*) iread
     num_g5_d_aa = iread
   endif
   if( tag == "num_g5_d_ab" )then
     read(cread,*) iread
     num_g5_d_ab = iread
   endif
   if( tag == "num_g5_d_ac" )then
     read(cread,*) iread
     num_g5_d_ac = iread
   endif
   if( tag == "num_g5_d_ad" )then
     read(cread,*) iread
     num_g5_d_ad = iread
   endif
   if( tag == "num_g5_d_bb" )then
     read(cread,*) iread
     num_g5_d_bb = iread
   endif
   if( tag == "num_g5_d_bc" )then
     read(cread,*) iread
     num_g5_d_bc = iread
   endif
   if( tag == "num_g5_d_bd" )then
     read(cread,*) iread
     num_g5_d_bd = iread
   endif
   if( tag == "num_g5_d_cc" )then
     read(cread,*) iread
     num_g5_d_cc = iread
   endif
   if( tag == "num_g5_d_cd" )then
     read(cread,*) iread
     num_g5_d_cd = iread
   endif
   if( tag == "num_g5_d_dd" )then
     read(cread,*) iread
     num_g5_d_dd = iread
   endif

   !other parameters
   if( tag == "rc" )then
     read(cread,*) dread
     rc = dread
   endif
   if( tag == "rn" )then
     read(cread,*) dread
     rn = dread
   endif
   if( tag == "neighbor" )then
     read(cread,*) iread
     neighbor = iread
   endif
   if( tag == "phi" )then
     phi = cread
   endif
   if( tag == "dir_data" )then
     dir_data = cread
   endif
   if( tag == "dir_sf" )then
     dir_sf = cread
   endif


 enddo ! i = 1, counter


 close(98) ! input_nnp.dat


 write(*,*) " "
 write(*,*) " Parameters for symmetry function"
 write(*,*) " ============================================================"

 !1 element
 !A
 if( elem_a /= "none" .and. elem_b == "none" &
     .and. elem_c == "none" .and. elem_d == "none" )then

   write(*,*) " atom_a: ", elem_a
   write(*,*) " input_node: ", node_in_a
   write(*,*) " num_g2"
   write(*,*) " a-a: ", num_g2_a_a
   write(*,*) " num_g5"
   write(*,*) " a-aa: ", num_g5_a_aa
   write(*,*) " "

 !B
 elseif( elem_a == "none" .and. elem_b /= "none" &
         .and. elem_c == "none" .and. elem_d == "none" )then

   write(*,*) " atom_b: ", elem_b
   write(*,*) " input_node: ", node_in_b
   write(*,*) " num_g2"
   write(*,*) " b-b: ", num_g2_b_b
   write(*,*) " num_g5"
   write(*,*) " b-bb: ", num_g5_b_bb
   write(*,*) " "

 !C
 elseif( elem_a == "none" .and. elem_b == "none" &
         .and. elem_c /= "none" .and. elem_d == "none" )then

   write(*,*) " atom_c: ", elem_c
   write(*,*) " input_node: ", node_in_c
   write(*,*) " num_g2"
   write(*,*) " c-c: ", num_g2_c_c
   write(*,*) " num_g5"
   write(*,*) " c-cc: ", num_g5_c_cc
   write(*,*) " "

 !D
 elseif( elem_a == "none" .and. elem_b == "none" &
         .and. elem_c == "none" .and. elem_d /= "none" )then

   write(*,*) " atom_d: ", elem_d
   write(*,*) " input_node: ", node_in_d
   write(*,*) " num_g2"
   write(*,*) " d-d: ", num_g2_d_d
   write(*,*) " num_g5"
   write(*,*) " d-dd: ", num_g5_d_dd
   write(*,*) " "

 endif


 !2 elements
 !A,B
 if( elem_a /= "none" .and. elem_b /= "none" &
     .and. elem_c == "none" .and. elem_d == "none" )then

   write(*,*) " atom_a: ", elem_a
   write(*,*) " input_node: ", node_in_a
   write(*,*) " num_g2"
   write(*,*) " a-a: ", num_g2_a_a, "a-b: ", num_g2_a_b
   write(*,*) " num_g5"
   write(*,*) " a-aa: ", num_g5_a_aa, "a-ab: ", num_g5_a_ab, "a-bb: ", num_g5_a_bb
   write(*,*) " "
   write(*,*) " atom_b: ", elem_b
   write(*,*) " input_node: ", node_in_b
   write(*,*) " num_g2"
   write(*,*) " b-a: ", num_g2_b_a, "b-b: ", num_g2_b_b
   write(*,*) " num_g5"
   write(*,*) " b-aa: ", num_g5_b_aa, "b-ab: ", num_g5_b_ab, "b-bb: ", num_g5_b_bb
   write(*,*) " "

 !A,C
 elseif( elem_a /= "none" .and. elem_b == "none" &
         .and. elem_c /= "none" .and. elem_d == "none" )then

   write(*,*) " atom_a: ", elem_a
   write(*,*) " input_node: ", node_in_a
   write(*,*) " num_g2"
   write(*,*) " a-a: ", num_g2_a_a, "a-c: ", num_g2_a_c
   write(*,*) " num_g5"
   write(*,*) " a-aa: ", num_g5_a_aa, "a-ac: ", num_g5_a_ac, "a-cc: ", num_g5_a_cc
   write(*,*) " "
   write(*,*) " atom_c: ", elem_c
   write(*,*) " input_node: ", node_in_c
   write(*,*) " num_g2"
   write(*,*) " c-a: ", num_g2_c_a, "c-c: ", num_g2_c_c
   write(*,*) " num_g5"
   write(*,*) " c-aa: ", num_g5_c_aa, "c-ac: ", num_g5_c_ac, "c-cc: ", num_g5_c_cc
   write(*,*) " "

 !A,D
 elseif( elem_a /= "none" .and. elem_b == "none" &
         .and. elem_c == "none" .and. elem_d /= "none" )then

   write(*,*) " atom_a: ", elem_a
   write(*,*) " input_node: ", node_in_a
   write(*,*) " num_g2"
   write(*,*) " a-a: ", num_g2_a_a, "a-d: ", num_g2_a_d
   write(*,*) " num_g5"
   write(*,*) " a-aa: ", num_g5_a_aa, "a-ad: ", num_g5_a_ad, "a-dd: ", num_g5_a_dd
   write(*,*) " "
   write(*,*) " atom_d: ", elem_d
   write(*,*) " input_node: ", node_in_d
   write(*,*) " num_g2"
   write(*,*) " d-a: ", num_g2_d_a, "d-d: ", num_g2_d_d
   write(*,*) " num_g5"
   write(*,*) " d-aa: ", num_g5_d_aa, "d-ad: ", num_g5_d_ad, "d-dd: ", num_g5_d_dd
   write(*,*) " "

 !B,C
 elseif( elem_a == "none" .and. elem_b /= "none" &
         .and. elem_c /= "none" .and. elem_d == "none" )then

   write(*,*) " atom_b: ", elem_b
   write(*,*) " input_node: ", node_in_b
   write(*,*) " num_g2"
   write(*,*) " b-b: ", num_g2_b_b, "b-c: ", num_g2_b_c
   write(*,*) " num_g5"
   write(*,*) " b-bb: ", num_g5_b_bb, "b-bc: ", num_g5_b_bc, "b-cc: ", num_g5_b_cc
   write(*,*) " "
   write(*,*) " atom_c: ", elem_c
   write(*,*) " input_node: ", node_in_c
   write(*,*) " num_g2"
   write(*,*) " c-b: ", num_g2_c_b, "c-c: ", num_g2_c_c
   write(*,*) " num_g5"
   write(*,*) " c-bb: ", num_g5_c_bb, "c-bc: ", num_g5_c_bc, "c-cc: ", num_g5_c_cc
   write(*,*) " "

 !B,D
 elseif( elem_a == "none" .and. elem_b /= "none" &
         .and. elem_c == "none" .and. elem_d /= "none" )then

   write(*,*) " atom_b: ", elem_b
   write(*,*) " input_node: ", node_in_b
   write(*,*) " num_g2"
   write(*,*) " b-b: ", num_g2_b_b, "b-d: ", num_g2_b_d
   write(*,*) " num_g5"
   write(*,*) " b-bb: ", num_g5_b_bb, "b-bd: ", num_g5_b_bd, "b-dd: ", num_g5_b_dd
   write(*,*) " "
   write(*,*) " atom_d: ", elem_d
   write(*,*) " input_node: ", node_in_d
   write(*,*) " num_g2"
   write(*,*) " d-b: ", num_g2_d_b, "d-d: ", num_g2_d_d
   write(*,*) " num_g5"
   write(*,*) " d-bb: ", num_g5_d_bb, "d-bd: ", num_g5_d_bd, "d-dd: ", num_g5_d_dd
   write(*,*) " "

 !C,D
 elseif( elem_a == "none" .and. elem_b == "none" &
         .and. elem_c /= "none" .and. elem_d /= "none" )then

   write(*,*) " atom_c: ", elem_c
   write(*,*) " input_node: ", node_in_c
   write(*,*) " num_g2"
   write(*,*) " c-c: ", num_g2_c_c, "c-d: ", num_g2_c_d
   write(*,*) " num_g5"
   write(*,*) " c-cc: ", num_g5_c_cc, "c-cd: ", num_g5_c_cd, "c-dd: ", num_g5_c_dd
   write(*,*) " "
   write(*,*) " atom_d: ", elem_d
   write(*,*) " input_node: ", node_in_d
   write(*,*) " num_g2"
   write(*,*) " d-c: ", num_g2_d_c, "d-d: ", num_g2_d_d
   write(*,*) " num_g5"
   write(*,*) " d-cc: ", num_g5_d_cc, "d-cd: ", num_g5_d_cd, "d-dd: ", num_g5_d_dd
   write(*,*) " "

 endif


 !3 elements
 !A,B,C
 if( elem_a /= "none" .and. elem_b /= "none" &
     .and. elem_c /= "none" .and. elem_d == "none" )then

   write(*,*) " atom_a: ", elem_a
   write(*,*) " input_node: ", node_in_a
   write(*,*) " num_g2"
   write(*,*) " a-a: ", num_g2_a_a, "a-b: ", num_g2_a_b, "a-c: ", num_g2_a_c
   write(*,*) " num_g5"
   write(*,*) " a-aa: ", num_g5_a_aa, "a-ab: ", num_g5_a_ab, "a-ac: ", num_g5_a_ac, &
              "a-bb: ", num_g5_a_bb, "a-bc: ", num_g5_a_bc, "a-cc: ", num_g5_a_cc
   write(*,*) " "
   write(*,*) " atom_b: ", elem_b
   write(*,*) " input_node: ", node_in_b
   write(*,*) " num_g2"
   write(*,*) " b-a: ", num_g2_b_a, "b-b: ", num_g2_b_b, "b-c: ", num_g2_b_c
   write(*,*) " num_g5"
   write(*,*) " b-aa: ", num_g5_b_aa, "b-ab: ", num_g5_b_ab, "b-ac: ", num_g5_b_ac, &
              "b-bb: ", num_g5_b_bb, "b-bc: ", num_g5_b_bc, "b-cc: ", num_g5_b_cc
   write(*,*) " "
   write(*,*) " atom_c: ", elem_c
   write(*,*) " input_node: ", node_in_c
   write(*,*) " num_g2"
   write(*,*) " c-a: ", num_g2_c_a, "c-b: ", num_g2_c_b, "c-c: ", num_g2_c_c
   write(*,*) " num_g5"
   write(*,*) " c-aa: ", num_g5_c_aa, "c-ab: ", num_g5_c_ab, "c-ac: ", num_g5_c_ac, &
              "c-bb: ", num_g5_c_bb, "c-bc: ", num_g5_c_bc, "c-cc: ", num_g5_c_cc
   write(*,*) " "

 !A,B,D
 elseif( elem_a /= "none" .and. elem_b /= "none" &
         .and. elem_c == "none" .and. elem_d /= "none" )then

   write(*,*) " atom_a: ", elem_a
   write(*,*) " input_node: ", node_in_a
   write(*,*) " num_g2"
   write(*,*) " a-a: ", num_g2_a_a, "a-b: ", num_g2_a_b, "a-d: ", num_g2_a_d
   write(*,*) " num_g5"
   write(*,*) " a-aa: ", num_g5_a_aa, "a-ab: ", num_g5_a_ab, "a-ad: ", num_g5_a_ad, &
              "a-bb: ", num_g5_a_bb, "a-bd: ", num_g5_a_bd, "a-dd: ", num_g5_a_dd
   write(*,*) " "
   write(*,*) " atom_b: ", elem_b
   write(*,*) " input_node: ", node_in_b
   write(*,*) " num_g2"
   write(*,*) " b-a: ", num_g2_b_a, "b-b: ", num_g2_b_b, "b-d: ", num_g2_b_d
   write(*,*) " num_g5"
   write(*,*) " b-aa: ", num_g5_b_aa, "b-ab: ", num_g5_b_ab, "b-ad: ", num_g5_b_ad, &
              "b-bb: ", num_g5_b_bb, "b-bd: ", num_g5_b_bd, "b-dd: ", num_g5_b_dd
   write(*,*) " "
   write(*,*) " atom_d: ", elem_d
   write(*,*) " input_node: ", node_in_d
   write(*,*) " num_g2"
   write(*,*) " d-a: ", num_g2_d_a, "d-b: ", num_g2_d_b, "d-d: ", num_g2_d_d
   write(*,*) " num_g5"
   write(*,*) " d-aa: ", num_g5_d_aa, "d-ab: ", num_g5_d_ab, "d-ad: ", num_g5_d_ad, &
              "d-bb: ", num_g5_d_bb, "d-bd: ", num_g5_d_bd, "d-dd: ", num_g5_d_dd
   write(*,*) " "

 !A,C,D
 elseif( elem_a /= "none" .and. elem_b == "none" &
         .and. elem_c /= "none" .and. elem_d /= "none" )then

   write(*,*) " atom_a: ", elem_a
   write(*,*) " input_node: ", node_in_a
   write(*,*) " num_g2"
   write(*,*) " a-a: ", num_g2_a_a, "a-c: ", num_g2_a_c, "a-d: ", num_g2_a_d
   write(*,*) " num_g5"
   write(*,*) " a-aa: ", num_g5_a_aa, "a-ac: ", num_g5_a_ac, "a-ad: ", num_g5_a_ad, &
              "a-cc: ", num_g5_a_cc, "a-cd: ", num_g5_a_cd, "a-dd: ", num_g5_a_dd
   write(*,*) " "
   write(*,*) " atom_c: ", elem_c
   write(*,*) " input_node: ", node_in_c
   write(*,*) " num_g2"
   write(*,*) " c-a: ", num_g2_c_a, "c-c: ", num_g2_c_c, "c-d: ", num_g2_c_d
   write(*,*) " num_g5"
   write(*,*) " c-aa: ", num_g5_c_aa, "c-ac: ", num_g5_c_ab, "c-ad: ", num_g5_c_ac, &
              "c-cc: ", num_g5_c_bb, "c-cd: ", num_g5_c_bc, "c-dd: ", num_g5_c_cc
   write(*,*) " "
   write(*,*) " atom_d: ", elem_d
   write(*,*) " input_node: ", node_in_d
   write(*,*) " num_g2"
   write(*,*) " d-a: ", num_g2_d_a, "d-c: ", num_g2_d_c, "d-d: ", num_g2_d_d
   write(*,*) " num_g5"
   write(*,*) " d-aa: ", num_g5_d_aa, "d-ac: ", num_g5_d_ac, "d-ad: ", num_g5_d_ad, &
              "d-cc: ", num_g5_d_cc, "d-cd: ", num_g5_d_cd, "d-dd: ", num_g5_d_dd
   write(*,*) " "

 !B,C,D
 elseif( elem_a == "none" .and. elem_b /= "none" &
         .and. elem_c /= "none" .and. elem_d /= "none" )then

   write(*,*) " atom_b: ", elem_b
   write(*,*) " input_node: ", node_in_b
   write(*,*) " num_g2"
   write(*,*) " b-b: ", num_g2_b_b, "b-c: ", num_g2_b_c, "b-d: ", num_g2_b_d
   write(*,*) " num_g5"
   write(*,*) " b-bb: ", num_g5_b_bb, "b-bc: ", num_g5_b_bc, "b-bd: ", num_g5_b_bd, &
              "b-cc: ", num_g5_b_cc, "b-cd: ", num_g5_b_cd, "b-dd: ", num_g5_b_dd
   write(*,*) " "
   write(*,*) " atom_c: ", elem_c
   write(*,*) " input_node: ", node_in_c
   write(*,*) " num_g2"
   write(*,*) " c-b: ", num_g2_c_b, "c-c: ", num_g2_c_c, "c-d: ", num_g2_c_d
   write(*,*) " num_g5"
   write(*,*) " c-bb: ", num_g5_c_bb, "c-bc: ", num_g5_c_bc, "c-bd: ", num_g5_c_bd, &
              "c-cc: ", num_g5_c_cc, "c-cd: ", num_g5_c_cd, "c-dd: ", num_g5_c_dd
   write(*,*) " "
   write(*,*) " atom_d: ", elem_d
   write(*,*) " input_node: ", node_in_d
   write(*,*) " num_g2"
   write(*,*) " d-b: ", num_g2_d_b, "d-c: ", num_g2_d_c, "d-d: ", num_g2_d_d
   write(*,*) " num_g5"
   write(*,*) " d-bb: ", num_g5_d_bb, "d-bc: ", num_g5_d_bc, "d-bd: ", num_g5_d_bd, &
              "d-cc: ", num_g5_d_cc, "d-cd: ", num_g5_d_cd, "d-dd: ", num_g5_d_dd
   write(*,*) " "

 endif


 !4 elements
 !A,B,C,D
 if( elem_a /= "none" .and. elem_b /= "none" &
     .and. elem_c /= "none" .and. elem_d /= "none" )then

   write(*,*) " atom_a: ", elem_a
   write(*,*) " input_node: ", node_in_a
   write(*,*) " num_g2"
   write(*,*) " a-a: ", num_g2_a_a, "a-b: ", num_g2_a_b, &
              "a-c: ", num_g2_a_c, "a-d: ", num_g2_a_d
   write(*,*) " num_g5"
   write(*,*) " a-aa: ", num_g5_a_aa, "a-ab: ", num_g5_a_ab, "a-ac: ", num_g5_a_ac, &
              "a-ad: ", num_g5_a_ad
   write(*,*) " a-bb: ", num_g5_a_bb, "a-bc: ", num_g5_a_bc, "a-bd: ", num_g5_a_bd
   write(*,*) " a-cc: ", num_g5_a_cc, "a-cd: ", num_g5_a_cd
   write(*,*) " a-dd: ", num_g5_a_dd
   write(*,*) " "
   write(*,*) " atom_b: ", elem_b
   write(*,*) " input_node: ", node_in_b
   write(*,*) " num_g2"
   write(*,*) " b-a: ", num_g2_b_a, "b-b: ", num_g2_b_b, &
              "b-c: ", num_g2_b_c, "b-d: ", num_g2_b_d
   write(*,*) " num_g5"
   write(*,*) " b-aa: ", num_g5_b_aa, "b-ab: ", num_g5_b_ab, "b-ac: ", num_g5_b_ac, &
              "b-ad: ", num_g5_b_ad
   write(*,*) " b-bb: ", num_g5_b_bb, "b-bc: ", num_g5_b_bc, "b-bd: ", num_g5_b_bd
   write(*,*) " b-cc: ", num_g5_b_cc, "b-cd: ", num_g5_b_cd
   write(*,*) " b-dd: ", num_g5_b_dd
   write(*,*) " "
   write(*,*) " atom_c: ", elem_c
   write(*,*) " input_node: ", node_in_c
   write(*,*) " num_g2"
   write(*,*) " c-a: ", num_g2_c_a, "c-b: ", num_g2_c_b, &
              "c-c: ", num_g2_c_c, "c-d: ", num_g2_c_d
   write(*,*) " num_g5"
   write(*,*) " c-aa: ", num_g5_c_aa, "c-ab: ", num_g5_c_ab, "c-ac: ", num_g5_c_ac, &
              "c-ad: ", num_g5_c_ad
   write(*,*) " c-bb: ", num_g5_c_bb, "c-bc: ", num_g5_c_bc, "c-bd: ", num_g5_c_bd
   write(*,*) " c-cc: ", num_g5_c_cc, "c-cd: ", num_g5_c_cd
   write(*,*) " c-dd: ", num_g5_c_dd
   write(*,*) " "
   write(*,*) " atom_d: ", elem_d
   write(*,*) " input_node: ", node_in_d
   write(*,*) " num_g2"
   write(*,*) " d-a: ", num_g2_d_a, "d-b: ", num_g2_d_b, &
              "d-c: ", num_g2_d_c, "d-d: ", num_g2_d_d
   write(*,*) " num_g5"
   write(*,*) " d-aa: ", num_g5_d_aa, "d-ab: ", num_g5_d_ab, "d-ac: ", num_g5_d_ac, &
              "d-ad: ", num_g5_d_ad
   write(*,*) " d-bb: ", num_g5_d_bb, "d-bc: ", num_g5_d_bc, "d-bd: ", num_g5_d_bd
   write(*,*) " d-cc: ", num_g5_d_cc, "d-cd: ", num_g5_d_cd
   write(*,*) " d-dd: ", num_g5_d_dd
   write(*,*) " "

 endif


 write(*,*) " "
 write(*,*) " Activation function: ", phi
 write(*,*) " Cutoff distance:", rc, " [Angst]"
 write(*,*) " Margin distance:", rn, " [Angst]"
 write(*,*) " Neighbor list size:", neighbor, " atoms"
 write(*,*) " "
 write(*,*) " Data: ", trim(adjustl(dir_data))
 write(*,*) " SF: ", trim(adjustl(dir_sf))
 write(*,*) " "
 write(*,*) " ============================================================"


 112 format( A9, A9 )
 113 format( A11, I5 )
 114 format( A9, A7, I5, A7, I5, A7, I5 )
 115 format( A9, A7, I5, A7, I5, A7, I5, A7, I5, A7, I5, A7, I5 )
 116 format( A17, F5.2, A9, A17, F5.2, A8 )
 117 format( A20, I6, A6 )


end subroutine read_input_nnp4
