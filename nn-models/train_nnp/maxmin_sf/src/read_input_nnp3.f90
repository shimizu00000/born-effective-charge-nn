subroutine read_input_nnp3( &
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

implicit none

character(len=100) tag
character(len=100) cread
integer            iread
double precision   dread

integer i, stat, counter

integer node_in_a, node_in_b, node_in_c

!A-*
integer num_g2_a_a,  num_g2_a_b,  num_g2_a_c, &
        num_g5_a_aa, num_g5_a_ab, num_g5_a_ac, &
        num_g5_a_bb, num_g5_a_bc, &
        num_g5_a_cc
!B-*
integer num_g2_b_a,  num_g2_b_b,  num_g2_b_c, &
        num_g5_b_aa, num_g5_b_ab, num_g5_b_ac, &
        num_g5_b_bb, num_g5_b_bc, &
        num_g5_b_cc
!C-*
integer num_g2_c_a,  num_g2_c_b,  num_g2_c_c, &
        num_g5_c_aa, num_g5_c_ab, num_g5_c_ac, &
        num_g5_c_bb, num_g5_c_bc, &
        num_g5_c_cc

integer   nelem, natom

integer   natom_a, natom_b, natom_c
integer   natom_1, natom_2, natom_3

character(len=4) elem_a, elem_b, elem_c
character(len=4) elem_1, elem_2, elem_3

integer md_steps

integer neighbor

double precision rc, rn

character(len=10) phi

character(len=120) dir_data, dir_sf

integer io_mix

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
 num_g2_a_a  = 0 ; num_g2_a_b  = 0 ; num_g2_a_c  = 0
 num_g5_a_aa = 0 ; num_g5_a_ab = 0 ; num_g5_a_ac = 0
 num_g5_a_bb = 0 ; num_g5_a_bc = 0
 num_g5_a_cc = 0
!B
 elem_b = "none"
 node_in_b   = 0 
 num_g2_b_a  = 0 ; num_g2_b_b  = 0 ; num_g2_b_c  = 0
 num_g5_b_aa = 0 ; num_g5_b_ab = 0 ; num_g5_b_ac = 0
 num_g5_b_bb = 0 ; num_g5_b_bc = 0
 num_g5_b_cc = 0
!C
 elem_c = "none"
 node_in_c   = 0 
 num_g2_c_a  = 0 ; num_g2_c_b  = 0 ; num_g2_c_c  = 0
 num_g5_c_aa = 0 ; num_g5_c_ab = 0 ; num_g5_c_ac = 0
 num_g5_c_bb = 0 ; num_g5_c_bc = 0
 num_g5_c_cc = 0

 io_mix = 0
 

print*,"COUNTER",counter

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
   if( tag == "num_g5_a_bb" )then
     read(cread,*) iread
     num_g5_a_bb = iread
   endif
   if( tag == "num_g5_a_bc" )then
     read(cread,*) iread
     num_g5_a_bc = iread
   endif
   if( tag == "num_g5_a_cc" )then
     read(cread,*) iread
     num_g5_a_cc = iread
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
   if( tag == "num_g5_b_bb" )then
     read(cread,*) iread
     num_g5_b_bb = iread
   endif
   if( tag == "num_g5_b_bc" )then
     read(cread,*) iread
     num_g5_b_bc = iread
   endif
   if( tag == "num_g5_b_cc" )then
     read(cread,*) iread
     num_g5_b_cc = iread
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
   if( tag == "num_g5_c_bb" )then
     read(cread,*) iread
     num_g5_c_bb = iread
   endif
   if( tag == "num_g5_c_bc" )then
     read(cread,*) iread
     num_g5_c_bc = iread
   endif
   if( tag == "num_g5_c_cc" )then
     read(cread,*) iread
     num_g5_c_cc = iread
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
   if( tag == "io_mix" )then
     read(cread,*) iread
     io_mix = iread
   endif


 enddo ! i = 1, counter


 close(98) ! input_nnp.dat


 write(*,*) " "
 write(*,*) " Parameters for symmetry function"
 write(*,*) " ============================================================"

 !1 element
 !A
 if( elem_a /= "none" .and. elem_b == "none" .and. elem_c == "none" )then

   write(*,*) " atom_a: ", elem_a
   write(*,*) " input_node: ", node_in_a
   write(*,*) " num_g2"
   write(*,*) " a-a: ", num_g2_a_a
   write(*,*) " num_g5"
   write(*,*) " a-aa: ", num_g5_a_aa
   write(*,*) " "

 !B
 elseif( elem_a == "none" .and. elem_b /= "none" .and. elem_c == "none" )then

   write(*,*) " atom_b: ", elem_b
   write(*,*) " input_node: ", node_in_b
   write(*,*) " num_g2"
   write(*,*) " b-b: ", num_g2_b_b
   write(*,*) " num_g5"
   write(*,*) " b-bb: ", num_g5_b_bb
   write(*,*) " "

 !C
 elseif( elem_a == "none" .and. elem_b == "none" .and. elem_c /= "none" )then

   write(*,*) " atom_c: ", elem_c
   write(*,*) " input_node: ", node_in_c
   write(*,*) " num_g2"
   write(*,*) " c-c: ", num_g2_c_c
   write(*,*) " num_g5"
   write(*,*) " c-cc: ", num_g5_c_cc
   write(*,*) " "

 endif


 !2 elements
 !A,B
 if( elem_a /= "none" .and. elem_b /= "none" .and. elem_c == "none" )then

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
 elseif( elem_a /= "none" .and. elem_b == "none" .and. elem_c /= "none" )then

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

 !B,C
 elseif( elem_a == "none" .and. elem_b /= "none" .and. elem_c /= "none" )then

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

 endif


 !3 elements
 !A,B,C
 if( elem_a /= "none" .and. elem_b /= "none" .and. elem_c /= "none" )then

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


end subroutine read_input_nnp3
