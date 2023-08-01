subroutine read_layer_node3( nlayer_a, nlayer_b, nlayer_c )

implicit none
character(len=120) tag
character(len=120) cread
integer            iread
!double precision   dread

integer i, stat, counter
character(len=4) elem_a, elem_b, elem_c
integer nlayer_a, nlayer_b, nlayer_c

character dummyc


!***********************
!Read input_nnp.dat file
!***********************

 open(unit=98,file="input_nnp.dat",action="read")
!Number of lines in input_nnp.dat
 counter = 0
 rewind(98)
 do
   read( 98, '(A1)', iostat=stat ) dummyc
   if( stat /= 0 ) exit
   if( trim( dummyc ) /= '' ) counter = counter + 1
 enddo

!A
 elem_a = "none"
 nlayer_a = 0
!B
 elem_b = "none"
 nlayer_b = 0
!C
 elem_c = "none"
 nlayer_c = 0
 
 rewind(98)
 do i = 1, counter

   read(98,*) tag, dummyc, cread

   !parameters for A
   if( tag == "element_a" )then
     elem_a = cread
   endif
   if( tag == "nlayer_a" )then
     read(cread,*) iread
     nlayer_a    = iread
   endif

   !parameters for B
   if( tag == "element_b" )then
     elem_b = cread
   endif
   if( tag == "nlayer_b" )then
     read(cread,*) iread
     nlayer_b    = iread
   endif

   !parameters for C
   if( tag == "element_c" )then
     elem_c = cread
   endif
   if( tag == "nlayer_c" )then
     read(cread,*) iread
     nlayer_c    = iread
   endif

 enddo ! i = 1, counter


 close(98) ! input_nnp.dat


end subroutine read_layer_node3
