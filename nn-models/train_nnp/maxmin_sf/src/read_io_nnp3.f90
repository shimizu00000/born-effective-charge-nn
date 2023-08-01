subroutine read_io_nnp3( io_nnp_a, io_nnp_b, io_nnp_c )

implicit none
character(len=100) tag
character(len=100) cread
integer            iread
double precision   dread

integer i, stat, counter
character(len=4) elem_a, elem_b, elem_c

integer io_nnp_a, io_nnp_b, io_nnp_c

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
 io_nnp_a = 0
!B
 elem_b = "none"
 io_nnp_b = 0
!C
 elem_c = "none"
 io_nnp_c = 0
 
 rewind(98)
 do i = 1, counter

   read(98,*) tag, dummyc, cread

   !parameters for A
   if( tag == "element_a" )then
     elem_a = cread
   endif
   if( tag == "nnp_io_a" )then
     read(cread,*) iread
     io_nnp_a    = iread
   endif

   !parameters for B
   if( tag == "element_b" )then
     elem_b = cread
   endif
   if( tag == "nnp_io_b" )then
     read(cread,*) iread
     io_nnp_b    = iread
   endif

   !parameters for C
   if( tag == "element_c" )then
     elem_c = cread
   endif
   if( tag == "nnp_io_c" )then
     read(cread,*) iread
     io_nnp_c    = iread
   endif

 enddo ! i = 1, counter


 close(98) ! input_nnp.dat


end subroutine read_io_nnp3
