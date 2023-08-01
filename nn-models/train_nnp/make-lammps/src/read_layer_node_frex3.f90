subroutine read_layer_node_frex3( &
           network_a, network_b, network_c, &  
           nlayer_a, nlayer_b, nlayer_c, &
           io_nnp_a, io_nnp_b, io_nnp_c )

implicit none

character(len=120) tag
character(len=120) cread
integer            iread
!double precision   dread

integer i, stat, counter
character(len=4) elem_a, elem_b, elem_c

integer nlayer_a, nlayer_b, nlayer_c
integer network_a( nlayer_a ), &
        network_b( nlayer_b ), &
        network_c( nlayer_c )
!A
integer node_in_a, node_out_a
integer node_h1_a, node_h2_a, node_h3_a, node_h4_a, node_h5_a, &
        node_h6_a, node_h7_a, node_h8_a, node_h9_a, node_h10_a
!B
integer node_in_b, node_out_b, &
        node_h1_b, node_h2_b, node_h3_b, node_h4_b, node_h5_b, &
        node_h6_b, node_h7_b, node_h8_b, node_h9_b, node_h10_b
!C
integer node_in_c, node_out_c, &
        node_h1_c, node_h2_c, node_h3_c, node_h4_c, node_h5_c, &
        node_h6_c, node_h7_c, node_h8_c, node_h9_c, node_h10_c

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
 network_a = 0
 node_in_a = 0 ; node_out_a = 0
 node_h1_a = 0 ; node_h2_a = 0 ; node_h3_a = 0 ; node_h4_a = 0 ; node_h5_a = 0
 node_h6_a = 0 ; node_h7_a = 0 ; node_h8_a = 0 ; node_h9_a = 0 ; node_h10_a = 0
!B
 elem_b = "none"
 io_nnp_b = 0
 network_b = 0
 node_in_b = 0 ; node_out_b = 0 
 node_h1_b = 0 ; node_h2_b = 0 ; node_h3_b = 0 ; node_h4_b = 0 ; node_h5_b = 0 
 node_h6_b = 0 ; node_h7_b = 0 ; node_h8_b = 0 ; node_h9_b = 0 ; node_h10_b = 0
!C
 elem_c = "none"
 io_nnp_c = 0
 network_c = 0
 node_in_c = 0 ; node_out_c = 0 
 node_h1_c = 0 ; node_h2_c = 0 ; node_h3_c = 0 ; node_h4_c = 0 ; node_h5_c = 0
 node_h6_c = 0 ; node_h7_c = 0 ; node_h8_c = 0 ; node_h9_c = 0 ; node_h10_c = 0


 rewind(98)
 do i = 1, counter

   read(98,*) tag, dummyc, cread

   ! Parameters for A
   !------------------
   if( tag == "element_a" )then
     elem_a = cread
   endif
   if( tag == "nlayer_a" )then
     read(cread,*) iread
     nlayer_a    = iread
   endif
   if( tag == "nnp_io_a" )then
     read(cread,*) iread
     io_nnp_a    = iread
   endif

   if( tag == "node_in_a" )then
     read(cread,*)  iread
     network_a(1) = iread
   endif
   if( tag == "node_out_a" )then
     read(cread,*) iread
     network_a(nlayer_a) = iread
   endif

   if( tag == "node_h1_a" )then
     read(cread,*)  iread
     network_a(2) = iread
   endif
   if( tag == "node_h2_a" )then
     read(cread,*)  iread
     network_a(3) = iread
   endif
   if( tag == "node_h3_a" )then
     read(cread,*)  iread
     network_a(4) = iread
   endif
   if( tag == "node_h4_a" )then
     read(cread,*)  iread
     network_a(5) = iread
   endif
   if( tag == "node_h5_a" )then
     read(cread,*)  iread
     network_a(6) = iread
   endif
   if( tag == "node_h6_a" )then
     read(cread,*)  iread
     network_a(7) = iread
   endif
   if( tag == "node_h7_a" )then
     read(cread,*)  iread
     network_a(8) = iread
   endif
   if( tag == "node_h8_a" )then
     read(cread,*)  iread
     network_a(9) = iread
   endif
   if( tag == "node_h9_a" )then
     read(cread,*)  iread
     network_a(10) = iread
   endif
   if( tag == "node_h10_a" )then
     read(cread,*)  iread
     network_a(11) = iread
   endif

   ! Parameters for B
   !------------------
   if( tag == "element_b" )then
     elem_b = cread
   endif
   if( tag == "nlayer_b" )then
     read(cread,*) iread
     nlayer_b    = iread
   endif
   if( tag == "nnp_io_b" )then
     read(cread,*) iread
     io_nnp_b    = iread
   endif

   if( tag == "node_in_b" )then
     read(cread,*)  iread
     network_b(1) = iread
   endif
   if( tag == "node_out_b" )then
     read(cread,*)  iread
     network_b(nlayer_b) = iread
   endif

   if( tag == "node_h1_b" )then
     read(cread,*)  iread
     network_b(2) = iread
   endif
   if( tag == "node_h2_b" )then
     read(cread,*)  iread
     network_b(3) = iread
   endif
   if( tag == "node_h3_b" )then
     read(cread,*)  iread
     network_b(4) = iread
   endif
   if( tag == "node_h4_b" )then
     read(cread,*)  iread
     network_b(5) = iread
   endif
   if( tag == "node_h5_b" )then
     read(cread,*)  iread
     network_b(6) = iread
   endif
   if( tag == "node_h6_b" )then
     read(cread,*)  iread
     network_b(7) = iread
   endif
   if( tag == "node_h7_b" )then
     read(cread,*)  iread
     network_b(8) = iread
   endif
   if( tag == "node_h8_b" )then
     read(cread,*)  iread
     network_b(9) = iread
   endif
   if( tag == "node_h9_b" )then
     read(cread,*)  iread
     network_b(10) = iread
   endif
   if( tag == "node_h10_b" )then
     read(cread,*)  iread
     network_b(11) = iread
   endif

   ! Parameters for C
   !------------------
   if( tag == "element_c" )then
     elem_c = cread
   endif
   if( tag == "nlayer_c" )then
     read(cread,*) iread
     nlayer_c    = iread
   endif
   if( tag == "nnp_io_c" )then
     read(cread,*) iread
     io_nnp_c    = iread
   endif

   if( tag == "node_in_c" )then
     read(cread,*)  iread
     network_c(1) = iread
   endif
   if( tag == "node_out_c" )then
     read(cread,*)  iread
     network_c(nlayer_c) = iread
   endif

   if( tag == "node_h1_c" )then
     read(cread,*)  iread
     network_c(2) = iread
   endif
   if( tag == "node_h2_c" )then
     read(cread,*)  iread
     network_c(3) = iread
   endif
   if( tag == "node_h3_c" )then
     read(cread,*)  iread
     network_c(4) = iread
   endif
   if( tag == "node_h4_c" )then
     read(cread,*)  iread
     network_c(5) = iread
   endif
   if( tag == "node_h5_c" )then
     read(cread,*)  iread
     network_c(6) = iread
   endif
   if( tag == "node_h6_c" )then
     read(cread,*)  iread
     network_c(7) = iread
   endif
   if( tag == "node_h7_c" )then 
     read(cread,*)  iread
     network_c(8) = iread
   endif
   if( tag == "node_h8_c" )then
     read(cread,*)  iread
     network_c(9) = iread
   endif
   if( tag == "node_h9_c" )then
     read(cread,*)  iread
     network_c(10) = iread
   endif
   if( tag == "node_h10_c" )then
     read(cread,*)  iread
     network_c(11) = iread
   endif


 enddo ! i = 1, counter

 close(98) ! input_nnp.dat


 write(*,*) " "
 write(*,*) " Parameters for neural network"
 write(*,*) " ===================================="

 !A
 if( elem_a /= "none" )then
   write(*,fmt='(A10)',advance='no') " NNP(a): "
   do i = 1, nlayer_a - 1
     write(*,fmt='(I4,A2)',advance='no') network_a(i),"-"
   enddo
   write(*,fmt='(I4)',advance='no') network_a(nlayer_a)
   if( io_nnp_a == 0 ) write(*,fmt='(A4)') "OFF"
   if( io_nnp_a == 1 ) write(*,fmt='(A4)') "ON"
 endif
 !B
 if( elem_b /= "none" )then
   write(*,fmt='(A10)',advance='no') " NNP(b): "
   do i = 1, nlayer_b - 1
     write(*,fmt='(I4,A2)',advance='no') network_b(i),"-"
   enddo
   write(*,fmt='(I4)',advance='no') network_b(nlayer_b)
   if( io_nnp_b == 0 ) write(*,fmt='(A4)') "OFF"
   if( io_nnp_b == 1 ) write(*,fmt='(A4)') "ON"
 endif
 !C
 if( elem_c /= "none" )then
   write(*,fmt='(A10)',advance='no') " NNP(c): "
   do i = 1, nlayer_c - 1
     write(*,fmt='(I4,A2)',advance='no') network_c(i),"-"
   enddo
   write(*,fmt='(I4)',advance='no') network_c(nlayer_c)
   if( io_nnp_c == 0 ) write(*,fmt='(A4)') "OFF"
   if( io_nnp_c == 1 ) write(*,fmt='(A4)') "ON"
 endif

 write(*,*) " ===================================="
 write(*,*) " "


end subroutine read_layer_node_frex3
