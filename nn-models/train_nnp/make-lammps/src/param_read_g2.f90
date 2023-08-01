subroutine param_read_g2( tag_in, g2_maxmin, eta2, rs2, num_g2 )
!use parameters
implicit none
integer num_g2
double precision g2_maxmin( num_g2, 2 ), eta2( num_g2 ), rs2( num_g2 )
integer i, stat, iparam, counter_g2
character(len=6) tag, tag_in
character(len=20) dummyc
character(len=100) filename

 filename = './data_maxmin/'//trim(adjustl(tag_in))//'_maxmin.dat'
 open(unit=97,file=filename,action="read",form="unformatted")
 read(97)g2_maxmin
 close(97)


 open(unit=97,file="./param_nnp.dat",action="read")
 !Number of lines in input_tag.dat
 counter_g2 = 0
 rewind(97)
 do
   read( 97, '(A20)', iostat=stat ) dummyc
   if( stat /= 0 ) exit
!   if( trim( dummyc ) == 'g' ) counter_g2 = counter_g2 + 1
!   if( dummyc(1:1) == 'g' ) counter_g2 = counter_g2 + 1
   if( dummyc(1:10) /= '' ) counter_g2 = counter_g2 + 1
 enddo

 rewind(97)
 do iparam = 1, counter_g2

   read(97,*) tag

   if( trim(adjustl(tag)) == trim(adjustl(tag_in)) )then

       read(97,*) dummyc

     do i = 1, num_g2

       read(97,*) eta2(i), rs2(i)

     enddo

     rewind(97)
     cycle

   endif

 enddo ! iparam

 close(97)


end subroutine param_read_g2
