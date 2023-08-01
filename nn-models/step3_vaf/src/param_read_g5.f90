subroutine param_read_g5( tag_in, eta5, theta5, zeta5, lambda5, num_g5 )
!use parameters
implicit none
integer num_g5
double precision eta5( num_g5 ), theta5( num_g5 ), zeta5( num_g5 ), lambda5( num_g5 )
integer i, stat, iparam, counter_g5
character(len=7) tag, tag_in
character(len=20) dummyc


 open(unit=97,file="./param_nnp.dat",action="read")
 !Number of lines in input_tag.dat
 counter_g5 = 0
 rewind(97)
 do
   read( 97, '(A20)', iostat=stat ) dummyc
   if( stat /= 0 ) exit
!   if( trim( dummyc ) == 'g' ) counter_g5 = counter_g5 + 1
!   if( dummyc(1:1) == 'g' ) counter_g5 = counter_g5 + 1
   if( dummyc(1:10) /= '' ) counter_g5 = counter_g5 + 1
 enddo

 rewind(97)
 do iparam = 1, counter_g5

   read(97,*) tag

   if( trim(adjustl(tag)) == trim(adjustl(tag_in)) )then

       read(97,*) dummyc

     do i = 1, num_g5

       read(97,*) eta5(i), theta5(i), zeta5(i), lambda5(i)

     enddo

     rewind(97)
     cycle

   endif

 enddo ! iparam


 close(97)


end subroutine param_read_g5
