program main_data_binary
!use parameters
implicit none
integer i, j, k, counter, stat, nfile, imd

integer   nelem, natom
integer   natom_a,  natom_b,  natom_c,  natom_d
character elem_a*4, elem_b*4, elem_c*4, elem_d*4
integer   md_steps

double precision lattice, x(3), y(3), z(3)
double precision,allocatable :: posi(:,:)
character dir*6, system*10, sel*9

double precision energy
double precision,allocatable :: born(:,:)

double precision t0, t1, t2, t_req
integer           dummyi
double precision  dummy
character(len=5)  dummyc

character(len=120) tag, filename

 1111 format( A120 )

 call cpu_time( t0 )


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
 write(*,*) counter, "tags are recognized"
 write(*,*) " "


!Tag loop
 rewind(99)
!************************************************************
! DO LOOP nfile
!************************************************************
 do nfile = 1, counter

   read(99,*) tag
   write(*,*) " "
   write(*,*) "#", nfile
   write(*,*) "------------------------------------------------------------"
   write(*,*) " "
   write(*,*) "Tag name: ", trim( tag )


   ! info READ
   write( filename, 1111 ) tag
   filename = '../step1_data/data_nnp/info_'//trim(adjustl(filename))//'.dat'
   open(unit=1,file=filename,action="read")
   ! positions READ
   write( filename, 1111 ) tag
   filename = '../step1_data/data_nnp/positions_'//trim(adjustl(filename))//'.dat'
   open(unit=2,file=filename,action="read")
   ! energies READ
   write( filename, 1111 ) tag
   filename = '../step1_data/data_nnp/energies_'//trim(adjustl(filename))//'.dat'
   open(unit=3,file=filename,action="read")
   ! born READ
   write( filename, 1111 ) tag
   filename = '../step1_data/data_nnp/born_'//trim(adjustl(filename))//'.dat'
   open(unit=4,file=filename,action="read")

   ! info (binary) WRITE
   write( filename, 1111 ) tag
   filename = './data/info_'//trim(adjustl(filename))//'.dat'
   open(unit=11,file=filename,action="write",form="unformatted")
   ! positions (binary) WRITE
   write( filename, 1111 ) tag
   filename = './data/positions_'//trim(adjustl(filename))//'.dat'
   open(unit=12,file=filename,action="write",form="unformatted")
   ! energies (binary) WRITE
   write( filename, 1111 ) tag
   filename = './data/energies_'//trim(adjustl(filename))//'.dat'
   open(unit=13,file=filename,action="write",form="unformatted")
   ! borns (binary) WRITE as force_*
   write( filename, 1111 ) tag
   filename = './data/forces_'//trim(adjustl(filename))//'.dat'
   open(unit=14,file=filename,action="write",form="unformatted")


   !******************************
   !Read parameters from info file
   !******************************
   !Total number of ions
   read(1,*) dummyc, natom
   write(*,*) " "
   write(*,*) "Total no. of atoms: ", natom
   write(11) natom
   !Ions per type
   read(1,*) dummyc, nelem
   write(*,*) "Total no. of elements: ", nelem
   write(11) nelem

   !Kind of element & number of atoms
   natom_a = 0     ; natom_b = 0     ; natom_c = 0     ; natom_d = 0
   elem_a = "none" ; elem_b = "none" ; elem_c = "none" ; elem_d = "none"

   if( nelem >= 1 ) read(1,*) elem_a, natom_a
   if( nelem >= 1 ) write(11) elem_a, natom_a

   if( nelem >= 2 ) read(1,*) elem_b, natom_b
   if( nelem >= 2 ) write(11) elem_b, natom_b

   if( nelem >= 3 ) read(1,*) elem_c, natom_c
   if( nelem >= 3 ) write(11) elem_c, natom_c

   if( nelem >= 4 ) read(1,*) elem_d, natom_d
   if( nelem >= 4 ) write(11) elem_d, natom_d

   if( nelem >= 5 ) write(*,*) "Error in main : revise program !!"
   if( nelem >= 5 ) stop

   write(*,*) " "
   write(*,*) " Ions per type: "
   write(*,*) elem_a, natom_a, elem_b, natom_b, elem_c, natom_c, elem_d, natom_d

   !MD steps
   read(1,*) dummyc, md_steps
   write(*,*) "MD steps: ", md_steps
   write(11) md_steps


   !*******************
   !Read positions file
   !*******************

   !Basis
   read(2,*) system
   read(2,*) lattice
   read(2,*) x(1), x(2), x(3)
   read(2,*) y(1), y(2), y(3)
   read(2,*) z(1), z(2), z(3)

   write(12) lattice
   write(12) x
   write(12) y
   write(12) z


   !Element & number of atoms
   if( nelem == 1 )then

     read(2,*) elem_a
     read(2,*) natom_a
     write(12) elem_a
     write(12) natom_a

   elseif( nelem == 2 )then

     read(2,*) elem_a,  elem_b
     read(2,*) natom_a, natom_b
     write(12) elem_a,  elem_b
     write(12) natom_a, natom_b

   elseif( nelem == 3 )then

     read(2,*) elem_a,  elem_b,  elem_c
     read(2,*) natom_a, natom_b, natom_c
     write(12) elem_a,  elem_b,  elem_c
     write(12) natom_a, natom_b, natom_c

   elseif( nelem == 4 )then

     read(2,*) elem_a,  elem_b,  elem_c,  elem_d
     read(2,*) natom_a, natom_b, natom_c, natom_d
     write(12) elem_a,  elem_b,  elem_c,  elem_d
     write(12) natom_a, natom_b, natom_c, natom_d

   else
     write(*,*) "Error in main : revise program !!" ; stop
     stop
   endif


   !Coordinations
   allocate( posi( natom, 3 ) )


   do imd = 1, md_steps

     !Direct configuration= ***
     read(2,*) dummyc

     !cartesian coordinates
     do i = 1, natom
       read(2,*) posi(i,1), posi(i,2), posi(i,3)
     enddo

     write(12) "Direct"
     write(12) posi

   enddo


   deallocate( posi )


   !******************
   !Read energies file
   !******************

   do i = 1, md_steps

     read(3,*) energy
     write(13) energy

   enddo


   !****************
   !Read born file
   !****************

   allocate( born( natom, 3 ) )

   do imd = 1, md_steps

     !--
     read(4,*) dummyc
     !POSITION, TOTAL-FORCE (eV/Angst)
     read(4,*) dummyc
     !--------------------------------
     read(4,*) dummyc
     !force_x,y,z for each atom
     do i = 1, natom
       read(4,*) dummy, dummy, dummy, born(i,1), born(i,2), born(i,3)
     enddo

     write(14) born

   enddo


   deallocate( born )


   close(1) ! info step1
   close(2) ! positions step1
   close(3) ! energies step1
   close(4) ! born step1

   close(11) ! info step2
   close(12) ! positions step2
   close(13) ! energies step2
   close(14) ! born step2


   write(*,*) " "


 enddo ! nfile, tag loop


 close(99) ! input_tag.dat


 write(*,*) "----------------------------"
 call cpu_time( t2 )
 t_req = t2 - t0
 write(*,*) "Total CPU time", t_req, "[s]"


end
