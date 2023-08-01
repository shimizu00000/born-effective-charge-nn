subroutine read_io3( natom, nelem, &
                     elem_a, elem_b, elem_c, &
                     natom_a, natom_b, natom_c, &
                     io_a, io_b, io_c )
implicit none
integer natom, nelem
integer natom_a, natom_b, natom_c
integer natom_1, natom_2, natom_3
integer io_a, io_b, io_c
character(len=4) elem_a, elem_b, elem_c
character(len=4) elem_1, elem_2, elem_3


!Number of atoms and elements
 read(1) natom
 read(1) nelem

 !write(*,*) " "
 !write(*,*) " Tatal atoms:",    natom
 !write(*,*) " Total elements:", nelem
 !write(*,*) " "


!Kind of element & number of atoms
 natom_1 = 0 ; elem_1 = "empt"
 natom_2 = 0 ; elem_2 = "empt"
 natom_3 = 0 ; elem_3 = "empt"

 natom_a = 0
 natom_b = 0
 natom_c = 0

 io_a = 0
 io_b = 0
 io_c = 0

! 1 element
 if( nelem >= 1 )then

   read(1) elem_1, natom_1

   if( elem_1 == elem_a )then
     io_a = 1
     natom_a = natom_1
   elseif( elem_1 == elem_b )then
     io_b = 1
     natom_b = natom_1
   elseif( elem_1 == elem_c )then
     io_c = 1
     natom_c = natom_1
   else
     write(*,*) "Error in main-1 : revise program !!"
     stop
   endif

 endif ! 1 element

! 2 elements
 if( nelem >= 2 )then

   read(1) elem_2, natom_2

   if( elem_2 == elem_a )then
     io_a = 1
     natom_a = natom_2
   elseif( elem_2 == elem_b )then
     io_b = 1
     natom_b = natom_2
   elseif( elem_2 == elem_c )then
     io_c = 1
     natom_c = natom_2
   else
     write(*,*) "Error in main-2 : revise program !!"
     stop
   endif

endif ! 2 elements

! 3 elements
 if( nelem >= 3 )then

   read(1) elem_3, natom_3

   if( elem_3 == elem_a )then
     io_a = 1
     natom_a = natom_3
   elseif( elem_3 == elem_b )then
     io_b = 1
     natom_b = natom_3
   elseif( elem_3 == elem_c )then
     io_c = 1
     natom_c = natom_3
   else
     write(*,*) "Error in main-3 : revise program !!"
     stop
   endif

 endif ! 3 elements

!others
 if( nelem >= 4 )then
   write(*,*) "Error in main-5 : revise program !!"
   stop
 endif


end subroutine read_io3
