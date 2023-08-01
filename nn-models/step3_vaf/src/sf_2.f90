subroutine sf_2( posi, xlattice, ylattice, zlattice, &
                 g2_x_x,  eta2_x_x,  rs2_x_x, &
                 g2_x_y,  eta2_x_y,  rs2_x_y, &
                 g5_x_xx, eta5_x_xx, theta5_x_xx, zeta5_x_xx, lambda5_x_xx, &
                 g5_x_xy, eta5_x_xy, theta5_x_xy, zeta5_x_xy, lambda5_x_xy, &
                 g5_x_yy, eta5_x_yy, theta5_x_yy, zeta5_x_yy, lambda5_x_yy, &
                 g2_y_x,  eta2_y_x,  rs2_y_x, &
                 g2_y_y,  eta2_y_y,  rs2_y_y, &
                 g5_y_xx, eta5_y_xx, theta5_y_xx, zeta5_y_xx, lambda5_y_xx, &
                 g5_y_xy, eta5_y_xy, theta5_y_xy, zeta5_y_xy, lambda5_y_xy, &
                 g5_y_yy, eta5_y_yy, theta5_y_yy, zeta5_y_yy, lambda5_y_yy, &
                 g2_deriv_x_x,  g2_deriv_x_y, &
                 g5_deriv_x_xx, g5_deriv_x_xy, g5_deriv_x_yy, &
                 g2_deriv_y_x,  g2_deriv_y_y, &
                 g5_deriv_y_xx, g5_deriv_y_xy, g5_deriv_y_yy, &
                 num_g2_x_x,  num_g2_x_y, &
                 num_g5_x_xx, num_g5_x_xy, num_g5_x_yy, &
                 num_g2_y_x,  num_g2_y_y, &
                 num_g5_y_xx, num_g5_y_xy, num_g5_y_yy, &
                 natom_x, natom_y,  &
                 natom, rc, neighbor, io_sf, num_mpi, myrank )
implicit none
include 'mpif.h'
integer ierr, num_mpi, myrank
integer i, j, k, l, m, n, ii, ig5, i0
integer ix1, iy1, iz1, ix2, lattice_box, boxsize, boxcenter, neighbor, center
integer natom, atom1, atom2 !, time, jamp
integer natom_x, natom_y
! natom : total number of atoms
double precision xlattice, ylattice, zlattice
double precision posi( natom, 3 ), posi1( natom, 3 )
double precision,allocatable :: posi2( : , : ) 
integer list(natom,neighbor)

integer num_g2_x_x,  num_g2_x_y, &
        num_g5_x_xx, num_g5_x_xy, num_g5_x_yy
integer num_g2_y_x,  num_g2_y_y, &
        num_g5_y_xx, num_g5_y_xy, num_g5_y_yy

!X
double precision &
  g2_x_x( natom_x, num_g2_x_x ), eta2_x_x( num_g2_x_x ), rs2_x_x( num_g2_x_x ), &
  g2_x_y( natom_x, num_g2_x_y ), eta2_x_y( num_g2_x_y ), rs2_x_y( num_g2_x_y )
double precision &
  g5_x_xx( natom_x, num_g5_x_xx ), &
   eta5_x_xx( num_g5_x_xx ),  theta5_x_xx( num_g5_x_xx ), &
  zeta5_x_xx( num_g5_x_xx ), lambda5_x_xx( num_g5_x_xx ), &
  g5_x_xy( natom_x, num_g5_x_xy ), &
   eta5_x_xy( num_g5_x_xy ),  theta5_x_xy( num_g5_x_xy ), &
  zeta5_x_xy( num_g5_x_xy ), lambda5_x_xy( num_g5_x_xy ), &
  g5_x_yy( natom_x, num_g5_x_yy ), &
   eta5_x_yy( num_g5_x_yy ),  theta5_x_yy( num_g5_x_yy ), &
  zeta5_x_yy( num_g5_x_yy ), lambda5_x_yy( num_g5_x_yy )
double precision &
  g2_x_x_mpi(  natom_x, num_g2_x_x ), &
  g2_x_y_mpi(  natom_x, num_g2_x_y ), &
  g5_x_xx_mpi( natom_x, num_g5_x_xx ), &
  g5_x_xy_mpi( natom_x, num_g5_x_xy ), &
  g5_x_yy_mpi( natom_x, num_g5_x_yy )
!Y
double precision &
  g2_y_x( natom_y, num_g2_y_x ), eta2_y_x( num_g2_y_x ), rs2_y_x( num_g2_y_x ), &
  g2_y_y( natom_y, num_g2_y_y ), eta2_y_y( num_g2_y_y ), rs2_y_y( num_g2_y_y )
double precision &
  g5_y_xx( natom_y, num_g5_y_xx ), &
   eta5_y_xx( num_g5_y_xx ),  theta5_y_xx( num_g5_y_xx ), &
  zeta5_y_xx( num_g5_y_xx ), lambda5_y_xx( num_g5_y_xx ), &
  g5_y_xy( natom_y, num_g5_y_xy ), &
   eta5_y_xy( num_g5_y_xy ),  theta5_y_xy( num_g5_y_xy ), &
  zeta5_y_xy( num_g5_y_xy ), lambda5_y_xy( num_g5_y_xy ), &
  g5_y_yy( natom_y, num_g5_y_yy ), &
   eta5_y_yy( num_g5_y_yy ),  theta5_y_yy( num_g5_y_yy ), &
  zeta5_y_yy( num_g5_y_yy ), lambda5_y_yy( num_g5_y_yy )
double precision &
  g2_y_x_mpi(  natom_y, num_g2_y_x ), &
  g2_y_y_mpi(  natom_y, num_g2_y_y ), &
  g5_y_xx_mpi( natom_y, num_g5_y_xx ), &
  g5_y_xy_mpi( natom_y, num_g5_y_xy ), &
  g5_y_yy_mpi( natom_y, num_g5_y_yy )
 
double precision &
  g2_deriv_x_x(  natom_x, natom, 3, num_g2_x_x  ), &
  g2_deriv_x_y(  natom_x, natom, 3, num_g2_x_y  ), &
  g5_deriv_x_xx( natom_x, natom, 3, num_g5_x_xx ), &
  g5_deriv_x_xy( natom_x, natom, 3, num_g5_x_xy ), &
  g5_deriv_x_yy( natom_x, natom, 3, num_g5_x_yy )
double precision &
  g2_deriv_x_x_mpi(  natom_x, natom, 3, num_g2_x_x  ), &
  g2_deriv_x_y_mpi(  natom_x, natom, 3, num_g2_x_y  ), &
  g5_deriv_x_xx_mpi( natom_x, natom, 3, num_g5_x_xx ), &
  g5_deriv_x_xy_mpi( natom_x, natom, 3, num_g5_x_xy ), &
  g5_deriv_x_yy_mpi( natom_x, natom, 3, num_g5_x_yy )
double precision &
  g2_deriv_y_x(  natom_y, natom, 3, num_g2_y_x  ), &
  g2_deriv_y_y(  natom_y, natom, 3, num_g2_y_y  ), &
  g5_deriv_y_xx( natom_y, natom, 3, num_g5_y_xx ), &
  g5_deriv_y_xy( natom_y, natom, 3, num_g5_y_xy ), &
  g5_deriv_y_yy( natom_y, natom, 3, num_g5_y_yy )
double precision &
  g2_deriv_y_x_mpi(  natom_y, natom, 3, num_g2_y_x  ), &
  g2_deriv_y_y_mpi(  natom_y, natom, 3, num_g2_y_y  ), &
  g5_deriv_y_xx_mpi( natom_y, natom, 3, num_g5_y_xx ), &
  g5_deriv_y_xy_mpi( natom_y, natom, 3, num_g5_y_xy ), &
  g5_deriv_y_yy_mpi( natom_y, natom, 3, num_g5_y_yy )

double precision Rij, Rik, g5_const_3
double precision rc
integer io_sf


 do i = 1, natom
   posi1(i,1) = posi(i,1) * xlattice
   posi1(i,2) = posi(i,2) * ylattice
   posi1(i,3) = posi(i,3) * zlattice
 enddo


!Repeating box
 if( min( xlattice, ylattice, zlattice ) >= rc )then
   lattice_box = 27
   boxsize     = 3
   boxcenter   = 2
   center      = natom * 13
 elseif( rc > min( xlattice, ylattice, zlattice ) &
         .and. 2.0d0 * min( xlattice, ylattice, zlattice ) >= rc )then
   lattice_box = 125
   boxsize     = 5
   boxcenter   = 3
   center      = natom * 62
 elseif( rc > 2.0d0 * min( xlattice, ylattice, zlattice ) &
         .and. 3.0d0 * min( xlattice, ylattice, zlattice ) >= rc )then
   lattice_box = 343
   boxsize     = 7
   boxcenter   = 4
   center      = natom * 171
 else
   write(*,*) "Error(sf_2): lattice_box"
   stop
 endif

 if( myrank == 0 .and. io_sf == 0 )then 
   write(*,*) " "
   write(*,*) " Repeating box size: ", lattice_box
   write(*,*) " "
   io_sf = 1
 endif


 allocate( posi2( lattice_box * natom, 3 ) )


 m = 0
 do ix1 = 1, boxsize
   do iy1 = 1, boxsize
     do iz1 = 1, boxsize
       do i = 1, natom
         m = m + 1
         posi2(m,1) = posi1(i,1) + xlattice*( ix1 - boxcenter )
         posi2(m,2) = posi1(i,2) + ylattice*( iy1 - boxcenter )
         posi2(m,3) = posi1(i,3) + zlattice*( iz1 - boxcenter )
       enddo
     enddo
   enddo
 enddo


!Neighboring list
!if( mod(time,neighbor_check)==1 .or. jamp == 1 )then
     list = 0
  do i = 1, natom
     m = 0
     n = 0
    do ix1 = 1, boxsize
      do iy1 = 1, boxsize
        do iz1 = 1, boxsize
          do j = 1, natom

            m = m + 1

            Rij = dsqrt(   ( posi1(i,1) - posi2(m,1) )**2 &
                         + ( posi1(i,2) - posi2(m,2) )**2 &
                         + ( posi1(i,3) - posi2(m,3) )**2 )

            if( Rij <= rc )then !+ Rn )then

              n = n + 1

              list(i,n) = m

              if( n == neighbor )then
                write(*,*) "Insufficient neighbor list !!"
                stop
              endif

            endif
          enddo
        enddo
      enddo
    enddo
  enddo
!endif


   g2_x_x  = 0.0d0 ; g2_x_y  = 0.0d0 
   g5_x_xx = 0.0d0 ; g5_x_xy = 0.0d0 ; g5_x_yy = 0.0d0
   g2_y_x  = 0.0d0 ; g2_y_y  = 0.0d0
   g5_y_xx = 0.0d0 ; g5_y_xy = 0.0d0 ; g5_y_yy = 0.0d0

   g2_deriv_x_x  = 0.0d0 ; g2_deriv_x_y  = 0.0d0
   g5_deriv_x_xx = 0.0d0 ; g5_deriv_x_xy = 0.0d0 ; g5_deriv_x_yy = 0.0d0
   g2_deriv_y_x  = 0.0d0 ; g2_deriv_y_y  = 0.0d0
   g5_deriv_y_xx = 0.0d0 ; g5_deriv_y_xy = 0.0d0 ; g5_deriv_y_yy = 0.0d0

! i : taget atom
! j : 1st neighbor atom
! k : 2nd neighbor atom
 do i0 = 1, natom, num_mpi
   i = i0 + myrank
   if( i <= natom )then
 !do i = 1, natom ! Target atom

   ! G2_X_... & G5_X_...
   if( i <= natom_x )then

     do ix1 = 1, neighbor

       atom1 = list(i,ix1)

       if( atom1 == 0 )cycle
       if( atom1 - center == i )cycle

       j = atom1 - natom * ( ( atom1 - 1 ) / natom )
       Rij = dsqrt(   ( posi1(i,1) - posi2(atom1,1) )**2 &
                    + ( posi1(i,2) - posi2(atom1,2) )**2 &
                    + ( posi1(i,3) - posi2(atom1,3) )**2 )

       if( Rij < rc )then

         ! G2_X_X
         if( j <= natom_x )then
           call sf_2_g2( g2_x_x, g2_deriv_x_x, eta2_x_x, rs2_x_x, num_g2_x_x, &
                         posi1(i,1),     posi1(i,2),     posi1(i,3), &
                         posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                         natom, natom_x, natom_y, &
                         Rij, i, j, rc, natom_x )
         ! G2_X_Y
         elseif( natom_x < j .and. j <= natom_x + natom_y )then
           call sf_2_g2( g2_x_y, g2_deriv_x_y, eta2_x_y, rs2_x_y, num_g2_x_y, &
                         posi1(i,1),     posi1(i,2),     posi1(i,3), &
                         posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                         natom, natom_x, natom_y, &
                         Rij, i, j, rc, natom_x )

         else
           write(*,*) "Error(sf_2): G2_X_*"
           stop
         endif ! j

         do ix2 = 1, neighbor

           atom2 = list(i,ix2)

           if( atom2 == 0 )cycle
           if( atom2 - center == i )cycle
           if( atom1 == atom2 )cycle

           k = atom2 - natom * ( ( atom2 - 1 ) / natom )
           Rik = dsqrt(   ( posi1(i,1) - posi2(atom2,1) )**2 &
                        + ( posi1(i,2) - posi2(atom2,2) )**2 &
                        + ( posi1(i,3) - posi2(atom2,3) )**2 )

           if( Rik <= rc )then
             
             ! G5_X_XX
             if( j <= natom_x .and. &
                 k <= natom_x )then
               call sf_2_g5( g5_x_xx, g5_deriv_x_xx, &
                             eta5_x_xx, lambda5_x_xx, zeta5_x_xx, num_g5_x_xx, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, &
                             Rij, Rik, i, j, k, rc, natom_x )
             ! G5_X_XY
             elseif( j <= natom_x .and. &
                     natom_x < k .and. k <= natom_x + natom_y )then
               call sf_2_g5( g5_x_xy, g5_deriv_x_xy, &
                             eta5_x_xy, lambda5_x_xy, zeta5_x_xy, num_g5_x_xy, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, &
                             Rij, Rik, i, j, k, rc, natom_x )
             ! G5_X_YY
             elseif( natom_x < j .and. j <= natom_x + natom_y .and. &
                     natom_x < k .and. k <= natom_x + natom_y )then
               call sf_2_g5( g5_x_yy, g5_deriv_x_yy, &
                             eta5_x_yy, lambda5_x_yy, zeta5_x_yy, num_g5_x_yy, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, &
                             Rij, Rik, i, j, k, rc, natom_x )
             ! G5_X_YX
             elseif( natom_x < j .and. j <= natom_x + natom_y .and. &
                     k <= natom_x )then
               cycle
             else
               write(*,*) "Error(sf_2): G5_X_*"
               stop
             endif ! j,k

           endif ! if( Rik <= rc )
         enddo ! do ix2
       endif ! if( Rij < rc )
     enddo ! do ix1 (g2_x... & g5_x...)

   ! G2_Y_... & G5_Y_...
   elseif( natom_x < i .and. i <= natom_x + natom_y )then

     do ix1 = 1, neighbor

       atom1 = list(i,ix1)

       if( atom1 == 0 )cycle
       if( atom1 - center == i )cycle

       j = atom1 - natom * ( ( atom1 - 1 ) / natom )
       Rij = dsqrt(   ( posi1(i,1) - posi2(atom1,1) )**2 &
                    + ( posi1(i,2) - posi2(atom1,2) )**2 &
                    + ( posi1(i,3) - posi2(atom1,3) )**2 )

       if( Rij < rc )then
         ! G2_Y_X
         if( j <= natom_x )then
           call sf_2_g2( g2_y_x, g2_deriv_y_x, eta2_y_x, rs2_y_x, num_g2_y_x, &
                         posi1(i,1),     posi1(i,2),     posi1(i,3), &
                         posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                         natom, natom_x, natom_y, &
                         Rij, i, j, rc, natom_y )
         ! G2_Y_Y
         elseif( natom_x < j .and. j <= natom_x + natom_y )then
           call sf_2_g2( g2_y_y, g2_deriv_y_y, eta2_y_y, rs2_y_y, num_g2_y_y, &
                         posi1(i,1),     posi1(i,2),     posi1(i,3), &
                         posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                         natom, natom_x, natom_y, &
                         Rij, i, j, rc, natom_y )
         else
           write(*,*) "Error(sf_2): G2_Y_*"
           stop
         endif

         do ix2 = 1, neighbor

           atom2 = list(i,ix2)

           if( atom2 == 0 )cycle
           if( atom2 - center == i )cycle
           if( atom1 == atom2 )cycle

           k = atom2 - natom * ( ( atom2 - 1 ) / natom )
           Rik = dsqrt(   ( posi1(i,1) - posi2(atom2,1) )**2 &
                        + ( posi1(i,2) - posi2(atom2,2) )**2 &
                        + ( posi1(i,3) - posi2(atom2,3) )**2 )

           if( Rik <= rc )then
             
             ! G5_Y_XX
             if( j <= natom_x .and. &
                 k <= natom_x )then
               call sf_2_g5( g5_y_xx, g5_deriv_y_xx, &
                             eta5_y_xx, lambda5_y_xx, zeta5_y_xx, num_g5_y_xx, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, &
                             Rij, Rik, i, j, k, rc, natom_y )
             ! G5_Y_XY
             elseif( j <= natom_x .and. &
                     natom_x < k .and. k <= natom_x + natom_y )then
               call sf_2_g5( g5_y_xy, g5_deriv_y_xy, &
                             eta5_y_xy, lambda5_y_xy, zeta5_y_xy, num_g5_y_xy, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, &
                             Rij, Rik, i, j, k, rc, natom_y )
             ! G5_Y_YY
             elseif( natom_x < j .and. j <= natom_x + natom_y .and. &
                     natom_x < k .and. k <= natom_x + natom_y )then
               call sf_2_g5( g5_y_yy, g5_deriv_y_yy, &
                             eta5_y_yy, lambda5_y_yy, zeta5_y_yy, num_g5_y_yy, &
                             posi1(i,1),     posi1(i,2),     posi1(i,3), &
                             posi2(atom1,1), posi2(atom1,2), posi2(atom1,3), &
                             posi2(atom2,1), posi2(atom2,2), posi2(atom2,3), &
                             natom, natom_x, natom_y, &
                             Rij, Rik, i, j, k, rc, natom_y )
             ! G5_Y_YX
             elseif( natom_x < j .and. j <= natom_x + natom_y .and. &
                     k <= natom_x )then
               cycle
             else
               write(*,*) "Error(sf_2): G5_Y_*"
               stop
             endif ! j,k

           endif ! if( Rik <= rc )
         enddo ! do ix2
       endif ! if( Rij < rc )
     enddo ! do ix1 (g2_y... & g5_y...)

   endif ! if( i <= natom_x )

   endif ! i <= natom
 enddo ! i0

 call mpi_barrier( mpi_comm_world, ierr )

 
!*****************
!Coefficient of g5
!*****************
 do i0 = 1, natom, num_mpi
   i = i0 + myrank
   if( i <= natom )then

!  A
   if( i <= natom_x )then

     do ig5 = 1, num_g5_x_xx
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_x_xx(ig5) )
         g5_x_xx(i,ig5) = g5_x_xx(i,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_x_xx(i,j,1,ig5) = g5_deriv_x_xx(i,j,1,ig5) * g5_const_3
         g5_deriv_x_xx(i,j,2,ig5) = g5_deriv_x_xx(i,j,2,ig5) * g5_const_3
         g5_deriv_x_xx(i,j,3,ig5) = g5_deriv_x_xx(i,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_x_xy
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_x_xy(ig5) )
         g5_x_xy(i,ig5) = g5_x_xy(i,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_x_xy(i,j,1,ig5) = g5_deriv_x_xy(i,j,1,ig5) * g5_const_3
         g5_deriv_x_xy(i,j,2,ig5) = g5_deriv_x_xy(i,j,2,ig5) * g5_const_3
         g5_deriv_x_xy(i,j,3,ig5) = g5_deriv_x_xy(i,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_x_yy
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_x_yy(ig5) )
         g5_x_yy(i,ig5) = g5_x_yy(i,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_x_yy(i,j,1,ig5) = g5_deriv_x_yy(i,j,1,ig5) * g5_const_3
         g5_deriv_x_yy(i,j,2,ig5) = g5_deriv_x_yy(i,j,2,ig5) * g5_const_3
         g5_deriv_x_yy(i,j,3,ig5) = g5_deriv_x_yy(i,j,3,ig5) * g5_const_3
       enddo
     enddo

!  B
   elseif( natom_x < i .and. i <= natom_x + natom_y )then

     ii = i - natom_x
  
     do ig5 = 1, num_g5_y_xx
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_y_xx(ig5) )
         g5_y_xx(ii,ig5) = g5_y_xx(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_y_xx(ii,j,1,ig5) = g5_deriv_y_xx(ii,j,1,ig5) * g5_const_3
         g5_deriv_y_xx(ii,j,2,ig5) = g5_deriv_y_xx(ii,j,2,ig5) * g5_const_3
         g5_deriv_y_xx(ii,j,3,ig5) = g5_deriv_y_xx(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_y_xy
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_y_xy(ig5) )
         g5_y_xy(ii,ig5) = g5_y_xy(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_y_xy(ii,j,1,ig5) = g5_deriv_y_xy(ii,j,1,ig5) * g5_const_3
         g5_deriv_y_xy(ii,j,2,ig5) = g5_deriv_y_xy(ii,j,2,ig5) * g5_const_3
         g5_deriv_y_xy(ii,j,3,ig5) = g5_deriv_y_xy(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

     do ig5 = 1, num_g5_y_yy
         g5_const_3 = 2.0d0**( 1.0d0 - zeta5_y_yy(ig5) )
         g5_y_yy(ii,ig5) = g5_y_yy(ii,ig5) * g5_const_3
       do j = 1, natom
         g5_deriv_y_yy(ii,j,1,ig5) = g5_deriv_y_yy(ii,j,1,ig5) * g5_const_3
         g5_deriv_y_yy(ii,j,2,ig5) = g5_deriv_y_yy(ii,j,2,ig5) * g5_const_3
         g5_deriv_y_yy(ii,j,3,ig5) = g5_deriv_y_yy(ii,j,3,ig5) * g5_const_3
       enddo
     enddo

   endif ! i

   endif ! i < natom

 enddo ! i0


 call mpi_barrier( mpi_comm_world, ierr )

 call mpi_reduce( g2_x_x, g2_x_x_mpi, natom_x*num_g2_x_x, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_x_y, g2_x_y_mpi, natom_x*num_g2_x_y, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_x_xx, g5_x_xx_mpi, natom_x*num_g5_x_xx, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_x_xy, g5_x_xy_mpi, natom_x*num_g5_x_xy, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_x_yy, g5_x_yy_mpi, natom_x*num_g5_x_yy, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )

 call mpi_reduce( g2_y_x, g2_y_x_mpi, natom_y*num_g2_y_x, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_y_y, g2_y_y_mpi, natom_y*num_g2_y_y, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_y_xx, g5_y_xx_mpi, natom_y*num_g5_y_xx, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_y_xy, g5_y_xy_mpi, natom_y*num_g5_y_xy, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_y_yy, g5_y_yy_mpi, natom_y*num_g5_y_yy, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 

 call mpi_reduce( g2_deriv_x_x, g2_deriv_x_x_mpi, natom_x*natom*3*num_g2_x_x, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_deriv_x_y, g2_deriv_x_y_mpi, natom_x*natom*3*num_g2_x_y, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_x_xx, g5_deriv_x_xx_mpi, natom_x*natom*3*num_g5_x_xx, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_x_xy, g5_deriv_x_xy_mpi, natom_x*natom*3*num_g5_x_xy, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_x_yy, g5_deriv_x_yy_mpi, natom_x*natom*3*num_g5_x_yy, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )

 call mpi_reduce( g2_deriv_y_x, g2_deriv_y_x_mpi, natom_y*natom*3*num_g2_y_x, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g2_deriv_y_y, g2_deriv_y_y_mpi, natom_y*natom*3*num_g2_y_y, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_y_xx, g5_deriv_y_xx_mpi, natom_y*natom*3*num_g5_y_xx, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_y_xy, g5_deriv_y_xy_mpi, natom_y*natom*3*num_g5_y_xy, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv_y_yy, g5_deriv_y_yy_mpi, natom_y*natom*3*num_g5_y_yy, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )

 call mpi_barrier( mpi_comm_world, ierr )


 deallocate( posi2 )


!--------
!Check SF
!--------

!X
do j = 1, num_g2_x_x
  do i = 1, natom_x
    if( isnan( g2_x_x(i,j) ) ) stop 'NaN g2_x_x'
  enddo
enddo

do j = 1, num_g2_x_y
  do i = 1, natom_x
    if( isnan( g2_x_y(i,j) ) ) stop 'NaN g2_x_y'
  enddo
enddo

do j = 1, num_g5_x_xx
  do i = 1, natom_x
    if( isnan( g5_x_xx(i,j) ) ) stop 'NaN g5_x_xx'
  enddo
enddo

do j = 1, num_g5_x_xy
  do i = 1, natom_x
    if( isnan( g5_x_xy(i,j) ) ) stop 'NaN g5_x_xy'
  enddo
enddo

do j = 1, num_g5_x_yy
  do i = 1, natom_x
    if( isnan( g5_x_yy(i,j) ) ) stop 'NaN g5_x_yy'
  enddo
enddo

!Y
do j = 1, num_g2_y_x
  do i = 1, natom_y
    if( isnan( g2_y_x(i,j) ) ) stop 'NaN g2_y_x'
  enddo
enddo

do j = 1, num_g2_y_y
  do i = 1, natom_y
    if( isnan( g2_y_y(i,j) ) ) stop 'NaN g2_y_y'
  enddo
enddo

do j = 1, num_g5_y_xx
  do i = 1, natom_y
    if( isnan( g5_y_xx(i,j) ) ) stop 'NaN g5_y_xx'
  enddo
enddo

do j = 1, num_g5_y_xy
  do i = 1, natom_y
    if( isnan( g5_y_xy(i,j) ) ) stop 'NaN g5_y_xy'
  enddo
enddo

do j = 1, num_g5_y_yy
  do i = 1, natom_y
    if( isnan( g5_y_yy(i,j) ) ) stop 'NaN g5_y_yy'
  enddo
enddo


!Deriv X
do l = 1, num_g2_x_x
do k = 1, 3
do j = 1, natom
do i = 1, natom_x
  if( isnan( g2_deriv_x_x(i,j,k,l) ) ) stop 'NaN g2_deriv_x_x'
enddo
enddo
enddo
enddo

do l = 1, num_g2_x_y
do k = 1, 3
do j = 1, natom
do i = 1, natom_x
  if( isnan( g2_deriv_x_y(i,j,k,l) ) ) stop 'NaN g2_deriv_x_y'
enddo
enddo
enddo
enddo

do l = 1, num_g5_x_xx
do k = 1, 3
do j = 1, natom
do i = 1, natom_x
  if( isnan( g5_deriv_x_xx(i,j,k,l) ) ) stop 'NaN g5_deriv_x_xx'
enddo
enddo
enddo
enddo

do l = 1, num_g5_x_xy
do k = 1, 3
do j = 1, natom
do i = 1, natom_x
  if( isnan( g5_deriv_x_xy(i,j,k,l) ) ) stop 'NaN g5_deriv_x_xy'
enddo
enddo
enddo
enddo

do l = 1, num_g5_x_yy
do k = 1, 3
do j = 1, natom
do i = 1, natom_x
  if( isnan( g5_deriv_x_yy(i,j,k,l) ) ) stop 'NaN g5_deriv_x_yy'
enddo
enddo
enddo
enddo

!Deriv Y
do l = 1, num_g2_y_x
do k = 1, 3
do j = 1, natom
do i = 1, natom_y
  if( isnan( g2_deriv_y_x(i,j,k,l) ) ) stop 'NaN g2_deriv_y_x'
enddo
enddo
enddo
enddo

do l = 1, num_g2_y_y
do k = 1, 3
do j = 1, natom
do i = 1, natom_y
  if( isnan( g2_deriv_y_y(i,j,k,l) ) ) stop 'NaN g2_deriv_y_y'
enddo
enddo
enddo
enddo

do l = 1, num_g5_y_xx
do k = 1, 3
do j = 1, natom
do i = 1, natom_y
  if( isnan( g5_deriv_y_xx(i,j,k,l) ) ) stop 'NaN g5_deriv_y_xx'
enddo
enddo
enddo
enddo

do l = 1, num_g5_y_xy
do k = 1, 3
do j = 1, natom
do i = 1, natom_y
  if( isnan( g5_deriv_y_xy(i,j,k,l) ) ) stop 'NaN g5_deriv_y_xy'
enddo
enddo
enddo
enddo

do l = 1, num_g5_y_yy
do k = 1, 3
do j = 1, natom
do i = 1, natom_y
  if( isnan( g5_deriv_y_yy(i,j,k,l) ) ) stop 'NaN g5_deriv_y_yy'
enddo
enddo
enddo
enddo


if( myrank == 0 )then

!--------
!Write sf
!--------

!g2_x-... & g5_x_...
 write(11)g2_x_x_mpi
 write(11)g2_x_y_mpi
 write(11)g5_x_xx_mpi
 write(11)g5_x_xy_mpi
 write(11)g5_x_yy_mpi

!g2_y-... & g5_y_...
 write(11)g2_y_x_mpi
 write(11)g2_y_y_mpi
 write(11)g5_y_xx_mpi
 write(11)g5_y_xy_mpi
 write(11)g5_y_yy_mpi


!--------------
!Write sf_deriv
!--------------

!g2_deriv_x-... & g5_deriv_x-...
 write(12)g2_deriv_x_x_mpi
 write(12)g2_deriv_x_y_mpi
 write(12)g5_deriv_x_xx_mpi
 write(12)g5_deriv_x_xy_mpi
 write(12)g5_deriv_x_yy_mpi

!g2_deriv_y-... & g5_deriv_y-...
 write(12)g2_deriv_y_x_mpi
 write(12)g2_deriv_y_y_mpi
 write(12)g5_deriv_y_xx_mpi
 write(12)g5_deriv_y_xy_mpi
 write(12)g5_deriv_y_yy_mpi

endif ! myrank == 0 


end subroutine sf_2
