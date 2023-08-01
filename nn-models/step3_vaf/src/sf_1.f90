subroutine sf_1( posi, xlattice, ylattice, zlattice, &
                 g2_mpi, eta2, rs2, &
                 g5_mpi, eta5, theta5, zeta5, lambda5, &
                 g2_deriv_mpi, g5_deriv_mpi, &
                 num_g2, &
                 num_g5, &
                 natom_x, &
                 natom, rc, neighbor, io_sf, num_mpi, myrank )
implicit none
include 'mpif.h'
integer ierr, num_mpi, myrank
integer i, j, k, l, m, n, ii, iii, i0
integer ig2, ig5
integer ix1, iy1, iz1, lattice_box, boxsize, boxcenter, neighbor, center
integer natom_x, natom, atom1, atom2
! natom_x : number of atom_x
! natom   : total number of atoms
double precision xlattice, ylattice, zlattice
double precision posi(  natom, 3 )
double precision posi1( natom, 3 )
double precision,allocatable :: posi2(:,:)
integer list( natom, neighbor )

integer num_g2, num_g5
double precision g2( natom_x, num_g2 ), eta2( num_g2 ), rs2( num_g2 ) 
double precision g2_mpi( natom_x, num_g2 ), g5_mpi( natom_x, num_g5 )
double precision g5( natom_x, num_g5 ), &
                 eta5( num_g5 ), theta5( num_g5 ), zeta5( num_g5 ), lambda5( num_g5 )
double precision g2_deriv( natom_x, natom, 3, num_g2 )
double precision g5_deriv( natom_x, natom, 3, num_g5 )
double precision g2_deriv_mpi( natom_x, natom, 3, num_g2 )
double precision g5_deriv_mpi( natom_x, natom, 3, num_g5 )
double precision g2_const_1, g2_deriv_1
double precision g5_const_1, g5_const_2, g5_const_3, g5_const_x, g5_const_y, g5_const_z
double precision Rij, Rik, Rij_v(3), Rik_v(3)
double precision exp_ij_ik, fc_ij_ik, theta_ijk
double precision theta_ijk_dx, theta_ijk_dy, theta_ijk_dz
double precision theta_ijk_dx1, theta_ijk_dy1, theta_ijk_dz1
double precision theta_ijk_dx2, theta_ijk_dy2, theta_ijk_dz2
double precision rc
double precision,parameter :: pi = 4.0d0*datan(1.0d0)
integer io_sf


 do i = 1, natom_x
   posi1(i,1) = posi(i,1) * xlattice
   posi1(i,2) = posi(i,2) * ylattice
   posi1(i,3) = posi(i,3) * zlattice
 enddo


!Repeating box
 if( min( xlattice, ylattice, zlattice ) >= rc )then
   lattice_box = 27
   boxsize     = 3
   boxcenter   = 2
   center      = natom_x * 13
 elseif( rc > min( xlattice, ylattice, zlattice ) &
         .and. 2.0d0 * min( xlattice, ylattice, zlattice ) >= rc )then
   lattice_box = 125
   boxsize     = 5
   boxcenter   = 3
   center      = natom_x * 62
 else
   write(*,*) "Error(sf_1): lattice_box"
   stop
 endif

 if( myrank == 0 .and. io_sf == 0 )then
   write(*,*) " "
   write(*,*) " Repeating box size: ", lattice_box
   write(*,*) " "
   io_sf = 1
 endif


 allocate( posi2( lattice_box * natom_x, 3 ) )

 m = 0
 do ix1 = 1, boxsize
   do iy1 = 1, boxsize
     do iz1 = 1, boxsize
       do i = 1, natom_x
         m = m + 1
         posi2(m,1) = posi1(i,1) + xlattice*( ix1 - boxcenter )
         posi2(m,2) = posi1(i,2) + ylattice*( iy1 - boxcenter )
         posi2(m,3) = posi1(i,3) + zlattice*( iz1 - boxcenter )
       enddo
     enddo
   enddo
 enddo


!Neighboring list
!if( mod(time,neighbor_check)==1 .or. jump == 1 )then
     list = 0
  do i = 1, natom_x
     m = 0
     n = 0
    do ix1 = 1, boxsize
      do iy1 = 1, boxsize
        do iz1 = 1, boxsize
          do j = 1, natom_x

            m = m + 1

            Rij = dsqrt(   ( posi1(i,1) - posi2(m,1) )**2 &
                         + ( posi1(i,2) - posi2(m,2) )**2 &
                         + ( posi1(i,3) - posi2(m,3) )**2 )

            if( Rij <= Rc )then !+ Rn )then

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


!  G2 & G5
   g2 = 0.0d0
   g5 = 0.0d0
!  G2_deriv & G5_deriv
   g2_deriv = 0.0d0
   g5_deriv = 0.0d0

 do i0 = 1, natom_x, num_mpi
   i = i0 + myrank
   if( i <= natom_x )then
 !do i = 1, natom_x ! Target atom

     do ii = 1, neighbor

       atom1 = list(i,ii)

       if( atom1 == 0 )cycle
       if( atom1 - center == i )cycle

       j = atom1 - natom_x*( (atom1-1)/natom_x )

       Rij = dsqrt(   ( posi1(i,1) - posi2(atom1,1) )**2 &
                    + ( posi1(i,2) - posi2(atom1,2) )**2 &
                    + ( posi1(i,3) - posi2(atom1,3) )**2 )

       if( Rij < Rc )then
         do ig2 = 1, num_g2

           g2_const_1 = dexp( -eta2(ig2)*( ( Rij - rs2(ig2) )**2 ) )
           g2(i,ig2) =   g2(i,ig2) + g2_const_1*fc(Rij)
           g2_deriv_1 = -2.0d0*eta2(ig2)*( Rij - rs2(ig2) )*fc(Rij) &
                        -0.50d0*(pi/Rc)*dsin( ( pi*Rij/Rc ) ) 

           !(i,i)
           g2_deriv(i,i,1,ig2) &
           =   g2_deriv(i,i,1,ig2)                                           &
             + g2_const_1*( ( posi1(i,1) - posi2(atom1,1) )/Rij )*g2_deriv_1

           g2_deriv(i,i,2,ig2)                                               &
           =   g2_deriv(i,i,2,ig2)                                           &
             + g2_const_1*( ( posi1(i,2) - posi2(atom1,2) )/Rij )*g2_deriv_1

           g2_deriv(i,i,3,ig2)                                               &
           =   g2_deriv(i,i,3,ig2)                                           &
             + g2_const_1*( ( posi1(i,3) - posi2(atom1,3) )/Rij )*g2_deriv_1

!          (i,j)
           g2_deriv(i,j,1,ig2)                                               &
           =   g2_deriv(i,j,1,ig2)                                           &
             - g2_const_1*( ( posi1(i,1) - posi2(atom1,1) )/Rij )*g2_deriv_1

           g2_deriv(i,j,2,ig2)                                               &
           =   g2_deriv(i,j,2,ig2)                                           &
             - g2_const_1*( ( posi1(i,2) - posi2(atom1,2) )/Rij )*g2_deriv_1

           g2_deriv(i,j,3,ig2)                                               &
           =   g2_deriv(i,j,3,ig2)                                           &
             - g2_const_1*( ( posi1(i,3) - posi2(atom1,3) )/Rij )*g2_deriv_1

         enddo

           Rij_v(1) = posi2(atom1,1) - posi1(i,1)
           Rij_v(2) = posi2(atom1,2) - posi1(i,2)
           Rij_v(3) = posi2(atom1,3) - posi1(i,3)

       endif


     do iii = 1, neighbor

       atom2 = list(i,iii)

       if( atom2 == 0 )cycle
       if( atom2-center == i )cycle
       if( atom1 == atom2 )cycle

       k = atom2 - natom_x*( (atom2-1)/natom_x )

         Rik = dsqrt(                                      &
                        ( posi1(i,1) - posi2(atom2,1) )**2 &
                      + ( posi1(i,2) - posi2(atom2,2) )**2 &
                      + ( posi1(i,3) - posi2(atom2,3) )**2 &
                    )

         if( Rij <= Rc .and. Rik <= Rc )then

           Rik_v(1) = posi2(atom2,1) - posi1(i,1)
           Rik_v(2) = posi2(atom2,2) - posi1(i,2)
           Rik_v(3) = posi2(atom2,3) - posi1(i,3)

           theta_ijk = (   Rij_v(1)*Rik_v(1) &
                         + Rij_v(2)*Rik_v(2) &
                         + Rij_v(3)*Rik_v(3) &
                       )/( Rij*Rik )

           theta_ijk_dx =  -( Rij_v(1) + Rik_v(1) )/( Rij*Rik ) &
                          + (                                   &
                             (  Rij_v(1)*Rik_v(1)               &
                              + Rij_v(2)*Rik_v(2)               &
                              + Rij_v(3)*Rik_v(3)               &
                             )*Rij_v(1)                         &
                            )/( (Rij**3)*Rik )                  &
                          + (                                   &
                             (  Rij_v(1)*Rik_v(1)               &
                              + Rij_v(2)*Rik_v(2)               &
                              + Rij_v(3)*Rik_v(3)               &
                             )*Rik_v(1)                         &
                            )/( Rij*(Rik**3) )

           theta_ijk_dy =  -( Rij_v(2) + Rik_v(2) )/( Rij*Rik ) &
                          + (                                   &
                             (  Rij_v(1)*Rik_v(1)               &
                              + Rij_v(2)*Rik_v(2)               &
                              + Rij_v(3)*Rik_v(3)               &
                             )*Rij_v(2)                         &
                            )/( (Rij**3)*Rik )                  &
                          + (                                   &
                             (  Rij_v(1)*Rik_v(1)               &
                              + Rij_v(2)*Rik_v(2)               &
                              + Rij_v(3)*Rik_v(3)               &
                             )*Rik_v(2)                         &
                            )/( Rij*(Rik**3) )

           theta_ijk_dz =  -( Rij_v(3) + Rik_v(3) )/( Rij*Rik ) &
                          + (                                   &
                             (  Rij_v(1)*Rik_v(1)               &
                              + Rij_v(2)*Rik_v(2)               &
                              + Rij_v(3)*Rik_v(3)               &
                             )*Rij_v(3)                         &
                            )/( (Rij**3)*Rik )                  &
                          + (                                   &
                             (  Rij_v(1)*Rik_v(1)               &
                              + Rij_v(2)*Rik_v(2)               &
                              + Rij_v(3)*Rik_v(3)               &
                              )*Rik_v(3)                        &
                            )/( Rij*(Rik**3) )

           theta_ijk_dx1 =   ( ( posi2(atom2,1) - posi1(i,1) )/(Rij*Rik) )     &
                           - (  Rij_v(1)*Rik_v(1)                              &
                              + Rij_v(2)*Rik_v(2)                              &
                              + Rij_v(3)*Rik_v(3)                              &
                             )*( posi2(atom1,1) - posi1(i,1) )/( (Rij**3)*Rik )

           theta_ijk_dy1 =   ( ( posi2(atom2,2) - posi1(i,2) )/(Rij*Rik) )     &
                           - (  Rij_v(1)*Rik_v(1)                              &
                              + Rij_v(2)*Rik_v(2)                              &
                              + Rij_v(3)*Rik_v(3)                              &
                             )*( posi2(atom1,2) - posi1(i,2) )/( (Rij**3)*Rik )

           theta_ijk_dz1 =   ( ( posi2(atom2,3) - posi1(i,3) )/(Rij*Rik) )     &
                           - (  Rij_v(1)*Rik_v(1)                              &
                              + Rij_v(2)*Rik_v(2)                              &
                              + Rij_v(3)*Rik_v(3)                              &
                             )*( posi2(atom1,3) - posi1(i,3) )/( (Rij**3)*Rik )

           theta_ijk_dx2 =   ( ( posi2(atom1,1) - posi1(i,1) )/(Rij*Rik) )     &
                           - (  Rij_v(1)*Rik_v(1)                              &
                              + Rij_v(2)*Rik_v(2)                              &
                              + Rij_v(3)*Rik_v(3)                              &
                             )*( posi2(atom2,1) - posi1(i,1) )/( (Rik**3)*Rij )

           theta_ijk_dy2 =   ( ( posi2(atom1,2) - posi1(i,2) )/(Rij*Rik) )     &
                           - (  Rij_v(1)*Rik_v(1)                              &
                              + Rij_v(2)*Rik_v(2)                              &
                              + Rij_v(3)*Rik_v(3)                              &
                             )*( posi2(atom2,2) - posi1(i,2) )/( (Rik**3)*Rij )

           theta_ijk_dz2 =   ( ( posi2(atom1,3) - posi1(i,3) )/(Rij*Rik) )     &
                           - (  Rij_v(1)*Rik_v(1)                              &
                              + Rij_v(2)*Rik_v(2)                              &
                              + Rij_v(3)*Rik_v(3)                              &
                             )*( posi2(atom2,3) - posi1(i,3) )/( (Rik**3)*Rij )

           g5_const_x =   fc_deriv( Rij, posi1(i,1), posi2(atom1,1) )*fc(Rik) &
                        + fc(Rij)*fc_deriv( Rik, posi1(i,1), posi2(atom2,1) )

           g5_const_y =   fc_deriv( Rij, posi1(i,2), posi2(atom1,2) )*fc(Rik) &
                        + fc(Rij)*fc_deriv( Rik, posi1(i,2), posi2(atom2,2) )

           g5_const_z =   fc_deriv( Rij, posi1(i,3), posi2(atom1,3) )*fc(Rik) &
                        + fc(Rij)*fc_deriv( Rik, posi1(i,3), posi2(atom2,3) )

           fc_ij_ik  = fc(Rij)*fc(Rik)

             do ig5=1,num_g5

               exp_ij_ik = dexp( -eta5(ig5)*( Rij**2 + Rik**2 ) )

               g5_const_1 = ( 1.0d0 + lambda5(ig5)*(theta_ijk) )**(zeta5(ig5))
               g5_const_2 = ( 1.0d0 + lambda5(ig5)*(theta_ijk) )**(zeta5(ig5)-1.0d0)

               g5(i,ig5) =   g5(i,ig5)                     &
                           + g5_const_1*exp_ij_ik*fc_ij_ik

!              (i,i)
               g5_deriv(i,i,1,ig5)                        &
               =   g5_deriv(i,i,1,ig5)                    &
                 + zeta5(ig5)*g5_const_2*lambda5(ig5)     &
                   *theta_ijk_dx*exp_ij_ik*fc_ij_ik       &
                 + g5_const_1                             &
                   *( -2.0d0*eta5(ig5)*exp_ij_ik          &
                      *(                                  &
                          ( posi1(i,1) - posi2(atom1,1) ) &
                        + ( posi1(i,1) - posi2(atom2,1) ) &
                       )                                  &
                      *fc_ij_ik                           &
                      +                                   &
                      exp_ij_ik*g5_const_x                &
                    )

               g5_deriv(i,i,2,ig5)                        &
               =   g5_deriv(i,i,2,ig5)                    &
                 + zeta5(ig5)*g5_const_2*lambda5(ig5)     &
                   *theta_ijk_dy*exp_ij_ik*fc_ij_ik       &
                 + g5_const_1                             &
                   *( -2.0d0*eta5(ig5)*exp_ij_ik          &
                      *(                                  &
                          ( posi1(i,2) - posi2(atom1,2) ) &
                        + ( posi1(i,2) - posi2(atom2,2) ) &
                       )                                  &
                      *fc_ij_ik                           &
                      +                                   &
                      exp_ij_ik*g5_const_y                &
                    )

               g5_deriv(i,i,3,ig5)                        &
               =   g5_deriv(i,i,3,ig5)                    &
                 + zeta5(ig5)*g5_const_2*lambda5(ig5)     &
                   *theta_ijk_dz*exp_ij_ik*fc_ij_ik       &
                 + g5_const_1                             &
                   *( -2.0d0*eta5(ig5)*exp_ij_ik          &
                      *(                                  &
                          ( posi1(i,3) - posi2(atom1,3) ) &
                        + ( posi1(i,3) - posi2(atom2,3) ) &
                       )                                  &
                      *fc_ij_ik                           &
                      +                                   &
                      exp_ij_ik*g5_const_z                &
                    )

!              (i,j)
               g5_deriv(i,j,1,ig5)                                         &
               =   g5_deriv(i,j,1,ig5)                                     &
                 + zeta5(ig5)*g5_const_2*lambda5(ig5)                      &
                   *theta_ijk_dx1*exp_ij_ik*fc_ij_ik                       &
                 + g5_const_1                                              &
                   *( -2.0d0*eta5(ig5)*exp_ij_ik                           &
                      *( posi2(atom1,1) - posi1(i,1) )*fc_ij_ik            &
                      -                                                    &
                      exp_ij_ik                                            &
                      *fc_deriv( Rij, posi1(i,1), posi2(atom1,1) )*fc(Rik) &
                    )

               g5_deriv(i,j,2,ig5)                                         &
               =   g5_deriv(i,j,2,ig5)                                     &
                 + zeta5(ig5)*g5_const_2*lambda5(ig5)                      &
                   *theta_ijk_dy1*exp_ij_ik*fc_ij_ik                       &
                 + g5_const_1                                              &
                   *( -2.0d0*eta5(ig5)*exp_ij_ik                           &
                      *( posi2(atom1,2) - posi1(i,2) )*fc_ij_ik            &
                      -                                                    &
                      exp_ij_ik                                            &
                      *fc_deriv( Rij, posi1(i,2), posi2(atom1,2) )*fc(Rik) &
                    )

               g5_deriv(i,j,3,ig5)                                         &
               =   g5_deriv(i,j,3,ig5)                                     &
                 + zeta5(ig5)*g5_const_2*lambda5(ig5)                      &
                   *theta_ijk_dz1*exp_ij_ik*fc_ij_ik                       &
                 + g5_const_1                                              &
                   *( -2.0d0*eta5(ig5)*exp_ij_ik                           &
                      *( posi2(atom1,3) - posi1(i,3) )*fc_ij_ik            &
                      -                                                    &
                      exp_ij_ik                                            &
                      *fc_deriv( Rij, posi1(i,3), posi2(atom1,3) )*fc(Rik) &
                    )

               g5_deriv(i,k,1,ig5)                                         &
               =   g5_deriv(i,k,1,ig5)                                     &
                 + zeta5(ig5)*g5_const_2*lambda5(ig5)                      &
                   *theta_ijk_dx2*exp_ij_ik*fc_ij_ik                       &
                 + g5_const_1                                              &
                   *( -2.0d0*eta5(ig5)*exp_ij_ik                           &
                      *( posi2(atom2,1) - posi1(i,1) )*fc_ij_ik            &
                      -                                                    &
                      exp_ij_ik                                            &
                      *fc_deriv( Rik, posi1(i,1), posi2(atom2,1) )*fc(Rij) &
                    )

               g5_deriv(i,k,2,ig5)                                         &
               =   g5_deriv(i,k,2,ig5)                                     &
                 + zeta5(ig5)*g5_const_2*lambda5(ig5)                      &
                   *theta_ijk_dy2*exp_ij_ik*fc_ij_ik                       &
                 + g5_const_1                                              &
                   *( -2.0d0*eta5(ig5)*exp_ij_ik                           &
                      *( posi2(atom2,2) - posi1(i,2) )*fc_ij_ik            &
                      -                                                    &
                      exp_ij_ik                                            &
                      *fc_deriv( Rik, posi1(i,2), posi2(atom2,2) )*fc(Rij) &
                    )

               g5_deriv(i,k,3,ig5)                                         &
               =   g5_deriv(i,k,3,ig5)                                     &
                 + zeta5(ig5)*g5_const_2*lambda5(ig5)                      &
                   *theta_ijk_dz2*exp_ij_ik*fc_ij_ik                       &
                 + g5_const_1                                              &
                   *( -2.0d0*eta5(ig5)*exp_ij_ik                           &
                      *( posi2(atom2,3) - posi1(i,3) )*fc_ij_ik            &
                      -                                                    &
                      exp_ij_ik                                            &
                      *fc_deriv( Rik, posi1(i,3), posi2(atom2,3) )*fc(Rij) &
                    )

             enddo

         endif ! Rij <= Rc .and. Rik <= Rc

     enddo ! iii
   enddo ! ii

   endif ! i <= natom_x
 enddo ! i0

 call mpi_barrier( mpi_comm_world, ierr )


!G5 coefficient
 do i0 = 1, natom, num_mpi
   i = i0 + myrank
   if( i <= natom )then

   do ig5 = 1, num_g5

       g5_const_3 = 2.0d0**( 1.0d0 - zeta5(ig5) )

       g5(i,ig5) = g5(i,ig5)*g5_const_3

     do j = 1, natom_x

       g5_deriv(i,j,1,ig5) = g5_deriv(i,j,1,ig5)*g5_const_3
       g5_deriv(i,j,2,ig5) = g5_deriv(i,j,2,ig5)*g5_const_3
       g5_deriv(i,j,3,ig5) = g5_deriv(i,j,3,ig5)*g5_const_3

     enddo
   enddo

 endif

 enddo

 call mpi_barrier( mpi_comm_world, ierr )

 call mpi_reduce( g2, g2_mpi, natom_x*num_g2, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5, g5_mpi, natom_x*num_g5, mpi_double_precision, &
                  mpi_sum, 0, mpi_comm_world, ierr )


 call mpi_reduce( g2_deriv, g2_deriv_mpi, natom_x*natom*3*num_g2, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_reduce( g5_deriv, g5_deriv_mpi, natom_x*natom*3*num_g5, &
                  mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr )
 call mpi_barrier( mpi_comm_world, ierr )

 deallocate( posi2 )


!--------
!Check SF
!--------

!X
do j = 1, num_g2
  do i = 1, natom_x
    if( isnan( g2(i,j) ) ) stop 'NaN g2_x_x'
  enddo
enddo

do j = 1, num_g5
  do i = 1, natom_x
    if( isnan( g5(i,j) ) ) stop 'NaN g5_x_xx'
  enddo
enddo


!Deriv X
do l = 1, num_g2
do k = 1, 3
do j = 1, natom
do i = 1, natom_x
  if( isnan( g2_deriv(i,j,k,l) ) ) stop 'NaN g2_deriv_x_x'
enddo
enddo
enddo
enddo

do l = 1, num_g5
do k = 1, 3
do j = 1, natom
do i = 1, natom_x
  if( isnan( g5_deriv(i,j,k,l) ) ) stop 'NaN g5_deriv_x_xx'
enddo
enddo
enddo
enddo


!--------
!Write sf
!--------
if( myrank == 0 )then
  write(11)g2_mpi
  write(11)g5_mpi
endif ! myrank == 0
!--------------
!Write sf_deriv
!--------------
if( myrank == 0 )then
  write(12)g2_deriv_mpi
  write(12)g5_deriv_mpi
endif ! myrank == 0 



contains
!***
double precision function fc( r_ij )
!use parameters
implicit none
double precision,parameter :: pi = 4.0d0*datan(1.0d0)
double precision r_ij
  fc = 0.50d0*( dcos( ( pi*r_ij ) / Rc ) + 1.0d0 )
end function fc
!***
double precision function fc_deriv( r_ij, x, x_j )
!use parameters
implicit none
double precision,parameter :: pi = 4.0d0*datan(1.0d0)
double precision r_ij, x, x_j
  fc_deriv = -(0.50d0*pi*( x - x_j )*dsin( ( pi*r_ij )/Rc ) )/( Rc*r_ij )
end function fc_deriv


end subroutine sf_1
