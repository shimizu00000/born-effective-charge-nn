subroutine sf_3_g5( g5, g5_deriv, eta5, R_s, zeta5, num_g5, &
                    posi1_x,   posi1_y,   posi1_z, &
                    posi2_j_x, posi2_j_y, posi2_j_z, &
                    posi2_k_x, posi2_k_y, posi2_k_z, &
                    natom, natom_x, natom_y, natom_z, &
                    Rij, Rik, i, j, k, rc, n )
implicit none
integer i, j, k, ig5, num_g5, n, ii
integer natom
integer natom_x, natom_y, natom_z
double precision g5( n, num_g5 ), g5_deriv( n, natom, 3, num_g5 )
double precision eta5( num_g5 ), R_s( num_g5 ), zeta5( num_g5 )
double precision Rij, Rik, exp_ij_ik, fc_ij_ik
double precision g5_const_x, g5_const_y, g5_const_z, g5_const_1, g5_const_2, g5_const_3
double precision posi1_x, posi1_y, posi1_z
double precision posi2_j_x, posi2_j_y, posi2_j_z
double precision posi2_k_x, posi2_k_y, posi2_k_z
double precision cos_theta_ijk
double precision cos_theta_ijk_dx,  cos_theta_ijk_dy,  cos_theta_ijk_dz
double precision cos_theta_ijk_dx1, cos_theta_ijk_dy1, cos_theta_ijk_dz1
double precision cos_theta_ijk_dx2, cos_theta_ijk_dy2, cos_theta_ijk_dz2
double precision Rij_v(3), Rik_v(3)
double precision rc
double precision,parameter :: pi = 4.0d0 * datan(1.0d0)


 if( i <= natom_x )then
   ii = i
 elseif( natom_x < i .and. i <= natom_x + natom_y )then
   ii = i - natom_x
 elseif( natom_x + natom_y < i .and. i <= natom_x + natom_y + natom_z )then
   ii = i - natom_x - natom_y
 else
   write(*,*) "Error(sf_3_g5): ii"
   stop
 endif


 Rij_v(1) = posi2_j_x - posi1_x !posi2(atom1,1) - posi1(i,1)
 Rij_v(2) = posi2_j_y - posi1_y !posi2(atom1,2) - posi1(i,2)
 Rij_v(3) = posi2_j_z - posi1_z !posi2(atom1,3) - posi1(i,3)

 Rik_v(1) = posi2_k_x - posi1_x !posi2(atom2,1) - posi1(i,1)
 Rik_v(2) = posi2_k_y - posi1_y !posi2(atom2,2) - posi1(i,2)
 Rik_v(3) = posi2_k_z - posi1_z !posi2(atom2,3) - posi1(i,3)

 cos_theta_ijk = (   Rij_v(1) * Rik_v(1) &
               + Rij_v(2) * Rik_v(2) &
               + Rij_v(3) * Rik_v(3) &
             )/( Rij * Rik )
             
 cos_theta_ijk_dx =  -( Rij_v(1) + Rik_v(1) )/( Rij * Rik ) &
                + ( &
                   (  Rij_v(1) * Rik_v(1) &
                    + Rij_v(2) * Rik_v(2) &
                    + Rij_v(3) * Rik_v(3) &
                   ) * Rij_v(1) &
                  )/( ( Rij**3 ) * Rik ) &
                + ( &
                   (  Rij_v(1) * Rik_v(1) &
                    + Rij_v(2) * Rik_v(2) &
                    + Rij_v(3) * Rik_v(3) &
                   ) * Rik_v(1) &
                  )/( Rij * ( Rik**3 ) )

 cos_theta_ijk_dy =  -( Rij_v(2) + Rik_v(2) )/( Rij * Rik ) &
                + ( &
                   (  Rij_v(1) * Rik_v(1) &
                    + Rij_v(2) * Rik_v(2) &
                    + Rij_v(3) * Rik_v(3) &
                   ) * Rij_v(2) &
                  )/( ( Rij**3 ) * Rik ) &
                + ( &
                   (  Rij_v(1) * Rik_v(1) &
                    + Rij_v(2) * Rik_v(2) &
                    + Rij_v(3) * Rik_v(3) &
                   ) * Rik_v(2) &
                  )/( Rij * ( Rik**3 ) )

 cos_theta_ijk_dz =  -( Rij_v(3) + Rik_v(3) )/( Rij * Rik ) &
                + ( &
                   (  Rij_v(1) * Rik_v(1) &
                    + Rij_v(2) * Rik_v(2) &
                    + Rij_v(3) * Rik_v(3) &
                   ) * Rij_v(3) &
                  )/( ( Rij**3 ) * Rik ) &
                + ( &
                   (  Rij_v(1) * Rik_v(1) &
                    + Rij_v(2) * Rik_v(2) &
                    + Rij_v(3) * Rik_v(3) &
                    ) * Rik_v(3) &
                  )/( Rij * ( Rik**3 ) )

 !theta_ijk_dx1 =   ( ( posi2(atom2,1) - posi1(i,1) )/(Rij*Rik) ) &
 cos_theta_ijk_dx1 =   ( ( posi2_k_x - posi1_x )/( Rij * Rik ) ) &
                 - (  Rij_v(1) * Rik_v(1) &
                    + Rij_v(2) * Rik_v(2) &
                    + Rij_v(3) * Rik_v(3) &
                   ) * ( posi2_j_x - posi1_x )/( ( Rij**3 ) * Rik )
                   !)*( posi2(atom1,1) - posi1(i,1) )/( (Rij**3)*Rik )
 !theta_ijk_dy1 =   ( ( posi2(atom2,2) - posi1(i,2) )/(Rij*Rik) ) &
 cos_theta_ijk_dy1 =   ( ( posi2_k_y - posi1_y )/( Rij * Rik ) ) &
                 - (  Rij_v(1) * Rik_v(1) &
                    + Rij_v(2) * Rik_v(2) &
                    + Rij_v(3) * Rik_v(3) &
                   ) * ( posi2_j_y - posi1_y )/( ( Rij**3 ) * Rik )
                   !)*( posi2(atom1,2) - posi1(i,2) )/( (Rij**3)*Rik )
 !theta_ijk_dz1 =   ( ( posi2(atom2,3) - posi1(i,3) )/(Rij*Rik) ) &
 cos_theta_ijk_dz1 =   ( ( posi2_k_z - posi1_z )/( Rij * Rik ) ) &
                 - (  Rij_v(1) * Rik_v(1) &
                    + Rij_v(2) * Rik_v(2) &
                    + Rij_v(3) * Rik_v(3) &
                   ) * ( posi2_j_z - posi1_z )/( ( Rij**3 ) * Rik )
                   !)*( posi2(atom1,3) - posi1(i,3) )/( (Rij**3)*Rik )

 !theta_ijk_dx2 =   ( ( posi2(atom1,1) - posi1(i,1) )/(Rij*Rik) ) &
 cos_theta_ijk_dx2 =   ( ( posi2_j_x - posi1_x )/( Rij * Rik ) ) &
                 - (  Rij_v(1) * Rik_v(1) &
                    + Rij_v(2) * Rik_v(2) &
                    + Rij_v(3) * Rik_v(3) &
                   ) * ( posi2_k_x - posi1_x )/( ( Rik**3 ) * Rij )
                   !)*( posi2(atom2,1) - posi1(i,1) )/( (Rik**3)*Rij )
 !theta_ijk_dy2 =   ( ( posi2(atom1,2) - posi1(i,2) )/(Rij*Rik) ) &
 cos_theta_ijk_dy2 =   ( ( posi2_j_y - posi1_y )/( Rij * Rik ) ) &
                 - (  Rij_v(1) * Rik_v(1) &
                    + Rij_v(2) * Rik_v(2) &
                    + Rij_v(3) * Rik_v(3) &
                   ) * ( posi2_k_y - posi1_y )/( ( Rik**3 ) * Rij )
                   !)*( posi2(atom2,2) - posi1(i,2) )/( (Rik**3)*Rij )
 !theta_ijk_dz2 =   ( ( posi2(atom1,3) - posi1(i,3) )/(Rij*Rik) ) &
 cos_theta_ijk_dz2 =   ( ( posi2_j_z - posi1_z )/( Rij * Rik ) ) &
                 - (  Rij_v(1) * Rik_v(1) &
                    + Rij_v(2) * Rik_v(2) &
                    + Rij_v(3) * Rik_v(3) &
                   )*( posi2_k_z - posi1_z )/( ( Rik**3 ) * Rij )
                   !)*( posi2(atom2,3) - posi1(i,3) )/( (Rik**3)*Rij )

 !g5_const_x =   fc_deriv( Rij, posi1(i,1), posi2(atom1,1) )*fc(Rik) &
 g5_const_x =   fc_deriv5( Rij, posi1_x, posi2_j_x ) * fc5(Rik) &
              + fc5(Rij) * fc_deriv5( Rik, posi1_x, posi2_k_x )
 g5_const_y =   fc_deriv5( Rij, posi1_y, posi2_j_y ) * fc5(Rik) &
              + fc5(Rij) * fc_deriv5( Rik, posi1_y, posi2_k_y )
 g5_const_z =   fc_deriv5( Rij, posi1_z, posi2_j_z ) * fc5(Rik) &
              + fc5(Rij) * fc_deriv5( Rik, posi1_z, posi2_k_z )

 fc_ij_ik  = fc5(Rij) * fc5(Rik)


 do ig5 = 1, num_g5

   exp_ij_ik = dexp( -eta5(ig5) * (( Rij + Rik ) * 0.5d0 - R_s(ig5) )**2 )

   g5_const_1 = ( 1.0d0 + cos_theta_ijk )**( zeta5(ig5) )
   g5_const_2 = ( 1.0d0 + cos_theta_ijk )**( zeta5(ig5) - 1.0d0 )
   g5_const_3 = ( posi1_z - posi2_j_z ) + ( posi1_z - posi2_k_z )
   g5(ii,ig5) = g5(ii,ig5) + g5_const_1 * g5_const_3 * exp_ij_ik * fc_ij_ik

!  (i,i)
   g5_deriv(ii,i,1,ig5) &
   =   g5_deriv(ii,i,1,ig5) &
     + zeta5(ig5) * g5_const_2 &
       * cos_theta_ijk_dx * exp_ij_ik * fc_ij_ik * g5_const_3 &
     + g5_const_1 * g5_const_3 &
       *( -eta5(ig5) * exp_ij_ik &
          *( &
              !( posi1(i,1) - posi2(atom1,1) ) &
              ( posi1_x - posi2_j_x ) * ( Rij + Rik - 2.0d0 * R_s(ig5) ) / Rij &
            + ( posi1_x - posi2_k_x ) * ( Rij + Rik - 2.0d0 * R_s(ig5) ) / Rik &
            !+( posi1(i,1) - posi2(atom2,1) ) &
           ) &
          * fc_ij_ik &
          + &
          exp_ij_ik * g5_const_x &
        )
   if( isnan( g5_deriv(ii,i,1,ig5) ) ) then
      write(*,*) zeta5(ig5), g5_const_2, cos_theta_ijk_dx, g5_const_3, cos_theta_ijk, exp_ij_ik, fc_ij_ik, posi1_x,   posi1_y,   posi1_z, &
                    posi2_j_x, posi2_j_y, posi2_j_z, &
                    posi2_k_x, posi2_k_y, posi2_k_z
      stop 'g5_deriv_ii_i_1_ig5 is NaN'
   end if

   g5_deriv(ii,i,2,ig5) &
   =   g5_deriv(ii,i,2,ig5) &
     + zeta5(ig5) * g5_const_2 &
       * cos_theta_ijk_dy * exp_ij_ik * fc_ij_ik * g5_const_3 &
     + g5_const_1 * g5_const_3 &
       *( -eta5(ig5) * exp_ij_ik &
          *( &
              !( posi1(i,2) - posi2(atom1,2) ) &
              ( posi1_y - posi2_j_y ) * (Rij + Rik - 2.0d0 * R_s(ig5) ) / Rij &
            + ( posi1_y - posi2_k_y ) * (Rij + Rik - 2.0d0 * R_s(ig5) ) / Rik &
            !+ ( posi1(i,2) - posi2(atom2,2) ) &
           ) &
          * fc_ij_ik &
          + &
          exp_ij_ik * g5_const_y &
        )
   g5_deriv(ii,i,3,ig5) &
   =   g5_deriv(ii,i,3,ig5) &
     + zeta5(ig5) * g5_const_2 &
       * cos_theta_ijk_dz * exp_ij_ik * fc_ij_ik * g5_const_3&
     + g5_const_1 * g5_const_3 &
       *( -eta5(ig5) * exp_ij_ik &
          *( &
              !( posi1(i,3) - posi2(atom1,3) ) &
              ( posi1_z - posi2_j_z ) * (Rij + Rik - 2.0d0 * R_s(ig5) ) / Rij &
            + ( posi1_z - posi2_k_z ) * (Rij + Rik - 2.0d0 * R_s(ig5) ) / Rik &
            !+ ( posi1(i,3) - posi2(atom2,3) ) &
           ) &
          * fc_ij_ik &
          + &
          exp_ij_ik * g5_const_z &
        ) &
     + 2.0d0 * g5_const_1 * exp_ij_ik * fc_ij_ik


!(i,j)
   g5_deriv(ii,j,1,ig5) &
   =   g5_deriv(ii,j,1,ig5) &
     + zeta5(ig5) * g5_const_2 &
       * cos_theta_ijk_dx1 * exp_ij_ik * fc_ij_ik * g5_const_3 &
     + g5_const_1 * g5_const_3 &
       *( -eta5(ig5) * exp_ij_ik &
          !*( posi2(atom1,1) - posi1(i,1) )*fc_ij_ik &
          *( posi2_j_x - posi1_x ) * (Rij + Rik - 2.0d0 * R_s(ig5) ) / Rij * fc_ij_ik &
          - &
          exp_ij_ik &
          * fc_deriv5( Rij, posi1_x, posi2_j_x ) * fc5(Rik) &
          !*fc_deriv( Rij, posi1(i,1), posi2(atom1,1) )*fc(Rik) &
        )
   g5_deriv(ii,j,2,ig5) &
   =   g5_deriv(ii,j,2,ig5) &
     + zeta5(ig5) * g5_const_2 &
       * cos_theta_ijk_dy1 * exp_ij_ik * fc_ij_ik * g5_const_3 &
     + g5_const_1 * g5_const_3 &
       *( -eta5(ig5) * exp_ij_ik &
          !*( posi2(atom1,2) - posi1(i,2) )*fc_ij_ik &
          *( posi2_j_y - posi1_y ) * (Rij + Rik - 2.0d0 * R_s(ig5) ) / Rij * fc_ij_ik &
          - &
          exp_ij_ik &
          * fc_deriv5( Rij, posi1_y, posi2_j_y ) * fc5(Rik) &
          !*fc_deriv( Rij, posi1(i,2), posi2(atom1,2) )*fc(Rik) &
        )
   g5_deriv(ii,j,3,ig5) &
   =   g5_deriv(ii,j,3,ig5) &
     + zeta5(ig5) * g5_const_2 &
       * cos_theta_ijk_dz1 * exp_ij_ik * fc_ij_ik * g5_const_3 &
     + g5_const_1 * g5_const_3 &
       *( -eta5(ig5) * exp_ij_ik &
          !*( posi2(atom1,3) - posi1(i,3) )*fc_ij_ik &
          *( posi2_j_z - posi1_z ) * (Rij + Rik - 2.0d0 * R_s(ig5) ) / Rij * fc_ij_ik &
          - &
          exp_ij_ik &
          * fc_deriv5( Rij, posi1_z, posi2_j_z ) * fc5(Rik) &
          !*fc_deriv( Rij, posi1(i,3), posi2(atom1,3) )*fc(Rik) &
        ) &
     - g5_const_1 * exp_ij_ik * fc_ij_ik
!(i,k)
   g5_deriv(ii,k,1,ig5) &
   =   g5_deriv(ii,k,1,ig5) &
     + zeta5(ig5) * g5_const_2 &
       * cos_theta_ijk_dx2 * exp_ij_ik * fc_ij_ik * g5_const_3 &
     + g5_const_1 * g5_const_3 &
       *( -eta5(ig5) * exp_ij_ik &
          !*( posi2(atom2,1) - posi1(i,1) )*fc_ij_ik &
          * ( posi2_k_x - posi1_x ) * (Rij + Rik - 2.0d0 * R_s(ig5) ) / Rik * fc_ij_ik &
          - &
          exp_ij_ik &
          * fc_deriv5( Rik, posi1_x, posi2_k_x ) * fc5(Rij) &
        )
   g5_deriv(ii,k,2,ig5) &
   =   g5_deriv(ii,k,2,ig5) &
     + zeta5(ig5) * g5_const_2 &
       * cos_theta_ijk_dy2 * exp_ij_ik * fc_ij_ik * g5_const_3 &
     + g5_const_1 * g5_const_3 &
       *( -eta5(ig5) * exp_ij_ik &
          !*( posi2(atom2,2) - posi1(i,2) )*fc_ij_ik &
          *( posi2_k_y - posi1_y ) * (Rij + Rik - 2.0d0 * R_s(ig5) ) / Rik * fc_ij_ik &
          - &
          exp_ij_ik &
          * fc_deriv5( Rik, posi1_y, posi2_k_y ) * fc5(Rij) &
          !*fc_deriv( Rik, posi1(i,2), posi2(atom2,2) )*fc(Rij) &
        )
   g5_deriv(ii,k,3,ig5) &
   =   g5_deriv(ii,k,3,ig5) &
     + zeta5(ig5) * g5_const_2 &
       * cos_theta_ijk_dz2 * exp_ij_ik * fc_ij_ik * g5_const_3 &
     + g5_const_1 * g5_const_3 &
       *( -eta5(ig5) * exp_ij_ik &
          !*( posi2(atom2,3) - posi1(i,3) )*fc_ij_ik &
          *( posi2_k_z - posi1_z ) * (Rij + Rik - 2.0d0 * R_s(ig5) ) / Rik * fc_ij_ik &
          - &
          exp_ij_ik &
          * fc_deriv5( Rik, posi1_z, posi2_k_z ) * fc5(Rij) &
          !*fc_deriv( Rik, posi1(i,3), posi2(atom2,3) )*fc(Rij) &
        ) &
     - g5_const_1 * exp_ij_ik * fc_ij_ik

enddo ! ig5


!*******
contains
!*******

double precision function fc5(r_ij)
!use parameters
implicit none
double precision r_ij

  fc5 = 0.50d0 * ( dcos( ( pi * r_ij ) / Rc ) + 1.0d0 )

end function fc5


double precision function fc_deriv5(r_ij,x,x_j)
!use parameters
implicit none
double precision r_ij, x, x_j

  fc_deriv5 = -( 0.50d0 * pi * ( x - x_j ) * dsin( ( pi * r_ij ) / Rc ) ) / ( Rc * r_ij )

end function fc_deriv5

!************
!end contains
!************


end subroutine sf_3_g5
