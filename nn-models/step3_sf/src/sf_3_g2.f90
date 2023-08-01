subroutine sf_3_g2( g2, g2_deriv, eta2, rs2, num_g2, &
                    posi1_x, posi1_y, posi1_z, &
                    posi2_x, posi2_y, posi2_z, &
                    natom, natom_x, natom_y, natom_z, &
                    Rij, i, j, rc, n )
implicit none
integer i, j, ig2, num_g2, n, ii
integer natom
integer natom_x, natom_y, natom_z
double precision g2( n, num_g2 ), g2_deriv( n, natom, 3, num_g2 )
double precision eta2( num_g2 ), rs2( num_g2 )
double precision Rij, g2_const_1, g2_deriv_1
double precision posi1_x, posi1_y, posi1_z
double precision posi2_x, posi2_y, posi2_z
double precision rc
double precision,parameter :: pi = 4.0d0 * datan(1.0d0)


 if( i <= natom_x )then
   ii = i
 elseif( natom_x < i .and. i <= natom_x + natom_y )then
   ii = i - natom_x
 elseif( natom_x + natom_y < i .and. i <= natom_x + natom_y + natom_z )then
   ii = i - natom_x - natom_y
 else
   write(*,*) "Error(sf_3_g2) ii"
   stop
 endif


 do ig2 = 1, num_g2

   g2_const_1 = dexp( -eta2(ig2) * ( ( Rij - rs2(ig2) )**2 ) )
   g2(ii,ig2) = g2(ii,ig2) + g2_const_1*fc(Rij)
   g2_deriv_1 = -2.0d0 * eta2(ig2) * ( Rij - rs2(ig2) ) * fc(Rij) &
                -0.50d0 * ( pi / Rc ) * dsin( ( pi * Rij / Rc ) )

!  (i,i)
   g2_deriv(ii,i,1,ig2)                                               &
   =   g2_deriv(ii,i,1,ig2)                                           &
     + g2_const_1 * ( ( posi1_x - posi2_x ) / Rij ) * g2_deriv_1
   g2_deriv(ii,i,2,ig2)                                               &
   =   g2_deriv(ii,i,2,ig2)                                           &
     + g2_const_1 * ( ( posi1_y - posi2_y ) / Rij ) * g2_deriv_1
   g2_deriv(ii,i,3,ig2)                                               &
   =   g2_deriv(ii,i,3,ig2)                                           &
     + g2_const_1 * ( ( posi1_z - posi2_z ) / Rij ) * g2_deriv_1

!  (i,j)
   g2_deriv(ii,j,1,ig2)                                               &
   =   g2_deriv(ii,j,1,ig2)                                           &
     - g2_const_1 * ( ( posi1_x - posi2_x ) / Rij ) * g2_deriv_1
   g2_deriv(ii,j,2,ig2)                                               &
   =   g2_deriv(ii,j,2,ig2)                                           &
     - g2_const_1 * ( ( posi1_y - posi2_y ) / Rij ) * g2_deriv_1
   g2_deriv(ii,j,3,ig2)                                               &
   =   g2_deriv(ii,j,3,ig2)                                           &
     - g2_const_1 * ( ( posi1_z - posi2_z ) / Rij ) * g2_deriv_1

 enddo ! do ig2=1,num_g2_A_A


!*******
contains
!*******

double precision function fc(r_ij)
!use parameters
implicit none
double precision r_ij

  fc = 0.50d0 * ( dcos( ( pi * r_ij ) / Rc ) + 1.0d0 )

end function fc

!************
!end contains
!************


end subroutine sf_3_g2
