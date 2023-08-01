subroutine f_ef_deriv3abc( &
           g2_x_x,  g2_deriv_x_x,  num_g2_x_x, &
           g2_x_y,  g2_deriv_x_y,  num_g2_x_y, &
           g2_x_z,  g2_deriv_x_z,  num_g2_x_z, &
           g5_x_xx, g5_deriv_x_xx, num_g5_x_xx, &
           g5_x_xy, g5_deriv_x_xy, num_g5_x_xy, &
           g5_x_xz, g5_deriv_x_xz, num_g5_x_xz, &
           g5_x_yy, g5_deriv_x_yy, num_g5_x_yy, &
           g5_x_yz, g5_deriv_x_yz, num_g5_x_yz, &
           g5_x_zz, g5_deriv_x_zz, num_g5_x_zz, &
           num_g_x, &
           network_x, nlayer_x, weight_x, nweight_x, hidden_x, nhidden_x, &
           beta, g_e_x, g_f_x, &
           natom, natom_x, energy, energy0, force, force0 )
implicit none
integer i, j, k, m
integer l1, l2, l3, l4
integer h1, h2
integer w1, w2, w3
integer h2_0
integer w2_0, w3_0
integer ll1, ll2, ll3, ll4
integer hh1, hh2
integer ww1, ww2, ww3
integer www2
!integer n ! nparam
integer natom_x, natom
integer num_g2_x_x,  num_g2_x_y,  num_g2_x_z, &
        num_g5_x_xx, num_g5_x_xy, num_g5_x_xz, & 
        num_g5_x_yy, num_g5_x_yz, &
        num_g5_x_zz
integer num_g_x(9)
double precision &
  g2_x_x(  natom_x, num_g2_x_x ), &
  g2_x_y(  natom_x, num_g2_x_y ), &
  g2_x_z(  natom_x, num_g2_x_z ), &
  g5_x_xx( natom_x, num_g5_x_xx ), &
  g5_x_xy( natom_x, num_g5_x_xy ), &
  g5_x_xz( natom_x, num_g5_x_xz ), &
  g5_x_yy( natom_x, num_g5_x_yy ), &
  g5_x_yz( natom_x, num_g5_x_yz ), &
  g5_x_zz( natom_x, num_g5_x_zz )
double precision &
  g2_deriv_x_x(  natom_x, natom, 3, num_g2_x_x ), &
  g2_deriv_x_y(  natom_x, natom, 3, num_g2_x_y ), &
  g2_deriv_x_z(  natom_x, natom, 3, num_g2_x_z ), &
  g5_deriv_x_xx( natom_x, natom, 3, num_g5_x_xx ), &
  g5_deriv_x_xy( natom_x, natom, 3, num_g5_x_xy ), &
  g5_deriv_x_xz( natom_x, natom, 3, num_g5_x_xz ), &
  g5_deriv_x_yy( natom_x, natom, 3, num_g5_x_yy ), &
  g5_deriv_x_yz( natom_x, natom, 3, num_g5_x_yz ), &
  g5_deriv_x_zz( natom_x, natom, 3, num_g5_x_zz )

integer nlayer_x
integer network_x( nlayer_x )
integer nweight_x, nhidden_x

double precision input_x( natom_x, 0 : network_x(1) )
double precision &
  dGdx( natom_x, natom, network_x(1) ), &
  dGdy( natom_x, natom, network_x(1) ), &
  dGdz( natom_x, natom, network_x(1) )

double precision weight_x( nweight_x ), hidden_x( natom_x, nhidden_x )
double precision dhidden_x( natom_x, nhidden_x )

double precision energy, energy0, denergy
double precision force( natom, 3 ), force0( natom, 3 )
double precision beta, const( natom_x ), const1

double precision g_e_x( nweight_x )
double precision g_f_x( nweight_x )
double precision,allocatable,dimension(:,:) :: dforce


 input_x = 0.0d0
 do i = 1, natom_x
   input_x(i,0) = 1.0d0
 enddo

 ! x_x 1
 do j = 1, num_g2_x_x
   do i = 1, natom_x
     input_x(i,j) = g2_x_x(i,j)
   enddo
 enddo
 ! x_y 2
 do j = 1, num_g2_x_y
   k = j + num_g_x(1)
   do i = 1, natom_x
     input_x(i,k) = g2_x_y(i,j)
   enddo
 enddo
 ! x_z 3
 do j = 1, num_g2_x_z
   k = j + sum( num_g_x(1:2) )
   do i = 1, natom_x
     input_x(i,k) = g2_x_z(i,j)
   enddo
 enddo
 ! x_xx 4
 do j = 1, num_g5_x_xx
   k = j + sum( num_g_x(1:3) )
   do i = 1, natom_x
     input_x(i,k) = g5_x_xx(i,j)
   enddo
 enddo
 ! x_xy 5
 do j = 1, num_g5_x_xy
   k = j + sum( num_g_x(1:4) )
   do i = 1, natom_x
     input_x(i,k) = g5_x_xy(i,j)
   enddo
 enddo
 ! x_xz 6
 do j = 1, num_g5_x_xz
   k = j + sum( num_g_x(1:5) )
   do i = 1, natom_x
     input_x(i,k) = g5_x_xz(i,j)
   enddo
 enddo
 ! x_yy 7
 do j = 1, num_g5_x_yy
   k = j + sum( num_g_x(1:6) )
   do i = 1, natom_x
     input_x(i,k) = g5_x_yy(i,j)
   enddo
 enddo
 ! x_yz 8
 do j = 1, num_g5_x_yz
   k = j + sum( num_g_x(1:7) )
   do i = 1, natom_x
     input_x(i,k) = g5_x_yz(i,j)
   enddo
 enddo
 ! x_zz 9
 do j = 1, num_g5_x_zz
   k = j + sum( num_g_x(1:8) )
   do i = 1, natom_x
     input_x(i,k) = g5_x_zz(i,j)
   enddo
 enddo

 
 dGdx = 0.0d0
 dGdy = 0.0d0
 dGdz = 0.0d0

 ! x_x 1
 do k = 1, num_g2_x_x
   do j = 1, natom
     do i = 1, natom_x
       dGdx(i,j,k) = g2_deriv_x_x(i,j,1,k)
       dGdy(i,j,k) = g2_deriv_x_x(i,j,2,k)
       dGdz(i,j,k) = g2_deriv_x_x(i,j,3,k)
     enddo ! i
   enddo ! j
 enddo ! k
 ! x_y 2
 do k = 1, num_g2_x_y
   m = k + num_g_x(1)
   do j = 1, natom
     do i = 1, natom_x
       dGdx(i,j,m) = g2_deriv_x_y(i,j,1,k)
       dGdy(i,j,m) = g2_deriv_x_y(i,j,2,k)
       dGdz(i,j,m) = g2_deriv_x_y(i,j,3,k)
     enddo
   enddo
 enddo
 ! x_z 3
 do k = 1, num_g2_x_z
   m = k + sum( num_g_x(1:2) )
   do j = 1, natom
     do i = 1, natom_x
       dGdx(i,j,m) = g2_deriv_x_z(i,j,1,k)
       dGdy(i,j,m) = g2_deriv_x_z(i,j,2,k)
       dGdz(i,j,m) = g2_deriv_x_z(i,j,3,k)
     enddo
   enddo
 enddo
 ! x_xx 4
 do k = 1, num_g5_x_xx
   m = k + sum( num_g_x(1:3) )
   do j = 1, natom
     do i = 1, natom_x
       dGdx(i,j,m) = g5_deriv_x_xx(i,j,1,k)
       dGdy(i,j,m) = g5_deriv_x_xx(i,j,2,k)
       dGdz(i,j,m) = g5_deriv_x_xx(i,j,3,k)
     enddo ! i
   enddo ! j
 enddo ! k
 ! x_xy 5
 do k = 1, num_g5_x_xy
   m = k + sum( num_g_x(1:4) )
   do j = 1, natom
     do i = 1, natom_x
       dGdx(i,j,m) = g5_deriv_x_xy(i,j,1,k)
       dGdy(i,j,m) = g5_deriv_x_xy(i,j,2,k)
       dGdz(i,j,m) = g5_deriv_x_xy(i,j,3,k)
     enddo ! i
   enddo ! j
 enddo
 ! x_xz 6
 do k = 1, num_g5_x_xz
   m = k + sum( num_g_x(1:5) )
   do j = 1, natom
     do i = 1, natom_x
       dGdx(i,j,m) = g5_deriv_x_xz(i,j,1,k)
       dGdy(i,j,m) = g5_deriv_x_xz(i,j,2,k)
       dGdz(i,j,m) = g5_deriv_x_xz(i,j,3,k)
     enddo ! i
   enddo ! j
 enddo ! k
 ! x_yy 7
 do k = 1, num_g5_x_yy
   m = k + sum( num_g_x(1:6) )
   do j = 1, natom
     do i = 1, natom_x
       dGdx(i,j,m) = g5_deriv_x_yy(i,j,1,k)
       dGdy(i,j,m) = g5_deriv_x_yy(i,j,2,k)
       dGdz(i,j,m) = g5_deriv_x_yy(i,j,3,k)
     enddo ! i
   enddo ! j
 enddo ! k
 ! x_yz 8
 do k = 1, num_g5_x_yz
   m = k + sum( num_g_x(1:7) )
   do j = 1, natom
     do i = 1, natom_x
       dGdx(i,j,m) = g5_deriv_x_yz(i,j,1,k)
       dGdy(i,j,m) = g5_deriv_x_yz(i,j,2,k)
       dGdz(i,j,m) = g5_deriv_x_yz(i,j,3,k)
     enddo ! i
   enddo ! j
 enddo ! k
 ! x_zz 9
 do k = 1, num_g5_x_zz
   m = k + sum( num_g_x(1:8) )
   do j = 1, natom
     do i = 1, natom_x
       dGdx(i,j,m) = g5_deriv_x_zz(i,j,1,k)
       dGdy(i,j,m) = g5_deriv_x_zz(i,j,2,k)
       dGdz(i,j,m) = g5_deriv_x_zz(i,j,3,k)
     enddo ! i
   enddo ! j
 enddo ! k


! ( 1.0d0 - hidden_x**2 )
 do j = 1, nhidden_x
   do i = 1, natom_x
     dhidden_x(i,j) = 1.0d0 - hidden_x(i,j)**2
   enddo
 enddo


 denergy = energy - energy0


!*****************
! 1 hidden layer
!*****************
!*****************
! 2 hidden layers
!*****************
 if( nlayer_x == 4 )then

   w2_0 = ( 1 + network_x(1) ) * network_x(2)
   w3_0 = w2_0 + ( 1 + network_x(2) ) * network_x(3)
   h2_0 = 1 + network_x(2)

     !-----------
     ! 1st layer
     !-----------
     do l2 = 1, network_x(2) !node_h1_x
       h1 = 1 + l2
       do l1 = 0, network_x(1) !node_in_x
         m = ( l1 + 1 ) + ( network_x(1) + 1 )*( l2 - 1 )
         do l3 = 1, network_x(3) !node_h2_x
           h2 = h2_0 + ( 1 + l3 )
           w2 = w2_0 + ( 1 + network_x(2) )*( l3 - 1 ) + ( 1 + l2 )
           do l4 = 1, network_x(4) !nodeout_x node_h3_x
             w3 = w3_0 + ( 1 + network_x(3) )*( l4 - 1 ) + ( 1 + l3 )
             do i = 1, natom_x
               g_e_x(m) =   g_e_x(m) &
                          + 2.0d0 * denergy &
                            * weight_x(w3) * dhidden_x(i,h2) &
                            * weight_x(w2) * dhidden_x(i,h1) &
                            * input_x(i,l1)
             enddo ! i
           enddo ! l4
         enddo ! l3
       enddo ! l1
     enddo ! l2
     !-----------
     ! 2nd layer
     !-----------
     do l3 = 1, network_x(3) !node_h2_x
       h2 = h2_0 + ( 1 + l3 )
       do l2 = 0, network_x(2) !node_h1_x
         h1 = 1 + l2
         m = w2_0 + ( l2 + 1 ) + ( network_x(2) + 1 )*( l3 - 1 )
         do l4 = 1, network_x(4) !node_out_x node_h3_x
           w3 = w3_0 + ( 1 + network_x(3) )*( l4 - 1 ) + ( 1 + l3 )
           do i = 1, natom_x
             g_e_x(m) =   g_e_x(m) &
                        + 2.0d0 * denergy &
                          * weight_x(w3) * dhidden_x(i,h2) &
                          * hidden_x(i,h1)
           enddo ! i
         enddo ! l4
       enddo ! l2
     enddo ! l3
     !-----------
     ! 3rd layer
     !-----------
     do l4 = 1, network_x(4) !node_out_x node_h3_x
       do l3 = 0, network_x(3) !node_h2_x
         m = w3_0 + ( l3 + 1 ) + ( network_x(3) + 1 )*( l4 - 1 )
         h2 = ( 1 + network_x(2) ) + ( 1 + l3 )
         do i = 1, natom_x
           g_e_x(m) =   g_e_x(m) &
                      + 2.0d0 * denergy &
                        * hidden_x(i,h2)
         enddo ! i
       enddo ! l3
     enddo ! l4


   !-------
   ! Force
   !-------
   if( beta /= 0.0d0 )then

     allocate( dforce( natom_x, network_x(1) ) )
     dforce = 0.0d0
     do l1 = 1, network_x(1)
     do i = 1, natom
     do j = 1, natom_x
       dforce(j,l1) = dforce(j,l1) + 2.0d0*( force(i,1) - force0(i,1) )*dGdx(j,i,l1)
       dforce(j,l1) = dforce(j,l1) + 2.0d0*( force(i,2) - force0(i,2) )*dGdy(j,i,l1)
       dforce(j,l1) = dforce(j,l1) + 2.0d0*( force(i,3) - force0(i,3) )*dGdz(j,i,l1)
     enddo
     enddo
     enddo


     w2_0 = ( 1 + network_x(1) ) * network_x(2)
     w3_0 = w2_0 + ( 1 + network_x(2) ) * network_x(3)
     h2_0 = 1 + network_x(2)

     do ll1 = 1, network_x(1)
     do ll2 = 1, network_x(2)
       hh1 = ( 1 + ll2 )
       ww1 = ( 1 + network_x(1) )*( ll2 - 1 ) + ( 1 + ll1 )
     do ll3 = 1, network_x(3)
       hh2 = h2_0 + ( 1 + ll3 )
       ww2 = w2_0 + ( 1 + network_x(2) )*( ll3 - 1 ) + ( 1 + ll2 )
     do ll4 = 1, network_x(4)
       ww3 = w3_0 + ( 1 + ll3 )
  
       const1 = weight_x(ww1) * weight_x(ww2) * weight_x(ww3)

       do j = 1, natom_x
         const(j) = const1 * dhidden_x(j,hh1) * dhidden_x(j,hh2) * dforce(j,ll1)
       enddo ! j
  
       !-----------
       ! 1st layer
       !-----------
       !bias l1 = 0
       do l2 = 1, network_x(2)
         h1 = ( 1 + l2 )
         w1 = ( 1 + network_x(1) )*( l2 - 1 ) + 1 !( 1 + l1 )
         www2 = w2_0 + ( 1 + network_x(2) )*( ll3 - 1 ) + ( 1 + l2 )
         do j = 1, natom_x
           g_f_x(w1) =   g_f_x(w1) &
                       + const(j) * 2.0d0 * hidden_x(j,hh2) &
                                  * weight_x(www2) &
                                  * dhidden_x(j,h1) 
         enddo ! j
         if( h1 == hh1 )then
           do j = 1, natom_x
             g_f_x(w1) =   g_f_x(w1) &
                         + const(j) * 2.0d0 * hidden_x(j,h1)
           enddo ! j
         endif
       enddo ! l2

       !weight
       do l2 = 1, network_x(2)
         h1 = ( 1 + l2 )
         www2 = w2_0 + ( 1 + network_x(2) )*( ll3 - 1 ) + ( 1 + l2 )
         do l1 = 1, network_x(1) !0, network_x(1)
           w1 = ( 1 + network_x(1) )*( l2 - 1 ) + ( 1 + l1 )
           do j = 1, natom_x
             g_f_x(w1) =   g_f_x(w1) &
                         + ( 2.0d0 * const(j) &
                             * hidden_x(j,hh2) &
                             * weight_x(www2) &
                             * dhidden_x(j,h1) &
                             * input_x(j,l1) )
           enddo ! j
           if( h1 == hh1 )then
             do j = 1, natom_x
               g_f_x(w1) =   g_f_x(w1) &
                           + ( 2.0d0 * const(j) * hidden_x(j,h1) * input_x(j,l1) )
             enddo ! j
           endif
           if( h1 == hh1 .and. l1 == ll1 )then
             do j = 1, natom_x
               g_f_x(w1) =   g_f_x(w1) &
                           - ( const(j) * ( 1.0d0 / weight_x(w1) ) ) !* input_x(j,l1)
             enddo ! j
           endif
         enddo ! l1
       enddo ! l2

       !-----------
       ! 2nd layer
       !-----------
       !bias l2 = 0
       do l3 = 1, network_x(3)
         if( l3 /= ll3 ) cycle
         h1 = 1 !( 1 + l2 )
         h2 = h2_0 + ( 1 + l3 )
         w2 = w2_0 + ( 1 + network_x(2) )*( l3 - 1 ) + 1 !( 1 + l2 )
         do j = 1, natom_x
           g_f_x(w2) =   g_f_x(w2) &
                       + ( 2.0d0 * const(j) * hidden_x(j,h2) )
         enddo!j
       enddo ! l3
  
       !weight
       do l3 = 1, network_x(3)
         if( l3 /= ll3 ) cycle
         do l2 = 1, network_x(2) !0, network_x(2)
           h1 = ( 1 + l2 )
           h2 = h2_0 + ( 1 + l3 )
           w2 = w2_0 + ( 1 + network_x(2) )*( l3 - 1 ) + ( 1 + l2 )
           do j = 1, natom_x
             g_f_x(w2) =   g_f_x(w2) &
                         + 2.0d0 * const(j) * hidden_x(j,h2) * hidden_x(j,h1) 
           enddo ! j
           if( h1 == hh1 .and. h2 == hh2 )then
             do j = 1, natom_x
               g_f_x(w2) =   g_f_x(w2) &
                           - ( const(j) * ( 1.0d0 / weight_x(w2) ) )
             enddo ! j
           endif
         enddo ! l2
       enddo ! l3
  
       !-----------
       ! 3rd layer
       !-----------
       do l4 = 1, network_x(4)
         if( l4 /= ll4 ) cycle
         do l3 = 1, network_x(3) !0, network_x(3)
           if( l3 /= ll3 ) cycle
           w3 = w3_0 + ( 1 + l3 )
             do j = 1, natom_x
               g_f_x(w3) =   g_f_x(w3) &
                           - ( const(j) * ( 1.0d0 / weight_x(w3) ) )
             enddo ! j
         enddo ! l3
       enddo ! l4
  
     enddo ! ll4
     enddo ! ll3
     enddo ! ll2
     enddo ! ll1
  
     deallocate( dforce )


   endif ! beta /= 0

 endif ! 2 hidden layers


!*****************
! 3 hidden layers
!*****************


end subroutine f_ef_deriv3abc
