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
           network_x, nlayer_x, weight_x, weight_inv_x, nweight_x, hidden_x, nhidden_x, &
           beta, g_e_x, g_f_x, &
           natom, natom_x, energy, energy0, force, force0 )
implicit none
integer i, j, k, l, m, m0, di
integer l1, l2, l3, l4
integer h1, h2, h1_1, h1_2
integer w1, w2, w3
integer w1_1, w1_2, w2_1, w2_2
integer h2_0
integer w2_0, w3_0
integer ll1, ll2, ll3, ll4
integer hh1, hh2
integer ww1, ww2, ww3, ww1_1, ww1_2
integer www2, www2_1, www2_2
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

double precision weight_x( nweight_x ), weight_inv_x( nweight_x )
double precision hidden_x( natom_x, nhidden_x )
double precision dhidden_x( natom_x, nhidden_x )

double precision dummy, dummy_wa

double precision energy, energy0, denergy, two_denergy
double precision force( natom, 3 ), force0( natom, 3 ), dforce0( natom, 3 )
double precision beta
double precision const_a, const_b, const_c, const_d, const_e
double precision inv_1, inv_2, inv_3
double precision const_w1( natom_x ), const_w2( natom_x ), const_w3( natom_x )
double precision const( natom_x ), const1( natom_x ), const2( natom_x ), &
                 const3( natom_x ), const4( natom_x )

double precision g_e_x( nweight_x )
double precision g_f_x( nweight_x )
double precision dforce( natom_x, network_x(1) )

double precision dhidden_x_t( nhidden_x, natom_x )
double precision array_i2( network_x(2), natom_x ), array_o2( network_x(2) )
double precision input_x_t( 0 : network_x(1), natom_x  )
double precision array_i1( network_x(1), natom_x ), array_o1( network_x(1) )
double precision hidden_x_t( nhidden_x, natom_x )
double precision array_i3( network_x(2), natom_x )

double precision d1, d2, d3, d4, d5, d6
double precision t0, t1, t2


!print*,"start f_deriv"
!call cpu_time(t0)


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
 di = num_g_x(1)
 do j = 1, num_g2_x_y
   k = j + di
   do i = 1, natom_x
     input_x(i,k) = g2_x_y(i,j)
   enddo
 enddo
 ! x_z 3
 di = sum( num_g_x(1:2) )
 do j = 1, num_g2_x_z
   k = j + di
   do i = 1, natom_x
     input_x(i,k) = g2_x_z(i,j)
   enddo
 enddo
 ! x_xx 4
 di = sum( num_g_x(1:3) )
 do j = 1, num_g5_x_xx
   k = j + di
   do i = 1, natom_x
     input_x(i,k) = g5_x_xx(i,j)
   enddo
 enddo
 ! x_xy 5
 di = sum( num_g_x(1:4) )
 do j = 1, num_g5_x_xy
   k = j + di
   do i = 1, natom_x
     input_x(i,k) = g5_x_xy(i,j)
   enddo
 enddo
 ! x_xz 6
 di = sum( num_g_x(1:5) )
 do j = 1, num_g5_x_xz
   k = j + di
   do i = 1, natom_x
     input_x(i,k) = g5_x_xz(i,j)
   enddo
 enddo
 ! x_yy 7
 di = sum( num_g_x(1:6) )
 do j = 1, num_g5_x_yy
   k = j + di
   do i = 1, natom_x
     input_x(i,k) = g5_x_yy(i,j)
   enddo
 enddo
 ! x_yz 8
 di = sum( num_g_x(1:7) )
 do j = 1, num_g5_x_yz
   k = j + di
   do i = 1, natom_x
     input_x(i,k) = g5_x_yz(i,j)
   enddo
 enddo
 ! x_zz 9
 di = sum( num_g_x(1:8) )
 do j = 1, num_g5_x_zz
   k = j + di
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
 di = num_g_x(1)
 do k = 1, num_g2_x_y
   m = k + di
   do j = 1, natom
     do i = 1, natom_x
       dGdx(i,j,m) = g2_deriv_x_y(i,j,1,k)
       dGdy(i,j,m) = g2_deriv_x_y(i,j,2,k)
       dGdz(i,j,m) = g2_deriv_x_y(i,j,3,k)
     enddo
   enddo
 enddo
 ! x_z 3
 di = sum( num_g_x(1:2) )
 do k = 1, num_g2_x_z
   m = k + di
   do j = 1, natom
     do i = 1, natom_x
       dGdx(i,j,m) = g2_deriv_x_z(i,j,1,k)
       dGdy(i,j,m) = g2_deriv_x_z(i,j,2,k)
       dGdz(i,j,m) = g2_deriv_x_z(i,j,3,k)
     enddo
   enddo
 enddo
 ! x_xx 4
 di = sum( num_g_x(1:3) )
 do k = 1, num_g5_x_xx
   m = k + di
   do j = 1, natom
     do i = 1, natom_x
       dGdx(i,j,m) = g5_deriv_x_xx(i,j,1,k)
       dGdy(i,j,m) = g5_deriv_x_xx(i,j,2,k)
       dGdz(i,j,m) = g5_deriv_x_xx(i,j,3,k)
     enddo ! i
   enddo ! j
 enddo ! k
 ! x_xy 5
 di = sum( num_g_x(1:4) )
 do k = 1, num_g5_x_xy
   m = k + di
   do j = 1, natom
     do i = 1, natom_x
       dGdx(i,j,m) = g5_deriv_x_xy(i,j,1,k)
       dGdy(i,j,m) = g5_deriv_x_xy(i,j,2,k)
       dGdz(i,j,m) = g5_deriv_x_xy(i,j,3,k)
     enddo ! i
   enddo ! j
 enddo
 ! x_xz 6
 di = sum( num_g_x(1:5) )
 do k = 1, num_g5_x_xz
   m = k + di
   do j = 1, natom
     do i = 1, natom_x
       dGdx(i,j,m) = g5_deriv_x_xz(i,j,1,k)
       dGdy(i,j,m) = g5_deriv_x_xz(i,j,2,k)
       dGdz(i,j,m) = g5_deriv_x_xz(i,j,3,k)
     enddo ! i
   enddo ! j
 enddo ! k
 ! x_yy 7
 di = sum( num_g_x(1:6) )
 do k = 1, num_g5_x_yy
   m = k + di
   do j = 1, natom
     do i = 1, natom_x
       dGdx(i,j,m) = g5_deriv_x_yy(i,j,1,k)
       dGdy(i,j,m) = g5_deriv_x_yy(i,j,2,k)
       dGdz(i,j,m) = g5_deriv_x_yy(i,j,3,k)
     enddo ! i
   enddo ! j
 enddo ! k
 ! x_yz 8
 di = sum( num_g_x(1:7) )
 do k = 1, num_g5_x_yz
   m = k + di
   do j = 1, natom
     do i = 1, natom_x
       dGdx(i,j,m) = g5_deriv_x_yz(i,j,1,k)
       dGdy(i,j,m) = g5_deriv_x_yz(i,j,2,k)
       dGdz(i,j,m) = g5_deriv_x_yz(i,j,3,k)
     enddo ! i
   enddo ! j
 enddo ! k
 ! x_zz 9
 di = sum( num_g_x(1:8) )
 do k = 1, num_g5_x_zz
   m = k + di
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
 two_denergy = 2.0d0 * denergy


 input_x_t   = transpose( input_x )
 hidden_x_t  = transpose( hidden_x )
 dhidden_x_t = transpose( dhidden_x )


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
     do l2 = 1, network_x(2)
       h1 = 1 + l2
       m0 = ( network_x(1) + 1 )*( l2 - 1 )
       do l1 = 0, network_x(1)
         m = ( l1 + 1 ) + m0 !( l1 + 1 ) + ( network_x(1) + 1 )*( l2 - 1 )
         do l3 = 1, network_x(3)
           h2 = h2_0 + ( 1 + l3 )
           w2 = w2_0 + ( 1 + network_x(2) )*( l3 - 1 ) + ( 1 + l2 )
           do l4 = 1, network_x(4)
             w3 = w3_0 + ( 1 + network_x(3) )*( l4 - 1 ) + ( 1 + l3 )

             dummy_wa = 0.0d0
             do j = 1, natom_x
               dummy_wa = dummy_wa + ( dhidden_x(j,h1) * dhidden_x(j,h2) * input_x(j,l1) )
             enddo
             dummy    = two_denergy * weight_x(w3) * weight_x(w2)
             g_e_x(m) = g_e_x(m) + ( dummy * dummy_wa )

             !do i = 1, natom_x
               !g_e_x(m) =     g_e_x(m) &
                          !+ ( 2.0d0 * denergy * weight_x(w3) * dhidden_x(i,h2) &
                          !                    * weight_x(w2) * dhidden_x(i,h1) * input_x(i,l1) )
             !enddo ! i

           enddo ! l4
         enddo ! l3
       enddo ! l1
     enddo ! l2

     !-----------
     ! 2nd layer
     !-----------
     do l3 = 1, network_x(3)
       h2 = h2_0 + ( 1 + l3 )
       m0 = w2_0 + ( network_x(2) + 1 )*( l3 - 1 )
       do l2 = 0, network_x(2)
         h1 = 1 + l2
         m = m0 + ( l2 + 1 ) !w2_0 + ( l2 + 1 ) + ( network_x(2) + 1 )*( l3 - 1 )
         do l4 = 1, network_x(4)
           w3 = w3_0 + ( 1 + network_x(3) )*( l4 - 1 ) + ( 1 + l3 )

           dummy_wa = 0.0d0
           do j = 1, natom_x
             dummy_wa = dummy_wa + ( dhidden_x(j,h2) * hidden_x(j,h1) )
           enddo
           dummy    = two_denergy * weight_x(w3)
           g_e_x(m) = g_e_x(m) + ( dummy * dummy_wa )

           !do i = 1, natom_x
           !  g_e_x(m) =   g_e_x(m) &
           !             !+ 2.0d0 * denergy &
           !             + ( two_denergy &
           !                 * weight_x(w3) * dhidden_x(i,h2) &
           !                 * hidden_x(i,h1) )
           !enddo ! i

         enddo ! l4
       enddo ! l2
     enddo ! l3

     !-----------
     ! 3rd layer
     !-----------
     do l4 = 1, network_x(4)
       do l3 = 0, network_x(3)
         m = w3_0 + ( l3 + 1 ) + ( network_x(3) + 1 )*( l4 - 1 )
         h2 = ( 1 + network_x(2) ) + ( 1 + l3 )

         dummy_wa = 0.0d0
         do j = 1, natom_x
           dummy_wa = dummy_wa + hidden_x(j,h2)
         enddo
         g_e_x(m) = g_e_x(m) + ( two_denergy * dummy_wa )

         !do i = 1, natom_x
         !  g_e_x(m) =   g_e_x(m) &
         !             !+ 2.0d0 * denergy &
         !             + ( two_denergy &
         !                 * hidden_x(i,h2) )
         !enddo ! i

       enddo ! l3
     enddo ! l4


   !-------
   ! Force
   !-------
   if( beta /= 0.0d0 )then


     do j = 1, natom_x
       do l2 = 1, network_x(2)
         array_i2(l2,j) = dhidden_x_t(l2+1,j)
       enddo
     enddo


     do j = 1, natom_x
       do l1 = 1, network_x(1)
         array_i1(l1,j) = input_x_t(l1,j)
       enddo
     enddo


     do j = 1, natom_x
       do l2 = 1, network_x(2)
         array_i3(l2,j) = hidden_x_t(l2+1,j)
       enddo
     enddo




     dforce = 0.0d0
     do l1 = 1, network_x(1)
       do i = 1, natom
           d1 = 4.0d0 * ( force(i,1) - force0(i,1) ) !2.0d0 * ( force(i,1) - force0(i,1) )
           d2 = 4.0d0 * ( force(i,2) - force0(i,2) ) !2.0d0 * ( force(i,2) - force0(i,2) )
           d3 = 4.0d0 * ( force(i,3) - force0(i,3) ) !2.0d0 * ( force(i,3) - force0(i,3) )

         do j = 1, natom_x
           dforce(j,l1) = dforce(j,l1) + ( d1*dGdx(j,i,l1) + d2*dGdy(j,i,l1) + d3*dGdz(j,i,l1) )
           !dforce(j,l1+1) = dforce(j,l1+1) + ( d1*dGdx(j,i,l1+1) + d2*dGdy(j,i,l1+1) + d3*dGdz(j,i,l1+1) )
           !dforce(j,l1) =   dforce(j,l1) &
           !               + 2.0d0*(   ( force(i,1) - force0(i,1) )*dGdx(j,i,l1) &
           !                         + ( force(i,2) - force0(i,2) )*dGdy(j,i,l1) &
           !                         + ( force(i,3) - force0(i,3) )*dGdz(j,i,l1) )
         enddo ! j

       enddo ! i
     enddo ! l1


     do ll1 = 1, network_x(1)
       do j = 1, natom_x
         const_w1(j) = dforce(j,ll1) !2.0d0 * dforce(j,ll1)
       enddo
     do ll2 = 1, network_x(2)
       hh1 = ( 1 + ll2 )
       ww1 = ( 1 + network_x(1) )*( ll2 - 1 ) + ( 1 + ll1 )
       d1 = weight_x(ww1)
       inv_1 = weight_inv_x(ww1) !0.50d0 / d1 !weight_x(ww1)
       do j = 1, natom_x
         const_w2(j) = const_w1(j) * dhidden_x(j,hh1) * d1 !weight_x(ww1)
       enddo
       ww1_1 = ( 1 + network_x(1) )*( ll2 - 1 )
     do ll3 = 1, network_x(3)
       hh2 = h2_0 + ( 1 + ll3 )
       ww2 = w2_0 + ( 1 + network_x(2) )*( ll3 - 1 ) + ( 1 + ll2 )
       d2 = weight_x(ww2)
       inv_2 = weight_inv_x(ww2) !0.50d0 / d2 !weight_x(ww2)
       do j = 1, natom_x
         const_w3(j) = const_w2(j) * dhidden_x(j,hh2) * d2 !weight_x(ww2)
       enddo
       w2_1 = w2_0 + ( 1 + network_x(2) )*( ll3 - 1 )
     do ll4 = 1, network_x(4)
       ww3 = w3_0 + ( 1 + ll3 )
       d3 = weight_x(ww3)
       inv_3 = weight_inv_x(ww3) !0.50d0 / d3 !weight_x(ww3)

       do j = 1, natom_x
         const(j) = const_w3(j) * d3 !weight_x(ww3)
       enddo

       !const =   weight_x(ww3) * ( 1.0d0 - hidden_x(j,hh2)**2 ) &
       !        * weight_x(ww2) * ( 1.0d0 - hidden_x(j,hh1)**2 ) &
       !        * weight_x(ww1) &
       !        * dforce(j,ll1)

       const_b = 0.0d0
       do j = 1, natom_x
         d1        = const(j) * hidden_x(j,hh1)
         const3(j) = d1 !const(j) * hidden_x(j,hh1)
         const_b   = const_b + d1
       enddo

       const_a = 0.0d0
       do j = 1, natom_x
         d1        = const(j) * hidden_x(j,hh2)
         const2(j) = d1 !const(j) * hidden_x(j,hh2)
         const_a   = const_a + d1 
       enddo ! j

       const_c = 0.0d0
       const_d = 0.0d0
       const_e = 0.0d0
       do j = 1, natom_x
         const_c = const_c + ( const(j) * inv_1 )
         const_d = const_d + ( const(j) * inv_2 )
         const_e = const_e + ( const(j) * inv_3 )
       enddo
       

       !-----------
       ! 1st layer
       !-----------
       !BIAS, l1 = 0
       CALL DGEMV('N', network_x(2), natom_x, 1.0d0, array_i2, network_x(2), &
                  const2, 1, 0.0d0, array_o2, 1)
       www2_1 = w2_0 + ( 1 + network_x(2) )*( ll3 - 1 )
       do l2 = 1, network_x(2) 
         h1 = ( 1 + l2 )
         w1 = ( 1 + network_x(1) )*( l2 - 1 ) + 1 !( 1 + l1 )
         www2 = www2_1 + ( 1 + l2 ) !w2_0 + ( 1 + network_x(2) )*( ll3 - 1 ) + ( 1 + l2 )
         g_f_x(w1) = g_f_x(w1) + ( array_o2(l2) * weight_x(www2) )
         !do j = 1, natom_x
           !g_f_x(w1) =   g_f_x(w1) &
           !            + 2.0*const(j)*hidden_x(j,hh2)*weight_x(www2)*dhidden_x(j,h1) 
         !enddo
       enddo ! l2
       !h1 == hh1 : l2 == ll2
       !w1 = ( 1 + network_x(1) )*( ll2 - 1 ) + 1 
       w1 = ww1_1 + 1 
       g_f_x(w1) = g_f_x(w1) + const_b

       !WEIGHT
       do l2 = 1, network_x(2)
         www2 = www2_1 + ( 1 + l2 ) !w2_0 + ( 1 + network_x(2) )*( ll3 - 1 ) + ( 1 + l2 )
         h1 = 1 + l2
         do j = 1, natom_x
           const1(j)  = const2(j) * dhidden_x(j,h1) ! * weight_x(www2)
         enddo ! j

         CALL DGEMV('N', network_x(1), natom_x, 1.0d0, array_i1, network_x(1), &
                    const1, 1, 0.0d0, array_o1, 1)
         w1_1 = ( 1 + network_x(1) )*( l2 - 1 )
         do l1 = 1, network_x(1)
           w1 = w1_1 + ( 1 + l1 ) !( 1 + network_x(1) )*( l2 - 1 ) + ( 1 + l1 )
           g_f_x(w1) = g_f_x(w1)  + ( array_o1(l1)  * weight_x(www2) )
           !g_f_x(w1) = g_f_x(w1) + ( 2.0*const(j)*hidden_x(j,hh2)*weight_x(www2)*dhidden_x(j,h1)*input_x(j,l1) )
         enddo ! l1

       enddo ! l2

       !h1 == hh1 : l2 == ll2
       CALL DGEMV('N', network_x(1), natom_x, 1.0d0, array_i1, network_x(1), &
                  const3, 1, 0.0d0, array_o1, 1)
       do l1 = 1, network_x(1) !0, network_x(1) !*****
         w1_1 = ( 1 + network_x(1) )*( ll2 - 1 ) + ( 1 + l1 )
         g_f_x(w1_1) = g_f_x(w1_1) + array_o1(l1)
       enddo ! l1

       !h1 == hh1 : l2 == ll2 .and. l1 == ll1
       g_f_x(ww1) = g_f_x(ww1) - const_c !dummy_wa
       !g_f_x(ww1) = g_f_x(ww1) - ( 0.5*const(j)/weight_x(ww1) ) !* input_x(j,l1)


       !-----------
       ! 2nd layer
       !-----------
       !BIAS, l2 = 0
       !do l3 = 1, network_x(3)
         !if( l3 /= ll3 ) cycle
         w2 = w2_1 + 1 !w2_0 + ( 1 + network_x(2) )*( ll3 - 1 ) + 1 !( 1 + l2 )
         g_f_x(w2) = g_f_x(w2) + const_a
         !do j = 1, natom_x
         !  g_f_x(w2) = g_f_x(w2) + ( const(j) * hidden_x(j,hh2) )
         !enddo!j
       !enddo ! l3
  
       !WEIGHT
       !do l3 = 1, network_x(3)
         !if( l3 /= ll3 ) cycle
         CALL DGEMV('N', network_x(2), natom_x, 1.0d0, array_i3, network_x(2), &
                    const2, 1, 0.0d0, array_o2, 1)
         do l2 = 1, network_x(2)
           h1 = ( 1 + l2 )
           w2 = w2_1 + ( 1 + l2 ) !w2_0 + ( 1 + network_x(2) )*( ll3 - 1 ) + ( 1 + l2 )
           g_f_x(w2) = g_f_x(w2) + array_o2(l2)
           !g_f_x(w2)   = g_f_x(w2) + ( const(j)*hidden_x(j,h2)*hidden_x(j,h1) )
         enddo ! l2
       !enddo ! l3
  
          !h1 == hh1 : l2 == ll2 .and. h2 == hh2 : l3 == ll3
          g_f_x(ww2) = g_f_x(ww2) - const_d !dummy_wa
          !g_f_x(ww2) = g_f_x(ww2) - ( 0.5*const(j)/weight_x(ww2) )


       !-----------
       ! 3rd layer
       !-----------
       !do l4 = 1, network_x(4)
         !if( l4 /= ll4 ) cycle
         !do l3 = 1, network_x(3) !0, network_x(3)
           !if( l3 /= ll3 ) cycle
             g_f_x(ww3) = g_f_x(ww3) - const_e !dummy_wa
             !g_f_x(w3)  = g_f_x(w3) - ( 0.5*const(j)/weight_x(w3) )
         !enddo ! l3
       !enddo ! l4
  
     enddo ! ll4
     enddo ! ll3
     enddo ! ll2
     enddo ! ll1
  

   endif ! beta /= 0

 endif ! 2 hidden layers


!*****************
! 3 hidden layers
!*****************


!call cpu_time(t1)
!print*,"g_f_x",t1-t0,"sec"


end subroutine f_ef_deriv3abc
