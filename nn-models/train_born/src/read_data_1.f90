subroutine read_data_1( num_atom_x, num_data_x, natom_x_total, natom_x_sum, &
                        counter_x, counter_natom_x, counter_deriv_x, &
                        imd, natom_x, natom, &
                        tag2_x_x,  g2_x_x,  g2_deriv_x_x,  g2_x_x_maxmin,  num_g2_x_x, &
                        tag5_x_xx, g5_x_xx, g5_deriv_x_xx, g5_x_xx_maxmin, num_g5_x_xx )
implicit none
integer i, j, k, m, imd
character(len=120) filename
character(len=6)   tag2_x_x
character(len=7)   tag5_x_xx
integer num_data_x
integer natom, natom_x
integer natom_x_total
integer natom_x_sum
integer num_atom_x( num_data_x )
integer counter_x
integer counter_natom_x
integer counter_deriv_x
integer num_g2_x_x, num_g5_x_xx
double precision,allocatable,dimension(:,:) :: g2, g5
double precision :: &
  g2_x_x(  natom_x_total, num_g2_x_x ), &
  g5_x_xx( natom_x_total, num_g5_x_xx )
double precision :: &
  g2_x_x_maxmin(  num_g2_x_x, 2 ), &
  g5_x_xx_maxmin( num_g5_x_xx, 2 )
double precision,allocatable,dimension(:,:,:,:) :: g2_deriv, g5_deriv
double precision :: &
  g2_deriv_x_x(  natom_x_sum, 3, num_g2_x_x ), &
  g5_deriv_x_xx( natom_x_sum, 3, num_g5_x_xx )


  num_atom_x( imd + counter_x ) = natom_x

  ! SF
  allocate( g2( natom_x, num_g2_x_x ) )
  read(11) g2
  do j = 1, num_g2_x_x
    do i = 1, natom_x
      g2_x_x( i + counter_natom_x, j ) = ( 2.0d0*( g2(i,j) - g2_x_x_maxmin(j,2) )/ &
                                         ( g2_x_x_maxmin(j,1) - g2_x_x_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g2 )

  allocate( g5( natom_x, num_g5_x_xx ) )
  read(11) g5
  do j = 1, num_g5_x_xx
    do i = 1, natom_x
      g5_x_xx( i + counter_natom_x, j ) = ( 2.0d0*( g5(i,j) - g5_x_xx_maxmin(j,2) )/ &
                                          ( g5_x_xx_maxmin(j,1) - g5_x_xx_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  counter_natom_x = counter_natom_x + natom_x


  ! SF_deriv
  allocate( g2_deriv( natom_x, natom, 3, num_g2_x_x ) )
  read(12) g2_deriv
  do k = 1, num_g2_x_x
    do j = 1, natom
      do i = 1, natom_x
        m = i + natom_x * ( j - 1 )
        g2_deriv_x_x( counter_deriv_x + m,1,k) = 2.0d0*g2_deriv(i,j,1,k)/ &
                                                 ( g2_x_x_maxmin(k,1) - g2_x_x_maxmin(k,2) )
        g2_deriv_x_x( counter_deriv_x + m,2,k) = 2.0d0*g2_deriv(i,j,2,k)/ &
                                                 ( g2_x_x_maxmin(k,1) - g2_x_x_maxmin(k,2) )
        g2_deriv_x_x( counter_deriv_x + m,3,k) = 2.0d0*g2_deriv(i,j,3,k)/ &
                                                 ( g2_x_x_maxmin(k,1) - g2_x_x_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g2_deriv )

  allocate( g5_deriv( natom_x, natom, 3, num_g5_x_xx ) )
  read(12) g5_deriv
  do k = 1, num_g5_x_xx
    do j = 1, natom
      do i = 1, natom_x
        m = i + natom_x * ( j - 1 )
        g5_deriv_x_xx( counter_deriv_x + m ,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                   ( g5_x_xx_maxmin(k,1) - g5_x_xx_maxmin(k,2) )
        g5_deriv_x_xx( counter_deriv_x + m ,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                   ( g5_x_xx_maxmin(k,1) - g5_x_xx_maxmin(k,2) )
        g5_deriv_x_xx( counter_deriv_x + m ,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                   ( g5_x_xx_maxmin(k,1) - g5_x_xx_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  counter_deriv_x = counter_deriv_x + natom_x * natom


end subroutine read_data_1
