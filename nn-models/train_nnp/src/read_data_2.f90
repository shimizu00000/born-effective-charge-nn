subroutine read_data_2( num_atom_x, num_atom_y, &
                        num_data_x, num_data_y, &
                        natom_x_total, natom_xy_total, &
                        natom_y_total, natom_yx_total, &
                        natom_x_sum, natom_xy_sum, &
                        natom_y_sum, natom_yx_sum, &
                        counter_x, counter_y, &
                        counter_natom_x, counter_natom_xy, &
                        counter_natom_y, counter_natom_yx, &
                        counter_deriv_x, counter_deriv_xy, &
                        counter_deriv_y, counter_deriv_yx, &
                        imd, natom_x, natom_y, natom, &
                        tag2_x_x,  g2_x_x,  g2_deriv_x_x,  g2_x_x_maxmin,  num_g2_x_x, &
                        tag2_x_y,  g2_x_y,  g2_deriv_x_y,  g2_x_y_maxmin,  num_g2_x_y, &
                        tag5_x_xx, g5_x_xx, g5_deriv_x_xx, g5_x_xx_maxmin, num_g5_x_xx, &
                        tag5_x_xy, g5_x_xy, g5_deriv_x_xy, g5_x_xy_maxmin, num_g5_x_xy, &
                        tag5_x_yy, g5_x_yy, g5_deriv_x_yy, g5_x_yy_maxmin, num_g5_x_yy, &
                        tag2_y_x,  g2_y_x,  g2_deriv_y_x,  g2_y_x_maxmin,  num_g2_y_x, &
                        tag2_y_y,  g2_y_y,  g2_deriv_y_y,  g2_y_y_maxmin,  num_g2_y_y, &
                        tag5_y_xx, g5_y_xx, g5_deriv_y_xx, g5_y_xx_maxmin, num_g5_y_xx, &
                        tag5_y_xy, g5_y_xy, g5_deriv_y_xy, g5_y_xy_maxmin, num_g5_y_xy, &
                        tag5_y_yy, g5_y_yy, g5_deriv_y_yy, g5_y_yy_maxmin, num_g5_y_yy )
implicit none
integer i, j, k, m, imd
character(len=120) filename
character(len=6) tag2_x_x,  tag2_x_y
character(len=7) tag5_x_xx, tag5_x_xy, tag5_x_yy
character(len=6) tag2_y_x,  tag2_y_y
character(len=7) tag5_y_xx, tag5_y_xy, tag5_y_yy
integer num_data_x, num_data_y
integer natom, natom_x, natom_y
integer natom_x_total, natom_xy_total, &
        natom_y_total, natom_yx_total
integer natom_x_sum, natom_xy_sum, &
        natom_y_sum, natom_yx_sum
integer num_atom_x( num_data_x ), &
        num_atom_y( num_data_y )
integer counter_x, counter_y
integer counter_natom_x, counter_natom_xy, &
        counter_natom_y, counter_natom_yx
integer counter_deriv_x, counter_deriv_xy, &
        counter_deriv_y, counter_deriv_yx
integer num_g2_x_x,  num_g2_x_y, &
        num_g5_x_xx, num_g5_x_xy, num_g5_x_yy
integer num_g2_y_x,  num_g2_y_y, &
        num_g5_y_xx, num_g5_y_xy, num_g5_y_yy
double precision,allocatable,dimension(:,:) :: g2, g5
double precision :: &
  g2_x_x(  natom_x_total,  num_g2_x_x ), &
  g2_x_y(  natom_xy_total, num_g2_x_y ), &
  g5_x_xx( natom_x_total,  num_g5_x_xx ), &
  g5_x_xy( natom_xy_total, num_g5_x_xy ), &
  g5_x_yy( natom_xy_total, num_g5_x_yy )
double precision :: &
  g2_y_x(  natom_yx_total, num_g2_y_x ), &
  g2_y_y(  natom_y_total,  num_g2_y_y ), &
  g5_y_xx( natom_yx_total, num_g5_y_xx ), &
  g5_y_xy( natom_yx_total, num_g5_y_xy ), &
  g5_y_yy( natom_y_total,  num_g5_y_yy )
double precision :: &
  g2_x_x_maxmin(  num_g2_x_x, 2 ), &
  g2_x_y_maxmin(  num_g2_x_y, 2 ), &
  g5_x_xx_maxmin( num_g5_x_xx, 2 ), &
  g5_x_xy_maxmin( num_g5_x_xy, 2 ), &
  g5_x_yy_maxmin( num_g5_x_yy, 2 )
double precision :: &
  g2_y_x_maxmin(  num_g2_y_x, 2 ), &
  g2_y_y_maxmin(  num_g2_y_y, 2 ), &
  g5_y_xx_maxmin( num_g5_y_xx, 2 ), &
  g5_y_xy_maxmin( num_g5_y_xy, 2 ), &
  g5_y_yy_maxmin( num_g5_y_yy, 2 )
double precision,allocatable,dimension(:,:,:,:) :: g2_deriv, g5_deriv
double precision :: &
  g2_deriv_x_x(  natom_x_sum,  3, num_g2_x_x ), &
  g2_deriv_x_y(  natom_xy_sum, 3, num_g2_x_y ), &
  g5_deriv_x_xx( natom_x_sum,  3, num_g5_x_xx ), &
  g5_deriv_x_xy( natom_xy_sum, 3, num_g5_x_xy ), &
  g5_deriv_x_yy( natom_xy_sum, 3, num_g5_x_yy )
double precision :: &
  g2_deriv_y_x(  natom_yx_sum, 3, num_g2_y_x ), &
  g2_deriv_y_y(  natom_y_sum,  3, num_g2_y_y ), &
  g5_deriv_y_xx( natom_yx_sum, 3, num_g5_y_xx ), &
  g5_deriv_y_xy( natom_yx_sum, 3, num_g5_y_xy ), &
  g5_deriv_y_yy( natom_y_sum,  3, num_g5_y_yy )
  

!X
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

  allocate( g2( natom_x, num_g2_x_y ) )
  read(11) g2
  do j = 1, num_g2_x_y
    do i = 1, natom_x
      g2_x_y( i + counter_natom_xy, j ) = ( 2.0d0*( g2(i,j) - g2_x_y_maxmin(j,2) )/ &
                                          ( g2_x_y_maxmin(j,1) - g2_x_y_maxmin(j,2) ) ) - 1.0d0
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

  allocate( g5( natom_x, num_g5_x_xy ) )
  read(11) g5
  do j = 1, num_g5_x_xy
    do i = 1, natom_x
      g5_x_xy( i + counter_natom_xy, j ) = ( 2.0d0*( g5(i,j) - g5_x_xy_maxmin(j,2) )/ &
                                           ( g5_x_xy_maxmin(j,1) - g5_x_xy_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  allocate( g5( natom_x, num_g5_x_yy ) )
  read(11) g5
  do j = 1, num_g5_x_yy
    do i = 1, natom_x
      g5_x_yy( i + counter_natom_xy, j ) = ( 2.0d0*( g5(i,j) - g5_x_yy_maxmin(j,2) )/ &
                                           ( g5_x_yy_maxmin(j,1) - g5_x_yy_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  counter_natom_x  = counter_natom_x  + natom_x
  counter_natom_xy = counter_natom_xy + natom_x


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

  allocate( g2_deriv( natom_x, natom, 3, num_g2_x_y ) )
  read(12) g2_deriv
  do k = 1, num_g2_x_y
    do j = 1, natom
      do i = 1, natom_x
        m = i + natom_x * ( j - 1 )
        g2_deriv_x_y( counter_deriv_xy + m,1,k) = 2.0d0*g2_deriv(i,j,1,k)/ &
                                                 ( g2_x_y_maxmin(k,1) - g2_x_y_maxmin(k,2) )
        g2_deriv_x_y( counter_deriv_xy + m,2,k) = 2.0d0*g2_deriv(i,j,2,k)/ &
                                                 ( g2_x_y_maxmin(k,1) - g2_x_y_maxmin(k,2) )
        g2_deriv_x_y( counter_deriv_xy + m,3,k) = 2.0d0*g2_deriv(i,j,3,k)/ &
                                                 ( g2_x_y_maxmin(k,1) - g2_x_y_maxmin(k,2) )
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
        g5_deriv_x_xx( counter_deriv_x + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                 ( g5_x_xx_maxmin(k,1) - g5_x_xx_maxmin(k,2) )
        g5_deriv_x_xx( counter_deriv_x + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                 ( g5_x_xx_maxmin(k,1) - g5_x_xx_maxmin(k,2) )
        g5_deriv_x_xx( counter_deriv_x + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                 ( g5_x_xx_maxmin(k,1) - g5_x_xx_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  allocate( g5_deriv( natom_x, natom, 3, num_g5_x_xy ) )
  read(12) g5_deriv
  do k = 1, num_g5_x_xy
    do j = 1, natom
      do i = 1, natom_x
        m = i + natom_x * ( j - 1 )
        g5_deriv_x_xy( counter_deriv_xy + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                 ( g5_x_xy_maxmin(k,1) - g5_x_xy_maxmin(k,2) )
        g5_deriv_x_xy( counter_deriv_xy + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                 ( g5_x_xy_maxmin(k,1) - g5_x_xy_maxmin(k,2) )
        g5_deriv_x_xy( counter_deriv_xy + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                 ( g5_x_xy_maxmin(k,1) - g5_x_xy_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  allocate( g5_deriv( natom_x, natom, 3, num_g5_x_yy ) )
  read(12) g5_deriv
  do k = 1, num_g5_x_yy
    do j = 1, natom
      do i = 1, natom_x
        m = i + natom_x * ( j - 1 )
        g5_deriv_x_yy( counter_deriv_xy + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                 ( g5_x_yy_maxmin(k,1) - g5_x_yy_maxmin(k,2) )
        g5_deriv_x_yy( counter_deriv_xy + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                 ( g5_x_yy_maxmin(k,1) - g5_x_yy_maxmin(k,2) )
        g5_deriv_x_yy( counter_deriv_xy + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                 ( g5_x_yy_maxmin(k,1) - g5_x_yy_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  counter_deriv_x  = counter_deriv_x  + natom_x * natom
  counter_deriv_xy = counter_deriv_xy + natom_x * natom


!Y
  num_atom_y( imd + counter_y ) = natom_y

  ! SF
  allocate( g2( natom_y, num_g2_y_x ) )
  read(11) g2
  do j = 1, num_g2_y_x
    do i = 1, natom_y
      g2_y_x( i + counter_natom_yx, j ) = ( 2.0d0*( g2(i,j) - g2_y_x_maxmin(j,2) )/ &
                                          ( g2_y_x_maxmin(j,1) - g2_y_x_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g2 )

  allocate( g2( natom_y, num_g2_y_y ) )
  read(11) g2
  do j = 1, num_g2_y_y
    do i = 1, natom_y
      g2_y_y( i + counter_natom_y, j ) = ( 2.0d0*( g2(i,j) - g2_y_y_maxmin(j,2) )/ &
                                         ( g2_y_y_maxmin(j,1) - g2_y_y_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g2 )

  allocate( g5( natom_y, num_g5_y_xx ) )
  read(11) g5
  do j = 1, num_g5_y_xx
    do i = 1, natom_y
      g5_y_xx( i + counter_natom_yx, j ) = ( 2.0d0*( g5(i,j) - g5_y_xx_maxmin(j,2) )/ &
                                           ( g5_y_xx_maxmin(j,1) - g5_y_xx_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  allocate( g5( natom_y, num_g5_y_xy ) )
  read(11) g5
  do j = 1, num_g5_y_xy
    do i = 1, natom_y
      g5_y_xy( i + counter_natom_yx, j ) = ( 2.0d0*( g5(i,j) - g5_y_xy_maxmin(j,2) )/ &
                                           ( g5_y_xy_maxmin(j,1) - g5_y_xy_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  allocate( g5( natom_y, num_g5_y_yy ) )
  read(11) g5
  do j = 1, num_g5_y_yy
    do i = 1, natom_y
      g5_y_yy( i + counter_natom_y, j ) = ( 2.0d0*( g5(i,j) - g5_y_yy_maxmin(j,2) )/ &
                                          ( g5_y_yy_maxmin(j,1) - g5_y_yy_maxmin(j,2) ) ) - 1.0d0
    enddo
  enddo
  deallocate( g5 )

  counter_natom_y  = counter_natom_y  + natom_y
  counter_natom_yx = counter_natom_yx + natom_y


  ! SF_deriv
  allocate( g2_deriv( natom_y, natom, 3, num_g2_y_x ) )
  read(12) g2_deriv
  do k = 1, num_g2_y_x
    do j = 1, natom
      do i = 1, natom_y
        m = i + natom_y * ( j - 1 )
        g2_deriv_y_x( counter_deriv_yx + m,1,k) = 2.0d0*g2_deriv(i,j,1,k)/ &
                                                 ( g2_y_x_maxmin(k,1) - g2_y_x_maxmin(k,2) )
        g2_deriv_y_x( counter_deriv_yx + m,2,k) = 2.0d0*g2_deriv(i,j,2,k)/ &
                                                 ( g2_y_x_maxmin(k,1) - g2_y_x_maxmin(k,2) )
        g2_deriv_y_x( counter_deriv_yx + m,3,k) = 2.0d0*g2_deriv(i,j,3,k)/ &
                                                 ( g2_y_x_maxmin(k,1) - g2_y_x_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g2_deriv )

  allocate( g2_deriv( natom_y, natom, 3, num_g2_y_y ) )
  read(12) g2_deriv
  do k = 1, num_g2_y_y
    do j = 1, natom
      do i = 1, natom_y
        m = i + natom_y * ( j - 1 )
        g2_deriv_y_y( counter_deriv_y + m,1,k) = 2.0d0*g2_deriv(i,j,1,k)/ &
                                                 ( g2_y_y_maxmin(k,1) - g2_y_y_maxmin(k,2) )
        g2_deriv_y_y( counter_deriv_y + m,2,k) = 2.0d0*g2_deriv(i,j,2,k)/ &
                                                 ( g2_y_y_maxmin(k,1) - g2_y_y_maxmin(k,2) )
        g2_deriv_y_y( counter_deriv_y + m,3,k) = 2.0d0*g2_deriv(i,j,3,k)/ &
                                                 ( g2_y_y_maxmin(k,1) - g2_y_y_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g2_deriv )

  allocate( g5_deriv( natom_y, natom, 3, num_g5_y_xx ) )
  read(12) g5_deriv
  do k = 1, num_g5_y_xx
    do j = 1, natom
      do i = 1, natom_y
        m = i + natom_y * ( j - 1 )
        g5_deriv_y_xx( counter_deriv_yx + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                 ( g5_y_xx_maxmin(k,1) - g5_y_xx_maxmin(k,2) )
        g5_deriv_y_xx( counter_deriv_yx + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                 ( g5_y_xx_maxmin(k,1) - g5_y_xx_maxmin(k,2) )
        g5_deriv_y_xx( counter_deriv_yx + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                 ( g5_y_xx_maxmin(k,1) - g5_y_xx_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  allocate( g5_deriv( natom_y, natom, 3, num_g5_y_xy ) )
  read(12) g5_deriv
  do k = 1, num_g5_y_xy
    do j = 1, natom
      do i = 1, natom_y
        m = i + natom_y * ( j - 1 )
        g5_deriv_y_xy( counter_deriv_yx + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                 ( g5_y_xy_maxmin(k,1) - g5_y_xy_maxmin(k,2) )
        g5_deriv_y_xy( counter_deriv_yx + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                 ( g5_y_xy_maxmin(k,1) - g5_y_xy_maxmin(k,2) )
        g5_deriv_y_xy( counter_deriv_yx + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                 ( g5_y_xy_maxmin(k,1) - g5_y_xy_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  allocate( g5_deriv( natom_y, natom, 3, num_g5_y_yy ) )
  read(12) g5_deriv
  do k = 1, num_g5_y_yy
    do j = 1, natom
      do i = 1, natom_y
        m = i + natom_y * ( j - 1 )
        g5_deriv_y_yy( counter_deriv_y + m,1,k) = 2.0d0*g5_deriv(i,j,1,k)/ &
                                                 ( g5_y_yy_maxmin(k,1) - g5_y_yy_maxmin(k,2) )
        g5_deriv_y_yy( counter_deriv_y + m,2,k) = 2.0d0*g5_deriv(i,j,2,k)/ &
                                                 ( g5_y_yy_maxmin(k,1) - g5_y_yy_maxmin(k,2) )
        g5_deriv_y_yy( counter_deriv_y + m,3,k) = 2.0d0*g5_deriv(i,j,3,k)/ &
                                                 ( g5_y_yy_maxmin(k,1) - g5_y_yy_maxmin(k,2) )
      enddo
    enddo
  enddo
  deallocate( g5_deriv )

  counter_deriv_y  = counter_deriv_y  + natom_y * natom
  counter_deriv_yx = counter_deriv_yx + natom_y * natom


end subroutine read_data_2
