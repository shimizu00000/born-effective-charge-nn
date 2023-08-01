subroutine pre_random
!This subroutine creats random number seed.
implicit none
integer::seedsize,c
integer,allocatable::seed(:)
 call random_seed(size=seedsize)
 allocate(seed(1:seedsize))
 call system_clock(count=c)
 seed=c
 call random_seed(put=seed)

return
end subroutine pre_random
