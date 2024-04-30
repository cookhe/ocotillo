module auxiliary
  implicit none
  private
  public :: tridag 
contains
subroutine tridag(a,b,c,r,u,err,msg)
!                                                                                                                
!  Solves a tridiagonal system.                                                                                  
!                                                                                                                
!  01-apr-03/tobi: from Numerical Recipes (p42-43).                                                              
!                                                                                                                
!  | b1 c1 0  ...            | | u1   |   | r1   |                                                               
!  | a2 b2 c2 ...            | | u2   |   | r2   |                                                               
!  | 0  a3 b3 c3             | | u3   | = | r3   |                                                               
!  |          ...            | | ...  |   | ...  |                                                               
!  |          an-1 bn-1 cn-1 | | un-1 |   | rn-1 |                                                               
!  |          0    a_n  b_n  | | un   |   | rn   |                                                               
!                                                                                                                

  real, dimension(:), intent(in) :: a,b,c,r
      real, dimension(:), intent(out) :: u
      real, dimension(size(b)) :: gam
      logical, intent(out), optional :: err
      character(len=*), intent(out), optional :: msg
      integer :: n,j
      real :: bet
!                                                                                                                
      if (present(err)) err=.false.
      n=size(b)
      bet=b(1)
      if (bet==0.0) then
        print*, ' ERROR tridiag stage 1: bet=b(1)) = ', bet
        if (present(msg)) msg = 'tridag: Error at code stage 1'
        if (present(err)) err = .true.
      endif
!                                                                                                                
      u(1)=r(1)/bet
      do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if (bet==0.0) then
          print*, 'ERROR tridiag stage 2: bet = b(j)-a(j)*gam(j) =', bet
          if (present(msg)) msg = 'tridag: Error at code stage 2'
          if (present(err)) err = .true.
          return
        endif
        u(j)=(r(j)-a(j)*u(j-1))/bet
      enddo
!                                                                                                                
      do j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
      enddo
!                                                                                                                
    endsubroutine tridag
endmodule auxiliary    
