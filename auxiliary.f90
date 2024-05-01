module auxiliary

  use grid

  implicit none
  private
  
  public :: tridag,update_ghosts,der

contains
!******************************************
  subroutine tridag(a,b,c,r,u)
!
!  Solves a tridiagonal system.
!  Imported from numerical recipes.    
!
    real, dimension(:), intent(in) :: a,b,c,r
    real, dimension(:), intent(out) :: u
    real, dimension(size(b)) :: gam
    integer :: n,j
    real :: bet
!
    n=size(b)
    bet=b(1)
    if (bet==0.0) then
      print*,"ERROR tridiag stage 1: bet=b(1)) = ", bet
      stop
    endif
!
    u(1)=r(1)/bet
    do j=2,n
      gam(j)=c(j-1)/bet
      bet=b(j)-a(j)*gam(j)
      if (bet==0.0) then
        print*,"ERROR tridiag stage 2: bet = b(j)-a(j)*gam(j) =", bet
        stop
      endif
      u(j)=(r(j)-a(j)*u(j-1))/bet
    enddo
!
    do j=n-1,1,-1
      u(j)=u(j)-gam(j+1)*u(j+1)
    enddo
!
  endsubroutine tridag
!******************************************
  subroutine update_ghosts(f)
!
!  Update the ghost zones of an array
!  using constant gradient. Used on U.     
!
    real, dimension(:), intent(inout) :: f
    integer :: i
!      
    do i=1,ng
      f(n1-i)=2*f(n1) - f(n1+i)
      f(n2+i)=2*f(n2) - f(n2-i)
    enddo
!
  endsubroutine update_ghosts
!******************************************
  subroutine der(f,df)
!
!  Sixth-order first derivative.
!    
    real, dimension(mz) :: f
    real, dimension(nz) :: df
!      
    intent(in) :: f
    intent(out) :: df
!
    df=1./60.*(+ 45.0*(f(n1+1:n2+1)-f(n1-1:n2-1)) &
               -  9.0*(f(n1+2:n2+2)-f(n1-2:n2-2)) &
               +      (f(n1+3:n2+3)-f(n1-3:n2-3)))
!    
  endsubroutine der
!******************************************         
  endmodule auxiliary
