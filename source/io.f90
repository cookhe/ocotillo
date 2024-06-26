module FileIO

  use Common
  
  implicit none
  private
  
  public :: output_ascii,output_grid,output_binary,itoa,file_exists
  
  !interface operator( .f. )
  !   module procedure file_exists
  !endinterface operator( .f. )

contains
!************************************************************************************
  function file_exists(filename) result(res)
    implicit none
    character(len=*),intent(in) :: filename
    logical                     :: res

  ! Check if the file exists
    inquire( file=trim(filename), exist=res )

  endfunction file_exists
!************************************************************************************
    character (len=21) function itoa(n)
! plucked from pencil
      integer, intent(in) :: n
!
      write (itoa, '(I21)') n   ! (64 bit integer plus sign)
      itoa = adjustl(itoa)
!
    endfunction itoa
!************************************************************************************
  subroutine output_ascii(U,absorp_coeff)

    real, dimension(mz,nw) :: U
    real, dimension(nz,nw) :: absorp_coeff
    integer :: i,iw
!
    open(10,file="output/diagnostics.dat",status="replace",action='write')
    do i=1,nz
      do iw=1,nw
         write(unit=10,FMT="(2I6,2e15.5)") i,iw,U(n1+i-1,iw),absorp_coeff(i,iw)
      enddo
    enddo
    close(10)

  endsubroutine output_ascii
!************************************************************************************
  subroutine output_grid(z)

    real, dimension(mz) :: z
    integer :: i
!
    open(15,file="output/zgrid.dat",status="replace",action='write')
    do i=1,nz
      write(unit=15,FMT=*) i,z(n1+i-1)
    enddo
    close(15)
!
    open(55, file='output/zgrid.fvar',form='unformatted',status='replace',action='write')
    write(55) z
    close(55)
!    
  endsubroutine output_grid
!************************************************************************************  
  subroutine output_binary(U,absorp_coeff,iprocx,iprocy)

    real, dimension(mz,nyloc,nxloc,nw) :: U
    real, dimension(nz,nyloc,nxloc,nw) :: absorp_coeff
    integer :: iprocx,iprocy
    character(len=90) :: proc_str

    proc_str = '_xproc'//trim(itoa(iprocx))//'_yproc'//trim(itoa(iprocy))
    
    open(35, file='output/mean_intensity'//trim(proc_str)//'.fvar', &
         form='unformatted',status='replace',action='write')
    write(35) U
    close(35)

    open(45, file='output/absorption_coefficients'//trim(proc_str)//'.fvar', &
         form='unformatted',status='replace',action='write')
    write(45) absorp_coeff
    close(45)

  endsubroutine output_binary
!************************************************************************************
end module FileIO
