! Module for disk model generation/manipulation.
module disk
    use grid
    implicit none
    
contains

    subroutine temperature_gaussian(T, z)
        real, dimension(nz), intent(inout) :: T
        real, dimension(mz), intent(in) :: z
        real, dimension(nz) :: Tgauss
        real :: eps,isoTemp,sigma,midTemp,floorTemp,switchTemp
        integer :: i
        
        eps = 0.2
        isoTemp = 64857.04298198191
        sigma = 2.54d14
        midTemp = 61341
        floorTemp = 3
        switchTemp = (1+eps) * isoTemp ! must be larger than mdiTemp
        
        ! Caluclate the gaussian profile
        do i=1, nz
            Tgauss(i) = midTemp * exp(-0.5 * z(n1+i-1)**2/sigma**2)
        enddo
        
        ! Choose between the gaussian profile, the floor, or the current value.
        do i = 1, nz
            if ((T(i) < switchTemp) .and. (Tgauss(i) > floorTemp)) then
                T(i) = Tgauss(i)
            elseif ((T(i) < switchTemp) .and. (Tgauss(i) < floorTemp)) then
                T(i) = floorTemp
            endif
        end do

    end subroutine temperature_gaussian
    
end module disk