! Grid module to interface with python
! Compile with "f2py -c -m fgrid grid.f90"
! Import and use in python with:
! >>>import fgrid
! >>>fgrid.grid.{function name}(arguments)
module grid

    implicit none

    private
    public volume

    contains

        ! Calculate the volume of grid boxes give the dimensions
        ! Computed using a spherical integral
        subroutine volume(rho,bounds,theta,phi,nz,ny,nx,boxes)
            !F2PY integer, intent(hide):: nz
            !F2PY integer, intent(hide):: ny
            !F2PY integer, intent(hide):: nx
            !F2PY real, intent(in):: rho(nz,ny,nx)
            !F2PY real, intent(in):: bounds(nz,ny,nx,2)
            !F2PY real, intent(in):: theta(nx)
            !F2PY real, intent(in):: phi(ny)
            !F2PY real, intent(out):: boxes(nz,ny,nx)

            real, intent(out) :: boxes(nz,ny,nx)

            !-Input Variables
            ! Grid dimensions
            integer, intent(in) :: nz, ny, nx

            ! Height (Earth Radius + Altitude)
            real, intent(in) :: rho(nz,ny,nx)

            ! Height at edges of each gridbox
            real, intent(in) :: bounds(nz,ny,nx,2)

            ! Longitude
            real, intent(in) :: theta(nx)

            ! 90 - Latitude
            real, intent(in) :: phi(ny)

            !-Local Variables
            real :: dtheta, dphi
            integer :: i, j, k

            ! Calculate lat/lon grid spacing (assumed to be even)
            dtheta = theta(2) - theta(1)
            dphi = phi(2) - phi(1)

            ! Loop over each gridbox
            do k=1,nz
                do j=1,ny
                    do i=1,nx
                        ! Calculate the volume in the gridbox as a spherical integral
                        boxes(k,j,i) = &
                                abs(rho(k,j,i)**2 * sin(phi(j)) * (bounds(k,j,i,2) - bounds(k,j,i,1)) * dtheta * dphi)
                    end do
                end do
            end do

        end subroutine volume
end module grid
