! Interpolation module to interface with python
! Compile with "f2py -c -m finterpolate interpolate.f90"
! Import and use in python with:
! >>>import finterpolate
! >>>finterpolate.interpolate.{function name}(arguments)
module interpolate
    implicit none

    public to_level

    contains

        ! Interpolate the input variable to levels defined by coord_out.
        ! Interpolate in the vertical to the highest point in coord_in with the desired value from coord out
        subroutine to_level(variable, coord_in, coord_out, type, order, np, nz, ny, nx, output, mask)
            !F2PY integer, intent(hide) :: np
            !F2PY integer, intent(hide) :: nz
            !F2PY integer, intent(hide) :: ny
            !F2PY integer, intent(hide) :: nx
            !F2PY real   , intent(in)   :: variable(nz,ny,nx)
            !F2PY real   , intent(in)   :: coord_in(nz,ny,nx)
            !F2PY real   , intent(in)   :: coord_out(np,ny,nx)
            !F2PY integer, intent(in)   :: type
            !F2PY integer, intent(in)   :: order
            !F2PY real   , intent(out)  :: output(np,ny,nx)
            !F2PY logical, intent(out)  :: mask(np,ny,nx)

            !-Input variables-
            ! Grid Dimensions
            integer, intent(in) :: np, nz, ny, nx

            ! Input Field
            real, intent(in) :: variable(nz, ny, nx)

            ! Vertical coord on same levels as input field
            real, intent(in) :: coord_in(nz, ny, nx)

            ! Vertical coord levels to interpolate to
            real, intent(in) :: coord_out(np, ny, nx)

            ! Interpolation type
            ! 0 = Linear
            ! 1 = Log-Linear
            ! 2 = Lagrange
            integer, intent(in) :: type

            ! Order of Lagrange polynomial interpolation
            integer, intent(in) :: order

            !-Output Variables-
            ! Output field interpolated to given levels
            real, intent(out) :: output(np, ny, nx)
            logical, intent(out) :: mask(np, ny, nx)

            !-Local Variables-
            integer :: i, j, k, l
            real :: z

            !-START-
            ! Loop over zcoord_out points
            do l=1,np
                do i=1,nx
                    do j=1,ny
                        z = coord_out(l,j,i)
                        if (out_of_bounds(z, coord_in(:,j,i), nz)) then
                            ! Set the mask to True if the requested value is not within the bounds of zcoord
                            mask(l, j, i) = .True.
                        else
                            ! Set the mask to False and do the interpolation if the requested value is within the
                            ! bounds of zcoord
                            mask(l, j, i) = .False.

                            ! Find the index of the first zcoord point that has crossed the output coordinate by
                            ! searching downwards
                            k = search_downwards(z, coord_in(:,j,i), nz)

                            ! Interpolate to the output point
                            if (type == 0) then
                                output(l,j,i) = linear(variable(:,j,i), coord_in(:,j,i), z, k, nz)
                            else if (type == 1) then
                                output(l,j,i) = log_linear(variable(:,j,i), coord_in(:,j,i), z, k, nz)
                                !else if (type == 2)
                                ! output(l,j,i) = lagrange(z, coord(:,j,i), k, nz)
                            end if
                        end if
                    end do
                end do
            end do

        end subroutine to_level
        !-------------------------------------------------------------------------------

        ! Search downwards to find the index where zcoord first crosses the value z.
        ! Before using this function it must be checked that z is within the bounds of zcoord
        function search_downwards(z, zcoord, nz)
            integer :: search_downwards

            !-Input Variables-
            real, intent(in) :: z, zcoord(nz)
            integer, intent(in) :: nz

            !-Local variables-
            integer :: k
            integer :: plusminus_at_top

            !-START-
            plusminus_at_top = plusminus(z - zcoord(nz))
            k = nz - 1
            do while (plusminus_at_top == plusminus(z-zcoord(k)))
                k = k-1
            end do

            search_downwards = k
        end function search_downwards

        ! Check whether the value z is contained within zcoord
        pure function out_of_bounds(z, zcoord, nz)
            logical :: out_of_bounds

            !-Input Variables-
            real, intent(in) :: z, zcoord(nz)
            integer, intent(in) :: nz

            !-Local variables-
            real :: zmax, zmin
            integer :: k

            !-START-
            zmax = zcoord(1)
            zmin = zcoord(1)
            do k=2, nz
                zmax = max(zmax, zcoord(k))
                zmin = min(zmin, zcoord(k))
            end do

            out_of_bounds = z > zmax.or.z < zmin
        end function out_of_bounds

        ! Calculate the interpolation weight assuming linear variation across two points
        function linear(variable, coord, z, k, nz)

            real :: linear

            !-Input Variables-
            real, intent(in) :: variable(nz), coord(nz), z
            integer, intent(in) :: k, nz

            !-Local Variables-
            ! Interpolation weight
            real :: alpha

            !-START-
            alpha = (z - coord(k)) / (coord(k+1) - coord(k))
            linear = alpha * variable(k+1) + (1 - alpha) * variable(k)
        end function linear

        ! Calculate the interpolation weight assuming linear logarithmic variation across two points
        pure function log_linear(variable, coord, z, k, nz)
            real :: log_linear

            !-Input Variables-
            real, intent(in) :: variable(nz), coord(nz), z
            integer, intent(in) :: k, nz

            !-Local Variables-
            ! Interpolation weight
            real :: alpha

            !-START-
            alpha = log(z / coord(k)) / log(coord(k+1) / coord(k))
            log_linear = alpha * variable(k+1) + (1 - alpha) * variable(k)
        end function log_linear

        ! Check whether x is positive or negative
        pure function plusminus(x)
            integer :: plusminus

            !-Input Variables-
            real, intent(in) :: x

            !-START-
            if (x < 0.0) then
                plusminus = -1
            else
                plusminus = 1
            end if
        end function plusminus

end module interpolate
