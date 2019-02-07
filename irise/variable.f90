! UM specific variable calculation module to interface with python
module variable

    implicit none

    real, parameter :: g=9.81, pi=3.14159, a=6378100

    contains

        subroutine calc_pv(u, v, theta, rho, r_theta, r_rho, r_u, r_v, sec_v_lat, tan_v_lat, sec_theta_lat, f3, &
                delta_lambda, delta_phi, nx, ny, nz, pv)
            !f2py integer, intent(hide) :: nx, ny, nz
            !f2py real, intent(in) :: u(nz, ny, nx)
            !f2py real, intent(in) :: v(nz, ny-1, nx)
            !f2py real, intent(in) :: theta(nz, ny, nx)
            !f2py real, intent(in) :: rho(nz, ny, nx)
            !f2py real, intent(in) :: r_theta(nz, ny, nx)
            !f2py real, intent(in) :: r_rho(nz, ny, nx)
            !f2py real, intent(in) :: r_u(nz, ny, nx)
            !f2py real, intent(in) :: r_v(nz, ny-1, nx)
            !f2py real, intent(in) :: sec_v_lat(ny-1, nx)
            !f2py real, intent(in) :: tan_v_lat(ny-1, nx)
            !f2py real, intent(in) :: sec_theta_lat(ny, nx)
            !f2py real, intent(in) :: f3(ny, nx)
            !f2py real, intent(in) :: delta_lambda, delta_phi
            !f2py real, intent(out) :: pv(nz, ny, nx)

            !-Input variables-
            ! Grid dimensions (ny-1 for v points having extra outer gridpoint)
            integer, intent(in) :: nx, ny, nz

            ! Wind field on staggered grid points
            real, intent(in) :: u(nz, ny, nx), v(nz, ny-1, nx)

            ! Potential temperature
            real, intent(in) :: theta(nz, ny, nx)

            ! Prognostic density (density * r^2)
            real, intent(in) :: rho(nz, ny, nx)

            ! Radius (a + z) at different grid points
            real, intent(in) :: r_theta(nz+1, ny, nx), r_rho(nz, ny, nx),                &
            r_u(nz, ny, nx), r_v(nz, ny-1, nx)

            ! Trigonometric functions at different grid points
            real, intent(in) :: sec_v_lat(ny-1, nx), tan_v_lat(ny-1, nx),                &
            sec_theta_lat(ny,nx)

            ! Vertical coriolis parameter at v-points
            real, intent(in) :: f3(ny, nx)

            ! Grid spacing
            real, intent(in) :: delta_lambda, delta_phi

            ! Potential Vorticity (output)
            real, intent(out) :: pv(nz, ny, nx)

            !-Local variables-
            integer :: i,j,k

            real :: recip_delta_lambda, recip_delta_phi, weight1(ny,nx), weight2(ny,nx)

            real :: r_pv(nz,ny,nx), r_pv_theta(nz,ny,nx)
            real :: dtheta_dr(ny,nx), du_dr(ny,nx) , dv_dr(ny,nx)
            real :: dtheta_dx(ny,nx), dtheta_dx_ave(nz,ny,nx)
            real :: du_dr_ave(ny,nx), dv_dr_ave(ny,nx)
            real :: dtheta_dy(ny,nx), dtheta_dy_ave(nz,ny,nx), dtheta_dr_ave(ny,nx)
            real :: x_term(ny,nx), y_term(ny,nx), z_term(ny,nx)
            real :: density(ny,nx), density_ave(ny,nx)

            !-START-
            ! Section 0. Initialisation
            recip_delta_lambda = 1.0/ delta_lambda
            recip_delta_phi = 1.0/ delta_phi


            ! Section 1. Calculate horizontal part of terms
            do k=1,nz
            ! calculate dtheta/dx and dtheta/dy
            dtheta_dx = d_dx(theta(k,:,:), recip_delta_lambda, sec_theta_lat, r_theta(k,:,:), ny, nx)
            dtheta_dx_ave(k, :, :) = average_y(dtheta_dx, ny, nx)

            dtheta_dy = d_dy(theta(k,:,:), recip_delta_phi, r_theta(k,:,:), ny, nx)
            dtheta_dy_ave(k, :, :) = average_x(dtheta_dy, ny, nx)

            ! Calculate altitude on theta levels
            r_pv(k,1:ny-1,:) = 0.5 * (r_u(k,2:ny,:) + r_u(k,1:ny-1,:))
            r_pv_theta(k,1:ny-1,:) = centre(r_theta(k,:,:), ny, nx)

            ! end loop over levels
            end do

            ! loop over model levels for rest of code.
            do k=1,nz
                ! Section 1.1 Calculate x_term and y_term.
                if (k == 1) then
                    x_term(:,:) = -dtheta_dx_ave(k,:,:)
                    y_term(:,:) = dtheta_dy_ave(k,:,:)

                else
                    weight1 = (r_pv_theta(k,:,:) - r_pv(k,:,:) ) / (r_pv_theta(k,:,:) - r_pv_theta(k-1,:,:))
                    weight2 = 1.0 - weight1

                    x_term = -(weight2 * dtheta_dx_ave(k-1,:,:) + weight1 * dtheta_dx_ave(k,:,:))
                    y_term = (weight2 * dtheta_dy_ave(k-1,:,:) + weight1 * dtheta_dy_ave(k,:,:))
                end if


                ! Section 1.2 Calculate z_term.
                do j=1,ny-1
                    do i=1,nx-1
                        z_term(j,i) = f3(j,i) + ((v(k,j,i+1)-v(k,j,i)) * recip_delta_lambda * sec_v_lat(j,i) + &
                                0.5 * (u(k,j+1,i) + u(k,j,i)) * tan_v_lat(j,i) - &
                                (u(k,j+1,i) - u(k,j,i)) * recip_delta_phi) / r_pv(k,j,i)
                    end do
                end do

                ! Section 2. Multiply horizontal terms by vertical terms and form full pv.
                if (k==1) then
                    ! Since theta gradient is un-defined use gradient from levels above.
                    dtheta_dr = (theta(k+1,:,:) - theta(k,:,:)) / (r_theta(k+1,:,:) - r_theta(k,:,:))
                else if (k/=2) then
                    ! no code required at level 2 since the value is the same as level 1.
                    ! calculate theta gradient.
                    dtheta_dr = (theta(k,:,:) - theta(k-1,:,:)) / (r_theta(k,:,:) - r_theta(k-1,:,:))
                end if
                ! Now average to required location.
                dtheta_dr_ave = centre(dtheta_dr, ny, nx)

                if (k==1) then
                    dv_dr(1:ny-1,:) = (v(k+1,:,:) - v(k,:,:)) / (r_v(k+1,:,:) - r_v(k,:,:))
                    du_dr = (u(k+1,:,:) - u(k,:,:)) / (r_u(k+1,:,:) - r_u(k,:,:))
                else if (k==nz) then
                    dv_dr(1:ny-1,:) = (v(k,:,:) - v(k-1,:,:)) / (r_v(k,:,:) - r_v(k-1,:,:))
                    du_dr = (u(k,:,:) - u(k-1,:,:)) / (r_u(k,:,:) - r_u(k-1,:,:))
                else
                    dv_dr(1:ny-1,:) = (v(k+1,:,:) - v(k-1,:,:)) / (r_v(k+1,:,:) - r_v(k-1,:,:))
                    du_dr = (u(k+1,:,:) - u(k-1,:,:)) / (r_u(k+1,:,:) - r_u(k-1,:,:))
                end if

                ! now average quantities
                du_dr_ave = average_y(du_dr, ny, nx)
                dv_dr_ave = average_x(dv_dr, ny, nx)

                ! convert rho to true density by removing factor of r squared.
                density = rho(k,:,:) / (r_rho(k,:,:) * r_rho(k,:,:))

                density_ave = centre(density, ny, nx)

                ! Calculate full PV.
                pv(k,:,:) = (x_term * dv_dr_ave + y_term * du_dr_ave + z_term * dtheta_dr_ave) / density_ave

                ! end loop over model levels
            end do
        end subroutine calc_PV


        subroutine hydrostatic_mass(P_surf, dlambda, phi, sigh, nx, ny, nz, mass)
            !f2py integer, intent(hide) :: nx, ny, nz
            !f2py real, intent(in) :: P_surf(ny, nx)
            !f2py real, intent(in) :: dlambda
            !f2py real, intent(in) :: phi(ny)
            !f2py real, intent(in) :: sigh(0:nz)
            !f2py real, intent(out) :: mass(nz, ny, nx)

            !-Input variables-
            ! Grid dimensions
            integer, intent(in) :: nx, ny, nz

            ! Surface pressure
            real, intent(in) :: P_surf(ny, nx)

            ! Longitude spacing (radians)
            real, intent(in) :: dlambda

            ! Latitude (radians)
            real, intent(in) :: phi(ny)

            ! Half sigma levels
            real, intent(in) :: sigh(0:nz)

            ! Mass in each gridbox
            real, intent(out) :: mass(nz, ny, nx)

            !-Local variables-
            integer :: i, j, k
            real :: dphi(ny)
            real :: rho_dz

            !-START-
            ! Calculate latitude spacing
            dphi(1) = (phi(2) + pi / 2.0) / 2.0
            do j = 2, ny - 1
                dphi(j) = (phi(j+1) - phi(j-1)) / 2.0
            end do
            dphi(ny) = (pi/2.0 - phi(ny-1)) / 2.0

            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        rho_dz = (sigh(k-1) - sigh(k)) * P_surf(j, i) / g
                        mass(k,j,i) = rho_dz * a**2 * cos(phi(j)) * dlambda * dphi(j)
                    end do
                end do
            end do
        end subroutine hydrostatic_mass


        ! da/dx = da/d(lambda) * (1/(r*cos(phi)))
        pure function d_dx(a, recip_dlambda, sec_lat, r, ny, nx)
            !-Input Variables-
            integer, intent(in) :: ny, nx
            ! Variable to differentiate
            real, intent(in) :: a(ny,nx)
            ! Grid spacing
            real, intent(in) :: recip_dlambda
            ! Grid spacing factor with latitude
            real, intent(in) :: sec_lat(ny,nx)
            ! Radius at grid points
            real, intent(in) :: r(ny,nx)

            !-Output Variables-
            real :: d_dx(ny, nx)

            !-START-
            d_dx(:, 1:nx-1) = (a(:, 2:nx) - a(:, 1:nx-1)) * recip_dlambda * sec_lat(:, 1:nx-1) * 2.0 / &
                    (r(:, 2:nx) + r(:, 1:nx-1))
        end function d_dx


        ! da/dy = da/d(phi) * (1/r)
        pure function d_dy(a, recip_dphi, r, ny, nx)
            !-Input Variables-
            integer, intent(in) :: ny, nx
            ! Variable to differentiate
            real, intent(in) :: a(ny,nx)
            ! Grid spacing
            real, intent(in) :: recip_dphi
            ! Radius at grid points
            real, intent(in) :: r(ny,nx)

            !-Output Variables-
            real :: d_dy(ny, nx)

            !-START-
            d_dy(1:ny-1, :) = (a(2:ny, :) - a(1:ny-1, :)) * recip_dphi * 2.0 / (r(2:ny, :) + r(1:ny-1, :))

        end function d_dy



        pure function average_x(a, ny, nx)
            !-Input Variables-
            integer, intent(in) :: ny, nx
            ! Variable to average
            real, intent(in) :: a(ny,nx)

            !-Output Variables-
            real :: average_x(ny, nx)

            !-START-
            average_x(:, 1:nx-1) = 0.5*(a(:, 2:nx) - a(:, 1:nx-1))

        end function average_x


        pure function average_y(a, ny, nx)
            !-Input Variables-
            integer, intent(in) :: ny, nx
            ! Variable to average
            real, intent(in) :: a(ny,nx)

            !-Output Variables-
            real :: average_y(ny, nx)

            !-START-
            average_y(1:ny-1, :) = 0.5*(a(2:ny, :) - a(1:ny-1, :))

        end function average_y


        pure function centre(x, ny, nx)
            ! Calculate values at the centre of grid boxes

            !-Input Variables-
            integer, intent(in) :: ny, nx
            ! Variable to differentiate
            real, intent(in) :: x(ny, nx)

            !-Output Variables-
            real :: centre(ny, nx)

            !-START-
            centre(1:ny-1, 1:nx-1) = 0.25 * (x(1:ny-1, 1:nx-1) + x(2:ny, 1:nx-1) + x(1:ny-1, 2:nx) + x(2:ny,2:nx))
        end function centre
end module variable

