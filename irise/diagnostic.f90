! Diagnostics module to interface with python
module diagnostic
    use interpolate, only: linear, search_downwards

    implicit none

    private
    public dynamical_tropopause, dynamical_tropopause_sc, tropopause, &
    mask_trop, thermal_tropopause, thermal_tropopause_sc

    contains

        ! Find the height of the dynamical tropopause for each grid point using a PV based definition
        subroutine dynamical_tropopause(pv, q, z, pvtrop, qtrop, nz, ny, nx, ztrop, fold_top, fold_bottom)
            !F2PY integer, intent(hide) :: nz
            !F2PY integer, intent(hide) :: ny
            !F2PY integer, intent(hide) :: nx
            !F2PY real, intent(in) :: pv(nz,ny,nx)
            !F2PY real, intent(in) :: q(nz,ny,nx)
            !F2PY real, intent(in) :: z(nz,ny,nx)
            !F2PY real, intent(in) :: pvtrop
            !F2PY real, intent(in) :: qtrop
            !F2PY real, intent(out) :: ztrop(ny,nx)
            !F2PY real, intent(out) :: fold_top(ny,nx)
            !F2PY real, intent(out) :: fold_bottom(ny,nx)

            !-Input variables-
            ! Grid Dimensions
            integer, intent(in) :: nz, ny, nx

            ! Potential vorticity
            real, intent(in) :: pv(nz,ny,nx)

            ! Specific humidity
            real, intent(in) :: q(nz,ny,nx)

            ! Height at each gridpoint
            real, intent(in) :: z(nz,ny,nx)

            ! Tropopause pv value. Typically 1.5-3.5 PVU
            real, intent(in) :: pvtrop

            ! Stratospheric humidity threshold
            real, intent(in) :: qtrop

            !-Output Variables-
            ! The altitude of the tropopause
            real, intent(out) :: ztrop(ny,nx)
            real, intent(out) :: fold_top(ny,nx)
            real, intent(out) :: fold_bottom(ny,nx)

            !-Local variables-
            integer :: i,j,k

            !-START-
            ztrop = 0.0
            fold_top = 0.0
            fold_bottom = 0.0
            ! Loop over grid columns
            do j=1,ny
                do i=1,nx
                    k = nz
                    call dynamical_tropopause_sc(pv(:,j,i), q(:,j,i), z(:,j,i), pvtrop, qtrop, ztrop(j,i), k)
                    if (k > 1) then
                        call dynamical_tropopause_sc(pv(:,j,i), q(:,j,i), z(:,j,i), pvtrop, qtrop, fold_top(j,i), k)
                    end if
                    if (k > 1) then
                        call dynamical_tropopause_sc(pv(:,j,i), q(:,j,i), z(:,j,i), pvtrop, qtrop, fold_bottom(j,i), k)
                    end if
                end do
            end do
        end subroutine dynamical_tropopause

        ! Find the height of the dynamical tropopause
        subroutine dynamical_tropopause_sc(pv, q, z, pvtrop, qtrop, ztrop, k)
            !-Input variables-
            ! Potential vorticity profile
            real, intent(in) :: pv(:)

            ! Humidity profile
            real, intent(in) :: q(:)

            ! Height of each point in the profile
            real, intent(in) :: z(:)

            ! Tropopause pv value. Typically 1.5-3.5 PVU
            real, intent(in) :: pvtrop

            ! Stratospheric humidity threshold
            real, intent(in) :: qtrop

            !-Output Variables-
            ! The altitude of the tropopause
            real, intent(out) :: ztrop

            ! Index of position in array
            integer, intent(inout) :: k

            !-Local variables-
            ! Input array size
            integer :: nz
            logical :: found

            !-START-
            found = .true.
            nz = size(pv)
            ztrop = 0.0
            ! Search downward for where PV crosses threshold
            do while (k > 1.and.found)
                k = k-1
                ! Check whether the pv crosses the threshold
                if ((pv(k+1) > pvtrop).neqv.(pv(k) > pvtrop)) then
                    ! Check the humidity threshold
                    if (q(k) < qtrop) then
                        ! Interpolate to the exact height where the threshold is crossed
                        ztrop = linear(z, pv, pvtrop, k, nz)

                        ! Exit while loop
                        found = .false.
                    end if
                end if
            end do
        end subroutine dynamical_tropopause_sc

        ! Find the height of the thermal tropopause for each grid point using the WMO temperature lapse rate definition
        subroutine thermal_tropopause(T, z, zmin, lapse_rate_trop, dz, nz, ny, nx, ztrop)
            !F2PY integer, intent(hide) :: nz
            !F2PY integer, intent(hide) :: ny
            !F2PY integer, intent(hide) :: nx
            !F2PY real, intent(in) :: T(nz,ny,nx)
            !F2PY real, intent(in) :: z(nz,ny,nx)
            !F2PY real, intent(in) :: lapse_rate_trop
            !F2PY real, intent(in) :: dz
            !F2PY real, intent(out) :: ztrop(ny,nx)

            !-Input variables-
            ! Grid Dimensions
            integer, intent(in) :: nz, ny, nx

            ! Temperature
            real, intent(in) :: T(nz,ny,nx)

            ! Height at each gridpoint
            real, intent(in) :: z(nz,ny,nx)

            ! Minimum starting height to diagnose tropopause to ignore boundary layer
            real, intent(in) :: zmin

            ! Stratospheric lapse rate. 0.002 K m-1 for WMO definition
            real, intent(in) :: lapse_rate_trop

            ! Height above tropopause for average lapse rate required to be stratospheric
            real, intent(in) :: dz

            !-Output Variables-
            ! The altitude of the tropopause
            real, intent(out) :: ztrop(ny,nx)

            !-Local variables-
            integer :: i, j

            ! Starting height index
            integer :: k

            !-START-
            ! Loop over grid columns
            do j=1,ny
                do i=1,nx
                    k = 1
                    do while (z(k,j,i) < zmin)
                        k = k+1
                    end do
                    ztrop(j,i) = thermal_tropopause_sc(T(k:,j,i), z(k:,j,i), lapse_rate_trop, dz)
                end do
            end do
        end subroutine thermal_tropopause

        ! Find the height of the thermal tropopause using the WMO temperature lapse
        ! rate definition
        function thermal_tropopause_sc(T, z, lapse_rate_trop, dz)
            !-Input variables-
            ! Temperature profile
            real, intent(in) :: T(:)

            ! Height of each temperature measurement
            real, intent(in) :: z(:)

            ! Stratospheric lapse rate. 0.002 K m-1 for WMO definition
            real, intent(in) :: lapse_rate_trop

            ! Height above tropopause for average lapse rate required to be
            ! stratospheric
            real, intent(in) :: dz

            !-Output Variables-
            ! The altitude of the tropopause
            real :: thermal_tropopause_sc

            !-Local variables-
            ! Input array size
            integer :: nz

            ! Counter for vertical ascent
            integer :: k

            ! Index for checking average lapse rate above tropopause
            integer :: k_above

            ! Temporary lapse rate and altitude at which it is calculated
            real :: lapse_rate(2), z_lr(2)

            ! Temperature at the diagnosed tropopause and a threshold distance
            ! above. Used to calculate the average lapse rate
            real :: T_trop, T_above, z_above, av_lapse_rate

            ! Output tropopause height
            real :: ztrop

            !-START-
            nz = size(T)
            lapse_rate(2) = 0.0
            z_lr(2) = 0.0
            ! Search upward for where lapse rate drops below threshold
            k = 0
            do while (k < nz)
                k = k+1
                ! Store the previous lapse rate and altitude
                lapse_rate(1) = lapse_rate(2)
                z_lr(1) = z_lr(2)

                ! Calculate the lapse rate between two grid boxes
                z_lr(2) = 0.5 * (z(k+1) + z(k))
                lapse_rate(2) = (T(k+1) - T(k)) / (z(k+1) - z(k))

                ! Check whether the lapse crosses the threshold
                if (lapse_rate(2) > lapse_rate_trop.and.lapse_rate(1) < lapse_rate_trop) then
                    ! Interpolate to the exact height where the threshold is crossed
                    ztrop = linear(z_lr, lapse_rate, lapse_rate_trop, 1, 2)

                    ! Interpolate to temperature at the tropopause
                    T_trop = linear(T, z, ztrop, k, nz)

                    ! Check that the average lapse rate is stratospheric above the tropopause
                    ! (next 2km for WMO)
                    z_above = ztrop + dz

                    ! If the top of the domain is within 2km use the average lapse rate to the
                    ! top of the domain
                    if (z(nz) <= z_above) then
                        z_above = z(nz)
                        T_above = T(nz)
                    else
                        k_above = search_downwards(z_above, z, nz)
                        ! Interpolate to temperature above the tropopause
                        T_above = linear(T, z, z_above, k_above, nz)
                    end if

                    ! Calculate the average lapse rate over the height
                    av_lapse_rate = (T_above - T_trop) / (z_above - ztrop)

                    if (av_lapse_rate > lapse_rate_trop) then
                        ! If tropopause is found then end while loop
                        k = nz+1
                    else
                    ! Reset tropopause height and continue
                    ztrop = 0.0
                    end if
                end if
            end do
            thermal_tropopause_sc = ztrop
        end function thermal_tropopause_sc

        subroutine mask_trop(pv, z, pvtrop, height, nz, ny, nx, k, j, i, mask)
            !F2PY      INTEGER, INTENT(HIDE)    :: nz
            !F2PY      INTEGER, INTENT(HIDE)    :: ny
            !F2PY      INTEGER, INTENT(HIDE)    :: nx
            !F2PY      INTEGER, INTENT(IN)      :: k
            !F2PY      INTEGER, INTENT(IN)      :: j
            !F2PY      INTEGER, INTENT(IN)      :: i
            !F2PY      REAL, INTENT(IN)         :: pv(nz, ny, nx)
            !F2PY      REAL, INTENT(IN)         :: z(nz, ny, nx)
            !F2PY      REAL, INTENT(IN)         :: pvtrop
            !F2PY      REAL, INTENT(IN)         :: height
            !F2PY      REAL, INTENT(INOUT)      :: mask(nz, ny, nx)
            ! Grid Dimensions
            integer nz, ny, nx, k, j, i

            ! Potential vorticity used to locate tropopause
            real    pv(nz, ny, nx)

            ! Height at each gridpoint
            real z(nz, ny, nx)

            ! Tropopause pv value
            real pvtrop

            ! Height above and below the tropopause to produce the mask
            real height

            ! A ones and zeros mask defining points close to the dynamical tropopause
            real mask(nz, ny, nx)



            ! Local variables
            real :: trop_height

            !-START-

            ! Linearly interpolate to find the height of the tropopause
            trop_height = z(k, j, i) + ((z(k+1,j,i) - z(k,j,i)) / (pv(k+1,j,i) - pv(k,j,i))) * (pvtrop - pv(k,j,i))

            ! Fill in the mask with ones where within a threshold distance
            ! of the tropopause
            do k=1,nz
                if (z(k,j,i) < (trop_height + height).and.z(k,j,i) > (trop_height - height)) then
                    mask(k, j, i) = 1.0
                end if
            end do
        end subroutine mask_trop



        ! Create a 1's and 0's mask for near tropopause points
        !
        ! Search downwards at each gridpoint to find the dynamical tropopause defined by pvtrop. Then produces a mask
        ! of 1's if the gridpoint is within the threshold distance of the tropopause
        subroutine tropopause(pv, q, z, pvtrop, qmax, height, nz, ny, nx, mask)
            !F2PY      INTEGER, INTENT(HIDE)    :: nz
            !F2PY      INTEGER, INTENT(HIDE)    :: ny
            !F2PY      INTEGER, INTENT(HIDE)    :: nx
            !F2PY      REAL, INTENT(IN)         :: pv(nz, ny, nx)
            !F2PY      REAL, INTENT(IN)         :: q(nz, ny, nx)
            !F2PY      REAL, INTENT(IN)         :: z(nz, ny, nx)
            !F2PY      REAL, INTENT(IN)         :: pvtrop
            !F2PY      REAL, INTENT(IN)         :: qmax
            !F2PY      REAL, INTENT(IN)         :: height
            !F2PY      REAL, INTENT(OUT)        :: mask(nz, ny, nx)

            ! Grid Dimensions
            integer nz, ny, nx

            ! Potential vorticity, used to locate tropopause
            real pv(nz, ny, nx)

            ! Specific humidity, used to mask diabatic pv anomalies
            real q(nz, ny, nx)

            ! Height at each gridpoint
            real z(nz, ny, nx)

            real pvtrop
            ! Tropopause pv value

            real qmax
            ! Maximum stratospheric specific humidity

            ! Height above and below the tropopause to produce the mask
            real height

            ! A ones and zeros mask defining points close to the dynamical tropopause
            real mask(nz, ny, nx)

            !-Local Variables-
            integer :: i, j, k

            !-START-
            mask = 0.0
            do j = 1, ny
                do i = 1, nx
                    ! Search downwards for first crossing of the dynamical
                    ! tropopause
                    k = nz - 1
                    do while (pv(k,j,i) > pvtrop)
                        k = k-1
                    end do
                    ! Mask points near the tropopause
                    call mask_trop(pv, z, pvtrop, height, nz, ny, nx, k, j, i, mask)

                    ! Search for folded tropopause points
                    do while (pv(k,j,i) < pvtrop.and.k > 0)
                        k = k-1
                    end do

                    ! Ignore diabatic PV anomalies
                    if (q(k,j,i) < qmax) k = 0

                    ! Mask within a threshold distance of the tropopause
                    if (k > 0) call mask_trop(pv, z, pvtrop, height, nz, ny, nx, k, j, i, mask)

                    ! Search for the bottom of the fold
                    do while (pv(k,j,i) > pvtrop.and.k > 0)
                        k = k-1
                    end do

                    ! Ignore diabatic PV anomalies
                    if (q(k,j,i) < qmax) k = 0

                    ! Mask within a threshold distance of the tropopause
                    if (k > 0) call mask_trop(pv, z, pvtrop, height, nz, ny, nx, k, j, i, mask)
                end do
            end do
        end subroutine tropopause
end module diagnostic
