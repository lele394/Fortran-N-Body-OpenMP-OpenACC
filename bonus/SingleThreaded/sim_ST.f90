program sim

    ! I hate fortran but I can't hate it too much because it works and it's fast...
    implicit none

    integer, parameter :: n_stars = 2000
    real(8), dimension(:,:), pointer :: positions, velocities, new_positions, new_velocities, accelerations  ! Arrays declared as pointers
    real(8) :: mass = 1.0/n_stars
    real(8) :: dt = 0.005  ! Time step for the integration
    
    real(8), dimension(3) :: dp, new_acceleration


    real(8) :: epsilon = 0.005
    real(8) :: G = 1.0

    integer, parameter :: number_of_steps = 400
    integer :: i, j, n

    ! Allocate memory for arrays
    allocate(positions(n_stars, 3))
    allocate(velocities(n_stars, 3))
    allocate(new_positions(n_stars, 3))
    allocate(new_velocities(n_stars, 3))
    allocate(accelerations(n_stars, 3))

    ! ============================ MY STUFF =================================================================================
    ! Load positions
    open(unit=10, file='positions.dat', form='unformatted', access='stream')
    read(10) positions
    close(10)

    ! Load velocities
    open(unit=10, file='velocities.dat', form='unformatted', access='stream')
    read(10) velocities
    close(10)
    ! =============================================================================================================


    ! ========================================= BARNABE'S STUFF ====================================================================
    ! do i=1, n_stars
    !     ! Initialize positions
    !     positions(i, :) = 2

    !     do while (sqrt(sum(positions(i, :)**2)) > 1.0)
    !         call random_number(positions(i, :))
    !         positions(i, :) = 1-2*positions(i, :)
    !     end do

    !     ! Initialize velocities based on solid rotation
    !     velocities(i, 1) = -positions(i, 2)
    !     velocities(i, 2) = positions(i, 1)
    !     velocities(i, 3) = 0

    ! end do
    ! =============================================================================================================



    ! Opening file for save    
    ! open(unit=11, file='data/sim.dat', form='unformatted', access='stream')
    open(unit = 11, file = 'position.dat', action = 'write')

    ! First write to save starting positions
    write(11, *) n_stars, number_of_steps, 8
    write(11, *) positions

    ! === Perform Simulation iterations here.
    do n = 1, number_of_steps
        ! Loop to evolve the system by one time step

        if (mod(n, 10) == 0) then
            print *, 'Iteration ', n, ': ', n, '/', number_of_steps, ' ', real(n)/real(number_of_steps)*100.0, '%'
        end if


        ! Now apply the leapfrog integrator to all stars
        do i = 1, n_stars
            ! Compute the new position using the leapfrog integration
            new_positions(i, :) = positions(i, :) + velocities(i, :) * dt + 0.5 * accelerations(i, :) * dt**2
        
            ! Compute the acceleration for the current position
            new_acceleration = 0.0
            do j = 1, n_stars
                if (i /= j) then
                    dp = positions(j, :) - positions(i, :)  ! Compute the distance vector between stars i and j
                    new_acceleration = new_acceleration + (G * mass / sqrt(sum(dp**2) + epsilon**2)**3) * dp  ! Gravitational acceleration
                end if
            end do
        
            ! Update the velocity using the newly computed acceleration
            new_velocities(i, :) = velocities(i, :) + 0.5 * (accelerations(i, :) + new_acceleration) * dt

            ! Store the new acceleration to use in the next step
            accelerations(i, :) = new_acceleration
        end do

        ! Dump the positions to a file in the "data" folder
        write(11, *) positions
        
        ! Swap pointers (no data copying)
        call swap_pointers(positions, new_positions)
        call swap_pointers(velocities, new_velocities)
        
    end do
    
    close(11)

    print *, "integration done"

contains

    subroutine swap_pointers(a, b)
        implicit none
        real(8), dimension(:,:), pointer :: a, b
        real(8), dimension(:,:), pointer :: temp

        ! Swap the pointers
        temp => a
        a => b
        b => temp
    end subroutine swap_pointers

end program sim




