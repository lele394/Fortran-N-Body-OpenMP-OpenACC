program sim

    ! I hate fortran but I can't hate it too much because it works and it's fast...
    implicit none

    integer, parameter :: n_stars = 4000
    real(8), dimension(:,:), pointer :: positions, velocities, new_positions, new_velocities, accelerations, accelerations_i  ! Arrays declared as pointers
    real(8) :: mass = 1.0/n_stars
    real(8) :: dt = 0.005  ! Time step for the integration
    
    real(8), dimension(3) :: dp


    real(8), dimension(3) :: f_vec    ! f_vec is a 3D vector used to compute the interaction force.

    real(8) :: epsilon = 0.005
    real(8) :: G = 1.0


    integer, parameter :: number_of_steps = 200
    integer :: i, j, n

    ! Allocate memory for arrays
    ! =========================== MEMORY ALLOCATION ================================================================
    print *, "Allocating arrays"
    allocate(positions(n_stars, 3))
    allocate(velocities(n_stars, 3))
    allocate(new_positions(n_stars, 3))
    allocate(new_velocities(n_stars, 3))
    allocate(accelerations(n_stars, 3))
    allocate(accelerations_i(n_stars, 3))





    ! ============================ LOAD INITIAL POSITIONS =================================================================================
    
    print *, "Loading initial state"
    
    ! Load positions
    open(unit=10, file='positions.dat', form='unformatted', access='stream')
    read(10) positions
    close(10)

    ! Load velocities
    open(unit=10, file='velocities.dat', form='unformatted', access='stream')
    read(10) velocities
    close(10)






    ! ============================ ARRAYS DEBUG PRINT =================================================================================
    print *, "========== Arrays debug =========="

    ! Positions array
    print *, 'positions array (in bytes): ', size(positions) * kind(positions)
    print '(A,Z16)', 'Memory offset of positions array: 0x', LOC(positions)
    
    ! Velocities array
    print *, 'velocities array (in bytes): ', size(velocities) * kind(velocities)
    print '(A,Z16)', 'Memory offset of velocities array: 0x', LOC(velocities)
    
    ! New positions array
    print *, 'new_positions array (in bytes): ', size(new_positions) * kind(new_positions)
    print '(A,Z16)', 'Memory offset of new_positions array: 0x', LOC(new_positions)
    
    ! New velocities array
    print *, 'new_velocities array (in bytes): ', size(new_velocities) * kind(new_velocities)
    print '(A,Z16)', 'Memory offset of new_velocities array: 0x', LOC(new_velocities)
    
    ! Accelerations array
    print *, 'accelerations array (in bytes): ', size(accelerations) * kind(accelerations)
    print '(A,Z16)', 'Memory offset of accelerations array: 0x', LOC(accelerations)
    
    ! Accelerations_i array
    print *, 'accelerations_i array (in bytes): ', size(accelerations_i) * kind(accelerations_i)
    print '(A,Z16)', 'Memory offset of accelerations_i array: 0x', LOC(accelerations_i)
    print *, "=================================="
    

    ! =============================================================================================================


    ! Opening file for save    
    ! open(unit=11, file='data/sim.dat', form='unformatted', access='stream')
    open(unit = 11, file = 'position.dat', action = 'write')

    ! First write to save starting positions
    write(11, *) n_stars, number_of_steps, 8
    write(11, *) positions

    print *, "=============== STARTING SIM ==================="

    accelerations = 0.0
    do i = 1, n_stars - 1
        ! Initialize the new acceleration for the current star
        ! Loop over all other stars to compute the gravitational acceleration
        do j = i + 1, n_stars

            if (i /= j) then
                ! Compute the distance vector between stars i and j
                dp = positions(j, :) - positions(i, :)

                f_vec = (G * mass / (sqrt(sum(dp**2) + epsilon**2)**3)) * dp ! compute the interaction force between the objects

                accelerations(i, :) = accelerations(i, :) + f_vec
                accelerations(j, :) = accelerations(j, :) - f_vec ! Takes advantage of the anti-symmetry
            end if

        end do
    end do

    ! === Perform Simulation iterations here.
    do n = 1, number_of_steps
        ! Loop to evolve the system by one time step

        if (mod(n, 10) == 0) then
            print *, 'Iteration ', n, ': ', n, '/', number_of_steps, ' ', real(n)/real(number_of_steps)*100.0, '%'
        end if

        do i = 1, n_stars
            ! Compute the new position using the leapfrog integration
            new_positions(i, :) = positions(i, :) + velocities(i, :) * dt + 0.5 * accelerations(i, :) * dt**2
        end do

        ! resetting accelerations_i array to compute new accelerations
        accelerations_i = 0.0
        do i = 1, n_stars - 1
            ! Initialize the new acceleration for the current star
            ! Loop over all other stars to compute the gravitational acceleration
            do j = i + 1, n_stars

                if (i /= j) then
                    ! Compute the distance vector between stars i and j
                    dp = positions(j, :) - positions(i, :)

                    f_vec = (G * mass / (sqrt(sum(dp**2) + epsilon**2)**3)) * dp ! compute the interaction force between the objects

                    accelerations_i(i, :) = accelerations_i(i, :) + f_vec
                    accelerations_i(j, :) = accelerations_i(j, :) - f_vec ! Takes advantage of the anti-symmetry
                end if

            end do
        end do



        ! Now apply the leapfrog integrator to all stars
        do i = 1, n_stars
            ! Update the velocity using the newly computed acceleration
            new_velocities(i, :) = velocities(i, :) + 0.5 * (accelerations(i, :) + accelerations_i(i, :)) * dt

        end do



        ! Dump the positions to a file in the "data" folder
        write(11, *) positions
        
        ! Swap pointers (no data copying)
        call swap_pointers(positions, new_positions)
        call swap_pointers(velocities, new_velocities)
        call swap_pointers(accelerations, accelerations_i)
        
    end do
    
    close(11)
    print *, "=============== SIM END ==================="

    print *, "Deallocating arrays"
    deallocate(positions)
    deallocate(velocities)
    deallocate(new_positions)
    deallocate(new_velocities)
    deallocate(accelerations)
    deallocate(accelerations_i)


    print *, "Quitting"

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




