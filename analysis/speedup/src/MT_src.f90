program sim

    implicit none


    !==============================================================
    ! Parameters
    !==============================================================
    ! Not using env vars or shell options. Praying
    ! the compiler does it jobs and shorts to right
    ! stuff.
    !==============================================================

    ! Main parameters ========
    integer, parameter :: n_stars = 4000
    integer, parameter :: number_of_steps = 200

    ! Mass ===================
    real(8) :: mass = 1.0/n_stars
    
    ! Time parameters ========
    real(8) :: dt = 0.005  ! Time step for the integration

    ! Energy Toggles =========
    logical :: do_V = .False.
    logical :: do_E = .False.
    
    ! Physics settings =======
    real(8) :: epsilon = 0.05
    real(8) :: G = 1.0
    !==============================================================





    !==============================================================
    ! Variables
    !==============================================================
    ! Arrays ==================
    real(8), dimension(:,:), allocatable :: positions, velocities, accelerations, accelerations_i
    ! note arrays are as allocatable, fortran has built-in pointer swap with move_alloc it seems.

    ! vars ====================
    real(8) :: ec, ep
    real(8), dimension(3) :: dp
    real(8), dimension(3) :: f_vec    ! f_vec is a 3D vector used to compute the interaction force.
    integer :: i, j, n ! loop indices
    !==============================================================





    ! =========================== MEMORY ALLOCATION ================================================================
    ! Allocate memory for arrays
    print *, "Allocating arrays"
    allocate(positions(3, n_stars))
    allocate(velocities(3, n_stars))
    allocate(accelerations(3, n_stars))
    allocate(accelerations_i(3, n_stars))





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
    
    ! Accelerations array
    print *, 'accelerations array (in bytes): ', size(accelerations) * kind(accelerations)
    print '(A,Z16)', 'Memory offset of accelerations array: 0x', LOC(accelerations)
    
    ! Accelerations_i array
    print *, 'accelerations_i array (in bytes): ', size(accelerations_i) * kind(accelerations_i)
    print '(A,Z16)', 'Memory offset of accelerations_i array: 0x', LOC(accelerations_i)
    print *, "=================================="
    

    ! =============================================================================================================


    ! Opening file for save    
    open(unit = 11, file = 'position.dat', action = 'write', asynchronous="yes", form="unformatted")
    
    if (do_E) open(unit = 12, file = 'energy.dat', action = 'write')
    if (do_V) open(unit = 13, file = 'velocity.dat', action = 'write', asynchronous="yes", form="unformatted")

    ! First write to save starting positions
    write(11, asynchronous='yes') n_stars, number_of_steps, 8
    ! write(11, asynchronous='yes') transpose(positions)
    write(11, asynchronous='yes') positions

    if (do_V) write(13, asynchronous='yes') n_stars, number_of_steps, 8
    if (do_V) write(13, asynchronous='yes') velocities
    print *, "=============== STARTING SIM ==================="

    

    ! Parallel zone declaration ========================================================
    !$omp parallel private(i, j, dp, f_vec) default(none) shared(accelerations_i, &
    !$omp& mass, ep, ec, do_E, positions, velocities,&
    !$omp& dt,  accelerations, epsilon, G, do_V)


    ! Initialization of the accelerations ==============================================
    !$OMP do schedule(runtime) reduction(+:ep)
    do i = 1, n_stars - 1
        ! Initialize the new acceleration for the current star
        ! Loop over all other stars to compute the gravitational acceleration
        do j = i + 1, n_stars

            ! Compute the distance vector between stars i and j
            dp = positions(:, j) - positions(:, i)

            f_vec = (G * mass / (sqrt(sum(dp**2) + epsilon**2)**3)) * dp ! compute the interaction force between the objects

            accelerations(:, i) = accelerations(:, i) + f_vec
            accelerations(:, j) = accelerations(:, j) - f_vec ! Takes advantage of the anti-symmetry

            ! Energy computation for potential energy (x2 because of symmetry)
            if (do_E) then
                ! Ensure only one thread updates `ep`
                !!$OMP atomic ! can get rid of atomic due to reduction
                ep = ep - 2.0/sqrt(sum(dp**2) + epsilon**2)
            end if

        end do
    end do
    !$omp end do
    !=================================================================================


    ! Simulation iterations loop =====================================================
    do n = 1, number_of_steps

        !$omp master
        if (mod(n, 10) == 0) then
            print *, 'Iteration ', n, ': ', n, '/', number_of_steps, ' ', real(n)/real(number_of_steps)*100.0, '%'
        end if
        
        ! energy reset
        if (do_E) ec = 0.0
        if (do_E) ep = 0.0

        ! resetting accelerations_i array to compute new accelerations
        accelerations_i = 0.0 !Using OMP loops has no noticeable effect
        !$omp end master


        ! First leapfrog step computation =============================================
        !$OMP do schedule(runtime)
        do i = 1, n_stars
            ! Compute the new position using the leapfrog integration
            positions(:, i) = positions(:, i) + velocities(:, i) * dt + 0.5 * accelerations(:, i) * dt**2
        end do
        !$omp end do


        ! Accelerations at the next step ==============================================
        !$OMP do schedule(runtime) reduction(+:accelerations_i,ep)
        do i = 1, n_stars - 1
            ! Initialize the new acceleration for the current star
            ! Loop over all other stars to compute the gravitational acceleration
            do j = i + 1, n_stars

                ! Compute the distance vector between stars i and j
                dp = positions(:, j) - positions(:, i)

                f_vec = (G * mass / (sqrt(sum(dp**2) + epsilon**2)**3)) * dp ! compute the interaction force between the objects

                accelerations_i(:, i) = accelerations_i(:, i) + f_vec
                accelerations_i(:, j) = accelerations_i(:, j) - f_vec ! Takes advantage of the anti-symmetry

                ! Energy computation for potential energy (x2 because of symmetry)
                if (do_E) then
                    ! Ensure only one thread updates `ep`
                    !!$OMP atomic ! can get rid of atomic due to reduction
                    ep = ep - 2.0/sqrt(sum(dp**2) + epsilon**2)
                end if

            end do
        end do
        !$omp end do
        !==============================================================================


        ! Second Leapfrog step ========================================================
        !$OMP do schedule(runtime)
        do i = 1, n_stars
            ! Update the velocity using the newly computed acceleration
            velocities(:, i) = velocities(:, i) + 0.5 * (accelerations(:, i) + accelerations_i(:, i)) * dt
        end do
        !$omp end do

        
        ! Dump the positions to a file 
        !$omp master
        if (mod(n, 1) == 0) then
            ! write(11, asynchronous="yes") transpose(positions)
            write(11, asynchronous="yes") positions
            if (do_V) write(13, asynchronous='yes') velocities
        end if
        
        ! kinetic energy
        if (do_E) ec = ec + sum(velocities**2)
        if (do_E) write(12, *) 0.5*ep*mass**2, ec*0.5*mass
        ! Swap pointers (no data copying)
        call swap_array(accelerations, accelerations_i)
        !$omp end master
    end do
    !$omp end parallel
    
    close(11)
    print *, "=============== SIM END ==================="

    ! Memory cleanup =========================================
    print *, "Deallocating arrays"
    deallocate(positions)
    deallocate(velocities)
    deallocate(accelerations)
    deallocate(accelerations_i)


    print *, "Quitting"

contains

    subroutine swap_array(old, new)
        real(kind=8), allocatable, intent(inout) :: old(:,:), new(:,:)
        real(kind=8), allocatable :: temp(:,:)

        ! Use move_alloc to swap without copying data
        call move_alloc(from=old, to=temp)
        call move_alloc(from=new, to=old)
        call move_alloc(from=temp, to=new)

    end subroutine swap_array

end program sim




