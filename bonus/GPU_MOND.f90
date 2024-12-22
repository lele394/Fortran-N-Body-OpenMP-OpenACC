program sim
    ! use openacc

    implicit none
    ! ! ==============================================================
    ! !   ~~CHECK MT FINAL FOR ARRAY INDEX SWITCH AND ADD IT~~
    ! !       ^  ARRAY SWITCH IS SLOWER ON GPU   ^
    ! ! ==============================================================

    !==============================================================
    ! Notes
    !==============================================================
    ! unit = 11, file = 'position.dat'
    ! unit = 12, file = 'energy.dat'
    ! unit = 13, file = 'velocity.dat'
    ! unit = 14, file = 'last_velocities.dat'



    !==============================================================
    ! Parameters
    !==============================================================
    ! Not using env vars or shell options. Praying
    ! the compiler does it jobs and shorts to right
    ! stuff.
    !==============================================================

    ! Main parameters ========
    integer, parameter :: n_stars = 100000 !SE
    integer, parameter :: number_of_steps = 5 !SE
    integer :: save_modulo = 1 !save state every X
    integer :: display_modulo = 1 ! display remaining steps every X

    ! MOND potential =========
    real(8) :: a0 = 1.2e-10

    ! Mass ===================
    real(8) :: mass = 1.0/n_stars
    
    ! Time parameters ========
    real(8) :: dt = 0.005  ! Time step for the integration

    ! Energy Toggles =========
    logical :: do_V = .False.
    logical :: do_E = .False.
    logical :: save_last_V = .True. ! Allows to load simulation from a previous state 
    !               ^  Saving works, change velocities files i guess [UNTESTED]
    
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
    real(8) :: ec, ep, ep_local
    real(8), dimension(3) :: dp
    real(8), dimension(3) :: f_vec    ! f_vec is a 3D vector used to compute the interaction force.
    integer :: i, j, n ! loop indices

    ! ETF computation ========
    real(8) :: start_time, current_time, remaining_time

    !==============================================================

    ! ETF start ==============
    call cpu_time(start_time)

    ! =========================== MEMORY ALLOCATION ================================================================
    ! Allocate memory for arrays
    print *, "Allocating arrays"
    allocate(positions(n_stars, 3))
    allocate(velocities(n_stars, 3))
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
    write(11, asynchronous='yes') positions
    
    if (do_V) write(13, asynchronous='yes') n_stars, number_of_steps, 8
    if (do_V) write(13, asynchronous='yes') velocities
    print *, "=============== STARTING SIM ==================="

    


    !$acc enter data copyin(mass, positions, velocities, dt, epsilon, G, do_V, do_E) &
    !$acc& create(accelerations_i, accelerations, ep, ec)

    print *, "Acceleration initialization..."

    ! Initialization of the accelerations ==============================================
        !$acc parallel loop  private(i, j, dp, f_vec, ep_local)
        do i = 1, n_stars
        ! Initialize the new acceleration for the current star
        ! Loop over all other stars to compute the gravitational acceleration
        f_vec = 0.0
        ep_local = 0.0
        do j = 1, n_stars

            ! Compute the distance vector between stars i and j
            dp = positions(j, :) - positions(i, :)

            f_vec = f_vec + (G * mass / ((sum(dp**2) + epsilon**2) * sqrt(sum(dp**2) + epsilon**2))) * &
            ((sum(dp**2) + epsilon**2) / (a0 * sqrt((sum(dp**2) + epsilon**2) / a0**2 + 1))) * dp ! compute the interaction force between the objects


            ! Energy computation for potential energy
            if (do_E) then
                ep_local = ep_local - 1.0/sqrt(sum(dp**2) + epsilon**2)
            end if
        end do

        !$acc atomic
        accelerations(i, 1) = accelerations(i, 1) + f_vec(1)
        !$acc atomic
        accelerations(i, 2) = accelerations(i, 2) + f_vec(2)
        !$acc atomic
        accelerations(i, 3) = accelerations(i, 3) + f_vec(3)

        !$acc atomic
        ep = ep + ep_local

    end do
    !$acc end parallel loop
    !!$acc end parallel loop
    !=================================================================================

    print *, "Acceleration initialization Done!"

    ! Simulation iterations loop =====================================================
    do n = 1, number_of_steps

        if (mod(n, display_modulo) == 0) then
            call cpu_time(current_time)
            
            ! Calculate the remaining time (ETF) in seconds
            remaining_time = current_time / real(n) * (number_of_steps - n)
            
            ! Print the formatted output in a single statement
            print '(A, I6.2,A, I6.2, ":", F6.2, A, I3.2, ":", I2.2, ":", I2.2)', &
            'Iteration ', n, '/', number_of_steps,  &
            real(n) / real(number_of_steps) * 100.0, '%  ETF :', &
            floor(remaining_time / 3600), mod(floor(remaining_time / 60), 60),  mod(floor(remaining_time), 60)
        
        end if

        ! create parallel zone for this step
        ! energy reset

        if (do_E) then 
            ec = 0
            ep = 0
            !$acc update device(ep, ec)
        end if







        ! resetting accelerations_i array to compute new accelerations
        !$acc parallel private(i)

        !$acc loop
        do i = 1, n_stars
            ! Compute the new position using the leapfrog integration
            accelerations_i(i, :)= 0.0
        end do



        ! First leapfrog step computation =============================================
        !$acc loop
        do i = 1, n_stars
            ! Compute the new position using the leapfrog integration
            positions(i, :) = positions(i, :) + velocities(i, :) * dt + 0.5 * accelerations(i, :) * dt**2
        end do
        
        !$acc end parallel 


        

        ! Accelerations at the next step ==============================================
        !$acc parallel loop  private(i, j, dp, f_vec, ep_local)
        do i = 1, n_stars
            ! Initialize the new acceleration for the current star
            ! Loop over all other stars to compute the gravitational acceleration
            f_vec = 0.0
            ep_local = 0.0
            do j = 1, n_stars
    
                ! Compute the distance vector between stars i and j
                dp = positions(j, :) - positions(i, :)
    
                f_vec = f_vec + (G * mass / ((sum(dp**2) + epsilon**2) * sqrt(sum(dp**2) + epsilon**2))) * &
                ((sum(dp**2) + epsilon**2) / (a0 * sqrt((sum(dp**2) + epsilon**2) / a0**2 + 1))) * dp ! compute the interaction force between the objects
    
    
                ! Energy computation for potential energy
                if (do_E) then
                    ep_local = ep_local - 1.0/sqrt(sum(dp**2) + epsilon**2)
                end if
            end do

            !$acc atomic
            accelerations_i(i, 1) = accelerations_i(i, 1) + f_vec(1)
            !$acc atomic
            accelerations_i(i, 2) = accelerations_i(i, 2) + f_vec(2)
            !$acc atomic
            accelerations_i(i, 3) = accelerations_i(i, 3) + f_vec(3)

            !$acc atomic
            ep = ep + ep_local

        end do
        !$acc end parallel loop
        !==============================================================================














        ! Second Leapfrog step ========================================================
        !$acc parallel private(i)
        
        !$acc loop
        do i = 1, n_stars
            ! Update the velocity using the newly computed acceleration
            velocities(i, :) = velocities(i, :) + 0.5 * (accelerations(i, :) + accelerations_i(i, :)) * dt
        end do

        !$acc wait

        !$acc loop
        do i = 1, n_stars
            accelerations(i, :) = accelerations_i(i, :)
        end do
        !$acc end parallel

        
        
        
        !$acc update host(positions)
        if (mod(n, save_modulo) == 0) then
            ! update data CPU side
            write(11, asynchronous="yes") positions
            if (do_E) then
                !$acc update host(velocities, ep, ec)
                if (do_V) write(13, asynchronous='yes') velocities
                ec = ec + sum(velocities**2)
                write(12, *) 0.5*ep*mass**2, ec*0.5*mass
            end if
        end if
        ! Swap pointers (no data copying)



    end do
    
    close(11)
    print *, "=============== SIM END ==================="


    if (save_last_V) then
        !$acc update host(velocities)
        open(unit = 14, file = 'last_velocities.dat', action = 'write', form="unformatted")
        write(14) velocities
    end if

    ! Memory cleanup =========================================
    print *, "Deallocating arrays"
    deallocate(positions)
    deallocate(velocities)
    deallocate(accelerations)
    deallocate(accelerations_i)

    !$acc exit data delete(mass, positions, velocities, dt, epsilon, G, do_V, do_E)
    !$acc exit data delete(accelerations_i, accelerations, ep, ec)

    print *, "Quitting"



end program sim




