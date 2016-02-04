PROGRAM driver
    USE dsmacc_global
    USE dsmacc_Parameters  !ONLY: IND_*
    USE dsmacc_Rates,       ONLY: Update_SUN, Update_RCONST, J
    USE dsmacc_integrator,  ONLY: integrate
    USE dsmacc_monitor,     ONLY: SPC_NAMES, MONITOR
    USE dsmacc_Util

    IMPLICIT NONE

    REAL(dp) :: TNOX, TNOX_OLD
    REAL(dp) :: CALCJO1D, CALCJNO2
    REAL(dp) :: RSTATE(20)

    REAL(dp) :: BASE_JDAY
    INTEGER  :: ERROR, IJ
    ! Photolysis calculation variables
    REAL(dp) :: Alta

    INTEGER  :: i, CONSTNOXSPEC, JK, counter
    INTEGER :: day_tsteps, STEPS_PER_DAY
    LOGICAL :: reached_steady_state

    REAL(dp) :: NOXRATIO
    REAL(dp), POINTER :: DIURNAL_OLD(:,:), DIURNAL_NEW(:,:)
    REAL(dp), ALLOCATABLE :: DIURNAL_RATES(:,:)

    ! DT - timestep var defined by KPP in the dsmacc_Global module.
    ! It is the timestep for output, rate constant and photolysis rates.
    DT = 1200. ! 20 minutes
    STEPS_PER_DAY = (24 * 60 * 60) / INT(DT)

    STEPMIN = 0.0_dp
    STEPMAX = 0.0_dp
    RTOL(1:NVAR) = 1.0e-5_dp
    ATOL(1:NVAR) = 1.0_dp

    ! Get run parameters from Init_cons.dat file
    CALL OpenInitCons()

    IF (CONSTRAIN_RUN) THEN
        ! allocate arrays for calculating daily difference for stop condition
        ALLOCATE(diurnal_old(NVAR, STEPS_PER_DAY))
        ALLOCATE(diurnal_new(NVAR, STEPS_PER_DAY))
        ALLOCATE(diurnal_rates(NREACT, STEPS_PER_DAY))
    END IF

!$OMP PARALLEL
!$OMP DO
! XXX this loop isn't going to be parallelizable without using separate
! XXX unit numbers for each Spec_*.dat and Rate_*.dat files

    !This is the loop of different points in the Init_cons.dat file
    DO counter = 1, LINECOUNT-3
!$OMP CRITICAL
        CALL NextInitCons(counter)
!$OMP END CRITICAL

        M  = CFACTOR
        O2 = 0.21 * CFACTOR
        N2 = 0.78 * CFACTOR

        WRITE(OUT_UNIT,*) 'Starting Jday:', jday

        ! tstart is the starting time, variations due to day of year are dealt with somewhere else 
        tstart = MOD(jday,1.)*24.*60.*60.

        ! convert tstart to local time
        tstart = tstart + LON/360.*24*60*60.
        BASE_JDAY = JDAY - tstart/24./60./60.
        tend = tstart + IntTime ! IntTime is from Init_cons.dat

        WRITE(OUT_UNIT,*) 'Starting time:', tstart
        WRITE(OUT_UNIT,*) 'Ending time:',  tend
        WRITE(OUT_UNIT,*) 'Time step:', DT

        ! Set up the photolysis rates
        ! First calculate pressure altitude from altitude
        WRITE(OUT_UNIT,*) ''
        WRITE(OUT_UNIT,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
        WRITE(OUT_UNIT,*) 'Using TUV to calculate photolysis rates as a function of SZA'
        alta=(1-(press/1013.25)**0.190263)*288.15/0.00198122*0.304800/1000.
        WRITE(OUT_UNIT,*) 'Aerosol surface area', SAREA
        WRITE(OUT_UNIT,*) 'Aerosol particle radius 1', RP1
        WRITE(OUT_UNIT,*) 'Altitude =', alta
        WRITE(OUT_UNIT,*) 'Pressure =', Press
        WRITE(OUT_UNIT,*) 'Temperature =', Temp
        WRITE(OUT_UNIT,*) 'Latitude =', Lat
        WRITE(OUT_UNIT,*) 'Lon =', Lon
        WRITE(OUT_UNIT,*) 'Local Time =', Tstart/(60.*60.)
        !WRITE(OUT_UNIT,*) 'SZA =',ZENANG(int(jday),Tstart/(60.*60.),lat)*180./(4*ATAN(1.))

        ! XXX this stuff belongs over in util.inc
        IF (o3col == 0) THEN 
            o3col = 260.
            WRITE(OUT_UNIT,*) 'Ozone column not specified; using 260 Dobsons'
        ELSE
            WRITE(OUT_UNIT,*) 'Ozone column =', o3col
        END IF
        IF (albedo == 0) THEN 
            albedo = 0.1
            WRITE(OUT_UNIT,*) 'Albedo not specified; using 0.1'
        ELSE
            WRITE(OUT_UNIT,*) 'Albedo =', albedo
        END IF

        ! Calculate the photolysis rates for the run
        !$OMP CRITICAL 
        !IF (JREPEAT == 0 .OR. counter == 1) THEN 
        !    CALL set_up_photol(O3col,Albedo, alta, temp, bs,cs,ds,szas,svj_tj)
        !ELSE
        !    WRITE(OUT_UNIT,*) 'Using previously calculated photolysis params'
        !END IF
        !$OMP END CRITICAL
        !WRITE(OUT_UNIT,*) 'Photolysis rates calculated'
        !WRITE(OUT_UNIT,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
        time = tstart

        IF (CONSTRAIN_NOX) THEN ! calculate the total NOx in the model 
            TNOX_OLD=0.
            DO JK=1,NVAR
                TNOX_OLD = TNOX_OLD + C(JK) * NOX(JK)
            END DO
        END IF

        ! Initialize model state
        IF (CONSTRAIN_RUN) THEN
            DO I = 1, NVAR
                DO IJ = 1, STEPS_PER_DAY
                    DIURNAL_OLD(I,IJ) = 0.
                END DO
            END DO
            reached_steady_state = .FALSE.
            day_tsteps = 0
        END IF

        ! Calculate clear sky photolysis rates
        JFACTNO2=1.
        JFACTO1D=1.

        ! Update the rate constants
        CALL Update_RCONST()

        CALCJO1D = J(1)
        CALCJNO2 = J(4)
        WRITE(OUT_UNIT,*) 'JO1D Calc=', CALCJO1D
        WRITE(OUT_UNIT,*) 'JO1D Measre =', JO1D

        ! Calculate correction factors for the model photolysis rates
        IF (JO1D /= 0. .AND. CALCJO1D > 0.) THEN
            JFACTO1D = JO1D / J(1)
        END IF
        IF (JNO2 /= 0. .AND. CALCJNO2 > 0.) THEN
            JFACTNO2 = JNO2 / J(4)
        END IF
        IF (JNO2 == 0. .AND. JO1D /= 0.) THEN 
            JFACTNO2 = JFACTO1D
        END IF
        IF (JO1D == 0. .AND. JNO2 /= 0.) THEN 
            JFACTO1D = JFACTNO2
        END IF

        WRITE(OUT_UNIT,*) 'Correction JO1D and JNO2 by', JFACTO1D, JFACTNO2

        ! Open next Spec_*.dat and Rate_*.dat files
        CALL OpenDataFiles(counter)

        ! If we are running a free running model output the initial condition
        ! so T=0 of the output file gives the initial condition
        IF (.NOT. CONSTRAIN_RUN) THEN 
            CALL WriteCurrentData(TIME)
        END IF

        WRITE(ERROR_UNIT,*) 'Concentrations in ppb'
        IF (NMONITOR > 0) THEN
            WRITE(ERROR_UNIT,'(100000(a25,"!"))') 'TIME', (SPC_NAMES(MONITOR(i)),i=1,NMONITOR)
            WRITE(ERROR_UNIT,'(100000(E25.16E3,"!"))') time, (C(MONITOR(i))/CFACTOR * 1e9,i=1,NMONITOR)
        END IF

        ! This is the main loop for integrations
        time_loop: DO WHILE (time < TEND)

            ! Update the rate constants
            CALL Update_RCONST()

            ! XXX this copying will be SLOW. Any way to avoid it?
            COLD(:) = C(:)

            ! Integrate model forward 1 timestep
            CALL INTEGRATE( TIN = time, TOUT = time+DT, RSTATUS_U = RSTATE, &
                ICNTRL_U = (/ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /),&
                IERR_U=ERROR)

            IF (ERROR /= 1) THEN
                WRITE(OUT_UNIT,*) 'Integration error at point', counter, &
                    'time', time
                IF (NMONITOR > 0) THEN
                    WRITE(OUT_UNIT,'(100000(E25.16E3,"!"))') time, &
                        (C(MONITOR(i))/CFACTOR*1e9, i=1,NMONITOR)
                    WRITE(OUT_UNIT,'(100000(E25.16E3,"!"))') time, &
                        ((C(MONITOR(i))-COLD(MONITOR(i))) / &
                         COLD(MONITOR(i))/CFACTOR*1e9, i=1,NMONITOR)
                END IF

                ! XXX why did orig. code zero out C here?!
                !DO I=1,NVAR
                !    C(I)=0.
                !END DO
                GOTO 1000
            END IF

            ! Traps for NaN
            DO I=1,NVAR
                IF (ISNAN(C(I))) THEN
                    WRITE(ERROR_UNIT, *) 'NaN found in - XXX -, exiting'
                    ! XXX why did orig. code zero out C here?!
                    !DO i=1,NVAR
                    !    C(I)=0.
                    !END DO
                    GOTO 1000
                END IF
            END DO

            ! Update the time to reflect the integration has taken place and 
            time = RSTATE(1) ! XXX why is rstate(1) different than tout?
            JDAY = BASE_JDAY + time / 24. / 60. / 60.

            IF (CONSTRAIN_NOX) THEN
                ! Calculate total NOx in the box
                TNOX = 0
                DO I=1,NVAR
                    IF (NOX(I) /= 0) THEN 
                        TNOX = TNOX + C(I) * NOX(I)
                    END IF
                END DO

                ! Update all NOx variables so that the total NOx in the box
                ! is the same as it was
                NOXRATIO = TNOX_OLD / TNOX
                DO I=1,NVAR
                    IF (NOX(I) /= 0) THEN
                        C(I)=C(I)*NOXRATIO
                    END IF
                END DO
            END IF

            ! If constrain species concentrations if necessary
            DO I=1,NVAR
                IF (CONSTRAIN(I) > 0) THEN             
                    C(I)=CONSTRAIN(I)
                END IF
            END DO

            IF (NMONITOR > 0) THEN
                WRITE(ERROR_UNIT,'(100000(E25.16E3,"!"))') time,&
                    (C(MONITOR(i))/CFACTOR * 1e9,i=1,NMONITOR)
            END IF

            IF (CONSTRAIN_RUN) THEN ! see if run is done
                CALL CONSTRAINED_STEP
                IF (reached_steady_state) EXIT
            ELSE
                CALL WriteCurrentData(TIME)
            END IF
        END DO time_loop


1000    print *, "out of time loop"

        CALL CloseDataFiles(counter)
    END DO ! each independent 'point' to run at
!$OMP END DO NOWAIT
!$OMP END PARALLEL 

    CALL CloseInitCons()

    IF (CONSTRAIN_RUN) THEN
        DEALLOCATE(diurnal_old)
        DEALLOCATE(diurnal_new)
        DEALLOCATE(diurnal_rates)
    END IF

CONTAINS

    ! Save timestep data and, if at end of day, compare with previous day
    ! to see if model has reached steady state
    SUBROUTINE CONSTRAINED_STEP()
        REAL(dp) :: fracdiff, when
        INTEGER  :: nfrac

        ! completed another constrained timestep
        day_tsteps = day_tsteps + 1

        ! save diurnal profile of all species
        DO I=1,NVAR
            DIURNAL_NEW(I,day_tsteps) = C(I)
        END DO
        DO I=1,NREACT
            DIURNAL_RATES(I,day_tsteps) = RCONST(I)
        END DO

        IF (day_tsteps < steps_per_day) RETURN ! not end of day

        ! Do end-of-day checks

        ! We need to fiddle with the NOX to ensure it has the right concentrations and see if we have reached a steady state

        ! Sort out the NOx. Need to increase the NOx concentration so that the constrained species is right

        IF (CONSTRAIN_NOX) THEN
            ! Calculate the ratio between the value we the constrained NOx species and what we have
            ! Remember the constrained NOx species is given by the negative constrained value
            NOXRATIO = -CONSTRAIN(CONST_NOX_SPC_IDX) / C(CONST_NOX_SPC_IDX)

            ! Multiply all the NOx species by the ratio so
            DO I=1,NVAR
                IF (NOX(I) /= 0) THEN
                    C(I) = C(I) * NOXRATIO
                END IF
            END DO
        END IF
        ! Update the total amount of NOx in box
        TNOX_OLD = TNOX_OLD * NOXRATIO

        ! See how much the diurnal ratios have changed since yesterday
        fracdiff=0.
        nfrac=0
        DO I=1,NVAR
            DO JK=1,day_tsteps
                ! Add up for all species and time point in the day
                ! if there is a concentration calculated
                IF (DIURNAL_NEW(I,JK) > 1.e2 .AND. &
                    TRIM(SPC_NAMES(I)) /= 'DUMMY') THEN
                    fracdiff = fracdiff + &
                        ABS(DIURNAL_OLD(I,JK) - DIURNAL_NEW(I,JK)) / &
                        DIURNAL_NEW(I,JK)
                    nfrac = nfrac + 1
                END IF
            END DO
        END DO
        fracdiff = fracdiff / nfrac

        ! Output the diagnostic
        WRITE(OUT_UNIT,*) 'Fractional difference in diurnal profile:', fracdiff

        ! check if model has converged to steady state
        IF (fracdiff <= 1e-3) THEN
            reached_steady_state = .TRUE. ! yay~ we're done integrating

            IF (OUTPUT_LAST_24) THEN
                DO I = 1, day_tsteps
                    when = I * DT
                    CALL WriteData(when, DIURNAL_NEW(:,I), DIURNAL_RATES(:,I))
                END DO
            ELSE
                ! output last timestep
                CALL WriteCurrentData(JDAY)
            END IF
        ELSE
            CALL SWAP_DIURNAL
            day_tsteps = 0 ! move on to next day
        END IF
    END SUBROUTINE CONSTRAINED_STEP

    ! Swap diurnal new and old arrays to quickly reuse the former as the latter
    SUBROUTINE SWAP_DIURNAL
        REAL(dp), POINTER :: dp(:,:)

        dp => DIURNAL_OLD
        DIURNAL_OLD => DIURNAL_NEW
        DIURNAL_NEW => dp
    END SUBROUTINE SWAP_DIURNAL
END PROGRAM driver
