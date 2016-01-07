PROGRAM driver
    USE dsmacc_global
    USE dsmacc_Parameters  !ONLY: IND_*
    USE dsmacc_Rates,       ONLY: Update_SUN, Update_RCONST, J
    USE dsmacc_integrator,  ONLY: integrate!, IERR_NAMES
    USE dsmacc_monitor,     ONLY: spc_names, MONITOR
    USE dsmacc_Util
    !USE constants

    IMPLICIT NONE
    REAL(dp) :: ENDSTATE(NVAR), total, RATIO, TNOX, TNOX_OLD
    REAL(dp) :: STARTSTATE(NVAR), TIMESCALE
    REAL(dp) :: RH, CALCJO1D, CALCJNO2
    REAL(dp) :: RSTATE(20)
    REAL(dp) :: DIURNAL_OLD(NVAR,3000), DIURNAL_NEW(NVAR,3000)
    REAL(dp) :: DIURNAL_RATES(NREACT, 3000)
    REAL(dp) :: BASE_JDAY
    INTEGER  :: ERROR, IJ
    LOGICAL :: SCREWED
    ! Photolysis calculation variables
    REAL(dp) :: Alta

    INTEGER  :: i, Daycounter, CONSTNOXSPEC, JK,counter

    REAL(dp) :: NOXRATIO, NEW_TIME
    REAL(dp) :: Fracdiff, SpeedRatio, oldfracdiff, FRACCOUNT

    STEPMIN = 0.0_dp
    STEPMAX = 0.0_dp
    RTOL(1:NVAR) = 1.0e-5_dp
    ATOL(1:NVAR) = 1.0_dp
    counter=0
    LAST_POINT=.FALSE.

    !  IF YOU WANT TO CONSTRAIN THE NOX THEN
    CONSTRAIN_NOX=.FALSE.

    !  If we are running a constrained run we want one file with the final points calculated
    IF ((CONSTRAIN_RUN .EQV. .TRUE.) .AND. (OUTPUT_LAST .EQV. .TRUE.)) THEN 
        CALL newinitsavedata(1)
    END IF

    CALL NewInitVal(0)
    !This is the loop of different points in the Init_cons.dat file

    !$OMP PARALLEL
    !$OMP DO

    DO counter=1,LINECOUNT-3
        !100 counter=counter+1

        ! Read in the next initial conditions
        WRITE(OUTPUT_UNIT,*) 'Reading in point', counter
        !$OMP CRITICAL
        CALL NewInitVal(counter)
        !$OMP END CRITICAL 
        ! Set up the output files file

        M   = CFACTOR
        O2 = 0.21 * CFACTOR
        N2 = 0.78 * CFACTOR

        WRITE(OUTPUT_UNIT,*) 'Starting Jday:',jday

        ! tstart is the starting time, variations due to day of year are dealt with somewhere else 
        tstart = (MOD(jday,1.))*24.*60.*60.              ! time start

        ! convert tstart to local time
        tstart = tstart+LON/360.*24*60*60.
        BASE_JDAY = JDAY - tstart / 24. /60./60.
        ! tend is the end time. IntTime is determined from the Init_cons.dat file
        tend = tstart + IntTime    

        !dt is the output timestep and the timestep between times rate constants and notably photolysis rates are calcualted
        dt = 1200.

        WRITE(OUTPUT_UNIT,*) 'Starting time:',tstart
        WRITE(OUTPUT_UNIT,*) 'Ending time:', tend
        WRITE(OUTPUT_UNIT,*) 'Time step:', dt

        !    Set up the photolysis rates
        !    First calculate pressure altitude from altitude
        WRITE(OUTPUT_UNIT,*) 'hvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhv'
        WRITE(OUTPUT_UNIT,*) 'Using TUV to calculate photolysis rates as a function of SZA'
        alta=(1-(press/1013.25)**0.190263)*288.15/0.00198122*0.304800/1000.
        WRITE(OUTPUT_UNIT,*) 'Aerosol surface area', SAREA
        WRITE(OUTPUT_UNIT,*) 'Aerosol particle radius 1', RP1
        WRITE(OUTPUT_UNIT,*) 'Altitude =', alta
        WRITE(OUTPUT_UNIT,*) 'Pressure =', Press
        WRITE(OUTPUT_UNIT,*) 'Temperature =', Temp
        WRITE(OUTPUT_UNIT,*) 'Latitude =', Lat
        WRITE(OUTPUT_UNIT,*) 'Lon =', Lon
        WRITE(OUTPUT_UNIT,*) 'Local Time =', Tstart/(60.*60.)
        !        WRITE(OUTPUT_UNIT,*) 'SZA =',ZENANG(int(jday),Tstart/(60.*60.),lat)*180./(4*ATAN(1.))
        IF (o3col .eq. 0) THEN 
            o3col=260.
            WRITE(OUTPUT_UNIT,*) 'Ozone column not specified using 260 Dobsons'
        ELSE
            WRITE(OUTPUT_UNIT,*) 'Ozone column =', o3col
        END IF

        IF (albedo .eq. 0) THEN 
            albedo=0.1
            WRITE(OUTPUT_UNIT,*) 'Albedo not specified using 0.1'
        ELSE
            WRITE(OUTPUT_UNIT,*) 'Albedo =', albedo
        END IF
        !    Calculate the photolysis rates for the run
        !$OMP CRITICAL 
        !        IF (JREPEAT .EQ. 0 .OR. COUNTER .EQ. 1) THEN 
        !            CALL set_up_photol(O3col,Albedo, alta, temp, bs,cs,ds,szas,svj_tj)
        !        ELSE
        !            WRITE(OUTPUT_UNIT,*) 'Using previously calculated photolysis params'
        !        END IF
        !$OMP END CRITICAL
        !        WRITE(OUTPUT_UNIT,*) 'Photolysis rates calculated'
        !        WRITE(OUTPUT_UNIT,*) 'hvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhv'
        time = tstart


        OLDFRACDIFF=0.

        ! If NOx is being constrained calculate the total NOx in the model 
        IF (CONSTRAIN_NOX) THEN 
            TNOX_OLD=0.
            DO JK=1,NVAR
                TNOX_OLD=TNOX_OLD+C(JK)*NOX(JK)
            END DO
        END IF

        ! Define the initial state of the model 
        DO I=1,NVAR
            STARTSTATE(I)=C(I)
            DO IJ=1,3000
                DIURNAL_OLD(I,IJ)=0.
            END DO
        END DO

        ! Calculate clear sky photolysis rates
        JFACTNO2=1.
        JFACTO1D=1.

        ! Update the rate constants
        CALL Update_RCONST()

        CALCJO1D = J(1)
        CALCJNO2 = J(4)

        WRITE(OUTPUT_UNIT,*) 'JO1D Calc=', CALCJO1D
        WRITE(OUTPUT_UNIT,*) 'JO1D Measre =', JO1D
        ! Calcualte correction factors for the model photolysis rates
        IF (JO1D .NE. 0. .AND. CALCJO1D .GT. 0.) THEN
            JFACTO1D=JO1D/J(1)
        END IF

        IF (JNO2 .NE. 0. .AND. CALCJNO2 .GT. 0.) THEN
            JFACTNO2=JNO2/J(4)
        END IF

        IF (JNO2 .EQ. 0. .AND. JO1D .NE. 0.) THEN 
            JFACTNO2=JFACTO1D
        END IF

        IF (JO1D .EQ. 0. .AND. JNO2 .NE. 0.) THEN 
            JFACTO1D=JFACTNO2
        END IF

        WRITE(OUTPUT_UNIT,*) 'Correction JO1D and JNO2 by', JFACTO1D,JFACTNO2

        ! If we are running a non-constrained run then we want one file per input in the Init_cons.dat file
        IF ((CONSTRAIN_RUN .EQV. .FALSE.) .OR. (OUTPUT_LAST .EQV. .FALSE.)) THEN 
            CALL NEWINITSAVEDATA(COUNTER)
        END IF


        ! If we are running a free running model output the initial condition so T=0 of the output file gives the initial condition
        IF (CONSTRAIN_RUN .EQV. .FALSE.) THEN 
            CALL NEWSAVEDATA()
        END IF

        ! Set up a counter to count the number time that the model has been run for 
        Daycounter=0
        WRITE(ERROR_UNIT,*)'Concentrations in ppb'
        IF (NMONITOR > 0) THEN
            WRITE(ERROR_UNIT,'(100000(a25,"!"))') 'TIME', (SPC_NAMES(MONITOR(i)),i=1,NMONITOR)
            WRITE(ERROR_UNIT,'(100000(E25.16E3,"!"))') time, (C(MONITOR(i))/CFACTOR * 1e9,i=1,NMONITOR)
        END IF
        ! This is the main loop for integrations
        time_loop: DO WHILE (time < TEND)

            ! Update the rate constants
            CALL Update_RCONST()
            COLD(:) = C(:)
            ! Integrate the model forwards 1 timestep
            CALL INTEGRATE( TIN = time, TOUT = time+DT, RSTATUS_U = RSTATE, &
                ICNTRL_U = (/ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /),&
                IERR_U=ERROR)

            IF (ERROR .NE. 1) THEN 
                WRITE(OUTPUT_UNIT,*) 'Integration error.'
                WRITE(OUTPUT_UNIT,*) 'Skipping point'
                IF (NMONITOR > 0) THEN
                    WRITE(OUTPUT_UNIT,'(100000(E25.16E3,"!"))') time,&
                        (C(MONITOR(i))/CFACTOR * 1e9,i=1,NMONITOR)
                    WRITE(OUTPUT_UNIT,'(100000(E25.16E3,"!"))') time,&
                        ((C(MONITOR(i))-COLD(MONITOR(i)))/COLD(MONITOR(i))/CFACTOR&
                        * 1e9,i=1,NMONITOR)
                END IF
                DO I=1,NVAR
                    C(I)=0.
                END DO
                GOTO 1000
            END IF
            ! Traps for NaN
            SCREWED=.FALSE.
            DO I=1,NVAR
                IF (ISNAN(C(I))) SCREWED=.TRUE.
            END DO

            IF (SCREWED) THEN 
                SCREWED=.FALSE.
                DO i=1,NVAR
                    C(I)=0.
                END DO
                GOTO 1000
            END IF
            ! Update the time to reflect the integration has taken place and 
            time = RSTATE(1)
            JDAY = BASE_JDAY + time / 24. /60./60.
            Daycounter=Daycounter+1

            ! If we are constraining NOx then:
            IF (CONSTRAIN_NOX) THEN 
                ! Calcualte the total NOx in the box
                TNOX=0
                DO I=1,NVAR
                    IF (NOX(I) .NE. 0) THEN 
                        TNOX=TNOX+C(I)*NOX(I)
                    END IF
                END DO

                ! Update all NOx variables so that the total NOx in the box is the same as it was
                DO I=1,NVAR
                    IF (NOX(I) .NE. 0) THEN 
                        C(I)=C(I)*TNOX_OLD/TNOX
                    END IF
                END DO
            END IF

            ! If constrain species concentALrations if necessary
            DO I=1,NVAR
                IF (CONSTRAIN(I) .GT. 0) THEN             
                    C(I)=CONSTRAIN(I)
                END IF
            END DO

            ! If we are not doing a constrained run then output the concentrations
            IF (CONSTRAIN_RUN .EQV. .FALSE.) THEN 
                CALL NEWSAVEDATA()
            END IF

            IF (NMONITOR > 0) THEN
                WRITE(ERROR_UNIT,'(100000(E25.16E3,"!"))') time,&
                    (C(MONITOR(i))/CFACTOR * 1e9,i=1,NMONITOR)
            END IF

            ! If we are doing a constrained run we need to store the diurnal profile of all the species
            IF (CONSTRAIN_RUN .EQV. .TRUE.) THEN
                DO I=1,NVAR
                    DIURNAL_NEW(I,DAYCOUNTER)=C(I)
                END DO

                DO I=1,NREACT
                    DIURNAL_RATES(I,DAYCOUNTER)=RCONST(I)
                END DO

                ! Are we at the end of a day?
                ! If so we need to 
                !   1) fiddle with the NOX to ensure it has the right concentrations see if we have reached a steady state
                IF (DAYCOUNTER*DT .GE. 24.*60.*60.) THEN
                    ! Sort out the NOx. Need to increase the NOx concentration so that the constrained species is right
                    ! What is  the constrained NOx species? Put result into CONSTNOXSPEC
                    DO I=1,NVAR
                        IF (NOX(I) .NE. 0) THEN 
                            IF (CONSTRAIN(I) .LT. 0) THEN 
                                CONSTNOXSPEC=I
                            END IF
                        END IF
                    END DO

                    ! If we are constraining NOx then:
                    IF (CONSTRAIN_NOX) THEN 
                        ! Calculate the ratio between the value we the constrained NOx species and what we have
                        ! Remember the constrained NOx species is given by the negative constrained value
                        NOXRATIO=-CONSTRAIN(CONSTNOXSPEC)/C(CONSTNOXSPEC)

                        ! Multiply all the NOx species by the ratio so 
                        DO I=1,NVAR
                            IF (NOX(I) .NE. 0) THEN 
                                C(I)=C(I)*NOXRATIO
                            END IF
                        END DO
                    END IF
                    ! Update the total amount of NOx in box
                    TNOX_OLD=TNOX_OLD*NOXRATIO

                    ! Lets see how much the diurnal ratios have changed since the last itteration

                    ! Frac diff is our metric for how much it has changed 
                    FRACDIFF=0.
                    FRACCOUNT=0.
                    ! Add up for all species and for each time point in the day
                    DO I=1,NVAR
                        DO JK=1,DAYCOUNTER
                            !If there is a concentration calculated
                            IF (DIURNAL_NEW(I,JK) .GT. 1.e2 .AND. &
                                TRIM(SPC_NAMES(I)) .NE. 'DUMMY') THEN 
                                !Calculate the absolute value of the fractional difference and add it on
                                ! Increment the counter to calculate the average
                                FRACDIFF=FRACDIFF+&
                                    ABS(DIURNAL_OLD(I,JK)-DIURNAL_NEW(I,JK))/&
                                    DIURNAL_NEW(I,JK)
                                FRACCOUNT=FRACCOUNT+1
                            END IF
                        END DO
                    END DO

                    ! Calculate the average fractional difference
                    FRACDIFF=FRACDIFF/FRACCOUNT

                    ! Output the diagnostic
                    WRITE(OUTPUT_UNIT,*)&
                        'Fraction difference in the diurnal profile:', FRACDIFF

                    ! Store the new diurnal profile as the old one so we can compare with the next day
                    DO I=1,NVAR
                        DO JK=1,DAYCOUNTER
                            DIURNAL_OLD(I,JK)=DIURNAL_NEW(I,JK)
                        END DO
                    END DO

                    ! reset the day counter to 0

                    ! if the system has converged then end the simulation for this point
                    IF (FRACDIFF .LE. 1e-3) THEN
                        GOTO 1000
                    END IF
                    DAYCOUNTER=0
                    OLDFRACDIFF=FRACDIFF
                END IF
            END IF
        END DO time_loop


1000    IF ((CONSTRAIN_RUN .EQV. .TRUE.) .AND. (OUTPUT_LAST .EQV. .FALSE.)) THEN 
            CALL newsavedata()
        END IF

        IF (OUTPUT_LAST .EQV. .TRUE.) THEN 
            DO I=1,DAYCOUNTER
                NEW_TIME=I*DT
                WRITE (SPEC_UNIT,999) NEW_TIME,LAT, LON, PRESS, TEMP,H2O, CFACTOR, RO2, &
                    (DIURNAL_NEW(JK,I),JK=1,NVAR)
                WRITE (RATE_UNIT,999) NEW_TIME,LAT, LON, PRESS, TEMP,H2O, CFACTOR,& 
                    (DIURNAL_RATES(JK,I),JK=1,NREACT)
            END DO
999         FORMAT(E24.16,100000(1X,E24.16))
        END IF

        IF (CONSTRAIN_RUN .EQV. .FALSE.) THEN
            ! close output file
            CALL newclosedata()
        END IF
        WRITE(OUTPUT_UNIT,*) 'Outputed point', counter

    END DO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL 

    !if (.NOT. LAST_POINT) goto 100

    IF (CONSTRAIN_RUN .EQV. .TRUE.) THEN 
        CALL Newclosedata()
    END IF
END PROGRAM driver
