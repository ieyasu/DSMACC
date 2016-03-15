#INLINE F90_UTIL

! Opens Init_cons.dat file, gets number of lines in the file (put into
! the LINECOUNT global var), then reads the header information and sets
! some globals controlling how the model is run.
SUBROUTINE OpenInitCons()
    USE dsmacc_Global
    USE dsmacc_Monitor, ONLY: SPC_NAMES

    INTEGER :: I, J, iostat
    LOGICAL :: SPEC_CH4, SPEC_H2
    CHARACTER(LEN=15) :: str(10000)

    OPEN(UNIT=CONS_UNIT, FILE='Init_cons.dat', IOSTAT=iostat)
    IF (iostat /= 0) THEN
        PRINT *, "Error opening Init_cons.dat"
        STOP
    END IF

    n_points = 0

    READ(CONS_UNIT,*,IOSTAT=iostat) ! integration time
    IF (iostat /= 0) GOTO 99
    READ(CONS_UNIT,15,IOSTAT=iostat) str ! species list
    IF (iostat /= 0) GOTO 99
    READ(CONS_UNIT,*,IOSTAT=iostat) ! constraints
    IF (iostat /= 0) GOTO 99

    ! count species
    n_file_specs = 0
    DO I = 1, 10000
        IF (LEN_TRIM(str(I)) == 0) EXIT
        n_file_specs = n_file_specs + 1
    END DO

    IF (n_file_specs < 1) THEN
        WRITE(ERROR_UNIT,*) 'Input file has no species defined?!'
        STOP
    ELSE IF (n_file_specs > 10000) THEN
        WRITE(ERROR_UNIT,*) 'No more than 10,000 species allowed'
        STOP
    END IF

    ALLOCATE(fconcs(n_file_specs))

    ! count points
    DO WHILE (.NOT. IS_IOSTAT_END(iostat))
        fconcs(1) = -999.99
        READ(CONS_UNIT,'(10000(e15.4,x))',IOSTAT=iostat) fconcs
        IF (fconcs(1) == -999.99) THEN
            print *, "blank line found"
            EXIT ! blank line, assume end of points
        END IF
        n_points = n_points + 1
    END DO
    WRITE(OUT_UNIT, *) 'Input file has', n_file_specs, 'species at', &
        n_points, 'points'

99  IF (n_points < 1) THEN
        PRINT *, "Init_cons.dat must have an integration time, " // &
            "list of species, list of constraints and at least " // &
            "one line of initial concentrations"
        STOP
    END IF

    ! reopen file, read header info
    CLOSE(CONS_UNIT)
    OPEN(UNIT=CONS_UNIT, FILE='Init_cons.dat')

    ALLOCATE(file_spec_names(n_file_specs))
    ALLOCATE(const_method(n_file_specs))

    READ (CONS_UNIT, '(i10)') IntTime
    READ (CONS_UNIT,15) file_spec_names
    READ (CONS_UNIT,16) const_method
15 FORMAT (10000(a15,x))
16 FORMAT (10000(i15,x))

    ! run configuration
    CONSTRAIN_NOX = .FALSE.
    SPEC_CH4 = .FALSE.
    SPEC_H2 = .FALSE.
    DO I = 1, n_file_specs
        file_spec_names(I) = ADJUSTL(file_spec_names(I))

        IF (TRIM(file_spec_names(I)) == 'NOx') THEN
            CONSTRAIN_NOX = .TRUE.
            WRITE (OUT_UNIT,*) 'Constraining total NOx concentation'
        END IF

        IF (const_method(I) == 1) THEN
            IF (TRIM(file_spec_names(I)) == 'CH4') SPEC_CH4 = .TRUE.
            IF (TRIM(file_spec_names(I)) == 'H2')  SPEC_H2  = .TRUE.
        END IF
    END DO

    ! if Init_cons.dat does not supply CH4 or H2 but the model computes it,
    ! record the index of that species so it may be set in NextInitCons.
    NEED_CH4_IDX = 0
    NEED_H2_IDX = 0
    DO J = 1, NVAR
        IF ((.NOT. SPEC_CH4) .AND. &
            TRIM(SPC_NAMES(J)) == 'CH4') NEED_CH4_IDX = J
        IF ((.NOT. SPEC_H2) .AND. &
            TRIM(SPC_NAMES(J)) == 'H2') NEED_H2_IDX = J
    END DO

    IF (CONSTRAIN_NOX) THEN ! record which indexes are NOx species
        DO J = 1, NVAR
            SELECT CASE (TRIM(SPC_NAMES(J)))
            CASE ('NO2', 'NO', 'NO3', 'HONO', 'HNO2', 'HO2NO2', 'HNO4', 'PNA')
                NOX(J) = 1
            CASE ('N2O5')
                NOX(J) = 2
            CASE DEFAULT
                ! MSB this wasn't in original code - it appears that NOX was
                ! MSB assumed to already be initialized to 0 - unsafe!
                NOX(J) = 0
            END SELECT

            IF (NOX(J) /= 0) THEN
                WRITE (OUT_UNIT,*) SPC_NAMES(J),' IN NOX FAMILY'

                ! remember which species index we're constraining
                CONST_NOX_SPC_IDX = J
            END IF
        END DO
    END IF

    ! store mapping of file species name to spc_names()
    ALLOCATE(file_spec_idx(n_file_specs))
    DO I = 1, n_file_specs
        file_spec_idx(I) = 0
        DO J = 1, NVAR
            IF (TRIM(file_spec_names(I)) == TRIM(SPC_NAMES(J))) THEN
                file_spec_idx(I) = J
                EXIT
            END IF
        END DO
    END DO

    OUTPUT_LAST_24 = .FALSE.

    IF (IntTime > 0) THEN ! normal # of seconds of integration time run
        WRITE (OUT_UNIT,*) 'Integrate for ',IntTime,' seconds'
        CONSTRAIN_RUN = .FALSE.
    ELSE
        IF (IntTime == -1) THEN
            WRITE (OUT_UNIT,*) 'Integration to convergence; output final timestep'
        ELSE IF (IntTime == -2) THEN
            WRITE (OUT_UNIT,*) 'Integration to convergence; output last 24 hours'
            OUTPUT_LAST_24 = .TRUE.
        ELSE ! IntTime < -2 or == 0
            WRITE (ERROR_UNIT,*) 'First line of Init_cons.dat must be '// &
                '-1, -2 or an integration time in seconds'
            STOP
        END IF
        CONSTRAIN_RUN = .TRUE.
        IntTime = 50*24.*60.*60. ! max run time of 50 days
    END IF
END SUBROUTINE OpenInitCons


! Reads the next line of initial concentrations from Init_cons.dat (opened
! in OpenInitCons(), above).  The counter argument is only for logging;
! it does not control what is read in from the file.
SUBROUTINE NextInitCons(counter)
    USE dsmacc_Global
    USE dsmacc_Monitor

    INTEGER :: I, J, counter
    LOGICAL :: FOUND
    INTEGER :: COUNT_NOX_CONSTRAINTS, iostat

    ! Open the file with the info
    WRITE (OUT_UNIT,*) 'Initializing model point', counter
    CALL flush(OUT_UNIT)

    ! Set everything to zero
    SAREA=0.
    ALBEDO=0.
    RP1=0.
    !VAR(:) = 0.    ! MSB should not need this,
    !C(1:NVAR) = 0. ! MSB or this - in fact, was wrong max index

    JDAY = -999

    READ(CONS_UNIT,'(10000(e15.4,x))',IOSTAT=iostat) fconcs

    DO I = 1, n_file_specs
        FOUND = .TRUE.
        SELECT CASE (TRIM(file_spec_names(I)))
        CASE ('DEPOS')
            DEPOSITION = fconcs(I)
        CASE ('H2O')
            H2O = fconcs(I)
        CASE ('PRESS')
            PRESS = fconcs(I)
        CASE ('LAT')
            LAT = fconcs(I)
        CASE ('LON')
            LON = fconcs(I)
        CASE ('TEMP')
            TEMP = fconcs(I)
        CASE ('JDAY')
            JDAY = fconcs(I)
        CASE ('NOx', 'TIME(GMTs)')
            ! found - but do nothing else
        CASE ('O3COL')
            O3COL = fconcs(I)
        CASE ('ALBEDO')
            ALBEDO = fconcs(I)
        CASE ('SAREA')
            SAREA = fconcs(I)
        CASE ('RP1')
            RP1 = fconcs(I)
        CASE ('JNO2')
            IF (CONST_METHOD(I) >= 1)  JNO2 = fconcs(I)
            IF (CONST_METHOD(I) == 2) THEN 
                JREPEAT = 1
            ELSE
                JREPEAT = 0
            END IF
        CASE ('JO1D')
            IF (CONST_METHOD(I) >= 1)  JO1D = fconcs(I)
            IF (CONST_METHOD(I) == 2) THEN 
                JREPEAT = 1
            ELSE
                JREPEAT = 0
            END IF
        CASE ('')
            print *, "BUG! counted too many species in Init_cons.dat"
            STOP
        CASE DEFAULT
            FOUND = .FALSE.
        END SELECT

        J = file_spec_idx(I)
        IF (J > 0) THEN
            VAR(J) = fconcs(I)

            IF (const_method(I) == 1) THEN
                CONSTRAIN(J) = fconcs(I)
            END IF

            FOUND = .TRUE.
        END IF

        IF (.NOT. FOUND) THEN
            WRITE (OUT_UNIT,*) FILE_SPEC_NAMES(I),' NOT FOUND'
            IF (FILE_SPEC_NAMES(I)(1:1) /= 'X') STOP
            WRITE (OUT_UNIT,*) 'Starts with an X so ignored and continued'   
        END IF
    END DO ! each concentration

    IF (JDAY == -999) THEN
        WRITE(ERROR_UNIT,*) "JDAY not specified; don't know when to start"
        STOP
    END IF

    CFACTOR=PRESS*1e2*1e-6/(8.314*TEMP)*6.022E23

    H2O=H2O*CFACTOR 
    DO I=1,NVAR
        VAR(I)=VAR(I)*CFACTOR
        CONSTRAIN(I)=CONSTRAIN(I)*CFACTOR
    END DO

    ! FIND NOX species
    IF (CONSTRAIN_NOX) THEN 
        COUNT_NOX_CONSTRAINTS = 0
        DO I = 1, NVAR
            IF (NOX(I) /= 0) THEN
                CONSTRAIN(I) = -CONSTRAIN(I)

                IF (CONSTRAIN(I) /= 0 .AND. TRIM(SPC_NAMES(I)) /= 'PNA') THEN
                    COUNT_NOX_CONSTRAINTS = COUNT_NOX_CONSTRAINTS + 1
                END IF
            END IF
        END DO

        IF (COUNT_NOX_CONSTRAINTS > 1) THEN 
            WRITE(OUT_UNIT,*) 'You cannot constrain multiple NOX species:'

            DO I = 1, NVAR
                IF (NOX(I) /= 0 .AND.CONSTRAIN(I) /= 0 .AND. &
                    TRIM(SPC_NAMES(I)) /= 'PNA') THEN
                    WRITE(OUT_UNIT,*) " -", SPC_NAMES(I)
                END IF
            END DO
            WRITE(OUT_UNIT,*) 'are currently constrained'

            STOP
        END IF
    END IF

    ! set some defaults
    IF (NEED_CH4_IDX > 0) THEN
        J = NEED_CH4_IDX
        WRITE(OUT_UNIT,*) 'No CH4 specified, assuming 1770 ppbv'
        VAR(J) = 1770e-9 * CFACTOR
        CONSTRAIN(J) = VAR(J)
    END IF
    IF (NEED_H2_IDX > 0) THEN
        J = NEED_H2_IDX
        WRITE(OUT_UNIT,*) 'No H2 specified, assuming 550 ppbv'
        VAR(J) = 550e-9 * CFACTOR
        CONSTRAIN(J) = VAR(J)
    END IF
    
    IF (O3COL == 0) THEN
        WRITE(OUT_UNIT,*) 'Ozone column not specified; using 260 Dobsons'
        O3COL = 260.
    END IF
    IF (ALBEDO == 0) THEN 
        WRITE(OUT_UNIT,*) 'Albedo not specified; using 0.1'
        ALBEDO = 0.1
    END IF
END SUBROUTINE NextInitCons


SUBROUTINE CloseInitCons()
    USE dsmacc_Global

    DEALLOCATE(fconcs)
    DEALLOCATE(file_spec_names)
    DEALLOCATE(const_method)
    DEALLOCATE(file_spec_idx)

    CLOSE(CONS_UNIT)
END SUBROUTINE CloseInitCons


! Open and write headers to Spec_<counter>.dat and Rate_<counter>.dat files.
SUBROUTINE OpenDataFiles(counter)
    USE dsmacc_Parameters
    USE dsmacc_monitor
    USE dsmacc_GLOBAL

    INTEGER :: counter
    CHARACTER*64 :: filename
    CHARACTER*16 :: fnum
    INTEGER :: i

    WRITE(fnum, *) counter
    fnum = ADJUSTL(fnum)

    filename = 'Spec_'//TRIM(fnum)//'.dat'
    OPEN(SPEC_UNIT, file=filename, iostat=i)
    IF (i /= 0) THEN
        WRITE(ERROR_UNIT, *) 'Error opening', TRIM(filename)
        STOP
    END IF

    filename = 'Rate_'//TRIM(fnum)//'.dat'
    OPEN(RATE_UNIT, file=filename, iostat=i)
    IF (i /= 0) THEN
        WRITE(ERROR_UNIT, *) 'Error opening ', TRIM(filename)
        STOP
    END IF

    ! XXX look into storing JNO2FACT, JO1DFACT (others?) only once,
    ! XXX as they do not appear to change in the time step loop

    ! write headers
    WRITE(SPEC_UNIT,'(100000(a25,"!"))') 'TIME', 'LAT', 'LON', 'PRESS', &
        'TEMP', 'H2O', 'CFACTOR', 'M', 'N2', 'O2', 'JNO2FACT', 'JO1DFACT', &
        'RO2', (SPC_NAMES(LOOKAT(i)), i=1,NLOOKAT)

    WRITE(RATE_UNIT,'(100000(a50,"!"))') 'TIME', 'LAT', 'LON', 'PRESS', &
        'TEMP','H2O', 'CFACTOR', 'M', 'N2', 'O2', 'JNO2FACT', 'JO1DFACT', &
        'RO2', EQN_NAMES
END SUBROUTINE OpenDataFiles


! Write current model values to Spec and Rate data files.
SUBROUTINE WriteCurrentData(when)
    USE dsmacc_Global, ONLY: C, RCONST
    REAL(dp), INTENT(IN) :: when
    CALL WriteData(when, C, RCONST)
END SUBROUTINE WriteCurrentData

SUBROUTINE WriteData(when, concs, rates)
    USE dsmacc_Global
    USE dsmacc_Monitor
    USE dsmacc_Integrator

    REAL(dp), INTENT(IN) :: when, concs(NVAR), rates(NREACT)
    INTEGER :: i

    WRITE(SPEC_UNIT, 125) when, LAT, LON, PRESS, TEMP, H2O, CFACTOR, &
        M, N2, O2, JFACTNO2, JFACTO1D, RO2, (concs(LOOKAT(i)), i=1,NLOOKAT)
    WRITE(RATE_UNIT, 150) when, LAT, LON, PRESS, TEMP, H2O, CFACTOR, &
        M, N2, O2, JFACTNO2, JFACTO1D, RO2, rates
125 FORMAT (100000(E25.16E3,"!"))
150 FORMAT (100000(E50.16E3,"!"))

END SUBROUTINE WriteData


! Close Spec and Rate data files.
SUBROUTINE CloseDataFiles(counter)
    USE dsmacc_Global, ONLY : OUT_UNIT, SPEC_UNIT, RATE_UNIT

    INTEGER :: counter

    CLOSE(SPEC_UNIT)
    CLOSE(RATE_UNIT)

    WRITE(OUT_UNIT,*) 'Output point', counter
END SUBROUTINE CloseDataFiles

#ENDINLINE
