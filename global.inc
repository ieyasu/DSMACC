#INLINE F90_GLOBAL
    REAL(dp) :: PRESS, LAT, LON, O3COL, JO1D, JNO2
    REAL(dp) :: JDAY, JREPEAT, ALBEDO, SAREA, RP1
    INTEGER :: INIT_TIME, NOX(NVAR)
    REAL(dp):: CONSTRAIN(NVAR)
    LOGICAL :: SPEC_CH4, SPEC_H2
    INTEGER :: IntTime
    LOGICAL :: CONSTRAIN_NOX, CONSTRAIN_RUN
    LOGICAL :: LAST_POINT = .false., OUTPUT_LAST_24
    CHARACTER(LEN=15) :: spec_name(10000) ! as listed in Init_cons.dat
    INTEGER :: const_method(10000), LINECOUNT
    ! C - Concentration of all species
    REAL(kind=dp) :: COLD(NSPEC)
    INTEGER, PARAMETER :: OUT_UNIT = 6, ERROR_UNIT = 0, &
        SPEC_UNIT = 8, RATE_UNIT = 7, CONS_UNIT = 21

    !Photolysis variables
    INCLUDE '../tuv/params'
    REAL(dp) :: bs(19,kj), cs(19,kj),ds(19,kj)
    REAL(dp) :: THETA
    REAL::svj_tj(kt,kj), szas(kt), jfactno2, jfacto1d
#ENDINLINE {above lines go into MODULE dsmacc_Global}
