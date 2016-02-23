#INLINE F90_GLOBAL
    REAL(dp) :: DEPOSITION
    REAL(dp) :: PRESS, LAT, LON, O3COL, JO1D, JNO2
    REAL(dp) :: JDAY, JREPEAT, ALBEDO, SAREA, RP1
    INTEGER :: NOX(NVAR)
    REAL(dp):: CONSTRAIN(NVAR)
    ! indexes to these species, non-zero when not in Init_cons.dat
    INTEGER :: NEED_CH4_IDX, NEED_H2_IDX
    INTEGER :: CONST_NOX_SPC_IDX
    INTEGER :: IntTime
    CHARACTER(LEN=15), ALLOCATABLE :: file_spec_names(:)
    ! at index of of species in file, stores index into spc_names
    INTEGER, ALLOCATABLE :: file_spec_idx(:), const_method(:)
    REAL(dp), ALLOCATABLE :: fconcs(:)
    INTEGER :: n_file_specs ! number of species given in Init_cons.dat
    INTEGER :: n_points
    LOGICAL :: CONSTRAIN_RUN, CONSTRAIN_NOX
    LOGICAL :: OUTPUT_LAST_24

    INTEGER, PARAMETER :: OUT_UNIT = 6, ERROR_UNIT = 0, &
        SPEC_UNIT = 8, RATE_UNIT = 7, CONS_UNIT = 21

    !Photolysis variables
    INCLUDE '../tuv/params'
    REAL(dp) :: bs(19,kj), cs(19,kj),ds(19,kj)
    REAL(dp) :: THETA
    REAL::svj_tj(kt,kj), szas(kt), jfactno2, jfacto1d
#ENDINLINE {above lines go into MODULE dsmacc_Global}
