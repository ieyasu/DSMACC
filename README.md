This was a fun fork to make, but no one I work with uses it anymore so I'm archiving it.

The Dynamically Simple Model for Atmospheric Chemical Complexity (DSMACC) is a tropospheric chemistry box model designed to help understand the composition of the troposphere in a flexible and friendly manner. It is written to address a range of problems ranging from calculating the expected concentrations of atmospheric radicals to comparing chemistry schemes.


Read the original [DSMACC Manual](https://sites.google.com/site/dsmaccmanual/home)

[Supplemental Notes](http://wiki.seas.harvard.edu/geos-chem/index.php/DSMACC_chemical_box_model#Generate_kpp_files_based_on_your_mechanism) can be found on the GEOS-Chem Wiki.


Prerequisites
=============

DSMACC is setup to run in a POSIX (UNIX) programming environment:

1.  POSIX environment: Linux, BSD, Solaris, [Cygwin](http://www.cygwin.com/) on Windows, etc.
2.  C and Fortran compilers (e.g., [GNU](http://gcc.gnu.org/wiki/GFortran), Intel, or Portland Group)
3.  [Bison](http://www.gnu.org/software/bison/)
4.  [FLEX](http://www.gnu.org/software/flex/)


Getting DSMACC
==============

In a terminal window run the command `git clone https://github.com/ieyasu/DSMACC.git`.  This will download a copy of the model source code from GitHub into a new directory `DSMACC`.


Compile and Run
===============

In a terminal window, cd into the `DSMACC` directory, and run

1. `./configure` - pick a compiler and other compile options with (`./configure --help` to see configuration options)
2. `make` - build the model
3. `bin/dsmacc` - run the model with the supplied initial conditions


Running Your Own Chemistry Scheme
=================================

Once you've got DSMACC compiled and running (see above), you can rebuild it to use your own chemistry scheme.

1. Generate a chemistry scheme using the [Master Chemical Mechanism](http://mcm.leeds.ac.uk/MCM/)
2. Rebuild the model with `make`.  This will
   - regenerate the deposition scheme (`depos.kpp` file) for the new `organic.kpp`,
   - generate new Fortran code with KPP, and
   - compile the code, producing an updated executable `bin/dsmacc`.
3. Edit the [initial conditions](#InitCond).
4. Run the simulations with `bin/dsmacc`


Generate a Chemistry Scheme
---------------------------

1. In your web browser open http://mcm.leeds.ac.uk/MCM/
2. Click the "Browse the mechanism" link.
3. Select your compound! (Click on checkbox(es) next to compound(s) of interest).
4. Choose "Extract" from top links.
5. Choose "KPP" as the format.  **DO NOT** include either inorganic reaction or generic rate constants.
6. Click Extract to download a file called `mcm_subset.kpp`.
7. Go to your `DSMACC` directory and delete the existing organic.kpp file (`rm organic.kpp`).
8. Rename the `mcm_subset.kpp` file you downloaded to `organic.kpp` and move it to the `DSMACC` directory.


<a name="InitCond">Initial Conditions
=====================================

The model initial condition and control information are contained in a file called Init\_cons.dat. It looks like a spreadsheet with columns representing different aspects of the initial conditions (time, pressure, latitude, concentrations) and the rows representing different independent simulations of the model.

    864000
               TEMP!                 LAT!                  LON!               JDAY!                 H2O!
             PRESS!         ALBEDO!              EMISS!                 NO2!                XNO!
               XNOX!                   O3!                   CO!                 CH4!               JO1D!
                JNO2!              C5H8!               XMVK!           XMACR!
                       0!                     0!                       0!                      0!                       1!
                       0!                     0!                       1!                      0!                       0!
                       0!                     0!                       1!                      1!                       1!
                       1!                     1!                       0!                      0!
     2.6824E+02!    4.0135E+01!    -1.0947E+02!     4.6792E+01!       3.6268E-03!
     8.3000E+02!     8.5000E-01!      1.0000E-09!      5.1910E-09!       1.0000E-09!
     1.0000E+00!     3.0000E-08!      1.8284E-07!      3.4669E-06!       1.0560E-05!
      9.1500E-03!     1.0000E-09!      5.0000E-10!      5.0000E-10!


### Line 1

If the first line of the file contains an integer that tells the model to run forward that number of seconds. The output of each independent simulation is written to the files Spec\_\*.dat and Rate\_\*.dat where the \* represents an integer value representing the simulation number.

If the first line contains `-1` the model is run forwards until a diurnal steady state has been reached with output for the final timestep of all the independent simulations being output into the files Spec\_1.dat and Spec\_2.dat

If the first line contains `-2` the model is run forwards until a diurnal steady state has been reached with output for the final 24 hours for each independent simulations being written to the files Spec\_\*.dat
and Rate\_\*.dat.


### Line 2

The second line in the file lists parameters your are constraining in the model. Each parameter name is 15 characters long, separated by an exclamation mark ("!") (it is read in with the Fortran `FORMAT (100000(a15,x))`.

The following parameters can be set (case sensitive):

  Key               | Value
  ----------------- | -------------------------------------
  PRESS             | Pressure hPa
  H2O               | H2O (v/v)
  LAT               | Latitude (decimal degrees)
  LON               | Longitude (decimal degrees)
  TEMP              | Temperature (K)
  JDAY              | Julian day fractional
  O3COL             | Ozone column (Dobsons)
  ALBEDO            | The surface albedo (fraction)
  SAREA             | Surface area of aerosols (m\^2/m\^3)
  RP1               | Radius of particles (m)
  *SPECIES NAME*    | Mixing ratio of species (v/v)

If a parameter is set which is not in the above list or is a species name as defined by the chemistry of the model compilation will stop (unless it starts with an X, XOH will not cause the compilation to crash).


### Line 3

The third line sets the constraint type for each parameter. If the constraint type is 1, the parameter will be fixed for the entire run; if 0, only the initial concentration is constrained. Each constraint field is 15 characters long, separated by an exclamation mark.

Where total NOx is to be constrained it is necessary to constrain either NO or NO<sub>2</sub>, but not both. While the parameter 'NOx' must be included in the Init\_cons.dat file for total NOx to be constrained, its values in the file can be set to zero.

In order to constrain NOx the model will calculate a number every 24 hours by which the NO (or NO<sub>2</sub> if NO<sub>2</sub> is constrained in preference to NO) must be multiplied so that its modelled value remains in agreement with its observed value input into the model. All NOx species will subsequently be multiplied by this value, and hence constrained by proxy to NO (or NO<sub>2</sub>).

If neither *J* (O(^1^D)) nor *J* (NO<sub>2</sub>) are included in the input file, clear-sky values will be calculated at the altitude in question (determined from the pressure input) using TUV cross-sections at solar zenith angles varying between 0 and 90 degrees in 5 degree steps. The solar zenith angle (SZA) at which the observations were made is then calculated from the observed latitude, longitude and time of day, and a spline fit to the calculated *J*-values as a function of SZA used to determine the appropriate *J*-value.

If *J* (O(^1^D)) or *J* (NO<sub>2</sub>) are present in the input file the model will compare calculated *J*-values to their observed values and scale all calculated values accordingly.

Unless otherwise stated in the input file the model assigns [CH~4~] = 1770 ppm, [H<sub>2</sub>] = 550 ppm, and an ozone column of 260 Dobsons.


### Subsequent Lines

The remaining lines in Init_cons.dat give the initial concentrations of each parameter. These are 15-character numbers separated by exclamation marks, read in with the Fortran `FORMAT (10000(e15.4,x))`.


Output Files
============

The output of each independent simulation is written to the files Spec\_\*.dat and Rate\_\*.dat where the \* represents an integer value representing the simulation number.

The format of Spec file data is `FORMAT (100000(E25.16E3,"!"))` and the format of the Rate file data is `FORMAT (100000(E50.16E3,"!"))`.

Other Topics
============

- changing timestep
- analysing output
- maintaining directories for separate chemical schemes and initial conditions


Credits
=======

The original code development was by [Mat Evans](http://www.env.leeds.ac.uk/people/m.evans) and [Kathryn
Emmerson](http://www.env.leeds.ac.uk/people/k.emmerson) at the University of Leeds. Subsquently developement / testing work has been undertaken by [Barron Henderson](mailto:barronh@ufl.edu) (UF), [Dylan
Millet](http://www.atmoschem.umn.edu/) , [Michael Berkley](http://www.geos.ed.ac.uk/qeo/postgraduate/PhD/Applications/people/person.html?indv=1476)
(U. Edinburgh), and [Daniel Stone](http://www.chem.leeds.ac.uk/Atmospheric/Field/fage/daniel.html)
(U. Leeds). The interface for GEOS-Chem was developed by [Jingqiu Mao](http://www.people.fas.harvard.edu/%7Emao/) (Harvard) and Mat Evans.

The DSMACC model code and testcase is available from [here](http://www.github.com/barronh/DSMACC). It is based on the [KPP chemistry integration code](http://people.cs.vt.edu/%7Easandu/Software/Kpp/) written at Virginia Tech by Adrian Sandu's group.

Matthew Bishop [forked](//github.com/ieyasu/DSMACC) Barron Henderson's repository and did major cleanup work on the code to reduce repeated computation and organize the code some.

If you use the code and wish to cite the model please use:

Emmerson, KM; Evans, MJ (2009) Comparison of tropospheric gas-phase
chemistry schemes for use within global models, *ATMOS CHEM PHYS*,
**9(5)**, pp1831-1845 [doi:
10.5194/acp-9-1831-2009](http://dx.doi.org/10.5194/acp-9-1831-2009) .

