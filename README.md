Introduction
============

Overview
--------

The Dynamically Simple Model for Atmospheric Chemical Complexity
(DSMACC) is a tropospheric chemistry box model designed to help
understand the composition of the troposphere in a flexible and
friendly manner. It is written to address a range of problems ranging
from calculating the expected concentrations of atmospheric radicals to
comparing chemistry schemes.

Credits
-------

The original code development was by [Mat
Evans](http://www.env.leeds.ac.uk/people/m.evans) and [Kathryn
Emmerson](http://www.env.leeds.ac.uk/people/k.emmerson) at the
University of Leeds. Subsquently developement / testing work has been
undertaken by [Barron Henderson](mailto:barronh@ufl.edu) (UF), [Dylan
Millet](http://www.atmoschem.umn.edu/) , [Michael
Berkley](http://www.geos.ed.ac.uk/qeo/postgraduate/PhD/Applications/people/person.html?indv=1476)
(U. Edinburgh), and [Daniel
Stone](http://www.chem.leeds.ac.uk/Atmospheric/Field/fage/daniel.html)
(U. Leeds). The interface for GEOS-Chem was developed by [Jingqiu
Mao](http://www.people.fas.harvard.edu/%7Emao/) (Harvard) and Mat Evans.

The DSMACC model code and testcase is available from
[here](http://www.github.com/barronh/DSMACC). It is based on the [KPP
chemistry integration
code](http://people.cs.vt.edu/%7Easandu/Software/Kpp/) written at
Virginia Tech by Adrian Sandu's group.

Matthew Bishop [forked](//github.com/ieyasu/DSMACC) Barron Henderson's
repository and made some major cleanup work on the code to reduce
repeated computation and organize the code some.

If you use the code and wish to cite the model please use:

Emmerson, KM; Evans, MJ (2009) Comparison of tropospheric gas-phase
chemistry schemes for use within global models, *ATMOS CHEM PHYS*,
**9(5)**, pp1831-1845 [doi:
10.5194/acp-9-1831-2009](http://dx.doi.org/10.5194/acp-9-1831-2009) .

Instructions
============

[Supplemental instructions](http://wiki.seas.harvard.edu/geos-chem/index.php/DSMACC_chemical_box_model#Generate_kpp_files_based_on_your_mechanism) can be found on the GEOS-Chem Wiki.

Overview
--------

This descibes the steps necessary to run the model. These instructions
assume basic knowledge of Linux Bash or C-Shell.

1.  [Install prerequisites](Install%20Prerequisites)
2.  [Run the test case](Test%20Case)
3.  [Modifying the model chemistry scheme
    (optional)](Chemistry%20Configuration)
4.  [Re-compiling and compiling the model
    (optional)](Compilation%20Instructions)
5.  [Creating your own the model
    inputs](Inputs%20and%20Initial%20Conditions)
6.  [Running the model](Running%20the%20model)

### <a id="Install%20Prerequisites"/>Install prerequisites

DSMACC has relatively few requirements and all can be obtained free of
charge. DSMACC itself requires a basic Linux development environment.
DSMACC also uses KPP (distributed with DSMACC), which requires a parser
and a lexer.

1.  Linux/Unix environment (including [cygwin](http://www.cygwin.com/)
    on windows)
2.  C and Fortran compilers (e.g.,
    [GNU](http://gcc.gnu.org/wiki/GFortran), Intel, or Portland Group)
3.  [Bison](http://www.gnu.org/software/bison/)
4.  [FLEX](www.gnu.org/software/flex/)

### <a id="Test%20Case"/>Test Case

The test case assumes that you have KPP installed with kpp in your PATH
and the KPP\_HOME environment variable set. (if not see [KPP
Configuration](KPP Configuration)) These instructions will build the
DSMACC model, and run it with example inputs. It will then compare your
results with a set of archived results.

1. Open a terminal and go to a directory where you want to put DSMACC.
2. Get a copy of the source code with `git clone https://github.com/ieyasu/DSMACC.git`
3. `cd DSMACC`
4. `./configure` (`./configure --help` to see configuration options)
5. `make`
6. (optional) `make check`

This may take some time, but will give you a list of pass/fails for a
diurnal constrained steady state run.

### <a id="Compilation%20Instructions"/>Compilation Instructions
#### Standard KPP input

Which chemistry scheme you are going to run in the model is described by
the series of files at the beginning of the dsmacc.kpp file. For example
if you have the line \`\#INCLUDE inorganic.kpp\` in your dsmacc.kpp file
the chemistry included in that file will be included in chemistry
scheme. For details about how to write these files see the kpp user
manual.

#### MCM input

If the MCM is being used the chemistry deck is composed of two files.
The first inorganic.kpp contains the inorganic chemistry. This is the
same for all simulations. A copy is included in the DSMACC distribution.
We make some efforts to keep this file up to date but make no guarantees
about this. The second organic.kpp is generated from the MCM web-site
and contains the organic aspects of the chemistry scheme.

#### To generate the organic.kpp file

Go to
[http://mcm.leeds.ac.uk/MCMv3.3.1/extract.htt](http://mcm.leeds.ac.uk/MCM/roots.http))
. Then select the base VOCs that you would like to simulate the
chemistry of by selecting the check boxes associated with each species.
Once all the VOC's have been selected click 'Added Selection to Marked
List.' Then click on 'Extract' from the menu at the top of the page.
Select KPP, experimental KPP format. Click on Extract Save the generated
page as organic.kpp

#### Using the GEOS-Chem globchem.dat file

**(this is under developement)**

This uses the GEOS-Chem globchem.dat file to generate a series of files
(globchem.eqn, globchem.spc and globchem.def.) in kpp format to use in
the model;

change all species in globchem.dat to Active if possible. there are bugs
in kpp code for some reactions with two lines. In order to make it work
right, we have to modify the globchem.dat for reactions, HNO4+M; N2O5+M;
HO2+MCO3; DMS+OH; change the reaction rate for photolysis reactions to
match up with TUV if you have 5 lines for each reaction in glochem.dat
(such as v8-02-05) change line 225 in Geos2kppbox\_parser.pl from "4" to
"5" type `Geos2kppbox_parser  globchem.dat` to generate globchem.eqn,
globchem.spc and globchem.def.

There are two steps to compiling the model. The first is to get kpp to
generate the appropriate FORTRAN code and the second is to compile the
FORTAN code into a program.

### <a id="Compilation%20Instructions"/> 

#### Generate the code

To construct the model type `kpp dsmacc.kpp dsmacc`. Information about
error messages etc can be found in the kpp user handbook. A series of
FORTRAN 90 (\*.f90) files will now be generated by KPP which contain the
chemistry scheme you have requested.

#### Making the model

To make the model as an executable type `make`. The model will be
constructed and named model.

### <a id="Inputs%20and%20Initial%20Conditions"/> Input and Initial Conditions

The model initial condition and control information are contained in a
file called *Init\_cons.dat*. It looks like a spreadsheet with columns
representing different aspects of the initial conditions (time,
pressure, latitude, concentrations) and the rows representing different
independent simulations of the model.

#### Line 2

If the first line of the file contains a ***positive integer*** this
tells the model to run forward that number of seconds. The
output of each independent simulation is written to the files
Spec\_\*.dat and Rate\_\*.dat where the \* represents an integer value
representing the simulation number.

If the first line contains `-1` the model is run forwards until a
diurnal steady state has been reached with output for the final timestep
of all the independent simulations being output into the files
Spec\_1.dat and Spec\_2.dat

If the first line contains `-2` the model is run forwards until a
diurnal steady state has been reached with output for the final 24 hours
for each independent simulations being written to the files Spec\_\*.dat
and Rate\_\*.dat.

#### Line 2

The second row in the file should contain the parameters to be input
into the model. Each parameter name is 15 characters long, separated by
an space (it is read in by the FORTRAN format statement
‘100000(a15,x))’. The following parameters can be set (case sensitive).

  Key               Value
  ----------------- -------------------------------------
  PRESS             Pressure hPa
  H2O               H2O (v/v)
  LAT               Latitude (decimal degrees)
  LON               Longitude (decimal degrees)
  TEMP              Temperature (K)
  JDAY              Julian day fractional
  O3COL             Ozone column (Dobsons)
  ALBEDO            The surface albedo (fraction)
  SAREA             Surface area of aerosols (m\^2/m\^3)
  RP1               Radius of particles (m)
  *SPECIES NAME*    Mixing ratio of species (v/v)

If a parameter is set which is not in the above list or is a species
name as defined by the chemistry of the model the will stop (unless it
starts with an X, XOH will not cause the model to crash).

#### Line 3

The third line gives information about which values should be
constrained and which ones allowed to run freely in the model
simulations. Each parameter is 15 characters long, separated by an space
(it is read in by the FORTRAN format statement ‘100000(a15,x))’. contain
either a "1" or a "0", indicating whether the parameter is to be
constrained in the model, with a ‘1’ indicating constraint on the
parameter and a ‘0’ indicating no constraint. Subsequent rows contain
the input data to the model, with concentrations in mixing ratio.

Where total NOx is to be constrained it is necessary to constrain either
NO or NO<sub>2</sub>, but not both. While the parameter ‘NOx’ must be included in
the Init\_cons.dat file for total NOx to be constrained, its values in
the file can be set to zero.

In order to constrain NOx the model will calculate a number every 24
hours by which the NO (or NO<sub>2</sub> if NO<sub>2</sub> is constrained in preference to
NO) must be multiplied so that its modelled value remains in agreement
with its observed value input into the model. All NOx species will
subsequently be multiplied by this value, and hence constrained by proxy
to NO (or NO<sub>2</sub>).

If neither *J* (O(^1^D)) nor *J* (NO<sub>2</sub>) are included in the input file,
clear-sky values will be calculated at the altitude in question
(determined from the pressure input) using TUV cross-sections at solar
zenith angles varying between 0 and 90 degrees in 5 degree steps. The
solar zenith angle (SZA) at which the observations were made is then
calculated from the observed latitude, longitude and time of day, and a
spline fit to the calculated *J*-values as a function of SZA used to
determine the appropriate *J*-value.

If *J* (O(^1^D)) or *J* (NO<sub>2</sub>) are present in the input file the model
will compare calculated *J*-values to their observed values and scale
all calculated values accordingly.

Unless otherwise stated in the input file the model assigns [CH~4~] =
1770 ppm, [H<sub>2</sub>] = 550 ppm, and an ozone column of 260 Dobsons.

### <a id="OutputFiles"/> Output Files


### Running the Model

To run the model type `model`

To log the diagnostic information produced by the model to a file type
`dsmacc > model_output`
