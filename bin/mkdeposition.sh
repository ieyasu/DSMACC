#!/bin/bash
# Reads (in)organic.kpp files and writes a corresponding deposition scheme
# to depos.kpp.

extract_species() {
    # pulls out species from $1 between #DEFVAR and the next '#' or '{' line
    sed -e '0,/#DEFVAR/ d' -e '/^[^#{]*[#{]/,$ d' < "$1" |
        grep '^ *[A-Za-z0-9]\+ *=[^;]\+;' | # remove anything not equation-like
        sed 's/^ *\([A-Za-z0-9]\+\) *=.\+; */\1/' # extract species name  
}

# ---

# name of deposition file
dep=depos.kpp

#echo '#DEFVAR'           > $dep
#echo 'DUMMY = IGNORE ;' >> $dep
echo '#EQUATIONS'       > $dep

specs=$( extract_species inorganic.kpp; extract_species organic.kpp; )

n=1
for spec in $specs; do
    echo "{${n}.} ${spec} = DUMMY : DEPOSITION ;" >>$dep
    n=$(expr $n + 1)
done
