#!/bin/bash
# Runs KPP on the DSMACC source, then checks each generated file for
# differences with the existing source (if any) before moving them down to
# src/ to reduce the amount of recompilation required when making small
# changes to the code.

# move the given file to .. if its counterpart doesn't exist or is different
move_changed() {
    newsrc="$1"
    dest="../src/$1"

    if [ -e "$dest" ]; then # see if file has changed
        diff -q <(sed '/^! Time.*/ D' "$newsrc") \
                <(sed '/^! Time.*/ D' "$dest")
        if [ $? -ne 0 ]; then
            echo "$newsrc changed"
            mv "$newsrc" "$dest"
        fi
    else # doesn't exist, must move
        mv -v "$newsrc" "$dest"
    fi
}

${KPP_HOME}/bin/kpp ../dsmacc.kpp dsmacc

if [ $? -ne 0 ]; then
    echo "KPP failed to process the DSMACC source"
    exit 1
fi

for f in *.f90; do
    move_changed "$f"
done
echo

# no need to keep remaining files
rm *.f90
