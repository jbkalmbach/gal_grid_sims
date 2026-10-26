#!/bin/bash
# Replicate what stack-prod.sh does, minus the exec into ipykernel.
# Edit the lines below to match the setup portion of stack-prod.sh.
set -eo pipefail
echo "wrapper: starting, args=$*"

source /opt/lsst/software/stack/loadLSST.bash ""
setup lsst_distrib

export OMP_NUM_THREADS=1

unset PYTHONSTARTUP

export PYTHONNOUSERSITE=' '

export DESC_GCR_SITE='nersc'

export RUBIN_SIM_DATA_DIR=/global/cfs/cdirs/lsst/groups/CO/rubin_sim

# skyCatalogs StarObject resolves its SED as $SIMS_SED_LIBRARY_DIR/<sed_filepath>.
# Our grid stars use the flat SED (flatSED/sed_flat.txt.gz); without this set,
# os.getenv returns None and os.path.join raises TypeError.  Galaxies use an
# inline tophat SED and never needed it, which is why galaxy-only runs worked.
export SIMS_SED_LIBRARY_DIR=/global/cfs/cdirs/descssim/imSim/lsst/data/sims_sed_library

if [ -n "$DESCPYTHONPATH" ]; then
    export PYTHONPATH="$DESCPYTHONPATH:$PYTHONPATH"
    echo "Including DESCPYTHONPATH: $DESCPYTHONPATH"
    echo "Wondering Why? DESCPYTHONPATH is likely set in your $HOME/.basrhc, $HOME/.bashrc.ext, or similar config script"
fi

if [ -n "$DESCSTACKUSERENV" ]; then
   conda activate $DESCSTACKUSERENV
   echo "Activated your DESCSTACKUSERENV: $DESCSTACKUSERENV"
    echo "Wondering Why? DESCSTACKUSERENV is likely set in your $HOME/.basrhc, $HOME/.bashrc.ext, or similar config script"
fi

python_ver_major=$(python -c 'import sys; print(sys.version_info.major)')
python_ver_minor=$(python -c 'import sys; print(sys.version_info.minor)')
export DESCPYTHONVER="python$python_ver_major.$python_ver_minor"

#export PYTHONPATH=$PYTHONPATH:$LSST_INST_DIR/$LSST_PYTHON_VER

if [ -n "$DESCSTACKUSERBASE" ]; then
    export PYTHONUSERBASE=$DESCSTACKUSERBASE
    unset PYTHONUSERSITE
    unset PYTHONNOUSERSITE
    export PATH=$PYTHONUSERBASE/bin:$PATH
    export PYTHONPATH="$PYTHONUSERBASE/lib/$DESCPYTHONVER/site-packages:$PYTHONPATH"
    echo "using DESCSTACKUSERBASE: $DESCSTACKUSERBASE"
fi


export PYTHONPATH=".:$PYTHONPATH"

echo "wrapper: setup done, exec-ing: $*"
exec "$@"
