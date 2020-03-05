#!/bin/bash

VENV="$HOME/seqtools-robertf-venv"
BASH="$VENV"/bash
SEQTOOLS="$VENV"/seqtools
SEQTOOLS_BASH="$SEQTOOLS"/bash

if [ "$1" == "clean" ]
then
    echo "Removing python virtual environment at $VENV"
    rm -R "$VENV"
fi
if [ ! -d "$VENV" ]
then
    echo "Creating python virtual environment at $VENV"
    python3 -m venv "$VENV"
fi
echo "Updating python libraries"
pip uninstall -y SeqTools
pip install git+https://git@github.com/francoisrobertlab/seqtools.git
echo "Updating bash scripts"
rm -R "$BASH"
mkdir "$BASH"
git clone https://github.com/francoisrobertlab/seqtools.git "$SEQTOOLS"
cp "$SEQTOOLS_BASH"/*.sh "$BASH"
rm -Rf "$SEQTOOLS"
