#!/bin/bash

VENV="$HOME/seqtools-robertf-venv"
BASH="$VENV"/bash
SEQTOOLS="$VENV"/seqtools
SEQTOOLS_BASH="$SEQTOOLS"/bash
EMAIL="$USER_EMAIL"

if [ "$1" == "clean" ]
then
    echo "Removing python virtual environment at $VENV"
    rm -R "$VENV"
    exit 0
fi
if [[ ! "$EMAIL" =~ ^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,4}$ ]]
then
    echo "You must supply your email address as the first argument"
    exit 1
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
git clone --depth 1 https://github.com/francoisrobertlab/seqtools.git "$SEQTOOLS"
cp "$SEQTOOLS_BASH"/*.sh "$BASH"
find "$BASH" -type f -name "*.sh" -exec sed -i "s/christian\.poitras@ircm\.qc\.ca/$EMAIL/g" {} \;
rm -Rf "$SEQTOOLS"
