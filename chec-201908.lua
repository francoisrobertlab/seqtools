help([[
For detailed instructions, go to:
   https://github.com/francoisrobertlab/chec-201908

This module loads the following modules and their requirements:
   - bwa/0.7.17
   - samtools/1.9
   - bedtools/2.27.1
   - kentutils/20180716
   - r/3.6.0
]])

whatis("Version: 1.0.0")
whatis("Keywords: ChEC-SEQ, Utility")
whatis("URL: https://github.com/francoisrobertlab/chec-201908")
whatis("Description: ChEC-SEQ analysis for Celia Jeronimo data")

always_load("nixpkgs/16.09")
always_load("gcc/7.3.0")
always_load("bwa/0.7.17")
always_load("samtools/1.9")
always_load("bedtools/2.27.1")
always_load("kentutils/20180716")
always_load("r/3.6.0")
always_load("vap")

local base = "~/projects/def-robertf/chec-201908"
prepend_path("PATH", pathJoin(base,"src/main/bash"))
setenv("CHEC_BASE", base)
setenv("CHEC_PATH", pathJoin(base,"src/main/python"))
setenv("CHEC_VENV", pathJoin(base,"venv/bin"))
