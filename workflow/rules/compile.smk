import os
from glob import glob

# Directories
SRC_DIR = "scripts"
BUILD_DIR = "build"
LOG_DIR = "logs/compile"

# Ensure the build directory exists
os.makedirs(BUILD_DIR, exist_ok=True)

# Find all CXX scripts
scripts = glob(os.path.join(SRC_DIR, "*.cxx"))

# # Rule to compile each script into an executable
rule compile_all:
    input:
        [os.path.join(BUILD_DIR, os.path.basename(script).replace(".cxx", "")) for script in scripts]

rule compile:
    input:
        SRC_DIR + "/{script}.cxx"
    output:
        BUILD_DIR + "/{script}"
    shell:
        """
        export PATH=/home/shihai/sw/root/root_install/bin:$PATH
        export LD_LIBRARY_PATH=/home/shihai/sw/root/root_install/lib:$LD_LIBRARY_PATH
        export ROOTSYS=/home/shihai/sw/root/root_install
        export PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH

        # Only re-run cmake if build dir doesn't exist
        if [ ! -d build ]; then
            mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release ..
        fi

        # Build only the modified target
        cd build && make {wildcards.script} -j
        """