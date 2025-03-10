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

# Rule to compile each script into an executable
rule all:
    input:
        [os.path.join(BUILD_DIR, os.path.basename(script).replace(".cxx", "")) for script in scripts]

rule compile:
    input:
        SRC_DIR + "/{script}.cxx"
    output:
        BUILD_DIR + "/{script}"
    shell:
        "cd build && cmake .. && make {wildcards.script}"