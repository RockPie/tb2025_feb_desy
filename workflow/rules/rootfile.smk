import os
from glob import glob
import re

# Directories
BUILD_DIR = "build"
DATA_DIR = "data/DESY_2025/data/beam"
DUMP_DIR = "dump/100_rootifier"
LOG_DIR = "logs/100_rootifier"

# Ensure directories exist
rule init:
    output:
        DUMP_DIR,
        LOG_DIR
    shell:
        "mkdir -p {output}"

# Find all RunXXX.h2g files
data_files = glob(os.path.join(DATA_DIR, "Run*.h2g"))

# Extract RunXXX identifiers
run_ids = [re.search(r'Run(\d+)\.h2g$', f).group(1) for f in data_files if re.search(r'Run(\d+)\.h2g$', f)]

# Rule to process data files using the compiled 100_rootifier script
rule process:
    input:
        BUILD_DIR + "/100_rootifier",
        DATA_DIR + "/Run{run_id}.h2g"
    output:
        DUMP_DIR + "/Run{run_id}_rootified.root"
    log:
        LOG_DIR + "/Run{run_id}.log"
    shell:
        """
        mkdir -p {LOG_DIR}
        {input[0]} -f {input[1]} -o {output} &> {log}
        """