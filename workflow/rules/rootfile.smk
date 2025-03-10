import os
from glob import glob
import re

# Directories
BUILD_DIR = "build"
DATA_DIR = "data/DESY_2025/data/beam"
DUMP_DIR = "dump"
LOG_DIR = "logs"

# Find all RunXXX.h2g files
data_files = glob(os.path.join(DATA_DIR, "Run*.h2g"))

# Extract RunXXX identifiers
run_ids = [re.search(r'Run(\d+)\.h2g$', f).group(1) for f in data_files if re.search(r'Run(\d+)\.h2g$', f)]

# Rule to process data files using the compiled 100_Rootifier script
rule Rootifier:
    input:
        BUILD_DIR + "/100_Rootifier",
        DATA_DIR + "/Run{run_id}.h2g"
    output:
        DUMP_DIR + "/100_Rootifier/Run{run_id}_rootified.root"
    log:
        LOG_DIR + "/100_Rootifier/Run{run_id}.log"
    shell:
        """
        mkdir -p {LOG_DIR}
        mkdir -p {DUMP_DIR}/100_Rootifier
        {input[0]} -f {input[1]} -o {output} &> {log}
        """

rule EventRecon:
    input:
        BUILD_DIR + "/101_EventRecon",
        DUMP_DIR + "/100_Rootifier/Run{run_id}_rootified.root"
    output:
        DUMP_DIR + "/101_EventRecon/Run{run_id}_recon.root"
    log:
        LOG_DIR + "/101_EventRecon/Run{run_id}_recon.log"
    shell:
        """
        mkdir -p {LOG_DIR}
        mkdir -p {DUMP_DIR}/101_EventRecon
        {input[0]} -f {input[1]} -o {output} &> {log}
        """

rule EventMatch:
    input:
        BUILD_DIR + "/102_EventMatch",
        DUMP_DIR + "/101_EventRecon/Run{run_id}_recon.root"
    output:
        DUMP_DIR + "/102_EventMatch/Run{run_id}_matched.root"
    log:
        LOG_DIR + "/102_EventMatch/Run{run_id}_match.log"
    shell:
        """
        mkdir -p {LOG_DIR}
        mkdir -p {DUMP_DIR}/102_EventMatch
        {input[0]} -f {input[1]} -o {output} &> {log}
        """

rule all:
    input:
        expand(DUMP_DIR + "/102_EventMatch/Run{run_id}_matched.root", run_id=run_ids)
    output:
        DUMP_DIR + "/all_matched.done"
    shell:
        """
        touch {output}
        """