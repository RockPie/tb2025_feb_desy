import glob
import json

config_files_eeemcal = glob.glob("data/DESY_2025/config/PhaseScan_*.json")
config_files_focal = glob.glob("data/SPS_2024/config/PhaseScan_*.json")

# Extract all samples from filenames
samples_eeemal = [file.split('/')[-1].replace('.json','') for file in config_files_eeemcal]
samples_focal = [file.split('/')[-1].replace('.json','') for file in config_files_focal]

rule phasescan_EEEMCal:
    input:
        config="data/DESY_2025/config/{sample}.json",
        runfiles=lambda wildcards: json.load(open(f"data/DESY_2025/config/{wildcards.sample}.json"))["Run Files"],
        script="build/200_PhaseScan"
    output:
        rootfile="dump/200_PhaseScan/EEEMCal/{sample}.root"
    params:
        event = lambda wildcards: config.get("event", -1)
    log:
        "logs/200_PhaseScan/EEEMCal/{sample}.log"
    shell:
        "build/200_PhaseScan -f {input.config} -o {output.rootfile} -e {params.event} &> {log}"

rule phasescan_FoCal:
    input:
        config="data/SPS_2024/config/{sample}.json",
        runfiles=lambda wildcards: json.load(open(f"data/SPS_2024/config/{wildcards.sample}.json"))["Run Files"],
        script="build/200_PhaseScan"
    output:
        rootfile="dump/200_PhaseScan/FoCal/{sample}.root"
    params:
        event = lambda wildcards: config.get("event", -1)
    log:
        "logs/200_PhaseScan/FoCal/{sample}.log"
    shell:
        "build/200_PhaseScan -f {input.config} -o {output.rootfile} -e {params.event} &> {log}"

rule phasescan_EEEMCal_all:
    input:
        expand("dump/201_PhaseScanToA/EEEMCal/{sample}.root", sample=samples_eeemal)

rule phasescan_FoCal_all:
    input:
        expand("dump/201_PhaseScanToA/FoCal/{sample}.root", sample=samples_focal)

rule phasescantoa_EEEMCal_all:
    input:
        expand("dump/201_PhaseScanToA/EEEMCal/{sample}.root", sample=samples_eeemal)

rule phasescantoa_FoCal_all:
    input:
        expand("dump/201_PhaseScanToA/FoCal/{sample}.root", sample=samples_focal)

rule phasescantoa_EEEMCal:
    input:
        config="data/DESY_2025/config/{sample}.json",
        runfiles=lambda wildcards: json.load(open(f"data/DESY_2025/config/{wildcards.sample}.json"))["Run Files"],
        script="build/201_PhaseScanToA"
    output:
        rootfile="dump/201_PhaseScanToA/EEEMCal/{sample}.root"
    params:
        event = lambda wildcards: config.get("event", -1)
    log:
        "logs/201_PhaseScanToA/EEEMCal/{sample}.log"
    shell:
        "build/201_PhaseScanToA -f {input.config} -o {output.rootfile} -e {params.event} &> {log}"

rule phasescantoa_FoCal:
    input:
        config="data/SPS_2024/config/{sample}.json",
        runfiles=lambda wildcards: json.load(open(f"data/SPS_2024/config/{wildcards.sample}.json"))["Run Files"],
        script="build/201_PhaseScanToA"
    output:
        rootfile="dump/201_PhaseScanToA/FoCal/{sample}.root"
    params:
        event = lambda wildcards: config.get("event", -1)
    log:
        "logs/201_PhaseScanToA/FoCal/{sample}.log"
    shell:
        "build/201_PhaseScanToA -f {input.config} -o {output.rootfile} -e {params.event} --focal &> {log}"