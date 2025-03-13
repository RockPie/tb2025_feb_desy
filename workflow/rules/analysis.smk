import glob
import json

config_files = glob.glob("data/config/PhaseScan_*.json")

# Extract all samples from filenames
samples = [file.split('/')[-1].replace('.json','') for file in config_files]

rule phasescan_all:
    input:
        expand("dump/200_PhaseScan/{sample}.root", sample=samples)

rule phasescantoa_all:
    input:
        expand("dump/201_PhaseScanToA/{sample}.root", sample=samples)

rule phasescan:
    input:
        config="data/config/{sample}.json",
        runfiles=lambda wildcards: json.load(open(f"data/config/{wildcards.sample}.json"))["Run Files"],
        script="build/200_PhaseScan"
    output:
        rootfile="dump/200_PhaseScan/{sample}.root"
    log:
        "logs/200_PhaseScan/{sample}.log"
    shell:
        "build/200_PhaseScan -f {input.config} -o {output.rootfile} &> {log}"

rule phasescantoa:
    input:
        config="data/config/{sample}.json",
        runfiles=lambda wildcards: json.load(open(f"data/config/{wildcards.sample}.json"))["Run Files"],
        script="build/201_PhaseScanToA"
    output:
        rootfile="dump/201_PhaseScanToA/{sample}.root"
    log:
        "logs/201_PhaseScanToA/{sample}.log"
    shell:
        "build/201_PhaseScanToA -f {input.config} -o {output.rootfile} &> {log}"