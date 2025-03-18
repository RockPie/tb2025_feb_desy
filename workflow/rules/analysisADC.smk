import glob
import json

config_files_eeemcal = glob.glob("data/DESY_2025/config/Pedestal_*.json")
config_files_focal = glob.glob("data/SPS_2024/config/Pedestal_*.json")

# Extract all samples from filenames
samples_eeemal = [file.split('/')[-1].replace('.json','') for file in config_files_eeemcal]
samples_focal = [file.split('/')[-1].replace('.json','') for file in config_files_focal]

rule Pedestal_EEEMCal:
    input:
        config="data/DESY_2025/config/{sample}.json",
        runfiles=lambda wildcards: json.load(open(f"data/DESY_2025/config/{wildcards.sample}.json"))["Run Files"],
        script="build/202_PedestalSub"
    output:
        rootfile="dump/202_PedestalSub/EEEMCal/{sample}.root"
    params:
        event = lambda wildcards: config.get("event", -1)
    log:
        "logs/202_PedestalSub/EEEMCal/{sample}.log"
    shell:
        "build/202_PedestalSub -f {input.config} -o {output.rootfile} -e {params.event} &> {log}"

rule Pedestal_FoCal:
    input:
        config="data/SPS_2024/config/{sample}.json",
        runfiles=lambda wildcards: json.load(open(f"data/SPS_2024/config/{wildcards.sample}.json"))["Run Files"],
        script="build/202_PedestalSub"
    output:
        rootfile="dump/202_PedestalSub/FoCal/{sample}.root"
    params:
        event = lambda wildcards: config.get("event", -1)
    log:
        "logs/202_PedestalSub/FoCal/{sample}.log"
    shell:
        "build/202_PedestalSub -f {input.config} -o {output.rootfile} -e {params.event} --focal &> {log}"

rule Pedestal_EEEMCal_all:
    input:
        expand("dump/202_PedestalSub/EEEMCal/{sample}.root", sample=samples_eeemal)

rule Pedestal_FoCal_all:
    input:
        expand("dump/202_PedestalSub/FoCal/{sample}.root", sample=samples_focal)