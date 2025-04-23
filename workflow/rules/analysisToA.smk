import glob
import json

rule ToA_DNL_Calibration_FoCal:
    input:
        root="dump/102_EventMatch/FoCal/Run{sample}_matched.root",
        script="build/207_ToACalib",
    output:
        rootfile="dump/207_ToACalib/FoCal/Run{sample}_toa_calib.root",
        jsonfile="dump/207_ToACalib/FoCal/Run{sample}_toa_calib_fineDNL.json"
    params:
        event = lambda wildcards: config.get("event", -1)
    log:
        "logs/204_ToA_DNL_Calibration/EEEMCal/{sample}.log"
    shell:
        "{input.script} -f {input.root} -o {output.rootfile} -e {params.event} --focal &> {log}"