import subprocess

remote_dir = "/eos/experiment/alice/focal/tb2024_Sep_SPSH2/hcal/FoCalH"
ssh_user = "sjia"
ssh_server = "lxplus.cern.ch"
ssh_pass = "workflow/rules/scp_pass.txt"

cmd = f"sshpass -f {ssh_pass} ssh {ssh_user}@{ssh_server} 'ls {remote_dir}/*.h2g'"
result = subprocess.check_output(cmd, shell=True, text=True)

RUN_IDS = [
    int(line.strip().split("Run")[1].split(".h2g")[0])
    for line in result.strip().split("\n") if line.strip()
]

rule download_focal_all:
    input:
        expand("data/SPS_2024/data/beam/Run{run_id}.h2g", run_id=RUN_IDS)
    output:
        "data/SPS_2024/data/beam/all.done"
    shell:
        """
        touch {output}
        """

rule download_focal_run_file:
    output:
        "data/SPS_2024/data/beam/Run{run_id}.h2g"
    params:
        user="sjia",
        server="lxplus.cern.ch",
        password="workflow/rules/scp_pass.txt",
        remote_path=lambda wildcards: f"/eos/experiment/alice/focal/tb2024_Sep_SPSH2/hcal/FoCalH{'_1000' if int(wildcards.run_id) < 500 else ''}/Run{wildcards.run_id:0>3}.h2g"
    shell:
        """
        sshpass -f {params.password} scp {params.user}@{params.server}:{params.remote_path} {output}
        """