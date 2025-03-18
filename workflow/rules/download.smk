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