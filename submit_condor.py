import os
import json
import argparse
import subprocess
from pathlib import Path
from analysis.utils import make_output_directory
from analysis.filesets.utils import divide_list, get_nano_version, fileset_checker


def move_proxy() -> str:
    """move x509 proxy file from /tmp to /afs/private. Returns the afs path"""
    try:
        subprocess.run("voms-proxy-info -exists -valid 0:20", shell=True, check=True)
    except subprocess.CalledProcessError:
        raise Exception(
            "VOMS proxy expired or non-existing: please run 'voms-proxy-init --voms cms'"
        )
    user = os.environ["USER"]
    x509_localpath = subprocess.check_output(
        "voms-proxy-info -path", shell=True, text=True
    ).strip()
    x509_path = (
        f"/afs/cern.ch/user/{user[0]}/{user}/private/{x509_localpath.split('/')[-1]}"
    )
    subprocess.run(["cp", x509_localpath, x509_path])
    return x509_path


def write_condor_file(template_file, condor_file, replacements):
    """Fill the condor submit template (condor/submit.sub): every key of `replacements`
    (CONDORDIR, BASEDIR, X509PATH, LOGDIR, JOBNAME, INPUTFILES, JOBNUM_FILE, JOBFLAVOR,
    IMAGE, MEMORY, NCPUS, MAXRETRIES) is replaced by its value. Keys are applied longest
    first so that no key is a substring hit of another one."""
    keys = sorted(replacements, key=len, reverse=True)
    with open(template_file) as template, open(condor_file, "w") as out:
        for line in template:
            for key in keys:
                line = line.replace(key, str(replacements[key]))
            out.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-w",
        "--workflow",
        dest="workflow",
        required=True,
        type=str,
        choices=[
            f.stem for f in (Path.cwd() / "analysis" / "workflows").glob("*.yaml")
        ],
        help="workflow to run",
    )
    parser.add_argument(
        "-y",
        "--year",
        dest="year",
        type=str,
        choices=[
            "2016preVFP",
            "2016postVFP",
            "2017",
            "2018",
            "2022preEE",
            "2022postEE",
            "2023preBPix",
            "2023postBPix",
            "2024",
        ],
        help="dataset year",
    )
    parser.add_argument(
        "-d",
        "--dataset",
        dest="dataset",
        type=str,
        help="dataset",
    )
    parser.add_argument(
        "--nfiles",
        dest="nfiles",
        type=int,
        default=5,
        help="number of root files per condor job (default 5)",
    )
    parser.add_argument(
        "--eos",
        action="store_true",
        help="Enable saving outputs to /eos",
    )
    parser.add_argument(
        "--submit",
        action="store_true",
        help="Enable Condor job submission. If not provided, it just builds condor files",
    )
    parser.add_argument(
        "--output_format",
        type=str,
        default="coffea",
        choices=["coffea", "root", "parquet"],
        help="format of output histogram",
    )
    parser.add_argument(
        "--image",
        dest="image",
        type=str,
        default="/cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-base-almalinux9:0.7.30-py3.10",
    )
    parser.add_argument("--jobflavor", dest="jobflavor", type=str, default="workday")
    parser.add_argument(
        "-m",
        "--memory",
        type=str,
        default="4000",
        help="Requested memory (in MB) for the condor job (default 4000)",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=50000,
        help="events per chunk (default 50000)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=2,
        help="worker processes per job = request_cpus (default 2)",
    )
    parser.add_argument(
        "--skipbadfiles",
        action="store_true",
        help="skip unreadable files/chunks instead of failing the job (see submit.py)",
    )
    parser.add_argument(
        "--max_retries",
        type=int,
        default=2,
        help="condor reruns of a job that exits with a non-zero code (default 2)",
    )
    args = parser.parse_args()

    # check if the fileset for the given year exists, generate it otherwise
    fileset_checker(samples=[args.dataset], year=args.year)

    print(f"Creating {args.workflow}-{args.year}-{args.dataset} condor file")
    jobname = f"{args.workflow}_{args.dataset}"

    # make condor and log directories
    condor_dir = Path.cwd() / "condor"
    job_dir = condor_dir / args.workflow / args.year / args.dataset
    if not job_dir.exists():
        job_dir.mkdir(parents=True, exist_ok=True)

    log_dir = condor_dir / "logs" / args.workflow / args.year / args.dataset
    if not log_dir.exists():
        log_dir.mkdir(parents=True, exist_ok=True)

    # save partitions json and jobnums to job directory
    jobnum_list = []
    partition_dataset = {}
    fileset_path = Path.cwd() / "analysis" / "filesets"

    nano_version = get_nano_version(args.year)
    with open(
        f"{fileset_path}/fileset_{args.year}_nanov{nano_version}_lxplus.json", "r"
    ) as f:
        root_files = json.load(f)[args.dataset]
    root_files_list = divide_list(root_files, args.nfiles)
    for i in range(len(root_files_list)):
        dataset_key = (
            f"{args.dataset}_{i+1}" if len(root_files_list) > 1 else args.dataset
        )
        partition_dataset[i + 1] = {dataset_key: root_files_list[i]}
        jobnum_list.append(i + 1)

    partition_file = job_dir / "partitions.json"
    with open(f"{partition_file}", "w") as json_file:
        json.dump(partition_dataset, json_file, indent=4)

    jobnum_file = job_dir / "jobnum.txt"
    with open(f"{jobnum_file}", "w") as f:
        print(*jobnum_list, sep="\n", file=f)

    # build and save arguments json
    args.output_path = make_output_directory(args)
    args.user = os.environ["USER"]
    args_file = job_dir / "arguments.json"
    with open(args_file, "w") as json_file:
        json.dump(vars(args), json_file, indent=4)

    # make condor file
    local_condor = f"{job_dir}/{jobname}.sub"
    write_condor_file(
        condor_dir / "submit.sub",
        local_condor,
        {
            "CONDORDIR": condor_dir,
            "BASEDIR": Path.cwd(),
            "X509PATH": move_proxy(),
            "LOGDIR": log_dir,
            "JOBNAME": jobname,
            "INPUTFILES": f"{partition_file},{jobnum_file},{args_file}",
            "JOBNUM_FILE": jobnum_file,
            "JOBFLAVOR": f'"{args.jobflavor}"',
            "IMAGE": f'"{args.image}"',
            "MEMORY": args.memory,
            "NCPUS": args.workers,
            "MAXRETRIES": args.max_retries,
        },
    )

    if args.submit:
        subprocess.run(["condor_submit", local_condor])
