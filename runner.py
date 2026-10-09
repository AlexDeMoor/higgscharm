import argparse
import subprocess
from pathlib import Path
from analysis.filesets.utils import (
    fileset_checker,
    get_datasets_to_run_over,
)


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
            "2024"
        ],
        help="dataset year",
    )
    parser.add_argument(
        "--nfiles",
        dest="nfiles",
        type=int,
        default=5,
        help="number of root files per condor job (default 5: ~45 min for the heaviest MC; a "
        "failed job is rerun as a whole, so short jobs lose less)",
    )
    parser.add_argument(
        "--submit",
        action="store_true",
        help="Enable Condor job submission. If not provided, it just builds condor files",
    )
    parser.add_argument(
        "--eos",
        action="store_true",
        help="Enable saving outputs to /eos",
    )
    parser.add_argument(
        "--output_format",
        type=str,
        default="coffea",
        choices=["coffea", "parquet"],
        help="format of output file",
    )
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
        help="events per chunk in the condor job (default 50000)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=2,
        help="worker processes per condor job, also its request_cpus (default 2)",
    )
    parser.add_argument(
        "--skipbadfiles",
        action="store_true",
        help="skip unreadable files/chunks instead of failing the job (see submit.py; off by default)",
    )
    parser.add_argument(
        "--max_retries",
        type=int,
        default=2,
        help="condor reruns of a job that exits with a non-zero code (default 2)",
    )
    args = parser.parse_args()

    # get datasets to run over for selected workflow and year
    datasets_to_run_over = get_datasets_to_run_over(args.workflow, args.year)

    # check if the input fileset for the given year exists, generate it otherwise
    fileset_checker(samples=datasets_to_run_over, year=args.year)

    # submit (or prepare) a job for each dataset using the given arguments
    cmd = ["python3", "submit_condor.py"]
    for dataset in datasets_to_run_over:
        cmd_args = [
            "--workflow",
            args.workflow,
            "--year",
            args.year,
            "--dataset",
            dataset,
            "--nfiles",
            str(args.nfiles),
            "--output_format",
            args.output_format,
            "--memory",
            args.memory,
            "--chunksize",
            str(args.chunksize),
            "--workers",
            str(args.workers),
            "--max_retries",
            str(args.max_retries),
        ]
        if args.skipbadfiles:
            cmd_args.append("--skipbadfiles")
        if args.submit:
            cmd_args.append("--submit")
        if args.eos:
            cmd_args.append("--eos")
        subprocess.run(cmd + cmd_args)
