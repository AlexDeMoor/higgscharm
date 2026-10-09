"""
This script inspects Condor job outputs, identifies missing jobs, analyzes recent xrootd-related failures, and optionally regenerates input filesets with problematic sites blacklisted and resubmits only the missing jobs
"""

import os
import sys
import yaml
import json
import shutil
import argparse
import subprocess
from pathlib import Path
from datetime import datetime, timedelta
from collections import defaultdict

from analysis.utils import make_output_directory
from analysis.filesets.xrootd_sites import xroot_to_site
from analysis.filesets.utils import (
    divide_list,
    modify_site_list,
    extract_xrootd_errors,
    get_nano_version,
    get_dataset_config,
    get_datasets_to_run_over,
)


def parse_args():
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
        help="Workflow name (must correspond to a YAML file in analysis/workflows)",
    )
    parser.add_argument(
        "-y",
        "--year",
        type=str,
        required=True,
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
        help="Dataset year",
    )
    parser.add_argument(
        "--eos",
        action="store_true",
        help="Read job outputs from EOS instead of local storage",
    )
    parser.add_argument(
        "--output_format",
        type=str,
        default="coffea",
        choices=["coffea", "parquet"],
        help="Output file format of the produced histograms",
    )
    parser.add_argument(
        "--hours_ago",
        type=int,
        default=8,
        help="Only consider .err files modified within the last N hours",
    )
    parser.add_argument(
        "--reset",
        action="store_true",
        help="Start this workflow/year from scratch: MOVE its job outputs, condor job files "
        "and logs to timestamped '<dir>_reset_<time>' backups next to them, back up the fileset "
        "JSON and the site list and drop the workflow's non-signal fileset entries (the "
        "private 'era: signal' entries are kept), then resubmit every job with "
        "'runner.py --submit' (asks for confirmation; with --dry-run only prints the plan)",
    )
    parser.add_argument(
        "--dry-run",
        dest="dry_run",
        action="store_true",
        help="with --reset: print what would be moved and the runner.py command, change nothing",
    )
    # forwarded to runner.py by --reset (same defaults as runner.py)
    parser.add_argument("--nfiles", type=int, default=5, help="--reset: files per job")
    parser.add_argument("-m", "--memory", type=str, default="4000", help="--reset: request_memory [MB]")
    parser.add_argument("--chunksize", type=int, default=50000, help="--reset: events per chunk")
    parser.add_argument("--workers", type=int, default=2, help="--reset: workers (= request_cpus) per job")
    parser.add_argument("--max_retries", type=int, default=2, help="--reset: condor max_retries")
    parser.add_argument(
        "--skipbadfiles", action="store_true", help="--reset: skip unreadable files (see submit.py)"
    )
    parser.add_argument(
        "--single_fetch",
        action="store_true",
        help="Run fetch.py only once by merging all failing sites and datasets into a single group",
    )
    return parser.parse_args()


def get_jobs_info(job_dir, output_dir, log_dir, output_format, hours_ago=3):
    """
    Inspect the job directory structure and determine:
      1. Which jobs were expected to run (from jobnum.txt)
      2. Which jobs successfully produced output files
      3. Which jobs produced recent error logs

    error logs are collected *per dataset*, which allows later
    analysis and remediation to remain dataset-local rather than global.

    Parameters
    ----------
    job_dir : Path
        Directory containing one subdirectory per dataset with Condor configs.
    output_dir : Path
        Directory containing produced output files.
    log_dir : Path
        Directory containing Condor log files (.err, .out, .log).
    output_format : str
        File extension of produced outputs (e.g. "coffea", "parquet").
    hours_ago : int
        Only .err files modified within this many hours are considered relevant.

    Returns
    -------
    tuple
        jobnum : dict
            dataset -> list of expected job numbers (strings)
        jobnum_done : dict
            dataset -> list of job numbers that produced output
        error_files : dict
            dataset -> list of Path objects pointing to recent .err files
    """
    jobnum, jobnum_done = {}, {}
    error_files = {}

    for dataset_dir in job_dir.iterdir():
        if not dataset_dir.is_dir():
            continue

        dataset = dataset_dir.name
        jobnum_path = dataset_dir / "jobnum.txt"
        if not jobnum_path.exists():
            raise FileNotFoundError(
                f"Missing jobnum.txt for dataset '{dataset}'. "
                f"Expected at: {jobnum_path}"
            )

        # Read expected job numbers
        jobnum[dataset] = jobnum_path.read_text().splitlines()

        # Discover completed jobs by looking for output files.
        #
        # Completion is keyed off the per-job ".coffea" END-marker, NOT the
        # presence of output shards. submit.py writes "<dataset>_<jobnum>.coffea"
        # via save(out, ...) ONLY after run_uproot_job fully completes (and it does
        # so unconditionally, even for jobs that select zero events). This is the
        # only reliable "job finished" signal:
        #   * an OOM/killed job that wrote some parquet shards before dying leaves
        #     NO .coffea -> correctly flagged missing (a shard-presence check would
        #     wrongly mark it done and silently keep truncated output);
        #   * a job that legitimately selects zero events writes a .coffea but no
        #     shards -> correctly marked done (a shard-presence check would wrongly
        #     flag low-efficiency samples as missing).
        # Both parquet and coffea layouts write the same "<dataset>/<dataset>_<n>"
        # .coffea path, so the discovery is identical for both.
        marker_dir = output_dir / dataset
        jobnum_done[dataset] = []
        for n in jobnum[dataset]:
            # split sample: "<dataset>/<dataset>_<n>.coffea";
            # unsplit single-job fallback: "<dataset>/<dataset>.coffea"
            if (marker_dir / f"{dataset}_{n}.coffea").exists() or (
                marker_dir / f"{dataset}.coffea"
            ).exists():
                jobnum_done[dataset].append(n)

        # Collect recent error logs for this dataset
        x_hours_ago = datetime.now() - timedelta(hours=hours_ago)
        dataset_errs = []
        for err_file in (log_dir / dataset).glob("*.err"):
            if datetime.fromtimestamp(err_file.stat().st_mtime) > x_hours_ago:
                dataset_errs.append(err_file)

        if dataset_errs:
            error_files[dataset] = dataset_errs

    return jobnum, jobnum_done, error_files


def print_job_status(jobnum, jobnum_done):
    """
    Print a concise but informative summary of the job execution state.

    The function reports:
      - Total expected jobs
      - Total completed jobs
      - Total missing jobs
      - A YAML-formatted list of datasets with missing jobs

    Parameters
    ----------
    jobnum : dict
        dataset -> list of expected job numbers
    jobnum_done : dict
        dataset -> list of completed job numbers

    Returns
    -------
    tuple
        jobnum_missing : dict
            dataset -> set of missing job numbers
        datasets_with_missing : list
            List of dataset names with at least one missing job
    """
    jobnum_missing = {d: set(jobnum[d]) - set(jobnum_done.get(d, [])) for d in jobnum}
    n_expected = sum(len(v) for v in jobnum.values())
    n_done = sum(len(v) for v in jobnum_done.values())
    n_missing = sum(len(v) for v in jobnum_missing.values())

    print("------------------------------------------------------------")
    print("JOB STATUS SUMMARY")
    print("------------------------------------------------------------")
    print(f"Total jobs expected : {n_expected}")
    print(f"Total jobs finished : {n_done}")
    print(f"Total jobs missing  : {n_missing}")

    datasets_with_missing = [d for d in jobnum_missing if jobnum_missing[d]]
    if n_missing:
        print("------------------------------------------------------------")
        print(
            f"Datasets with missing jobs ({len(datasets_with_missing)}/{len(jobnum)}):"
        )
        print("------------------------------------------------------------")
        print(
            yaml.dump(
                datasets_with_missing,
                default_flow_style=False,
                sort_keys=False,
                indent=2,
            )
        )

    return jobnum_missing, datasets_with_missing


def analyze_xrootd_errors_by_dataset(error_files):
    """
    Analyze xrootd-related failures on a per-dataset basis.

    For each dataset, this function:
      1. Parses its recent .err files
      2. Extracts xrootd endpoints or error patterns
      3. Maps them to physical sites using `xroot_to_site`
      4. Produces a set of failing sites per dataset

    This dataset-level granularity avoids the pathological behavior
    of blacklisting a site globally when it only affects one dataset.

    Parameters
    ----------
    error_files : dict
        dataset -> list of Path objects pointing to .err files

    Returns
    -------
    dict
        dataset -> set of sites that exhibited xrootd errors
    """
    dataset_sites = {}

    for dataset, err_files in error_files.items():
        xrootd_errs = extract_xrootd_errors(err_files)

        # Translate xrootd endpoints to site names where possible
        sites = {xroot_to_site[err] for err in xrootd_errs if err in xroot_to_site}
        if sites:
            dataset_sites[dataset] = sites

        # Warn explicitly if some endpoints could not be mapped
        for err in xrootd_errs:
            if err not in xroot_to_site:
                print(
                    f"[{dataset}] Could not map xrootd endpoint '{err}' to a site name"
                )

    if not dataset_sites:
        print(f"No xrootd site failures detected in recent error logs.\n")
    return dataset_sites


def group_datasets_by_sites(dataset_sites):
    """
    Group datasets by identical sets of failing sites.

    Many datasets often fail on exactly the same sites.
    Rather than regenerating input filesets separately
    for each dataset, we group datasets by their failing-site
    signature and regenerate filesets once per group.

    Example:
        Input:
            {
              "A": {"T2_US_Florida", "T2_FR_GRIF"},
              "B": {"T2_US_Florida", "T2_FR_GRIF"},
              "C": {"T2_IT_Pisa"},
            }

        Output:
            {
              frozenset({"T2_US_Florida", "T2_FR_GRIF"}): ["A", "B"],
              frozenset({"T2_IT_Pisa"}): ["C"],
            }

    Parameters
    ----------
    dataset_sites : dict
        dataset -> set of failing sites

    Returns
    -------
    dict
        frozenset(failing sites) -> list of datasets
    """
    groups = defaultdict(list)
    for dataset, sites in dataset_sites.items():
        groups[frozenset(sites)].append(dataset)
    return groups


def update_input_filesets_for_group(
    sites_to_blacklist, year, fileset_dir, job_dir, datasets
):
    """
    Regenerate input filesets for a group of datasets that share the same
    failing sites.

    The procedure is:
      1. Reset all sites to "white" (allowed)
      2. Blacklist only the sites known to be problematic for this group
      3. Run fetch.py once for all datasets in the group
      4. Regenerate partition files (partitions.json) for each dataset

    This minimizes expensive calls to fetch.py while still preserving
    dataset-level correctness.

    Parameters
    ----------
    sites_to_blacklist : iterable
        Collection of site names to blacklist for this group
    year : str
        Dataset year
    fileset_dir : Path
        Directory containing JSON filesets
    job_dir : Path
        Directory containing Condor job files
    datasets : list
        List of dataset names to regenerate
    """
    print("------------------------------------------------------------")
    print("REGENERATING FILESETS FOR DATASET GROUP")
    print("------------------------------------------------------------")
    print(f"Datasets           : {datasets}")
    print(f"Blacklisted sites  : {sorted(sites_to_blacklist)}")

    # First, reset all sites to "white" to ensure no cross-group contamination.
    # This guarantees that each group is regenerated with exactly its own
    # blacklist, and nothing else.
    for site in xroot_to_site.values():
        modify_site_list(year, site, "white")

    # Apply the blacklist for this specific group.
    for site in sites_to_blacklist:
        modify_site_list(year, site, "black")

    # Run fetch.py once for the entire dataset group.
    samples_str = " ".join(datasets)
    print("Running fetch.py for this group...")
    subprocess.run(
        ["python3", "fetch.py", "--year", year, "--samples", samples_str],
        check=True,
    )

    # Load the regenerated filesets
    nano_version = get_nano_version(year)
    fileset_path = fileset_dir / f"fileset_{year}_nanov{nano_version}_lxplus.json"
    all_filesets = json.loads(fileset_path.read_text())

    # Regenerate partitions.json for each dataset in the group
    for dataset in datasets:
        if dataset not in all_filesets:
            print(f"[{dataset}] Not found in fileset JSON — skipping")
            continue

        root_files = all_filesets[dataset]
        args_json = job_dir / dataset / "arguments.json"
        if not args_json.exists():
            print(f"[{dataset}] Missing arguments.json — cannot repartition")
            continue

        nfiles = json.loads(args_json.read_text())["nfiles"]
        root_files_list = divide_list(root_files, nfiles)

        partition_dataset = {
            i
            + 1: {(f"{dataset}_{i+1}" if len(root_files_list) > 1 else dataset): chunk}
            for i, chunk in enumerate(root_files_list)
        }

        partition_file = job_dir / dataset / "partitions.json"
        with open(partition_file, "w") as json_file:
            json.dump(partition_dataset, json_file, indent=4)

    print("Fileset regeneration for this group completed.\n")


def resubmit_jobs(job_dir, jobnum_missing, datasets_with_missing_jobs, workflow, year):
    """
    Prepare and resubmit only the jobs that are missing output files.

    For each affected dataset:
      1. Write a missing.txt file listing the missing job numbers
      2. Patch the Condor submission file to use missing.txt instead of jobnum.txt
      3. Submit the modified job description to Condor

    A backup of the original submission file (*_all.sub) is created once
    per dataset to avoid destructive overwrites.

    Parameters
    ----------
    job_dir : Path
        Directory containing Condor job files
    jobnum_missing : dict
        dataset -> set of missing job numbers
    datasets_with_missing_jobs : list
        List of dataset names with missing jobs
    workflow : str
        Workflow name
    year : str
        Dataset year
    """
    print("------------------------------------------------------------")
    print("RESUBMITTING MISSING JOBS")
    print("------------------------------------------------------------")

    to_resubmit = []
    for dataset in datasets_with_missing_jobs:
        missing_jobs = sorted(jobnum_missing[dataset])
        missing_file = job_dir / dataset / "missing.txt"
        with open(missing_file, "w") as f:
            print(*missing_jobs, sep="\n", file=f)

        condor_file = job_dir / dataset / f"{workflow}_{dataset}.sub"
        condor_backup = condor_file.with_name(condor_file.stem + "_all.sub")
        if not condor_backup.exists():
            subprocess.run(["cp", str(condor_file), str(condor_backup)], check=True)

        submit_text = condor_file.read_text().replace("jobnum.txt", "missing.txt")
        condor_file.write_text(submit_text)
        to_resubmit.append(str(condor_file))

    for submit_file in to_resubmit:
        print(f"Submitting {submit_file}")
        subprocess.run(["condor_submit", submit_file], check=True)

    print("Job resubmission completed.\n")


def reset_workflow(
    workflow,
    year,
    base_dir,
    output_dirs,
    datasets,
    dataset_config,
    runner_cmd,
    timestamp,
    run=subprocess.run,
    dry_run=False,
):
    """Reset one workflow/year (jobs_status.py --reset). Nothing is deleted:

      1. output directories (`output_dirs`, i.e. outputs/<workflow>/<year>: per-job .coffea,
         parquet shards, sumw_records and the postprocessed files, which would be stale),
         condor/<workflow>/<year> and condor/logs/<workflow>/<year> are RENAMED to
         '<dir>_reset_<timestamp>' next to them (outside every glob of run_postprocess.py
         and jobs_status.py);
      2. analysis/filesets/<year>_sites.yaml (the site black/white list that the xrootd
         remediation edits; shared by all workflows of the year) is renamed to
         '.bak_reset_<timestamp>' so that fetch.py builds a fresh one -- only when some
         datasets are refetched (a signal-only workflow never runs fetch.py);
      3. the fileset JSON is copied to '.bak_reset_<timestamp>' and the workflow's datasets
         are removed from it EXCEPT the 'era: signal' ones (private samples: plain fetch.py
         skips them, they are only rebuilt by `fetch.py --onlysignal`); runner.py's
         fileset_checker then refetches the removed ones (make_filesets updates the JSON in
         place, other workflows' entries are untouched);
      4. `runner_cmd` (runner.py ... --submit) is run.

    Returns the list of (source, destination) moves."""
    base_dir = Path(base_dir)
    nano_version = get_nano_version(year)
    fileset_file = base_dir / "analysis" / "filesets" / f"fileset_{year}_nanov{nano_version}_lxplus.json"
    sites_file = base_dir / "analysis" / "filesets" / f"{year}_sites.yaml"

    signal = [d for d in datasets if (dataset_config.get(d) or {}).get("era") == "signal"]
    to_refetch = [d for d in datasets if d not in signal]

    moves = []
    for path in [Path(p) for p in output_dirs] + [
        base_dir / "condor" / workflow / year,
        base_dir / "condor" / "logs" / workflow / year,
    ]:
        if path.exists():
            moves.append((path, path.with_name(f"{path.name}_reset_{timestamp}")))
    if sites_file.exists() and to_refetch:
        # only when fetch.py will run (it rebuilds the list; the xrootd remediation needs it).
        # NB: the site list is shared by all the workflows of this year
        moves.append((sites_file, sites_file.with_name(f"{sites_file.name}.bak_reset_{timestamp}")))
    filesets = json.loads(fileset_file.read_text()) if fileset_file.exists() else {}
    missing_signal = [d for d in signal if d not in filesets]
    if missing_signal:
        # runner.py's fetch.py call cannot build them (make_filesets skips era: signal)
        raise RuntimeError(
            f"signal datasets {missing_signal} are not in {fileset_file}: run "
            f"'python3 fetch.py -y {year} --onlysignal' first"
        )

    print("------------------------------------------------------------")
    print(f"RESET {workflow} / {year}" + ("  (dry run: nothing is changed)" if dry_run else ""))
    print("------------------------------------------------------------")
    for src, dst in moves:
        print(f"move   {src}\n    -> {dst}")
    if fileset_file.exists():
        print(f"backup {fileset_file} -> {fileset_file.name}.bak_reset_{timestamp}")
        print(f"fileset entries kept (era: signal): {signal}")
        print(f"fileset entries dropped (refetched by runner.py): {to_refetch}")
    print("resubmit: " + " ".join(runner_cmd))
    if dry_run:
        return moves

    for src, dst in moves:
        if dst.exists():
            raise FileExistsError(f"reset backup {dst} already exists")
        shutil.move(str(src), str(dst))

    if fileset_file.exists():
        shutil.copy2(fileset_file, fileset_file.with_name(f"{fileset_file.name}.bak_reset_{timestamp}"))
        for dataset in to_refetch:
            filesets.pop(dataset, None)
        fileset_file.write_text(json.dumps(filesets, indent=4, sort_keys=True))

    run(runner_cmd, check=True)
    return moves


def main():
    args = parse_args()

    # Optional full reset (backups, no deletion; see reset_workflow). Ends the script: the
    # job status of a freshly submitted workflow is meaningless (every job "missing")
    if args.reset:
        base_dir = Path.cwd()
        output_dirs = [Path(make_output_directory(args))]
        if args.output_format == "parquet" and not args.eos:
            # submit.py writes the parquet shards to AFS when --eos is not given
            user = os.environ["USER"]
            output_dirs.append(
                Path(f"/afs/cern.ch/user/{user[0]}/{user}/public/higgscharm/outputs")
                / args.workflow
                / args.year
            )
        runner_cmd = [
            sys.executable, "runner.py",
            "-w", args.workflow,
            "-y", args.year,
            "--output_format", args.output_format,
            "--nfiles", str(args.nfiles),
            "--memory", args.memory,
            "--chunksize", str(args.chunksize),
            "--workers", str(args.workers),
            "--max_retries", str(args.max_retries),
            "--submit",
        ]
        if args.eos:
            runner_cmd.append("--eos")
        if args.skipbadfiles:
            runner_cmd.append("--skipbadfiles")
        reset_args = dict(
            workflow=args.workflow,
            year=args.year,
            base_dir=base_dir,
            output_dirs=output_dirs,
            datasets=get_datasets_to_run_over(args.workflow, args.year),
            dataset_config=get_dataset_config(args.year),
            runner_cmd=runner_cmd,
            timestamp=datetime.now().strftime("%Y%m%d_%H%M%S"),
        )
        reset_workflow(**reset_args, dry_run=True)
        if args.dry_run:
            return
        if input("Proceed with the reset and resubmit? (y/n): ").lower() not in ["y", "yes"]:
            print("Reset aborted, nothing changed.")
            return
        reset_workflow(**reset_args)
        print("Reset complete: jobs resubmitted.\n")
        return

    # Directory setup
    output_dir = Path(make_output_directory(args))
    base_dir = Path.cwd()
    condor_dir = base_dir / "condor"
    job_dir = condor_dir / args.workflow / args.year
    log_dir = condor_dir / "logs" / args.workflow / args.year
    fileset_dir = base_dir / "analysis" / "filesets"

    # Discover job state
    jobnum, jobnum_done, error_files = get_jobs_info(
        job_dir, output_dir, log_dir, args.output_format, args.hours_ago
    )

    jobnum_missing, datasets_with_missing_jobs = print_job_status(jobnum, jobnum_done)

    if not datasets_with_missing_jobs:
        return

    # Analyze xrootd errors per dataset
    dataset_sites = analyze_xrootd_errors_by_dataset(error_files)

    # Restrict to datasets that actually have missing jobs
    dataset_sites = {
        d: sites
        for d, sites in dataset_sites.items()
        if d in datasets_with_missing_jobs
    }

    if not dataset_sites:
        print("No xrootd-related failures detected for datasets with missing jobs.")
        print("You may want to inspect other error modes or resubmit directly.\n")
    else:
        # Group datasets by identical failing-site sets
        groups = group_datasets_by_sites(dataset_sites)
        if args.single_fetch:
            all_sites = set().union(*dataset_sites.values())
            all_datasets = sorted(dataset_sites.keys())
            groups = {frozenset(all_sites): all_datasets}
        else:
            groups = group_datasets_by_sites(dataset_sites)

        # Optional fileset regeneration
        if input("Regenerate input filesets for these groups? (y/n): ").lower() in [
            "y",
            "yes",
        ]:
            for sites, datasets in groups.items():
                update_input_filesets_for_group(
                    sites, args.year, fileset_dir, job_dir, datasets
                )

    # Optional job resubmission
    if input("Resubmit missing jobs now? (y/n): ").lower() in ["y", "yes"]:
        resubmit_jobs(
            job_dir,
            jobnum_missing,
            datasets_with_missing_jobs,
            args.workflow,
            args.year,
        )


if __name__ == "__main__":
    main()
