import os

# XRootD client settings BEFORE coffea/uproot are imported (submit.sh exports the same;
# kept here so that a direct `python3 submit.py` behaves like a condor job).
# Retries: coffea 0.7.30 cannot retry a chunk with the futures executor (`retries` is not
# a FuturesExecutor/Runner field -> TypeError, and Runner.automatic_retries only retries
# with skipbadfiles=True, which skips OSError, i.e. the usual xrootd failure, without
# retrying). A failed job is therefore rerun as a whole by condor (max_retries in
# condor/submit.sub); the .coffea is written only at the end, so a rerun is idempotent.
os.environ.setdefault("XRD_RUNFORKHANDLER", "1")
os.environ.setdefault("XRD_NETWORKSTACK", "IPv4")

import sys
import json
import concurrent.futures
import argparse
import threading
import traceback
import faulthandler
from pathlib import Path
from coffea import processor
from coffea.util import save
from coffea.nanoevents import NanoAODSchema
from analysis.utils import write_root
from analysis.processors.base import BaseProcessor


def _run_and_save(args, output_location):
    """load the partition fileset, run the processor and save the output"""
    with open(args.partition_json) as f:
        partition_fileset = json.load(f)
    # our own pool instance: coffea then uses it as is, without the `with` block whose
    # shutdown(wait=True) deadlocks when a worker died, so the exception reaches __main__
    pool = concurrent.futures.ProcessPoolExecutor(max_workers=args.workers)
    out = processor.run_uproot_job(
        partition_fileset,
        treename="Events",
        processor_instance=BaseProcessor(
            workflow=args.workflow,
            year=args.year,
            output_format=args.output_format,
            output_location=output_location,
        ),
        executor=processor.futures_executor,
        # workers=2 (was 4): each worker holds a NanoAOD chunk in RAM; 4 pushed peak
        # to ~5.9 GB and OOM-truncated heavy samples at the 2-3 GB condor request.
        # 2 roughly halves peak (~3 GB) so jobs fit a modest request; ~2x slower/job.
        # skipbadfiles / xrootdtimeout go to coffea's Runner (supported in 0.7.30)
        executor_args={
            "schema": NanoAODSchema,
            "pool": pool,
            "workers": args.workers,
            "skipbadfiles": args.skipbadfiles,
            "xrootdtimeout": args.xrootdtimeout,
        },
        chunksize=args.chunksize,
    )
    pool.shutdown(wait=True)
    savepath = f"{args.output_path}/{args.dataset}"
    if args.output_format in ["coffea", "parquet"]:
        save(out, f"{savepath}.coffea")
    elif args.output_format == "root":
        write_root(out, savepath, args)


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
        "-d",
        "--dataset",
        dest="dataset",
        type=str,
        help="dataset",
    )
    parser.add_argument(
        "--partition_json",
        dest="partition_json",
        type=str,
        help="json with partition dataset",
    )
    parser.add_argument(
        "--output_path",
        dest="output_path",
        type=str,
        help="output path",
    )
    parser.add_argument(
        "--output_format",
        type=str,
        default="coffea",
        choices=["coffea", "root", "parquet"],
        help="format of output histogram",
    )
    parser.add_argument(
        "--eos",
        action="store_true",
        help="Enable saving outputs to /eos",
    )
    parser.add_argument(
        "--user",
        type=str,
        help="User name",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=50000,
        help="events per chunk (default 50000; coffea's own default is 100000)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=2,
        help="futures-executor worker processes (= request_cpus of the condor job)",
    )
    parser.add_argument(
        "--skipbadfiles",
        action="store_true",
        help="skip chunks whose file cannot be read (OSError / missing tree) instead of "
        "failing the job. The job then still writes its .coffea (counted as done by "
        "jobs_status.py) with fewer events: MC stays normalised (sumw is summed over the "
        "processed chunks only), data silently loses luminosity. Off by default.",
    )
    parser.add_argument(
        "--xrootdtimeout",
        type=int,
        default=600,
        help="uproot xrootd timeout per request in seconds (default 600; coffea's 60 made one "
        "slow remote read ('Operation expired') abort a whole job)",
    )
    parser.add_argument(
        "--watchdog",
        type=int,
        default=4 * 3600,
        help="abort the job (exit 3, thread dump in stderr) if it has not finished after this "
        "many seconds (default 4 h; a 5-file job takes < 1.5 h)",
    )
    args = parser.parse_args()

    # set output location (used when --output_format parquet)
    if args.eos:
        output_location = f"root://eosuser.cern.ch//eos/user/{args.user[0]}/{args.user}/higgscharm/outputs/"
    else:
        output_location = (
            f"/afs/cern.ch/user/{args.user[0]}/{args.user}/public/higgscharm/outputs/"
        )

    # Fail fast instead of hanging. When a futures-executor worker dies (seen 2026-10-08:
    # 'std::bad_alloc' in a worker after xrootd errors, cgroup far below its limit), coffea
    # raises BrokenProcessPool and then deadlocks in the pool shutdown of its `with` block
    # (a queue feeder thread blocked in pipe_write): the job never exits and condor never
    # retries it. _run_and_save passes its own pool, so the exception reaches here; any
    # exception -> traceback + os._exit(1) (skips the hanging shutdown); a watchdog
    # thread ends a job stuck anywhere else (exit 3). Condor reruns non-zero exits.
    def _watchdog_abort():
        print(f"submit.py watchdog: not finished after {args.watchdog} s, aborting (exit 3)",
              file=sys.stderr, flush=True)
        faulthandler.dump_traceback(file=sys.stderr, all_threads=True)
        sys.stderr.flush()
        os._exit(3)

    watchdog = threading.Timer(args.watchdog, _watchdog_abort)
    watchdog.daemon = True
    watchdog.start()
    try:
        _run_and_save(args, output_location)
    except BaseException:
        traceback.print_exc()
        sys.stdout.flush()
        sys.stderr.flush()
        os._exit(1)
    watchdog.cancel()
