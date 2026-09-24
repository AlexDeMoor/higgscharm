"""Run a workflow locally (no condor) over one dataset of the input fileset.

Writes exactly what a condor job writes -- <outputs>/<workflow>/<year>/<dataset>/<dataset>.coffea
(plus the parquet shards when --output_format parquet) -- so run_postprocess.py and
jobs_status.py work unchanged. Meant for the private signal samples and quick tests:

    python3 run_local.py -w htautau_mutau -y 2024 -d HplusBottom_HtoTauTau --eos --workers 8
    python3 run_local.py -w htautau_mutau -y 2024 -d HplusBottom_HtoTauTau --files /path/to/file.root --max-chunks 1
"""
import os

# XRootD client settings BEFORE coffea/uproot are imported: the futures executor forks
# the workers after the main process has already opened the files (preprocessing), and
# without the fork handler the forked xrootd client stalls (same as condor/submit.sh)
os.environ.setdefault("XRD_RUNFORKHANDLER", "1")
os.environ.setdefault("XRD_NETWORKSTACK", "IPv4")

import json
import time
import argparse
from pathlib import Path
from coffea import processor
from coffea.util import save
from coffea.nanoevents import NanoAODSchema
from analysis.utils import make_output_directory
from analysis.filesets.utils import get_nano_version
from analysis.processors.base import BaseProcessor


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-w", "--workflow", required=True, type=str,
        choices=[f.stem for f in (Path.cwd() / "analysis" / "workflows").glob("*.yaml")],
    )
    parser.add_argument("-y", "--year", required=True, type=str)
    parser.add_argument("-d", "--dataset", required=True, type=str, help="dataset name (fileset key)")
    parser.add_argument("--files", nargs="*", default=None, help="explicit ROOT files (bypass the fileset)")
    parser.add_argument("--max-files", type=int, default=None, help="use only the first N files")
    parser.add_argument("--max-chunks", type=int, default=None, help="process only N chunks (test)")
    parser.add_argument("--chunksize", type=int, default=50000)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--output_format", type=str, default="coffea", choices=["coffea", "parquet"])
    parser.add_argument("--eos", action="store_true", help="write outputs to /eos/user/<u>/higgscharm/outputs")
    parser.add_argument("--tag", type=str, default="", help="optional suffix for the output file name")
    parser.add_argument("--iterative", action="store_true", help="single-process executor (full tracebacks, debugging)")
    parser.add_argument("--stage-dir", type=str, default=None,
                        help="copy the input files to this local directory (xrdcp) and read them from there; "
                             "columnar reads over xrootd are latency-bound (~10 MB/s per node), local disk is not")
    parser.add_argument("--stage-only", action="store_true", help="only stage the files, do not process")
    parser.add_argument("--stage-streams", type=int, default=4, help="parallel xrdcp copies")
    parser.add_argument("--clean-staged", action="store_true", help="delete the staged copies at the end")
    return parser.parse_args()


def stage_files(files, stage_dir, streams=4):
    """copy (xrdcp) the files into stage_dir unless already there; returns local paths"""
    import subprocess
    from concurrent.futures import ThreadPoolExecutor

    stage_dir = Path(stage_dir)
    stage_dir.mkdir(parents=True, exist_ok=True)

    def to_url(f):
        return f if f.startswith("root://") else f"root://eosuser.cern.ch/{f}" if f.startswith("/eos/") else f

    def copy(f):
        # basename + short hash of the full source path: DAS files of different blocks
        # share basenames (nano_1.root, ...)
        import hashlib
        digest = hashlib.sha1(f.encode()).hexdigest()[:8]
        local = stage_dir / f"{Path(f).stem}_{digest}{Path(f).suffix}"
        if local.exists():
            return str(local)
        tmp = local.with_suffix(".part")
        cmd = ["xrdcp", "-f", "-s", to_url(f), str(tmp)]
        for attempt in range(3):
            if subprocess.run(cmd).returncode == 0:
                tmp.rename(local)
                return str(local)
        raise RuntimeError(f"staging failed for {f}")

    t0 = time.time()
    with ThreadPoolExecutor(max_workers=streams) as pool:
        local_files = list(pool.map(copy, files))
    size = sum(Path(f).stat().st_size for f in local_files) / 1e9
    print(f"staged {len(local_files)} files ({size:.1f} GB) in {stage_dir} [{time.time() - t0:.0f} s]")
    return local_files


if __name__ == "__main__":
    args = parse_args()
    args.user = os.environ["USER"]

    if args.files:
        files = args.files
    else:
        nano_version = get_nano_version(args.year)
        fileset_file = (
            Path.cwd() / "analysis" / "filesets"
            / f"fileset_{args.year}_nanov{nano_version}_lxplus.json"
        )
        with open(fileset_file) as f:
            files = json.load(f)[args.dataset]
    if args.max_files:
        files = files[: args.max_files]
    if args.stage_dir:
        files = stage_files(files, args.stage_dir, args.stage_streams)
        if args.stage_only:
            raise SystemExit(0)
    print(f"{args.workflow} / {args.year} / {args.dataset}: {len(files)} files, "
          f"{args.workers} workers, chunksize {args.chunksize}")

    # same output roots as submit.py / submit_condor.py
    output_path = make_output_directory(args)
    if args.eos:
        output_location = f"root://eosuser.cern.ch//eos/user/{args.user[0]}/{args.user}/higgscharm/outputs/"
    else:
        output_location = f"/afs/cern.ch/user/{args.user[0]}/{args.user}/public/higgscharm/outputs/"

    t0 = time.time()
    out = processor.run_uproot_job(
        {args.dataset: files},
        treename="Events",
        processor_instance=BaseProcessor(
            workflow=args.workflow,
            year=args.year,
            output_format=args.output_format,
            output_location=output_location,
        ),
        executor=processor.iterative_executor if args.iterative else processor.futures_executor,
        executor_args=(
            {"schema": NanoAODSchema}
            if args.iterative
            else {"schema": NanoAODSchema, "workers": args.workers}
        ),
        chunksize=args.chunksize,
        maxchunks=args.max_chunks,
    )
    savepath = f"{output_path}/{args.dataset}{args.tag}.coffea"
    save(out, savepath)
    if args.stage_dir and args.clean_staged:
        for f in files:
            Path(f).unlink(missing_ok=True)
    print(f"saved {savepath}  ({time.time() - t0:.0f} s)")
    if "metadata" in out:
        print("sumw:", out["metadata"].get("sumw"))
        for key, val in out["metadata"].items():
            if isinstance(val, dict) and "cutflow" in val:
                print(f"[{key}]", {k: round(float(v), 3) for k, v in val["cutflow"].items()})
