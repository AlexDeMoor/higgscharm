import glob
import json
import argparse
import subprocess
from pathlib import Path
from analysis.filesets.utils import get_nano_version, get_dataset_config


def get_signal_fileset(query):
    query = f"file dataset=/{query}"
    farray = subprocess.run(
        ["dasgoclient", f"-query={query} instance=prod/phys03"],
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    stdout = farray.stdout
    stdout_array = stdout.split("\n")
    stdout_array = stdout_array[:-1]
    stdout_array[-1] = stdout_array[-1].replace(",", "")
    return stdout_array


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
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
    )
    parser.add_argument(
        "--site",
        dest="site",
        default="root://maite.iihe.ac.be:1094",
        type=str,
        help="site to read the signal samples from",
    )
    args = parser.parse_args()

    # open dataset configs
    dataset_configs = get_dataset_config(args.year)

    # load already generated fileset
    nano_version = get_nano_version(args.year)
    filesets_dir = Path.cwd() / "analysis" / "filesets"
    fileset_file = filesets_dir / f"fileset_{args.year}_nanov{nano_version}_lxplus.json"
    built_dataset = {}
    if fileset_file.exists():
        with open(fileset_file, "r") as f:
            built_dataset = json.load(f)

    for sample in list(dataset_configs.keys()):
        if dataset_configs[sample]["era"] == "signal":
            query = dataset_configs[sample].get("query")
            path = dataset_configs[sample].get("path")
            if path:
                # private production on EOS: glob the local path(s) and expose them
                # through the EOS xrootd endpoint (readable from lxplus and condor)
                patterns = path if isinstance(path, list) else [path]
                files = sorted(f for pattern in patterns for f in glob.glob(pattern))
                if not files:
                    raise FileNotFoundError(f"no ROOT files found for {sample}: {patterns}")
                files = [
                    f"root://eosuser.cern.ch/{f}" if f.startswith("/eos/") else f
                    for f in files
                ]
                print(f"Adding sample {sample} ({len(files)} files from {patterns})")
                built_dataset[sample] = files
            elif len(query.split("/")) != 1:
                # from DAS queries
                files = get_signal_fileset(query)
                files = [f"{args.site}/{f}" for f in files]
                print(f"Adding sample {sample}")
                built_dataset[sample] = files
            else:
                # from already built signal samples
                signal_file = filesets_dir / f"{query}.json"
                print(f"Looking for file {signal_file}")
                if signal_file.exists():
                    with open(signal_file, "r") as json_file:
                        signal_json = json.load(json_file)
                    if sample not in signal_json:
                        print(f"sample key '{sample}' was not found. Available keys are {list(signal_json.keys())}")
                    else:
                        print(
                            f"Adding signal sample '{sample}'"
                        )
                        built_dataset.update({sample: signal_json[sample]})
                else:
                    print(f"File '{signal_file}_signal.json' was not found!")

    with open(fileset_file, "w") as json_file:
        json.dump(built_dataset, json_file, indent=4, sort_keys=True)