import re
import glob
import yaml
import json
import subprocess
import numpy as np
from pathlib import Path
from collections import defaultdict
from analysis.workflows.config import WorkflowConfigBuilder


def get_rootfiles(year: str, dataset: str):
    main_dir = Path.cwd()
    fileset_path = Path(f"{main_dir}/analysis/filesets")
    with open(f"{fileset_path}/{year}_fileset.yaml", "r") as f:
        dataset_config = yaml.safe_load(f)[dataset]

    # check for .root files in the specified path
    root_files = glob.glob(f"{dataset_config['path']}/*.root")
    if not root_files:
        # if no files found, check in the subdirectories
        root_files = glob.glob(f"{dataset_config['path']}/*/*.root")
    if not root_files:
        raise FileNotFoundError(
            f"No .root files found in {dataset_config['path']} or its subdirectories."
        )
    return root_files


def divide_list(lst: list, nfiles: int = 20) -> list:
    """Divide a list into sublists such that each sublist has at least 20 elements."""
    if len(lst) < nfiles:
        return [lst]

    # Dynamically calculate the number of sublists such that each has at least 20 elements
    n = len(lst) // nfiles  # This gives the number of groups with at least 20 elements
    if len(lst) % nfiles != 0:
        n += 1  # Increase n by 1 if there is a remainder, to accommodate extra elements

    # Divide the list into 'n' sublists
    size = len(lst) // n
    remainder = len(lst) % n
    result = []
    start = 0

    for i in range(n):
        if i < remainder:
            end = start + size + 1
        else:
            end = start + size
        result.append(lst[start:end])
        start = end
    return result


def get_nano_version(year: str):
    if year.startswith("201"):
        nano_version = "9"
    elif year == "2024":
        nano_version = "15"
    else:
        nano_version = "12"
    return nano_version


def get_dataset_config(year: str):
    aux_year_map = {
        "2016": "2016preVFP",
        "2022": "2022preEE",
        "2023": "2023preBPix",
    }
    aux_year = aux_year_map.get(year, year)
    fileset_path = Path.cwd() / "analysis" / "filesets"
    nano_version = get_nano_version(year)
    fileset_file = f"{fileset_path}/{aux_year}_nanov{nano_version}.yaml"
    with open(fileset_file, "r") as f:
        dataset_config = yaml.safe_load(f)
    return dataset_config


def get_dataset_name(dataset):
    datasets = ["MuonEG", "Muon", "EGamma", "SingleMuon", "DoubleMuon", "Tau"]
    for dataset_key in datasets:
        if dataset.startswith(dataset_key):
            return dataset_key
    return "mc"


def get_dataset_era(dataset: str, year: str):
    dataset_config = get_dataset_config(year)
    for dataset_key in dataset_config:
        if dataset.startswith(dataset_key):
            return dataset_config[dataset_key]["era"]


def modify_site_list(year: str, site: str, status: str) -> None:
    """
    Move a given site to the specified status list ("white" or "black")

    Parameters:
    -----------
        site: The site identifier to modify.
        status: Desired list to move the site to ("white" or "black").
    """
    yaml_file = Path.cwd() / "analysis" / "filesets" / f"{year}_sites.yaml"
    with open(yaml_file, "r") as f:
        data = yaml.safe_load(f)

    if status == "white":
        if site not in data["white"]:
            data["white"].append(site)
        if site in data["black"]:
            data["black"].remove(site)
    else:
        if site not in data["black"]:
            data["black"].append(site)
        if site in data["white"]:
            data["white"].remove(site)

    with open(yaml_file, "w") as f:
        yaml.dump(data, f, default_flow_style=False)


def extract_xrootd_errors(error_files: list) -> set:
    """
    Extract xrootd URLs from a list of Condor error log files.

    Parameters:
    -----------
        error_files: (list of str or Path): Paths to error log files.

    Returns:
    --------
        set: A set of unique xrootd endpoints found in the logs.
    """
    xrootd_errs = set()
    for error_file in error_files:
        with open(error_file, "r") as f:
            err_content = f.read()
        xrootd_matches = re.findall(r"root://[a-zA-Z0-9\-.]+(?:[:]\d+)?", err_content)
        xrootd_errs.update(xrootd_matches)
    return xrootd_errs


def fileset_checker(samples: list, year: str):
    """check if the fileset for the given year exists, generate it otherwise"""
    nano_version = get_nano_version(year)
    filesets_path = Path.cwd() / "analysis" / "filesets"
    fileset_file = filesets_path / f"fileset_{year}_nanov{nano_version}_lxplus.json"

    build_input_fileset = False
    if not fileset_file.exists():
        build_input_fileset = True
    else:
        with open(fileset_file, "r") as f:
            filesets = json.load(f)
        build_input_fileset = any(ds not in filesets for ds in samples)

    if build_input_fileset:
        print("\nBuilding input filesets for:")
        print(yaml.dump(samples, default_flow_style=False, sort_keys=False, indent=2))
        cmd = f"python3 fetch.py --year {year} --samples {' '.join(samples)}"
        subprocess.run(cmd, shell=True)


def get_datasets_map(year: str):
    dataset_config = get_dataset_config(year)
    datasets_map = {}
    for sample in dataset_config:
        keys = dataset_config[sample]["key"]
        if isinstance(keys, str):
            if keys in datasets_map:
                datasets_map[keys] += [sample]
            else:
                datasets_map[keys] = [sample]
        elif isinstance(keys, list):
            for key in keys:
                if key in datasets_map:
                    datasets_map[key] += [sample]
                else:
                    datasets_map[key] = [sample]
    return datasets_map


# `datasets:` groups of a workflow yaml that are processed by default, and groups that are only
# declared (validated, never run): `data_blinded` keeps the primary datasets of a blinded analysis
# out of every default run and keeps run_postprocess.py / the plotter blinded (they test `data`).
DATASET_GROUPS_TO_RUN = ("data", "mc", "signal")
DATASET_GROUPS_DECLARED_ONLY = ("data_blinded",)
ALL_YEARS = (
    "2016preVFP",
    "2016postVFP",
    "2017",
    "2018",
    "2022preEE",
    "2022postEE",
    "2023preBPix",
    "2023postBPix",
    "2024",
)


def get_datasets_to_run_over(workflow: str, year: str):
    """fileset samples of the `data` / `mc` / `signal` keys of a workflow for one year

    Fails loudly on a key that no fileset yaml of any year defines (a typo would otherwise
    silently drop a sample) and on an unknown `datasets:` group. Keys that exist for other
    years only are skipped with one summary warning: several Run 2 / Run 3 workflows
    (zzto4l, hplusc, zplusl*, zplusll*, ctag_eff) list the union of the per-year keys.
    """
    config_builder = WorkflowConfigBuilder(workflow=workflow)
    workflow_config = config_builder.build_workflow_config()
    datasets = workflow_config.datasets

    unknown_groups = [
        g for g in datasets
        if g not in DATASET_GROUPS_TO_RUN + DATASET_GROUPS_DECLARED_ONLY
    ]
    if unknown_groups:
        raise ValueError(
            f"workflow '{workflow}': unknown `datasets:` group(s) {unknown_groups}; "
            f"expected {list(DATASET_GROUPS_TO_RUN + DATASET_GROUPS_DECLARED_ONLY)}"
        )

    samples_keys_to_run_over = []
    for group in DATASET_GROUPS_TO_RUN:
        samples_keys_to_run_over += list(datasets.get(group) or [])
    declared_keys = []
    for group in DATASET_GROUPS_DECLARED_ONLY:
        declared_keys += list(datasets.get(group) or [])

    datasets_map = get_datasets_map(year)
    all_keys = set(datasets_map)
    for other_year in ALL_YEARS:
        if other_year != year:
            all_keys |= set(get_datasets_map(other_year))
    unknown = [k for k in samples_keys_to_run_over + declared_keys if k not in all_keys]
    if unknown:
        raise ValueError(
            f"workflow '{workflow}': dataset key(s) {unknown} match no `key:` of any fileset "
            f"yaml in analysis/filesets/ (typo?)"
        )

    datasets_to_run_over = []
    missing_this_year = []
    for key in samples_keys_to_run_over:
        if key not in datasets_map:
            missing_this_year.append(key)
            continue
        datasets_to_run_over += datasets_map[key]
    if missing_this_year:
        print(
            f"\nWARNING: workflow '{workflow}': key(s) {missing_this_year} have no {year} "
            f"sample (defined for other years only): skipped"
        )
    return datasets_to_run_over


def get_workflow_key_process_map(workflow_config, year: str):
    datasets = workflow_config.datasets
    dataset_configs = get_dataset_config(year)
    sample_keys = np.concatenate(list(datasets.values())).tolist()

    processes = []
    key_process_map = {}

    for sample in dataset_configs:
        sample_process = dataset_configs[sample]["process"]
        if sample_process == "Data":
            continue
        sample_key = dataset_configs[sample]["key"]
        if isinstance(sample_key, str):
            if sample_key in sample_keys:
                if sample_process not in processes:
                    processes.append(sample_process)
                    key_process_map[sample_key] = sample_process
        elif isinstance(sample_key, list):
            for sk in sample_key:
                if sk in sample_keys:
                    if sample_process not in processes:
                        processes.append(sample_process)
                        key_process_map[sk] = sample_process

    return key_process_map


def get_workflow_processes(workflow_config, year: str):
    """every (non-Data) process of the workflow's dataset keys, in fileset order. Unlike
    get_workflow_key_process_map (one process per key), a key shared by several processes
    (e.g. 2024 `dy_nlo` = DYll + DYtautau) contributes all of them"""
    sample_keys = set(np.concatenate(list(workflow_config.datasets.values())).tolist())
    processes = []
    for sample, config in get_dataset_config(year).items():
        if config["process"] == "Data":
            continue
        keys = config["key"] if isinstance(config["key"], list) else [config["key"]]
        if sample_keys.intersection(keys) and config["process"] not in processes:
            processes.append(config["process"])
    return processes


def get_process_era_map(year: str):
    dataset_configs = get_dataset_config(year)
    process_era_map = {}
    for sample in dataset_configs:
        sample_process = dataset_configs[sample]["process"]
        if sample_process not in process_era_map:
            process_era_map[sample_process] = dataset_configs[sample]["era"]

    return process_era_map


def get_process_sample_map(datasets: list[str], year: str) -> dict[str, list[str]]:
    """map processes to their corresponding samples based on dataset config"""
    dataset_configs = get_dataset_config(year)
    process_sample_map = defaultdict(list)
    for sample in datasets:
        config = dataset_configs[sample]
        process_sample_map[config["process"]].append(sample)
    return dict(process_sample_map)
