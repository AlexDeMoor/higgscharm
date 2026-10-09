import awkward as ak


def add_cutflow(events, output, selection_manager, workflow_config, initial=None, skip_categories=()):
    """Add cutflow to output metadata

    `cutflow`: generator-weighted (genWeight only, no scale factors; scaled by
    lumi * xsec / sumw at postprocess); `cutflow_raw`: unweighted event counts.
    Rows: "initial" = the chunk BEFORE the jet veto map (the same events as the
    metadata sumw; `initial` = (sumw, nevents) of the pre-veto chunk, data: sumw =
    nevents), "jet_veto" = after the jet veto (the events the selections run on),
    then one row per cut of the category. Categories in `skip_categories` (the
    mc_only_categories on data) get no cutflow."""
    is_mc = hasattr(events, "genWeight")
    sumw = ak.sum(events.genWeight) if is_mc else len(events)
    if initial is None:
        initial = (sumw, len(events))
    for category, category_cuts in workflow_config.event_selection[
        "categories"
    ].items():
        if category in skip_categories:
            continue
        output["metadata"].update(
            {
                category: {
                    "cutflow": {"initial": initial[0], "jet_veto": sumw},
                    "cutflow_raw": {"initial": initial[1], "jet_veto": len(events)},
                }
            }
        )
        selections = []
        for cut_name in category_cuts:
            selections.append(cut_name)
            current_selection = selection_manager.all(*selections)
            n_raw = int(ak.sum(current_selection))
            if n_raw != 0:
                sumw_cutflow = (
                    ak.sum(events.genWeight[current_selection])
                    if is_mc
                    else n_raw
                )
                output["metadata"][category]["cutflow"][cut_name] = sumw_cutflow
            else:
                output["metadata"][category]["cutflow"][cut_name] = 0
            output["metadata"][category]["cutflow_raw"][cut_name] = n_raw


def update(events, collections):
    """Return a shallow copy of events array with some collections swapped out"""
    out = events
    for name, value in collections.items():
        out = ak.with_field(out, value, name)
    return out
