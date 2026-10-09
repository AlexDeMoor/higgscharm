import warnings
import numpy as np


def add_scalevar_weight(events, weights_container, shift):
    """
    Twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopSystematics#Factorization_and_renormalizatio

    __doc__:
    ['LHE scale variation weights (w_var / w_nominal)',
    ' [0] is renscfact=0.5d0 facscfact=0.5d0 ',
    ' [1] is renscfact=0.5d0 facscfact=1d0 ',
    ' [2] is renscfact=0.5d0 facscfact=2d0 ',
    ' [3] is renscfact=1d0 facscfact=0.5d0 ',
    ' [4] is renscfact=1d0 facscfact=1d0 ',
    ' [5] is renscfact=1d0 facscfact=2d0 ',
    ' [6] is renscfact=2d0 facscfact=0.5d0 ',
    ' [7] is renscfact=2d0 facscfact=1d0 ',
    ' [8] is renscfact=2d0 facscfact=2d0 ']
    """
    lhe_weights = events.LHEScaleWeight
    nom = np.ones(len(weights_container.weight()))
    if shift is None:
        if len(lhe_weights) > 0:
            if len(lhe_weights[0]) == 9:
                # the nominal weight is NOT modified (each `add` multiplies it): [4]
                # (muR = muF = 1) is 1 in most samples but not in the 4FS t-channel
                # single top (up to ~7 or negative), where using it as nominal x3
                # distorted the yields. Up/down = variation / [4], relative to the
                # central scale; coffea divides them by the (unit) nominal passed here.
                central = lhe_weights[:, 4]
                weights_container.add(
                    "scalevar_muR",
                    nom,
                    lhe_weights[:, 1] / central,
                    lhe_weights[:, 7] / central,
                )
                weights_container.add(
                    "scalevar_muF",
                    nom,
                    lhe_weights[:, 3] / central,
                    lhe_weights[:, 5] / central,
                )
                weights_container.add(
                    "scalevar_muR_muF",
                    nom,
                    lhe_weights[:, 0] / central,
                    lhe_weights[:, 8] / central,
                )
            elif len(lhe_weights[0]) > 1:
                print("Scale variation vector has length ", len(lhe_weights[0]))
        else:
            warnings.warn(
                "LHE scale variation weights are not available, put nominal weights"
            )
            weights_container.add("scalevar_muR", nom, nom, nom)
            weights_container.add("scalevar_muF", nom, nom, nom)
            weights_container.add("scalevar_muR_muF", nom, nom, nom)

    else:
        weights_container.add("scalevar_3pt", nom)
