import pickle
from pathlib import Path


def generate_experiment_name(use_cell_frequencies=False, use_mon_ratio=False, use_mye_pathways=False,
                             use_TNFA_SIGNALING_VIA_NFKB_CD16_Mono=False, use_metadata_as_features=False,
                             keep_only_NOV=False, use_manualy_chosen_feats=False):
    features_names = []
    if use_cell_frequencies:
        features_names.append('comp')
    if use_mye_pathways:
        features_names.append('mye_pathways')
    elif use_TNFA_SIGNALING_VIA_NFKB_CD16_Mono:
        features_names += ["tnfa_via_nfkb_CD16_Mono"]
    if use_metadata_as_features:
        features_names.append('metadata')
    if use_mon_ratio:
        features_names.append('mon_ratio')
    if keep_only_NOV:
        features_names.append('NOV')

    if use_manualy_chosen_feats:
        features_names = ["manualy_chosen_feats"]
    experiment_name = '_'.join(features_names + ['results'])
    return experiment_name


def load_results_of_exp_name(exp_name: str, from_dir: Path):
    exp_results_path = Path(from_dir, f"{exp_name}.pkl")
    with open(exp_results_path, 'rb') as handle:
        results = pickle.load(handle)
    return results
