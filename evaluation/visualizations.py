import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap
from sklearn.metrics import ConfusionMatrixDisplay, confusion_matrix
from sklearn.metrics import RocCurveDisplay


def plot_ROC_PRauc_CM_stem(y_true, y_pred, y_proba, pos_label=1, flip_stem=False, save_figure_to_path=None,
                           use_all_score_range=False):
    fig_height = 8 if flip_stem else 6
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24, fig_height))

    RocCurveDisplay.from_predictions(y_true=y_true, y_pred=y_proba, pos_label=pos_label,
                                     plot_chance_level=True, ax=ax1)

    cm = confusion_matrix(y_true, y_pred)
    cm_display = ConfusionMatrixDisplay(cm)

    cm_display.plot(ax=ax2)
    plot_stemplot(y_true, y_proba, ax=ax3, rotate=flip_stem, use_all_score_range=use_all_score_range)

    if save_figure_to_path is not None:
        plt.savefig(save_figure_to_path, format="svg")
    plt.show()


def print_shap_plots(model, X):
    try:  # tree
        explainer = shap.TreeExplainer(model)
    except:
        try:  # kernel
            explainer = shap.Explainer(model, X)
        except:
            explainer = shap.KernelExplainer(model.predict, X)
    shap_values = explainer(X)
    if len(shap_values.shape) > 2:
        shap_values = shap_values[:, :, 1]
    # shap.plots.beeswarm(shap_values)
    shap.plots.bar(shap_values, max_display=5)
    return shap_values


def plot_stemplot(y_true, y_proba, ax=None, rotate=False, plot_sample_ind=True, use_all_score_range=False,
                  classifier_threshold=0.5):
    y_df = pd.DataFrame(y_true)
    y_df["score"] = y_proba
    y_df = y_df.sort_values(by="score", ascending=rotate)
    if use_all_score_range:
        thresh = classifier_threshold
        max_score = max(y_df["score"])
        min_score = min(y_df["score"])

        def _above_thresh_extrapolation(x, threshold):
            return threshold + threshold * ((x - threshold) / (max_score - threshold))

        def _bellow_thresh_extrapolation(x, threshold):
            return threshold - threshold * ((x - threshold) / (min_score - threshold))

        y_df["score"] = y_df["score"].apply(
            lambda x: _bellow_thresh_extrapolation(x, thresh) if x < thresh else _above_thresh_extrapolation(x, thresh))
    y_df = y_df.reset_index().reset_index()

    y_df_response = y_df[y_df[y_true.name] == 1]
    y_df_non_response = y_df[y_df[y_true.name] == 0]

    R_NR_colors = {
        "R": "#71cceb",  # "blue",
        "NR": "#b30000"  # "red"
    }
    stem_orientation = 'horizontal' if rotate else 'vertical'
    stem_kwrgs = {"bottom": classifier_threshold,
                  "orientation": stem_orientation,
                  "basefmt": "k-"}
    if ax is None:
        (markers1, stemlines1, baseline1) = plt.stem(y_df_response['index'], y_df_response["score"],
                                                     **stem_kwrgs)
        (markers2, stemlines2, baseline2) = plt.stem(y_df_non_response['index'], y_df_non_response["score"],
                                                     **stem_kwrgs)
    else:
        (markers1, stemlines1, baseline1) = ax.stem(y_df_response['index'], y_df_response["score"], **stem_kwrgs)
        (markers2, stemlines2, baseline2) = ax.stem(y_df_non_response['index'], y_df_non_response["score"],
                                                    **stem_kwrgs)
    plt.setp(markers1, marker='o', markersize=3, markeredgecolor=R_NR_colors["R"], markeredgewidth=2)
    plt.setp(markers2, marker='o', markersize=3, markeredgecolor=R_NR_colors["NR"], markeredgewidth=2)
    plt.setp(stemlines1, linestyle="-", color=R_NR_colors["R"], linewidth=2)
    plt.setp(stemlines2, linestyle="-", color=R_NR_colors["NR"], linewidth=2)

    if not rotate:
        plt.ylabel('Prediction scores')
    else:
        plt.xlabel('Prediction scores')

    buffer = 0.1
    if ax is not None and plot_sample_ind:
        if not rotate:
            ax.set_xticks(np.arange(len(y_df['patient'])))
            ax.set_xticklabels(y_df['patient'])
            plt.xticks(rotation=90)
            plt.ylim([0 - buffer, 1 + buffer])
        else:
            ax.set_yticks(np.arange(len(y_df['patient'])))
            ax.set_yticklabels(y_df['patient'])
            plt.xlim([0 - buffer, 1 + buffer])


def plot_loocv_roc_curves(all_experiments_dict, exp_names_to_plot, fig_path=None):
    fig, ax = plt.subplots(figsize=(6, 6))
    for i, (exp_name, new_exp_name) in enumerate(exp_names_to_plot.items()):
        results = all_experiments_dict[exp_name]
        exp_name = new_exp_name
        y_true = results['loocv']["y"]
        y_proba = np.array(results['loocv']["y_proba"])
        viz = RocCurveDisplay.from_predictions(
            y_true,
            y_proba,
            name=f"{exp_name}",
            alpha=1,
            lw=1,
            ax=ax,
            plot_chance_level=(i == len(exp_names_to_plot) - 1),
        )

        xlabel = "False Positive Rate"
        ylabel = "True Positive Rate"
        ax.set(xlabel=xlabel, ylabel=ylabel)
    if fig_path is not None:
        fig.savefig(fig_path, bbox_inches='tight', format="svg")
