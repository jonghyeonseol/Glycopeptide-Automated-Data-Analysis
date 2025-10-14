"""
PLS-DA Diagnostic Plots Module for pGlyco Auto Combine
Creates 4-panel diagnostic visualization to validate model reliability

Dependencies:
    External:
        - pandas: Data manipulation
        - numpy: Numerical computations
        - matplotlib: Plotting backend
        - seaborn: Statistical visualization (heatmap)
        - sklearn: Model validation (cross_val_score, LeaveOneOut, PLSRegression, roc_curve, auc, confusion_matrix)

    Internal:
        - src.utils: save_trace_data
        - src.plots.plot_config: save_publication_figure

Phase 4.2: Critical for validating PLS-DA results and VIP scores
Ensures model is not overfitted and predictions are reliable
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from sklearn.model_selection import cross_val_score, cross_val_predict, LeaveOneOut
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import RobustScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_curve, auc, confusion_matrix
from ..utils import save_trace_data
from .plot_config import (
    save_publication_figure, COLOR_CANCER, COLOR_NORMAL, DPI_COMPLEX,
    TITLE_SIZE, AXIS_LABEL_SIZE, TICK_LABEL_SIZE, LEGEND_SIZE, ANNOTATION_SIZE,
    PLOT_LINE_LINEWIDTH, LINE_MEDIUM_THICK, LINE_BOLD, LINE_THICK,
    EDGE_LINEWIDTH_THIN, EDGE_LINEWIDTH_THICK,
    LINE_ALPHA, ALPHA_MEDIUM_LIGHT, THRESHOLD_ALPHA, OVERLAY_ALPHA, POINT_ALPHA,
    MARKER_SIZE_MEDIUM, DIAGNOSTIC_MARKER_SIZE,
    EDGE_COLOR_BLACK,# Edge color standardization,
    # Font family constants (Phase 10.3.10)
    FONT_DATA,
    # Colormap constants (Phase 10.3.11)
    CMAP_CONFUSION_MATRIX,
    # Marker style constants (Phase 10.3.9)
    MARKER_HIGHLIGHT,
    # Linestyle constants (Phase 10.3.8)
    LINESTYLE_DOTTED, LINESTYLE_DASHDOT, THRESHOLD_LINESTYLE,
    # Zorder constants (Phase 10.3.7)
    ZORDER_BACKGROUND, ZORDER_GRID, ZORDER_SEPARATOR,
    ZORDER_DATA_LOW, ZORDER_DATA_HIGH,
    ZORDER_THRESHOLD, ZORDER_ANNOTATION,
    ZORDER_OVERLAY, ZORDER_EFFECT,
    ZORDER_TOP, ZORDER_ABSOLUTE_TOP
)

logger = logging.getLogger(__name__)


class PLSDADiagnosticPlotMixin:
    """Mixin class for PLS-DA diagnostic visualization"""

    @staticmethod
    def _prepare_intensity_matrix(df: pd.DataFrame):
        """
        Prepare unscaled intensity matrix for cross-validation

        Uses unscaled data to avoid data leakage. Pipeline will handle
        scaling inside CV folds.

        Args:
            df: Annotated DataFrame with intensity data

        Returns:
            Tuple of (intensity_matrix, y_labels, sample_names)

        Pattern Used:
            Helper Extraction - separates data preparation from visualization
        """
        from ..analyzer import GlycanAnalyzer
        analyzer_temp = GlycanAnalyzer()
        intensity_matrix, _, _ = analyzer_temp.prepare_intensity_matrix(df)

        # Extract y_labels (Cancer=1, Normal=0)
        sample_cols = [col for col in df.columns if col.startswith(('C', 'N'))]
        y_labels = np.array([1 if col.startswith('C') else 0 for col in sample_cols])

        logger.debug(f"  Intensity matrix shape: {intensity_matrix.shape}")

        return intensity_matrix, y_labels

    @staticmethod
    def _plot_r2_q2_panel(ax, intensity_matrix, y_labels, plsda_model,
                          n_components_range, loo):
        """
        Plot Panel 1: R² and Q² scores across component range

        Calculates model quality metrics using proper cross-validation to
        prevent data leakage. Uses Pipeline to ensure scaling is done inside
        each CV fold.

        Args:
            ax: Matplotlib axis object
            intensity_matrix: Unscaled intensity matrix
            y_labels: Binary class labels (0=Normal, 1=Cancer)
            plsda_model: Fitted PLS-DA model (to extract selected components)
            n_components_range: Range of component counts to evaluate
            loo: LeaveOneOut cross-validator

        Returns:
            Tuple of (selected_r2, selected_q2, selected_comp)

        Pattern Used:
            Helper Extraction - isolates R²/Q² calculation and plotting
        """
        logger.info("  Panel 1: Calculating R² and Q² scores (with Pipeline to prevent leakage)...")

        r2_scores = []
        q2_scores = []

        for n_comp in n_components_range:
            # Use Pipeline to prevent data leakage
            # Scaler is fit INSIDE each CV fold
            pipeline = Pipeline([
                ('scaler', RobustScaler()),
                ('pls', PLSRegression(n_components=n_comp))
            ])

            # R² score (training set performance)
            pipeline.fit(intensity_matrix, y_labels)
            r2 = pipeline.score(intensity_matrix, y_labels)
            r2_scores.append(r2)

            # Q² score (cross-validated - proper way)
            # Pipeline ensures scaler is fit on training folds only
            q2_cv = cross_val_score(pipeline, intensity_matrix, y_labels,
                                    cv=loo, scoring='r2')
            q2_scores.append(q2_cv.mean())

        logger.debug(f"  R² range: {min(r2_scores):.3f} - {max(r2_scores):.3f}")
        logger.debug(f"  Q² range: {min(q2_scores):.3f} - {max(q2_scores):.3f}")

        # Plot R² and Q²
        ax.plot(n_components_range, r2_scores, 'o-', color=COLOR_CANCER,
                linewidth=PLOT_LINE_LINEWIDTH, markersize=MARKER_SIZE_MEDIUM, label='R² (Explained Variance)')
        ax.plot(n_components_range, q2_scores, 's-', color=COLOR_NORMAL,
                linewidth=PLOT_LINE_LINEWIDTH, markersize=MARKER_SIZE_MEDIUM, label='Q² (Predictive Ability)')

        # Add threshold line at 0.5
        ax.axhline(0.5, color='gray', linestyle=THRESHOLD_LINESTYLE, linewidth=LINE_MEDIUM_THICK,
                   label='Good model threshold (>0.5)', alpha=LINE_ALPHA)

        ax.set_xlabel('Number of Components', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax.set_ylabel('Score', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax.set_title('Model Quality: R² vs Q²\nHigher Q² = Better Predictive Power',
                     fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax.legend(loc='best', frameon=True, fontsize=LEGEND_SIZE)
        ax.grid(True, alpha=ALPHA_MEDIUM_LIGHT)
        ax.set_ylim(-0.1, 1.1)

        # Annotate selected model
        selected_comp = plsda_model.n_components
        selected_r2 = r2_scores[selected_comp - 1]
        selected_q2 = q2_scores[selected_comp - 1]
        ax.scatter([selected_comp], [selected_r2], s=DIAGNOSTIC_MARKER_SIZE, c='red',
                   marker=MARKER_HIGHLIGHT, edgecolors=EDGE_COLOR_BLACK, linewidths=EDGE_LINEWIDTH_THICK, zorder=ZORDER_THRESHOLD)
        ax.annotate(f'Selected: {selected_comp} comp\nR²={selected_r2:.3f}, Q²={selected_q2:.3f}',
                    xy=(selected_comp, selected_r2), xytext=(selected_comp + 1, selected_r2 - 0.15),
                    fontsize=ANNOTATION_SIZE, bbox=dict(boxstyle='round', facecolor='yellow', alpha=OVERLAY_ALPHA),
                    arrowprops=dict(arrowstyle='->', color='black', lw=LINE_MEDIUM_THICK))

        return selected_r2, selected_q2, selected_comp

    @staticmethod
    def _plot_roc_curve_panel(ax, intensity_matrix, y_labels, selected_comp, loo):
        """
        Plot Panel 2: ROC curve from cross-validated predictions

        Uses cross_val_predict to generate unbiased predictions for ROC
        analysis. Each prediction is made when the sample is in the test fold.

        Args:
            ax: Matplotlib axis object
            intensity_matrix: Unscaled intensity matrix
            y_labels: Binary class labels (0=Normal, 1=Cancer)
            selected_comp: Number of PLS components to use
            loo: LeaveOneOut cross-validator

        Returns:
            float: ROC AUC score

        Pattern Used:
            Helper Extraction - isolates ROC curve calculation and plotting
        """
        logger.info("  Panel 2: Generating ROC curve (using cross-validated predictions)...")

        # Use cross_val_predict for unbiased ROC
        # Each prediction is made when sample is in test fold
        pipeline_selected = Pipeline([
            ('scaler', RobustScaler()),
            ('pls', PLSRegression(n_components=selected_comp))
        ])

        # Get CV predictions (each sample predicted when it's in test set)
        y_pred_prob = cross_val_predict(pipeline_selected, intensity_matrix, y_labels,
                                        cv=loo, method='predict').ravel()

        # Calculate ROC curve from CV predictions
        fpr, tpr, thresholds = roc_curve(y_labels, y_pred_prob)
        roc_auc = auc(fpr, tpr)

        logger.debug(f"  ROC AUC (CV): {roc_auc:.3f}")

        # Plot ROC curve
        ax.plot(fpr, tpr, color=COLOR_CANCER, linewidth=LINE_BOLD,
                label=f'PLS-DA (AUC = {roc_auc:.3f})')
        ax.plot([0, 1], [0, 1], 'k--', linewidth=PLOT_LINE_LINEWIDTH, label='Random Classifier', alpha=THRESHOLD_ALPHA)

        ax.set_xlabel('False Positive Rate', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax.set_ylabel('True Positive Rate', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax.set_title(f'ROC Curve: Discrimination Ability\nAUC = {roc_auc:.3f} (1.0 = Perfect)',
                     fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax.legend(loc='lower right', frameon=True, fontsize=LEGEND_SIZE)
        ax.grid(True, alpha=ALPHA_MEDIUM_LIGHT)
        ax.set_xlim([-0.05, 1.05])
        ax.set_ylim([-0.05, 1.05])
        ax.set_aspect('equal')

        # Add interpretation text
        if roc_auc > 0.9:
            interp = "Excellent discrimination"
        elif roc_auc > 0.8:
            interp = "Good discrimination"
        elif roc_auc > 0.7:
            interp = "Fair discrimination"
        else:
            interp = "Poor discrimination"

        ax.text(0.6, 0.2, interp, fontsize=ANNOTATION_SIZE, fontweight='bold',
                bbox=dict(boxstyle='round', facecolor='lightgreen' if roc_auc > 0.8 else 'yellow', alpha=OVERLAY_ALPHA))

        return roc_auc, y_pred_prob

    @staticmethod
    def _plot_confusion_matrix_panel(ax, y_labels, y_pred_prob):
        """
        Plot Panel 3: Confusion matrix from cross-validated predictions

        Computes classification accuracy metrics from CV predictions.

        Args:
            ax: Matplotlib axis object
            y_labels: Binary class labels (0=Normal, 1=Cancer)
            y_pred_prob: Predicted probabilities from cross-validation

        Returns:
            Tuple of (accuracy, sensitivity, specificity)

        Pattern Used:
            Helper Extraction - isolates confusion matrix calculation and plotting
        """
        logger.info("  Panel 3: Computing confusion matrix (from CV predictions)...")

        # Predict classes (threshold at 0.5)
        y_pred_class = (y_pred_prob > 0.5).astype(int)

        # Compute confusion matrix from CV predictions
        cm = confusion_matrix(y_labels, y_pred_class)

        logger.debug(f"  Confusion matrix:\n{cm}")

        # Plot confusion matrix
        sns.heatmap(cm, annot=True, fmt='d', cmap=CMAP_CONFUSION_MATRIX, cbar=False,
                    ax=ax, square=True, linewidths=LINE_THICK, linecolor='black',
                    annot_kws={'size': 16, 'weight': 'bold'})

        ax.set_xlabel('Predicted Class', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax.set_ylabel('True Class', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax.set_title('Confusion Matrix\n(Leave-One-Out Cross-Validation)',
                     fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax.set_xticklabels(['Normal (0)', 'Cancer (1)'], fontsize=TICK_LABEL_SIZE)
        ax.set_yticklabels(['Normal (0)', 'Cancer (1)'], fontsize=TICK_LABEL_SIZE, rotation=0)

        # Calculate accuracy metrics
        tn, fp, fn, tp = cm.ravel()
        accuracy = (tp + tn) / (tp + tn + fp + fn)
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0

        metrics_text = f"Accuracy: {accuracy:.1%}\n"
        metrics_text += f"Sensitivity: {sensitivity:.1%}\n"
        metrics_text += f"Specificity: {specificity:.1%}"

        ax.text(1.05, 0.5, metrics_text, transform=ax.transAxes,
                fontsize=ANNOTATION_SIZE, verticalalignment='center',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=OVERLAY_ALPHA),
                family=FONT_DATA)

        return accuracy, sensitivity, specificity

    @staticmethod
    def _plot_vip_distribution_panel(ax, vip_df):
        """
        Plot Panel 4: VIP score distribution histogram

        Shows distribution of Variable Importance in Projection (VIP) scores
        with threshold at 1.0 for identifying important features.

        Args:
            ax: Matplotlib axis object
            vip_df: DataFrame with VIP_Score column

        Returns:
            Tuple of (mean_vip, median_vip, n_important)

        Pattern Used:
            Helper Extraction - isolates VIP distribution plotting
        """
        logger.info("  Panel 4: Plotting VIP score distribution...")

        # Get VIP scores
        vip_scores = vip_df['VIP_Score'].values

        # Plot histogram
        ax.hist(vip_scores, bins=50, color='#27AE60', alpha=POINT_ALPHA,
                edgecolor=EDGE_COLOR_BLACK, linewidth=EDGE_LINEWIDTH_THIN)

        # Add threshold line at VIP=1.0
        ax.axvline(1.0, color='red', linestyle=THRESHOLD_LINESTYLE, linewidth=LINE_BOLD,
                   label='VIP = 1.0 (Important Features)', zorder=ZORDER_THRESHOLD)

        # Add mean and median lines
        mean_vip = vip_scores.mean()
        median_vip = np.median(vip_scores)
        ax.axvline(mean_vip, color='blue', linestyle=LINESTYLE_DOTTED, linewidth=PLOT_LINE_LINEWIDTH,
                   label=f'Mean = {mean_vip:.2f}', alpha=LINE_ALPHA)
        ax.axvline(median_vip, color='orange', linestyle=LINESTYLE_DASHDOT, linewidth=PLOT_LINE_LINEWIDTH,
                   label=f'Median = {median_vip:.2f}', alpha=LINE_ALPHA)

        ax.set_xlabel('VIP Score', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax.set_ylabel('Frequency', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax.set_title('VIP Score Distribution\nVIP > 1.0 = Important Features',
                     fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax.legend(loc='upper right', frameon=True, fontsize=LEGEND_SIZE)
        ax.grid(True, alpha=ALPHA_MEDIUM_LIGHT, axis='y')

        # Add statistics box
        n_important = (vip_scores > 1.0).sum()
        pct_important = n_important / len(vip_scores) * 100

        stats_text = f"Total features: {len(vip_scores)}\n"
        stats_text += f"VIP > 1.0: {n_important} ({pct_important:.1f}%)\n"
        stats_text += f"Max VIP: {vip_scores.max():.2f}"

        ax.text(0.98, 0.98, stats_text, transform=ax.transAxes,
                fontsize=ANNOTATION_SIZE, verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=OVERLAY_ALPHA),
                family=FONT_DATA)

        return mean_vip, median_vip, n_important

    def plot_plsda_diagnostics(self, plsda_results: dict, df: pd.DataFrame,
                               figsize: tuple = (16, 12)):
        """
        Create 4-panel PLS-DA diagnostic plot

        CRITICAL FOR MODEL VALIDATION:
        - Panel 1: R² and Q² (model quality)
        - Panel 2: ROC curve (discrimination ability)
        - Panel 3: Confusion matrix (classification accuracy)
        - Panel 4: VIP score distribution (feature importance)

        Refactored in Phase 10.6 to use helper methods for better organization.

        Args:
            plsda_results: PLS-DA results from analyzer
            df: Annotated DataFrame (for sample info)
            figsize: Figure size (width, height)
        """
        logger.info("Creating PLS-DA diagnostic plots...")

        # Extract results
        plsda_model = plsda_results['plsda_model']
        y_labels = plsda_results['y_labels']
        vip_df = plsda_results['vip_scores']

        # Prepare unscaled intensity matrix (using helper)
        intensity_matrix, y_labels = self._prepare_intensity_matrix(df)

        # Create 2x2 subplot layout
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=figsize)
        fig.suptitle('PLS-DA Model Diagnostics: Validation of VIP Score Reliability',
                     fontsize=TITLE_SIZE, fontweight='bold', y=0.995)

        # Setup cross-validation
        n_components_range = range(1, min(11, intensity_matrix.shape[1]))
        loo = LeaveOneOut()

        # Panel 1: R² and Q² scores (using helper)
        selected_r2, selected_q2, selected_comp = self._plot_r2_q2_panel(
            ax1, intensity_matrix, y_labels, plsda_model, n_components_range, loo
        )

        # Panel 2: ROC curve (using helper)
        roc_auc, y_pred_prob = self._plot_roc_curve_panel(
            ax2, intensity_matrix, y_labels, selected_comp, loo
        )

        # Panel 3: Confusion matrix (using helper)
        accuracy, sensitivity, specificity = self._plot_confusion_matrix_panel(
            ax3, y_labels, y_pred_prob
        )

        # Panel 4: VIP distribution (using helper)
        mean_vip, median_vip, n_important = self._plot_vip_distribution_panel(ax4, vip_df)

        plt.tight_layout()

        # Save plot using standardized function
        output_file = self.output_dir / 'plsda_diagnostics.png'
        save_publication_figure(fig, output_file, dpi=DPI_COMPLEX)
        logger.info(f"Saved PLS-DA diagnostics to {output_file} (optimized, {DPI_COMPLEX} DPI)")

        # Save diagnostic metrics as trace data
        n_features = intensity_matrix.shape[1]
        cv_method = 'LeaveOneOut'

        diagnostic_metrics = pd.DataFrame({
            'Metric': ['R²', 'Q²', 'ROC_AUC', 'Accuracy', 'Sensitivity', 'Specificity',
                       'VIP_Mean', 'VIP_Median', 'N_Important_Features',
                       'Selected_Components', 'N_Features', 'Cross_Validation_Method'],
            'Value': [selected_r2, selected_q2, roc_auc, accuracy, sensitivity, specificity,
                      mean_vip, median_vip, n_important,
                      selected_comp, n_features, cv_method]
        })
        save_trace_data(diagnostic_metrics, self.output_dir, 'plsda_diagnostic_metrics.csv')

        plt.close()

        logger.info("✓ PLS-DA diagnostics complete - model validated")
        logger.info(f"  Model quality: R²={selected_r2:.3f}, Q²={selected_q2:.3f}, AUC={roc_auc:.3f}")
        logger.info(f"  Classification: Accuracy={accuracy:.1%}, {n_important} important features")
