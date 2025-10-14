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

    def plot_plsda_diagnostics(self, plsda_results: dict, df: pd.DataFrame,
                               figsize: tuple = (16, 12)):
        """
        Create 4-panel PLS-DA diagnostic plot

        CRITICAL FOR MODEL VALIDATION:
        - Panel 1: R² and Q² (model quality)
        - Panel 2: ROC curve (discrimination ability)
        - Panel 3: Confusion matrix (classification accuracy)
        - Panel 4: VIP score distribution (feature importance)

        Args:
            plsda_results: PLS-DA results from analyzer
            df: Annotated DataFrame (for sample info)
            figsize: Figure size (width, height)
        """
        logger.info("Creating PLS-DA diagnostic plots...")

        # Extract results
        plsda_model = plsda_results['plsda_model']
        plsda_results['X_scores']
        y_labels = plsda_results['y_labels']
        plsda_results['sample_names']
        vip_df = plsda_results['vip_scores']

        # ========================================
        # FIX: Use unscaled data to avoid leakage
        # Pipeline will handle scaling inside CV
        # ========================================
        from ..analyzer import GlycanAnalyzer
        analyzer_temp = GlycanAnalyzer()
        intensity_matrix, _, _ = analyzer_temp.prepare_intensity_matrix(df)

        logger.debug(f"  Intensity matrix shape: {intensity_matrix.shape}")

        # Create 2x2 subplot layout
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=figsize)
        fig.suptitle('PLS-DA Model Diagnostics: Validation of VIP Score Reliability',
                     fontsize=TITLE_SIZE, fontweight='bold', y=0.995)

        # =====================================================================
        # PANEL 1: R² and Q² across components (NO LEAKAGE)
        # =====================================================================
        logger.info("  Panel 1: Calculating R² and Q² scores (with Pipeline to prevent leakage)...")

        n_components_range = range(1, min(11, intensity_matrix.shape[1]))
        r2_scores = []
        q2_scores = []
        loo = LeaveOneOut()

        for n_comp in n_components_range:
            # ========================================
            # FIX: Use Pipeline to prevent data leakage
            # Scaler is fit INSIDE each CV fold
            # ========================================
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
        ax1.plot(n_components_range, r2_scores, 'o-', color=COLOR_CANCER,
                 linewidth=PLOT_LINE_LINEWIDTH, markersize=MARKER_SIZE_MEDIUM, label='R² (Explained Variance)')
        ax1.plot(n_components_range, q2_scores, 's-', color=COLOR_NORMAL,
                 linewidth=PLOT_LINE_LINEWIDTH, markersize=MARKER_SIZE_MEDIUM, label='Q² (Predictive Ability)')

        # Add threshold line at 0.5
        ax1.axhline(0.5, color='gray', linestyle=THRESHOLD_LINESTYLE, linewidth=LINE_MEDIUM_THICK,
                    label='Good model threshold (>0.5)', alpha=LINE_ALPHA)

        ax1.set_xlabel('Number of Components', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax1.set_ylabel('Score', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax1.set_title('Model Quality: R² vs Q²\nHigher Q² = Better Predictive Power',
                      fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax1.legend(loc='best', frameon=True, fontsize=LEGEND_SIZE)
        ax1.grid(True, alpha=ALPHA_MEDIUM_LIGHT)
        ax1.set_ylim(-0.1, 1.1)

        # Annotate selected model
        selected_comp = plsda_model.n_components
        selected_r2 = r2_scores[selected_comp - 1]
        selected_q2 = q2_scores[selected_comp - 1]
        ax1.scatter([selected_comp], [selected_r2], s=DIAGNOSTIC_MARKER_SIZE, c='red',
                    marker='*', edgecolors=EDGE_COLOR_BLACK, linewidths=EDGE_LINEWIDTH_THICK, zorder=ZORDER_THRESHOLD)
        ax1.annotate(f'Selected: {selected_comp} comp\nR²={selected_r2:.3f}, Q²={selected_q2:.3f}',
                     xy =(selected_comp, selected_r2), xytext =(selected_comp + 1, selected_r2 - 0.15),
                     fontsize=ANNOTATION_SIZE, bbox=dict(boxstyle='round', facecolor='yellow', alpha=OVERLAY_ALPHA),
                     arrowprops=dict(arrowstyle='->', color='black', lw=LINE_MEDIUM_THICK))

        # =====================================================================
        # PANEL 2: ROC Curve (NO LEAKAGE)
        # =====================================================================
        logger.info("  Panel 2: Generating ROC curve (using cross-validated predictions)...")

        # ========================================
        # FIX: Use cross_val_predict for unbiased ROC
        # Each prediction is made when sample is in test fold
        # ========================================
        selected_comp = plsda_model.n_components
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
        ax2.plot(fpr, tpr, color=COLOR_CANCER, linewidth=LINE_BOLD,
                 label=f'PLS-DA (AUC = {roc_auc:.3f})')
        ax2.plot([0, 1], [0, 1], 'k--', linewidth=PLOT_LINE_LINEWIDTH, label='Random Classifier', alpha=THRESHOLD_ALPHA)

        ax2.set_xlabel('False Positive Rate', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax2.set_ylabel('True Positive Rate', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax2.set_title(f'ROC Curve: Discrimination Ability\nAUC = {roc_auc:.3f} (1.0 = Perfect)',
                      fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax2.legend(loc='lower right', frameon=True, fontsize=LEGEND_SIZE)
        ax2.grid(True, alpha=ALPHA_MEDIUM_LIGHT)
        ax2.set_xlim([-0.05, 1.05])
        ax2.set_ylim([-0.05, 1.05])
        ax2.set_aspect('equal')

        # Add interpretation text
        if roc_auc > 0.9:
            interp = "Excellent discrimination"
        elif roc_auc > 0.8:
            interp = "Good discrimination"
        elif roc_auc > 0.7:
            interp = "Fair discrimination"
        else:
            interp = "Poor discrimination"

        ax2.text(0.6, 0.2, interp, fontsize=ANNOTATION_SIZE, fontweight='bold',
                 bbox=dict(boxstyle='round', facecolor='lightgreen' if roc_auc > 0.8 else 'yellow', alpha=OVERLAY_ALPHA))

        # =====================================================================
        # PANEL 3: Confusion Matrix (NO LEAKAGE)
        # =====================================================================
        logger.info("  Panel 3: Computing confusion matrix (from CV predictions)...")

        # ========================================
        # FIX: Use CV predictions (already computed above)
        # These are unbiased predictions from LOO-CV
        # ========================================
        # Predict classes (threshold at 0.5)
        y_pred_class = (y_pred_prob > 0.5).astype(int)

        # Compute confusion matrix from CV predictions
        cm = confusion_matrix(y_labels, y_pred_class)

        logger.debug(f"  Confusion matrix:\n{cm}")

        # Plot confusion matrix
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', cbar=False,
                    ax=ax3, square=True, linewidths=LINE_THICK, linecolor='black',
                    annot_kws={'size': 16, 'weight': 'bold'})

        ax3.set_xlabel('Predicted Class', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax3.set_ylabel('True Class', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax3.set_title('Confusion Matrix\n(Leave-One-Out Cross-Validation)',
                      fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax3.set_xticklabels(['Normal (0)', 'Cancer (1)'], fontsize=TICK_LABEL_SIZE)
        ax3.set_yticklabels(['Normal (0)', 'Cancer (1)'], fontsize=TICK_LABEL_SIZE, rotation=0)

        # Calculate accuracy metrics
        tn, fp, fn, tp = cm.ravel()
        accuracy = (tp + tn) / (tp + tn + fp + fn)
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0

        metrics_text = f"Accuracy: {accuracy:.1%}\n"
        metrics_text += f"Sensitivity: {sensitivity:.1%}\n"
        metrics_text += f"Specificity: {specificity:.1%}"

        ax3.text(1.05, 0.5, metrics_text, transform=ax3.transAxes,
                 fontsize=ANNOTATION_SIZE, verticalalignment='center',
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=OVERLAY_ALPHA),
                 family='monospace')

        # =====================================================================
        # PANEL 4: VIP Score Distribution
        # =====================================================================
        logger.info("  Panel 4: Plotting VIP score distribution...")

        # Get VIP scores
        vip_scores = vip_df['VIP_Score'].values

        # Plot histogram
        ax4.hist(vip_scores, bins=50, color='#27AE60', alpha=POINT_ALPHA,
                 edgecolor=EDGE_COLOR_BLACK, linewidth=EDGE_LINEWIDTH_THIN)

        # Add threshold line at VIP=1.0
        ax4.axvline(1.0, color='red', linestyle=THRESHOLD_LINESTYLE, linewidth=LINE_BOLD,
                    label='VIP = 1.0 (Important Features)', zorder=ZORDER_THRESHOLD)

        # Add mean and median lines
        mean_vip = vip_scores.mean()
        median_vip = np.median(vip_scores)
        ax4.axvline(mean_vip, color='blue', linestyle=LINESTYLE_DOTTED, linewidth=PLOT_LINE_LINEWIDTH,
                    label=f'Mean = {mean_vip:.2f}', alpha=LINE_ALPHA)
        ax4.axvline(median_vip, color='orange', linestyle=LINESTYLE_DASHDOT, linewidth=PLOT_LINE_LINEWIDTH,
                    label=f'Median = {median_vip:.2f}', alpha=LINE_ALPHA)

        ax4.set_xlabel('VIP Score', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax4.set_ylabel('Frequency', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax4.set_title('VIP Score Distribution\nVIP > 1.0 = Important Features',
                      fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax4.legend(loc='upper right', frameon=True, fontsize=LEGEND_SIZE)
        ax4.grid(True, alpha=ALPHA_MEDIUM_LIGHT, axis='y')

        # Add statistics box
        n_important = (vip_scores > 1.0).sum()
        pct_important = n_important / len(vip_scores) * 100

        stats_text = f"Total features: {len(vip_scores)}\n"
        stats_text += f"VIP > 1.0: {n_important} ({pct_important:.1f}%)\n"
        stats_text += f"Max VIP: {vip_scores.max():.2f}"

        ax4.text(0.98, 0.98, stats_text, transform=ax4.transAxes,
                 fontsize=ANNOTATION_SIZE, verticalalignment='top', horizontalalignment='right',
                 bbox=dict(boxstyle='round', facecolor='lightblue', alpha=OVERLAY_ALPHA),
                 family='monospace')

        plt.tight_layout()

        # Save plot using standardized function
        output_file = self.output_dir / 'plsda_diagnostics.png'
        save_publication_figure(fig, output_file, dpi=DPI_COMPLEX)
        logger.info(f"Saved PLS-DA diagnostics to {output_file} (optimized, {DPI_COMPLEX} DPI)")

        # Save diagnostic metrics as trace data
        # Enhanced with pipeline metadata (Phase 10.7)
        n_features = intensity_matrix.shape[1]  # Total number of features
        cv_method = 'LeaveOneOut'  # Cross-validation method used

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
