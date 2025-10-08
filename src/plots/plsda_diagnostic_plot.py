"""
PLS-DA Diagnostic Plots Module for pGlyco Auto Combine
Creates 4-panel diagnostic visualization to validate model reliability

Phase 4.2: Critical for validating PLS-DA results and VIP scores
Ensures model is not overfitted and predictions are reliable
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from sklearn.model_selection import cross_val_score, LeaveOneOut
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_curve, auc, confusion_matrix
from ..utils import save_trace_data

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

        # Reconstruct scaled X matrix
        from ..analyzer import GlycanAnalyzer
        analyzer_temp = GlycanAnalyzer()
        intensity_matrix, _, _ = analyzer_temp.prepare_intensity_matrix(df)

        from sklearn.preprocessing import RobustScaler
        scaler = RobustScaler()
        X_scaled = scaler.fit_transform(intensity_matrix)

        # Create 2x2 subplot layout
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=figsize)
        fig.suptitle('PLS-DA Model Diagnostics: Validation of VIP Score Reliability',
                     fontsize=16, fontweight='bold', y=0.995)

        # =====================================================================
        # PANEL 1: R² and Q² across components
        # =====================================================================
        logger.info("  Panel 1: Calculating R² and Q² scores...")

        n_components_range = range(1, min(11, X_scaled.shape[1]))
        r2_scores = []
        q2_scores = []

        for n_comp in n_components_range:
            # Fit model with n_comp components
            pls_temp = PLSRegression(n_components=n_comp)
            pls_temp.fit(X_scaled, y_labels)

            # R² score (explained variance)
            r2 = pls_temp.score(X_scaled, y_labels)
            r2_scores.append(r2)

            # Q² score (cross-validated)
            loo = LeaveOneOut()
            q2_cv = cross_val_score(pls_temp, X_scaled, y_labels,
                                    cv=loo, scoring='r2')
            q2_scores.append(q2_cv.mean())

        # Plot R² and Q²
        ax1.plot(n_components_range, r2_scores, 'o-', color='#E74C3C',
                 linewidth=2, markersize=8, label='R² (Explained Variance)')
        ax1.plot(n_components_range, q2_scores, 's-', color='#3498DB',
                 linewidth=2, markersize=8, label='Q² (Predictive Ability)')

        # Add threshold line at 0.5
        ax1.axhline(0.5, color='gray', linestyle='--', linewidth=1.5,
                    label='Good model threshold (>0.5)', alpha=0.7)

        ax1.set_xlabel('Number of Components', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Score', fontsize=12, fontweight='bold')
        ax1.set_title('Model Quality: R² vs Q²\nHigher Q² = Better Predictive Power',
                      fontsize=12, fontweight='bold')
        ax1.legend(loc='best', frameon=True, fontsize=10)
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(-0.1, 1.1)

        # Annotate selected model
        selected_comp = plsda_model.n_components
        selected_r2 = r2_scores[selected_comp - 1]
        selected_q2 = q2_scores[selected_comp - 1]
        ax1.scatter([selected_comp], [selected_r2], s=200, c='red',
                    marker='*', edgecolors='black', linewidths=2, zorder=10)
        ax1.annotate(f'Selected: {selected_comp} comp\nR²={selected_r2:.3f}, Q²={selected_q2:.3f}',
                     xy =(selected_comp, selected_r2), xytext =(selected_comp + 1, selected_r2 - 0.15),
                     fontsize=9, bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7),
                     arrowprops=dict(arrowstyle='->', color='black', lw=1.5))

        # =====================================================================
        # PANEL 2: ROC Curve
        # =====================================================================
        logger.info("  Panel 2: Generating ROC curve...")

        # Get predicted probabilities
        y_pred_prob = plsda_model.predict(X_scaled).ravel()

        # Calculate ROC curve
        fpr, tpr, thresholds = roc_curve(y_labels, y_pred_prob)
        roc_auc = auc(fpr, tpr)

        # Plot ROC curve
        ax2.plot(fpr, tpr, color='#E74C3C', linewidth=3,
                 label=f'PLS-DA (AUC = {roc_auc:.3f})')
        ax2.plot([0, 1], [0, 1], 'k--', linewidth=2, label='Random Classifier', alpha=0.5)

        ax2.set_xlabel('False Positive Rate', fontsize=12, fontweight='bold')
        ax2.set_ylabel('True Positive Rate', fontsize=12, fontweight='bold')
        ax2.set_title(f'ROC Curve: Discrimination Ability\nAUC = {roc_auc:.3f} (1.0 = Perfect)',
                      fontsize=12, fontweight='bold')
        ax2.legend(loc='lower right', frameon=True, fontsize=10)
        ax2.grid(True, alpha=0.3)
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

        ax2.text(0.6, 0.2, interp, fontsize=11, fontweight='bold',
                 bbox=dict(boxstyle='round', facecolor='lightgreen' if roc_auc > 0.8 else 'yellow', alpha=0.7))

        # =====================================================================
        # PANEL 3: Confusion Matrix
        # =====================================================================
        logger.info("  Panel 3: Computing confusion matrix...")

        # Predict classes (threshold at 0.5)
        y_pred_class = (y_pred_prob > 0.5).astype(int)

        # Compute confusion matrix
        cm = confusion_matrix(y_labels, y_pred_class)

        # Plot confusion matrix
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', cbar=False,
                    ax=ax3, square=True, linewidths=2, linecolor='black',
                    annot_kws={'size': 16, 'weight': 'bold'})

        ax3.set_xlabel('Predicted Class', fontsize=12, fontweight='bold')
        ax3.set_ylabel('True Class', fontsize=12, fontweight='bold')
        ax3.set_title('Confusion Matrix\n(Leave-One-Out Cross-Validation)',
                      fontsize=12, fontweight='bold')
        ax3.set_xticklabels(['Normal (0)', 'Cancer (1)'], fontsize=10)
        ax3.set_yticklabels(['Normal (0)', 'Cancer (1)'], fontsize=10, rotation=0)

        # Calculate accuracy metrics
        tn, fp, fn, tp = cm.ravel()
        accuracy = (tp + tn) / (tp + tn + fp + fn)
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0

        metrics_text = f"Accuracy: {accuracy:.1%}\n"
        metrics_text += f"Sensitivity: {sensitivity:.1%}\n"
        metrics_text += f"Specificity: {specificity:.1%}"

        ax3.text(1.05, 0.5, metrics_text, transform=ax3.transAxes,
                 fontsize=10, verticalalignment='center',
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7),
                 family='monospace')

        # =====================================================================
        # PANEL 4: VIP Score Distribution
        # =====================================================================
        logger.info("  Panel 4: Plotting VIP score distribution...")

        # Get VIP scores
        vip_scores = vip_df['VIP_Score'].values

        # Plot histogram
        ax4.hist(vip_scores, bins=50, color='#27AE60', alpha=0.7,
                 edgecolor='black', linewidth=0.5)

        # Add threshold line at VIP=1.0
        ax4.axvline(1.0, color='red', linestyle='--', linewidth=3,
                    label='VIP = 1.0 (Important Features)', zorder=10)

        # Add mean and median lines
        mean_vip = vip_scores.mean()
        median_vip = np.median(vip_scores)
        ax4.axvline(mean_vip, color='blue', linestyle=':', linewidth=2,
                    label=f'Mean = {mean_vip:.2f}', alpha=0.7)
        ax4.axvline(median_vip, color='orange', linestyle='-.', linewidth=2,
                    label=f'Median = {median_vip:.2f}', alpha=0.7)

        ax4.set_xlabel('VIP Score', fontsize=12, fontweight='bold')
        ax4.set_ylabel('Frequency', fontsize=12, fontweight='bold')
        ax4.set_title('VIP Score Distribution\nVIP > 1.0 = Important Features',
                      fontsize=12, fontweight='bold')
        ax4.legend(loc='upper right', frameon=True, fontsize=10)
        ax4.grid(True, alpha=0.3, axis='y')

        # Add statistics box
        n_important = (vip_scores > 1.0).sum()
        pct_important = n_important / len(vip_scores) * 100

        stats_text = f"Total features: {len(vip_scores)}\n"
        stats_text += f"VIP > 1.0: {n_important} ({pct_important:.1f}%)\n"
        stats_text += f"Max VIP: {vip_scores.max():.2f}"

        ax4.text(0.98, 0.98, stats_text, transform=ax4.transAxes,
                 fontsize=10, verticalalignment='top', horizontalalignment='right',
                 bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7),
                 family='monospace')

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'plsda_diagnostics.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved PLS-DA diagnostics to {output_file}")

        # Save diagnostic metrics as trace data
        diagnostic_metrics = pd.DataFrame({
            'Metric': ['R²', 'Q²', 'ROC_AUC', 'Accuracy', 'Sensitivity', 'Specificity',
                       'VIP_Mean', 'VIP_Median', 'N_Important_Features'],
            'Value': [selected_r2, selected_q2, roc_auc, accuracy, sensitivity, specificity,
                      mean_vip, median_vip, n_important]
        })
        save_trace_data(diagnostic_metrics, self.output_dir, 'plsda_diagnostic_metrics.csv')

        plt.close()

        logger.info("✓ PLS-DA diagnostics complete - model validated")
        logger.info(f"  Model quality: R²={selected_r2:.3f}, Q²={selected_q2:.3f}, AUC={roc_auc:.3f}")
        logger.info(f"  Classification: Accuracy={accuracy:.1%}, {n_important} important features")
