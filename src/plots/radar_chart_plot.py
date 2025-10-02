"""
Radar Chart Plot Module for pGlyco Auto Combine
Visualizes glycan type profiles as radar charts
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import logging
from utils import replace_empty_with_zero, save_trace_data

logger = logging.getLogger(__name__)


class RadarChartPlotMixin:
    """Mixin class for radar/spider chart visualization"""

    def plot_radar_chart(self, df: pd.DataFrame, figsize: tuple = (12, 10)):
        """
        Create radar chart comparing glycan profiles between Cancer and Normal

        Args:
            df: Annotated DataFrame with intensity data
            figsize: Figure size (width, height)
        """
        # Get sample columns
        cancer_samples = [col for col in df.columns if col.startswith('C') and col[1:].isdigit()]
        normal_samples = [col for col in df.columns if col.startswith('N') and col[1:].isdigit()]

        # Calculate percentages for each category
        categories = []
        cancer_values = []
        normal_values = []

        # 1. Sialylation percentage (based on intensity)
        if 'IsSialylated' in df.columns:
            sial_mask = df['IsSialylated']
        else:
            sial_mask = df['Sialylation'] == 'Sialylated'

        cancer_sial_total = replace_empty_with_zero(df[sial_mask][cancer_samples]).sum().sum()
        cancer_total = replace_empty_with_zero(df[cancer_samples]).sum().sum()
        cancer_sial_pct = (cancer_sial_total / cancer_total * 100) if cancer_total > 0 else 0

        normal_sial_total = replace_empty_with_zero(df[sial_mask][normal_samples]).sum().sum()
        normal_total = replace_empty_with_zero(df[normal_samples]).sum().sum()
        normal_sial_pct = (normal_sial_total / normal_total * 100) if normal_total > 0 else 0

        categories.append('Sialylated\n(%)')
        cancer_values.append(cancer_sial_pct)
        normal_values.append(normal_sial_pct)

        # 2. Fucosylation percentage
        if 'IsFucosylated' in df.columns:
            fuc_mask = df['IsFucosylated']
        else:
            fuc_mask = df['Fucosylation'] == 'Fucosylated'

        cancer_fuc_total = replace_empty_with_zero(df[fuc_mask][cancer_samples]).sum().sum()
        cancer_fuc_pct = (cancer_fuc_total / cancer_total * 100) if cancer_total > 0 else 0

        normal_fuc_total = replace_empty_with_zero(df[fuc_mask][normal_samples]).sum().sum()
        normal_fuc_pct = (normal_fuc_total / normal_total * 100) if normal_total > 0 else 0

        categories.append('Fucosylated\n(%)')
        cancer_values.append(cancer_fuc_pct)
        normal_values.append(normal_fuc_pct)

        # 3. High Mannose percentage
        hm_mask = df['PrimaryClassification'] == 'High Mannose'
        cancer_hm_total = replace_empty_with_zero(df[hm_mask][cancer_samples]).sum().sum()
        cancer_hm_pct = (cancer_hm_total / cancer_total * 100) if cancer_total > 0 else 0

        normal_hm_total = replace_empty_with_zero(df[hm_mask][normal_samples]).sum().sum()
        normal_hm_pct = (normal_hm_total / normal_total * 100) if normal_total > 0 else 0

        categories.append('High Mannose\n(%)')
        cancer_values.append(cancer_hm_pct)
        normal_values.append(normal_hm_pct)

        # 4. Complex/Hybrid percentage
        ch_mask = df['PrimaryClassification'] == 'ComplexHybrid'
        cancer_ch_total = replace_empty_with_zero(df[ch_mask][cancer_samples]).sum().sum()
        cancer_ch_pct = (cancer_ch_total / cancer_total * 100) if cancer_total > 0 else 0

        normal_ch_total = replace_empty_with_zero(df[ch_mask][normal_samples]).sum().sum()
        normal_ch_pct = (normal_ch_total / normal_total * 100) if normal_total > 0 else 0

        categories.append('Complex/Hybrid\n(%)')
        cancer_values.append(cancer_ch_pct)
        normal_values.append(normal_ch_pct)

        # 5. Sialofucosylated percentage (both modifications)
        sialofuc_mask = df['SecondaryClassification'] == 'Sialofucosylated'
        cancer_sf_total = replace_empty_with_zero(df[sialofuc_mask][cancer_samples]).sum().sum()
        cancer_sf_pct = (cancer_sf_total / cancer_total * 100) if cancer_total > 0 else 0

        normal_sf_total = replace_empty_with_zero(df[sialofuc_mask][normal_samples]).sum().sum()
        normal_sf_pct = (normal_sf_total / normal_total * 100) if normal_total > 0 else 0

        categories.append('Sialofucosylated\n(%)')
        cancer_values.append(cancer_sf_pct)
        normal_values.append(normal_sf_pct)

        # Number of variables
        num_vars = len(categories)

        # Compute angle for each axis
        angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()

        # Complete the circle
        cancer_values += cancer_values[:1]
        normal_values += normal_values[:1]
        angles += angles[:1]

        # Create plot
        fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(projection='polar'))

        # Plot data
        ax.plot(angles, cancer_values, 'o-', linewidth=2, label='Cancer', color='#E74C3C')
        ax.fill(angles, cancer_values, alpha=0.25, color='#E74C3C')

        ax.plot(angles, normal_values, 'o-', linewidth=2, label='Normal', color='#3498DB')
        ax.fill(angles, normal_values, alpha=0.25, color='#3498DB')

        # Fix axis to go in the right order and start at 12 o'clock
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)

        # Set category labels
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(categories, fontsize=11, fontweight='bold')

        # Set y-axis limits (0-100%)
        ax.set_ylim(0, 100)
        ax.set_yticks([20, 40, 60, 80, 100])
        ax.set_yticklabels(['20%', '40%', '60%', '80%', '100%'], fontsize=9)

        # Add grid
        ax.grid(True, alpha=0.3)

        # Add legend
        ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize=12)

        # Title
        plt.title('Glycan Profile Comparison: Cancer vs Normal\n(Radar Chart)',
                 fontsize=14, fontweight='bold', pad=20)

        plt.tight_layout()

        # Save plot
        output_file = self.output_dir / 'radar_chart_glycan_profile.png'
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"Saved radar chart to {output_file}")

        # Prepare trace data (before completing the circle)
        radar_data = pd.DataFrame({
            'Category': categories,
            'Cancer_Percentage': cancer_values[:-1],  # Remove duplicate first element
            'Normal_Percentage': normal_values[:-1]   # Remove duplicate first element
        })

        save_trace_data(radar_data, self.output_dir, 'radar_chart_glycan_profile_data.csv')

        plt.close()

        return radar_data
