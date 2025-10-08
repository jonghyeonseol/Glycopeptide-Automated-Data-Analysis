"""
Interactive Dashboard Module for pGlyco Auto Combine
Generates interactive Plotly visualizations with drill-down capabilities

Phase 2.1: Post-doctoral level enhancement for publication-quality interactive figures
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from pathlib import Path
import logging
from typing import Tuple
import base64

logger = logging.getLogger(__name__)


class InteractiveDashboard:
    """Create interactive Plotly visualizations for glycoproteomics data"""

    def __init__(self, output_dir: str, colors: dict = None):
        """
        Initialize InteractiveDashboard

        Args:
            output_dir: Directory to save interactive HTML plots
            colors: Color scheme for glycan types
        """
        self.output_dir = Path(output_dir) / 'interactive'
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Match static plot color scheme
        self.colors = colors or {
            'HM': '#27AE60',        # Green - High Mannose
            'F': '#E74C3C',         # Red - Fucosylated
            'S': '#E91E63',         # Pink - Sialylated
            'SF': '#E67E22',        # Orange - Sialofucosylated
            'C/H': '#3498DB'        # Blue - Complex/Hybrid
        }

        self.group_colors = {
            'Cancer': '#E74C3C',    # Red
            'Normal': '#3498DB'     # Blue
        }

        logger.info(f"Interactive dashboard output: {self.output_dir}")

    def create_interactive_pca(self, pca_results: dict, df: pd.DataFrame = None) -> go.Figure:
        """
        Create interactive PCA plot with hover tooltips and confidence ellipses

        Features:
        - Hover: Sample ID, PC scores, group, glycan count
        - Click: Highlight sample
        - Zoom/pan: Explore sample clusters
        - 95% confidence ellipses

        Args:
            pca_results: PCA results from analyzer
            df: Optional DataFrame for additional metadata (glycan counts)

        Returns:
            Plotly Figure object
        """
        pca_df = pca_results['pca_df'].copy()
        explained_var = pca_results['explained_variance']

        # Add glycan count if df provided
        if df is not None:
            from src.utils import get_sample_columns
            cancer_samples, normal_samples = get_sample_columns(df)
            all_samples = cancer_samples + normal_samples

            # Count detected glycopeptides per sample
            glycan_counts = {}
            for sample in all_samples:
                if sample in df.columns:
                    glycan_counts[sample] = (df[sample] > 0).sum()

            pca_df['Glycan_Count'] = pca_df.index.map(glycan_counts)
        else:
            pca_df['Glycan_Count'] = 'N/A'

        # Create scatter plot
        fig = go.Figure()

        # Plot Cancer samples
        cancer_mask = pca_df['Group'] == 'Cancer'
        cancer_data = pca_df[cancer_mask]

        fig.add_trace(go.Scatter(
            x=cancer_data['PC1'],
            y=cancer_data['PC2'],
            mode='markers+text',
            name='Cancer',
            text=cancer_data.index,
            textposition='top center',
            textfont=dict(size=10, color=self.group_colors['Cancer']),
            marker=dict(
                size=12,
                color=self.group_colors['Cancer'],
                line=dict(color='black', width=1),
                opacity=0.8
            ),
            customdata=np.column_stack((
                cancer_data.index,
                cancer_data['PC1'].round(3),
                cancer_data['PC2'].round(3),
                cancer_data['Glycan_Count']
            )),
            hovertemplate='<b>%{customdata[0]}</b><br>' +
            'PC1: %{customdata[1]}<br>' +
            'PC2: %{customdata[2]}<br>' +
            'Glycans: %{customdata[3]}<br>' +
            '<extra></extra>'
        ))

        # Plot Normal samples
        normal_mask = pca_df['Group'] == 'Normal'
        normal_data = pca_df[normal_mask]

        fig.add_trace(go.Scatter(
            x=normal_data['PC1'],
            y=normal_data['PC2'],
            mode='markers+text',
            name='Normal',
            text=normal_data.index,
            textposition='top center',
            textfont=dict(size=10, color=self.group_colors['Normal']),
            marker=dict(
                size=12,
                color=self.group_colors['Normal'],
                line=dict(color='black', width=1),
                opacity=0.8,
                symbol='square'
            ),
            customdata=np.column_stack((
                normal_data.index,
                normal_data['PC1'].round(3),
                normal_data['PC2'].round(3),
                normal_data['Glycan_Count']
            )),
            hovertemplate='<b>%{customdata[0]}</b><br>' +
            'PC1: %{customdata[1]}<br>' +
            'PC2: %{customdata[2]}<br>' +
            'Glycans: %{customdata[3]}<br>' +
            '<extra></extra>'
        ))

        # Add 95% confidence ellipses
        for group_name, mask, color in [
            ('Cancer', cancer_mask, self.group_colors['Cancer']),
            ('Normal', normal_mask, self.group_colors['Normal'])
        ]:
            group_data = pca_df[mask]
            if len(group_data) >= 3:
                ellipse_x, ellipse_y = self._calculate_confidence_ellipse(
                    group_data['PC1'].values,
                    group_data['PC2'].values,
                    n_std=1.96  # 95% confidence
                )

                fig.add_trace(go.Scatter(
                    x=ellipse_x,
                    y=ellipse_y,
                    mode='lines',
                    name=f'{group_name} 95% CI',
                    line=dict(color=color, width=2, dash='dash'),
                    showlegend=True,
                    hoverinfo='skip'
                ))

        # Update layout
        fig.update_layout(
            title=dict(
                text='<b>Interactive PCA: Cancer vs Normal</b><br>' +
                     f'<sub>PC1: {explained_var[0] * 100:.2f}% | PC2: {explained_var[1] * 100:.2f}%</sub>',
                x=0.5,
                xanchor='center',
                font=dict(size=18)
            ),
            xaxis_title =f'PC1 ({explained_var[0] * 100:.2f}%)',
            yaxis_title =f'PC2 ({explained_var[1] * 100:.2f}%)',
            hovermode='closest',
            plot_bgcolor='white',
            width=1000,
            height=800,
            font=dict(size=12),
            legend=dict(
                x=1.02,
                y=1,
                xanchor='left',
                yanchor='top',
                bgcolor='rgba(255,255,255,0.9)',
                bordercolor='black',
                borderwidth=1
            )
        )

        # Add grid
        fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='lightgray', zeroline=True)
        fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgray', zeroline=True)

        # Save to HTML
        output_file = self.output_dir / 'interactive_pca.html'
        fig.write_html(output_file, include_plotlyjs='cdn')
        logger.info(f"Saved interactive PCA to {output_file}")

        return fig

    def create_interactive_volcano(self, volcano_df: pd.DataFrame,
                                   fdr_threshold: float = 0.05,
                                   fc_threshold: float = 1.5) -> go.Figure:
        """
        Create interactive volcano plot with drill-down capabilities

        Features:
        - Hover: Glycopeptide ID, FC, p-value, FDR, VIP score, glycan type
        - Click: Select glycopeptide for detailed view
        - Filter by glycan type using legend
        - Zoom to significant regions

        Args:
            volcano_df: Volcano plot data from VolcanoPlotMixin
            fdr_threshold: FDR significance threshold
            fc_threshold: Fold change threshold

        Returns:
            Plotly Figure object
        """
        # Ensure required columns exist
        if 'Glycopeptide' not in volcano_df.columns:
            volcano_df['Glycopeptide'] = volcano_df['Peptide'] + '_' + volcano_df['GlycanComposition']

        if '-Log10FDR' not in volcano_df.columns:
            volcano_df['-Log10FDR'] = -np.log10(volcano_df['FDR'])

        # Create figure
        fig = go.Figure()

        # Define regulation categories
        categories = {
            'Up in Cancer': {'color': self.group_colors['Cancer'], 'symbol': 'circle'},
            'Down in Cancer': {'color': self.group_colors['Normal'], 'symbol': 'square'},
            'Non-significant': {'color': '#95A5A6', 'symbol': 'diamond'}
        }

        # Plot by regulation status
        for regulation, style in categories.items():
            mask = volcano_df['Regulation'] == regulation
            subset = volcano_df[mask]

            if len(subset) == 0:
                continue

            # Prepare hover data
            hover_text = []
            for _, row in subset.iterrows():
                text = f"<b>{row['Glycopeptide']}</b><br>"
                text += f"Peptide: {row['Peptide']}<br>"
                text += f"Glycan: {row['GlycanComposition']}<br>"
                text += f"Type: {row.get('GlycanTypeCategory', 'N/A')}<br>"
                text += f"Log2FC: {row['Log2FC']:.3f}<br>"

                # Add CI if available (Phase 2.2)
                if 'Log2FC_CI_Lower' in row and not pd.isna(row['Log2FC_CI_Lower']):
                    text += f"95% CI: [{row['Log2FC_CI_Lower']:.3f}, {row['Log2FC_CI_Upper']:.3f}]<br>"

                text += f"FDR: {row['FDR']:.2e}<br>"
                text += f"-Log10FDR: {row['-Log10FDR']:.2f}<br>"
                text += f"VIP: {row.get('VIP_Score', 0):.2f}"
                hover_text.append(text)

            # Size by VIP score
            if 'VIP_Score' in subset.columns:
                sizes = 8 + (subset['VIP_Score'] / subset['VIP_Score'].max() * 12)
            else:
                sizes = 10

            fig.add_trace(go.Scatter(
                x=subset['Log2FC'],
                y=subset['-Log10FDR'],
                mode='markers',
                name=f"{regulation} (n={len(subset)})",
                marker=dict(
                    size=sizes,
                    color=style['color'],
                    symbol=style['symbol'],
                    line=dict(color='black', width=0.5),
                    opacity=0.7
                ),
                text=hover_text,
                hovertemplate='%{text}<extra></extra>'
            ))

        # Add threshold lines
        fig.add_hline(
            y=-np.log10(fdr_threshold),
            line=dict(color='gray', width=2, dash='dash'),
            annotation_text=f'FDR = {fdr_threshold}',
            annotation_position='right'
        )

        fig.add_vline(
            x=np.log2(fc_threshold),
            line=dict(color='gray', width=2, dash='dash'),
            annotation_text=f'FC = {fc_threshold}x',
            annotation_position='top'
        )

        fig.add_vline(
            x=-np.log2(fc_threshold),
            line=dict(color='gray', width=2, dash='dash'),
            annotation_text=f'FC = {fc_threshold}x',
            annotation_position='top'
        )

        # Calculate statistics
        n_up = len(volcano_df[volcano_df['Regulation'] == 'Up in Cancer'])
        n_down = len(volcano_df[volcano_df['Regulation'] == 'Down in Cancer'])
        n_ns = len(volcano_df[volcano_df['Regulation'] == 'Non-significant'])

        # Update layout
        fig.update_layout(
            title=dict(
                text='<b>Interactive Volcano Plot: Differential Expression</b><br>' +
                     f'<sub>Up: {n_up} | Down: {n_down} | NS: {n_ns} | ' +
                     f'Thresholds: FC>{fc_threshold}x, FDR<{fdr_threshold}</sub>',
                x=0.5,
                xanchor='center',
                font=dict(size=18)
            ),
            xaxis_title='Log2 Fold Change (Cancer / Normal)',
            yaxis_title='-Log10 FDR',
            hovermode='closest',
            plot_bgcolor='white',
            width=1200,
            height=900,
            font=dict(size=12),
            legend=dict(
                x=1.02,
                y=1,
                xanchor='left',
                yanchor='top',
                bgcolor='rgba(255,255,255,0.9)',
                bordercolor='black',
                borderwidth=1
            )
        )

        # Add grid
        fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='lightgray', zeroline=True)
        fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgray', zeroline=True)

        # Save to HTML
        output_file = self.output_dir / 'interactive_volcano.html'
        fig.write_html(output_file, include_plotlyjs='cdn')
        logger.info(f"Saved interactive volcano plot to {output_file}")

        return fig

    def create_interactive_vip_scores(self, vip_df: pd.DataFrame,
                                      df: pd.DataFrame,
                                      top_n: int = 30) -> go.Figure:
        """
        Create interactive VIP score plot with filtering

        Features:
        - Hover: Full glycopeptide info, VIP score, intensities
        - Filter by glycan type
        - Sort by VIP score or glycan type

        Args:
            vip_df: VIP scores DataFrame
            df: Annotated DataFrame with intensity data
            top_n: Number of top features to show

        Returns:
            Plotly Figure object
        """
        # Get top N features
        top_vip = vip_df.nlargest(top_n, 'VIP_Score').copy()

        # Merge with annotations
        top_vip = top_vip.merge(
            df[['Peptide', 'GlycanComposition', 'GlycanTypeCategory']].drop_duplicates(),
            on=['Peptide', 'GlycanComposition'],
            how='left'
        )

        # Calculate mean intensities for Cancer and Normal
        from src.utils import get_sample_columns
        cancer_samples, normal_samples = get_sample_columns(df)

        top_vip['Cancer_Mean'] = df.loc[
            top_vip.index, cancer_samples
        ].replace('', 0).astype(float).mean(axis=1)

        top_vip['Normal_Mean'] = df.loc[
            top_vip.index, normal_samples
        ].replace('', 0).astype(float).mean(axis=1)

        top_vip['FC'] = top_vip['Cancer_Mean'] / (top_vip['Normal_Mean'] + 1)

        # Create glycopeptide labels
        top_vip['Label'] = top_vip['Peptide'].str[:10] + '..._' + top_vip['GlycanComposition']

        # Sort by VIP score (descending)
        top_vip = top_vip.sort_values('VIP_Score', ascending=True)  # Ascending for horizontal bar

        # Create figure
        fig = go.Figure()

        # Plot by glycan type
        for glycan_type in ['HM', 'F', 'S', 'SF', 'C/H']:
            mask = top_vip['GlycanTypeCategory'] == glycan_type
            subset = top_vip[mask]

            if len(subset) == 0:
                continue

            # Prepare hover data
            hover_text = []
            for _, row in subset.iterrows():
                text = f"<b>{row['Label']}</b><br>"
                text += f"VIP: {row['VIP_Score']:.3f}<br>"
                text += f"Type: {glycan_type}<br>"
                text += f"Cancer: {row['Cancer_Mean']:.1e}<br>"
                text += f"Normal: {row['Normal_Mean']:.1e}<br>"
                text += f"FC: {row['FC']:.2f}x"
                hover_text.append(text)

            fig.add_trace(go.Bar(
                y=subset['Label'],
                x=subset['VIP_Score'],
                orientation='h',
                name=glycan_type,
                marker=dict(
                    color=self.colors.get(glycan_type, '#333333'),
                    line=dict(color='black', width=1)
                ),
                text=hover_text,
                hovertemplate='%{text}<extra></extra>'
            ))

        # Add VIP threshold line
        fig.add_vline(
            x=1.0,
            line=dict(color='red', width=2, dash='dash'),
            annotation_text='VIP = 1.0 (Important)',
            annotation_position='top right'
        )

        # Update layout
        fig.update_layout(
            title=dict(
                text=f'<b>Top {top_n} Features by VIP Score (PLS-DA)</b><br>' +
                '<sub>Variable Importance in Projection - Higher = More discriminative</sub>',
                x=0.5,
                xanchor='center',
                font=dict(size=18)
            ),
            xaxis_title='VIP Score',
            yaxis_title='Glycopeptide',
            hovermode='closest',
            plot_bgcolor='white',
            width=1200,
            height=max(800, top_n * 25),  # Scale height with number of features
            font=dict(size=11),
            barmode='stack',
            legend=dict(
                x=1.02,
                y=1,
                xanchor='left',
                yanchor='top',
                bgcolor='rgba(255,255,255,0.9)',
                bordercolor='black',
                borderwidth=1,
                title='Glycan Type'
            )
        )

        # Add grid
        fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')
        fig.update_yaxes(showgrid=False)

        # Save to HTML
        output_file = self.output_dir / 'interactive_vip_scores.html'
        fig.write_html(output_file, include_plotlyjs='cdn')
        logger.info(f"Saved interactive VIP scores to {output_file}")

        return fig

    def create_interactive_boxplot(self, df: pd.DataFrame,
                                   stat_results: pd.DataFrame = None) -> go.Figure:
        """
        Create interactive boxplot comparing Cancer vs Normal by glycan type

        Features:
        - Hover: Sample count, quartiles, outliers
        - Toggle glycan types
        - Show/hide individual points

        Args:
            df: Annotated DataFrame
            stat_results: Optional statistical test results

        Returns:
            Plotly Figure object
        """
        from src.utils import get_sample_columns, replace_empty_with_zero

        cancer_samples, normal_samples = get_sample_columns(df)

        # Prepare data for plotting
        plot_data = []

        for glycan_type in ['HM', 'F', 'S', 'SF', 'C/H']:
            mask = df['GlycanTypeCategory'] == glycan_type
            subset = df[mask]

            if len(subset) == 0:
                continue

            # Get intensities
            cancer_intensities = replace_empty_with_zero(subset[cancer_samples])
            normal_intensities = replace_empty_with_zero(subset[normal_samples])

            # Flatten and log transform
            cancer_values = np.log10(cancer_intensities.values.flatten() + 1)
            normal_values = np.log10(normal_intensities.values.flatten() + 1)

            # Remove zeros
            cancer_values = cancer_values[cancer_values > 0]
            normal_values = normal_values[normal_values > 0]

            # Add to plot data
            plot_data.append({
                'values': cancer_values,
                'group': 'Cancer',
                'glycan_type': glycan_type,
                'color': self.group_colors['Cancer']
            })

            plot_data.append({
                'values': normal_values,
                'group': 'Normal',
                'glycan_type': glycan_type,
                'color': self.group_colors['Normal']
            })

        # Create figure
        fig = go.Figure()

        # Plot boxplots
        for data in plot_data:
            fig.add_trace(go.Box(
                y=data['values'],
                name=f"{data['glycan_type']} - {data['group']}",
                marker=dict(color=data['color']),
                boxmean='sd',  # Show mean and std dev
                legendgroup=data['glycan_type'],
                showlegend=True
            ))

        # Update layout
        fig.update_layout(
            title=dict(
                text='<b>Interactive Boxplot: Glycan Type Intensities</b><br>' +
                     '<sub>Cancer vs Normal comparison by glycan type category</sub>',
                x=0.5,
                xanchor='center',
                font=dict(size=18)
            ),
            yaxis_title='Log10(Intensity + 1)',
            hovermode='closest',
            plot_bgcolor='white',
            width=1400,
            height=700,
            font=dict(size=12),
            boxmode='group',
            legend=dict(
                x=1.02,
                y=1,
                xanchor='left',
                yanchor='top',
                bgcolor='rgba(255,255,255,0.9)',
                bordercolor='black',
                borderwidth=1
            )
        )

        # Add grid
        fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgray')

        # Save to HTML
        output_file = self.output_dir / 'interactive_boxplot.html'
        fig.write_html(output_file, include_plotlyjs='cdn')
        logger.info(f"Saved interactive boxplot to {output_file}")

        return fig

    def _calculate_confidence_ellipse(self, x: np.ndarray, y: np.ndarray,
                                      n_std: float = 1.96, n_points: int = 100) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate 95% confidence ellipse coordinates

        Args:
            x, y: Data points
            n_std: Number of standard deviations (1.96 for 95% CI)
            n_points: Number of points to draw ellipse

        Returns:
            Tuple of (ellipse_x, ellipse_y) coordinates
        """
        # Calculate covariance
        cov = np.cov(x, y)

        # Eigendecomposition
        eigenvalues, eigenvectors = np.linalg.eig(cov)

        # Sort by eigenvalue
        order = eigenvalues.argsort()[::-1]
        eigenvalues = eigenvalues[order]
        eigenvectors = eigenvectors[:, order]

        # Calculate angle
        angle = np.degrees(np.arctan2(*eigenvectors[:, 0][::-1]))

        # Width and height
        width, height = 2 * n_std * np.sqrt(eigenvalues)

        # Generate ellipse points
        t = np.linspace(0, 2 * np.pi, n_points)
        ellipse_x_r = width / 2 * np.cos(t)
        ellipse_y_r = height / 2 * np.sin(t)

        # Rotation matrix
        angle_rad = np.radians(angle)
        R = np.array([
            [np.cos(angle_rad), -np.sin(angle_rad)],
            [np.sin(angle_rad), np.cos(angle_rad)]
        ])

        # Rotate and translate
        ellipse_coords = R @ np.array([ellipse_x_r, ellipse_y_r])
        ellipse_x = ellipse_coords[0, :] + np.mean(x)
        ellipse_y = ellipse_coords[1, :] + np.mean(y)

        # Close the ellipse
        ellipse_x = np.append(ellipse_x, ellipse_x[0])
        ellipse_y = np.append(ellipse_y, ellipse_y[0])

        return ellipse_x, ellipse_y

    def generate_all_interactive_plots(self, pca_results: dict, df: pd.DataFrame,
                                       vip_df: pd.DataFrame = None,
                                       volcano_df: pd.DataFrame = None,
                                       results_dir: Path = None):
        """
        Generate all interactive plots

        Phase 2.3: Embeds static PNG visualizations in index page

        Args:
            pca_results: PCA results from analyzer
            df: Annotated DataFrame
            vip_df: VIP scores DataFrame (optional)
            volcano_df: Volcano plot data (optional)
            results_dir: Path to Results directory for embedding static PNGs (optional)
        """
        logger.info("Generating all interactive visualizations...")

        # PCA
        self.create_interactive_pca(pca_results, df)

        # Boxplot
        self.create_interactive_boxplot(df)

        # VIP scores (if available)
        if vip_df is not None:
            self.create_interactive_vip_scores(vip_df, df, top_n=30)

        # Volcano plot (if available)
        if volcano_df is not None:
            self.create_interactive_volcano(volcano_df)

        # Create index page with embedded static PNGs (Phase 2.3)
        self._create_index_page(results_dir=results_dir)

        if results_dir is not None:
            logger.info("✓ Static PNG visualizations embedded as base64 for full portability")

        logger.info(f"All interactive visualizations saved to {self.output_dir}")

    def _png_to_base64(self, png_path: Path) -> str:
        """
        Convert PNG file to base64 string for embedding

        Args:
            png_path: Path to PNG file

        Returns:
            Base64 encoded string
        """
        try:
            with open(png_path, 'rb') as f:
                png_data = f.read()
                base64_str = base64.b64encode(png_data).decode('utf-8')
                return f"data:image/png;base64,{base64_str}"
        except Exception as e:
            logger.warning(f"Failed to encode {png_path}: {e}")
            return ""

    def _create_index_page(self, results_dir: Path = None):
        """
        Create an index.html landing page with links to all interactive plots

        Phase 2.3: Embeds static PNG visualizations as base64 for full portability

        Args:
            results_dir: Path to Results directory containing static PNGs
        """
        # Build static visualization gallery with base64 embedding
        static_gallery_html = ""

        if results_dir is not None:
            # Key static plots to embed
            key_static_plots = [
                ('pca_plot.png', 'PCA Analysis', 'Sample clustering with 95% confidence ellipses'),
                ('volcano_plot.png', 'Volcano Plot', 'Differential expression with bootstrap CIs'),
                ('heatmap_top20_main.png', 'Top 20 Heatmap', 'Hierarchical clustering of top glycopeptides'),
                ('boxplot_glycan_types.png', 'Boxplot Analysis', 'Glycan type intensities with effect sizes'),
                ('vip_score_glycopeptide_r.png', 'VIP Scores', 'Top discriminative features from PLS-DA'),
                ('pie_chart_glycan_types_enhanced.png', 'Glycan Distribution', 'Cancer vs Normal with fold changes')
            ]

            static_cards = []
            for filename, title, description in key_static_plots:
                png_path = results_dir / filename
                if png_path.exists():
                    # Embed as base64 for portability
                    base64_img = self._png_to_base64(png_path)
                    if base64_img:
                        card_html = """
        <div class="plot-card">
            <h3>{title}</h3>
            <p>{description}</p>
            <img src="{base64_img}" alt="{title}" style="width:100%; border-radius:4px; margin:10px 0;">
            <span class="static-label">Static PNG (Base64 embedded)</span>
        </div>"""
                        static_cards.append(card_html)

            if static_cards:
                static_gallery_html = """
    <h2>📸 Key Static Visualizations (Embedded)</h2>
    <p class="info-text">High-resolution static plots embedded as base64 for full portability.
    All original PNG files available in parent Results/ directory.</p>
    <div class="plot-grid">
        {''.join(static_cards)}
    </div>
"""

        html_content = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>pGlyco Interactive Dashboard</title>
    <style>
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }
        h1 {
            color: #2c3e50;
            text-align: center;
            border-bottom: 3px solid #3498db;
            padding-bottom: 15px;
        }
        h2 {
            color: #34495e;
            margin-top: 30px;
        }
        .plot-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(350px, 1fr));
            gap: 20px;
            margin: 30px 0;
        }
        .plot-card {
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            transition: transform 0.2s, box-shadow 0.2s;
        }
        .plot-card:hover {
            transform: translateY(-5px);
            box-shadow: 0 4px 8px rgba(0,0,0,0.2);
        }
        .plot-card h3 {
            color: #2980b9;
            margin-top: 0;
        }
        .plot-card p {
            color: #7f8c8d;
            font-size: 14px;
            line-height: 1.6;
        }
        .plot-card a {
            display: inline-block;
            margin-top: 10px;
            padding: 10px 20px;
            background-color: #3498db;
            color: white;
            text-decoration: none;
            border-radius: 5px;
            transition: background-color 0.2s;
        }
        .plot-card a:hover {
            background-color: #2980b9;
        }
        .info-box {
            background-color: #e8f4f8;
            border-left: 4px solid #3498db;
            padding: 15px;
            margin: 20px 0;
            border-radius: 4px;
        }
        .info-text {
            color: #7f8c8d;
            font-size: 14px;
            margin: 10px 0 20px 0;
        }
        .static-label {
            display: inline-block;
            background-color: #27ae60;
            color: white;
            padding: 4px 8px;
            border-radius: 3px;
            font-size: 11px;
            font-weight: bold;
            margin-top: 5px;
        }
        footer {
            text-align: center;
            margin-top: 50px;
            padding-top: 20px;
            border-top: 2px solid #bdc3c7;
            color: #7f8c8d;
        }
    </style>
</head>
<body>
    <h1>🔬 pGlyco Interactive Dashboard</h1>

    <div class="info-box">
        <strong>Note:</strong> All interactive plots are standalone HTML files using Plotly.js.
        Features include hover tooltips, zoom/pan, legend filtering, and data export.
    </div>

    <h2>📊 Available Interactive Visualizations</h2>

    <div class="plot-grid">
        <div class="plot-card">
            <h3>Interactive PCA</h3>
            <p>Explore sample clustering with hover tooltips showing sample ID, PC scores, and glycan counts.
            Includes 95% confidence ellipses for Cancer vs Normal groups.</p>
            <a href="interactive_pca.html" target="_blank">Open PCA Plot →</a>
        </div>

        <div class="plot-card">
            <h3>Interactive Volcano Plot</h3>
            <p>Differential expression analysis with drill-down capabilities. Hover to see glycopeptide details,
            FC, FDR, and VIP scores. Click legend to filter by regulation status.</p>
            <a href="interactive_volcano.html" target="_blank">Open Volcano Plot →</a>
        </div>

        <div class="plot-card">
            <h3>Interactive VIP Scores</h3>
            <p>Top discriminative features from PLS-DA analysis. Filter by glycan type using the legend.
            Hover for detailed intensity comparisons.</p>
            <a href="interactive_vip_scores.html" target="_blank">Open VIP Scores →</a>
        </div>

        <div class="plot-card">
            <h3>Interactive Boxplot</h3>
            <p>Compare glycan type intensities between Cancer and Normal groups. Toggle categories
            and view statistical distributions with outliers.</p>
            <a href="interactive_boxplot.html" target="_blank">Open Boxplot →</a>
        </div>
    </div>

    {static_gallery_html}

    <h2>💡 How to Use</h2>
    <ul>
        <li><strong>Hover:</strong> Move your mouse over data points to see detailed information</li>
        <li><strong>Zoom:</strong> Click and drag to zoom into regions of interest. Double-click to reset.</li>
        <li><strong>Pan:</strong> Hold shift and drag to pan across the plot</li>
        <li><strong>Filter:</strong> Click legend items to show/hide categories</li>
        <li><strong>Export:</strong> Use the camera icon in the top-right to download as PNG</li>
    </ul>

    <footer>
        <p>Generated by pGlyco Auto Combine Pipeline - Interactive Dashboard Module</p>
        <p>Post-doctoral Enhancement Phase 2.1</p>
    </footer>
</body>
</html>
"""

        index_file = self.output_dir / 'index.html'
        with open(index_file, 'w') as f:
            f.write(html_content.format(static_gallery_html=static_gallery_html))

        logger.info(f"Created dashboard index at {index_file}")


if __name__ == "__main__":
    # Test module
    pass
