"""
VIP Score Plot Module using R/ggplot2 for pGlyco Auto Combine
Handles VIP score visualizations with R graphics
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
from ..utils import replace_empty_with_zero, get_sample_columns, save_trace_data

logger = logging.getLogger(__name__)


class VIPScorePlotRMixin:
    """Mixin class for R-based VIP score plots using ggplot2"""

    def _create_vip_plot_r(self, vip_data: pd.DataFrame, heatmap_data: pd.DataFrame,
                           title: str, ylabel: str, output_file: str):
        """
        Create VIP score plot with ggplot2

        Args:
            vip_data: DataFrame with columns [Feature, VIP_Score]
            heatmap_data: DataFrame with columns [Feature, Cancer, Normal]
            title: Plot title
            ylabel: Y-axis label
            output_file: Output file path
        """
        # Prepare combined data
        plot_data = vip_data.copy()
        plot_data['Cancer_High'] = (heatmap_data['Cancer'] > heatmap_data['Normal']).astype(int)
        plot_data['Normal_High'] = (heatmap_data['Normal'] > heatmap_data['Cancer']).astype(int)

        # Reverse order for plotting (top to bottom)
        plot_data = plot_data.iloc[::-1].reset_index(drop=True)
        plot_data['y_pos'] = range(len(plot_data))

        # Create R script for plotting
        r_script = f"""
        library(ggplot2)
        library(grid)

        # Prepare data
        df <- data.frame(
            Feature = c({', '.join([f'"{x}"' for x in plot_data['Feature']])}),
            VIP_Score = c({', '.join(map(str, plot_data['VIP_Score']))}),
            Cancer_High = c({', '.join(map(str, plot_data['Cancer_High']))}),
            Normal_High = c({', '.join(map(str, plot_data['Normal_High']))}),
            y_pos = c({', '.join(map(str, plot_data['y_pos']))})
        )

        # Create heatmap data in long format
        heatmap_long <- data.frame(
            y_pos = rep(df$y_pos, 2),
            Group = rep(c("Cancer", "Normal"), each = nrow(df)),
            Value = c(df$Cancer_High, df$Normal_High),
            x_pos = rep(c(0, 1), each = nrow(df))
        )

        # Create plot
        p <- ggplot() +
            # VIP scores (left panel)
            geom_point(data = df, aes(x = VIP_Score, y = y_pos),
                      size = 3, color = "#555555", stroke = 1.5, shape = 21, fill = "#555555") +

            # Feature names as text on the left side of plot body
            geom_text(data = df, aes(x = min(VIP_Score) - 0.15, y = y_pos, label = Feature),
                     hjust = 1, size = 2.5, fontface = "plain") +

            # Heatmap tiles (right panel)
            geom_tile(data = heatmap_long,
                     aes(x = max(df$VIP_Score) + 0.3 + x_pos * 0.15, y = y_pos, fill = factor(Value)),
                     width = 0.14, height = 0.9, color = "white", linewidth = 0.5) +

            # Color scale for heatmap
            scale_fill_manual(values = c("0" = "#3498DB", "1" = "#E74C3C"),
                            labels = c("Low", "High"),
                            name = "Relative\\nIntensity") +

            # Add Cancer/Normal labels above heatmap
            annotate("text", x = max(df$VIP_Score) + 0.3 + 0 * 0.15, y = max(df$y_pos) + 1.5,
                    label = "Cancer", size = 3) +
            annotate("text", x = max(df$VIP_Score) + 0.3 + 1 * 0.15, y = max(df$y_pos) + 1.5,
                    label = "Normal", size = 3) +

            # Scales and labels
            scale_x_continuous(limits = c(min(df$VIP_Score) - 0.8, max(df$VIP_Score) + 0.6),
                             expand = c(0, 0)) +
            scale_y_continuous(limits = c(-0.5, max(df$y_pos) + 2), expand = c(0, 0)) +

            labs(title = "{title}",
                 x = "VIP Score",
                 y = "") +

            # Theme with frame (X and Y axis lines)
            theme_minimal() +
            theme(
                plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                axis.title.x = element_text(size = 12),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.x = element_line(color = "black", linewidth = 0.5),
                axis.line.y = element_line(color = "black", linewidth = 0.5),
                panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
                panel.grid.major.y = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                legend.position = "right",
                legend.justification = "center",
                legend.title = element_text(size = 10, face = "bold"),
                legend.text = element_text(size = 9),
                legend.key.height = unit(1, "cm"),
                plot.margin = margin(10, 15, 10, 10)
            )

        # Save plot
        ggsave("{output_file}", plot = p, width = 10, height = 6, dpi = 300, bg = "white")
        """

        # Execute R script
        ro.r(r_script)
        logger.info(f"Saved R-based VIP score plot to {output_file}")

    def plot_vip_scores_glycopeptide_r(self, df: pd.DataFrame, vip_df: pd.DataFrame, top_n: int = 10):
        """
        Plot top VIP scores by glycopeptide using R/ggplot2

        Args:
            df: Annotated DataFrame
            vip_df: DataFrame with VIP scores by glycopeptide
            top_n: Number of top glycopeptides to show
        """
        top_n_data = vip_df.head(top_n).copy()

        # Create labels
        top_n_data['Feature'] = top_n_data['Peptide'] + ' | ' + top_n_data['GlycanComposition']

        # Get sample columns
        # Get sample columns (C1-C24, N1-N24)
        cancer_samples, normal_samples = get_sample_columns(df)
        sample_cols = cancer_samples + normal_samples

        # Prepare heatmap data
        heatmap_data = []
        for _, row in top_n_data.iterrows():
            mask = (df['Peptide'] == row['Peptide']) & (df['GlycanComposition'] == row['GlycanComposition'])
            if mask.sum() > 0:
                glycopeptide_row = df[mask].iloc[0]
                cancer_mean = replace_empty_with_zero(glycopeptide_row[cancer_samples]).mean()
                normal_mean = replace_empty_with_zero(glycopeptide_row[normal_samples]).mean()
                heatmap_data.append({'Feature': row['Feature'], 'Cancer': cancer_mean, 'Normal': normal_mean})
            else:
                heatmap_data.append({'Feature': row['Feature'], 'Cancer': 0, 'Normal': 0})

        heatmap_df = pd.DataFrame(heatmap_data)
        vip_plot_data = top_n_data[['Feature', 'VIP_Score']].copy()

        output_file = self.output_dir / 'vip_score_glycopeptide_r.png'
        self._create_vip_plot_r(vip_plot_data, heatmap_df,
                                f'Top {top_n} Glycopeptides by VIP Score',
                                'Glycopeptide', str(output_file))

    def plot_vip_scores_glycan_composition_r(self, df: pd.DataFrame, vip_df: pd.DataFrame, top_n: int = 10):
        """
        Plot VIP scores by GlycanComposition using R/ggplot2

        Args:
            df: Annotated DataFrame
            vip_df: DataFrame with all VIP scores by glycopeptide
            top_n: Number of top GlycanCompositions to show
        """
        # Group by GlycanComposition and get the max VIP score for each
        glycan_vip = vip_df.groupby('GlycanComposition')['VIP_Score'].max().nlargest(top_n).reset_index()
        glycan_vip['Feature'] = glycan_vip['GlycanComposition']

        # Get sample columns
        # Get sample columns (C1-C24, N1-N24)
        cancer_samples, normal_samples = get_sample_columns(df)
        sample_cols = cancer_samples + normal_samples

        # Prepare heatmap data
        heatmap_data = []
        for _, row in glycan_vip.iterrows():
            glycan_comp = row['GlycanComposition']
            mask = df['GlycanComposition'] == glycan_comp

            if mask.sum() > 0:
                glycan_rows = df[mask]
                cancer_total = 0
                normal_total = 0
                for _, gly_row in glycan_rows.iterrows():
                    cancer_total += replace_empty_with_zero(gly_row[cancer_samples]).sum()
                    normal_total += replace_empty_with_zero(gly_row[normal_samples]).sum()
                heatmap_data.append({'Feature': row['Feature'], 'Cancer': cancer_total, 'Normal': normal_total})
            else:
                heatmap_data.append({'Feature': row['Feature'], 'Cancer': 0, 'Normal': 0})

        heatmap_df = pd.DataFrame(heatmap_data)
        vip_plot_data = glycan_vip[['Feature', 'VIP_Score']].copy()

        output_file = self.output_dir / 'vip_score_glycan_composition_r.png'
        self._create_vip_plot_r(vip_plot_data, heatmap_df,
                                f'Top {top_n} Glycan Compositions by VIP Score',
                                'Glycan Composition', str(output_file))

    def plot_vip_scores_peptide_r(self, df: pd.DataFrame, vip_df: pd.DataFrame, top_n: int = 10):
        """
        Plot VIP scores by Peptide using R/ggplot2

        Args:
            df: Annotated DataFrame
            vip_df: DataFrame with all VIP scores
            top_n: Number of top peptides to show
        """
        # Group by Peptide and get the max VIP score for each
        peptide_vip = vip_df.groupby('Peptide')['VIP_Score'].max().nlargest(top_n).reset_index()
        peptide_vip['Feature'] = peptide_vip['Peptide']

        # Get sample columns
        # Get sample columns (C1-C24, N1-N24)
        cancer_samples, normal_samples = get_sample_columns(df)
        sample_cols = cancer_samples + normal_samples

        # Prepare heatmap data
        heatmap_data = []
        for _, row in peptide_vip.iterrows():
            peptide = row['Peptide']
            mask = df['Peptide'] == peptide

            if mask.sum() > 0:
                peptide_rows = df[mask]
                cancer_total = 0
                normal_total = 0
                for _, pep_row in peptide_rows.iterrows():
                    cancer_total += replace_empty_with_zero(pep_row[cancer_samples]).sum()
                    normal_total += replace_empty_with_zero(pep_row[normal_samples]).sum()
                heatmap_data.append({'Feature': row['Feature'], 'Cancer': cancer_total, 'Normal': normal_total})
            else:
                heatmap_data.append({'Feature': row['Feature'], 'Cancer': 0, 'Normal': 0})

        heatmap_df = pd.DataFrame(heatmap_data)
        vip_plot_data = peptide_vip[['Feature', 'VIP_Score']].copy()

        output_file = self.output_dir / 'vip_score_peptide_r.png'
        self._create_vip_plot_r(vip_plot_data, heatmap_df,
                                f'Top {top_n} Peptides by VIP Score',
                                'Peptide', str(output_file))

    def plot_vip_scores_peptide_grouped_r(self, df: pd.DataFrame, vip_df: pd.DataFrame, top_n: int = 10):
        """
        Plot VIP scores with peptide grouping (bracket notation for multiple glycoforms)

        Args:
            df: Annotated DataFrame
            vip_df: DataFrame with all VIP scores by glycopeptide
            top_n: Number of top peptides to show
        """
        # Get top peptides by max VIP score
        top_peptides = vip_df.groupby('Peptide')['VIP_Score'].max().nlargest(top_n).index.tolist()

        # Build plot data with grouping structure
        plot_rows = []
        y_pos = 0

        for peptide in reversed(top_peptides):  # Reverse for top-to-bottom plotting
            # Get all glycoforms for this peptide, sorted by VIP score descending
            peptide_glycoforms = vip_df[vip_df['Peptide'] == peptide].sort_values('VIP_Score', ascending=False).reset_index(drop=True)

            # Get sample columns (C1-C24, N1-N24)
            cancer_samples, normal_samples = get_sample_columns(df)

            # Add each glycoform
            for idx, row in peptide_glycoforms.iterrows():
                glycan_comp = row['GlycanComposition']
                vip_score = row['VIP_Score']

                # Get heatmap data
                mask = (df['Peptide'] == peptide) & (df['GlycanComposition'] == glycan_comp)
                if mask.sum() > 0:
                    glycopeptide_row = df[mask].iloc[0]
                    cancer_mean = replace_empty_with_zero(glycopeptide_row[cancer_samples]).mean()
                    normal_mean = replace_empty_with_zero(glycopeptide_row[normal_samples]).mean()
                else:
                    cancer_mean = 0
                    normal_mean = 0

                cancer_high = 1 if cancer_mean > normal_mean else 0
                normal_high = 1 if normal_mean > cancer_mean else 0

                plot_rows.append({
                    'Peptide': peptide,
                    'GlycanComposition': glycan_comp,
                    'VIP_Score': vip_score,
                    'y_pos': y_pos,
                    'Cancer_High': cancer_high,
                    'Normal_High': normal_high,
                    'is_first': idx == 0,
                    'is_last': idx == len(peptide_glycoforms) - 1,
                    'group_size': len(peptide_glycoforms)
                })
                y_pos += 1

            # Add spacing between peptide groups
            y_pos += 0.5

        plot_data = pd.DataFrame(plot_rows)

        # Create R script
        r_script = f"""
        library(ggplot2)
        library(grid)

        # Prepare main data
        df <- data.frame(
            Peptide = c({', '.join([f'"{x}"' for x in plot_data['Peptide']])}),
            GlycanComposition = c({', '.join([f'"{x}"' for x in plot_data['GlycanComposition']])}),
            VIP_Score = c({', '.join(map(str, plot_data['VIP_Score']))}),
            y_pos = c({', '.join(map(str, plot_data['y_pos']))}),
            Cancer_High = c({', '.join(map(str, plot_data['Cancer_High']))}),
            Normal_High = c({', '.join(map(str, plot_data['Normal_High']))}),
            is_first = c({', '.join(map(lambda x: 'TRUE' if x else 'FALSE', plot_data['is_first']))}),
            is_last = c({', '.join(map(lambda x: 'TRUE' if x else 'FALSE', plot_data['is_last']))}),
            group_size = c({', '.join(map(str, plot_data['group_size']))})
        )

        # Create bracket segments (only for peptides with multiple glycoforms)
        bracket_segments <- data.frame()
        for (peptide in unique(df$Peptide)) {{
            peptide_rows <- df[df$Peptide == peptide, ]
            if (nrow(peptide_rows) > 1) {{
                y_start <- min(peptide_rows$y_pos)
                y_end <- max(peptide_rows$y_pos)
                bracket_x <- min(df$VIP_Score) - 0.5

                # Vertical line connecting all glycoforms
                bracket_segments <- rbind(bracket_segments, data.frame(
                    x = bracket_x,
                    xend = bracket_x,
                    y = y_start,
                    yend = y_end
                ))

                # Horizontal lines to each glycoform
                for (i in 1:nrow(peptide_rows)) {{
                    bracket_segments <- rbind(bracket_segments, data.frame(
                        x = bracket_x,
                        xend = bracket_x + 0.15,
                        y = peptide_rows$y_pos[i],
                        yend = peptide_rows$y_pos[i]
                    ))
                }}
            }}
        }}

        # Create heatmap data in long format
        heatmap_long <- data.frame(
            y_pos = rep(df$y_pos, 2),
            Group = rep(c("Cancer", "Normal"), each = nrow(df)),
            Value = c(df$Cancer_High, df$Normal_High),
            x_pos = rep(c(0, 1), each = nrow(df))
        )

        # Create plot
        p <- ggplot()

        # Add bracket segments only if there are any (conditional)
        if (nrow(bracket_segments) > 0) {{
            p <- p + geom_segment(data = bracket_segments,
                                 aes(x = x, xend = xend, y = y, yend = yend),
                                 color = "#666666", linewidth = 0.4)
        }}

        p <- p +
            # VIP scores (dots)
            geom_point(data = df, aes(x = VIP_Score, y = y_pos),
                      size = 3, color = "#555555", stroke = 1.5, shape = 21, fill = "#555555") +

            # Peptide names (only for groups with multiple glycoforms or centered for single)
            geom_text(data = df[df$group_size > 1 & df$is_first, ],
                     aes(x = min(VIP_Score) - 0.8, y = y_pos, label = Peptide),
                     hjust = 1, size = 2.5, fontface = "bold", color = "#000000") +

            # Peptide names for single glycoform entries
            geom_text(data = df[df$group_size == 1, ],
                     aes(x = min(VIP_Score) - 0.8, y = y_pos, label = Peptide),
                     hjust = 1, size = 2.5, fontface = "bold", color = "#000000") +

            # GlycanComposition labels (for multiple glycoform groups only)
            geom_text(data = df[df$group_size > 1, ],
                     aes(x = min(VIP_Score) - 0.15, y = y_pos, label = GlycanComposition),
                     hjust = 1, size = 2.0, fontface = "plain", color = "#555555") +

            # Heatmap tiles (right panel)
            geom_tile(data = heatmap_long,
                     aes(x = max(df$VIP_Score) + 0.3 + x_pos * 0.15, y = y_pos, fill = factor(Value)),
                     width = 0.14, height = 0.9, color = "white", linewidth = 0.5) +

            # Color scale for heatmap
            scale_fill_manual(values = c("0" = "#3498DB", "1" = "#E74C3C"),
                            labels = c("Low", "High"),
                            name = "Relative\\nIntensity") +

            # Add Cancer/Normal labels above heatmap
            annotate("text", x = max(df$VIP_Score) + 0.3 + 0 * 0.15, y = max(df$y_pos) + 1.5,
                    label = "Cancer", size = 3) +
            annotate("text", x = max(df$VIP_Score) + 0.3 + 1 * 0.15, y = max(df$y_pos) + 1.5,
                    label = "Normal", size = 3) +

            # Scales and labels
            scale_x_continuous(limits = c(min(df$VIP_Score) - 1.5, max(df$VIP_Score) + 0.6),
                             expand = c(0, 0)) +
            scale_y_continuous(limits = c(-0.5, max(df$y_pos) + 2), expand = c(0, 0)) +

            labs(title = "Top {top_n} Peptides by VIP Score (Grouped by Glycoforms)",
                 x = "VIP Score",
                 y = "") +

            # Theme with frame (X and Y axis lines)
            theme_minimal() +
            theme(
                plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                axis.title.x = element_text(size = 12),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.x = element_line(color = "black", linewidth = 0.5),
                axis.line.y = element_line(color = "black", linewidth = 0.5),
                panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
                panel.grid.major.y = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                legend.position = "right",
                legend.justification = "center",
                legend.title = element_text(size = 10, face = "bold"),
                legend.text = element_text(size = 9),
                legend.key.height = unit(1, "cm"),
                plot.margin = margin(10, 15, 10, 10)
            )

        # Save plot
        ggsave("{self.output_dir / 'vip_score_peptide_grouped_r.png'}", plot = p, width = 12, height = 8, dpi = 300, bg = "white")
        """

        # Execute R script
        ro.r(r_script)
        logger.info(f"Saved R-based peptide grouped VIP score plot to {self.output_dir / 'vip_score_peptide_grouped_r.png'}")
