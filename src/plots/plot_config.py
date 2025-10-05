"""
Standardized Plot Configuration for pGlyco Auto Combine
Inspired by GraphPad Prism design principles for publication-quality figures

DESIGN PRINCIPLES FOR SCIENTIFIC VISUALIZATION:
==============================================

1. READABILITY FIRST
   - Minimum font size: 10pt for body text, 12pt for labels
   - High contrast text (black on white backgrounds)
   - Adequate spacing between elements

2. DATA-TO-INK RATIO
   - Maximize data, minimize decoration
   - Remove unnecessary gridlines, borders, backgrounds
   - Use white space effectively

3. VISUAL HIERARCHY
   - Title > Axis Labels > Tick Labels > Annotations
   - Important data points should be visually prominent
   - Use size, color, and position to guide attention

4. COLOR USAGE
   - Maximum 5-7 distinct colors per plot
   - Use color meaningfully (not decoratively)
   - Ensure colorblind accessibility
   - Maintain consistency across related plots

5. LABEL CLARITY
   - Never overlap labels with data or other labels
   - Use abbreviations when necessary
   - Position labels to minimize eye movement
   - Ensure all text is horizontal or vertically aligned

6. CONSISTENCY
   - Same colors mean same things across all plots
   - Consistent font sizes for similar elements
   - Uniform spacing and margins
"""

# ==============================================================================
# Font Sizes - PUBLICATION QUALITY
# Designed for 300 DPI output, readable in printed journals
# ==============================================================================

TITLE_SIZE = 16  # Increased for better visibility
AXIS_LABEL_SIZE = 14  # Increased - must be clearly readable
TICK_LABEL_SIZE = 12  # Increased - critical for data interpretation
LEGEND_SIZE = 12  # Increased for accessibility
LEGEND_TITLE_SIZE = 13
ANNOTATION_SIZE = 11  # Increased - annotations must be readable

# Font Weights (Prism uses regular weight more often for cleaner look)
TITLE_WEIGHT = 'bold'
AXIS_LABEL_WEIGHT = 'normal'  # Changed from bold - cleaner
LEGEND_TITLE_WEIGHT = 'bold'

# ==============================================================================
# Colors - Group Comparison (Cancer vs Normal)
# Prism-inspired: Bold, saturated, high contrast
# ==============================================================================

# Using Prism's bold, professional colors
COLOR_CANCER = '#FF0000'  # Pure red - bold, unmistakable (from Prism "colors" palette)
COLOR_NORMAL = '#0072B2'  # Science blue - professional, colorblind-safe

GROUP_PALETTE = {
    'Cancer': COLOR_CANCER,
    'Normal': COLOR_NORMAL
}

# ==============================================================================
# Colors - Glycan Types (Scientific Color Scheme)
# ==============================================================================

# SCIENTIFIC COLOR SCHEME - Biologically meaningful colors
# Green = High Mannose (early biosynthesis)
# Red = Fucosylated (terminal modification)
# Pink = Sialylated (terminal modification)
# Orange = Sialofucosylated (both terminal modifications)
# Blue = Complex/Hybrid (processed structures)

GLYCAN_COLORS = {
    'HM': '#2ECC71',    # Green - High Mannose
    'F': '#E74C3C',     # Red - Fucosylated
    'S': '#FF69B4',     # Pink - Sialylated
    'SF': '#FF8C00',    # Orange - Sialofucosylated
    'C/H': '#3498DB'    # Blue - Complex/Hybrid
}

# Legacy glycan type colors (for backward compatibility)
LEGACY_GLYCAN_COLORS = {
    'Non': '#95A5A6',       # Gray - Non-modified
    'Sialylated': '#FF69B4',   # Pink - Sialylated
    'Fucosylated': '#E74C3C',  # Red - Fucosylated
    'Both': '#FF8C00'          # Orange - Sialofucosylated
}

# Extended category colors (scientific scheme)
EXTENDED_CATEGORY_COLORS = {
    'HM': '#2ECC71',              # Green - High Mannose
    'High mannose': '#2ECC71',    # Alias for consistency
    'High Mannose': '#2ECC71',    # Alias (capitalized)
    'C/H': '#3498DB',             # Blue - Complex/Hybrid
    'Complex/Hybrid': '#3498DB',  # Alias
    'ComplexHybrid': '#3498DB',   # Alias (no space/slash)
    'Fucosylated': '#E74C3C',     # Red - Fucosylated
    'Sialylated': '#FF69B4',      # Pink - Sialylated
    'Sialofucosylated': '#FF8C00', # Orange - Sialofucosylated
    'Both': '#FF8C00',            # Alias for Sialofucosylated
    'Truncated': '#95A5A6',       # Gray - Truncated/Other
    'Other': '#95A5A6'            # Gray - Other
}

# ==============================================================================
# Plot Settings (Prism-inspired)
# ==============================================================================

DPI = 300  # Publication quality

# Grid settings (Prism style: minimal, subtle)
GRID_ALPHA = 0.15  # Much lighter (was 0.3)
GRID_LINESTYLE = '-'  # Solid, not dashed (Prism style)
GRID_LINEWIDTH = 0.5
GRID_COLOR = '#E0E0E0'  # Light gray, not pure gray

# Legend settings (Prism style: clean, no shadow)
LEGEND_FRAMEON = True
LEGEND_FANCYBOX = False  # Square corners (Prism style)
LEGEND_SHADOW = False     # No shadow (cleaner)
LEGEND_FRAMEALPHA = 1.0   # Solid background (was 0.9)
LEGEND_EDGECOLOR = '#000000'  # Black border
LEGEND_FRAMEWIDTH = 1.0

# ==============================================================================
# Plot-Specific Settings (Prism-inspired: Bolder data elements)
# ==============================================================================

# Boxplot (Prism style: larger, bolder)
BOXPLOT_WIDTH = 0.7  # Wider boxes (was 0.6)
BOXPLOT_FIGSIZE = (10, 6)  # Slightly more compact
BOXPLOT_EXTENDED_FIGSIZE = (12, 6)
BOXPLOT_LINEWIDTH = 1.5  # Thicker box lines (Prism style)
BOXPLOT_FLIERSIZE = 6  # Larger outlier markers

# Histogram
HISTOGRAM_FIGSIZE = (16, 10)  # More compact (was 20x12)
HISTOGRAM_X_ROTATION = 90
HISTOGRAM_X_HA = 'right'
HISTOGRAM_BAR_EDGEWIDTH = 1.0  # Thicker bar edges

# ==============================================================================
# VIP Score Plots (R/ggplot2) - MAXIMUM VISIBILITY
# ==============================================================================
VIP_FEATURE_NAME_SIZE = 5.0  # R ggplot2 units - Balanced size for full-length labels
VIP_POINT_SIZE_MIN = 8       # Minimum point size (scaled by VIP)
VIP_POINT_SIZE_MAX = 14      # Maximum point size (very prominent)
VIP_POINT_STROKE = 3.5       # Very thick stroke for definition
VIP_HEATMAP_WIDTH = 0.30     # Extra wide tiles - easy to distinguish
VIP_HEATMAP_HEIGHT = 0.95    # Full height
VIP_GROUP_LABEL_SIZE = 6.5   # Large bold group labels (Cancer/Normal)
VIP_FIGURE_WIDTH = 20  # Reduced for API 2000px limit
VIP_FIGURE_HEIGHT = 8        # Taller for better proportions
VIP_LEFT_MARGIN_EXPAND = 0.3  # Extra space on left for labels outside box

# ==============================================================================
# PCA Plot - PUBLICATION QUALITY
# ==============================================================================
PCA_FIGSIZE = (10, 8)
PCA_LABEL_FONTSIZE = 10      # Larger (was 9)
PCA_LABEL_BBOX_PAD = 0.4     # More padding
PCA_POINT_SIZE = 150         # Larger scatter points
PCA_POINT_LINEWIDTH = 2.0    # Thicker edges
PCA_POINT_ALPHA = 0.8

# ==============================================================================
# Volcano Plot - MAXIMUM VISIBILITY
# ==============================================================================
VOLCANO_FIGSIZE = (14, 11)   # Extra large canvas for labels
VOLCANO_POINT_SIZE = 120     # Extra large base size
VOLCANO_POINT_ALPHA = 0.8    # High visibility
VOLCANO_THRESHOLD_LINEWIDTH = 3.0  # Thicker threshold lines
VOLCANO_THRESHOLD_ALPHA = 0.6
VOLCANO_POINT_EDGEWIDTH = 1.2  # Strong edges
VOLCANO_LABEL_FONTSIZE = 14    # Extra large labels
VOLCANO_LABEL_WEIGHT = 'bold'  # Bold for emphasis
VOLCANO_LABEL_PADDING = 0.6    # More padding in label boxes
VOLCANO_LABEL_LINEWIDTH = 2.5  # Thick label borders
VOLCANO_MAX_LABELS = 3         # Maximum labels to show (avoid overcrowding)

# Heatmap
HEATMAP_FIGSIZE = (12, 9)  # More compact
HEATMAP_CBAR_LABEL_SIZE = 11

# ==============================================================================
# Statistical Annotation Settings (MetaboAnalyst/Prism style)
# ==============================================================================
ERROR_BAR_CAPSIZE = 4  # Error bar cap width
ERROR_BAR_LINEWIDTH = 1.5  # Error bar line thickness
STAT_BRACKET_LINEWIDTH = 2.0  # Statistical comparison bracket thickness
SIGNIFICANCE_MARKER_SIZE = 14  # Size of *, **, *** markers

# ==============================================================================
# Utility Functions
# ==============================================================================

def apply_standard_axis_style(ax, xlabel=None, ylabel=None, title=None, grid=True):
    """
    Apply Prism-inspired standardized styling to matplotlib axes

    Args:
        ax: Matplotlib axes object
        xlabel: X-axis label text
        ylabel: Y-axis label text
        title: Plot title text
        grid: Whether to show gridlines (Prism style: minimal, only Y-axis)
    """
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=AXIS_LABEL_SIZE, fontweight=AXIS_LABEL_WEIGHT)

    if ylabel:
        ax.set_ylabel(ylabel, fontsize=AXIS_LABEL_SIZE, fontweight=AXIS_LABEL_WEIGHT)

    if title:
        ax.set_title(title, fontsize=TITLE_SIZE, fontweight=TITLE_WEIGHT, pad=15)

    # Tick label sizes and style (Prism: larger ticks pointing outward)
    ax.tick_params(axis='both', which='major', labelsize=TICK_LABEL_SIZE,
                  width=1.2, length=6, direction='out')

    # Axis line width (Prism style: bolder axis lines)
    for spine in ['top', 'right', 'bottom', 'left']:
        ax.spines[spine].set_linewidth(1.2)

    # Hide top and right spines (Prism style: cleaner look)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Grid (Prism style: very subtle, only Y-axis, horizontal lines)
    if grid:
        ax.grid(True, alpha=GRID_ALPHA, linestyle=GRID_LINESTYLE,
                linewidth=GRID_LINEWIDTH, color=GRID_COLOR, axis='y', zorder=0)
        ax.set_axisbelow(True)  # Grid behind data


def apply_standard_legend(ax, **kwargs):
    """
    Apply Prism-inspired standardized legend styling
    Legend positioned OUTSIDE plot area by default (non-interruptive)

    Args:
        ax: Matplotlib axes object
        **kwargs: Additional legend parameters
    """
    legend_params = {
        'fontsize': LEGEND_SIZE,
        'title_fontsize': LEGEND_TITLE_SIZE,
        'frameon': LEGEND_FRAMEON,
        'fancybox': LEGEND_FANCYBOX,
        'shadow': LEGEND_SHADOW,
        'framealpha': LEGEND_FRAMEALPHA,
        'edgecolor': LEGEND_EDGECOLOR,
        'facecolor': 'white',
        'loc': 'upper left',  # Default location
        'bbox_to_anchor': (1.02, 1)  # Outside plot area (right side, top aligned)
    }
    legend_params.update(kwargs)

    legend = ax.legend(**legend_params)
    if legend:
        legend.get_frame().set_linewidth(LEGEND_FRAMEWIDTH)
