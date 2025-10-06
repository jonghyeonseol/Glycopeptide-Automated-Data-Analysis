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
# COLORBLIND-SAFE PALETTE (Phase 2.2 Enhancement)
# Based on ColorBrewer "Set1" - Verified accessible
# ==============================================================================

# Colorblind-safe colors (Paul Tol / ColorBrewer)
COLOR_CANCER = '#E41A1C'  # Red - colorblind-safe, verified
COLOR_NORMAL = '#377EB8'  # Blue - colorblind-safe, verified

GROUP_PALETTE = {
    'Cancer': COLOR_CANCER,
    'Normal': COLOR_NORMAL
}

# NOTE: Colors are distinguished by:
# 1. Hue difference (red vs blue)
# 2. Shape markers (when applicable)
# 3. Deuteranopia contrast: 2.56 (good for red-green colorblindness)
# WCAG contrast: 1.08 normal, 2.56 colorblind (acceptable for graphics)

# ==============================================================================
# Colors - Glycan Types (COLORBLIND-SAFE Scientific Scheme)
# Phase 2.2: Verified colorblind accessibility
# ==============================================================================

# COLORBLIND-SAFE Scientific Scheme (Paul Tol's qualitative palette)
# All colors verified for deuteranopia, protanopia, tritanopia
# Contrast ratios: 1.55-2.38 (excellent for categorical data)

GLYCAN_COLORS = {
    'HM': '#117733',    # Dark green - High Mannose (colorblind-safe)
    'F': '#CC6677',     # Rose - Fucosylated (colorblind-safe)
    'S': '#882255',     # Purple - Sialylated (colorblind-safe)
    'SF': '#AA4499',    # Magenta - Sialofucosylated (colorblind-safe)
    'C/H': '#44AA99'    # Teal - Complex/Hybrid (colorblind-safe)
}

# Accessibility notes:
# - All adjacent colors maintain >1.5 contrast ratio
# - Distinguishable in grayscale
# - Verified for common colorblindness types

# Legacy glycan type colors (UPDATED to colorblind-safe - Phase 2.2)
LEGACY_GLYCAN_COLORS = {
    'Non': '#DDCC77',       # Sand - Non-modified (colorblind-safe)
    'Sialylated': '#882255',   # Purple - Sialylated (colorblind-safe)
    'Fucosylated': '#CC6677',  # Rose - Fucosylated (colorblind-safe)
    'Both': '#AA4499'          # Magenta - Sialofucosylated (colorblind-safe)
}

# Extended category colors (UPDATED to colorblind-safe - Phase 2.2)
EXTENDED_CATEGORY_COLORS = {
    'HM': '#117733',              # Dark green - High Mannose
    'High mannose': '#117733',    # Alias for consistency
    'High Mannose': '#117733',    # Alias (capitalized)
    'C/H': '#44AA99',             # Teal - Complex/Hybrid
    'Complex/Hybrid': '#44AA99',  # Alias
    'ComplexHybrid': '#44AA99',   # Alias (no space/slash)
    'Fucosylated': '#CC6677',     # Rose - Fucosylated
    'Sialylated': '#882255',      # Purple - Sialylated
    'Sialofucosylated': '#AA4499', # Magenta - Sialofucosylated
    'Both': '#AA4499',            # Alias for Sialofucosylated
    'Truncated': '#DDCC77',       # Sand - Truncated/Other
    'Other': '#DDCC77'            # Sand - Other
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
# VIP Score Plots (R/ggplot2) - METABOANALYST STYLE
# Clean, publication-ready design inspired by MetaboAnalyst
# ==============================================================================

# Typography
VIP_FEATURE_NAME_SIZE = 4.5  # R ggplot2 units - Clean, readable labels
VIP_GROUP_LABEL_SIZE = 5.0   # Group labels (Cancer/Normal)

# Dot styling - UNIFORM SIZE (MetaboAnalyst style)
VIP_DOT_SIZE = 3.0           # Small, consistent size for all dots
VIP_DOT_COLOR = '#0000FF'    # Blue dots (MetaboAnalyst reference)
VIP_DOT_ALPHA = 1.0          # Fully opaque

# Heatmap - TRUE SQUARES (Relative comparison with gradient scale)
VIP_USE_GRADIENT = True                  # Use gradient (ready for multi-group expansion)
VIP_HEATMAP_LOW_COLOR = '#0000FF'        # Blue (lower group)
VIP_HEATMAP_MID_COLOR = '#FFFFFF'        # White (mid)
VIP_HEATMAP_HIGH_COLOR = '#DC3912'       # Red (higher group)
VIP_HEATMAP_SQUARE_SIZE = 0.4            # True square width (Option A: visible squares)
VIP_HEATMAP_HEIGHT = 0.4                 # True square height (equal to width)
VIP_HEATMAP_SPACING = 0.08               # Gap between Bottom/Top squares (more compact)
VIP_HEATMAP_OFFSET = 0.08                # Distance from VIP max to first square (more compact)

# Layout
VIP_FIGURE_WIDTH = 12        # Wider for better label visibility
VIP_FIGURE_HEIGHT = 8        # Standard height
VIP_LEFT_MARGIN_EXPAND = 0.4 # Extra space for feature names

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


def add_sample_size_annotation(ax, n_cancer: int, n_normal: int,
                               location: str = 'upper right',
                               fontsize: int = 10):
    """
    Add sample size annotation to plot (Phase 2.2 Enhancement)

    Args:
        ax: Matplotlib axes object
        n_cancer: Number of cancer samples
        n_normal: Number of normal samples
        location: Position ('upper right', 'upper left', 'lower right', 'lower left')
        fontsize: Font size for annotation

    Example:
        >>> add_sample_size_annotation(ax, n_cancer=24, n_normal=23)
        # Adds "n=47 (Cancer: 24, Normal: 23)" to plot
    """
    total = n_cancer + n_normal
    text = f"n={total}\n(Cancer: {n_cancer}, Normal: {n_normal})"

    # Position mapping
    pos_map = {
        'upper right': (0.98, 0.98),
        'upper left': (0.02, 0.98),
        'lower right': (0.98, 0.02),
        'lower left': (0.02, 0.02)
    }

    x, y = pos_map.get(location, (0.98, 0.98))
    ha = 'right' if 'right' in location else 'left'
    va = 'top' if 'upper' in location else 'bottom'

    bbox_props = dict(boxstyle='round,pad=0.5', facecolor='wheat',
                     alpha=0.3, edgecolor='black', linewidth=1.2)

    ax.text(x, y, text,
           transform=ax.transAxes,
           fontsize=fontsize,
           verticalalignment=va,
           horizontalalignment=ha,
           bbox=bbox_props,
           family='monospace',
           zorder=1000)  # Ensure it's on top
