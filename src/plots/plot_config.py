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
# PREMIUM DESIGN SYSTEM INTEGRATION
# ==============================================================================
try:
    from .design_system import (
        TypographySystem, ColorSystem, VisualEffects,
        LayoutSystem, SmartAnnotation, PanelSystem
    )
    DESIGN_SYSTEM_AVAILABLE = True
except ImportError:
    DESIGN_SYSTEM_AVAILABLE = False

# ==============================================================================
# Font Sizes - ENHANCED TYPOGRAPHY (Publication Quality)
# Now with premium font stack and enhanced hierarchy
# ==============================================================================

# Enhanced font sizes using typography system
TITLE_SIZE = 18  # Enhanced from 16 (bolder presence)
AXIS_LABEL_SIZE = 15  # Enhanced from 14 (improved readability)
TICK_LABEL_SIZE = 13  # Enhanced from 12 (clearer at print resolution)
LEGEND_SIZE = 12  # Maintained for balance
LEGEND_TITLE_SIZE = 13
ANNOTATION_SIZE = 11  # Maintained

# Font Weights (stronger hierarchy)
TITLE_WEIGHT = 'bold'
AXIS_LABEL_WEIGHT = 'medium'  # Enhanced from 'normal' (better presence)
LEGEND_TITLE_WEIGHT = 'bold'

# Premium font preferences (fallback chain - macOS optimized)
# SF Pro Display is the macOS system font, excellent alternative to Inter
FONT_FAMILY_DISPLAY = 'SF Pro Display'  # Titles (macOS native, similar to Inter)
FONT_FAMILY_TEXT = 'SF Pro Text'        # Body text (macOS native, similar to Roboto)
FONT_FAMILY_MONO = 'SF Mono'            # Data/code (macOS native monospace)

# ==============================================================================
# Colors - Group Comparison (Cancer vs Normal)
# COLORBLIND-SAFE PALETTE (Phase 2.2 Enhancement)
# Based on ColorBrewer "Set1" - Verified accessible
# ==============================================================================

# Group colors - PREMIUM ENHANCED (Cancer vs Normal)
# Warmer, bolder colors with better emotional resonance
COLOR_CANCER = '#E63946'  # Cinnabar Red - Enhanced (warmer, bolder than #E41A1C)
COLOR_CANCER_LIGHT = '#FF758F'  # Light variant for gradients
COLOR_CANCER_DARK = '#B8283A'   # Dark variant for emphasis

COLOR_NORMAL = '#457B9D'  # Sapphire Blue - Enhanced (cooler, more stable than #377EB8)
COLOR_NORMAL_LIGHT = '#7DA3C6'  # Light variant for gradients
COLOR_NORMAL_DARK = '#2D5670'   # Dark variant for emphasis

GROUP_PALETTE = {
    'Cancer': COLOR_CANCER,
    'Normal': COLOR_NORMAL
}

# Extended palette with gradients
GROUP_PALETTE_PREMIUM = {
    'Cancer': COLOR_CANCER,
    'Cancer_Light': COLOR_CANCER_LIGHT,
    'Cancer_Dark': COLOR_CANCER_DARK,
    'Normal': COLOR_NORMAL,
    'Normal_Light': COLOR_NORMAL_LIGHT,
    'Normal_Dark': COLOR_NORMAL_DARK
}

# NOTE: Group colors are colorblind-safe:
# - Deuteranopia contrast: 2.56 (good for red-green colorblindness)
# - Shape markers used for additional distinction

# ==============================================================================
# Colors - Glycan Types (SEMANTIC COLOR SCHEME)
# Colors reflect biological/functional meaning of glycan types
# ==============================================================================

# SEMANTIC GLYCAN COLORS - Biologically meaningful color assignments
# These colors are chosen to match conventional glycan type representations:
# - Green: High Mannose (early/simple glycosylation)
# - Blue: Complex/Hybrid (mature/complex structures)
# - Pink: Sialylated (charged, acidic modifications)
# - Orange: Sialofucosylated (dual modifications)
# - Red: Fucosylated (core modifications)

GLYCAN_COLORS = {
    'HM': '#10B981',    # Emerald (Material 3.0) - Enhanced green, more vibrant
    'F': '#F59E0B',     # Amber (Material 3.0) - Warm energy color
    'S': '#8B5CF6',     # Violet (Material 3.0) - Transformation/charge
    'SF': '#EC4899',    # Pink (Material 3.0) - Complexity indicator
    'C/H': '#3B82F6'    # Blue (Material 3.0) - Stability/maturity
}

# PREMIUM: Material Design 3.0 palette (better color harmony + accessibility)
GLYCAN_COLORS_PREMIUM = {
    'HM': '#10B981',    # Emerald - High Mannose (growth, foundation)
    'F': '#F59E0B',     # Amber - Fucosylated (energy, core modification)
    'S': '#8B5CF6',     # Violet - Sialylated (transformation, charged)
    'SF': '#EC4899',    # Pink - Sialofucosylated (complexity, dual)
    'C/H': '#3B82F6'    # Blue - Complex/Hybrid (stability, mature)
}

# Color meaning rationale:
# - HM (Green): Simple, foundational glycan structures
# - F (Red): Core fucosylation, high visibility for important modification
# - S (Pink): Sialic acid (negatively charged), distinctive from others
# - SF (Orange): Combination color (red + yellow tones)
# - C/H (Blue): Complex structures, stable/mature

# Legacy glycan type colors - UPDATED TO MATERIAL DESIGN 3.0
LEGACY_GLYCAN_COLORS = {
    'Non': '#95A5A6',           # Gray - Non-modified (neutral, unchanged)
    'Sialylated': '#8B5CF6',    # Violet (Material 3.0) - Sialylated (charged, transformation)
    'Fucosylated': '#F59E0B',   # Amber (Material 3.0) - Fucosylated (energy, core modification)
    'Both': '#EC4899'           # Pink (Material 3.0) - Sialofucosylated (complexity, dual)
}

# Extended category colors - USER-REQUESTED SCHEME (aligned with constants.py)
# Biological meaning: HM=Green, C/H=Blue, F=Red, S=Pink, SF=Orange
EXTENDED_CATEGORY_COLORS = {
    'HM': '#00CC00',              # Green - High Mannose (early/simple glycosylation)
    'High mannose': '#00CC00',    # Alias for consistency
    'High Mannose': '#00CC00',    # Alias (capitalized)
    'C/H': '#0000FF',             # Blue - Complex/Hybrid (mature/complex structures)
    'Complex/Hybrid': '#0000FF',  # Alias
    'ComplexHybrid': '#0000FF',   # Alias (no space/slash)
    'F': '#FF0000',               # Red - Fucosylated (core modifications)
    'Fucosylated': '#FF0000',     # Alias
    'S': '#FF69B4',               # Pink - Sialylated (charged, acidic modifications)
    'Sialylated': '#FF69B4',      # Alias
    'SF': '#FFA500',              # Orange - Sialofucosylated (dual modifications)
    'Sialofucosylated': '#FFA500', # Alias
    'Both': '#FFA500',            # Alias for Sialofucosylated
    'Truncated': '#95A5A6',       # Gray - Truncated/Other (unchanged)
    'Other': '#95A5A6'            # Gray - Other (unchanged)
}

# ==============================================================================
# Shape Markers - COLORBLIND ACCESSIBILITY (Phase 3)
# Matplotlib marker styles for distinguishing categories without color
# ==============================================================================

# Group markers - Cancer vs Normal (colorblind-safe distinction)
GROUP_MARKERS = {
    'Cancer': 'o',      # Circle - standard marker
    'Normal': 's'       # Square - easily distinguished from circle
}

# Glycan type markers - Distinct shapes for each category
# Chosen for maximum visual distinction and clarity
GLYCAN_TYPE_MARKERS = {
    'HM': 'o',          # Circle - High Mannose (simple structure → simple shape)
    'F': '^',           # Triangle up - Fucosylated (pointing up → core modification)
    'S': 'v',           # Triangle down - Sialylated (pointing down → terminal modification)
    'SF': 'D',          # Diamond - Sialofucosylated (combination → combined shape)
    'C/H': 's'          # Square - Complex/Hybrid (structured → geometric)
}

# Legacy glycan type markers
LEGACY_GLYCAN_MARKERS = {
    'Non': 'o',         # Circle - Non-modified
    'Sialylated': 'v',  # Triangle down - Sialylated
    'Fucosylated': '^',  # Triangle up - Fucosylated
    'Both': 'D'         # Diamond - Both
}

# Regulation markers (for volcano plots)
REGULATION_MARKERS = {
    'Up in Cancer': '^',        # Triangle up - increased
    'Down in Cancer': 'v',      # Triangle down - decreased
    'Non-significant': 'o'      # Circle - no change
}

# All available matplotlib markers (for reference)
# Filled: 'o', 's', '^', 'v', 'D', 'p', 'h', '*', 'X', 'P', '<', '>'
# Unfilled: shapes + '' suffix (e.g., 'o', 'sf')
# Note: We use filled markers for better visibility at 300 DPI

# Marker size scaling (relative to base size)
MARKER_SIZE_SCALE = {
    'small': 0.7,
    'normal': 1.0,
    'large': 1.3
}

# ==============================================================================
# Plot Settings - PUBLICATION OPTIMIZED (Nature/Science Standards)
# ==============================================================================

# DPI TIERED SYSTEM - Optimizes file size while maintaining quality
# Based on plot complexity and typical usage:
# - Main figures: 200 DPI (sufficient for publication, ~60% file size reduction)
# - Supplementary: 150 DPI (~75% file size reduction)
# - Complex plots: 150 DPI (prevents multi-MB files)
DPI_MAIN = 200           # Main publication figures (standard plots)
DPI_SUPPLEMENTARY = 150  # Supplementary figures (detailed plots)
DPI_COMPLEX = 150        # Complex visualizations (heatmaps, volcano plots)
DPI = DPI_MAIN          # Default DPI for backward compatibility

# Grid settings (Prism style: minimal, subtle)
GRID_ALPHA = 0.15  # Much lighter (was 0.3)
GRID_LINESTYLE = '-'  # Solid, not dashed (Prism style)
GRID_LINEWIDTH = 0.5
GRID_COLOR = '#E0E0E0'  # Light gray, not pure gray

# PNG Compression Settings (File Size Optimization)
# Matplotlib supports automatic PNG optimization
SAVE_KWARGS = {
    'bbox_inches': 'tight',    # Tight bounding box (removes whitespace)
    'pad_inches': 0.1,         # Minimal padding
    'facecolor': 'white',      # White background (no transparency)
    'edgecolor': 'none'        # No edge color
}

# Legend settings (Prism style: clean, no shadow)
LEGEND_FRAMEON = True
LEGEND_FANCYBOX = False  # Square corners (Prism style)
LEGEND_SHADOW = False     # No shadow (cleaner)
LEGEND_FRAMEALPHA = 1.0   # Solid background (was 0.9)
LEGEND_EDGECOLOR = '#000000'  # Black border
LEGEND_FRAMEWIDTH = 1.0

# ==============================================================================
# Plot-Specific Settings - OPTIMIZED FOR PUBLICATION
# Nature/Science standards: Compact, high information density
# ==============================================================================

# Boxplot (Publication-optimized: compact, clear)
BOXPLOT_WIDTH = 0.65          # Optimized box width
BOXPLOT_FIGSIZE = (8, 5)      # More compact (was 10x6) - reduces file size
BOXPLOT_EXTENDED_FIGSIZE = (10, 5)  # More compact (was 12x6)
BOXPLOT_LINEWIDTH = 1.2       # Optimized line width (was 1.5)
BOXPLOT_FLIERSIZE = 4         # Smaller outlier markers (was 6) - reduces complexity
BOXPLOT_DPI = DPI_MAIN        # Main figure quality (200 DPI)

# Histogram (Optimized for readability)
HISTOGRAM_FIGSIZE = (12, 7)   # Compact (was 16x10) - significant size reduction
HISTOGRAM_X_ROTATION = 90
HISTOGRAM_X_HA = 'right'
HISTOGRAM_BAR_EDGEWIDTH = 0.8  # Slightly thinner (was 1.0)
HISTOGRAM_DPI = DPI_MAIN      # Main figure quality (200 DPI)

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
VIP_LEFT_MARGIN_EXPAND = 0.4  # Extra space for feature names

# VIP Plot Enhancements (Publication Quality)
VIP_SIGNIFICANCE_THRESHOLD = 1.0            # Standard cutoff for VIP significance
VIP_SHOW_THRESHOLD_LINE = True              # Add vertical line at VIP=1.0
VIP_THRESHOLD_LINE_COLOR = '#808080'        # Gray
VIP_THRESHOLD_LINE_TYPE = 'dashed'          # Dashed line
VIP_THRESHOLD_LINE_WIDTH = 0.8              # Line width
VIP_HEATMAP_SQUARE_SIZE_ENHANCED = 8        # Larger squares for better visibility (ggplot2 units)
VIP_SHOW_FOLD_CHANGE = True                 # Add fold change column
VIP_SHOW_SIGNIFICANCE = True                # Add significance stars column
VIP_FC_TEXT_SIZE = 3.5                      # Fold change text size
VIP_SIG_TEXT_SIZE = 4.0                     # Significance star size

# ==============================================================================
# PCA Plot - PUBLICATION OPTIMIZED
# ==============================================================================
PCA_FIGSIZE = (8, 6)         # Compact (was 10x8) - reduces file size
PCA_LABEL_FONTSIZE = 9       # Balanced (was 10)
PCA_LABEL_BBOX_PAD = 0.3     # Tighter padding (was 0.4)
PCA_POINT_SIZE = 120         # Optimized (was 150)
PCA_POINT_LINEWIDTH = 1.5    # Lighter edges (was 2.0)
PCA_POINT_ALPHA = 0.8
PCA_DPI = DPI_MAIN          # Main figure quality (200 DPI)

# ==============================================================================
# Volcano Plot - OPTIMIZED FOR CLARITY AND FILE SIZE
# ==============================================================================
VOLCANO_FIGSIZE = (10, 8)    # Compact (was 14x11) - major size reduction
VOLCANO_POINT_SIZE = 80      # Optimized (was 120) - reduces complexity
VOLCANO_POINT_ALPHA = 0.7    # Slightly transparent (was 0.8)
VOLCANO_THRESHOLD_LINEWIDTH = 2.0  # Optimized (was 3.0)
VOLCANO_THRESHOLD_ALPHA = 0.5
VOLCANO_POINT_EDGEWIDTH = 0.8  # Lighter edges (was 1.2)
VOLCANO_LABEL_FONTSIZE = 11    # Balanced (was 14)
VOLCANO_LABEL_WEIGHT = 'bold'
VOLCANO_LABEL_PADDING = 0.4    # Tighter (was 0.6)
VOLCANO_LABEL_LINEWIDTH = 1.5  # Lighter (was 2.5)
VOLCANO_MAX_LABELS = 3
VOLCANO_DPI = DPI_COMPLEX   # Complex plot (150 DPI) - prevents >1MB files

# Heatmap (File size critical - often >500K)
HEATMAP_FIGSIZE = (10, 7)   # Compact (was 12x9)
HEATMAP_CBAR_LABEL_SIZE = 10  # Slightly smaller (was 11)
HEATMAP_DPI = DPI_COMPLEX   # Complex plot (150 DPI) - major size reduction

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


# ==============================================================================
# Colorblind Accessibility Helpers (Phase 3)
# ==============================================================================

def get_group_style(group: str) -> tuple:
    """
    Get color and marker for a group (Cancer/Normal)

    Phase 3 Enhancement: Returns both color and marker for dual encoding

    Args:
        group: Group name ('Cancer' or 'Normal')

    Returns:
        Tuple of (color, marker)

    Example:
        >>> color, marker = get_group_style('Cancer')
        >>> # color = '#E41A1C', marker = 'o'
    """
    color = GROUP_PALETTE.get(group, '#333333')
    marker = GROUP_MARKERS.get(group, 'o')
    return color, marker


def get_glycan_type_style(glycan_type: str) -> tuple:
    """
    Get color and marker for a glycan type

    Phase 3 Enhancement: Returns both color and marker for dual encoding

    Args:
        glycan_type: Glycan type ('HM', 'F', 'S', 'SF', 'C/H')

    Returns:
        Tuple of (color, marker)

    Example:
        >>> color, marker = get_glycan_type_style('HM')
        >>> # color = '#27AE60', marker = 'o'
    """
    color = GLYCAN_COLORS.get(glycan_type, '#333333')
    marker = GLYCAN_TYPE_MARKERS.get(glycan_type, 'o')
    return color, marker


def get_regulation_style(regulation: str) -> tuple:
    """
    Get color and marker for regulation status (volcano plots)

    Phase 3 Enhancement: Returns both color and marker for dual encoding

    Args:
        regulation: Regulation status ('Up in Cancer', 'Down in Cancer', 'Non-significant')

    Returns:
        Tuple of (color, marker)

    Example:
        >>> color, marker = get_regulation_style('Up in Cancer')
        >>> # color = '#E41A1C', marker = '^'
    """
    color_map = {
        'Up in Cancer': GROUP_PALETTE['Cancer'],
        'Down in Cancer': GROUP_PALETTE['Normal'],
        'Non-significant': '#95A5A6'
    }

    color = color_map.get(regulation, '#333333')
    marker = REGULATION_MARKERS.get(regulation, 'o')
    return color, marker


# ==============================================================================
# Optimized Saving Function - PUBLICATION QUALITY + FILE SIZE OPTIMIZATION
# ==============================================================================

def save_publication_figure(fig, filepath, dpi=None, **kwargs):
    """
    Save figure with publication-quality settings and file size optimization

    Automatically applies:
    - Appropriate DPI based on plot type
    - Tight bounding box (removes whitespace)
    - White background (no transparency)
    - Optimal compression

    Args:
        fig: Matplotlib figure object
        filepath: Output file path (Path object or string)
        dpi: DPI override (if None, uses DPI_MAIN=200)
        **kwargs: Additional arguments passed to savefig()

    Example:
        >>> fig, ax = plt.subplots()
        >>> ax.plot([1, 2, 3])
        >>> save_publication_figure(fig, 'output.png')
        # Saves with 200 DPI, tight layout, optimized

    File Size Reduction:
        - 200 DPI vs 300 DPI: ~60% reduction
        - 150 DPI vs 300 DPI: ~75% reduction
        - bbox_inches='tight': Removes whitespace padding
    """
    if dpi is None:
        dpi = DPI_MAIN

    # Merge default save kwargs with user-provided kwargs
    save_params = SAVE_KWARGS.copy()
    save_params.update(kwargs)
    save_params['dpi'] = dpi

    # Save with optimized settings
    fig.savefig(filepath, **save_params)


def get_plot_dpi(plot_type: str) -> int:
    """
    Get appropriate DPI for plot type (tiered system)

    Plot Types:
        - 'main': Standard plots (boxplot, PCA, etc.) - 200 DPI
        - 'complex': Complex plots (heatmap, volcano) - 150 DPI
        - 'supplementary': Supplementary figures - 150 DPI

    Args:
        plot_type: Type of plot ('main', 'complex', 'supplementary')

    Returns:
        Appropriate DPI value

    Example:
        >>> dpi = get_plot_dpi('complex')
        >>> # Returns 150 (for heatmaps, volcano plots)
    """
    dpi_map = {
        'main': DPI_MAIN,              # 200 DPI
        'complex': DPI_COMPLEX,        # 150 DPI
        'supplementary': DPI_SUPPLEMENTARY  # 150 DPI
    }
    return dpi_map.get(plot_type, DPI_MAIN)
"""
Enhanced Plot Configuration for pGlyco Auto Combine
Publication-quality styling with advanced visual effects

NEW ENHANCEMENTS:
- Perceptually uniform colormaps (viridis, RdBu_r)
- Enhanced typography (Arial/Helvetica)
- Subtle shadows and rounded corners
- Gradient fills and professional polish
"""

import numpy as np
import matplotlib as mpl
from matplotlib.patches import FancyBboxPatch

# ==============================================================================
# Enhanced Colormaps - PERCEPTUALLY UNIFORM (Publication Quality)
# ==============================================================================

# Sequential colormaps for intensity data (perceptually uniform)
CMAP_SEQUENTIAL = 'viridis'  # Best for general intensity data
CMAP_SEQUENTIAL_ALT = 'plasma'  # Alternative for variety
CMAP_SEQUENTIAL_BLUE = 'Blues'  # For single-color sequential

# Diverging colormaps for fold change / differential data
CMAP_DIVERGING = 'RdBu_r'  # Red-Blue reversed (Red=high, Blue=low)
CMAP_DIVERGING_ALT = 'coolwarm'  # Alternative red-blue
CMAP_DIVERGING_GREEN = 'PiYG'  # Purple-Yellow-Green

# Colorblind-safe sequential
CMAP_COLORBLIND = 'cividis'  # Perceptually uniform + colorblind-safe

# Heatmap colormap defaults
HEATMAP_CMAP_INTENSITY = 'viridis'  # For intensity heatmaps
HEATMAP_CMAP_FOLDCHANGE = 'RdBu_r'  # For fold change heatmaps
HEATMAP_CMAP_CORRELATION = 'RdBu_r'  # For correlation matrices

# ==============================================================================
# Enhanced Visual Effects - PUBLICATION POLISH
# ==============================================================================

# Font Family Preferences (Nature/Science journal standards)
FONT_FAMILY = 'sans-serif'
FONT_NAME = 'Arial'  # Preferred: Arial, Helvetica, DejaVu Sans

# Shadow effects (subtle, professional)
SHADOW_ENABLED = True
SHADOW_ALPHA = 0.15  # Very subtle shadow
SHADOW_OFFSET = (2, -2)  # Offset in points (x, y)
SHADOW_COLOR = '#000000'  # Black shadow

# Rounded corners for boxes and annotations
BBOX_CORNER_RADIUS = 0.3  # Rounded corner radius (pad units)
BBOX_EDGE_WIDTH = 1.2  # Border width for annotation boxes
BBOX_ALPHA = 0.9  # Slight transparency for overlay boxes

# Gradient fills for visual appeal
GRADIENT_ENABLED = True  # Enable gradient fills where appropriate
GRADIENT_ALPHA_START = 0.7  # Starting alpha for gradient
GRADIENT_ALPHA_END = 0.3  # Ending alpha for gradient

# Enhanced marker properties
MARKER_EDGE_ALPHA = 0.8  # Marker edge transparency
MARKER_GLOW_ENABLED = False  # Enable glow effect (use sparingly)

# ==============================================================================
# Enhanced Styling Helpers - PUBLICATION POLISH
# ==============================================================================

def create_fancy_bbox(facecolor='white', edgecolor='black', alpha=0.9, linewidth=1.2):
    """
    Create a fancy bounding box for annotations with rounded corners

    Args:
        facecolor: Background color
        edgecolor: Border color
        alpha: Transparency
        linewidth: Border width

    Returns:
        Dictionary of bbox properties for matplotlib text
    """
    bbox_props = {
        'boxstyle': f'round,pad={BBOX_CORNER_RADIUS}',
        'facecolor': facecolor,
        'edgecolor': edgecolor,
        'alpha': alpha,
        'linewidth': linewidth
    }
    return bbox_props


def add_shadow_effect(ax, artist, offset=(2, -2), alpha=0.15):
    """
    Add subtle shadow effect to a matplotlib artist

    Args:
        ax: Matplotlib axes
        artist: Artist object (patch, line, etc.)
        offset: Shadow offset in points (x, y)
        alpha: Shadow transparency
    """
    from matplotlib.patheffects import withSimplePatchShadow

    if SHADOW_ENABLED:
        shadow = withSimplePatchShadow(
            offset=offset,
            shadow_rgbFace='black',
            alpha=alpha
        )
        artist.set_path_effects([shadow])


def create_gradient_fill(ax, x, y1, y2, color, alpha_start=0.7, alpha_end=0.3):
    """
    Create gradient fill between two curves

    Args:
        ax: Matplotlib axes
        x: X coordinates
        y1, y2: Y coordinates for fill boundaries
        color: Base color
        alpha_start: Starting alpha (bottom)
        alpha_end: Ending alpha (top)

    Returns:
        Collection of gradient patches
    """
    if not GRADIENT_ENABLED:
        # Simple fill without gradient
        ax.fill_between(x, y1, y2, color=color, alpha=alpha_start)
        return

    # Create gradient using multiple alpha levels
    n_gradients = 20
    alphas = np.linspace(alpha_start, alpha_end, n_gradients)

    for i in range(n_gradients - 1):
        y_bottom = y1 + (y2 - y1) * i / n_gradients
        y_top = y1 + (y2 - y1) * (i + 1) / n_gradients
        ax.fill_between(x, y_bottom, y_top, color=color, alpha=alphas[i],
                       linewidth=0, zorder=0)


def apply_publication_theme(fig, ax_list=None):
    """
    Apply comprehensive publication theme to figure and axes

    Args:
        fig: Matplotlib figure object
        ax_list: List of axes (if None, applies to all axes in figure)
    """
    # Set figure background
    fig.patch.set_facecolor('white')
    fig.patch.set_alpha(1.0)

    # Get all axes if not provided
    if ax_list is None:
        ax_list = fig.get_axes()

    # Apply theme to each axis
    for ax in ax_list:
        # Set font properties
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                     ax.get_xticklabels() + ax.get_yticklabels()):
            try:
                item.set_fontfamily(FONT_FAMILY)
            except:
                pass

        # Enhanced spine styling
        for spine in ax.spines.values():
            spine.set_linewidth(1.2)
            spine.set_edgecolor('#333333')

        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)


def enhance_statistical_bracket(ax, x1, x2, y, text, color='black', fontsize=14):
    """
    Add enhanced statistical significance bracket with rounded ends

    Args:
        ax: Matplotlib axes
        x1, x2: X positions for bracket ends
        y: Y position for bracket
        text: Significance text (e.g., '***', 'p<0.001')
        color: Bracket and text color
        fontsize: Font size for text
    """
    # Calculate bracket height
    y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
    bracket_height = 0.02 * y_range

    # Draw horizontal line with rounded caps
    ax.plot([x1, x2], [y, y], color=color, linewidth=2.5,
            solid_capstyle='round', zorder=10)

    # Vertical caps
    ax.plot([x1, x1], [y, y - bracket_height], color=color, linewidth=2.5,
            solid_capstyle='round', zorder=10)
    ax.plot([x2, x2], [y, y - bracket_height], color=color, linewidth=2.5,
            solid_capstyle='round', zorder=10)

    # Text annotation with fancy box
    bbox_props = create_fancy_bbox(facecolor='white', edgecolor=color, alpha=0.95)
    ax.text((x1 + x2) / 2, y, text, ha='center', va='bottom',
            fontsize=fontsize, fontweight='bold', color=color,
            bbox=bbox_props, zorder=11, fontfamily=FONT_FAMILY)


def enhance_heatmap_colorbar(cbar, label='', fontsize=12):
    """
    Enhance colorbar styling for heatmaps

    Args:
        cbar: Matplotlib colorbar object
        label: Colorbar label text
        fontsize: Font size for label
    """
    # Set label with enhanced styling (handle both Colorbar and Axes objects)
    if label:
        cbar.set_label(label)
        # Handle different colorbar types
        if hasattr(cbar, 'ax'):
            # Standard matplotlib Colorbar object
            label_obj = cbar.ax.yaxis.label
            axes_obj = cbar.ax
        else:
            # Seaborn clustermap - cbar IS the axes
            label_obj = cbar.yaxis.label
            axes_obj = cbar

        label_obj.set_fontsize(fontsize)
        label_obj.set_fontweight('bold')
        label_obj.set_fontfamily(FONT_FAMILY)

    # Enhanced tick styling
    if hasattr(cbar, 'ax'):
        axes_obj = cbar.ax
    else:
        axes_obj = cbar
    axes_obj.tick_params(labelsize=fontsize - 2, width=1.2, length=4)

    # Colorbar outline (only for matplotlib Colorbar objects)
    if hasattr(cbar, 'outline'):
        cbar.outline.set_linewidth(1.2)
        cbar.outline.set_edgecolor('#333333')
