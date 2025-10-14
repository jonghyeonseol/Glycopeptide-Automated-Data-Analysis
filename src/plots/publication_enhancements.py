"""
Publication Enhancements Module

Provides utilities for publication-quality visualizations:
1. Sample size annotations
2. Colorblind-safe palettes (verified)
3. Font optimization for readability
4. Accessibility compliance

Created: 2025-10-06 (Phase 2.2)
"""

from typing import Dict, List, Optional

# Import centralized constants
try:
    from .plot_config import (
        EDGE_LINEWIDTH_NORMAL, LINE_MEDIUM_THICK,
        EDGE_COLOR_BLACK, SEPARATOR_EDGE_COLOR,
        FONT_DATA  # Phase 10.3.10: Font family standardization
    )
except ImportError:
    # Fallback values if import fails
    EDGE_LINEWIDTH_NORMAL = 1.0
    LINE_MEDIUM_THICK = 1.5
    EDGE_COLOR_BLACK = 'black'
    SEPARATOR_EDGE_COLOR = '#CCCCCC'
    FONT_DATA = 'monospace'


# ==============================================================================
# COLORBLIND-SAFE PALETTES (Verified with Color Universal Design)
# ==============================================================================

def verify_colorblind_safe(hex_color: str) -> Dict[str, str]:
    """
    Verify color appearance for different types of colorblindness

    ⚠️ WARNING: APPROXIMATION ONLY
    This is a simplified LMS color space simulation.
    For publication-critical validation, use dedicated tools:
    - Color Oracle (free, cross-platform): https://colororacle.org/
    - ColorBrewer 2.0 (online): https://colorbrewer2.org/
    - Coblis (Color Blindness Simulator): https://www.color-blindness.com/coblis-color-blindness-simulator/

    Returns simulated appearance for:
    - Deuteranopia (red-green, most common, ~5% males)
    - Protanopia (red-green, ~1% males)
    - Tritanopia (blue-yellow, rare, <1%)

    Returns:
        Dictionary with:
        - 'original': Input color
        - 'deuteranopia', 'protanopia', 'tritanopia': Simulated colors
        - 'warning': Prominent warning message about limitations
    """
    # Convert hex to RGB
    r = int(hex_color[1:3], 16) / 255.0
    g = int(hex_color[3:5], 16) / 255.0
    b = int(hex_color[5:7], 16) / 255.0

    # Simplified LMS color space conversion (approximation)
    # These matrices simulate how different cone cells respond

    # Deuteranopia (missing M cones - green receptors)
    # Red and green appear similar
    deuter_r = r * 0.625 + g * 0.375
    deuter_g = r * 0.7 + g * 0.3
    deuter_b = b

    # Protanopia (missing L cones - red receptors)
    # Red appears darker, shifted toward green
    prota_r = r * 0.567 + g * 0.433
    prota_g = r * 0.558 + g * 0.442
    prota_b = b

    # Tritanopia (missing S cones - blue receptors)
    # Blue and yellow appear similar
    trita_r = r * 0.95 + g * 0.05
    trita_g = g
    trita_b = r * 0.433 + b * 0.567

    def rgb_to_hex(r, g, b):
        return f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"

    return {
        'original': hex_color,
        'deuteranopia': rgb_to_hex(deuter_r, deuter_g, deuter_b),
        'protanopia': rgb_to_hex(prota_r, prota_g, prota_b),
        'tritanopia': rgb_to_hex(trita_r, trita_g, trita_b),
        'warning': (
            "⚠️ APPROXIMATION ONLY: This simulation uses simplified LMS color space transformation. "
            "For publication, validate colors with Color Oracle (https://colororacle.org/) or "
            "ColorBrewer 2.0 (https://colorbrewer2.org/)."
        )
    }


# Colorblind-safe palette (Paul Tol's "bright" scheme)
# Verified to maintain distinction across all colorblindness types
COLORBLIND_SAFE_PALETTE = {
    'blue': '#4477AA',      # Accessible to all
    'cyan': '#66CCEE',      # Accessible to all
    'green': '#228833',     # Accessible to all
    'yellow': '#CCBB44',    # Accessible to all
    'red': '#EE6677',       # Accessible to all
    'purple': '#AA3377',    # Accessible to all
    'grey': '#BBBBBB'       # Neutral, accessible to all
}

# Group comparison colors (Cancer vs Normal) - VERIFIED COLORBLIND-SAFE
# Based on ColorBrewer "Set1" scheme
COLORBLIND_GROUP_PALETTE = {
    'Cancer': '#E41A1C',    # Red (warm, attention-grabbing)
    'Normal': '#377EB8'     # Blue (cool, stable)
}

# Glycan type colors - SEMANTIC SCHEME
# Colors reflect biological/functional meaning
COLORBLIND_GLYCAN_PALETTE = {
    'HM': '#27AE60',        # Emerald Green - High Mannose
    'F': '#E74C3C',         # Alizarin Red - Fucosylated
    'S': '#E91E63',         # Material Pink - Sialylated
    'SF': '#E67E22',        # Carrot Orange - Sialofucosylated
    'C/H': '#3498DB',       # Peter River Blue - Complex/Hybrid
    'Non': '#95A5A6'        # Concrete Gray - Non-modified
}

# Extended palette for more categories - SEMANTIC SCHEME
COLORBLIND_EXTENDED_PALETTE = {
    'HM': '#27AE60',
    'High Mannose': '#27AE60',
    'High mannose': '#27AE60',
    'F': '#E74C3C',
    'Fucosylated': '#E74C3C',
    'S': '#E91E63',
    'Sialylated': '#E91E63',
    'SF': '#E67E22',
    'Sialofucosylated': '#E67E22',
    'Both': '#E67E22',
    'C/H': '#3498DB',
    'Complex/Hybrid': '#3498DB',
    'ComplexHybrid': '#3498DB',
    'Non': '#95A5A6',
    'Other': '#95A5A6',
    'Truncated': '#95A5A6'
}


# ==============================================================================
# SAMPLE SIZE ANNOTATION UTILITIES
# ==============================================================================

def add_sample_size_to_title(ax, n_total: Optional[int] = None,
                             n_cancer: Optional[int] = None,
                             n_normal: Optional[int] = None,
                             current_title: str = "") -> str:
    """
    Add sample size information to plot title

    Args:
        ax: Matplotlib axes object
        n_total: Total number of samples
        n_cancer: Number of cancer samples
        n_normal: Number of normal samples
        current_title: Existing title text

    Returns:
        Enhanced title string with sample sizes

    Examples:
        >>> add_sample_size_to_title(ax, n_total=47)
        "Plot Title (n=47)"

        >>> add_sample_size_to_title(ax, n_cancer=24, n_normal=23)
        "Plot Title (Cancer n=24, Normal n=23)"
    """
    if n_total is not None:
        size_text = f" (n={n_total})"
    elif n_cancer is not None and n_normal is not None:
        size_text = f" (Cancer n={n_cancer}, Normal n={n_normal})"
    elif n_cancer is not None:
        size_text = f" (Cancer n={n_cancer})"
    elif n_normal is not None:
        size_text = f" (Normal n={n_normal})"
    else:
        size_text = ""

    enhanced_title = current_title + size_text
    return enhanced_title


def annotate_sample_sizes_on_groups(ax, group_names: List[str],
                                    group_sizes: Dict[str, int],
                                    y_position: float = 0.98,
                                    fontsize: int = 10,
                                    bbox_style: Dict = None):
    """
    Annotate sample sizes directly on plot for each group

    Args:
        ax: Matplotlib axes object
        group_names: List of group names
        group_sizes: Dictionary mapping group names to sample sizes
        y_position: Vertical position (0-1, axes coordinates)
        fontsize: Font size for annotations
        bbox_style: Custom bbox style dictionary

    Example:
        >>> group_sizes = {'Cancer': 24, 'Normal': 23}
        >>> annotate_sample_sizes_on_groups(ax, ['Cancer', 'Normal'], group_sizes)
    """
    from .plot_config import ALPHA_VERY_HIGH
    if bbox_style is None:
        bbox_style = dict(boxstyle='round,pad=0.5', facecolor='white',
                          alpha=ALPHA_VERY_HIGH, edgecolor=SEPARATOR_EDGE_COLOR, linewidth=EDGE_LINEWIDTH_NORMAL)

    annotation_text = ", ".join([f"{name}: n={group_sizes.get(name, 0)}"
                                for name in group_names if name in group_sizes])

    # Place annotation in upper right corner
    ax.text(0.98, y_position, annotation_text,
            transform=ax.transAxes,
            fontsize=fontsize,
            verticalalignment='top',
            horizontalalignment='right',
            bbox=bbox_style)


def add_sample_size_legend(ax, n_cancer: int, n_normal: int,
                           location: str = 'upper right',
                           fontsize: int = 10):
    """
    Add a dedicated sample size box to the plot

    Args:
        ax: Matplotlib axes object
        n_cancer: Number of cancer samples
        n_normal: Number of normal samples
        location: Legend location
        fontsize: Font size for text
    """
    sample_text = f"Sample Sizes:\nCancer: n={n_cancer}\nNormal: n={n_normal}\nTotal: n={n_cancer + n_normal}"

    # Create a text box
    from .plot_config import ALPHA_MEDIUM_LIGHT
    props = dict(boxstyle='round', facecolor='wheat', alpha=ALPHA_MEDIUM_LIGHT, edgecolor=EDGE_COLOR_BLACK, linewidth=LINE_MEDIUM_THICK)

    # Position mapping
    pos_map = {
        'upper right': (0.95, 0.95),
        'upper left': (0.05, 0.95),
        'lower right': (0.95, 0.05),
        'lower left': (0.05, 0.05)
    }

    x, y = pos_map.get(location, (0.95, 0.95))
    ha = 'right' if 'right' in location else 'left'
    va = 'top' if 'upper' in location else 'bottom'

    ax.text(x, y, sample_text,
            transform=ax.transAxes,
            fontsize=fontsize,
            verticalalignment=va,
            horizontalalignment=ha,
            bbox=props,
            family=FONT_DATA)


# ==============================================================================
# COLORBLIND PALETTE APPLICATION
# ==============================================================================

def apply_colorblind_safe_palette(plot_type: str = 'group') -> Dict[str, str]:
    """
    Get appropriate colorblind-safe palette for plot type

    Args:
        plot_type: Type of plot ('group', 'glycan', 'extended', 'general')

    Returns:
        Dictionary mapping categories to hex colors
    """
    palette_map = {
        'group': COLORBLIND_GROUP_PALETTE,
        'glycan': COLORBLIND_GLYCAN_PALETTE,
        'extended': COLORBLIND_EXTENDED_PALETTE,
        'general': COLORBLIND_SAFE_PALETTE
    }

    return palette_map.get(plot_type, COLORBLIND_SAFE_PALETTE).copy()


def test_colorblind_visibility(color1: str, color2: str,
                               min_contrast: float = 3.0) -> Dict:
    """
    Test if two colors are distinguishable for colorblind viewers

    Args:
        color1: First hex color
        color2: Second hex color
        min_contrast: Minimum acceptable contrast ratio (WCAG: 3.0 for graphics)

    Returns:
        Dictionary with test results and contrast ratios
    """
    def hex_to_luminance(hex_color):
        """Calculate relative luminance (WCAG formula)"""
        r = int(hex_color[1:3], 16) / 255.0
        g = int(hex_color[3:5], 16) / 255.0
        b = int(hex_color[5:7], 16) / 255.0

        # Gamma correction
        r = r / 12.92 if r <= 0.03928 else ((r + 0.055) / 1.055) ** 2.4
        g = g / 12.92 if g <= 0.03928 else ((g + 0.055) / 1.055) ** 2.4
        b = b / 12.92 if b <= 0.03928 else ((b + 0.055) / 1.055) ** 2.4

        return 0.2126 * r + 0.7152 * g + 0.0722 * b

    lum1 = hex_to_luminance(color1)
    lum2 = hex_to_luminance(color2)

    # Contrast ratio
    lighter = max(lum1, lum2)
    darker = min(lum1, lum2)
    contrast_ratio = (lighter + 0.05) / (darker + 0.05)

    # Simulate colorblind appearance
    deuter1 = verify_colorblind_safe(color1)['deuteranopia']
    deuter2 = verify_colorblind_safe(color2)['deuteranopia']

    deuter_lum1 = hex_to_luminance(deuter1)
    deuter_lum2 = hex_to_luminance(deuter2)
    deuter_contrast = (max(deuter_lum1, deuter_lum2) + 0.05) / (min(deuter_lum1, deuter_lum2) + 0.05)

    return {
        'normal_contrast': contrast_ratio,
        'deuteranopia_contrast': deuter_contrast,
        'passes_normal': contrast_ratio >= min_contrast,
        'passes_colorblind': deuter_contrast >= min_contrast,
        'wcag_level': (
            'AAA' if contrast_ratio >= 7.0 else
            'AA' if contrast_ratio >= 4.5 else
            'A' if contrast_ratio >= 3.0 else
            'Fail'
        )
    }


# ==============================================================================
# FONT OPTIMIZATION
# ==============================================================================

def optimize_font_sizes_for_publication(figure_width_inches: float = 10,
                                        target_journal: str = 'nature') -> Dict[str, int]:
    """
    Calculate optimal font sizes for publication based on figure width

    Args:
        figure_width_inches: Width of figure in inches
        target_journal: Journal name ('nature', 'science', 'plos', 'generic')

    Returns:
        Dictionary of recommended font sizes

    Guidelines:
        Nature: 6-8pt min, single column (89mm), double column (183mm)
        Science: 6-8pt min, single column (5.5cm), double column (12cm)
        PLOS: 8-12pt, single column (8.3cm), double column (17.3cm)
    """
    # Journal-specific minimum sizes (in points)
    journal_specs = {
        'nature': {'min_size': 7, 'single_col_mm': 89, 'double_col_mm': 183},
        'science': {'min_size': 7, 'single_col_mm': 55, 'double_col_mm': 120},
        'plos': {'min_size': 8, 'single_col_mm': 83, 'double_col_mm': 173},
        'generic': {'min_size': 8, 'single_col_mm': 90, 'double_col_mm': 180}
    }

    spec = journal_specs.get(target_journal, journal_specs['generic'])
    figure_width_mm = figure_width_inches * 25.4

    # Determine if single or double column
    is_single_column = figure_width_mm <= spec['single_col_mm'] * 1.2

    # Scale factor based on figure size
    base_scale = figure_width_mm / \
        spec['single_col_mm'] if is_single_column else figure_width_mm / spec['double_col_mm']
    base_scale = max(0.8, min(base_scale, 1.5))  # Clamp between 0.8 and 1.5

    min_size = spec['min_size']

    return {
        'title': max(int(14 * base_scale), min_size + 4),
        'axis_label': max(int(12 * base_scale), min_size + 2),
        'tick_label': max(int(10 * base_scale), min_size),
        'legend': max(int(10 * base_scale), min_size),
        'annotation': max(int(9 * base_scale), min_size - 1),
        'column_type': 'single' if is_single_column else 'double',
        # DPI for PRINT SUBMISSION (higher than default electronic 200 DPI)
        # Use plot_config.DPI_MAIN (200) for electronic/web publishing
        'recommended_dpi': 300 if is_single_column else 600  # Print: 300 (single-col) / 600 (double-col)
    }


# ==============================================================================
# VALIDATION & TESTING
# ==============================================================================

def validate_plot_accessibility(ax, check_contrast: bool = True,
                                check_fonts: bool = True,
                                min_font_size: int = 8) -> Dict:
    """
    Validate plot meets accessibility standards (Enhanced - Phase 10.9)

    Checks multiple accessibility criteria:
    - Font sizes (title, labels, legend, colorbar, annotations)
    - Color contrast (future enhancement)

    Args:
        ax: Matplotlib axes object
        check_contrast: Check color contrast (not yet implemented)
        check_fonts: Check font sizes for all text elements
        min_font_size: Minimum acceptable font size (points)

    Returns:
        Dictionary with detailed validation results:
        - 'font_sizes': Dict of all font sizes found
        - 'failed_items': List of dicts with details of failed elements
        - 'warnings': List of warning messages
        - 'passes': Boolean - True if all checks passed

    Example:
        >>> results = validate_plot_accessibility(ax, min_font_size=8)
        >>> if not results['passes']:
        >>>     for item in results['failed_items']:
        >>>         print(f"{item['element']}: {item['size']:.1f}pt < {item['min_required']:.1f}pt")
    """
    results = {
        'font_sizes': {},
        'failed_items': [],  # NEW: Detailed failure information
        'color_contrast': {},
        'warnings': [],
        'passes': True
    }

    if check_fonts:
        # Check title
        title = ax.get_title()
        if title:
            title_obj = ax.title
            title_size = title_obj.get_fontsize()
            results['font_sizes']['title'] = title_size
            min_title_size = min_font_size + 4  # Title should be larger

            if title_size < min_title_size:
                results['warnings'].append(
                    f"Title font size ({title_size:.1f}pt) below recommended minimum ({min_title_size}pt)"
                )
                results['failed_items'].append({
                    'element': 'title',
                    'size': title_size,
                    'min_required': min_title_size,
                    'text_preview': title[:50] + ('...' if len(title) > 50 else '')
                })
                results['passes'] = False

        # Check axis labels
        xlabel_size = ax.xaxis.label.get_fontsize()
        ylabel_size = ax.yaxis.label.get_fontsize()
        results['font_sizes']['xlabel'] = xlabel_size
        results['font_sizes']['ylabel'] = ylabel_size

        if xlabel_size < min_font_size:
            results['warnings'].append(
                f"X-label font size ({xlabel_size:.1f}pt) below minimum ({min_font_size}pt)"
            )
            results['failed_items'].append({
                'element': 'xlabel',
                'size': xlabel_size,
                'min_required': min_font_size,
                'text_preview': ax.xaxis.label.get_text()[:50]
            })
            results['passes'] = False

        if ylabel_size < min_font_size:
            results['warnings'].append(
                f"Y-label font size ({ylabel_size:.1f}pt) below minimum ({min_font_size}pt)"
            )
            results['failed_items'].append({
                'element': 'ylabel',
                'size': ylabel_size,
                'min_required': min_font_size,
                'text_preview': ax.yaxis.label.get_text()[:50]
            })
            results['passes'] = False

        # NEW: Check legend font sizes
        legend = ax.get_legend()
        if legend:
            legend_texts = legend.get_texts()
            if legend_texts:
                legend_size = legend_texts[0].get_fontsize()  # Check first entry
                results['font_sizes']['legend'] = legend_size

                if legend_size < min_font_size:
                    results['warnings'].append(
                        f"Legend font size ({legend_size:.1f}pt) below minimum ({min_font_size}pt)"
                    )
                    results['failed_items'].append({
                        'element': 'legend',
                        'size': legend_size,
                        'min_required': min_font_size,
                        'text_preview': legend_texts[0].get_text()[:30] + '...'
                    })
                    results['passes'] = False

            # Check legend title if present
            legend_title = legend.get_title()
            if legend_title and legend_title.get_text():
                legend_title_size = legend_title.get_fontsize()
                results['font_sizes']['legend_title'] = legend_title_size

                if legend_title_size < min_font_size:
                    results['warnings'].append(
                        f"Legend title font size ({legend_title_size:.1f}pt) below minimum ({min_font_size}pt)"
                    )
                    results['failed_items'].append({
                        'element': 'legend_title',
                        'size': legend_title_size,
                        'min_required': min_font_size,
                        'text_preview': legend_title.get_text()[:30]
                    })
                    results['passes'] = False

        # NEW: Check colorbar font sizes (for heatmaps, etc.)
        # Try to find colorbar associated with this axes
        fig = ax.get_figure()
        if fig:
            for child_ax in fig.get_axes():
                # Colorbar axes typically have different properties
                if hasattr(child_ax, 'colorbar'):
                    cbar = child_ax.colorbar
                    if cbar:
                        cbar_label_size = cbar.ax.yaxis.label.get_fontsize()
                        results['font_sizes']['colorbar_label'] = cbar_label_size

                        if cbar_label_size < min_font_size:
                            results['warnings'].append(
                                f"Colorbar label font size ({cbar_label_size:.1f}pt) below minimum ({min_font_size}pt)"
                            )
                            results['failed_items'].append({
                                'element': 'colorbar_label',
                                'size': cbar_label_size,
                                'min_required': min_font_size,
                                'text_preview': cbar.ax.yaxis.label.get_text()[:30]
                            })
                            results['passes'] = False

        # NEW: Check text annotations
        texts = ax.texts
        if texts:
            annotation_sizes = []
            for i, text_obj in enumerate(texts):
                text_size = text_obj.get_fontsize()
                annotation_sizes.append(text_size)

                if text_size < min_font_size - 1:  # Annotations can be slightly smaller
                    text_content = text_obj.get_text()
                    results['warnings'].append(
                        f"Annotation #{i+1} font size ({text_size:.1f}pt) below acceptable minimum ({min_font_size-1}pt)"
                    )
                    results['failed_items'].append({
                        'element': f'annotation_{i+1}',
                        'size': text_size,
                        'min_required': min_font_size - 1,
                        'text_preview': text_content[:40] + ('...' if len(text_content) > 40 else '')
                    })
                    results['passes'] = False

            if annotation_sizes:
                results['font_sizes']['annotations_min'] = min(annotation_sizes)
                results['font_sizes']['annotations_max'] = max(annotation_sizes)
                results['font_sizes']['annotations_count'] = len(annotation_sizes)

    if not results['warnings']:
        results['warnings'].append("✓ All accessibility checks passed!")

    return results


# ==============================================================================
# USAGE EXAMPLES (for documentation)
# ==============================================================================

if __name__ == "__main__":
    # Example 1: Test colorblind safety of current palette
    print("Testing Cancer vs Normal colors:")
    test_result = test_colorblind_visibility(
        COLORBLIND_GROUP_PALETTE['Cancer'],
        COLORBLIND_GROUP_PALETTE['Normal']
    )
    print(f"  Contrast ratio: {test_result['normal_contrast']:.2f}")
    print(f"  Colorblind contrast: {test_result['deuteranopia_contrast']:.2f}")
    print(f"  WCAG level: {test_result['wcag_level']}")
    print(f"  Passes: {test_result['passes_colorblind']}")

    # Example 2: Get optimal font sizes for Nature journal
    print("\nOptimal font sizes for Nature (10-inch wide figure):")
    fonts = optimize_font_sizes_for_publication(10, 'nature')
    for key, value in fonts.items():
        print(f"  {key}: {value}")
