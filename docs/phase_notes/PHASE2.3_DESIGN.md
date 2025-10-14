# Phase 2.3: Unified Base Method Design Document

## Objective

Create a single unified base method `_plot_comparison_heatmap_base()` that consolidates the three public methods:
1. `plot_glycopeptide_comparison_heatmap` (standard, top N)
2. `plot_glycopeptide_comparison_heatmap_full` (full-scale, ALL)
3. `plot_glycopeptide_comparison_heatmap_by_type` (type-specific, ALL of one type)

**Expected Impact**: ~200-300 lines reduction by eliminating duplicate plotting logic.

---

## Method Comparison Analysis

### Method 1: Standard Heatmap (Top N)

**Purpose**: Show top discriminative glycopeptides in a readable format

**Data Selection**:
- Top `max_peptides` unique peptides by VIP score (default: 50)
- Top `max_glycans_per_type` glycans per type (default: 15)
- Multi-type: HM, F, S, SF, C/H

**Figure Properties**:
- Size: Fixed `(24, 16)`
- Layout: 3-panel (top + colorbar + main), height_ratios=[1, 0.3, 5.5], hspace=0.05
- Marker size: `'large'`
- Linewidth: bold=5.0, normal=PLOT_LINE_LINEWIDTH_THICK

**Legend**:
- Glycan type color patches (5 types)
- Symbol indicators with linewidth distinction (4 items)
- Total: 9 legend elements

**Title**: `"Glycopeptide Comparison Heatmap: Cancer vs Normal\n(Top {N} peptides by VIP score)"`

**Output**: `glycopeptide_comparison_heatmap.png`

---

### Method 2: Full-Scale Heatmap (ALL)

**Purpose**: Show complete glycopeptide landscape without filtering

**Data Selection**:
- ALL glycopeptides (no filtering)
- Multi-type: HM, F, S, SF, C/H

**Figure Properties**:
- Size: Dynamic, `max(30, n_peptides * 0.2 + 5)` × `max(30, n_glycans * 0.23 + 10)`
- Layout: 3-panel (top + colorbar + main), height_ratios=[1, 0.3, 5.5], hspace=0.05
- Marker size: `'small'`
- Linewidth: bold=3.0, normal=EDGE_LINEWIDTH_THICK

**Legend**:
- Glycan type color patches (5 types)
- Symbol indicators with linewidth distinction (4 items)
- Total: 9 legend elements

**Title**: `"Full-Scale Glycopeptide Comparison Heatmap: Cancer vs Normal\n(ALL {N} peptides × {M} glycans = {total} glycopeptides)"`

**Output**: `glycopeptide_comparison_heatmap_full.png`

---

### Method 3: Type-Specific Heatmap (ALL of one type)

**Purpose**: Focus on single glycan type for detailed comparison

**Data Selection**:
- ALL glycopeptides filtered to specific `glycan_type` (HM, F, S, SF, or C/H)
- Single-type only

**Figure Properties**:
- Size: Dynamic, `max(20, n_peptides * 0.2 + 5)` × `max(20, n_glycans * 0.23 + 10)`
- Layout: 2-panel (top + main), height_ratios=[1, 6], hspace=0.08, **NO colorbar**
- Marker size: `'small'`
- Linewidth: bold=3.0, normal=EDGE_LINEWIDTH_THICK

**Legend**:
- NO glycan type color patches (single color only)
- Symbol indicators with single color (4 items)
- Total: 4 legend elements

**Title**: `"Glycopeptide Comparison: {Type Name} Type\nCancer vs Normal ({N} peptides × {M} glycans = {total} glycopeptides)"`

**Output**: `glycopeptide_comparison_heatmap_{type}.png` (e.g., `_HM`, `_F`, `_C_H`)

**Additional Parameters**:
- `glycan_type`: Specific type to filter (str)
- `glycan_type_color`: Color for that type
- `glycan_type_names`: Dict mapping type code to full name

---

## Key Differences Summary

| Feature | Method 1 (Standard) | Method 2 (Full) | Method 3 (By Type) |
|---------|---------------------|-----------------|---------------------|
| **Data Selection** | Top N peptides + glycans | ALL | ALL of one type |
| **Figure Size** | Fixed (24, 16) | Dynamic (min 30×30) | Dynamic (min 20×20) |
| **Panels** | 3 (top + colorbar + main) | 3 (top + colorbar + main) | 2 (top + main) |
| **Height Ratios** | [1, 0.3, 5.5] | [1, 0.3, 5.5] | [1, 6] |
| **Hspace** | 0.05 | 0.05 | 0.08 |
| **Marker Size** | 'large' | 'small' | 'small' |
| **Linewidth Bold** | 5.0 | 3.0 | 3.0 |
| **Linewidth Normal** | PLOT_LINE_LINEWIDTH_THICK | EDGE_LINEWIDTH_THICK | EDGE_LINEWIDTH_THICK |
| **Glycan Types** | Multiple (5) | Multiple (5) | Single (1) |
| **Colorbar** | Yes | Yes | No |
| **Legend Elements** | 9 (patches + symbols) | 9 (patches + symbols) | 4 (symbols only) |
| **Title Prefix** | "Glycopeptide Comparison Heatmap" | "Full-Scale Glycopeptide Comparison Heatmap" | "Glycopeptide Comparison: {Type Name} Type" |
| **Filename** | glycopeptide_comparison_heatmap | glycopeptide_comparison_heatmap_full | glycopeptide_comparison_heatmap_{type} |

---

## Unified Base Method Design

### Method Signature

```python
def _plot_comparison_heatmap_base(
    self,
    df_with_vip: pd.DataFrame,
    heatmap_variant: str,  # 'standard', 'full', or 'by_type'
    config: DataPreparationConfig,
    # Variant-specific parameters
    max_peptides: int = None,  # For 'standard' only
    max_glycans_per_type: int = None,  # For 'standard' only
    glycan_type: str = None,  # For 'by_type' only
    figsize: tuple = None  # For 'standard' only (others use dynamic)
):
    """
    Unified base method for all glycopeptide comparison heatmaps

    Args:
        df_with_vip: DataFrame with VIP scores already merged
        heatmap_variant: 'standard' (top N), 'full' (all), or 'by_type' (one type)
        config: Data preparation configuration
        max_peptides: Max peptides for standard variant (default: 50)
        max_glycans_per_type: Max glycans per type for standard variant (default: 15)
        glycan_type: Specific type for by_type variant ('HM', 'F', 'S', 'SF', 'C/H')
        figsize: Fixed figure size for standard variant (default: (24, 16))
    """
```

### Conditional Logic Flow

```python
# 1. Data Selection Strategy
if heatmap_variant == 'standard':
    # Top N peptides + top N glycans per type
    df_plot = select_top_glycopeptides(df_with_vip, max_peptides, max_glycans_per_type)
    use_multi_type = True
elif heatmap_variant == 'full':
    # ALL glycopeptides
    df_plot = df_with_vip.copy()
    use_multi_type = True
elif heatmap_variant == 'by_type':
    # ALL glycopeptides of specific type
    df_plot = df_with_vip[df_with_vip['GlycanTypeCategory'] == glycan_type].copy()
    use_multi_type = False

# 2. Glycan Ordering
if use_multi_type:
    glycan_type_order = ['HM', 'F', 'S', 'SF', 'C/H']
    glycan_type_colors = EXTENDED_CATEGORY_COLORS
    glycan_order, glycan_type_positions = create_glycan_order(df_plot, glycan_type_order, include_positions=True)
else:
    # Single type - no type grouping
    glycan_compositions = df_plot['GlycanComposition'].unique()
    glycan_order = sorted(glycan_compositions, key=glycan_sort_key)
    glycan_type_positions = None
    glycan_type_color = EXTENDED_CATEGORY_COLORS[glycan_type]

# 3. Figure Sizing
if heatmap_variant == 'standard':
    # Fixed size
    fig = plt.figure(figsize=figsize or (24, 16))
else:
    # Dynamic sizing
    n_peptides = len(peptide_order)
    n_glycans = len(glycan_order)
    if heatmap_variant == 'full':
        min_size = 30
    else:  # by_type
        min_size = 20
    fig_height = max(min_size, n_peptides * 0.2 + 5)
    fig_width = max(min_size, n_glycans * 0.23 + 10)
    fig = plt.figure(figsize=(fig_width, fig_height))

# 4. Panel Layout
if use_multi_type:
    # 3-panel layout with colorbar
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 0.3, 5.5], hspace=0.05)
    ax_top = fig.add_subplot(gs[0])
    ax_colorbar = fig.add_subplot(gs[1])
    ax_main = fig.add_subplot(gs[2])
else:
    # 2-panel layout without colorbar
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 6], hspace=0.08)
    ax_top = fig.add_subplot(gs[0])
    ax_main = fig.add_subplot(gs[1])

# 5. Styling Parameters
if heatmap_variant == 'standard':
    marker_size = 'large'
    linewidth_bold = 5.0
    linewidth_normal = PLOT_LINE_LINEWIDTH_THICK
    label_size = 'large'
else:  # 'full' or 'by_type'
    marker_size = 'small'
    linewidth_bold = 3.0
    linewidth_normal = EDGE_LINEWIDTH_THICK
    label_size = 'small'

# 6. Top Panel
self._plot_top_panel(ax_top, glycan_order, df_plot, marker_size=marker_size)

# 7. Colorbar Panel (conditional)
if use_multi_type:
    self._plot_colorbar_panel(ax_colorbar, glycan_order, glycan_type_positions, glycan_type_colors, label_size=label_size)

# 8. Symbol Heatmap
if use_multi_type:
    self._plot_symbol_heatmap(
        ax_main, df_plot, peptide_to_idx, glycan_to_idx, peptide_order, glycan_order,
        marker_size=marker_size, linewidth_bold=linewidth_bold, linewidth_normal=linewidth_normal,
        glycan_type_colors=glycan_type_colors, glycan_type_positions=glycan_type_positions
    )
else:
    self._plot_symbol_heatmap(
        ax_main, df_plot, peptide_to_idx, glycan_to_idx, peptide_order, glycan_order,
        marker_size=marker_size, linewidth_bold=linewidth_bold, linewidth_normal=linewidth_normal,
        fixed_color=glycan_type_color
    )

# 9. Grid (same for all)
# ... grid code

# 10. Legend (conditional)
if use_multi_type:
    # Create legend with glycan type patches + symbols
    legend_elements = []
    for gt in glycan_type_order:
        if gt in glycan_type_positions:
            legend_elements.append(Patch(...))
    # Add symbol indicators
    legend_elements.extend([...])
else:
    # Create legend with symbols only (single color)
    legend_elements = [
        Line2D([0], [0], marker=MARKER_ANNOTATION, color='w',
               markerfacecolor=glycan_type_color, ...)
    ]

# 11. Title (conditional)
if heatmap_variant == 'standard':
    title = f"Glycopeptide Comparison Heatmap: Cancer vs Normal\n(Top {len(peptide_order)} peptides by VIP score)"
elif heatmap_variant == 'full':
    title = f"Full-Scale Glycopeptide Comparison Heatmap: Cancer vs Normal\n(ALL {len(peptide_order)} peptides × {len(glycan_order)} glycans = {len(df_plot)} glycopeptides)"
else:  # by_type
    glycan_type_names = {'HM': 'High-Mannose', 'F': 'Fucosylated', ...}
    title = f"Glycopeptide Comparison: {glycan_type_names[glycan_type]} Type\nCancer vs Normal ({len(peptide_order)} peptides × {len(glycan_order)} glycans = {len(df_plot)} glycopeptides)"

# 12. Save
if heatmap_variant == 'standard':
    filename = 'glycopeptide_comparison_heatmap.png'
elif heatmap_variant == 'full':
    filename = 'glycopeptide_comparison_heatmap_full.png'
else:  # by_type
    filename = f'glycopeptide_comparison_heatmap_{glycan_type.replace("/", "_")}.png'

# 13. Trace Data
self._prepare_and_save_trace_data(
    df_plot, peptide_to_idx, glycan_to_idx,
    filename.replace('.png', ''),
    linewidth_bold=linewidth_bold,
    linewidth_normal=linewidth_normal,
    include_symbol_color=(not use_multi_type),
    symbol_color=(glycan_type_color if not use_multi_type else None)
)
```

---

## Helper Method Needed

### `_select_top_glycopeptides()` (for standard variant)

```python
def _select_top_glycopeptides(
    self,
    df_with_vip: pd.DataFrame,
    max_peptides: int,
    max_glycans_per_type: int,
    glycan_type_order: list
) -> pd.DataFrame:
    """
    Select top N peptides and top N glycans per type

    Logic from current method 1:
    - Sort by VIP score descending
    - Get first max_peptides unique peptides
    - Include all glycoforms of those peptides
    - For each glycan type, get top max_glycans_per_type glycans
    """
```

---

## Refactored Public Methods (Phase 2.4)

Each public method becomes a thin wrapper:

### Method 1: Standard

```python
def plot_glycopeptide_comparison_heatmap(
    self, df: pd.DataFrame, vip_scores: pd.DataFrame,
    config: DataPreparationConfig = None,
    figsize: tuple = (24, 16),
    max_peptides: int = 50,
    max_glycans_per_type: int = 15
):
    """Create dot-based heatmap comparing Cancer vs Normal groups (top N)"""
    logger.info("Creating glycopeptide comparison heatmap (Cancer vs Normal)...")

    # Use default config if not provided
    if config is None:
        config = DataPreparationConfig(min_detection_pct=0.30, min_samples=5, missing_data_method='skipna')

    # Prepare data with VIP scores
    df_with_vip = prepare_visualization_data(
        df=df, config=config, vip_scores=vip_scores,
        merge_method='left', apply_detection_filter=False,
        log_prefix="[Comparison Heatmap] "
    )

    if len(df_with_vip) == 0:
        logger.error("No glycopeptides available!")
        return

    # Call unified base method
    self._plot_comparison_heatmap_base(
        df_with_vip=df_with_vip,
        heatmap_variant='standard',
        config=config,
        max_peptides=max_peptides,
        max_glycans_per_type=max_glycans_per_type,
        figsize=figsize
    )
```

**Estimated size**: ~25 lines (down from ~200+ lines)

### Method 2: Full

```python
def plot_glycopeptide_comparison_heatmap_full(
    self, df: pd.DataFrame, vip_scores: pd.DataFrame,
    config: DataPreparationConfig = None
):
    """Create FULL-SCALE heatmap comparing Cancer vs Normal groups (ALL)"""
    logger.info("Creating FULL-SCALE glycopeptide comparison heatmap (ALL glycopeptides)...")

    # Use default config if not provided
    if config is None:
        config = DataPreparationConfig(min_detection_pct=0.30, min_samples=5, missing_data_method='skipna')

    # Prepare data with VIP scores
    df_with_vip = prepare_visualization_data(
        df=df, config=config, vip_scores=vip_scores,
        merge_method='left', apply_detection_filter=False,
        log_prefix="[Full Comparison Heatmap] "
    )

    if len(df_with_vip) == 0:
        logger.error("No glycopeptides available!")
        return

    logger.info(f"Processing {len(df_with_vip)} total glycopeptides (complete dataset)")

    # Call unified base method
    self._plot_comparison_heatmap_base(
        df_with_vip=df_with_vip,
        heatmap_variant='full',
        config=config
    )
```

**Estimated size**: ~30 lines (down from ~200+ lines)

### Method 3: By Type

```python
def plot_glycopeptide_comparison_heatmap_by_type(
    self, df: pd.DataFrame, vip_scores: pd.DataFrame,
    glycan_type: str,
    config: DataPreparationConfig = None
):
    """Create heatmap for a SPECIFIC glycan type (ALL of one type)"""
    logger.info(f"Creating glycopeptide comparison heatmap for {glycan_type} type...")

    # Use default config if not provided
    if config is None:
        config = DataPreparationConfig(min_detection_pct=0.30, min_samples=5, missing_data_method='skipna')

    # Prepare data with VIP scores
    df_with_vip = prepare_visualization_data(
        df=df, config=config, vip_scores=vip_scores,
        merge_method='left', apply_detection_filter=False,
        log_prefix=f"[{glycan_type} Heatmap] "
    )

    if len(df_with_vip) == 0:
        logger.error("No glycopeptides available!")
        return

    # Call unified base method
    self._plot_comparison_heatmap_base(
        df_with_vip=df_with_vip,
        heatmap_variant='by_type',
        config=config,
        glycan_type=glycan_type
    )
```

**Estimated size**: ~30 lines (down from ~180+ lines)

---

## Expected Line Count Reduction

**Current State** (after Phase 2.2b): 1,011 lines

**Breakdown**:
- Helper methods (4): ~418 lines
- Method 1 (standard): ~200 lines
- Method 2 (full): ~200 lines
- Method 3 (by_type): ~180 lines

**After Phase 2.3-2.4**:
- Helper methods (4): ~418 lines (unchanged)
- Helper method (1 new): ~50 lines (for data selection logic)
- Unified base method: ~250 lines (consolidates common plotting logic)
- Method 1 wrapper: ~25 lines
- Method 2 wrapper: ~30 lines
- Method 3 wrapper: ~30 lines

**Estimated total**: ~803 lines (208 lines reduction from 1,011)

**Final size**: ~803 lines (37% reduction from original 1,273 lines)

**Progress towards 600 line target**: 203 lines above target (acceptable given complexity)

---

## Implementation Plan

### Phase 2.3.1: Extract Data Selection Helper
Create `_select_top_glycopeptides()` method to handle standard variant's data selection logic.

### Phase 2.3.2: Create Unified Base Method
Implement `_plot_comparison_heatmap_base()` with all conditional logic branches.

### Phase 2.3.3: Test Base Method
Run pipeline with base method to ensure correctness before refactoring public methods.

### Phase 2.4: Refactor Public Methods
Convert all 3 public methods to thin wrappers calling the base method.

### Phase 2.5: Final Testing
- Full pipeline execution
- Visual comparison of all 7 PNG files
- CSV trace data validation
- Performance benchmarking

---

## Risks and Mitigation

**Risk 1**: Increased complexity in base method
- **Mitigation**: Clear documentation, well-structured conditional logic, comprehensive comments

**Risk 2**: Subtle behavioral differences between variants
- **Mitigation**: Thorough testing, visual comparison of before/after outputs

**Risk 3**: Harder to debug issues
- **Mitigation**: Detailed logging at each conditional branch

**Risk 4**: May not reach 600 line target
- **Mitigation**: Accept ~800 lines as success (37% reduction is significant)

---

## Status

**Phase 2.3**: Design complete, ready for implementation ✅
**Next**: Begin implementation of Phase 2.3.1
