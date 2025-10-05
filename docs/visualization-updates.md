# Glycopeptide Comparison Heatmap - Latest Updates

## Visual Comparison

### Before (Split-Cell)
```
┌─────────┬─────────┐
│ Cancer  │ Normal  │
│ (Red)   │ (Blue)  │
│ Left    │ Right   │
└─────────┴─────────┘
```
❌ Less intuitive
❌ Overlapping not clear

### After (Symbol-Based) ✓
```
Grid intersections:
     ×  = Cancer (Red)
     +  = Normal (Blue)

Example:
─────┼─────┼─────┼─────
     │  ×+ │     │  ×
─────┼─────┼─────┼─────
     │     │  +  │  ×+
─────┼─────┼─────┼─────
```
✓ Clear symbols (× vs +)
✓ Overlapping visible (both shown)
✓ Aligned to grid

---

## Key Features

### 1. Symbol Markers
- **× (cross)** = Cancer group (red #E74C3C)
- **+ (plus)** = Normal/Control group (blue #3498DB)
- Placed exactly on grid intersections
- Size: 300 points, linewidth: 2.5
- Zorder: × on top (4) when overlapping with + (3)

### 2. Transparency Encoding
- **Darkest** symbol → highest intensity in that glycan type
- **Lightest** symbol → lowest intensity in that glycan type
- Range: alpha 0.3 to 1.0
- **Relative within each glycan type** (HM, F, S, SF, C/H)

### 3. Glycan Sorting
Glycan compositions sorted **ascending** within each type:

**HM group**:
- H(12)N(2)
- H(5)N(2)
- H(6)N(2)
- H(7)N(2)
- H(8)N(2)

**F group**:
- H(3)N(3)F(1)
- H(3)N(5)F(1)
- H(4)N(4)F(1)
- H(5)N(4)F(1)
- ...

### 4. Grid Lines
- Visible grid (alpha 0.4, gray #CCCCCC)
- Solid lines, width 0.8
- Helps track symbols across rows/columns

### 5. Color Bar
- **Solid color blocks** per glycan type
- No gradient (clear boundaries)
- HM (green) → F (red) → S (pink) → SF (orange) → C/H (blue)

---

## Interpretation Guide

### Reading the Heatmap

1. **Find your peptide** (Y-axis, sorted by VIP score)
2. **Find your glycan** (X-axis, grouped by type, sorted alphabetically)
3. **Look at the symbol**:
   - **×** present? → Detected in Cancer
   - **+** present? → Detected in Normal
   - **Both ×+**? → Detected in both groups
   - **Darkness** → Relative intensity within that glycan type

### Example Patterns

**Cancer-specific**:
```
× only (no +) → Only in cancer samples
```

**Normal-specific**:
```
+ only (no ×) → Only in normal samples
```

**Both groups, higher in cancer**:
```
Dark × + Light + → Higher in cancer
```

**Both groups, higher in normal**:
```
Light × + Dark + → Higher in normal
```

---

## Technical Details

### Code Changes

**Symbol rendering** (`glycopeptide_comparison_heatmap.py` lines 229-246):
```python
# Cancer symbol (× cross) - Red
ax_main.scatter(x_pos, y_pos, s=300, c='#E74C3C', alpha=alpha,
              marker='x', linewidths=2.5, zorder=4)

# Normal symbol (+ plus) - Blue
ax_main.scatter(x_pos, y_pos, s=300, c='#3498DB', alpha=alpha,
              marker='+', linewidths=2.5, zorder=3)
```

**Glycan sorting** (lines 118-119):
```python
# Sort glycan compositions in ascending order within each type group
type_glycans_sorted = sorted(type_glycans)
```

**Grid** (line 278):
```python
ax_main.grid(True, alpha=0.4, linestyle='-', linewidth=0.8,
             color='#CCCCCC', zorder=0)
```

---

## Legend

**Glycan Types**:
- HM (green) - High-mannose
- F (red) - Fucosylated
- S (pink) - Sialylated
- SF (orange) - Sialofucosylated
- C/H (blue) - Complex/Hybrid

**Group Markers**:
- × Cancer (red cross)
- + Normal/Control (blue plus)

**Intensity**:
- Symbol darkness = Intensity (relative to glycan type max)

---

## Benefits

1. **Clear Distinction**: × vs + immediately recognizable
2. **Perfect Alignment**: Symbols on grid intersections (not offset)
3. **Overlapping Visible**: When both groups present, both symbols shown
4. **Organized**: Glycans sorted alphabetically within each type
5. **Traceable**: All data in `Results/Trace/` for verification

---

## Files Generated

1. **Visualization**: `Results/glycopeptide_comparison_heatmap.png`
2. **Summary data**: `Results/Trace/glycopeptide_comparison_heatmap_summary.csv` (22 columns)
3. **Full data**: `Results/Trace/glycopeptide_comparison_heatmap_data.csv` (69 columns)

---

**Last Updated**: 2025-10-05
**Status**: All updates implemented and tested ✓
