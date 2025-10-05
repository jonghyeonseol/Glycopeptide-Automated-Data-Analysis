# Glycan Composition Sorting - Technical Guide

## Overview

Glycan compositions are sorted **numerically within each glycan type group** using natural sorting algorithm.

---

## Sorting Algorithm

### Logic

1. Extract monosaccharide components and their counts
2. Sort by monosaccharide order: **H → N → A → F → G**
3. Within each monosaccharide, sort **numerically** (not alphabetically)

### Example: Why Numeric Sorting Matters

**Alphabetical (Wrong)**:
- H(12)N(2) comes before H(5)N(2) (because "1" < "5" in strings)

**Numeric (Correct)**:
- H(5)N(2) comes before H(12)N(2) (because 5 < 12 numerically)

---

## Sorting Results by Type

### HM (High-Mannose)
```
Position 0: H(5)N(2)
Position 1: H(6)N(2)
Position 2: H(7)N(2)
Position 3: H(8)N(2)
Position 4: H(12)N(2)
```
✓ Sorted by H count: 5 < 6 < 7 < 8 < 12

### F (Fucosylated)
```
Position 5:  H(3)N(3)F(1)
Position 6:  H(3)N(4)F(1)
Position 7:  H(3)N(5)F(1)
Position 8:  H(4)N(3)F(1)
Position 9:  H(4)N(4)F(1)
Position 10: H(4)N(5)F(1)
Position 11: H(5)N(4)F(1)
Position 12: H(5)N(5)F(1)
```
✓ First by H: 3,3,3,4,4,4,5,5
✓ Then by N: 3,4,5,3,4,5,4,5
✓ F is constant (all =1)

### S (Sialylated)
```
Position 13: H(5)N(4)A(1)
Position 14: H(5)N(4)A(2)
Position 15: H(6)N(3)A(1)
Position 16: H(6)N(5)A(3)
```
✓ First by H: 5,5,6,6
✓ Then by N: 4,4,3,5
✓ Then by A: 1,2,1,3

### SF (Sialofucosylated)
```
Position 17: H(4)N(4)A(1)F(1)
Position 18: H(5)N(4)A(1)F(1)
Position 19: H(5)N(4)A(1)F(2)
Position 20: H(5)N(4)A(1)F(3)
Position 21: H(5)N(4)A(2)F(1)
Position 22: H(5)N(5)A(1)F(1)
Position 23: H(6)N(5)A(3)F(1)
```
✓ H: 4,5,5,5,5,5,6
✓ N: 4,4,4,4,4,5,5
✓ A: 1,1,1,1,2,1,3
✓ F: 1,1,2,3,1,1,1

### C/H (Complex/Hybrid)
```
Position 24: H(3)N(3)
Position 25: H(3)N(4)
Position 26: H(4)N(3)
Position 27: H(4)N(4)
Position 28: H(4)N(5)
Position 29: H(5)N(3)
Position 30: H(5)N(4)
Position 31: H(5)N(5)
Position 32: H(6)N(3)
Position 33: H(6)N(4)
Position 34: H(6)N(5)
```
✓ H: 3,3,4,4,4,5,5,5,6,6,6
✓ N: 3,4,3,4,5,3,4,5,3,4,5

---

## Code Implementation

### Sort Key Function

```python
def glycan_sort_key(glycan_comp):
    """Extract numbers from glycan composition for natural sorting"""
    import re
    # Extract monosaccharide types and counts: e.g., H(5)N(4)A(2)
    parts = re.findall(r'([A-Z]+)\((\d+)\)', glycan_comp)

    # Define monosaccharide order
    monosaccharide_order = {'H': 0, 'N': 1, 'A': 2, 'F': 3, 'G': 4}

    # Create sort tuple: (mono_order, numeric_count)
    sort_tuple = []
    for mono, count in parts:
        order = monosaccharide_order.get(mono, 99)
        sort_tuple.append((order, int(count)))  # int() for numeric sort

    return tuple(sort_tuple)
```

### Usage

```python
glycans = ['H(12)N(2)', 'H(5)N(2)', 'H(6)N(2)', 'H(7)N(2)']
sorted_glycans = sorted(glycans, key=glycan_sort_key)
# Result: ['H(5)N(2)', 'H(6)N(2)', 'H(7)N(2)', 'H(12)N(2)']
```

---

## Verification

### Manual Check

Open `Results/Trace/glycopeptide_comparison_heatmap_summary.csv`:

```csv
GlycanComposition,GlycanTypeCategory,Plot_X_Position
H(5)N(2),HM,0
H(6)N(2),HM,1
H(7)N(2),HM,2
H(8)N(2),HM,3
H(12)N(2),HM,4
```

✓ X positions are in ascending order
✓ Glycan compositions are numerically sorted

### Python Verification

```python
import pandas as pd

df = pd.read_csv('Results/Trace/glycopeptide_comparison_heatmap_summary.csv')

for gtype in ['HM', 'F', 'S', 'SF', 'C/H']:
    glycans = df[df['GlycanTypeCategory']==gtype][
        ['GlycanComposition', 'Plot_X_Position']
    ].drop_duplicates().sort_values('Plot_X_Position')

    print(f'\n{gtype}:')
    for _, row in glycans.iterrows():
        print(f"  {row['GlycanComposition']} at position {int(row['Plot_X_Position'])}")
```

---

## Benefits

1. **Intuitive Order**: Smaller numbers appear first
2. **Consistent**: Same sorting logic across all glycan types
3. **Traceable**: Easy to find specific compositions on x-axis
4. **Biological**: Follows typical glycan biosynthesis progression

---

## File Location

**Implementation**: `src/plots/glycopeptide_comparison_heatmap.py` (lines 115-142)

**Key Code**:
- Line 116-129: `glycan_sort_key()` function
- Line 135: Apply sorting with `sorted(type_glycans, key=glycan_sort_key)`

---

**Last Updated**: 2025-10-05
**Status**: Implemented and verified ✓
