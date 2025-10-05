# Scripts Directory

Utility scripts for data verification and maintenance.

## Available Scripts

### verify_trace_data.py

**Purpose**: Comprehensive verification of glycopeptide comparison heatmap trace data.

**Usage**:
```bash
python3 scripts/verify_trace_data.py
```

**What it checks**:
1. ✓ Cancer_Mean calculations (across C1-C24 samples)
2. ✓ Normal_Mean calculations (across N1-N23 samples)
3. ✓ VIP score sorting (descending order)
4. ✓ Glycan type grouping (contiguous positions)
5. ✓ Fold change calculations
6. ✓ Plot flags (Cancer_Dot_Plotted, Normal_Dot_Plotted)
7. ✓ Alpha transparency values (0-1 range)
8. ✓ Sample counts (≤ max samples)

**Output**:
- Detailed verification report
- PASS/FAIL status for each check
- Exit code 0 if all checks pass, 1 if any fail

**Requirements**:
- Must run after generating comparison heatmap
- Requires trace data files in `Results/Trace/`:
  - `glycopeptide_comparison_heatmap_summary.csv`
  - `glycopeptide_comparison_heatmap_data.csv`

---

## Adding New Scripts

When adding utility scripts:

1. Place in `scripts/` directory
2. Add shebang: `#!/usr/bin/env python3`
3. Include docstring explaining purpose
4. Update this README
5. Make executable: `chmod +x scripts/your_script.py`

---

## Directory Structure

```
scripts/
├── README.md              # This file
└── verify_trace_data.py   # Trace data verification
```

---

**Last Updated**: 2025-10-05
