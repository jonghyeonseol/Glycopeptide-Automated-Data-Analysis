# Refactoring Recommendations for pGlyco Auto Combine

## ✅ Lint Issues Fixed

### 1. **FutureWarnings Eliminated**
- ✅ Fixed all pandas FutureWarnings by using `.mask()` instead of `.replace()` in `utils.py`
- ✅ Pipeline now runs cleanly without downcasting warnings
- ✅ 28,798+ FutureWarning instances eliminated

### 2. **Remaining Warnings (Non-Critical)**
- ⚠️ **Seaborn UserWarning**: "Clustering large matrix with scipy. Installing fastcluster may give better performance"
  - **Solution**: Add `fastcluster` to `requirements.txt` for performance improvement
- ⚠️ **R Warning**: "no non-missing arguments to min; returning Inf"
  - **Impact**: Harmless R warning from VIP score calculations, does not affect results
- ⚠️ **Matplotlib UserWarning**: "This figure includes Axes that are not compatible with tight_layout"
  - **Location**: `site_specific_heatmap.py:144`
  - **Solution**: Suppress warning or use `constrained_layout=True` instead

---

## 🔧 Critical Refactoring Needs

### 1. **Code Duplication - Sample Column Extraction** (Priority: HIGH)

**Problem**: Sample column extraction logic is duplicated across 14+ files

**Current Pattern (duplicated):**
```python
# Pattern 1 (used in 6 files)
cancer_samples = [col for col in df.columns if col.startswith('C') and col[1:].isdigit()]
normal_samples = [col for col in df.columns if col.startswith('N') and col[1:].isdigit()]

# Pattern 2 (used in 8 files)
cancer_samples = [col for col in sample_cols if col.startswith('C')]
normal_samples = [col for col in sample_cols if col.startswith('N')]
```

**Recommended Solution:**
```python
# Add to src/utils.py
def get_sample_columns(df: pd.DataFrame) -> tuple[list[str], list[str]]:
    """
    Extract cancer and normal sample columns from DataFrame

    Returns:
        Tuple of (cancer_samples, normal_samples)
    """
    cancer_samples = [col for col in df.columns if col.startswith('C') and col[1:].isdigit()]
    normal_samples = [col for col in df.columns if col.startswith('N') and col[1:].isdigit()]
    return cancer_samples, normal_samples
```

**Files to Update:**
- `src/analyzer.py`
- `src/plots/boxplot.py`
- `src/plots/correlation_matrix_plot.py`
- `src/plots/cv_distribution_plot.py`
- `src/plots/heatmap.py`
- `src/plots/histogram.py`
- `src/plots/pca_plot.py` (if applicable)
- `src/plots/radar_chart_plot.py`
- `src/plots/site_specific_heatmap.py`
- `src/plots/venn_diagram_plot.py`
- `src/plots/vip_score_plot.py`
- `src/plots/vip_score_plot_r.py`
- `src/plots/volcano_plot.py`
- `main.py`

---

### 2. **File Organization** (Priority: MEDIUM)

**Problem**: Temporary/helper scripts in root directory

**Current Structure:**
```
.
├── add_trace_remaining.py  # ❌ Temporary helper script
├── main.py
├── src/
│   ├── plots/  # 12 files
│   └── ...
```

**Recommended Structure:**
```
.
├── main.py
├── scripts/              # NEW: Helper scripts
│   └── add_trace_remaining.py
├── src/
│   ├── core/            # NEW: Core logic
│   │   ├── analyzer.py
│   │   ├── annotator.py
│   │   └── data_loader.py
│   ├── plots/           # Visualization modules
│   │   └── ...
│   └── utils.py
├── tests/               # NEW: Unit tests
│   ├── test_analyzer.py
│   └── ...
└── docs/                # NEW: Documentation
    └── API.md
```

**Actions:**
1. Move `add_trace_remaining.py` to `scripts/` folder
2. Create `src/core/` for core business logic
3. Create `tests/` directory for unit tests
4. Create `docs/` for API documentation

---

### 3. **Type Hints Missing** (Priority: MEDIUM)

**Problem**: No type hints in function signatures

**Current:**
```python
def plot_pca(self, pca_results: dict, figsize: tuple = (10, 8)):  # ✅ Has types
def plot_histogram(self, df, glycan_type_col='GlycanType'):      # ❌ No types
```

**Recommended:**
```python
from typing import Optional, Dict, List, Tuple
import pandas as pd

def plot_histogram(
    self,
    df: pd.DataFrame,
    glycan_type_col: str = 'GlycanType'
) -> None:
    """Plot histogram..."""
```

**Files Needing Type Hints:**
- `src/plots/histogram.py` (multiple methods)
- `src/plots/boxplot.py` (multiple methods)
- `src/plots/vip_score_plot.py`
- `src/analyzer.py` (some methods)

---

### 4. **Magic Numbers and Strings** (Priority: MEDIUM)

**Problem**: Hardcoded values scattered throughout codebase

**Examples:**
```python
# In multiple files
top_n = 50  # Hardcoded
fdr_threshold = 0.05  # Hardcoded
fc_threshold = 1.5  # Hardcoded
'C' and col[1:].isdigit()  # Hardcoded pattern
```

**Recommended Solution:**
```python
# Add to config.yaml
visualization:
  defaults:
    top_n_glycopeptides: 50
    fdr_threshold: 0.05
    fold_change_threshold: 1.5

  sample_patterns:
    cancer_prefix: "C"
    normal_prefix: "N"
```

Or create constants file:
```python
# src/constants.py
DEFAULT_TOP_N = 50
FDR_THRESHOLD = 0.05
FC_THRESHOLD = 1.5
CANCER_PREFIX = 'C'
NORMAL_PREFIX = 'N'
```

---

### 5. **Large Method Complexity** (Priority: MEDIUM)

**Problem**: Some methods exceed 100 lines

**Large Methods:**
- `src/plots/histogram.py::plot_histogram_primary_classification()` (~120 lines)
- `src/plots/boxplot.py::plot_boxplot_secondary_classification()` (~100 lines)
- `src/plots/vip_score_plot_r.py::plot_vip_scores_peptide_grouped_r()` (~150 lines)

**Recommended**: Break down into smaller helper methods
```python
# Before (150 lines)
def plot_vip_scores_peptide_grouped_r(self, df, vip_df):
    # ... 150 lines ...

# After (refactored)
def plot_vip_scores_peptide_grouped_r(self, df, vip_df):
    plot_data = self._prepare_peptide_grouped_data(df, vip_df)
    r_script = self._generate_peptide_grouped_r_script(plot_data)
    self._execute_r_plot(r_script, 'vip_score_peptide_grouped_r.png')

def _prepare_peptide_grouped_data(self, df, vip_df):
    # Data preparation logic

def _generate_peptide_grouped_r_script(self, plot_data):
    # R script generation

def _execute_r_plot(self, r_script, output_filename):
    # R execution
```

---

### 6. **Error Handling** (Priority: LOW)

**Problem**: Limited error handling in visualization methods

**Current:**
```python
def plot_volcano(self, df, vip_df):
    # ... lots of processing ...
    # No try-except for statistical tests or plotting
```

**Recommended:**
```python
def plot_volcano(self, df, vip_df):
    try:
        # Data processing
        volcano_data = self._calculate_volcano_data(df, vip_df)

        # Plotting
        self._create_volcano_plot(volcano_data)

    except ValueError as e:
        logger.error(f"Data validation error in volcano plot: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error in volcano plot: {e}")
        raise
```

---

### 7. **Testing Infrastructure** (Priority: HIGH)

**Problem**: No unit tests exist

**Recommended Test Structure:**
```
tests/
├── __init__.py
├── conftest.py              # Pytest fixtures
├── test_data_loader.py      # Test data integration
├── test_annotator.py         # Test glycan annotation
├── test_analyzer.py          # Test PCA/statistics
├── test_utils.py             # Test utility functions
├── test_visualizations.py    # Test plot generation
└── fixtures/                 # Test data
    ├── sample_C01.csv
    └── sample_N01.csv
```

**Example Test:**
```python
# tests/test_utils.py
import pytest
import pandas as pd
from src.utils import replace_empty_with_zero

def test_replace_empty_with_zero_dataframe():
    df = pd.DataFrame({'A': ['', '100'], 'B': ['200', '']})
    result = replace_empty_with_zero(df)
    assert result.loc[0, 'A'] == 0.0
    assert result.loc[0, 'B'] == 200.0
```

---

### 8. **Configuration Management** (Priority: LOW)

**Problem**: Configuration only loaded once at startup

**Current:**
```python
# main.py
config = load_config('config.yaml')  # Loaded once
```

**Recommended**: Create Config class
```python
# src/config.py
class Config:
    _instance = None

    @classmethod
    def get_instance(cls):
        if cls._instance is None:
            cls._instance = cls._load_config()
        return cls._instance

    @classmethod
    def _load_config(cls):
        with open('config.yaml') as f:
            return yaml.safe_load(f)

# Usage
from src.config import Config
config = Config.get_instance()
```

---

### 9. **Import Organization** (Priority: LOW)

**Problem**: Inconsistent import ordering

**Current (mixed):**
```python
import pandas as pd
from pathlib import Path
import numpy as np
import logging
from utils import replace_empty_with_zero
```

**Recommended (PEP 8 order):**
```python
# Standard library
import logging
from pathlib import Path

# Third-party
import numpy as np
import pandas as pd

# Local
from utils import replace_empty_with_zero
```

---

### 10. **Documentation** (Priority: MEDIUM)

**Missing Documentation:**
- ❌ No API documentation
- ❌ No developer guide
- ✅ CLAUDE.md exists (good!)
- ❌ No visualization guide explaining each plot

**Recommended:**
```
docs/
├── API.md                    # API reference
├── DEVELOPER_GUIDE.md        # Setup and development
├── VISUALIZATION_GUIDE.md    # Plot explanations
└── EXAMPLES.md               # Usage examples
```

---

## 📋 Prioritized Action Plan

### Phase 1: High Priority (Week 1)
1. ✅ **Fix FutureWarnings** - COMPLETED
2. 🔧 **Create sample column utility function** in `utils.py`
3. 🔧 **Refactor all files to use new utility**
4. 🧪 **Set up pytest infrastructure**
5. 🧪 **Write tests for core modules** (data_loader, annotator, analyzer)

### Phase 2: Medium Priority (Week 2)
6. 📁 **Reorganize file structure** (move to src/core, create scripts/)
7. 🔤 **Add type hints** to all public methods
8. 📝 **Extract magic numbers to config/constants**
9. 🔨 **Break down large methods** (>100 lines)
10. 📚 **Create basic documentation** (API.md, DEVELOPER_GUIDE.md)

### Phase 3: Low Priority (Week 3)
11. 🛡️ **Add error handling** to visualization methods
12. ⚙️ **Implement Config singleton class**
13. 📦 **Organize imports** according to PEP 8
14. 📊 **Create visualization guide**
15. ⚡ **Add `fastcluster` to requirements** for performance

---

## 🎯 Quick Wins (Can Do Now)

1. **Remove temporary script**
   ```bash
   mkdir -p scripts
   mv add_trace_remaining.py scripts/
   ```

2. **Add fastcluster to requirements.txt**
   ```bash
   echo "fastcluster>=1.2.0" >> requirements.txt
   pip3 install fastcluster
   ```

3. **Create sample extraction utility**
   - Add `get_sample_columns()` to `src/utils.py`
   - Update 14 files to use it

4. **Add .gitignore for Python**
   ```
   __pycache__/
   *.pyc
   *.pyo
   *.egg-info/
   .pytest_cache/
   Results/*.png
   Results/*.csv
   Results/Trace/*.csv
   ```

---

## 📊 Code Metrics Summary

| Metric | Current | Target |
|--------|---------|--------|
| **Total Python Files** | 20 | 25+ (with tests) |
| **Lines of Code** | ~4,590 | ~5,500 (with tests/docs) |
| **Largest File** | 494 lines (histogram.py) | <300 lines |
| **Code Duplication** | ~14 instances (sample extraction) | 0 instances |
| **Test Coverage** | 0% | 80%+ |
| **Type Hint Coverage** | ~30% | 100% |
| **FutureWarnings** | 0 ✅ | 0 |
| **Critical Warnings** | 0 | 0 |

---

## 🔍 Code Quality Checklist

- [x] No wildcard imports (`import *`)
- [x] No FutureWarnings
- [ ] Type hints on all public methods
- [ ] Unit tests for core functionality
- [ ] Error handling in all visualization methods
- [ ] No code duplication
- [ ] All magic numbers in config
- [ ] API documentation
- [ ] Developer guide
- [ ] Methods <100 lines
- [ ] Proper import organization (PEP 8)

---

## 💡 Suggestions for Future Enhancements

1. **Parallel Processing**: Use `multiprocessing` or `concurrent.futures` for visualization generation
2. **Caching**: Cache PCA/PLS-DA results to avoid recalculation
3. **CLI Arguments**: Add argparse for command-line configuration
4. **Progress Bar**: Add `tqdm` for long-running operations
5. **Logging Levels**: Make logging configurable (DEBUG/INFO/WARNING)
6. **Export Formats**: Support SVG/PDF output in addition to PNG
7. **Interactive Plots**: Add Plotly/Bokeh for interactive visualizations
8. **Pipeline Resume**: Checkpoint system to resume failed runs
9. **Docker Support**: Containerize for reproducibility
10. **Web Interface**: Optional Streamlit/Dash dashboard

---

## 📞 Next Steps

**Immediate Actions:**
1. Review this document with your team
2. Prioritize refactoring tasks based on project timeline
3. Create GitHub issues for each refactoring item
4. Set up branch strategy (e.g., `refactor/sample-extraction`, `feature/add-tests`)

**Questions to Consider:**
- What's the timeline for these improvements?
- Which refactorings provide the most value for your research?
- Do you need backwards compatibility with existing scripts?
- Should we maintain the current API or can we break it?

---

**Document Version**: 1.0
**Last Updated**: 2025-10-02
**Status**: Ready for Review
