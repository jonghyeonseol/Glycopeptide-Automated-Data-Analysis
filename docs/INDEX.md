# Documentation Index

This directory contains all project documentation organized by category.

---

## 📁 Directory Structure

```
docs/
├── INDEX.md (this file)          # Documentation index
│
├── Current Documentation          # Active project docs
│   ├── CHANGELOG.md              # Version history
│   ├── PROJECT_STATUS.md         # Current project status (v2.1.0)
│   ├── RESTRUCTURE_SUMMARY.md    # v2.1 restructure summary
│   └── DATA_INTEGRITY_REPORT.md  # Data integrity verification
│
├── Guides                        # User/developer guides
│   ├── README.md                 # Docs overview
│   ├── glycan-sorting-guide.md
│   ├── normalization.md
│   ├── trace-data-reference.md
│   ├── verification-guide.md
│   ├── visualization-guide.md
│   ├── visualization-updates.md
│   └── visualization-enhancements.md
│
├── redesign/                     # v3.0 Redesign documentation
│   ├── PIPELINE_REDESIGN_COMPLETE.md    # Complete summary
│   ├── PIPELINE_REDESIGN_PHASE1.md      # Phase 1: Core Pipeline
│   ├── PIPELINE_REDESIGN_PHASE2.md      # Phase 2: Analysis Split
│   └── PIPELINE_REDESIGN_PHASE3-5.md    # Phase 3-5: Viz/Report/Config
│
└── archive/                      # Historical development reports
    ├── CHANGELOG_old.md
    ├── DETECTION_FILTER_IMPACT.md
    ├── IMPLEMENTATION_COMPLETE.md
    ├── PERFORMANCE_NOTES.md
    ├── REFACTORING_REPORT.md
    ├── REFACTORING_SUMMARY.md
    ├── RELIABILITY_ANALYSIS.md
    ├── SCIENTIFIC_VALIDITY_REVIEW.md
    ├── VISUALIZATION_APPEARANCE_REVIEW.md
    └── VISUALIZATION_RELIABILITY_ISSUES.md
```

---

## 📚 Quick Links

### Getting Started
- [README](../README.md) - Main project readme
- [Architecture](../ARCHITECTURE.md) - System architecture
- [Visualization Guide](visualization-guide.md) - All plot types explained

### Current Status (v2.1.0)
- [Project Status](PROJECT_STATUS.md) - Current version status
- [Data Integrity Report](DATA_INTEGRITY_REPORT.md) - Comprehensive integrity review
- [Restructure Summary](RESTRUCTURE_SUMMARY.md) - v2.1 changes

### v3.0 Redesign
- [Redesign Complete](redesign/PIPELINE_REDESIGN_COMPLETE.md) - Full summary
- [Phase 1: Core Pipeline](redesign/PIPELINE_REDESIGN_PHASE1.md)
- [Phase 2: Analysis Split](redesign/PIPELINE_REDESIGN_PHASE2.md)
- [Phase 3-5: Complete Architecture](redesign/PIPELINE_REDESIGN_PHASE3-5.md)

### Development History
- [Changelog](CHANGELOG.md) - Version history
- [Archive](archive/) - Historical development reports

---

## 🔍 Find Documentation By Topic

### Data Processing
- [Normalization Guide](normalization.md)
- [Data Preparation](../ARCHITECTURE.md#data-preparation)
- [Filtering](PROJECT_STATUS.md#filtering)

### Visualization
- [Visualization Guide](visualization-guide.md)
- [Visualization Updates](visualization-updates.md)
- [Glycan Sorting](glycan-sorting-guide.md)

### Data Integrity
- [Data Integrity Report](DATA_INTEGRITY_REPORT.md)
- [Verification Guide](verification-guide.md)
- [Trace Data Reference](trace-data-reference.md)

### Architecture
- [Main Architecture Doc](../ARCHITECTURE.md)
- [v3.0 Redesign](redesign/PIPELINE_REDESIGN_COMPLETE.md)
- [Phase Summaries](redesign/)

---

## 📝 Documentation Types

### Active Documentation
Current, maintained documentation for the latest version.

**Location**: `docs/` (root level)

### Guides
User and developer guides for using the pipeline.

**Location**: `docs/` (mixed with active docs)

### Redesign Documentation
Complete v3.0 architectural redesign documentation.

**Location**: `docs/redesign/`

### Historical Archive
Old reports and analysis from development process.

**Location**: `docs/archive/`

---

*Last Updated*: 2025-10-06
