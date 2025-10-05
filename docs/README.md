# pGlyco Auto Combine - Documentation

This directory contains comprehensive documentation for the pGlyco Auto Combine pipeline.

## 📚 Documentation Index

### Getting Started
- [Main README](../README.md) - Project overview, installation, and quick start
- [Architecture](../ARCHITECTURE.md) - System design and code organization
- [Normalization Pipeline](normalization.md) - TIC normalization and data preprocessing

### Visualization
- [Visualization Guide](visualization-guide.md) - Complete guide to all visualizations
- [Visualization Updates](visualization-updates.md) - Latest symbol-based design
- [Visualization Enhancements](visualization-enhancements.md) - User-friendly improvements
- [Glycan Sorting Guide](glycan-sorting-guide.md) - Numeric sorting algorithm

### Data Verification
- [Verification Guide](verification-guide.md) - Manual verification walkthrough (Excel-friendly)
- [Trace Data Reference](trace-data-reference.md) - Technical reference for trace data files

### Change History
- [CHANGELOG](CHANGELOG.md) - All recent improvements and fixes

## 🎯 Quick Navigation

**For Users:**
- First time? Start with [Main README](../README.md)
- Want to verify results? See [Verification Guide](verification-guide.md)
- Understanding visualizations? Check [Visualization Guide](visualization-guide.md)

**For Developers:**
- Pipeline architecture: [Main README](../README.md#architecture)
- Data processing: [Normalization](normalization.md)
- Trace data format: [Trace Data Reference](trace-data-reference.md)

**For AI Assistants:**
- See [CLAUDE.md](../CLAUDE.md) for project-specific instructions

## 📁 File Organization

```
pGlyco_auto_combine/
├── README.md                    # Main project documentation
├── CLAUDE.md                    # AI assistant instructions
├── docs/                        # All user documentation (you are here)
│   ├── README.md               # This file (documentation index)
│   ├── visualization-guide.md  # Complete visualization reference
│   ├── verification-guide.md   # Manual verification guide
│   ├── trace-data-reference.md # Technical trace data details
│   └── normalization.md        # Normalization pipeline docs
├── src/                        # Source code
├── Dataset/                    # Input CSV files
└── Results/                    # Output files and visualizations
    └── Trace/                  # Trace data for verification
```

## 🔍 Finding What You Need

| I want to... | Go to... |
|--------------|----------|
| Install and run the pipeline | [Main README](../README.md) |
| Understand PCA/heatmap/boxplot outputs | [Visualization Guide](visualization-guide.md) |
| Verify a specific dot in the heatmap | [Verification Guide](verification-guide.md) |
| Understand trace data columns | [Trace Data Reference](trace-data-reference.md) |
| Learn about TIC normalization | [Normalization](normalization.md) |
| Modify the pipeline (developers) | [CLAUDE.md](../CLAUDE.md) + source code |

## 💡 Common Questions

**Q: How do I verify my results?**
A: See [Verification Guide](verification-guide.md) for step-by-step Excel-based verification.

**Q: What do the trace data files contain?**
A: See [Trace Data Reference](trace-data-reference.md) for complete column descriptions.

**Q: How does normalization work?**
A: See [Normalization](normalization.md) for TIC normalization details.

**Q: What visualizations are generated?**
A: See [Visualization Guide](visualization-guide.md) for all plot types and interpretations.

---

**Last Updated:** 2025-10-05
**Version:** 2.0 (Documentation consolidated)
