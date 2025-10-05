# Glycopeptide Comparison Heatmap - Visualization Guide

## 🎯 Understanding the Visualization

### What are Circles and Squares?

**The visualization compares TWO groups side-by-side:**

| Symbol | Meaning | Position | Color |
|--------|---------|----------|-------|
| **○ Circle** | **Cancer group** | Left dot | Red tint |
| **□ Square** | **Normal group** | Right dot | Blue tint |

### How to Read Each Dot Pair

For **EVERY peptide-glycan combination**, you see **TWO dots**:

```
     ○ □           ○ = Cancer (averaged across C1-C24)
    Cancer Normal  □ = Normal (averaged across N1-N23)
```

**Example Interpretation:**

```
○ ■■              Dark Cancer, Light Normal → Higher in Cancer
  ■■ □            Light Cancer, Dark Normal → Higher in Normal
○○ □□             Both dark → Present in both groups
○ .               Only Cancer → Only detected in Cancer
. □               Only Normal → Only detected in Normal
```

**Dot Darkness** = Intensity (darker = higher abundance)

---

## 📐 Three-Panel Layout

### Complete Layout (UPDATED)

```
┌─────────────────────────────────────────────────────────┐
│  TOP PANEL                                              │
│  Aggregated Intensity Line Plots                       │
│  • Red line = Cancer samples                           │
│  • Blue line = Normal samples                          │
├─────────────────────────────────────────────────────────┤
│  MIDDLE PANEL - Gradient Colored Bar                   │
│  [HM (Green)] [F (Red)] [S (Pink)] [SF (Orange)] [C/H (Blue)] │
├─────────────────────────────────────────────────────────┤
│                                                         │
│  MAIN PANEL - Heatmap-like Rectangular Box             │
│  • Y-axis: Peptides (sorted by VIP score)             │
│  • Dots: ○ (Cancer left) □ (Normal right)             │
│                                                         │
│  H  H  H  H  H  H  H  H  H  H  H  ...                 │  ← Glycan composition
│  (  (  (  (  (  (  (  (  (  (  (                      │    labels (rotated 90°)
│  5  5  6  7  5  5  5  6  6  7  5                       │    vertical text
│  )  )  )  )  )  )  )  )  )  )  )                       │
│  N  N  N  N  N  N  N  N  N  N  N                       │
│  (  (  (  (  (  (  (  (  (  (  (                       │
│  4  4  5  6  4  4  4  5  5  6  4                       │
│  )  )  )  )  )  )  )  )  )  )  )                       │
│  A  A  A  A  F  F  F  F  F  F  A                       │
│  (  (  (  (  (  (  (  (  (  (  (                       │
│  1  2  2  4  1  1  1  2  2  4  1                       │
│  )  )  )  )  )  )  )  )  )  )  )                       │
└─────────────────────────────────────────────────────────┘
       Glycan Composition (detailed labels below)
```

**Four-Level Information Hierarchy:**

1. **Top Panel**: Overall group intensity trends
2. **Colored Bar**: Glycan type grouping (HM, F, S, SF, C/H)
3. **Main Heatmap**: Detailed peptide-glycan comparisons
4. **Bottom Labels**: Specific glycan compositions (vertical text)

**Why this layout?**
- The colored bar now sits **between** the line graph and heatmap
- This provides **immediate visual context** before looking at detailed dots
- You can quickly see which columns belong to which glycan type

---

## 🎨 Color Coding

### Glycan Type Colors (on colored bar)

| Type | Full Name | Color | Hex | Biology |
|------|-----------|-------|-----|---------|
| **HM** | High-mannose | 🟢 Green | #00CC00 | Early glycosylation, ER-associated |
| **F** | Fucosylated | 🔴 Red | #FF0000 | Core fucosylation, immune signaling |
| **S** | Sialylated | 🩷 Pink | #FF69B4 | Terminal sialylation, cell recognition |
| **SF** | Sialofucosylated | 🟠 Orange | #FFA500 | Both modifications, complex |
| **C/H** | Complex/Hybrid | 🔵 Blue | #0000FF | Processed glycans, complex branching |

### Dot Colors

Each dot inherits the color of its glycan type:
- A dot in the **HM region** = Green
- A dot in the **F region** = Red
- etc.

**Transparency** encodes intensity:
- **Darker/Opaque** = High intensity
- **Lighter/Transparent** = Low intensity

---

## 📊 Step-by-Step Reading Guide

### Step 1: Check the Top Panel
**Question**: Which group has higher overall intensity?

Look at the line plots:
- **Red line higher** → Cancer group has more of those glycans
- **Blue line higher** → Normal group has more

### Step 2: Identify Glycan Types (Middle Panel)
**Question**: Which glycan types are present?

The colored bar shows:
- Width of each color region = number of glycans of that type
- Labels (HM, F, S, SF, C/H) identify each region

### Step 3: Compare Specific Glycopeptides (Main Panel)
**Question**: Which peptide-glycan combinations differ between groups?

Look for patterns:
- **Dark ○, Light □** → Higher in Cancer
- **Light ○, Dark □** → Higher in Normal
- **○ only, no □** → Cancer-specific
- **□ only, no ○** → Normal-specific

### Step 4: Find Important Biomarkers
**Question**: Which peptides are most discriminative?

Peptides are sorted by **VIP score** (top = highest):
- **Top rows** = Most important for distinguishing Cancer vs Normal
- **Bottom rows** = Less discriminative

---

## 🔍 Example Analyses

### Finding Cancer-Specific Markers

**Look for:**
1. **Top rows** (high VIP score)
2. **Dark ○ circles** (high Cancer intensity)
3. **Light or missing □ squares** (low/absent Normal intensity)

**Interpretation**: These glycopeptides are **upregulated in cancer** and have **high biomarker potential**.

### Finding Type-Specific Patterns

**Look for:**
1. **Entire columns** within a color region (e.g., all HM glycans)
2. Compare **left dots (○)** vs **right dots (□)** in that region

**Interpretation**: Does one group show higher levels of a specific glycan type?

Example:
- If **all HM dots** show **dark ○, light □** → High-mannose glycans are **upregulated in Cancer**
- If **all S dots** show **light ○, dark □** → Sialylated glycans are **upregulated in Normal**

### Finding Peptide-Specific Patterns

**Look for:**
1. **Entire rows** (one peptide with different glycan modifications)
2. Compare patterns **across glycan types** for the same peptide

**Interpretation**: How does glycosylation of this specific peptide differ?

Example:
- Peptide X shows **dark ○** in F region, **dark □** in S region
- → This peptide is **fucosylated in Cancer**, **sialylated in Normal**

---

## 💡 Visual Interpretation Tips

### Quick Patterns to Spot

| Visual Pattern | Meaning |
|----------------|---------|
| **Vertical line of dark ○** | Cancer-enriched glycan type |
| **Vertical line of dark □** | Normal-enriched glycan type |
| **Horizontal line of dark ○** | Cancer-enriched peptide (all glycoforms) |
| **Alternating ○ and □** | Mixed regulation across glycan types |
| **Sparse dots** | Low-abundance glycopeptides |
| **Dense dots** | High-abundance glycopeptides |

### Common Biological Interpretations

**High-mannose (HM) enriched in Cancer**:
- Suggests incomplete glycosylation processing
- May indicate ER stress or rapid cell proliferation

**Sialylated (S/SF) enriched in Cancer**:
- Suggests immune evasion mechanisms
- May indicate cancer immune microenvironment changes

**Fucosylated (F/SF) differences**:
- Core fucosylation often cancer-associated
- May affect cell adhesion and signaling

---

## 📈 Statistical Considerations

### VIP Score (Y-axis Sorting)

**VIP Score** = Variable Importance in Projection (from PLS-DA)

- **VIP > 1.0**: Important for group discrimination
- **VIP > 1.5**: Highly important
- **VIP > 2.0**: Extremely important (potential biomarkers)

Top peptides in the heatmap typically have **VIP > 2.0**.

### Aggregated Data

**Important**: The dots represent **averaged data**:
- **Cancer (○)**: Mean across C1-C24 (24 samples)
- **Normal (□)**: Mean across N1-N23 (23 samples)

**Benefits**:
- Reduces noise from individual sample variability
- Highlights consistent group differences
- More robust for biomarker discovery

**Caveat**:
- Individual sample heterogeneity is hidden
- Outlier samples may be masked

---

## 🎓 Advanced Reading

### Multi-Dimensional Comparison

You can compare **three dimensions simultaneously**:

1. **Peptide dimension** (Y-axis): Which protein/sequence?
2. **Glycan dimension** (X-axis): Which modification type?
3. **Group dimension** (○ vs □): Cancer vs Normal?

This allows questions like:
- "Is peptide A more sialylated in cancer?"
- "Are high-mannose glycans generally higher in normal?"
- "Which peptide shows the strongest cancer-specific fucosylation?"

### Combining with Other Data

**Cross-reference with:**
- **VIP scores table**: Get exact VIP values
- **Statistical tests**: Confirm significance of visual patterns
- **Protein information**: Understand biological context
- **Pathway analysis**: Link to functional consequences

---

## 🛠️ Customization Tips

### Adjusting for Your Data

Edit `config.yaml`:

```yaml
glycopeptide_comparison:
  max_peptides: 50          # More peptides = more detail
  max_glycans_per_type: 15  # More glycans = wider plot
  figsize: [24, 16]         # Larger = more readable
```

**Trade-offs:**
- **More features** = More comprehensive but crowded
- **Fewer features** = Cleaner but may miss patterns

### Recommended Settings by Goal

| Goal | max_peptides | max_glycans_per_type |
|------|--------------|---------------------|
| **Biomarker discovery** | 20-30 | 8-10 |
| **Comprehensive overview** | 50 | 15 |
| **Publication figure** | 30 | 10 |
| **Detailed analysis** | 50+ | 15+ |

---

## ❓ FAQ

### Q: Why are some combinations missing dots?

**A**: The glycopeptide was not detected (or below detection threshold) in that group.

### Q: What if both dots have similar darkness?

**A**: Similar intensity between groups → Not discriminative for this comparison.

### Q: Why is the colored bar important?

**A**: It provides **instant visual grouping** - you can immediately see where different glycan types are, without reading individual labels.

### Q: How do I find the best biomarkers?

**A**: Look at **top rows** (high VIP) with **large intensity difference** between ○ and □.

### Q: Can I trust dots that appear only in one group?

**A**: Yes, but verify with statistics. These may be **group-specific glycopeptides** (strong biomarker candidates).

---

## 📚 Related Files

- **Data Source**: `Results/glycopeptide_comparison_heatmap_data.csv`
- **VIP Scores**: `Results/vip_scores_all.csv`
- **Configuration**: `config.yaml`
- **Code**: `src/plots/glycopeptide_comparison_heatmap.py`

---

## 🎯 Summary

**Key Takeaways:**

1. **○ = Cancer (left), □ = Normal (right)** - compare side-by-side
2. **Colored bar (middle)** - shows glycan type groupings
3. **Darker dots** - higher intensity
4. **Top rows** - most important peptides (VIP score)
5. **Grouped columns** - glycans organized by biological type

**This visualization enables:**
- Quick group comparison (Cancer vs Normal)
- Glycan type pattern recognition
- Biomarker candidate identification
- Multi-dimensional glycoproteomics analysis

---

**For more details, see `VISUALIZATION_IMPROVEMENTS.md`**
