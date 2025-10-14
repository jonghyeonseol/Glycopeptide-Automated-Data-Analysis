#!/usr/bin/env python3
"""Quick test script to regenerate Group → Glycan Type Sankey plot"""

import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent / 'src'))

from src.visualizer import GlycanVisualizer
from src.data_preparation import DataPreparationConfig

# Load filtered data
results_dir = Path('Results')
df = pd.read_csv(results_dir / 'integrated_filtered.csv')
vip_scores = pd.read_csv(results_dir / 'vip_scores_all.csv')

# Create visualizer
visualizer = GlycanVisualizer(
    output_dir=results_dir,
    dpi=300,
    colors={}
)

# Create config
config = DataPreparationConfig(
    min_detection_pct=0.30,
    min_samples=5,
    missing_data_method='skipna'
)

# Generate Sankey plot
print("Generating Group → Glycan Type Sankey diagram...")
visualizer.plot_group_to_glycan_sankey(
    df=df,
    vip_scores=vip_scores,
    config=config,
    log2fc_threshold=1.0,
    fdr_threshold=0.05
)

print("✓ Sankey plot generated successfully!")
