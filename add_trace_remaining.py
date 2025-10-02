#!/usr/bin/env python3
"""
Helper script to add trace data exports to remaining visualization methods
"""

import re
from pathlib import Path

# Define the files and their save patterns
updates = [
    # Boxplot.py - 6 methods
    {
        'file': 'src/plots/boxplot.py',
        'replacements': [
            {
                'old': '        output_file = self.output_dir / \'boxplot_glycan_types.png\'\n        plt.savefig(output_file, dpi=self.dpi, bbox_inches=\'tight\')\n        logger.info(f"Saved boxplot to {output_file}")\n\n        plt.close()',
                'new': '        # Save plot\n        output_file = self.output_dir / \'boxplot_glycan_types.png\'\n        plt.savefig(output_file, dpi=self.dpi, bbox_inches=\'tight\')\n        logger.info(f"Saved boxplot to {output_file}")\n\n        # Save trace data\n        save_trace_data(boxplot_data, self.output_dir, \'boxplot_glycan_types_data.csv\')\n\n        plt.close()'
            },
            {
                'old': '        output_file = self.output_dir / \'boxplot_extended_categories.png\'\n        plt.savefig(output_file, dpi=self.dpi, bbox_inches=\'tight\')\n        logger.info(f"Saved extended boxplot to {output_file}")\n\n        plt.close()',
                'new': '        # Save plot\n        output_file = self.output_dir / \'boxplot_extended_categories.png\'\n        plt.savefig(output_file, dpi=self.dpi, bbox_inches=\'tight\')\n        logger.info(f"Saved extended boxplot to {output_file}")\n\n        # Save trace data\n        save_trace_data(boxplot_data, self.output_dir, \'boxplot_extended_categories_data.csv\')\n\n        plt.close()'
            }
        ]
    },
    # Distribution plot
    {
        'file': 'src/plots/distribution_plot.py',
        'replacements': [
            {
                'old': '        output_file = self.output_dir / \'glycan_type_distribution.png\'\n        plt.savefig(output_file, dpi=self.dpi, bbox_inches=\'tight\')\n        logger.info(f"Saved glycan type distribution to {output_file}")\n\n        plt.close()',
                'new': '        # Save plot\n        output_file = self.output_dir / \'glycan_type_distribution.png\'\n        plt.savefig(output_file, dpi=self.dpi, bbox_inches=\'tight\')\n        logger.info(f"Saved glycan type distribution to {output_file}")\n\n        # Save trace data\n        trace_data = count_df.copy()\n        save_trace_data(trace_data, self.output_dir, \'glycan_type_distribution_data.csv\')\n\n        plt.close()'
            }
        ]
    }
]

# Apply updates
for update in updates:
    file_path = Path(update['file'])
    if file_path.exists():
        content = file_path.read_text()
        for repl in update['replacements']:
            content = content.replace(repl['old'], repl['new'])
        file_path.write_text(content)
        print(f"Updated {file_path}")

print("All trace data exports added successfully!")
