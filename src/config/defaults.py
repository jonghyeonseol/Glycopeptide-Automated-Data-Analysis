"""
Default Configuration
Provides default configuration values
"""

DEFAULT_CONFIG = {
    'paths': {
        'dataset_dir': 'Dataset',
        'results_dir': 'Results',
        'output_file': 'integrated.csv'
    },

    'processing': {
        'required_columns': [
            'Peptide',
            'GlycanComposition',
            'IsotopeArea',
            'Proteins'
        ],
        'qc_filters': {
            'min_isotope_area': 0
        }
    },

    'annotation': {
        'sialylation_marker': 'A',
        'fucosylation_marker': 'F'
    },

    'analysis': {
        'pca': {
            'n_components': 2,
            'log_transform': True
        },
        'detection_filter': {
            'min_detection_pct': 0.30,
            'min_samples': 5
        },
        'missing_data_handling': {
            'method': 'skipna'
        },
        'statistical_tests': {
            'alpha': 0.05,
            'fdr_correction': True,
            'method': 'mannwhitneyu'
        }
    },

    'visualization': {
        'figsize': {
            'pca': [10, 8],
            'boxplot': [12, 6],
            'heatmap': [14, 10]
        },
        'dpi': 300,
        'colors': {
            'Non': '#CCCCCC',
            'Sialylated': '#E74C3C',
            'Fucosylated': '#3498DB',
            'Both': '#9B59B6'
        },
        'glycopeptide_comparison': {
            'enabled': True,
            'max_peptides': 20,
            'max_glycans_per_type': 15,
            'figsize': [24, 16]
        }
    }
}
