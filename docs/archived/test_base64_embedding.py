#!/usr/bin/env python3
"""Test base64 PNG embedding functionality"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from src.interactive_dashboard import InteractiveDashboard

# Test base64 encoding
dashboard = InteractiveDashboard(output_dir='Results')

# Find a test PNG file
test_png = Path('Results/pca_plot.png')

if test_png.exists():
    print(f"Testing base64 encoding of {test_png}...")
    base64_str = dashboard._png_to_base64(test_png)

    if base64_str and base64_str.startswith('data:image/png;base64,'):
        encoded_length = len(base64_str)
        print("✓ Successfully encoded PNG to base64")
        print(f"  Original file size: {test_png.stat().st_size:,} bytes")
        print(f"  Base64 string length: {encoded_length:,} characters")
        print(f"  Sample (first 100 chars): {base64_str[:100]}...")
        print("\n✓ Base64 embedding functionality working correctly!")
    else:
        print("✗ Failed to encode PNG")
        sys.exit(1)
else:
    print(f"✗ Test PNG not found at {test_png}")
    print("  Run the main pipeline first to generate visualizations")
    sys.exit(1)
