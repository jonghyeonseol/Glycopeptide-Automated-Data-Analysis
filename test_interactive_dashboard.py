#!/usr/bin/env python3
"""Quick test for interactive dashboard module"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

# Test imports
print("Testing interactive dashboard imports...")
try:
    from src.interactive_dashboard import InteractiveDashboard
    print("✓ InteractiveDashboard imported successfully")
except Exception as e:
    print(f"✗ Failed to import InteractiveDashboard: {e}")
    sys.exit(1)

# Test plotly imports
try:
    import plotly.graph_objects as go
    import plotly.express as px
    print("✓ Plotly imported successfully")
except Exception as e:
    print(f"✗ Failed to import Plotly: {e}")
    sys.exit(1)

# Test instantiation
try:
    dashboard = InteractiveDashboard(output_dir='Results')
    print("✓ Dashboard instantiated successfully")
    print(f"  Output directory: {dashboard.output_dir}")
except Exception as e:
    print(f"✗ Failed to instantiate dashboard: {e}")
    sys.exit(1)

print("\n✓ All tests passed! Interactive dashboard module is ready.")
