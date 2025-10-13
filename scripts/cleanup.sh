#!/bin/bash
# pGlyco Auto Combine - Automated Cleanup Script
# Removes temporary files and archives logs
# Version: 3.1.0

set -e  # Exit on error

# Color codes for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'  # No Color

echo "=================================="
echo "pGlyco Auto Combine Cleanup Script"
echo "=================================="
echo ""

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_ROOT="$( cd "$SCRIPT_DIR/.." && pwd )"

cd "$REPO_ROOT"

# 1. Remove Python cache files
echo -n "Cleaning Python cache files... "
find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
find . -type f \( -name "*.pyc" -o -name "*.pyo" \) -delete 2>/dev/null || true
echo -e "${GREEN}✓${NC}"

# 2. Remove system files
echo -n "Removing system files (.DS_Store)... "
find . -name ".DS_Store" -delete 2>/dev/null || true
echo -e "${GREEN}✓${NC}"

# 3. Archive audit logs from Results/ to Logs/audit_logs/
echo -n "Archiving audit logs... "
mkdir -p Logs/audit_logs
if ls Results/audit_log_*.jsonl 1> /dev/null 2>&1; then
    mv Results/audit_log_*.jsonl Logs/audit_logs/ 2>/dev/null || true
    echo -e "${GREEN}✓${NC} (moved to Logs/audit_logs/)"
else
    echo -e "${YELLOW}⊘${NC} (no audit logs found)"
fi

# 4. Summary
echo ""
echo "=================================="
echo -e "${GREEN}Cleanup complete!${NC}"
echo "=================================="
echo ""
echo "Cleaned:"
echo "  ✓ Python cache (__pycache__, *.pyc, *.pyo)"
echo "  ✓ System files (.DS_Store)"
echo "  ✓ Audit logs (archived to Logs/audit_logs/)"
echo ""
echo "Tip: Run this script periodically to keep your repository clean."
echo ""
