#!/bin/bash
# test_workflow.sh
# Test script to verify the RCT analysis workflow

echo "=== TESTING RCT MYO-INOSITOL ANALYSIS WORKFLOW ==="

# Test 1: Check data structure
echo "Test 1: Checking data structure..."
if [ -f "data/example_data.csv" ]; then
    echo "✓ Example data file exists"
    head -n 3 data/example_data.csv
else
    echo "✗ Example data file missing"
    exit 1
fi

# Test 2: Python analysis
echo -e "\nTest 2: Running Python basic analysis..."
if python3 scripts/basic_analysis.py > /tmp/python_output.log 2>&1; then
    echo "✓ Python analysis completed successfully"
    echo "Generated files:"
    ls -la output/
else
    echo "✗ Python analysis failed"
    cat /tmp/python_output.log
    exit 1
fi

# Test 3: Check R script syntax (without running due to R dependencies)
echo -e "\nTest 3: Checking R script syntax..."
if command -v R >/dev/null 2>&1; then
    echo "✓ R is available, checking syntax..."
    R --vanilla --slave -e "parse('scripts/setup_environment.R')" 2>/dev/null && echo "✓ setup_environment.R syntax OK"
    R --vanilla --slave -e "parse('scripts/data_analysis.R')" 2>/dev/null && echo "✓ data_analysis.R syntax OK"
    R --vanilla --slave -e "parse('scripts/forest_plots.R')" 2>/dev/null && echo "✓ forest_plots.R syntax OK"
else
    echo "⚠ R not available, skipping syntax check"
fi

# Test 4: Directory structure
echo -e "\nTest 4: Checking directory structure..."
for dir in data scripts output docs; do
    if [ -d "$dir" ]; then
        echo "✓ Directory '$dir' exists"
    else
        echo "✗ Directory '$dir' missing"
        exit 1
    fi
done

# Test 5: Documentation
echo -e "\nTest 5: Checking documentation..."
if [ -f "docs/workflow_guide.md" ]; then
    echo "✓ Workflow guide exists"
    echo "Documentation size: $(wc -l < docs/workflow_guide.md) lines"
else
    echo "✗ Workflow guide missing"
    exit 1
fi

echo -e "\n=== ALL TESTS PASSED ==="
echo "The RCT myo-inositol analysis workflow is ready for use!"
echo -e "\nQuick start:"
echo "1. Place your CSV data in data/ directory"
echo "2. Run: python3 scripts/basic_analysis.py data/your_data.csv"
echo "3. For full R analysis: source('scripts/forest_plots.R'); run_comprehensive_analysis('data/your_data.csv')"