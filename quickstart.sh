#!/bin/bash
# Quick Start Commands for RecJf-RPA Primer Design

echo "======================================================================"
echo "  RecJf-RPA lncRNA Primer Design - Quick Start"
echo "======================================================================"
echo ""

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 is not installed. Please install Python 3.8+"
    exit 1
fi

echo "✓ Python version:"
python3 --version
echo ""

# Check if dependencies are installed
echo "📦 Checking dependencies..."
pip list | grep -i "biopython\|pandas\|numpy" > /dev/null
if [ $? -ne 0 ]; then
    echo "⚠️  Some dependencies might be missing."
    echo "Installing dependencies..."
    pip install -r requirements.txt
fi

echo ""
echo "======================================================================"
echo "  Available Commands"
echo "======================================================================"
echo ""
echo "1️⃣  Design primers for HOTAIR:"
echo "    python3 examples.py 1"
echo ""
echo "2️⃣  Batch design (HOTAIR + MALAT1):"
echo "    python3 examples.py 2"
echo ""
echo "3️⃣  Design with validation:"
echo "    python3 examples.py 3"
echo ""
echo "4️⃣  Custom sequence:"
echo "    python3 examples.py 4"
echo ""
echo "5️⃣  All targets:"
echo "    python3 examples.py 5"
echo ""
echo "6️⃣  Comparison:"
echo "    python3 examples.py 6"
echo ""
echo "📊 Full pipeline (recommended):"
echo "    python3 run_pipeline.py"
echo ""
echo "🔍 See available targets:"
echo "    python3 run_pipeline.py --list-targets"
echo ""
echo "======================================================================"
echo ""

# Ask user which example to run
read -p "Enter example number (1-6) or press Enter for Example 1: " example_num
example_num=${example_num:-1}

echo ""
echo "Running Example $example_num..."
echo ""

python3 examples.py $example_num

echo ""
echo "======================================================================"
echo "✓ Done! Check the output files in your current directory."
echo "======================================================================"
