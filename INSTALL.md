# Installation Guide for abagpdb

This guide provides detailed instructions for installing the abagpdb package.

## Quick Install

### Option 1: Install from Source (Recommended)

```bash
# Clone or download the repository
git clone https://github.com/yourusername/abagpdb.git
cd abagpdb

# Install the package
pip install -e .
```

The `-e` flag installs in "editable" mode, allowing you to modify the code without reinstalling.

### Option 2: Install with All Optional Dependencies

```bash
pip install -e ".[all]"
```

This includes visualization tools (plotly, seaborn) and development tools.

### Option 3: Install Specific Extras

```bash
# Just visualization tools
pip install -e ".[viz]"

# Just development tools
pip install -e ".[dev]"
```

## Requirements

### Python Version
- **Python 3.8 or higher** is required
- Tested on Python 3.8, 3.9, 3.10, and 3.11

### Core Dependencies
These are installed automatically:
- `numpy >= 1.20.0`
- `pandas >= 1.3.0`
- `matplotlib >= 3.4.0`

### Optional Dependencies

#### For SASA Calculations (Highly Recommended)
```bash
pip install freesasa
```

FreeSASA provides fast and accurate solvent accessible surface area calculations. If not installed, the package will use a fallback method.

#### For Interactive Visualizations
```bash
pip install plotly seaborn
```

#### For Development
```bash
pip install pytest pytest-cov black flake8 mypy
```

## Verify Installation

After installation, verify it works:

```python
import abagpdb
print(abagpdb.__version__)

# Test basic functionality
from abagpdb import parse_pdb
cx = parse_pdb("examples/5GGS_wt.pdb")  # Assuming you have example data
print(f"Loaded structure with {len(cx.chains)} chains")
```

## Installation in Different Environments

### Virtual Environment (Recommended)

```bash
# Create virtual environment
python -m venv abagpdb-env

# Activate it
# On Linux/Mac:
source abagpdb-env/bin/activate
# On Windows:
abagpdb-env\Scripts\activate

# Install the package
cd abagpdb-package
pip install -e .
```

### Conda Environment

```bash
# Create conda environment
conda create -n abagpdb python=3.10
conda activate abagpdb

# Install the package
cd abagpdb-package
pip install -e .
```

### Jupyter Notebook

If you want to use the package in Jupyter notebooks:

```bash
# Install jupyter
pip install jupyter notebook

# Make the package available in Jupyter
python -m ipykernel install --user --name=abagpdb

# Start Jupyter
jupyter notebook examples/
```

## Troubleshooting

### Issue: "No module named 'abagpdb'"

**Solution**: Make sure you installed the package and are in the correct Python environment:
```bash
pip list | grep abagpdb
```

### Issue: FreeSASA not installing

**Solution**: FreeSASA requires compilation. Install build tools:

**Linux (Ubuntu/Debian)**:
```bash
sudo apt-get install build-essential python3-dev
pip install freesasa
```

**macOS**:
```bash
xcode-select --install
pip install freesasa
```

**Windows**: FreeSASA may be difficult to install. The package will work without it using fallback methods.

### Issue: Import errors

**Solution**: Reinstall with dependencies:
```bash
pip install -e . --force-reinstall
```

### Issue: Example notebooks not working

**Solution**: Make sure you're in the package root directory:
```bash
cd abagpdb-package
jupyter notebook examples/
```

## Development Installation

For developers who want to contribute:

```bash
# Clone the repository
git clone https://github.com/yourusername/abagpdb.git
cd abagpdb-package

# Install in development mode with all extras
pip install -e ".[dev,all]"

# Run tests
pytest

# Format code
black abagpdb/

# Lint code
flake8 abagpdb/
```

## Uninstallation

To remove the package:

```bash
pip uninstall abagpdb
```

## Getting Help

If you encounter installation issues:

1. Check Python version: `python --version`
2. Check pip version: `pip --version`
3. Update pip: `pip install --upgrade pip`
4. Open an issue: https://github.com/yourusername/abagpdb/issues

---

**After successful installation, check out the [examples/](examples/) directory for comprehensive tutorials!**
