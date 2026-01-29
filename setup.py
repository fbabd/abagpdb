"""
Setup configuration for abagpdb package
"""
from setuptools import setup, find_packages
from pathlib import Path

# Read the README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding='utf-8')

# Read requirements
requirements = []
requirements_path = this_directory / "requirements.txt"
if requirements_path.exists():
    requirements = requirements_path.read_text().strip().split('\n')
    requirements = [req.strip() for req in requirements if req.strip() and not req.startswith('#')]

setup(
    name="abagpdb",
    version="0.1.0",
    author="Faisal B Ashraf",  
    author_email="faisal.b.ashraf@gmail.com",   
    description="Comprehensive toolkit for protein-protein complex structure analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/fbabd/abagpdb",   
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        'dev': [
            'pytest>=7.0',
            'pytest-cov>=3.0',
            'black>=22.0',
            'flake8>=4.0',
            'mypy>=0.950',
        ],
        'viz': [
            'plotly>=5.0',
            'seaborn>=0.11',
        ],
        'all': [
            'pytest>=7.0',
            'pytest-cov>=3.0',
            'plotly>=5.0',
            'seaborn>=0.11',
        ],
    },
    include_package_data=True,
    zip_safe=False,
    keywords=[
        'bioinformatics',
        'protein structure',
        'antibody',
        'antigen',
        'PDB',
        'molecular interactions',
        'structural biology',
        'protein analysis',
    ],
    project_urls={
        'Documentation': 'https://github.com/fbabd/abagpdb#readme',
        'Source': 'https://github.com/fbabd/abagpdb',
        'Bug Reports': 'https://github.com/fbabd/abagpdb/issues',
    },
)
