"""
PDB complex Analysis Toolkit

A comprehensive Python package for analyzing protein structures,
with focus on antibody-antigen complexes.
"""

__version__ = "0.1.0"
__author__ = "Faisal B. Ashraf"
__email__ = "fashr003@ucr.edu"

# Core functionality
from .pdbparser import parse_pdb
from .models import Complex, Chain, Residue, Atom
from .selection import select, selection_residues

# Analysis modules
from .contacts import analyze_contacts, ContactParameters, ContactType
from .interface import compute_interface
from .distances import residue_distance_matrix, pairwise_min_distance
from .geometry import analyze_geometry
from .sasa_analysis import compare_bound_unbound_sasa

# Feature extraction
from .residue_features import (
    extract_residue_features,
    features_to_dataframe,
    filter_interface_residues,
    identify_hotspots,
    get_feature_summary,
    export_features_csv,
    export_features_json,
)

from .single_chain_features import (
    extract_single_chain_features,
    features_to_dataframe as single_chain_features_to_dataframe,
    export_features_csv as export_single_chain_csv,
    export_features_json as export_single_chain_json,
    get_feature_summary as get_single_chain_summary,
)

# Multi-complex analysis
from .multicomplex.multi_analysis import MultiComplexAnalyzer, AnalysisConfig

# Define what gets imported with "from abagpdb import *"
__all__ = [
    # Version
    '__version__',

    # Core
    'parse_pdb',
    'Complex',
    'Chain',
    'Residue',
    'Atom',
    'select',
    'selection_residues',

    # Analysis
    'analyze_contacts',
    'ContactParameters',
    'ContactType',
    'compute_interface',
    'residue_distance_matrix',
    'pairwise_min_distance',
    'analyze_geometry',
    'compare_bound_unbound_sasa',

    # Feature extraction
    'extract_residue_features',
    'features_to_dataframe',
    'filter_interface_residues',
    'identify_hotspots',
    'get_feature_summary',
    'export_features_csv',
    'export_features_json',

    # Single chain
    'extract_single_chain_features',
    'single_chain_features_to_dataframe',
    'export_single_chain_csv',
    'export_single_chain_json',
    'get_single_chain_summary',

    # Multi-complex
    'MultiComplexAnalyzer',
    'AnalysisConfig',
]
