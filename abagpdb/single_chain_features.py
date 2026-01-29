"""
Single Chain Residue Feature Extraction

This module extracts comprehensive features for residues in a single protein chain
or multiple chains treated independently (not as a complex interface).

Features extracted:
- Basic properties (name, position, sequence)
- Structural geometry (backbone and sidechain angles)
- Surface accessibility (SASA)
- Local environment (sequence neighbors, buried/exposed)
- Chemical properties (charge, polarity, hydrophobicity)
- B-factors and occupancy
- Secondary structure prediction from phi/psi

Use cases:
- Protein structure quality assessment
- Sequence-structure analysis
- Stability predictions
- Fold classification features
- Single protein ML applications
"""

from __future__ import annotations
from typing import Dict, List, Optional, Union
from dataclasses import dataclass, field, asdict
import csv
import json
import math

try:
    import pandas as pd
    _HAS_PANDAS = True
except ImportError:
    _HAS_PANDAS = False

try:
    import numpy as np
    _HAS_NUMPY = True
except ImportError:
    _HAS_NUMPY = False

from .models import Complex, Residue, Atom
from .pdbparser import parse_pdb
from .selection import select, selection_residues
from .geometry import analyze_geometry, GeometryAnalysis

# Optional SASA import
try:
    from .sasa import _sasa_fs
    _HAS_FREESASA_WRAPPER = True
except (ImportError, AttributeError):
    _HAS_FREESASA_WRAPPER = False


# ==================== Feature Data Structures ====================

@dataclass
class SingleChainResidueFeatures:
    """
    Feature representation for a single residue in a protein chain.
    Focuses on intrinsic properties rather than interface features.
    """
    # ===== Basic Identification =====
    chain_id: str
    resname: str
    resseq: int
    icode: str = ""
    residue_id: str = ""  # Format: "TYR A:32"

    # ===== Position in Sequence =====
    position_in_chain: int = 0  # 0-based position in chain
    sequence_length: int = 0  # Total length of chain
    relative_position: float = 0.0  # Position / length (0.0 to 1.0)
    is_n_terminal: bool = False  # First 5 residues
    is_c_terminal: bool = False  # Last 5 residues

    # ===== Geometry Features =====
    phi: Optional[float] = None  # Backbone dihedral
    psi: Optional[float] = None  # Backbone dihedral
    omega: Optional[float] = None  # Backbone dihedral
    bend_angle: Optional[float] = None  # Bend angle at this center
    chi1: Optional[float] = None  # Side-chain rotamer
    chi2: Optional[float] = None
    chi3: Optional[float] = None
    chi4: Optional[float] = None

    # ===== Secondary Structure (from phi/psi) =====
    secondary_structure: Optional[str] = None  # "alpha_helix", "beta_sheet", "loop", etc.

    # ===== Surface Accessibility =====
    sasa: Optional[float] = None  # Solvent accessible surface area (Ų)
    relative_sasa: Optional[float] = None  # SASA / max_theoretical_SASA
    is_buried: bool = False  # relative_sasa < 0.2
    is_exposed: bool = False  # relative_sasa > 0.5

    # ===== Local Environment =====
    # CA-CA distances to sequence neighbors
    ca_dist_to_prev: Optional[float] = None  # i-1
    ca_dist_to_next: Optional[float] = None  # i+1
    ca_dist_to_prev2: Optional[float] = None  # i-2
    ca_dist_to_next2: Optional[float] = None  # i+2

    # Average distance to local neighbors (i-2, i-1, i+1, i+2)
    avg_local_ca_distance: Optional[float] = None

    # Centroid coordinates
    centroid_x: Optional[float] = None
    centroid_y: Optional[float] = None
    centroid_z: Optional[float] = None

    # Distance from chain centroid
    dist_from_chain_center: Optional[float] = None

    # ===== Chemical Properties =====
    is_charged: bool = False  # LYS, ARG, ASP, GLU, HIS
    is_polar: bool = False  # SER, THR, ASN, GLN, TYR, CYS
    is_hydrophobic: bool = False  # ALA, VAL, ILE, LEU, MET, PHE, TRP, PRO
    is_aromatic: bool = False  # PHE, TYR, TRP, HIS
    is_positive: bool = False  # LYS, ARG
    is_negative: bool = False  # ASP, GLU
    is_small: bool = False  # GLY, ALA, SER
    is_tiny: bool = False  # GLY, ALA

    # ===== Atom Counts =====
    num_heavy_atoms: int = 0
    num_all_atoms: int = 0
    num_carbon_atoms: int = 0
    num_nitrogen_atoms: int = 0
    num_oxygen_atoms: int = 0
    num_sulfur_atoms: int = 0

    # ===== B-factors (dynamics/uncertainty) =====
    avg_bfactor: Optional[float] = None
    max_bfactor: Optional[float] = None
    min_bfactor: Optional[float] = None
    std_bfactor: Optional[float] = None

    # ===== Occupancy =====
    avg_occupancy: Optional[float] = None
    min_occupancy: Optional[float] = None

    # ===== Custom Features =====
    custom_features: Dict[str, Union[float, int, str, bool]] = field(default_factory=dict)

    def to_dict(self) -> Dict:
        """Convert to dictionary for export."""
        return asdict(self)

    @property
    def is_regular_helix(self) -> bool:
        """Check if in alpha helix region of Ramachandran plot."""
        if self.phi is None or self.psi is None:
            return False
        return -90 <= self.phi <= -30 and -70 <= self.psi <= 10

    @property
    def is_beta_strand(self) -> bool:
        """Check if in beta sheet region of Ramachandran plot."""
        if self.phi is None or self.psi is None:
            return False
        return -180 <= self.phi <= -90 and 90 <= self.psi <= 180


# ==================== Residue Classification ====================

CHARGED_RES = {"LYS", "ARG", "ASP", "GLU", "HIS"}
POSITIVE_RES = {"LYS", "ARG"}
NEGATIVE_RES = {"ASP", "GLU"}
POLAR_RES = {"SER", "THR", "ASN", "GLN", "TYR", "CYS"}
HYDROPHOBIC_RES = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO"}
AROMATIC_RES = {"PHE", "TYR", "TRP", "HIS"}
SMALL_RES = {"GLY", "ALA", "SER"}
TINY_RES = {"GLY", "ALA"}

# Maximum SASA values by residue type (Ų) - from Miller et al. 1987
MAX_SASA = {
    "ALA": 129.0, "CYS": 167.0, "ASP": 193.0, "GLU": 223.0, "PHE": 240.0,
    "GLY": 104.0, "HIS": 224.0, "ILE": 197.0, "LYS": 236.0, "LEU": 201.0,
    "MET": 224.0, "ASN": 195.0, "PRO": 159.0, "GLN": 225.0, "ARG": 274.0,
    "SER": 155.0, "THR": 172.0, "VAL": 174.0, "TRP": 285.0, "TYR": 263.0,
}


def classify_residue_chemistry(resname: str) -> Dict[str, bool]:
    """Return chemical classification flags for a residue."""
    resname = resname.upper()
    return {
        'is_charged': resname in CHARGED_RES,
        'is_positive': resname in POSITIVE_RES,
        'is_negative': resname in NEGATIVE_RES,
        'is_polar': resname in POLAR_RES,
        'is_hydrophobic': resname in HYDROPHOBIC_RES,
        'is_aromatic': resname in AROMATIC_RES,
        'is_small': resname in SMALL_RES,
        'is_tiny': resname in TINY_RES,
    }


# ==================== Main Feature Extraction ====================

def extract_single_chain_features(
    cx: Complex,
    chain_ids: Union[str, List[str]],
    compute_sasa: bool = True,
    compute_geometry: bool = True,
    verbose: bool = True,
) -> Dict[str, SingleChainResidueFeatures]:
    """
    Extract comprehensive residue-level features for single chain(s).

    Args:
        cx: Complex structure
        chain_ids: Chain ID(s) to analyze (string or list of strings)
        compute_sasa: Calculate solvent accessible surface area
        compute_geometry: Calculate backbone/sidechain angles
        verbose: Print progress messages

    Returns:
        Dictionary mapping residue_id -> SingleChainResidueFeatures

    Example:
        >>> cx = parse_pdb("protein.pdb")
        >>> features = extract_single_chain_features(cx, "A")
        >>> df = features_to_dataframe(features)
    """

    if verbose:
        print("Starting single chain residue feature extraction...")

    # Normalize chain IDs
    if isinstance(chain_ids, str):
        chain_ids = [chain_ids] if ',' not in chain_ids else [c.strip() for c in chain_ids.split(',')]

    # Initialize feature storage
    features: Dict[str, SingleChainResidueFeatures] = {}

    # Get all residues for selected chains
    all_residues: List[Residue] = []
    chain_residue_lists: Dict[str, List[Residue]] = {}

    for chain_id in chain_ids:
        if chain_id not in cx.chains:
            if verbose:
                print(f"  Warning: Chain '{chain_id}' not found in complex")
            continue

        chain = cx.chains[chain_id]
        residues = list(chain.iter_residues())
        all_residues.extend(residues)
        chain_residue_lists[chain_id] = residues

    if verbose:
        print(f"  Total residues to analyze: {len(all_residues)}")
        print(f"  Chains: {', '.join(chain_ids)}")

    # ===== STEP 1: Initialize basic features =====
    if verbose:
        print("  [1/5] Initializing basic features...")

    for chain_id, residues in chain_residue_lists.items():
        chain_length = len(residues)

        for idx, residue in enumerate(residues):
            res_id = residue.id_str

            # Basic identification
            feat = SingleChainResidueFeatures(
                chain_id=residue.chain_id,
                resname=residue.resname,
                resseq=residue.resseq,
                icode=residue.icode,
                residue_id=res_id
            )

            # Position in sequence
            feat.position_in_chain = idx
            feat.sequence_length = chain_length
            feat.relative_position = idx / chain_length if chain_length > 0 else 0.0
            feat.is_n_terminal = idx < 5
            feat.is_c_terminal = idx >= chain_length - 5

            # Chemical classification
            chem_props = classify_residue_chemistry(residue.resname)
            feat.is_charged = chem_props['is_charged']
            feat.is_positive = chem_props['is_positive']
            feat.is_negative = chem_props['is_negative']
            feat.is_polar = chem_props['is_polar']
            feat.is_hydrophobic = chem_props['is_hydrophobic']
            feat.is_aromatic = chem_props['is_aromatic']
            feat.is_small = chem_props['is_small']
            feat.is_tiny = chem_props['is_tiny']

            # Atom counts
            atoms = list(residue.iter_atoms(ignore_h=True))
            feat.num_heavy_atoms = len(atoms)
            feat.num_all_atoms = len(residue.atoms)

            # Count by element
            for atom in atoms:
                if atom.element == 'C':
                    feat.num_carbon_atoms += 1
                elif atom.element == 'N':
                    feat.num_nitrogen_atoms += 1
                elif atom.element == 'O':
                    feat.num_oxygen_atoms += 1
                elif atom.element == 'S':
                    feat.num_sulfur_atoms += 1

            # B-factors and occupancy
            if atoms:
                bfactors = [a.bfactor for a in atoms]
                occupancies = [a.occupancy for a in atoms]

                feat.avg_bfactor = sum(bfactors) / len(bfactors)
                feat.max_bfactor = max(bfactors)
                feat.min_bfactor = min(bfactors)

                if len(bfactors) > 1:
                    mean_bf = feat.avg_bfactor
                    variance = sum((b - mean_bf) ** 2 for b in bfactors) / len(bfactors)
                    feat.std_bfactor = math.sqrt(variance)

                feat.avg_occupancy = sum(occupancies) / len(occupancies)
                feat.min_occupancy = min(occupancies)

            # Centroid calculation
            if atoms:
                x_coords = [a.x for a in atoms]
                y_coords = [a.y for a in atoms]
                z_coords = [a.z for a in atoms]

                feat.centroid_x = sum(x_coords) / len(x_coords)
                feat.centroid_y = sum(y_coords) / len(y_coords)
                feat.centroid_z = sum(z_coords) / len(z_coords)

            features[res_id] = feat

    # ===== STEP 2: Calculate chain centroids and distances =====
    if verbose:
        print("  [2/5] Computing distances from chain centers...")

    for chain_id, residues in chain_residue_lists.items():
        # Calculate chain centroid
        all_centroids = []
        for res in residues:
            res_id = res.id_str
            if res_id in features:
                feat = features[res_id]
                if feat.centroid_x is not None:
                    all_centroids.append((feat.centroid_x, feat.centroid_y, feat.centroid_z))

        if all_centroids:
            chain_center_x = sum(c[0] for c in all_centroids) / len(all_centroids)
            chain_center_y = sum(c[1] for c in all_centroids) / len(all_centroids)
            chain_center_z = sum(c[2] for c in all_centroids) / len(all_centroids)

            # Calculate distance from chain center for each residue
            for res in residues:
                res_id = res.id_str
                if res_id in features:
                    feat = features[res_id]
                    if feat.centroid_x is not None:
                        dx = feat.centroid_x - chain_center_x
                        dy = feat.centroid_y - chain_center_y
                        dz = feat.centroid_z - chain_center_z
                        feat.dist_from_chain_center = math.sqrt(dx*dx + dy*dy + dz*dz)

    # ===== STEP 3: Local CA-CA distances =====
    if verbose:
        print("  [3/5] Computing local CA-CA distances...")

    for chain_id, residues in chain_residue_lists.items():
        for idx, residue in enumerate(residues):
            res_id = residue.id_str
            if res_id not in features:
                continue

            feat = features[res_id]

            # Get CA atom
            ca = next((a for a in residue.atoms if a.name.strip().upper() == "CA"), None)
            if ca is None:
                continue

            local_dists = []

            # Distance to i-2
            if idx >= 2:
                res_prev2 = residues[idx - 2]
                ca_prev2 = next((a for a in res_prev2.atoms if a.name.strip().upper() == "CA"), None)
                if ca_prev2:
                    dist = math.sqrt((ca.x - ca_prev2.x)**2 + (ca.y - ca_prev2.y)**2 + (ca.z - ca_prev2.z)**2)
                    feat.ca_dist_to_prev2 = dist
                    local_dists.append(dist)

            # Distance to i-1
            if idx >= 1:
                res_prev = residues[idx - 1]
                ca_prev = next((a for a in res_prev.atoms if a.name.strip().upper() == "CA"), None)
                if ca_prev:
                    dist = math.sqrt((ca.x - ca_prev.x)**2 + (ca.y - ca_prev.y)**2 + (ca.z - ca_prev.z)**2)
                    feat.ca_dist_to_prev = dist
                    local_dists.append(dist)

            # Distance to i+1
            if idx < len(residues) - 1:
                res_next = residues[idx + 1]
                ca_next = next((a for a in res_next.atoms if a.name.strip().upper() == "CA"), None)
                if ca_next:
                    dist = math.sqrt((ca.x - ca_next.x)**2 + (ca.y - ca_next.y)**2 + (ca.z - ca_next.z)**2)
                    feat.ca_dist_to_next = dist
                    local_dists.append(dist)

            # Distance to i+2
            if idx < len(residues) - 2:
                res_next2 = residues[idx + 2]
                ca_next2 = next((a for a in res_next2.atoms if a.name.strip().upper() == "CA"), None)
                if ca_next2:
                    dist = math.sqrt((ca.x - ca_next2.x)**2 + (ca.y - ca_next2.y)**2 + (ca.z - ca_next2.z)**2)
                    feat.ca_dist_to_next2 = dist
                    local_dists.append(dist)

            # Average local distance
            if local_dists:
                feat.avg_local_ca_distance = sum(local_dists) / len(local_dists)

    # ===== STEP 4: Geometry analysis =====
    if compute_geometry and verbose:
        print("  [4/5] Computing backbone and sidechain angles...")

    if compute_geometry:
        try:
            geom_analysis = analyze_geometry(
                cx,
                chain_ids,
                compute_backbone=True,
                compute_sidechains=True,
                compute_bends=True
            )

            # Extract backbone angles (phi/psi/omega)
            if geom_analysis.backbone:
                for backbone_angles in geom_analysis.backbone:
                    res_id = backbone_angles.residue_id
                    if res_id in features:
                        features[res_id].phi = backbone_angles.phi
                        features[res_id].psi = backbone_angles.psi
                        features[res_id].omega = backbone_angles.omega

                        # Assign secondary structure
                        if backbone_angles.phi and backbone_angles.psi:
                            phi, psi = backbone_angles.phi, backbone_angles.psi

                            # Alpha helix
                            if -90 <= phi <= -30 and -70 <= psi <= 10:
                                features[res_id].secondary_structure = "alpha_helix"
                            # Beta sheet
                            elif -180 <= phi <= -90 and 90 <= psi <= 180:
                                features[res_id].secondary_structure = "beta_sheet"
                            # Left-handed helix
                            elif 30 <= phi <= 100 and 0 <= psi <= 80:
                                features[res_id].secondary_structure = "left_helix"
                            else:
                                features[res_id].secondary_structure = "loop"

            # Extract sidechain chi angles
            if geom_analysis.sidechains:
                for sidechain_chi in geom_analysis.sidechains:
                    res_id = sidechain_chi.residue_id
                    if res_id in features:
                        chis = sidechain_chi.chis
                        features[res_id].chi1 = chis.get('chi1')
                        features[res_id].chi2 = chis.get('chi2')
                        features[res_id].chi3 = chis.get('chi3')
                        features[res_id].chi4 = chis.get('chi4')

            # Extract bend angles
            if geom_analysis.bends:
                for bent in geom_analysis.bends:
                    res_id = bent.residue_id
                    if res_id in features:
                        features[res_id].bend_angle = bent.angle_deg

        except Exception as e:
            if verbose:
                print(f"    Warning: Geometry calculation failed: {e}")

    # ===== STEP 5: SASA calculation =====
    if compute_sasa and verbose:
        print("  [5/5] Computing solvent accessible surface area...")

    if compute_sasa and _HAS_FREESASA_WRAPPER:
        try:
            # Compute SASA for the chains
            sasa_results = {}

            for chain_id in chain_ids:
                if chain_id not in cx.chains:
                    continue

                chain = cx.chains[chain_id]
                residues = list(chain.iter_residues())

                # Prepare atoms for FreeSASA
                atoms_list = []
                for res in residues:
                    for atom in res.iter_atoms(ignore_h=True):
                        atoms_list.append(atom)

                if atoms_list:
                    # Call FreeSASA
                    result = _sasa_fs(atoms_list, probe_radius=1.4)

                    # Extract per-residue SASA
                    for res in residues:
                        res_sasa = 0.0
                        for atom in res.iter_atoms(ignore_h=True):
                            atom_key = (atom.chain_id, atom.resseq, atom.icode, atom.name.strip())
                            res_sasa += result.get(atom_key, 0.0)

                        res_id = res.id_str
                        if res_id in features:
                            features[res_id].sasa = res_sasa

                            # Calculate relative SASA
                            max_sasa = MAX_SASA.get(res.resname.upper(), 200.0)
                            features[res_id].relative_sasa = min(res_sasa / max_sasa, 1.0)

                            # Classify as buried/exposed
                            features[res_id].is_buried = features[res_id].relative_sasa < 0.2
                            features[res_id].is_exposed = features[res_id].relative_sasa > 0.5

        except Exception as e:
            if verbose:
                print(f"    Warning: SASA calculation failed: {e}")
    elif compute_sasa and not _HAS_FREESASA_WRAPPER:
        if verbose:
            print("    Warning: FreeSASA not available, skipping SASA calculation")

    if verbose:
        print(f"  Feature extraction complete! Generated features for {len(features)} residues.")

        # Summary statistics
        if features:
            helix_count = sum(1 for f in features.values() if f.secondary_structure == "alpha_helix")
            sheet_count = sum(1 for f in features.values() if f.secondary_structure == "beta_sheet")
            loop_count = sum(1 for f in features.values() if f.secondary_structure == "loop")

            print(f"  - Alpha helix: {helix_count}")
            print(f"  - Beta sheet:  {sheet_count}")
            print(f"  - Loop/other:  {loop_count}")

            if compute_sasa and _HAS_FREESASA_WRAPPER:
                buried_count = sum(1 for f in features.values() if f.is_buried)
                exposed_count = sum(1 for f in features.values() if f.is_exposed)
                print(f"  - Buried residues:  {buried_count}")
                print(f"  - Exposed residues: {exposed_count}")

    return features


# ==================== Export Functions ====================

def features_to_dataframe(features: Dict[str, SingleChainResidueFeatures]) -> 'pd.DataFrame':
    """
    Convert feature dictionary to pandas DataFrame.

    Args:
        features: Dictionary from extract_single_chain_features()

    Returns:
        DataFrame with one row per residue
    """
    if not _HAS_PANDAS:
        raise ImportError("pandas required for features_to_dataframe()")

    rows = []
    for res_id, feat in features.items():
        row = feat.to_dict()
        # Exclude custom_features dict for simplicity (can be flattened if needed)
        if 'custom_features' in row and not row['custom_features']:
            row.pop('custom_features')
        rows.append(row)

    df = pd.DataFrame(rows)

    # Sort by chain and position
    df = df.sort_values(['chain_id', 'resseq', 'icode'], ignore_index=True)

    return df


def export_features_csv(
    features: Dict[str, SingleChainResidueFeatures],
    output_path: str,
    include_empty_columns: bool = False
) -> None:
    """Export features to CSV file."""
    if not features:
        raise ValueError("No features to export")

    sample_feat = next(iter(features.values()))
    all_columns = list(sample_feat.to_dict().keys())

    # Remove custom_features
    if 'custom_features' in all_columns:
        all_columns.remove('custom_features')

    # Filter columns if requested
    if not include_empty_columns:
        active_columns = []
        for col in all_columns:
            has_value = False
            for feat in features.values():
                val = getattr(feat, col, None)
                if val is not None and val != '' and val != {} and val != []:
                    has_value = True
                    break
            if has_value:
                active_columns.append(col)
        columns = active_columns
    else:
        columns = all_columns

    # Write CSV
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()

        for res_id in sorted(features.keys()):
            feat = features[res_id]
            row = feat.to_dict()
            row = {k: v for k, v in row.items() if k in columns}
            writer.writerow(row)

    print(f"✓ Exported {len(features)} residues to {output_path}")


def export_features_json(
    features: Dict[str, SingleChainResidueFeatures],
    output_path: str,
    pretty: bool = True
) -> None:
    """Export features to JSON file."""
    data = {res_id: feat.to_dict() for res_id, feat in features.items()}

    with open(output_path, 'w') as f:
        if pretty:
            json.dump(data, f, indent=2)
        else:
            json.dump(data, f)

    print(f"✓ Exported {len(features)} residues to {output_path}")


# ==================== Analysis Utilities ====================

def get_feature_summary(features: Dict[str, SingleChainResidueFeatures]) -> Dict:
    """Generate summary statistics."""
    summary = {
        'total_residues': len(features),
        'chains': len(set(f.chain_id for f in features.values())),
    }

    # Secondary structure
    summary['alpha_helix'] = sum(1 for f in features.values() if f.secondary_structure == "alpha_helix")
    summary['beta_sheet'] = sum(1 for f in features.values() if f.secondary_structure == "beta_sheet")
    summary['loop'] = sum(1 for f in features.values() if f.secondary_structure == "loop")

    # Chemical properties
    summary['charged_residues'] = sum(1 for f in features.values() if f.is_charged)
    summary['hydrophobic_residues'] = sum(1 for f in features.values() if f.is_hydrophobic)
    summary['aromatic_residues'] = sum(1 for f in features.values() if f.is_aromatic)
    summary['polar_residues'] = sum(1 for f in features.values() if f.is_polar)

    # SASA
    summary['buried_residues'] = sum(1 for f in features.values() if f.is_buried)
    summary['exposed_residues'] = sum(1 for f in features.values() if f.is_exposed)

    # B-factors
    bfactors = [f.avg_bfactor for f in features.values() if f.avg_bfactor is not None]
    if bfactors:
        summary['avg_bfactor'] = sum(bfactors) / len(bfactors)
        summary['max_bfactor'] = max(bfactors)
        summary['min_bfactor'] = min(bfactors)

    return summary
