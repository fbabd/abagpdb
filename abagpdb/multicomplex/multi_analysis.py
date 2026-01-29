"""
Modular Multi-Complex Analysis System
Allows users to select which properties to calculate and compare
"""

from __future__ import annotations
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field
import traceback

# Import analysis modules
from ..pdbparser import parse_pdb
from ..sasa_analysis import compare_bound_unbound_sasa
from ..contacts import analyze_contacts
from ..interface import compute_interface
from ..vdw import per_residue_LJ_decomposition
# geometry module not used yet - will add RMSD later if needed


@dataclass
class AnalysisConfig:
    """Configuration for which analyses to run"""
    calculate_sasa: bool = True
    calculate_contacts: bool = True
    calculate_interface: bool = True
    calculate_vdw: bool = False
    # RMSD removed - not implemented yet
    
    # Advanced options
    sasa_probe_radius: float = 1.4
    contact_distance: float = 5.0
    vdw_cutoff: float = 5.0
    interface_distance: float = 5.0


@dataclass
class ComplexData:
    """Data container for a single complex"""
    name: str
    filepath: str
    structure: object = None
    
    # Analysis results (populated as needed)
    sasa_data: Dict = field(default_factory=dict)
    contact_data: Dict = field(default_factory=dict)
    interface_data: Dict = field(default_factory=dict)
    vdw_data: Dict = field(default_factory=dict)
    # rmsd_data removed - not implemented yet
    
    # Summary statistics
    summary: Dict = field(default_factory=dict)


class MultiComplexAnalyzer:
    """
    Main class for multi-complex comparison with modular analysis selection
    """
    
    def __init__(self, pdb_files: List[Tuple[str, str]]):
        """
        Initialize analyzer with PDB files
        
        Args:
            pdb_files: List of (name, filepath) tuples
        """
        self.complexes: Dict[str, ComplexData] = {}
        self.wt_name: Optional[str] = None
        self.chain_groups: Dict[str, List[str]] = {}
        
        # Load all structures
        print(f"\nLoading {len(pdb_files)} structures...")
        for name, filepath in pdb_files:
            try:
                structure = parse_pdb(filepath)
                self.complexes[name] = ComplexData(
                    name=name,
                    filepath=filepath,
                    structure=structure
                )
                print(f"  ✓ Loaded {name}")
            except Exception as e:
                print(f"  ✗ Failed to load {name}: {e}")
                raise
    
    
    def set_reference(self, wt_name: str):
        """Set the wild-type reference structure"""
        if wt_name not in self.complexes:
            raise ValueError(f"Reference '{wt_name}' not found in loaded structures")
        self.wt_name = wt_name
        print(f"Reference set to: {wt_name}")
    
    
    def set_chain_groups(self, chain_groups: Dict[str, List[str]]):
        """Set chain groupings for analysis"""
        self.chain_groups = chain_groups
        print(f"Chain groups configured: {list(chain_groups.keys())}")
    
    
    def run_analysis(self, config: AnalysisConfig) -> Dict[str, ComplexData]:
        """
        Run selected analyses on all complexes
        
        Args:
            config: Configuration specifying which analyses to run
            
        Returns:
            Dictionary of complex name -> ComplexData with results
        """
        print("\n" + "="*60)
        print("STARTING MULTI-COMPLEX ANALYSIS")
        print("="*60)
        
        analyses_enabled = []
        if config.calculate_sasa:
            analyses_enabled.append("SASA")
        if config.calculate_contacts:
            analyses_enabled.append("Contacts")
        if config.calculate_interface:
            analyses_enabled.append("Interface")
        if config.calculate_vdw:
            analyses_enabled.append("VDW")
        
        print(f"Analyses enabled: {', '.join(analyses_enabled)}")
        print(f"Processing {len(self.complexes)} structures...")
        print()
        
        for name, complex_data in self.complexes.items():
            print(f"Analyzing: {name}")
            print("-" * 40)
            
            # Run each selected analysis
            if config.calculate_sasa:
                self._calculate_sasa(complex_data, config)
            
            if config.calculate_contacts:
                self._calculate_contacts(complex_data, config)
            
            if config.calculate_interface:
                self._calculate_interface(complex_data, config)
            
            if config.calculate_vdw:
                self._calculate_vdw(complex_data, config)
            
            print()
        
        print("="*60)
        print("ANALYSIS COMPLETE")
        print("="*60)
        
        return self.complexes
    
    
    def _calculate_sasa(self, complex_data: ComplexData, config: AnalysisConfig):
        """Calculate SASA for a complex"""
        try:
            print("  → Calculating SASA...")
            chain_ids = list(complex_data.structure.chains.keys())
            
            sasa_results = compare_bound_unbound_sasa(
                complex_data.structure, 
                chain_ids,
                probe=config.sasa_probe_radius 
            )
            
            complex_data.sasa_data = sasa_results
            
            # Calculate summary stats
            total_buried = sum(
                res_data.get('buried_sasa', 0) 
                for res_data in sasa_results.values()
            )
            complex_data.summary['total_buried_sasa'] = total_buried
            complex_data.summary['num_residues_with_sasa'] = len(sasa_results)
            
            print(f"     ✓ SASA calculated ({len(sasa_results)} residues)")
            
        except Exception as e:
            print(f"     ✗ SASA calculation failed: {e}")
            traceback.print_exc()
            complex_data.sasa_data = {}
    
    
    def _calculate_contacts(self, complex_data: ComplexData, config: AnalysisConfig):
        """Calculate inter-chain contacts"""
        try:
            print("  → Calculating contacts...")
            chain_ids = list(complex_data.structure.chains.keys())
            
            if len(chain_ids) < 2:
                print("     ⚠ Skipped (only 1 chain)")
                return
            
            # Analyze contacts between first chain and others
            contact_analysis = analyze_contacts( 
                complex_data.structure,
                selection_A=chain_ids[0],
                selection_B=chain_ids[1:]
            )
            
            complex_data.contact_data = {
                'residue_summary': contact_analysis.get_residue_summary(),
                'bond_counts': contact_analysis.get_contact_counts()
            }
            
            # Store summary
            bond_counts = contact_analysis.get_contact_counts()
            for bond_type, count in bond_counts.items():
                complex_data.summary[f'contacts_{bond_type.name}'] = count
            
            total_contacts = sum(bond_counts.values())
            complex_data.summary['total_contacts'] = total_contacts
            
            print(f"     ✓ Contacts calculated ({total_contacts} total)")
            
        except Exception as e:
            print(f"     ✗ Contact calculation failed: {e}")
            traceback.print_exc()
            complex_data.contact_data = {}
    
    
    def _calculate_interface(self, complex_data: ComplexData, config: AnalysisConfig):
        """Calculate interface properties"""
        try:
            print("  → Calculating interface...")
            chain_ids = list(complex_data.structure.chains.keys())
            
            if len(chain_ids) < 2:
                print("     ⚠ Skipped (only 1 chain)")
                return
            
            interface_analysis = compute_interface(
                complex_data.structure,
                selection_A=chain_ids[0],
                selection_B=chain_ids[1:],
                cutoff=config.interface_distance
            )
            
            # Convert InterfaceAnalysis object to dictionary format expected by dashboard
            interface_residues = interface_analysis.get_interface_residues()
            residues_dict = {}
            
            for res_data in interface_residues:
                # res_data is ResidueInterfaceData object
                res_id = res_data.id_str
                residues_dict[res_id] = {
                    'is_interface': res_data.is_interface,
                    'contact_count': res_data.contact_count,
                    'partner_residues': res_data.partner_residues
                }
            
            # Store as dictionary
            complex_data.interface_data = {
                'residues': residues_dict,
                'total_contacts': interface_analysis.total_contacts,
                'bsa_total': interface_analysis.bsa_total,
                'cutoff': interface_analysis.cutoff,
                'group_A_chains': interface_analysis.group_A_chains,
                'group_B_chains': interface_analysis.group_B_chains
            }
            
            # Store summary
            complex_data.summary['interface_residues'] = interface_analysis.total_contacts
            complex_data.summary['buried_surface_area'] = interface_analysis.bsa_total or 0.0
            
            
            print(f"     ✓ Interface calculated ({interface_analysis.total_contacts} contacts)")
            
        except Exception as e:
            print(f"     ✗ Interface calculation failed: {e}")
            traceback.print_exc()
            complex_data.interface_data = {}
    
    
    def _calculate_vdw(self, complex_data: ComplexData, config: AnalysisConfig):
        """Calculate VDW energies"""
        try:
            print("  → Calculating VDW energies...")
            chain_ids = list(complex_data.structure.chains.keys())
            
            if len(chain_ids) < 2:
                print("     ⚠ Skipped (only 1 chain)")
                return
            
            # Calculate VDW for each chain against others
            all_vdw = {}
            for i, chain in enumerate(chain_ids):
                group_a = [chain]
                group_b = [ch for j, ch in enumerate(chain_ids) if j != i]
                
                vdw_energies = per_residue_LJ_decomposition(
                    complex_data.structure,
                    group_a,
                    group_b,
                    cutoff=config.vdw_cutoff
                )
                all_vdw.update(vdw_energies)
            
            complex_data.vdw_data = all_vdw
            
            # Calculate summary
            total_vdw = sum(all_vdw.values())
            favorable_vdw = sum(e for e in all_vdw.values() if e < 0)
            
            complex_data.summary['total_vdw_energy'] = total_vdw
            complex_data.summary['favorable_vdw_energy'] = favorable_vdw
            complex_data.summary['num_residues_with_vdw'] = len(all_vdw)
            
            print(f"     ✓ VDW calculated ({len(all_vdw)} residues, total: {total_vdw:.2f} kcal/mol)")
            
        except Exception as e:
            print(f"     ✗ VDW calculation failed: {e}")
            traceback.print_exc()
            complex_data.vdw_data = {}
    
    
    def get_summary_table(self) -> Dict[str, Dict]:
        """
        Get summary statistics for all complexes
        
        Returns:
            Dictionary of complex name -> summary statistics
        """
        return {
            name: complex_data.summary
            for name, complex_data in self.complexes.items()
        }
    
    
    def export_results(self) -> Dict:
        """
        Export all results in format compatible with existing multi_dashboard.py
        
        Returns:
            Dictionary containing all analysis results in legacy format
        """
        # Convert to format expected by existing dashboard
        results = {}
        
        for name, complex_data in self.complexes.items():
            # Build residue-level data structure
            residues = {}
            
            # Collect all unique residues
            all_residue_ids = set()
            if complex_data.sasa_data:
                all_residue_ids.update(complex_data.sasa_data.keys())
            if complex_data.contact_data and 'residue_summary' in complex_data.contact_data:
                all_residue_ids.update(complex_data.contact_data['residue_summary'].keys())
            if complex_data.interface_data and 'residues' in complex_data.interface_data:
                all_residue_ids.update(complex_data.interface_data['residues'].keys())
            if complex_data.vdw_data:
                all_residue_ids.update(complex_data.vdw_data.keys())
            
            # Build per-residue data
            for res_id in all_residue_ids:
                residues[res_id] = {
                    'sasa': complex_data.sasa_data.get(res_id, {}),
                    'bonds': complex_data.contact_data.get('residue_summary', {}).get(res_id, {}),
                    'interface': complex_data.interface_data.get('residues', {}).get(res_id, 0),
                    'vdw_energy': complex_data.vdw_data.get(res_id, 0.0)
                }
            
            # Build summary data
            summary = complex_data.summary.copy()
            
            # Add contact counts if available
            if complex_data.contact_data and 'bond_counts' in complex_data.contact_data:
                for bond_type, count in complex_data.contact_data['bond_counts'].items():
                    summary[f'contacts_{bond_type}'] = count
            
            # Build result entry
            results[name] = {
                'residues': residues,
                'summary': summary
            }
        
        return results
    
    