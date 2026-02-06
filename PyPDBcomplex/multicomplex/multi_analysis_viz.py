"""
Multi-Complex Comparison Visualizations

This module provides interactive HTML visualizations for comparing multiple protein
complexes (WT and variants), showing mutations, property changes, and structural differences.

"""

from __future__ import annotations
from typing import Dict, List, Optional, Any, Tuple
import json
from collections import defaultdict


def sort_key(position_key: str) -> Tuple:
    """
    Sort key for position strings like 'A:10', 'A:100', 'B:5'
    Returns tuple of (chain, numeric_position) for proper sorting
    """
    try:
        chain, pos = position_key.split(':')
        # Try to extract numeric part from position (handle insertion codes)
        pos_num = int(''.join(filter(str.isdigit, pos)))
        return (chain, pos_num, pos)
    except:
        return (position_key, 0, '')


def _parse_interface_data(interface_data):
    """Helper to parse interface data consistently"""
    if isinstance(interface_data, dict):
        return {
            'is_interface': True,
            'partners': interface_data.get('partners', []),
            'contact_count': interface_data.get('contact_count', 0)
        }
    else:
        return {
            'is_interface': False,
            'partners': [],
            'contact_count': 0
        }


def _parse_bonds(bonds_data):
    """Helper to parse bonds data, handling empty dicts"""
    if not bonds_data:  # Handle empty dict or None
        return {
            'counts': {
                'hydrogen_bond': 0,
                'salt_bridge': 0,
                'disulfide': 0,
                'hydrophobic': 0,
                'pi_stacking': 0,
                'total': 0
            },
            'details': {}
        }
    
    bond_summary = {}
    bond_types = ["hydrogen_bond",  "salt_bridge", "disulfide", "hydrophobic", "pi_stacking"]
    for bond_type, bond_count in bonds_data.items():
        if bond_type in bond_types:
            bond_summary[bond_type] = bond_count 
    
    bond_summary['total'] = sum(bond_summary.values())
    
    return {
        'counts': bond_summary,
        'details': bonds_data
    }
    


def compare_residues(results):
    """
    Compare individual residue properties across variants.
    Handles mutations, insertions, and deletions.
    """
    variants = list(results.keys())
    
    # Get all unique residue positions (chain:position) across all variants
    all_positions = set()
    residue_identity_map = {}  # Maps position to residue identities in each variant
    
    for variant in variants:
        for res_id in results[variant]['residues'].keys():
            # Parse residue ID: e.g., "SER A:32" -> residue_type="SER", chain="A", position="32"
            parts = res_id.split()
            if len(parts) == 2:
                res_type = parts[0]
                chain_pos = parts[1]  # "A:32"
                chain, position = chain_pos.split(':')
                
                position_key = f"{chain}:{position}"
                all_positions.add(position_key)
                
                if position_key not in residue_identity_map:
                    residue_identity_map[position_key] = {}
                residue_identity_map[position_key][variant] = {
                    'full_id': res_id,
                    'residue_type': res_type,
                    'chain': chain,
                    'position': position
                }

    residue_comparisons = {}
    
    for position_key in sorted(all_positions, key=sort_key):
        residue_comparisons[position_key] = {
            'position': position_key,
            'variants': {},
            'mutation_status': {},
            'comparison_type': None  # Will be: 'identical', 'mutation', 'insertion', 'deletion'
        } 
        
        # Determine comparison type and mutation status
        present_variants = list(residue_identity_map[position_key].keys())
        residue_types = [residue_identity_map[position_key][v]['residue_type'] 
                        for v in present_variants]
        
        # Check if position exists in all variants
        if len(present_variants) == len(variants):
            # Position exists in all variants
            if len(set(residue_types)) == 1:
                residue_comparisons[position_key]['comparison_type'] = 'identical'
            else:
                residue_comparisons[position_key]['comparison_type'] = 'mutation'
                # Store what mutated to what
                if 'WT' in present_variants:
                    wt_type = residue_identity_map[position_key]['WT']['residue_type']
                    for variant in present_variants:
                        if variant != 'WT':
                            var_type = residue_identity_map[position_key][variant]['residue_type']
                            if wt_type != var_type:
                                residue_comparisons[position_key]['mutation_status'][variant] = {
                                    'from': wt_type,
                                    'to': var_type,
                                    'mutation': f"{wt_type}->{var_type}"
                                }
        else:
            # Position doesn't exist in all variants (insertion/deletion)
            if 'WT' in present_variants:
                residue_comparisons[position_key]['comparison_type'] = 'deletion'
                wt_type = residue_identity_map[position_key]['WT']['residue_type']
                for variant in variants:
                    if variant not in present_variants:
                        residue_comparisons[position_key]['mutation_status'][variant] = {
                            'status': 'deleted',
                            'wt_residue': wt_type
                        }
            else:
                residue_comparisons[position_key]['comparison_type'] = 'insertion'
                for variant in present_variants:
                    var_type = residue_identity_map[position_key][variant]['residue_type']
                    residue_comparisons[position_key]['mutation_status'][variant] = {
                        'status': 'inserted',
                        'residue_type': var_type
                    }
        
        # Store data for each variant
        for variant in variants:
            if variant in residue_identity_map[position_key]:
                res_id = residue_identity_map[position_key][variant]['full_id']
                res_type = residue_identity_map[position_key][variant]['residue_type']
                res_data = results[variant]['residues'][res_id]
    
                residue_comparisons[position_key]['variants'][variant] = {
                    'present': True,
                    'residue_id': res_id,
                    'residue_type': res_type,
                    'sasa': {
                        'buried_fraction': res_data['sasa']['buried_fraction'],
                        'delta': res_data['sasa']['delta'],
                        'bound': res_data['sasa']['bound'],
                        'unbound': res_data['sasa']['unbound']
                    },
                    'interface': _parse_interface_data(res_data['interface']),
                    'vdw_energy': res_data['vdw_energy'],
                    'bonds': _parse_bonds(res_data['bonds'])
                }
            else:
                # Position doesn't exist in this variant
                residue_comparisons[position_key]['variants'][variant] = {
                    'present': False,
                    'residue_id': None,
                    'residue_type': None,
                    'sasa': None,
                    'interface': None,
                    'vdw_energy': None,
                    'bonds': None
                }
    residue_comparisons = dict(sorted(residue_comparisons.items(), key=lambda x: sort_key(x[0]))) 
    return residue_comparisons


def generate_summary_comparison_html(results: Dict[str, Dict], title: str = "Summary Comparison") -> str:
    """
    Generate HTML visualization comparing summary statistics across variants.
    
    Args:
        results: Dictionary with variant names as keys, each containing 'summary' dict
        title: Title for the visualization
        
    Returns:
        HTML string with interactive summary comparison
    """
    variants = list(results.keys())
    
    # Extract summary data
    metrics = ['HYDROGEN_BOND', 'SALT_BRIDGE', 'HYDROPHOBIC', 'n_res_interface', 
               'buried_surface_area', 'total_contacts']
    
    summary_data = {}
    for metric in metrics:
        summary_data[metric] = {v: results[v]['summary'].get(metric, 0) for v in variants}
    
    # Calculate deltas from WT if available
    has_wt = 'WT' in variants
    
    html = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>{title}</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f5f5f5;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2c3e50;
            margin-bottom: 10px;
        }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .metric-card {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #3498db;
        }}
        .metric-card h3 {{
            margin: 0 0 15px 0;
            color: #2c3e50;
            font-size: 16px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        .metric-values {{
            display: flex;
            flex-direction: column;
            gap: 10px;
        }}
        .value-row {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 8px;
            background: white;
            border-radius: 4px;
        }}
        .value-row .variant-name {{
            font-weight: 600;
            color: #34495e;
        }}
        .value-row .value {{
            font-size: 18px;
            font-weight: bold;
            color: #3498db;
        }}
        .value-row .delta {{
            font-size: 14px;
            margin-left: 10px;
            padding: 2px 8px;
            border-radius: 12px;
        }}
        .delta.positive {{
            background: #e8f5e9;
            color: #2e7d32;
        }}
        .delta.negative {{
            background: #ffebee;
            color: #c62828;
        }}
        .delta.neutral {{
            background: #f5f5f5;
            color: #757575;
        }}
        .comparison-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .comparison-table th {{
            background: #34495e;
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: 600;
        }}
        .comparison-table td {{
            padding: 12px;
            border-bottom: 1px solid #ecf0f1;
        }}
        .comparison-table tr:hover {{
            background: #f8f9fa;
        }}
        .variant-badge {{
            display: inline-block;
            padding: 4px 12px;
            border-radius: 12px;
            font-size: 12px;
            font-weight: 600;
            background: #3498db;
            color: white;
        }}
        .wt-badge {{
            background: #2ecc71;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>{title}</h1>
        
        <h2>Summary Metrics Comparison</h2>
        <div class="summary-grid">
"""
    
    # Generate metric cards
    for metric in metrics:
        values = summary_data[metric]
        
        html += f"""
            <div class="metric-card">
                <h3>{metric.replace('_', ' ').title()}</h3>
                <div class="metric-values">
"""
        
        wt_value = values.get('WT', 0) if has_wt else None
        
        for variant in variants:
            value = values[variant]
            
            # Calculate delta if WT exists
            delta_html = ""
            if has_wt and variant != 'WT' and wt_value is not None:
                delta = value - wt_value
                pct_change = ((value - wt_value) / wt_value * 100) if wt_value != 0 else 0
                
                delta_class = "positive" if delta > 0 else ("negative" if delta < 0 else "neutral")
                delta_sign = "+" if delta > 0 else ""
                
                if isinstance(value, float):
                    delta_html = f'<span class="delta {delta_class}">{delta_sign}{delta:.2f} ({delta_sign}{pct_change:.1f}%)</span>'
                else:
                    delta_html = f'<span class="delta {delta_class}">{delta_sign}{delta} ({delta_sign}{pct_change:.1f}%)</span>'
            
            value_str = f"{value:.2f}" if isinstance(value, float) else str(value)
            badge_class = "wt-badge" if variant == 'WT' else ""
            
            html += f"""
                    <div class="value-row">
                        <span class="variant-name">
                            <span class="variant-badge {badge_class}">{variant}</span>
                        </span>
                        <div>
                            <span class="value">{value_str}</span>
                            {delta_html}
                        </div>
                    </div>
"""
        
        html += """
                </div>
            </div>
"""
    
    html += """
        </div>
        
        <h2>Detailed Comparison Table</h2>
        <table class="comparison-table">
            <thead>
                <tr>
                    <th>Metric</th>
"""
    
    for variant in variants:
        badge_class = "wt-badge" if variant == 'WT' else ""
        html += f'<th><span class="variant-badge {badge_class}">{variant}</span></th>'
    
    if has_wt and len(variants) > 1:
        html += '<th>Max Δ from WT</th>'
    
    html += """
                </tr>
            </thead>
            <tbody>
"""
    
    for metric in metrics:
        values = summary_data[metric]
        html += f"""
                <tr>
                    <td><strong>{metric.replace('_', ' ').title()}</strong></td>
"""
        
        for variant in variants:
            value = values[variant]
            value_str = f"{value:.2f}" if isinstance(value, float) else str(value)
            html += f'<td>{value_str}</td>'
        
        # Add max delta column if WT exists
        if has_wt and len(variants) > 1:
            wt_value = values.get('WT', 0)
            deltas = [abs(values[v] - wt_value) for v in variants if v != 'WT']
            max_delta = max(deltas) if deltas else 0
            max_delta_str = f"{max_delta:.2f}" if isinstance(max_delta, float) else str(max_delta)
            html += f'<td>{max_delta_str}</td>'
        
        html += """
                </tr>
"""
    
    html += """
            </tbody>
        </table>
    </div>
</body>
</html>
"""
    
    return html


def generate_residue_comparison_table(
    results: Dict[str, Dict],
    title: str = "Residue-by-Residue Comparison",
    show_identical: bool = False,
    filter_interface_only: bool = False
) -> str:
    """
    Generate interactive HTML table comparing all residues across variants.
    
    Shows overall differences in the table, with detailed per-variant information
    in interactive tooltips on hover.
    
    Args:
        results: Dictionary with variant data
        title: Title for the visualization
        show_identical: Show residues with no changes
        filter_interface_only: Only show interface residues
        
    Returns:
        HTML string with interactive comparison table
    """
    residue_comp = compare_residues(results)
    variants = list(results.keys())
    
    html = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>{title}</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f5f5f5;
        }}
        .container {{
            max-width: 1600px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2c3e50;
            margin-bottom: 10px;
        }}
        .controls {{
            margin: 20px 0;
            padding: 15px;
            background: #f8f9fa;
            border-radius: 6px;
            display: flex;
            gap: 20px;
            align-items: center;
            flex-wrap: wrap;
        }}
        .controls label {{
            display: flex;
            align-items: center;
            gap: 8px;
            font-size: 14px;
        }}
        .controls input[type="checkbox"] {{
            width: 18px;
            height: 18px;
        }}
        .controls input[type="text"] {{
            padding: 8px 12px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 14px;
        }}
        .comparison-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            font-size: 13px;
        }}
        .comparison-table th {{
            background: #34495e;
            color: white;
            padding: 12px 8px;
            text-align: left;
            font-weight: 600;
            position: sticky;
            top: 0;
            z-index: 10;
            cursor: pointer;
            user-select: none;
        }}
        .comparison-table th:hover {{
            background: #2c3e50;
        }}
        .comparison-table td {{
            padding: 10px 8px;
            border-bottom: 1px solid #ecf0f1;
        }}
        .comparison-table tbody tr:hover {{
            background: #f8f9fa;
        }}
        .residue-cell {{
            display: flex;
            align-items: center;
            gap: 5px;
        }}
        .res-type {{
            font-weight: 600;
            font-family: monospace;
        }}
        .mutation-badge {{
            display: inline-block;
            padding: 2px 6px;
            border-radius: 3px;
            font-size: 11px;
            font-weight: 600;
            background: #f39c12;
            color: white;
        }}
        .deletion-badge {{
            background: #e74c3c;
        }}
        .insertion-badge {{
            background: #3498db;
        }}
        .interface-badge {{
            display: inline-block;
            padding: 2px 6px;
            border-radius: 3px;
            font-size: 10px;
            background: #2ecc71;
            color: white;
        }}
        .value-cell {{
            text-align: right;
            font-family: monospace;
        }}
        .value-positive {{
            color: #27ae60;
        }}
        .value-negative {{
            color: #c0392b;
        }}
        .row-mutation {{
            background: #fff9e6 !important;
        }}
        .row-deletion {{
            background: #ffe6e6 !important;
        }}
        .row-insertion {{
            background: #e6f3ff !important;
        }}
        .missing-cell {{
            color: #bdc3c7;
            font-style: italic;
        }}
        .sort-icon {{
            margin-left: 5px;
            font-size: 10px;
        }}
        .stats-bar {{
            display: flex;
            gap: 20px;
            margin: 20px 0;
            padding: 15px;
            background: #ecf0f1;
            border-radius: 6px;
        }}
        .stat-item {{
            display: flex;
            flex-direction: column;
        }}
        .stat-value {{
            font-size: 24px;
            font-weight: bold;
            color: #2c3e50;
        }}
        .stat-label {{
            font-size: 12px;
            color: #7f8c8d;
            text-transform: uppercase;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>{title}</h1>
        
        <div class="stats-bar">
            <div class="stat-item">
                <span class="stat-value" id="total-residues">{len(residue_comp)}</span>
                <span class="stat-label">Total Positions</span>
            </div>
            <div class="stat-item">
                <span class="stat-value" id="mutations-count">
                    {sum(1 for r in residue_comp.values() if r['comparison_type'] == 'mutation')}
                </span>
                <span class="stat-label">Mutations</span>
            </div>
            <div class="stat-item">
                <span class="stat-value" id="insertions-count">
                    {sum(1 for r in residue_comp.values() if r['comparison_type'] == 'insertion')}
                </span>
                <span class="stat-label">Insertions</span>
            </div>
            <div class="stat-item">
                <span class="stat-value" id="deletions-count">
                    {sum(1 for r in residue_comp.values() if r['comparison_type'] == 'deletion')}
                </span>
                <span class="stat-label">Deletions</span>
            </div>
        </div>
        
        <div class="controls">
            <label>
                <input type="checkbox" id="show-identical" {'checked' if show_identical else ''}>
                Show Identical Residues
            </label>
            <label>
                <input type="checkbox" id="interface-only" {'checked' if filter_interface_only else ''}>
                Interface Residues Only
            </label>
            <label>
                Search:
                <input type="text" id="search-box" placeholder="Position, residue...">
            </label>
        </div>
        
        <table class="comparison-table" id="comparison-table">
            <thead>
                <tr>
                    <th onclick="sortTable(0)">Position <span class="sort-icon">▼</span></th>
                    <th onclick="sortTable(1)">Type <span class="sort-icon"></span></th>
"""
    
    # Add columns for each variant
    for i, variant in enumerate(variants, start=2):
        html += f'<th onclick="sortTable({i})">{"🔵 " if variant == "WT" else ""}{variant} <span class="sort-icon"></span></th>'
    
    # Add property columns
    prop_col_start = len(variants) + 2
    html += f"""
                    <th onclick="sortTable({prop_col_start})">SASA Δ <span class="sort-icon"></span></th>
                    <th onclick="sortTable({prop_col_start + 1})">Interface (atom-atom contact) <span class="sort-icon"></span></th>
                    <th onclick="sortTable({prop_col_start + 2})">VDW Energy <span class="sort-icon"></span></th>
                    <th onclick="sortTable({prop_col_start + 3})">#Bonds <span class="sort-icon"></span></th>
                </tr>
            </thead>
            <tbody>
"""
    
    # Generate table rows
    for position, data in residue_comp.items():
        comp_type = data['comparison_type']
        
        # Skip identical if not showing
        if not show_identical and comp_type == 'identical':
            # Check if any property changed
            has_property_change = False
            wt_data = data['variants'].get('WT')
            if wt_data and wt_data['present']:
                for variant in variants:
                    if variant == 'WT':
                        continue
                    var_data = data['variants'].get(variant)
                    if var_data and var_data['present']:
                        if abs(var_data['sasa']['buried_fraction'] - wt_data['sasa']['buried_fraction']) > 0.05:
                            has_property_change = True
                            break
                        if abs(var_data['vdw_energy'] - wt_data['vdw_energy']) > 1.0:
                            has_property_change = True
                            break
            
            if not has_property_change:
                continue
        
        # Filter interface if requested
        if filter_interface_only:
            is_interface = any(
                v['interface']['is_interface'] 
                for v in data['variants'].values() 
                if v['present'] and v['interface']
            )
            if not is_interface:
                continue
        
        # Row class based on type
        row_class = f"row-{comp_type}" if comp_type != 'identical' else ""
        
        # Type badge
        type_badge = ""
        if comp_type == 'mutation':
            type_badge = '<span class="mutation-badge">MUT</span>'
        elif comp_type == 'deletion':
            type_badge = '<span class="deletion-badge">DEL</span>'
        elif comp_type == 'insertion':
            type_badge = '<span class="insertion-badge">INS</span>'
        
        html += f"""
                <tr class="{row_class}">
                    <td><strong>{position}</strong></td>
                    <td>{type_badge}</td>
"""
        
        # Variant columns
        for variant in variants:
            var_data = data['variants'][variant]
            if var_data['present']:
                res_type = var_data['residue_type']
                interface_badge = '<span class="interface-badge">I</span>' if var_data['interface']['is_interface'] else ''
                
                # Check for mutation
                mutation_info = data.get('mutation_status', {}).get(variant, {})
                if 'mutation' in mutation_info:
                    mut_text = f" ({mutation_info['mutation']})"
                    html += f'<td><div class="residue-cell"><span class="res-type">{res_type}</span>{mut_text} {interface_badge}</div></td>'
                else:
                    html += f'<td><div class="residue-cell"><span class="res-type">{res_type}</span> {interface_badge}</div></td>'
            else:
                html += '<td class="missing-cell">—</td>'
        
        # Property columns - show range or specific values
        # SASA
        sasa_values = [v['sasa']['buried_fraction'] for v in data['variants'].values() if v['present']]
        if sasa_values:
            min_sasa = min(sasa_values)
            max_sasa = max(sasa_values)
            if abs(max_sasa - min_sasa) > 0.01:
                html += f'<td class="value-cell">{min_sasa:.2f} - {max_sasa:.2f}</td>'
            else:
                html += f'<td class="value-cell">{min_sasa:.2f}</td>'
        else:
            html += '<td class="missing-cell">—</td>'
        
        # Interface contacts
        contact_values = [v['interface']['contact_count'] for v in data['variants'].values() if v['present'] and v['interface']['is_interface']]
        if contact_values:
            min_contacts = min(contact_values)
            max_contacts = max(contact_values)
            if min_contacts != max_contacts:
                html += f'<td class="value-cell">{min_contacts} - {max_contacts}</td>'
            else:
                html += f'<td class="value-cell">{min_contacts}</td>'
        else:
            html += '<td class="missing-cell">—</td>'
        
        # VDW Energy
        energy_values = [v['vdw_energy'] for v in data['variants'].values() if v['present']]
        if energy_values:
            min_energy = min(energy_values)
            max_energy = max(energy_values)
            if abs(max_energy - min_energy) > 0.1:
                html += f'<td class="value-cell">{min_energy:.2f} - {max_energy:.2f}</td>'
            else:
                html += f'<td class="value-cell">{min_energy:.2f}</td>'
        else:
            html += '<td class="missing-cell">—</td>'
        
        # Bonds
        bond_values = [v['bonds']['counts']['total'] for v in data['variants'].values() if v['present']]
        if bond_values:
            min_bonds = min(bond_values)
            max_bonds = max(bond_values)
            if min_bonds != max_bonds:
                html += f'<td class="value-cell">{min_bonds} - {max_bonds}</td>'
            else:
                html += f'<td class="value-cell">{min_bonds}</td>'
        else:
            html += '<td class="missing-cell">—</td>'
        
        html += """
                </tr>
"""
    
    html += """
            </tbody>
        </table>
    </div>
    
    <script>
        // Filter controls
        document.getElementById('show-identical').addEventListener('change', filterTable);
        document.getElementById('interface-only').addEventListener('change', filterTable);
        document.getElementById('search-box').addEventListener('input', filterTable);
        
        function filterTable() {
            const showIdentical = document.getElementById('show-identical').checked;
            const interfaceOnly = document.getElementById('interface-only').checked;
            const searchTerm = document.getElementById('search-box').value.toLowerCase();
            
            const rows = document.querySelectorAll('#comparison-table tbody tr');
            let visibleCount = 0;
            
            rows.forEach(row => {
                let show = true;
                
                // Check if row has mutation/insertion/deletion class
                const hasChange = row.classList.contains('row-mutation') || 
                                 row.classList.contains('row-deletion') || 
                                 row.classList.contains('row-insertion');
                
                if (!showIdentical && !hasChange) {
                    show = false;
                }
                
                // Check interface filter
                if (interfaceOnly) {
                    const hasInterface = row.textContent.includes('I');
                    if (!hasInterface) show = false;
                }
                
                // Check search
                if (searchTerm && !row.textContent.toLowerCase().includes(searchTerm)) {
                    show = false;
                }
                
                row.style.display = show ? '' : 'none';
                if (show) visibleCount++;
            });
        }
        
        // Table sorting
        let sortDirection = {};
        
        function sortTable(columnIndex) {
            const table = document.getElementById('comparison-table');
            const tbody = table.querySelector('tbody');
            const rows = Array.from(tbody.querySelectorAll('tr'));
            
            // Toggle sort direction
            sortDirection[columnIndex] = !sortDirection[columnIndex];
            const ascending = sortDirection[columnIndex];
            
            // Update sort icons
            document.querySelectorAll('.sort-icon').forEach(icon => icon.textContent = '');
            const header = table.querySelectorAll('th')[columnIndex];
            header.querySelector('.sort-icon').textContent = ascending ? '▲' : '▼';
            
            rows.sort((a, b) => {
                const aCell = a.cells[columnIndex].textContent.trim();
                const bCell = b.cells[columnIndex].textContent.trim();
                
                // Try numeric comparison
                const aNum = parseFloat(aCell.replace(/[^0-9.-]/g, ''));
                const bNum = parseFloat(bCell.replace(/[^0-9.-]/g, ''));
                
                if (!isNaN(aNum) && !isNaN(bNum)) {
                    return ascending ? aNum - bNum : bNum - aNum;
                }
                
                // String comparison
                return ascending ? 
                    aCell.localeCompare(bCell) : 
                    bCell.localeCompare(aCell);
            });
            
            rows.forEach(row => tbody.appendChild(row));
        }
        
        // Initialize
        filterTable();
    </script>
</body>
</html>
"""
    
    return html


def generate_comparison_dashboard(
    results: Dict[str, Dict],
    title: str = "Protein Complex Comparison Dashboard",
    output_prefix: Optional[str] = None,
    include_3d: bool = True
) -> Dict[str, str]:
    """
    Generate a comprehensive comparison dashboard with all visualizations.
    
    Args:
        results: Dictionary with variant data
        title: Dashboard title
        output_prefix: If provided, saves files with this prefix
        include_3d: Whether to include 3D structure visualization
        
    Returns:
        Dictionary with HTML strings for each visualization type
    """
    
    # Generate all visualizations
    visualizations = {
        'summary': generate_summary_comparison_html(results, title=f"{title} - Summary"),
        'residues': generate_residue_comparison_table(results, title=f"{title} - Residues"),
    }
    
    # Add 3D visualization if requested
    if include_3d:
        visualizations['structure_3d'] = generate_3d_structure_comparison(
            results, 
            title=f"{title} - 3D Structure"
        )
    
    # Save files if output_prefix provided
    if output_prefix:
        for viz_type, html in visualizations.items():
            filename = f"{output_prefix}_{viz_type}.html"
            with open(filename, 'w') as f:
                f.write(html)
            print(f"Saved {viz_type} visualization to {filename}")
    
    return visualizations 



def visualize_comparison(
    results: Dict[str, Dict],
    viz_type: str = 'all',
    output_file: Optional[str] = None,
    **kwargs
) -> Dict[str, str]:
    """
    Main function to generate comparison visualizations.
    
    Args:
        results: Dictionary with variant comparison data
        viz_type: Type of visualization ('summary', 'residues', '3d', 'all')
        output_file: Output file path (without extension)
        **kwargs: Additional arguments passed to visualization functions
        
    Returns:
        Dictionary mapping visualization types to HTML strings or file paths
    """
    
    if viz_type == 'summary':
        html = generate_summary_comparison_html(results, **kwargs)
        if output_file:
            with open(f"{output_file}.html", 'w') as f:
                f.write(html)
            return {'summary': f"{output_file}.html"}
        return {'summary': html}
    
    elif viz_type == 'residues':
        html = generate_residue_comparison_table(results, **kwargs)
        if output_file:
            with open(f"{output_file}.html", 'w') as f:
                f.write(html)
            return {'residues': f"{output_file}.html"}
        return {'residues': html}
    
    elif viz_type == '3d' or viz_type == 'structure_3d':
        html = generate_3d_structure_comparison(results, output_file=output_file, **kwargs)
        if output_file:
            return {'structure_3d': output_file}
        return {'structure_3d': html}
    
    elif viz_type == 'all':
        return generate_comparison_dashboard(results, output_prefix=output_file, **kwargs)
    
    else:
        raise ValueError(f"Unknown viz_type: {viz_type}. Use 'summary', 'residues', '3d', or 'all'")



def generate_3d_structure_comparison(
    results: Dict[str, Dict],
    pdb_files: Tuple[str, str],
    output_file: Optional[str] = None
) -> str:
    """
    Generate interactive 3D visualization of multiple protein complexes using 3Dmol.js.
    
    This creates a real molecular viewer showing actual protein structures with toggleable
    property coloring. Each structure is loaded from PDB files and can be colored by
    different analytical properties.
    
    Args:
        results: Dictionary with variant names as keys, each containing:
            - 'residues': Dict of residue data with properties (SASA, interface, VdW, etc.)
            - 'pdb_path': Optional path to PDB file (can also use pdb_files parameter)
        pdb_files: Optional dict mapping variant names to PDB file paths
            If not provided, will try to use 'pdb_path' from results
        title: Title for the visualization
        output_file: Optional file path to save the HTML
        
    Returns:
        HTML string with interactive 3D molecular visualization
        
    Example:
        >>> results = comparator.analyze()
        >>> pdb_files = [('WT', 'wt.pdb'), ('Variant1', 'var1.pdb')] 
        >>> html = generate_3d_structure_comparison(
        ...     results, 
        ...     pdb_files=pdb_files,
        ...     output_file='comparison_3d.html'
        ... )
    """
    
    variants = list(results.keys())
    
    # Read PDB file contents
    pdb_contents = {}
    for variant, pdb_path in pdb_files:
        try:
            with open(pdb_path, 'r') as f:
                pdb_contents[variant] = f.read()
        except Exception as e:
            print(f"Warning: Could not read PDB file for {variant}: {e}")
            pdb_contents[variant] = ""
    
    # Prepare residue property data
    property_data = {}
    for variant in variants:
        residues = results[variant]['residues']
        variant_props = {}
        
        for res_id, res_data in residues.items():
            # Parse residue ID: "SER A:32" -> chain, resseq
            parts = res_id.split()
            if len(parts) != 2:
                continue
                
            chain_pos = parts[1]  # "A:32"
            chain, position = chain_pos.split(':')
            
            # Extract numeric part of position (handle insertion codes)
            resseq = int(''.join(filter(str.isdigit, position)))
            
            # Get properties
            sasa_buried = res_data['sasa']['buried_fraction'] if res_data['sasa'] else 0
            interface_data = _parse_interface_data(res_data['interface'])
            is_interface = 1.0 if interface_data['is_interface'] else 0.0
            vdw_energy = res_data['vdw_energy'] if res_data['vdw_energy'] is not None else 0
            bonds_data = _parse_bonds(res_data['bonds'])
            total_contacts = bonds_data['counts']['total']
            
            # Store by chain:resseq key for easy lookup
            key = f"{chain}:{resseq}"
            variant_props[key] = {
                'sasa': sasa_buried,
                'interface': is_interface,
                'vdw': vdw_energy,
                'contacts': total_contacts,
                'res_id': res_id
            }
        
        property_data[variant] = variant_props
    
    # Generate HTML with 3Dmol.js
    html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    
    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            padding: 20px;
        }}
        
        .container {{
            max-width: 1800px;
            margin: 0 auto;
            background: white;
            border-radius: 16px;
            box-shadow: 0 20px 60px rgba(0,0,0,0.3);
            overflow: hidden;
        }}
        
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px 40px;
            text-align: center;
        }}
        
        h1 {{
            font-size: 32px;
            margin-bottom: 10px;
            font-weight: 600;
        }}
        
        .subtitle {{
            font-size: 16px;
            opacity: 0.95;
            font-weight: 300;
        }}
        
        .main-content {{
            display: flex;
            height: calc(100vh - 100px);
            min-height: 600px;
        }}
        
        .controls-panel {{
            width: 320px;
            background: #f8f9fa;
            padding: 25px;
            overflow-y: auto;
            border-right: 1px solid #e0e0e0;
        }}
        
        .viewer-container {{
            flex: 1;
            position: relative;
            background: #000;
        }}
        
        #viewport {{
            width: 100%;
            height: 100%;
            position: relative;
        }}
        
        .control-section {{
            margin-bottom: 25px;
        }}
        
        .control-label {{
            font-size: 13px;
            font-weight: 600;
            color: #333;
            margin-bottom: 10px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        
        select, button {{
            width: 100%;
            padding: 12px;
            border: 2px solid #e0e0e0;
            border-radius: 8px;
            font-size: 14px;
            background: white;
            cursor: pointer;
            transition: all 0.3s;
            font-family: inherit;
        }}
        
        select:hover, select:focus {{
            border-color: #667eea;
            outline: none;
        }}
        
        button {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border: none;
            font-weight: 600;
            margin-top: 10px;
        }}
        
        button:hover {{
            transform: translateY(-2px);
            box-shadow: 0 6px 20px rgba(102, 126, 234, 0.4);
        }}
        
        button:active {{
            transform: translateY(0);
        }}
        
        .variant-checkboxes {{
            display: flex;
            flex-direction: column;
            gap: 10px;
        }}
        
        .checkbox-item {{
            display: flex;
            align-items: center;
            padding: 10px;
            background: white;
            border: 2px solid #e0e0e0;
            border-radius: 8px;
            cursor: pointer;
            transition: all 0.3s;
        }}
        
        .checkbox-item:hover {{
            border-color: #667eea;
            background: #f0f4ff;
        }}
        
        .checkbox-item input[type="checkbox"] {{
            width: 20px;
            height: 20px;
            margin-right: 10px;
            cursor: pointer;
        }}
        
        .checkbox-item label {{
            cursor: pointer;
            font-size: 14px;
            font-weight: 500;
            flex: 1;
        }}
        
        .legend {{
            background: white;
            padding: 15px;
            border-radius: 8px;
            border: 2px solid #e0e0e0;
            margin-top: 15px;
        }}
        
        .legend-title {{
            font-size: 12px;
            font-weight: 600;
            color: #666;
            margin-bottom: 10px;
            text-transform: uppercase;
        }}
        
        .legend-gradient {{
            height: 20px;
            border-radius: 4px;
            margin: 10px 0;
        }}
        
        .legend-labels {{
            display: flex;
            justify-content: space-between;
            font-size: 11px;
            color: #666;
        }}
        
        .info-badge {{
            background: #e3f2fd;
            color: #1976d2;
            padding: 8px 12px;
            border-radius: 6px;
            font-size: 12px;
            margin-top: 15px;
            line-height: 1.6;
        }}
        
        .residue-info-panel {{
            position: absolute;
            bottom: 20px;
            right: 20px;
            width: 350px;
            max-height: 400px;
            background: white;
            border-radius: 10px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.3);
            padding: 20px;
            display: none;
            overflow-y: auto;
            z-index: 1000;
        }}
        
        .residue-info-panel.visible {{
            display: block;
        }}
        
        .residue-info-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 15px;
            padding-bottom: 10px;
            border-bottom: 2px solid #667eea;
        }}
        
        .residue-info-title {{
            font-size: 16px;
            font-weight: 700;
            color: #333;
        }}
        
        .close-btn {{
            background: #ff6b6b;
            color: white;
            border: none;
            width: 24px;
            height: 24px;
            border-radius: 50%;
            cursor: pointer;
            font-size: 16px;
            line-height: 1;
            padding: 0;
        }}
        
        .close-btn:hover {{
            background: #ff5252;
        }}
        
        .residue-info-content {{
            font-size: 13px;
        }}
        
        .info-row {{
            display: flex;
            justify-content: space-between;
            padding: 8px 0;
            border-bottom: 1px solid #f0f0f0;
        }}
        
        .info-label {{
            font-weight: 600;
            color: #666;
        }}
        
        .info-value {{
            color: #333;
            text-align: right;
        }}
        
        .info-section {{
            margin-top: 15px;
            padding-top: 10px;
            border-top: 2px solid #f0f0f0;
        }}
        
        .info-section-title {{
            font-weight: 700;
            color: #667eea;
            margin-bottom: 8px;
            font-size: 12px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        
        .loading {{
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            color: white;
            font-size: 18px;
            text-align: center;
        }}
        
        .loading-spinner {{
            border: 4px solid rgba(255,255,255,0.3);
            border-top: 4px solid white;
            border-radius: 50%;
            width: 40px;
            height: 40px;
            animation: spin 1s linear infinite;
            margin: 0 auto 15px;
        }}
        
        @keyframes spin {{
            0% {{ transform: rotate(0deg); }}
            100% {{ transform: rotate(360deg); }}
        }}
        
        .style-buttons {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 8px;
            margin-top: 10px;
        }}
        
        .style-button {{
            padding: 8px;
            font-size: 12px;
            margin: 0;
        }}
        
        .property-info {{
            font-size: 12px;
            color: #666;
            margin-top: 8px;
            padding: 8px;
            background: #f0f4ff;
            border-radius: 6px;
            line-height: 1.5;
        }}
    </style>
</head>
<body>
    <div class="container">
        <!--div class="header">
            <h1>Interactive Molecular Visualization with Property Coloring</h1>
        </div-->
        
        <div class="main-content">
            <div class="controls-panel">
                <div class="control-section">
                    <div class="control-label">🎨 Color Property</div>
                    <select id="property-select">
                        <option value="default">By Variant (distinct colors)</option>
                        <option value="sasa">SASA (Buried Fraction)</option>
                        <option value="interface">Interface Residues</option>
                        <option value="vdw">VdW Energy</option>
                        <option value="contacts">Total Contacts</option>
                        <option value="bfactor">B-factor</option>
                    </select>
                    <div id="property-info" class="property-info"></div>
                </div>
                
                <div class="control-section">
                    <div class="control-label">🔬 Representation</div>
                    <select id="style-select">
                        <option value="cartoon">Cartoon</option>
                        <option value="stick">Stick</option>
                        <option value="sphere">Sphere</option>
                        <option value="line">Line</option>
                        <option value="surface">Surface</option>
                    </select>
                    <div class="style-buttons">
                        <button class="style-button" onclick="toggleSideChains()">Toggle Side Chains</button>
                        <button class="style-button" onclick="toggleChainOutlines()">Chain Outlines</button>
                    </div>
                </div>
                
                <div class="control-section">
                    <div class="control-label">📊 Display Options</div>
                    <select id="display-select">
                        <option value="all">Show All Structures</option>
                        <option value="overlay">Overlay (aligned)</option>
                        <option value="grid">Grid Layout</option>
                    </select>
                </div>
                
                <div class="control-section">
                    <div class="control-label">🔍 Difference Highlighting</div>
                    <select id="difference-mode">
                        <option value="none">None</option>
                        <option value="mutations">Highlight Mutations Only</option>
                        <option value="rmsd">Color by RMSD/Distance</option>
                        <option value="changed">Show Changed Regions Only</option>
                    </select>
                    <div style="margin-top: 10px;">
                        <label style="display: flex; align-items: center; gap: 5px; font-size: 12px;">
                            <input type="checkbox" id="fade-similar" onchange="updateVisualization()">
                            <span>Fade Similar Regions</span>
                        </label>
                    </div>
                </div>
                
                <div class="control-section">
                    <div class="control-label">🧬 Variants</div>
                    <div class="variant-checkboxes" id="variant-checkboxes">
"""
    
    # Add variant checkboxes
    for i, variant in enumerate(variants):
        checked = "checked" if i == 0 else ""
        html += f"""
                        <div class="checkbox-item">
                            <input type="checkbox" id="variant-{variant}" value="{variant}" {checked} onchange="updateVisualization()">
                            <label for="variant-{variant}">{variant}</label>
                        </div>
"""
    
    html += """
                    </div>
                </div>
                
                <div class="control-section">
                    <button onclick="updateVisualization()">🔄 Update Visualization</button>
                    <button onclick="resetView()">🔄 Reset View</button>
                    <button onclick="downloadImage()">📸 Save Image</button>
                </div>
                
                <div class="legend" id="legend">
                    <div class="legend-title">Color Scale</div>
                    <div id="legend-content">Select a property to see color scale</div>
                </div>
                
                <div class="info-badge">
                    💡 <strong>Tip:</strong> Click and drag to rotate, scroll to zoom, right-click and drag to pan. <strong>Click on any residue to see detailed information!</strong> Use "Chain Outlines" button to distinguish chains when using property coloring.
                </div>
            </div>
            
            <div class="viewer-container">
                <div id="viewport"></div>
                <div class="loading" id="loading">
                    <div class="loading-spinner"></div>
                    Loading structures...
                </div>
                
                <!-- Residue Information Panel -->
                <div class="residue-info-panel" id="residue-info-panel">
                    <div class="residue-info-header">
                        <div class="residue-info-title" id="residue-title">Residue Info</div>
                        <button class="close-btn" onclick="closeResidueInfo()">×</button>
                    </div>
                    <div class="residue-info-content" id="residue-info-content">
                        Click on a residue to see details
                    </div>
                </div>
            </div>
        </div>
    </div>
    
    <script>
"""
    
    # Embed PDB data and property data as JavaScript
    html += f"""
        // PDB data
        const pdbData = {json.dumps(pdb_contents, indent=2)};
        
        // Property data for each variant
        const propertyData = {json.dumps(property_data, indent=2)};
        
        // Variant list
        const variants = {json.dumps(variants)};
        
        // Initialize 3Dmol viewer
        let viewer = null;
        let showSideChains = false;
        let showChainOutlines = false;
        
        // Cache for difference detection (avoid recalculating)
        let cachedMutations = null;
        let cachedChangedResidues = null;
        let lastDifferenceMode = 'none';
        
        // Define distinct colors for each variant
        const variantColors = [
            '#FF6B6B', // Red
            '#52B788',  // Green
            '#45B7D1', // Blue
            '#BB8FCE', // Purple
            '#F7DC6F', // Yellow
            '#4ECDC4', // Teal
            '#FFA07A', // Light Salmon
            '#98D8C8', // Mint
            '#85C1E2', // Sky Blue
            '#F8B739' // Orange
            
        ];
        
        // Map each variant to a color
        const variantColorMap = {{}};
        variants.forEach((variant, idx) => {{
            variantColorMap[variant] = variantColors[idx % variantColors.length];
        }});
        
        // Color scales configuration
        const colorScales = {{
            sasa: {{
                name: 'SASA Buried Fraction',
                min: 0,
                max: 1,
                colors: ['#2166ac', '#67a9cf', '#d1e5f0', '#fddbc7', '#ef8a62', '#b2182b'],
                info: 'Blue = buried (0), Red = exposed (1)',
                getLabel: (v) => v.toFixed(2)
            }},
            interface: {{
                name: 'Interface Residues',
                min: 0,
                max: 1,
                colors: ['#cccccc', '#ff6b6b'],
                info: 'Gray = non-interface, Red = interface',
                getLabel: (v) => v === 1 ? 'Interface' : 'Core'
            }},
            vdw: {{
                name: 'VdW Energy (kcal/mol)',
                min: -15,
                max: 0,
                colors: ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0', '#f7f7f7', '#fddbc7'],
                info: 'Blue = favorable (negative), White = neutral',
                getLabel: (v) => v.toFixed(1)
            }},
            contacts: {{
                name: 'Total Contacts',
                min: 0,
                max: 20,
                colors: ['#ffffcc', '#ffeda0', '#fed976', '#feb24c', '#fd8d3c', '#fc4e2a', '#e31a1c', '#b10026'],
                info: 'Yellow = few contacts, Red = many contacts',
                getLabel: (v) => Math.round(v)
            }},
            bfactor: {{
                name: 'B-factor',
                min: 0,
                max: 100,
                colors: ['#0571b0', '#92c5de', '#f7f7f7', '#f4a582', '#ca0020'],
                info: 'Blue = low B-factor (rigid), Red = high B-factor (flexible)',
                getLabel: (v) => Math.round(v)
            }}
        }};
        
        // Initialize on load
        window.addEventListener('DOMContentLoaded', function() {{
            const element = document.getElementById('viewport');
            const config = {{ backgroundColor: 'white' }};
            viewer = $3Dmol.createViewer(element, config);
            
            // Initial visualization (click handler will be set up in updateVisualization)
            updateVisualization();
            
            // Event listeners
            document.getElementById('property-select').addEventListener('change', updateVisualization);
            document.getElementById('style-select').addEventListener('change', updateVisualization);
            document.getElementById('display-select').addEventListener('change', updateVisualization);
            document.getElementById('difference-mode').addEventListener('change', updateVisualization);
        }});
        
        function getSelectedVariants() {{
            const selected = [];
            variants.forEach(variant => {{
                const checkbox = document.getElementById(`variant-${{variant}}`);
                if (checkbox && checkbox.checked) {{
                    selected.push(variant);
                }}
            }});
            return selected;
        }}
        
        function identifyMutations(selectedVariants) {{
            // Identify positions where residue type differs between variants
            const mutations = new Map(); // key: chain:position, value: {{variants: {{variant: resType}}}}
            
            selectedVariants.forEach(variant => {{
                const props = propertyData[variant];
                if (!props) return;
                
                Object.keys(props).forEach(posKey => {{
                    if (!mutations.has(posKey)) {{
                        mutations.set(posKey, {{ variants: {{}} }});
                    }}
                    // We'd need residue type from property data - for now mark all as potential differences
                    mutations.get(posKey).variants[variant] = props[posKey].res_id;
                }});
            }});
            
            // Filter to only positions that differ
            const realMutations = new Map();
            mutations.forEach((value, key) => {{
                const resTypes = new Set(Object.values(value.variants));
                if (resTypes.size > 1) {{
                    realMutations.set(key, value);
                }}
            }});
            
            return realMutations;
        }}
        
        function identifyChangedResidues(selectedVariants) {{
            // Identify residues with significant property changes
            const changed = new Map();
            
            if (selectedVariants.length < 2) return changed;
            
            // Compare first variant to others
            const reference = selectedVariants[0];
            const refProps = propertyData[reference];
            if (!refProps) return changed;
            
            Object.keys(refProps).forEach(posKey => {{
                let maxDiff = 0;
                
                selectedVariants.slice(1).forEach(variant => {{
                    const varProps = propertyData[variant];
                    if (!varProps || !varProps[posKey]) return;
                    
                    // Calculate property differences
                    const sasaDiff = Math.abs(refProps[posKey].sasa - varProps[posKey].sasa);
                    const vdwDiff = Math.abs(refProps[posKey].vdw - varProps[posKey].vdw);
                    const contactDiff = Math.abs(refProps[posKey].contacts - varProps[posKey].contacts);
                    
                    // Weighted difference score
                    const diff = sasaDiff * 2 + vdwDiff * 0.5 + contactDiff * 0.3;
                    maxDiff = Math.max(maxDiff, diff);
                }});
                
                if (maxDiff > 0.3) {{ // Threshold for "significant" change
                    changed.set(posKey, maxDiff);
                }}
            }});
            
            return changed;
        }}
        
        function applyVisualizationStyle(variant, modelIdx, property, style, differenceMode, fadeSimilar, mutations, changedResidues) {{
            let styleSpec = {{}};
            const variantColor = variantColorMap[variant];
            const props = propertyData[variant];
            
            // Determine which residues to highlight based on difference mode
            let highlightPositions = new Set();
            
            if (differenceMode === 'mutations') {{
                mutations.forEach((value, key) => highlightPositions.add(key));
            }} else if (differenceMode === 'changed') {{
                changedResidues.forEach((value, key) => highlightPositions.add(key));
            }}
            
            // Base style
            if (property === 'default') {{
                if (style === 'cartoon') {{
                    styleSpec.cartoon = {{ color: variantColor, clickable: true }};
                }} else if (style === 'stick') {{
                    styleSpec.stick = {{ radius: 0.15, color: variantColor, clickable: true }};
                }} else if (style === 'sphere') {{
                    styleSpec.sphere = {{ radius: 1.0, color: variantColor, clickable: true }};
                }} else if (style === 'line') {{
                    styleSpec.line = {{ color: variantColor, clickable: true }};
                }} else if (style === 'surface') {{
                    styleSpec.surface = {{ opacity: 0.8, color: variantColor, clickable: true }};
                }}
            }} else if (property === 'bfactor') {{
                if (style === 'cartoon') {{
                    styleSpec.cartoon = {{ 
                        colorscheme: {{ prop: 'b', gradient: 'roygb', min: 0, max: 100 }},
                        clickable: true
                    }};
                }} else {{
                    styleSpec[style] = {{ 
                        colorscheme: {{ prop: 'b', gradient: 'roygb', min: 0, max: 100 }},
                        clickable: true
                    }};
                }}
            }} else {{
                if (style === 'cartoon') {{
                    styleSpec.cartoon = {{ style: 'trace', clickable: true }};
                }} else if (style === 'stick') {{
                    styleSpec.stick = {{ radius: 0.15, clickable: true }};
                }} else if (style === 'sphere') {{
                    styleSpec.sphere = {{ radius: 1.0, clickable: true }};
                }} else if (style === 'line') {{
                    styleSpec.line = {{ clickable: true }};
                }} else if (style === 'surface') {{
                    styleSpec.surface = {{ opacity: 0.8, clickable: true }};
                }}
            }}
            
            // Apply base style to all
            viewer.setStyle({{ model: modelIdx }}, styleSpec);
            
            // Apply property coloring if not default
            if (property !== 'default' && property !== 'bfactor' && props) {{
                Object.keys(props).forEach(key => {{
                    const [chain, resseq] = key.split(':');
                    const value = props[key][property];
                    
                    if (value !== undefined) {{
                        const color = getColor(property, value, variant);
                        if (color) {{
                            const sel = {{ chain: chain, resi: parseInt(resseq), model: modelIdx }};
                            const coloredStyle = {{ ...styleSpec }};
                            coloredStyle[style] = {{ ...styleSpec[style], color: color, clickable: true }};
                            viewer.setStyle(sel, coloredStyle);
                        }}
                    }}
                }});
            }}
            
            // Apply difference highlighting - RESPECTING COLOR SCHEME
            if (differenceMode !== 'none' && props) {{
                if (differenceMode === 'changed') {{
                    // Fade or highlight based on whether residue changed
                    Object.keys(props).forEach(key => {{
                        const [chain, resseq] = key.split(':');
                        const isHighlighted = highlightPositions.has(key);
                        const sel = {{ chain: chain, resi: parseInt(resseq), model: modelIdx }};
                        
                        if (fadeSimilar && !isHighlighted) {{
                            // Fade non-changed residues
                            const fadedStyle = {{}};
                            fadedStyle[style] = {{ opacity: 0.15, clickable: true }};
                            viewer.setStyle(sel, fadedStyle);
                        }} else if (isHighlighted) {{
                            // Highlight changed residues - keep their property color
                            const highlightStyle = {{}};
                            
                            // Get the current color for this residue from property scheme
                            let currentColor = variantColor; // default fallback
                            if (property !== 'default' && property !== 'bfactor') {{
                                const value = props[key][property];
                                if (value !== undefined) {{
                                    const propColor = getColor(property, value, variant);
                                    if (propColor) currentColor = propColor;
                                }}
                            }}
                            
                            // Emphasize with thicker representation
                            if (style === 'cartoon') {{
                                highlightStyle.cartoon = {{ color: currentColor, thickness: 1.5, clickable: true }};
                            }} else if (style === 'sphere') {{
                                highlightStyle.sphere = {{ radius: 1.3, color: currentColor, clickable: true }};
                            }} else {{
                                highlightStyle[style] = {{ color: currentColor, clickable: true }};
                            }}
                            
                            // Add stick outline for extra emphasis
                            highlightStyle.stick = {{ radius: 0.3, color: currentColor, clickable: true }};
                            
                            viewer.setStyle(sel, highlightStyle);
                        }}
                    }});
                }} else if (differenceMode === 'mutations') {{
                    // Highlight mutation positions - keep property colors
                    highlightPositions.forEach(posKey => {{
                        const [chain, resseq] = posKey.split(':');
                        const sel = {{ chain: chain, resi: parseInt(resseq), model: modelIdx }};
                        
                        // Get the current color for this residue from property scheme
                        let currentColor = variantColor; // default fallback
                        if (property !== 'default' && property !== 'bfactor' && props[posKey]) {{
                            const value = props[posKey][property];
                            if (value !== undefined) {{
                                const propColor = getColor(property, value, variant);
                                if (propColor) currentColor = propColor;
                            }}
                        }}
                        
                        // Add emphasis while keeping the property color
                        const mutStyle = {{}};
                        mutStyle.stick = {{ radius: 0.4, color: currentColor, clickable: true }};
                        mutStyle.sphere = {{ radius: 1.8, color: currentColor, opacity: 0.8, clickable: true }};
                        viewer.setStyle(sel, mutStyle, {{ keepStyle: true }});
                    }});
                    
                    // Fade everything else if requested
                    if (fadeSimilar) {{
                        Object.keys(props).forEach(key => {{
                            if (!highlightPositions.has(key)) {{
                                const [chain, resseq] = key.split(':');
                                const sel = {{ chain: chain, resi: parseInt(resseq), model: modelIdx }};
                                const fadedStyle = {{}};
                                fadedStyle[style] = {{ opacity: 0.15, clickable: true }};
                                viewer.setStyle(sel, fadedStyle);
                            }}
                        }});
                    }}
                }}
            }}
        }}
        
        function getColor(property, value, variant) {{
            if (property === 'default' || !property) {{
                return null; // Use default coloring
            }}
            
            const scale = colorScales[property];
            if (!scale) return null;
            
            // Normalize value to [0, 1]
            const normalized = Math.max(0, Math.min(1, (value - scale.min) / (scale.max - scale.min)));
            
            // Get color from gradient
            const colors = scale.colors;
            const idx = normalized * (colors.length - 1);
            const lower = Math.floor(idx);
            const upper = Math.ceil(idx);
            const frac = idx - lower;
            
            if (lower === upper) {{
                return colors[lower];
            }}
            
            // Interpolate between colors
            return colors[lower]; // Simple version - just return lower color
        }}
        
        function updateVisualization() {{
            if (!viewer) return;
            
            const property = document.getElementById('property-select').value;
            const style = document.getElementById('style-select').value;
            const displayMode = document.getElementById('display-select').value;
            const differenceMode = document.getElementById('difference-mode').value;
            const fadeSimilar = document.getElementById('fade-similar').checked;
            const selectedVariants = getSelectedVariants();
            
            if (selectedVariants.length === 0) {{
                alert('Please select at least one variant to display');
                return;
            }}
            
            // Show loading
            document.getElementById('loading').style.display = 'block';
            
            // Clear viewer
            viewer.clear();
            
            // Only calculate differences if difference mode is active AND it changed
            let mutations = new Map();
            let changedResidues = new Map();
            
            if (differenceMode !== 'none') {{
                // Use cache if mode hasn't changed
                if (differenceMode === lastDifferenceMode && cachedMutations && cachedChangedResidues) {{
                    mutations = cachedMutations;
                    changedResidues = cachedChangedResidues;
                }} else {{
                    // Calculate fresh
                    if (differenceMode === 'mutations' || differenceMode === 'rmsd') {{
                        mutations = identifyMutations(selectedVariants);
                    }}
                    if (differenceMode === 'changed') {{
                        changedResidues = identifyChangedResidues(selectedVariants);
                    }}
                    // Cache results
                    cachedMutations = mutations;
                    cachedChangedResidues = changedResidues;
                    lastDifferenceMode = differenceMode;
                }}
            }} else {{
                // Clear cache when not using difference mode
                cachedMutations = null;
                cachedChangedResidues = null;
                lastDifferenceMode = 'none';
            }}
            
            // Calculate positions for grid layout
            const spacing = 60;
            const gridSize = Math.ceil(Math.sqrt(selectedVariants.length));
            
            selectedVariants.forEach((variant, idx) => {{
                if (!pdbData[variant]) {{
                    console.warn(`No PDB data for ${{variant}}`);
                    return;
                }}
                
                // Add model
                const model = viewer.addModel(pdbData[variant], 'pdb');
                
                // Position for grid layout
                if (displayMode === 'grid') {{
                    const row = Math.floor(idx / gridSize);
                    const col = idx % gridSize;
                    const offsetX = (col - gridSize / 2) * spacing;
                    const offsetY = (row - gridSize / 2) * spacing;
                    model.setCoordinates(model.getCoordinates().map(coord => ({{
                        x: coord.x + offsetX,
                        y: coord.y + offsetY,
                        z: coord.z
                    }})));
                }}
                
                // Apply styling - single call, no duplication
                applyVisualizationStyle(variant, idx, property, style, differenceMode, fadeSimilar, mutations, changedResidues);
            }});
            
            // Add side chains if requested
            if (showSideChains && property === 'default') {{
                viewer.setStyle({{ atom: ['CA', 'C', 'N', 'O'] }}, {{ stick: {{ radius: 0.2, clickable: true }} }}, {{ byres: false }});
            }}
            
            // Add chain outlines if requested (helps distinguish chains with property coloring)
            if (showChainOutlines) {{
                selectedVariants.forEach((variant, idx) => {{
                    if (!pdbData[variant]) return;
                    
                    // Get all unique chains in this variant
                    const modelData = viewer.getModel(idx);
                    if (!modelData) return;
                    
                    const atoms = modelData.selectedAtoms({{}});
                    const chains = [...new Set(atoms.map(a => a.chain))];
                    
                    // Define colors for chain outlines (darker, more saturated colors)
                    const chainOutlineColors = [
                        '#8B0000', // Dark Red
                        '#006064', // Dark Cyan  
                        '#1A237E', // Dark Blue
                        '#E65100', // Dark Orange
                        '#004D40', // Dark Teal
                        '#F57F17', // Dark Yellow
                        '#4A148C', // Dark Purple
                        '#01579B', // Dark Sky Blue
                        '#BF360C', // Dark Deep Orange
                        '#1B5E20'  // Dark Green
                    ];
                    
                    chains.forEach((chain, chainIdx) => {{
                        const outlineColor = chainOutlineColors[chainIdx % chainOutlineColors.length];
                        
                        // Add thin cartoon trace for chain outline
                        viewer.setStyle(
                            {{ chain: chain, model: idx }},
                            {{ cartoon: {{ color: outlineColor, thickness: 0.4, opacity: 0.7 }} }},
                            {{ keepStyle: true }}
                        );
                    }});
                }});
            }}
            
            // Set up click handler after styles are applied
            viewer.setClickable({{}}, true, function(atom, viewer, event, container) {{
                if (atom) {{
                    showResidueInfo(atom);
                    event.stopPropagation();
                }}
            }});
            
            viewer.zoomTo();
            viewer.render();
            
            // Hide loading
            document.getElementById('loading').style.display = 'none';
            
            // Update legend
            updateLegend(property, selectedVariants);
            
            // Update property info
            updatePropertyInfo(property);
        }}
        
        function updateLegend(property, selectedVariants) {{
            const legendContent = document.getElementById('legend-content');
            
            if (property === 'default') {{
                // Show variant color assignments
                let html = '<div style="font-size: 11px; font-weight: 600; margin-bottom: 8px;">Variant Colors</div>';
                selectedVariants.forEach(variant => {{
                    const color = variantColorMap[variant];
                    html += `
                        <div style="display: flex; align-items: center; margin: 6px 0;">
                            <div style="width: 24px; height: 16px; background: ${{color}}; border-radius: 3px; margin-right: 8px; border: 1px solid #ccc;"></div>
                            <span style="font-size: 12px;">${{variant}}</span>
                        </div>
                    `;
                }});
                legendContent.innerHTML = html;
                return;
            }}
            
            const scale = colorScales[property];
            if (!scale) return;
            
            const gradient = scale.colors.join(', ');
            
            legendContent.innerHTML = `
                <div style="font-size: 11px; font-weight: 600; margin-bottom: 8px;">${{scale.name}}</div>
                <div class="legend-gradient" style="background: linear-gradient(to right, ${{gradient}});"></div>
                <div class="legend-labels">
                    <span>${{scale.getLabel(scale.min)}}</span>
                    <span>${{scale.getLabel(scale.max)}}</span>
                </div>
                <div style="font-size: 11px; color: #666; margin-top: 8px;">${{scale.info}}</div>
            `;
        }}
        
        function updatePropertyInfo(property) {{
            const infoDiv = document.getElementById('property-info');
            
            const descriptions = {{
                'default': 'Each variant/complex is colored with a distinct color for easy identification.',
                'sasa': 'Solvent Accessible Surface Area - measures how buried each residue is in the structure.',
                'interface': 'Highlights residues at protein-protein interfaces (contacts with other chains).',
                'vdw': 'Van der Waals energy - negative values indicate favorable interactions.',
                'contacts': 'Number of atomic contacts each residue makes with neighboring residues.',
                'bfactor': 'Temperature/B-factor from PDB - indicates structural flexibility or disorder.'
            }};
            
            infoDiv.textContent = descriptions[property] || '';
        }}
        
        function toggleSideChains() {{
            showSideChains = !showSideChains;
            updateVisualization();
        }}
        
        function toggleChainOutlines() {{
            showChainOutlines = !showChainOutlines;
            updateVisualization();
        }}
        
        function resetView() {{
            if (viewer) {{
                viewer.zoomTo();
                viewer.render();
            }}
        }}
        
        function downloadImage() {{
            if (viewer) {{
                viewer.render();
                const imgData = viewer.pngURI();
                const link = document.createElement('a');
                link.href = imgData;
                link.download = '3d_structure_comparison.png';
                link.click();
            }}
        }}
        
        function showResidueInfo(atom) {{
            const panel = document.getElementById('residue-info-panel');
            const title = document.getElementById('residue-title');
            const content = document.getElementById('residue-info-content');
            
            // Get atom information
            const resname = atom.resn || 'Unknown';
            const chain = atom.chain || '?';
            const resi = atom.resi || '?';
            const atomName = atom.atom || '?';
            const elem = atom.elem || '?';
            const bfactor = atom.b !== undefined ? atom.b.toFixed(2) : 'N/A';
            
            // Get model index to identify variant
            const modelIndex = atom.model || 0;
            const selectedVariants = getSelectedVariants();
            const variant = selectedVariants[modelIndex] || 'Unknown';
            
            // Find residue properties from our data
            const resKey = `${{chain}}:${{resi}}`;
            const variantProps = propertyData[variant];
            const resProps = variantProps ? variantProps[resKey] : null;
            
            // Build title
            title.textContent = `${{resname}} ${{chain}}:${{resi}}`;
            
            // Build content
            let html = `
                <div class="info-row">
                    <span class="info-label">Variant:</span>
                    <span class="info-value" style="color: ${{variantColorMap[variant]}}; font-weight: 700;">${{variant}}</span>
                </div>
                <div class="info-row">
                    <span class="info-label">Residue Type:</span>
                    <span class="info-value">${{resname}}</span>
                </div>
                <div class="info-row">
                    <span class="info-label">Chain:</span>
                    <span class="info-value">${{chain}}</span>
                </div>
                <div class="info-row">
                    <span class="info-label">Position:</span>
                    <span class="info-value">${{resi}}</span>
                </div>
                <div class="info-row">
                    <span class="info-label">Atom Clicked:</span>
                    <span class="info-value">${{atomName}} (${{elem}})</span>
                </div>
                <div class="info-row">
                    <span class="info-label">B-factor:</span>
                    <span class="info-value">${{bfactor}}</span>
                </div>
            `;
            
            // Add analysis properties if available
            if (resProps) {{
                html += `
                    <div class="info-section">
                        <div class="info-section-title">Analysis Properties</div>
                        
                        <div class="info-row">
                            <span class="info-label">SASA (Buried):</span>
                            <span class="info-value">${{resProps.sasa.toFixed(2)}}</span>
                        </div>
                        
                        <div class="info-row">
                            <span class="info-label">Interface:</span>
                            <span class="info-value">${{resProps.interface === 1 ? '✓ Yes' : '✗ No'}}</span>
                        </div>
                        
                        <div class="info-row">
                            <span class="info-label">VdW Energy:</span>
                            <span class="info-value">${{resProps.vdw.toFixed(2)}} kcal/mol</span>
                        </div>
                        
                        <div class="info-row">
                            <span class="info-label">Total Contacts:</span>
                            <span class="info-value">${{resProps.contacts}}</span>
                        </div>
                    </div>
                `;
            }} else {{
                html += `
                    <div class="info-section">
                        <div class="info-section-title">Analysis Properties</div>
                        <div style="color: #999; font-style: italic; font-size: 12px;">
                            No analysis data available for this residue
                        </div>
                    </div>
                `;
            }}
            
            content.innerHTML = html;
            panel.classList.add('visible');
        }}
        
        function closeResidueInfo() {{
            const panel = document.getElementById('residue-info-panel');
            panel.classList.remove('visible');
        }}
    </script>
</body>
</html>
"""
    
    # Save to file if requested
    if output_file:
        with open(output_file, 'w') as f:
            f.write(html)
        print(f"Saved 3D visualization to {output_file}")
    
    return html




def generate_3d_sidebyside_comparison(
    results: Dict[str, Dict],
    pdb_files: Tuple[str, str],
    output_file: Optional[str] = None
) -> str:
    """
    Generate side-by-side 3D visualization of multiple protein complexes with synchronized controls.
    
    Each structure is displayed in its own viewer panel, but all viewers share the same
    styling, coloring, and viewing options for easy comparison.
    
    Args:
        results: Dictionary with variant names as keys, each containing:
            - 'residues': Dict of residue data with properties (SASA, interface, VdW, etc.)
        pdb_files: List of tuples mapping variant names to PDB file paths
        title: Title for the visualization
        output_file: Optional file path to save the HTML
        
    Returns:
        HTML string with side-by-side 3D molecular visualization
    """
    
    variants = list(results.keys())
    
    # Read PDB file contents
    pdb_contents = {}
    for variant, pdb_path in pdb_files:
        try:
            with open(pdb_path, 'r') as f:
                pdb_contents[variant] = f.read()
        except Exception as e:
            print(f"Warning: Could not read PDB file for {variant}: {e}")
            pdb_contents[variant] = ""
    
    # Prepare residue property data
    property_data = {}
    for variant in variants:
        residues = results[variant]['residues']
        variant_props = {}
        
        for res_id, res_data in residues.items():
            # Parse residue ID: "SER A:32" -> chain, resseq
            parts = res_id.split()
            if len(parts) != 2:
                continue
                
            chain_pos = parts[1]  # "A:32"
            chain, position = chain_pos.split(':')
            
            # Extract numeric part of position (handle insertion codes)
            resseq = int(''.join(filter(str.isdigit, position)))
            
            # Get properties
            sasa_buried = res_data['sasa']['buried_fraction'] if res_data['sasa'] else 0
            interface_data = _parse_interface_data(res_data['interface'])
            is_interface = 1.0 if interface_data['is_interface'] else 0.0
            vdw_energy = res_data['vdw_energy'] if res_data['vdw_energy'] is not None else 0
            bonds_data = _parse_bonds(res_data['bonds'])
            total_contacts = bonds_data['counts']['total']
            
            # Store by chain:resseq key for easy lookup
            key = f"{chain}:{resseq}"
            variant_props[key] = {
                'sasa': sasa_buried,
                'interface': is_interface,
                'vdw': vdw_energy,
                'contacts': total_contacts,
                'res_id': res_id
            }
        
        property_data[variant] = variant_props
    
    # Generate HTML with 3Dmol.js
    html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">

    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            padding: 20px;
        }}
        
        .container {{
            max-width: 2000px;
            margin: 0 auto;
            background: white;
            border-radius: 16px;
            box-shadow: 0 20px 60px rgba(0,0,0,0.3);
            overflow: hidden;
        }}
        
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px 40px;
            text-align: center;
        }}
        
        h1 {{
            font-size: 32px;
            margin-bottom: 10px;
            font-weight: 600;
        }}
        
        .subtitle {{
            font-size: 16px;
            opacity: 0.95;
            font-weight: 300;
        }}
        
        .main-content {{
            display: flex;
           
            min-height: 600px;
        }}
        
        .controls-panel {{
            width: 300px;
            background: #f8f9fa;
            padding: 25px;
            overflow-y: auto;
            border-right: 1px solid #e0e0e0;
        }}
        
        .viewers-container {{
            flex: 1;
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 2px;
            background: #e0e0e0;
            padding: 2px;
            overflow-y: auto;
        }}
        
        .viewer-panel {{
            background: #000;
            position: relative;
            min-height: 400px;
            display: flex;
            flex-direction: column;
        }}
        
        .viewer-header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 12px 16px;
            font-weight: 600;
            font-size: 14px;
            text-align: center;
            letter-spacing: 0.5px;
        }}
        
        .viewer-viewport {{
            flex: 1;
            position: relative;
            min-height: 350px;
        }}
        
        .control-section {{
            margin-bottom: 25px;
        }}
        
        .control-label {{
            font-size: 13px;
            font-weight: 600;
            color: #333;
            margin-bottom: 10px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        
        select, button {{
            width: 100%;
            padding: 12px;
            border: 2px solid #e0e0e0;
            border-radius: 8px;
            font-size: 14px;
            background: white;
            cursor: pointer;
            transition: all 0.3s;
            font-family: inherit;
        }}
        
        select:hover, select:focus {{
            border-color: #667eea;
            outline: none;
        }}
        
        button {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border: none;
            font-weight: 600;
            margin-top: 10px;
        }}
        
        button:hover {{
            transform: translateY(-2px);
            box-shadow: 0 6px 20px rgba(102, 126, 234, 0.4);
        }}
        
        button:active {{
            transform: translateY(0);
        }}
        
        .legend {{
            background: white;
            padding: 15px;
            border-radius: 8px;
            border: 2px solid #e0e0e0;
            margin-top: 15px;
        }}
        
        .legend-title {{
            font-size: 12px;
            font-weight: 600;
            color: #666;
            margin-bottom: 10px;
            text-transform: uppercase;
        }}
        
        .legend-gradient {{
            height: 20px;
            border-radius: 4px;
            margin: 10px 0;
        }}
        
        .legend-labels {{
            display: flex;
            justify-content: space-between;
            font-size: 11px;
            color: #666;
        }}
        
        .info-badge {{
            background: #e3f2fd;
            color: #1976d2;
            padding: 8px 12px;
            border-radius: 6px;
            font-size: 12px;
            margin-top: 15px;
            line-height: 1.6;
        }}
        
        .residue-info-panel {{
            position: fixed;
            bottom: 20px;
            right: 20px;
            width: 350px;
            max-height: 500px;
            background: white;
            border-radius: 10px;
            box-shadow: 0 8px 32px rgba(0,0,0,0.4);
            padding: 20px;
            display: none;
            overflow-y: auto;
            z-index: 2000;
        }}
        
        .residue-info-panel.visible {{
            display: block;
        }}
        
        .residue-info-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 15px;
            padding-bottom: 10px;
            border-bottom: 2px solid #667eea;
        }}
        
        .residue-info-title {{
            font-size: 16px;
            font-weight: 700;
            color: #333;
        }}
        
        .close-btn {{
            background: #ff6b6b;
            color: white;
            border: none;
            width: 24px;
            height: 24px;
            border-radius: 50%;
            cursor: pointer;
            font-size: 16px;
            line-height: 1;
            padding: 0;
        }}
        
        .close-btn:hover {{
            background: #ff5252;
        }}
        
        .residue-info-content {{
            font-size: 13px;
        }}
        
        .info-row {{
            display: flex;
            justify-content: space-between;
            padding: 8px 0;
            border-bottom: 1px solid #f0f0f0;
        }}
        
        .info-label {{
            font-weight: 600;
            color: #666;
        }}
        
        .info-value {{
            color: #333;
            text-align: right;
        }}
        
        .info-section {{
            margin-top: 15px;
            padding-top: 10px;
            border-top: 2px solid #f0f0f0;
        }}
        
        .info-section-title {{
            font-weight: 700;
            color: #667eea;
            margin-bottom: 8px;
            font-size: 12px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        
        .loading {{
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            color: white;
            font-size: 14px;
            text-align: center;
            z-index: 10;
        }}
        
        .loading-spinner {{
            border: 3px solid rgba(255,255,255,0.3);
            border-top: 3px solid white;
            border-radius: 50%;
            width: 30px;
            height: 30px;
            animation: spin 1s linear infinite;
            margin: 0 auto 10px;
        }}
        
        @keyframes spin {{
            0% {{ transform: rotate(0deg); }}
            100% {{ transform: rotate(360deg); }}
        }}
        
        .style-buttons {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 8px;
            margin-top: 10px;
        }}
        
        .style-button {{
            padding: 8px;
            font-size: 12px;
            margin: 0;
        }}
        
        .property-info {{
            font-size: 12px;
            color: #666;
            margin-top: 8px;
            padding: 8px;
            background: #f0f4ff;
            border-radius: 6px;
            line-height: 1.5;
        }}
        
        .sync-checkbox {{
            display: flex;
            align-items: center;
            gap: 8px;
            padding: 12px;
            background: #fff3cd;
            border: 2px solid #ffc107;
            border-radius: 8px;
            margin-bottom: 20px;
        }}
        
        .sync-checkbox input[type="checkbox"] {{
            width: 18px;
            height: 18px;
            cursor: pointer;
        }}
        
        .sync-checkbox label {{
            font-size: 13px;
            font-weight: 600;
            color: #856404;
            cursor: pointer;
        }}
    </style>
</head>
<body>
    <div class="container">
        <!--div class="header">
            <h1>Side-by-Side Molecular Visualization with Synchronized Controls</h1>
        </div-->
        
        <div class="main-content">
            <div class="controls-panel">
                <div class="sync-checkbox">
                    <input type="checkbox" id="sync-rotation" checked onchange="toggleSync()">
                    <label for="sync-rotation">🔗 Sync Rotation/Zoom</label>
                </div>
                
                <div class="control-section">
                    <div class="control-label">🎨 Color Property</div>
                    <select id="property-select" onchange="updateAllViewers()">
                        <option value="default">By Variant (distinct colors)</option>
                        <option value="sasa">SASA (Buried Fraction)</option>
                        <option value="interface">Interface Residues</option>
                        <option value="vdw">VdW Energy</option>
                        <option value="contacts">Total Contacts</option>
                        <option value="bfactor">B-factor</option>
                    </select>
                    <div id="property-info" class="property-info"></div>
                </div>
                
                <div class="control-section">
                    <div class="control-label">🔬 Representation</div>
                    <select id="style-select" onchange="updateAllViewers()">
                        <option value="cartoon">Cartoon</option>
                        <option value="stick">Stick</option>
                        <option value="sphere">Sphere</option>
                        <option value="line">Line</option>
                        <option value="surface">Surface</option>
                    </select>
                    <div class="style-buttons">
                        <button class="style-button" onclick="toggleSideChains()">Toggle Side Chains</button>
                        <button class="style-button" onclick="toggleChainOutlines()">Chain Outlines</button>
                    </div>
                </div>
                
                <div class="control-section">
                    <button onclick="updateAllViewers()">🔄 Update All</button>
                    <button onclick="resetAllViews()">🔄 Reset All Views</button>
                    <button onclick="downloadAllImages()">📸 Save All Images</button>
                </div>
                
                <div class="legend" id="legend">
                    <div class="legend-title">Color Scale</div>
                    <div id="legend-content">Select a property to see color scale</div>
                </div>
                
                <div class="info-badge">
                    💡 <strong>Tip:</strong> All viewers share the same styling. Click and drag to rotate, scroll to zoom. Enable "Sync Rotation/Zoom" to move all viewers together. Click any residue for details!
                </div>
            </div>
            
            <div class="viewers-container" id="viewers-container">
"""
    
    # Create viewer panels for each variant
    for variant in variants:
        html += f"""
                <div class="viewer-panel">
                    <div class="viewer-header">{variant}</div>
                    <div class="viewer-viewport" id="viewport-{variant}">
                        <div class="loading" id="loading-{variant}">
                            <div class="loading-spinner"></div>
                            Loading...
                        </div>
                    </div>
                </div>
"""
    
    html += """
            </div>
        </div>
        
        <!-- Residue Information Panel -->
        <div class="residue-info-panel" id="residue-info-panel">
            <div class="residue-info-header">
                <div class="residue-info-title" id="residue-title">Residue Info</div>
                <button class="close-btn" onclick="closeResidueInfo()">×</button>
            </div>
            <div class="residue-info-content" id="residue-info-content">
                Click on a residue to see details
            </div>
        </div>
    </div>
    
    <script>
"""
    
    # Embed PDB data and property data as JavaScript
    html += f"""
        // PDB data
        const pdbData = {json.dumps(pdb_contents, indent=2)};
        
        // Property data for each variant
        const propertyData = {json.dumps(property_data, indent=2)};
        
        // Variant list
        const variants = {json.dumps(variants)};
        
        // Viewer instances
        const viewers = {{}};
        let showSideChains = false;
        let showChainOutlines = false;
        let syncRotation = true;
        
        // Define distinct colors for each variant
        const variantColors = [
            '#FF6B6B', // Red
            '#52B788',  // Green
            '#45B7D1', // Blue
            '#BB8FCE', // Purple
            '#F7DC6F', // Yellow
            '#4ECDC4', // Teal
            '#FFA07A', // Light Salmon
            '#98D8C8', // Mint
            '#85C1E2', // Sky Blue
            '#F8B739' // Orange
        ];
        
        // Map each variant to a color
        const variantColorMap = {{}};
        variants.forEach((variant, idx) => {{
            variantColorMap[variant] = variantColors[idx % variantColors.length];
        }});
        
        // Color scales configuration
        const colorScales = {{
            sasa: {{
                name: 'SASA Buried Fraction',
                min: 0,
                max: 1,
                colors: ['#2166ac', '#67a9cf', '#d1e5f0', '#fddbc7', '#ef8a62', '#b2182b'],
                info: 'Blue = buried (0), Red = exposed (1)',
                getLabel: (v) => v.toFixed(2)
            }},
            interface: {{
                name: 'Interface Residues',
                min: 0,
                max: 1,
                colors: ['#cccccc', '#ff6b6b'],
                info: 'Gray = non-interface, Red = interface',
                getLabel: (v) => v === 1 ? 'Interface' : 'Core'
            }},
            vdw: {{
                name: 'VdW Energy (kcal/mol)',
                min: -15,
                max: 0,
                colors: ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0', '#f7f7f7', '#fddbc7'],
                info: 'Blue = favorable (negative), White = neutral',
                getLabel: (v) => v.toFixed(1)
            }},
            contacts: {{
                name: 'Total Contacts',
                min: 0,
                max: 20,
                colors: ['#ffffcc', '#ffeda0', '#fed976', '#feb24c', '#fd8d3c', '#fc4e2a', '#e31a1c', '#b10026'],
                info: 'Yellow = few contacts, Red = many contacts',
                getLabel: (v) => Math.round(v)
            }},
            bfactor: {{
                name: 'B-factor',
                min: 0,
                max: 100,
                colors: ['#0571b0', '#92c5de', '#f7f7f7', '#f4a582', '#ca0020'],
                info: 'Blue = low B-factor (rigid), Red = high B-factor (flexible)',
                getLabel: (v) => Math.round(v)
            }}
        }};
        
        // Initialize on load
        window.addEventListener('DOMContentLoaded', function() {{
            // Create viewer for each variant
            variants.forEach(variant => {{
                const element = document.getElementById(`viewport-${{variant}}`);
                const config = {{ backgroundColor: 'white' }};
                viewers[variant] = $3Dmol.createViewer(element, config);
                
                // Set up synchronized rotation if enabled
                if (syncRotation) {{
                    setupSyncedRotation(variant);
                }}
            }});
            
            // Initial visualization
            updateAllViewers();
        }});
        
        function setupSyncedRotation(sourceVariant) {{
            const sourceViewer = viewers[sourceVariant];
            
            // Listen for rotation changes
            sourceViewer.spin(false); // Ensure spin is off
            
            // Note: 3Dmol.js doesn't have built-in rotation sync events
            // We'll implement a simple sync on mouse up
            const viewportElement = document.getElementById(`viewport-${{sourceVariant}}`);
            
            viewportElement.addEventListener('mouseup', function() {{
                if (!syncRotation) return;
                
                // Get current view matrix from source
                const view = sourceViewer.getView();
                
                // Apply to all other viewers
                variants.forEach(variant => {{
                    if (variant !== sourceVariant) {{
                        viewers[variant].setView(view);
                        viewers[variant].render();
                    }}
                }});
            }});
            
            // Also sync on zoom (wheel event)
            viewportElement.addEventListener('wheel', function() {{
                if (!syncRotation) return;
                
                setTimeout(() => {{
                    const view = sourceViewer.getView();
                    variants.forEach(variant => {{
                        if (variant !== sourceVariant) {{
                            viewers[variant].setView(view);
                            viewers[variant].render();
                        }}
                    }});
                }}, 50);
            }});
        }}
        
        function toggleSync() {{
            syncRotation = document.getElementById('sync-rotation').checked;
            
            if (syncRotation) {{
                // Re-setup sync for all viewers
                variants.forEach(variant => setupSyncedRotation(variant));
            }}
        }}
        
        function getColor(property, value, variant) {{
            if (property === 'default' || !property) {{
                return null;
            }}
            
            const scale = colorScales[property];
            if (!scale) return null;
            
            const normalized = Math.max(0, Math.min(1, (value - scale.min) / (scale.max - scale.min)));
            
            const colors = scale.colors;
            const idx = normalized * (colors.length - 1);
            const lower = Math.floor(idx);
            
            return colors[lower];
        }}
        
        function applyVisualizationStyle(variant, property, style) {{
            const viewer = viewers[variant];
            if (!viewer) return;
            
            let styleSpec = {{}};
            const variantColor = variantColorMap[variant];
            const props = propertyData[variant];
            
            // Base style
            if (property === 'default') {{
                if (style === 'cartoon') {{
                    styleSpec.cartoon = {{ color: variantColor, clickable: true }};
                }} else if (style === 'stick') {{
                    styleSpec.stick = {{ radius: 0.15, color: variantColor, clickable: true }};
                }} else if (style === 'sphere') {{
                    styleSpec.sphere = {{ radius: 1.0, color: variantColor, clickable: true }};
                }} else if (style === 'line') {{
                    styleSpec.line = {{ color: variantColor, clickable: true }};
                }} else if (style === 'surface') {{
                    styleSpec.surface = {{ opacity: 0.8, color: variantColor, clickable: true }};
                }}
            }} else if (property === 'bfactor') {{
                if (style === 'cartoon') {{
                    styleSpec.cartoon = {{ 
                        colorscheme: {{ prop: 'b', gradient: 'roygb', min: 0, max: 100 }},
                        clickable: true
                    }};
                }} else {{
                    styleSpec[style] = {{ 
                        colorscheme: {{ prop: 'b', gradient: 'roygb', min: 0, max: 100 }},
                        clickable: true
                    }};
                }}
            }} else {{
                if (style === 'cartoon') {{
                    styleSpec.cartoon = {{ style: 'trace', clickable: true }};
                }} else if (style === 'stick') {{
                    styleSpec.stick = {{ radius: 0.15, clickable: true }};
                }} else if (style === 'sphere') {{
                    styleSpec.sphere = {{ radius: 1.0, clickable: true }};
                }} else if (style === 'line') {{
                    styleSpec.line = {{ clickable: true }};
                }} else if (style === 'surface') {{
                    styleSpec.surface = {{ opacity: 0.8, clickable: true }};
                }}
            }}
            
            // Apply base style to all
            viewer.setStyle({{}}, styleSpec);
            
            // Apply property coloring if not default
            if (property !== 'default' && property !== 'bfactor' && props) {{
                Object.keys(props).forEach(key => {{
                    const [chain, resseq] = key.split(':');
                    const value = props[key][property];
                    
                    if (value !== undefined) {{
                        const color = getColor(property, value, variant);
                        if (color) {{
                            const sel = {{ chain: chain, resi: parseInt(resseq) }};
                            const coloredStyle = {{ ...styleSpec }};
                            coloredStyle[style] = {{ ...styleSpec[style], color: color, clickable: true }};
                            viewer.setStyle(sel, coloredStyle);
                        }}
                    }}
                }});
            }}
            
            // Add side chains if requested
            if (showSideChains && property === 'default') {{
                viewer.setStyle({{ atom: ['CA', 'C', 'N', 'O'] }}, {{ stick: {{ radius: 0.2, clickable: true }} }}, {{ byres: false }});
            }}
            
            // Add chain outlines if requested
            if (showChainOutlines) {{
                const modelData = viewer.getModel(0);
                if (modelData) {{
                    const atoms = modelData.selectedAtoms({{}});
                    const chains = [...new Set(atoms.map(a => a.chain))];
                    
                    const chainOutlineColors = [
                        '#8B0000', '#006064', '#1A237E', '#E65100', '#004D40',
                        '#F57F17', '#4A148C', '#01579B', '#BF360C', '#1B5E20'
                    ];
                    
                    chains.forEach((chain, chainIdx) => {{
                        const outlineColor = chainOutlineColors[chainIdx % chainOutlineColors.length];
                        viewer.setStyle(
                            {{ chain: chain }},
                            {{ cartoon: {{ color: outlineColor, thickness: 0.4, opacity: 0.7 }} }},
                            {{ keepStyle: true }}
                        );
                    }});
                }}
            }}
        }}
        
        function updateAllViewers() {{
            const property = document.getElementById('property-select').value;
            const style = document.getElementById('style-select').value;
            
            variants.forEach(variant => {{
                const viewer = viewers[variant];
                if (!viewer) return;
                
                // Show loading
                const loading = document.getElementById(`loading-${{variant}}`);
                if (loading) loading.style.display = 'block';
                
                // Clear and reload
                viewer.clear();
                
                if (!pdbData[variant]) {{
                    console.warn(`No PDB data for ${{variant}}`);
                    if (loading) loading.style.display = 'none';
                    return;
                }}
                
                // Add model
                viewer.addModel(pdbData[variant], 'pdb');
                
                // Apply styling
                applyVisualizationStyle(variant, property, style);
                
                // Set up click handler
                viewer.setClickable({{}}, true, function(atom, viewer, event, container) {{
                    if (atom) {{
                        showResidueInfo(atom, variant);
                        event.stopPropagation();
                    }}
                }});
                
                viewer.zoomTo();
                viewer.render();
                
                // Hide loading
                if (loading) loading.style.display = 'none';
            }});
            
            // Update legend and property info
            updateLegend(property);
            updatePropertyInfo(property);
        }}
        
        function updateLegend(property) {{
            const legendContent = document.getElementById('legend-content');
            
            if (property === 'default') {{
                let html = '<div style="font-size: 11px; font-weight: 600; margin-bottom: 8px;">Variant Colors</div>';
                variants.forEach(variant => {{
                    const color = variantColorMap[variant];
                    html += `
                        <div style="display: flex; align-items: center; margin: 6px 0;">
                            <div style="width: 24px; height: 16px; background: ${{color}}; border-radius: 3px; margin-right: 8px; border: 1px solid #ccc;"></div>
                            <span style="font-size: 12px;">${{variant}}</span>
                        </div>
                    `;
                }});
                legendContent.innerHTML = html;
                return;
            }}
            
            const scale = colorScales[property];
            if (!scale) return;
            
            const gradient = scale.colors.join(', ');
            
            legendContent.innerHTML = `
                <div style="font-size: 11px; font-weight: 600; margin-bottom: 8px;">${{scale.name}}</div>
                <div class="legend-gradient" style="background: linear-gradient(to right, ${{gradient}});"></div>
                <div class="legend-labels">
                    <span>${{scale.getLabel(scale.min)}}</span>
                    <span>${{scale.getLabel(scale.max)}}</span>
                </div>
                <div style="font-size: 11px; color: #666; margin-top: 8px;">${{scale.info}}</div>
            `;
        }}
        
        function updatePropertyInfo(property) {{
            const infoDiv = document.getElementById('property-info');
            
            const descriptions = {{
                'default': 'Each variant is colored with a distinct color for easy identification.',
                'sasa': 'Solvent Accessible Surface Area - measures how buried each residue is.',
                'interface': 'Highlights residues at protein-protein interfaces.',
                'vdw': 'Van der Waals energy - negative values indicate favorable interactions.',
                'contacts': 'Number of atomic contacts each residue makes.',
                'bfactor': 'Temperature/B-factor - indicates structural flexibility or disorder.'
            }};
            
            infoDiv.textContent = descriptions[property] || '';
        }}
        
        function toggleSideChains() {{
            showSideChains = !showSideChains;
            updateAllViewers();
        }}
        
        function toggleChainOutlines() {{
            showChainOutlines = !showChainOutlines;
            updateAllViewers();
        }}
        
        function resetAllViews() {{
            variants.forEach(variant => {{
                const viewer = viewers[variant];
                if (viewer) {{
                    viewer.zoomTo();
                    viewer.render();
                }}
            }});
        }}
        
        function downloadAllImages() {{
            variants.forEach(variant => {{
                const viewer = viewers[variant];
                if (viewer) {{
                    viewer.render();
                    const imgData = viewer.pngURI();
                    const link = document.createElement('a');
                    link.href = imgData;
                    link.download = `${{variant}}_structure.png`;
                    link.click();
                }}
            }});
        }}
        
        function showResidueInfo(atom, variant) {{
            const panel = document.getElementById('residue-info-panel');
            const title = document.getElementById('residue-title');
            const content = document.getElementById('residue-info-content');
            
            const resname = atom.resn || 'Unknown';
            const chain = atom.chain || '?';
            const resi = atom.resi || '?';
            const atomName = atom.atom || '?';
            const elem = atom.elem || '?';
            const bfactor = atom.b !== undefined ? atom.b.toFixed(2) : 'N/A';
            
            const resKey = `${{chain}}:${{resi}}`;
            const variantProps = propertyData[variant];
            const resProps = variantProps ? variantProps[resKey] : null;
            
            title.textContent = `${{resname}} ${{chain}}:${{resi}} (${{variant}})`;
            
            let html = `
                <div class="info-row">
                    <span class="info-label">Variant:</span>
                    <span class="info-value" style="color: ${{variantColorMap[variant]}}; font-weight: 700;">${{variant}}</span>
                </div>
                <div class="info-row">
                    <span class="info-label">Residue Type:</span>
                    <span class="info-value">${{resname}}</span>
                </div>
                <div class="info-row">
                    <span class="info-label">Chain:</span>
                    <span class="info-value">${{chain}}</span>
                </div>
                <div class="info-row">
                    <span class="info-label">Position:</span>
                    <span class="info-value">${{resi}}</span>
                </div>
                <div class="info-row">
                    <span class="info-label">Atom Clicked:</span>
                    <span class="info-value">${{atomName}} (${{elem}})</span>
                </div>
                <div class="info-row">
                    <span class="info-label">B-factor:</span>
                    <span class="info-value">${{bfactor}}</span>
                </div>
            `;
            
            if (resProps) {{
                html += `
                    <div class="info-section">
                        <div class="info-section-title">Analysis Properties</div>
                        
                        <div class="info-row">
                            <span class="info-label">SASA (Buried):</span>
                            <span class="info-value">${{resProps.sasa.toFixed(2)}}</span>
                        </div>
                        
                        <div class="info-row">
                            <span class="info-label">Interface:</span>
                            <span class="info-value">${{resProps.interface === 1 ? '✓ Yes' : '✗ No'}}</span>
                        </div>
                        
                        <div class="info-row">
                            <span class="info-label">VdW Energy:</span>
                            <span class="info-value">${{resProps.vdw.toFixed(2)}} kcal/mol</span>
                        </div>
                        
                        <div class="info-row">
                            <span class="info-label">Total Contacts:</span>
                            <span class="info-value">${{resProps.contacts}}</span>
                        </div>
                    </div>
                `;
            }} else {{
                html += `
                    <div class="info-section">
                        <div class="info-section-title">Analysis Properties</div>
                        <div style="color: #999; font-style: italic; font-size: 12px;">
                            No analysis data available
                        </div>
                    </div>
                `;
            }}
            
            content.innerHTML = html;
            panel.classList.add('visible');
        }}
        
        function closeResidueInfo() {{
            const panel = document.getElementById('residue-info-panel');
            panel.classList.remove('visible');
        }}
    </script>
</body>
</html>
"""
    
    # Save to file if requested
    if output_file:
        with open(output_file, 'w') as f:
            f.write(html)
        print(f"Saved side-by-side 3D visualization to {output_file}")
    
    return html

