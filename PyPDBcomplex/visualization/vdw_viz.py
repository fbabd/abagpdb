from __future__ import annotations
from typing import Dict, List, Tuple, Optional
import json
import os

from ..models import Complex, Residue
from ..vdw import (
    per_residue_LJ_decomposition,
    per_residue_pair_LJ,
    identify_energetic_hotspots,
    compare_energetic_contributions,
)


# ---------------------------
# Color Mapping for Energy Visualization
# ---------------------------
def energy_to_rgb(energy: float, vmin: float = -5.0, vmax: float = 1.0) -> Tuple[int, int, int]:
    """
    Map energy value to RGB tuple.
    
    Blue (favorable) -> White (neutral) -> Red (unfavorable)
    
    Args:
        energy: VDW energy value in kcal/mol
        vmin: Minimum energy for color scale
        vmax: Maximum energy for color scale
    
    Returns:
        RGB tuple (r, g, b) with values 0-255
    """
    # Normalize energy to [0, 1]
    if energy <= vmin:
        norm = 0.0
    elif energy >= vmax:
        norm = 1.0
    else:
        norm = (energy - vmin) / (vmax - vmin)
    
    # Blue -> White -> Red gradient
    if norm < 0.5:
        # Blue to White
        t = norm * 2.0
        r = int(0 + t * 255)
        g = int(0 + t * 255)
        b = 255
    else:
        # White to Red
        t = (norm - 0.5) * 2.0
        r = 255
        g = int(255 - t * 255)
        b = int(255 - t * 255)
    
    return (r, g, b)


def energy_to_hex(energy: float, vmin: float = -5.0, vmax: float = 1.0) -> str:
    """Convert energy to hex color string."""
    r, g, b = energy_to_rgb(energy, vmin, vmax)
    return f"#{r:02x}{g:02x}{b:02x}"


# ---------------------------
# 3D Molecular Visualization with 3Dmol.js
# ---------------------------
def generate_3dmol_viewer(
    pdb_content: str,
    per_res_energy: Dict[str, float],
    viewer_id: str = "vdw_viewer",
    width: int = 800,
    height: int = 600,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
) -> str:
    """
    Generate HTML div with 3Dmol.js viewer showing VDW energies.
    
    Args:
        pdb_content: PDB file content as string
        per_res_energy: Dictionary of residue_id -> energy
        viewer_id: HTML element ID for the viewer
        width: Viewer width in pixels
        height: Viewer height in pixels
        vmin: Minimum energy for color scale
        vmax: Maximum energy for color scale
    
    Returns:
        HTML string with embedded 3Dmol.js viewer
    """
    if not per_res_energy:
        return "<p>No energy data to visualize</p>"
    
    # Auto-detect energy range
    energies = list(per_res_energy.values())
    if vmin is None:
        vmin = min(energies)
    if vmax is None:
        vmax = max(energies)
    
    # Build residue color mapping
    residue_colors = {}
    for res_id, energy in per_res_energy.items():
        hex_color = energy_to_hex(energy, vmin, vmax)
        # Parse residue ID (format: "ASP A_123")
        parts = res_id.split()
        if len(parts) >= 2:
            resname = parts[0]  # Extract residue name (e.g., "ASP")
            chain_res = parts[1]
            if '_' in chain_res:
                chain, resnum = chain_res.split('_')
            elif ':' in chain_res:
                chain, resnum = chain_res.split(':')
            else:
                continue
            
            key = f"{chain}:{resnum}"
            residue_colors[key] = {
                'color': hex_color, 
                'energy': energy,
                'resname': resname  # Store residue name
            }
    
    # Escape PDB content for JavaScript
    pdb_escaped = json.dumps(pdb_content)
    colors_json = json.dumps(residue_colors)
    
    html = f"""
<div id="{viewer_id}" style="width: {width}px; height: {height}px; position: relative; margin: 20px auto;"></div>
<div id="{viewer_id}_legend" style="width: {width}px; margin: 10px auto; text-align: center;">
    <div style="display: inline-block; background: linear-gradient(to right, #0000ff, #ffffff, #ff0000); width: 300px; height: 20px; border: 1px solid #ccc;"></div>
    <div style="margin-top: 5px;">
        <span style="float: left; margin-left: calc(50% - 150px);">{vmin:.2f} kcal/mol (Favorable)</span>
        <span style="float: right; margin-right: calc(50% - 150px);">{vmax:.2f} kcal/mol (Unfavorable)</span>
    </div>
    <div style="clear: both; margin-top: 10px; font-size: 0.9em; color: #666;">
        <strong>Click</strong> on residues to show/hide details
    </div>
</div>

<script>
// Wait for both jQuery and 3Dmol to be loaded
(function initViewer() {{
    if (typeof $ === 'undefined' || typeof $3Dmol === 'undefined') {{
        console.log('Waiting for libraries to load...');
        setTimeout(initViewer, 100);
        return;
    }}
    
    console.log('Initializing 3D viewer: {viewer_id}');
    
    $(document).ready(function() {{
        try {{
            var viewer = $3Dmol.createViewer("{viewer_id}", {{
                backgroundColor: 'white'
            }});
            
            if (!viewer) {{
                console.error('Failed to create viewer');
                return;
            }}
            
            var pdbData = {pdb_escaped};
            var residueColors = {colors_json};
            
            console.log('Loading PDB data...');
            viewer.addModel(pdbData, "pdb");
            
            // Set cartoon style for entire structure
            viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray'}}}});
            
            // Color interface residues by energy
            console.log('Coloring', Object.keys(residueColors).length, 'residues');
            for (var resKey in residueColors) {{
                var parts = resKey.split(':');
                var chain = parts[0];
                var resi = parts[1];
                var color = residueColors[resKey].color;
                
                // Color cartoon and add sticks for interface residues
                // IMPORTANT: Add clickable: true to make these atoms clickable!
                viewer.setStyle(
                    {{chain: chain, resi: resi}},
                    {{
                        cartoon: {{color: color}},
                        stick: {{color: color, radius: 0.3, clickable: true}}
                    }}
                );
            }}
            
            viewer.zoomTo();
            viewer.render();
            
            console.log('Viewer rendered successfully');
            
            // Track persistent click labels
            var persistentLabels = [];
            
            // Add click handler to show residue details
            console.log('Setting up click handler');
            
            // IMPORTANT: Use simplified signature - just (atom, viewer)
            viewer.setClickable({{}}, true, function(atom, viewer) {{
                console.log('>>> CLICK EVENT <<<');
                console.log('Atom:', atom);
                
                if (atom) {{
                    var resKey = atom.chain + ':' + atom.resi;
                    console.log('Clicked residue key:', resKey);
                    var info = residueColors[resKey];
                    
                    if (info) {{
                        console.log('Found residue info:', info);
                        
                        // Check if this residue already has a label
                        var existingIndex = -1;
                        for (var i = 0; i < persistentLabels.length; i++) {{
                            if (persistentLabels[i].resKey === resKey) {{
                                existingIndex = i;
                                break;
                            }}
                        }}
                        
                        if (existingIndex >= 0) {{
                            // Remove existing label for this residue (toggle off)
                            console.log('Removing existing label');
                            viewer.removeLabel(persistentLabels[existingIndex].label);
                            persistentLabels.splice(existingIndex, 1);
                        }} else {{
                            // Add new label
                            console.log('Adding new label');
                            var labelText = info.resname + ' ' + atom.chain + ':' + atom.resi + 
                                       '\\nEnergy: ' + info.energy.toFixed(2) + ' kcal/mol' +
                                       '\\n(Click to toggle)';
                            
                            var label = viewer.addLabel(labelText, {{
                                position: atom,
                                backgroundColor: 'white',
                                fontColor: 'black',
                                fontSize: 14,
                                borderThickness: 2,
                                borderColor: '#333333',
                                inFront: true
                            }});
                            
                            persistentLabels.push({{
                                resKey: resKey,
                                label: label
                            }});
                            console.log('Label added! Total labels:', persistentLabels.length);
                        }}
                        
                        viewer.render();
                    }} else {{
                        console.log('No energy info for residue:', resKey);
                        console.log('Available residues:', Object.keys(residueColors));
                    }}
                }} else {{
                    // Click on empty space removes all labels
                    console.log('Clicked empty space, clearing', persistentLabels.length, 'labels');
                    for (var i = 0; i < persistentLabels.length; i++) {{
                        viewer.removeLabel(persistentLabels[i].label);
                    }}
                    persistentLabels = [];
                    viewer.render();
                }}
            }});
            
            console.log('Click handler set up successfully');
            console.log('Try clicking on the colored residues!');
            
        }} catch (error) {{
            console.error('Error initializing viewer:', error);
        }}
    }});
}})();
</script>
"""
    
    return html


# ---------------------------
# Interactive Plotly Charts
# ---------------------------
def create_plotly_distribution(
    per_res_energy: Dict[str, float],
    title: str = "VDW Energy Distribution",
) -> str:
    """
    Create interactive Plotly histogram of energy distribution.
    
    Args:
        per_res_energy: Dictionary of residue_id -> energy
        title: Chart title
    
    Returns:
        HTML string with embedded Plotly chart
    """
    energies = list(per_res_energy.values())
    
    data = [{
        'x': energies,
        'type': 'histogram',
        'nbinsx': 50,
        'marker': {
            'color': energies,
            'colorscale': [[0, 'blue'], [0.5, 'white'], [1, 'red']],
            'colorbar': {
                'title': 'Energy<br>(kcal/mol)',
                'thickness': 15,
                'len': 0.7
            },
            'line': {'color': 'black', 'width': 1}
        },
        'name': 'Residues'
    }]
    
    layout = {
        'title': {
            'text': title,
            'font': {'size': 18, 'family': 'Arial, sans-serif'}
        },
        'xaxis': {
            'title': 'VDW Energy (kcal/mol)',
            'gridcolor': '#e0e0e0'
        },
        'yaxis': {
            'title': 'Number of Residues',
            'gridcolor': '#e0e0e0'
        },
        'plot_bgcolor': 'white',
        'hovermode': 'closest',
        'showlegend': False
    }
    
    # Add statistics annotation
    import statistics
    mean_e = statistics.mean(energies)
    median_e = statistics.median(energies)
    std_e = statistics.stdev(energies) if len(energies) > 1 else 0
    
    annotation = {
        'text': f'Total: {len(energies)} residues<br>' +
                f'Mean: {mean_e:.2f} kcal/mol<br>' +
                f'Median: {median_e:.2f} kcal/mol<br>' +
                f'Std: {std_e:.2f} kcal/mol',
        'xref': 'paper',
        'yref': 'paper',
        'x': 0.02,
        'y': 0.98,
        'xanchor': 'left',
        'yanchor': 'top',
        'showarrow': False,
        'bgcolor': 'wheat',
        'bordercolor': 'black',
        'borderwidth': 1,
        'font': {'size': 11}
    }
    layout['annotations'] = [annotation]
    
    config = {
        'responsive': True,
        'displayModeBar': True,
        'displaylogo': False
    }
    
    plot_json = {
        'data': data,
        'layout': layout,
        'config': config
    }
    
    div_id = "plotly_distribution_" + str(hash(title))[:8]
    
    html = f"""
<div id="{div_id}" style="width: 100%; height: 500px;"></div>
<script>
    Plotly.newPlot('{div_id}', {json.dumps(plot_json['data'])}, 
                   {json.dumps(plot_json['layout'])}, 
                   {json.dumps(plot_json['config'])});
</script>
"""
    
    return html


def create_plotly_chain_comparison(
    per_res_energy: Dict[str, float],
    title: str = "VDW Energy by Chain",
) -> str:
    """
    Create interactive Plotly bar chart comparing chains.
    
    Args:
        per_res_energy: Dictionary of residue_id -> energy
        title: Chart title
    
    Returns:
        HTML string with embedded Plotly chart
    """
    summary = compare_energetic_contributions(per_res_energy)
    
    if not summary:
        return "<p>No chain data available</p>"
    
    chains = sorted(summary.keys())
    
    # Prepare data for grouped bar chart
    data = [
        {
            'x': chains,
            'y': [summary[c]['total_energy'] for c in chains],
            'name': 'Total Energy',
            'type': 'bar',
            'marker': {'color': 'steelblue'}
        },
        {
            'x': chains,
            'y': [summary[c]['favorable_energy'] for c in chains],
            'name': 'Favorable (< 0)',
            'type': 'bar',
            'marker': {'color': 'green'}
        },
        {
            'x': chains,
            'y': [summary[c]['unfavorable_energy'] for c in chains],
            'name': 'Unfavorable (≥ 0)',
            'type': 'bar',
            'marker': {'color': 'red'}
        }
    ]
    
    layout = {
        'title': {
            'text': title,
            'font': {'size': 18}
        },
        'xaxis': {'title': 'Chain'},
        'yaxis': {'title': 'VDW Energy (kcal/mol)'},
        'barmode': 'group',
        'plot_bgcolor': 'white',
        'hovermode': 'closest',
        'showlegend': True,
        'legend': {'x': 1.02, 'y': 1}
    }
    
    config = {'responsive': True, 'displaylogo': False}
    
    div_id = "plotly_chains_" + str(hash(title))[:8]
    
    html = f"""
<div id="{div_id}" style="width: 100%; height: 500px;"></div>
<script>
    Plotly.newPlot('{div_id}', {json.dumps(data)}, 
                   {json.dumps(layout)}, 
                   {json.dumps(config)});
</script>
"""
    
    return html


def create_plotly_hotspot_ranking(
    per_res_energy: Dict[str, float],
    energy_threshold: float = -2.0,
    top_n: int = 20,
    title: str = "Top Energetic Hotspots",
) -> str:
    """
    Create interactive Plotly horizontal bar chart of hotspots.
    
    Args:
        per_res_energy: Dictionary of residue_id -> energy
        top_n: Number of top hotspots to display
        title: Chart title
    
    Returns:
        HTML string with embedded Plotly chart
    """
    hotspots = identify_energetic_hotspots(per_res_energy, energy_threshold, top_n=top_n)
    
    if not hotspots:
        return "<p>No hotspots found</p>"
    
    res_ids, energies = zip(*hotspots)
    
    # Shorten residue IDs for display
    display_ids = [rid.split()[1] if ' ' in rid else rid for rid in res_ids]
    
    data = [{
        'y': display_ids[::-1],  # Reverse to show strongest at top
        'x': list(energies)[::-1],
        'type': 'bar',
        'orientation': 'h',
        'marker': {
            'color': list(energies)[::-1],
            'colorscale': [[0, 'darkblue'], [0.5, 'blue'], [1, 'lightblue']],
            'colorbar': {
                'title': 'Energy<br>(kcal/mol)',
                'thickness': 15,
                'len': 0.7
            },
            'line': {'color': 'black', 'width': 1}
        },
        'text': [f'{e:.2f}' for e in energies][::-1],
        'textposition': 'inside',
        'textfont': {'color': 'white', 'size': 10},
        'hovertemplate': '<b>%{y}</b><br>Energy: %{x:.2f} kcal/mol<extra></extra>'
    }]
    
    layout = {
        'title': {
            'text': title,
            'font': {'size': 18}
        },
        'xaxis': {
            'title': 'VDW Energy (kcal/mol)',
            'gridcolor': '#e0e0e0'
        },
        'yaxis': {
            'title': '',
            'tickfont': {'size': 10}
        },
        'plot_bgcolor': 'white',
        'height': max(400, len(hotspots) * 50),
        'margin': {'l': 100, 'r': 100, 't': 80, 'b': 60},
        'showlegend': False
    }
    
    config = {'responsive': True, 'displaylogo': False}
    
    div_id = "plotly_hotspots_" + str(hash(title))[:8]
    
    html = f"""
<div id="{div_id}" style="width: 100%; height: {layout['height']}px;"></div>
<script>
    Plotly.newPlot('{div_id}', {json.dumps(data)}, 
                   {json.dumps(layout)}, 
                   {json.dumps(config)});
</script>
"""
    
    return html


def create_plotly_residue_profile(
    per_res_energy: Dict[str, float],
    chain_id: str,
    title: Optional[str] = None,
    highlight_threshold: float = -2.0,
) -> str:
    """
    Create interactive Plotly line chart of energy along residue sequence.
    
    Args:
        per_res_energy: Dictionary of residue_id -> energy
        chain_id: Chain ID to plot
        title: Chart title
        highlight_threshold: Threshold for highlighting hotspots
    
    Returns:
        HTML string with embedded Plotly chart
    """
    # Filter residues for this chain
    chain_data = []
    for res_id, energy in per_res_energy.items():
        parts = res_id.split()
        if len(parts) >= 2:
            chain_res = parts[1]
            if '_' in chain_res:
                chain, resnum = chain_res.split('_')
            elif ':' in chain_res:
                chain, resnum = chain_res.split(':')
            else:
                continue
            
            if chain == chain_id:
                try:
                    resnum_int = int(resnum)
                    chain_data.append((resnum_int, energy, res_id))
                except ValueError:
                    continue
    
    if not chain_data:
        return f"<p>No data found for chain {chain_id}</p>"
    
    # Sort by residue number
    chain_data.sort(key=lambda x: x[0])
    resnums, energies, res_ids = zip(*chain_data)
    
    # Main energy trace
    data = [
        {
            'x': list(resnums),
            'y': list(energies),
            'type': 'scatter',
            'mode': 'lines+markers',
            'name': 'VDW Energy',
            'line': {'color': 'steelblue', 'width': 2},
            'marker': {'size': 6, 'color': 'steelblue'},
            'hovertemplate': '<b>Residue %{x}</b><br>Energy: %{y:.2f} kcal/mol<extra></extra>'
        }
    ]
    
    # Add hotspot markers
    hotspot_x = [r for r, e in zip(resnums, energies) if e <= highlight_threshold]
    hotspot_y = [e for e in energies if e <= highlight_threshold]
    
    if hotspot_x:
        data.append({
            'x': hotspot_x,
            'y': hotspot_y,
            'type': 'scatter',
            'mode': 'markers',
            'name': f'Hotspots (E ≤ {highlight_threshold})',
            'marker': {
                'size': 12,
                'color': 'red',
                'symbol': 'star',
                'line': {'color': 'darkred', 'width': 1}
            },
            'hovertemplate': '<b>Hotspot: Res %{x}</b><br>Energy: %{y:.2f} kcal/mol<extra></extra>'
        })
    
    if title is None:
        title = f"VDW Energy Profile - Chain {chain_id}"
    
    layout = {
        'title': {
            'text': title,
            'font': {'size': 18}
        },
        'xaxis': {
            'title': 'Residue Number',
            'gridcolor': '#e0e0e0'
        },
        'yaxis': {
            'title': 'VDW Energy (kcal/mol)',
            'gridcolor': '#e0e0e0',
            'zeroline': True,
            'zerolinecolor': 'black',
            'zerolinewidth': 2
        },
        'plot_bgcolor': 'white',
        'hovermode': 'closest',
        'showlegend': True,
        'shapes': [
            {
                'type': 'line',
                'x0': min(resnums),
                'x1': max(resnums),
                'y0': highlight_threshold,
                'y1': highlight_threshold,
                'line': {
                    'color': 'red',
                    'dash': 'dash',
                    'width': 2
                }
            }
        ]
    }
    
    # Add statistics annotation
    import statistics
    mean_e = statistics.mean(energies)
    
    annotation = {
        'text': f'Chain {chain_id}<br>' +
                f'Residues: {len(energies)}<br>' +
                f'Mean: {mean_e:.2f}<br>' +
                f'Hotspots: {len(hotspot_x)}',
        'xref': 'paper',
        'yref': 'paper',
        'x': 0.02,
        'y': 0.98,
        'xanchor': 'left',
        'yanchor': 'top',
        'showarrow': False,
        'bgcolor': 'wheat',
        'bordercolor': 'black',
        'borderwidth': 1
    }
    layout['annotations'] = [annotation]
    
    config = {'responsive': True, 'displaylogo': False}
    
    div_id = f"plotly_profile_{chain_id}_" + str(hash(title or ""))[:8]
    
    html = f"""
<div id="{div_id}" style="width: 100%; height: 500px;"></div>
<script>
    Plotly.newPlot('{div_id}', {json.dumps(data)}, 
                   {json.dumps(layout)}, 
                   {json.dumps(config)});
</script>
"""
    
    return html


def create_plotly_contact_matrix(
    cx: Complex,
    chain_exprs: List[str],
    cutoff: float = 5.0,
    title: str = "VDW Contact Matrix",
) -> str:
    """
    Create interactive Plotly heatmap of chain-chain VDW energies.
    
    Args:
        cx: Complex object
        chain_exprs: List of chain selection expressions
        cutoff: Distance cutoff for contacts
        title: Chart title
    
    Returns:
        HTML string with embedded Plotly heatmap
    """
    # Calculate pairwise energies
    pair_energies = per_residue_pair_LJ(cx, chain_exprs, cutoff=cutoff)
    
    # Build matrix
    n = len(chain_exprs)
    matrix = [[0.0 for _ in range(n)] for _ in range(n)]
    
    for i, chain1 in enumerate(chain_exprs):
        for j, chain2 in enumerate(chain_exprs):
            if i == j:
                continue
            key = (chain1, chain2)
            if key in pair_energies:
                total_energy = sum(pair_energies[key].values())
                matrix[i][j] = total_energy
    
    # Create annotations for cell values
    annotations = []
    for i in range(n):
        for j in range(n):
            if i != j:
                annotations.append({
                    'x': j,
                    'y': i,
                    'text': f'{matrix[i][j]:.1f}',
                    'showarrow': False,
                    'font': {'size': 12, 'color': 'black'}
                })
    
    data = [{
        'z': matrix,
        'x': chain_exprs,
        'y': chain_exprs,
        'type': 'heatmap',
        'colorscale': [[0, 'blue'], [0.5, 'white'], [1, 'red']],
        'colorbar': {
            'title': 'Total VDW<br>Energy<br>(kcal/mol)',
            'thickness': 20,
            'len': 0.7
        },
        'hovertemplate': '<b>%{y} → %{x}</b><br>Energy: %{z:.2f} kcal/mol<extra></extra>'
    }]
    
    layout = {
        'title': {
            'text': title,
            'font': {'size': 18}
        },
        'xaxis': {
            'title': 'Chain',
            'side': 'bottom'
        },
        'yaxis': {
            'title': 'Chain',
            'autorange': 'reversed'
        },
        'annotations': annotations,
        'width': 600,
        'height': 600
    }
    
    config = {'responsive': True, 'displaylogo': False}
    
    div_id = "plotly_matrix_" + str(hash(title))[:8]
    
    html = f"""
<div id="{div_id}" style="width: 600px; height: 600px; margin: 20px auto;"></div>
<script>
    Plotly.newPlot('{div_id}', {json.dumps(data)}, 
                   {json.dumps(layout)}, 
                   {json.dumps(config)});
</script>
"""
    
    return html


# ---------------------------
# Complete HTML Report Generator
# ---------------------------
def generate_html_report(
    cx: Complex,
    left_exprs: List[str],
    right_exprs: List[str],
    pdb_content: str,
    output_file: str = "vdw_report.html",
    cutoff: float = 5.0,
    energy_threshold: float = -2.0,
    top_n: int = 20,
    include_3d: bool = True,
) -> str:
    """
    Generate complete HTML report with interactive visualizations.
    
    Args:
        cx: Complex object
        left_exprs: Selection expression for left side
        right_exprs: Selection expression for right side
        pdb_content: PDB file content as string
        output_file: Output HTML filename
        cutoff: Distance cutoff (Angstroms)
        energy_threshold: Hotspot threshold (kcal/mol)
        top_n: Number of top hotspots
        include_3d: Include 3D molecular viewer
    
    Returns:
        Path to generated HTML file
    """
    print("=" * 80)
    print("GENERATING HTML VDW VISUALIZATION REPORT")
    print("=" * 80)
    
    # Calculate energies 
    per_res_energy = per_residue_LJ_decomposition(
        cx, left_exprs, right_exprs, cutoff=cutoff
    )
    
    if not per_res_energy:
        print("Error: No energies calculated")
        return ""
    
    # print(f"   Calculated energies for {len(per_res_energy)} residues")
    
    # Generate visualizations
    
    viewer_html = ""
    if include_3d:
        viewer_html = generate_3dmol_viewer(pdb_content, per_res_energy)
    
    distribution_html = create_plotly_distribution(per_res_energy)
    chain_comparison_html = create_plotly_chain_comparison(per_res_energy)
    hotspot_ranking_html = create_plotly_hotspot_ranking(per_res_energy, energy_threshold, top_n=top_n)
    
    # Generate per-chain profiles
    summary = compare_energetic_contributions(per_res_energy)
    chain_profiles_html = ""
    for chain_id in sorted(summary.keys()):
        profile = create_plotly_residue_profile(
            per_res_energy, 
            chain_id, 
            highlight_threshold=energy_threshold
        )
        chain_profiles_html += f"""
        <div class="chart-container">
            <h3>Chain {chain_id} Energy Profile</h3>
            {profile}
        </div>
        """
    
    # Calculate summary statistics
    all_energies = list(per_res_energy.values())
    total_energy = sum(all_energies)
    favorable = sum(e for e in all_energies if e < 0)
    unfavorable = sum(e for e in all_energies if e >= 0)
    hotspots = identify_energetic_hotspots(per_res_energy, energy_threshold=energy_threshold)
    
    import statistics
    mean_energy = statistics.mean(all_energies)
    median_energy = statistics.median(all_energies)
    
    # Build HTML document
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>VDW Energy Analysis Report</title>
    
    <!-- External Libraries -->
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: white;
            padding: 20px;
            color: #333;
        }}
        
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            border-radius: 10px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
            overflow: hidden;
        }}
        
        header {{
            background: linear-gradient(135deg, #fffffff 0%, #ffffff 100%);
            color: white;
            padding: 40px;
            text-align: center;
        }}
        
        header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
            color:black;
        }}
        
        header p {{
            font-size: 1.1em;
            opacity: 0.9;
        }}
        
        .summary {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            padding: 30px;
            background: #f8f9fa;
            border-bottom: 2px solid #e0e0e0;
        }}
        
        .stat-card {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            text-align: center;
            transition: transform 0.3s ease;
        }}
        
        .stat-card:hover {{
            transform: translateY(-5px);
        }}
        
        .stat-card h3 {{
            color: #667eea;
            font-size: 0.9em;
            text-transform: uppercase;
            letter-spacing: 1px;
            margin-bottom: 10px;
        }}
        
        .stat-card .value {{
            font-size: 2em;
            font-weight: bold;
            color: #333;
        }}
        
        .stat-card .unit {{
            font-size: 0.8em;
            color: #666;
        }}
        
        .content {{
            padding: 40px;
        }}
        
        .section {{
            margin-bottom: 50px;
        }}
        
        .section h2 {{
            color: #667eea;
            font-size: 2em;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 3px solid #667eea;
        }}
        
        .chart-container {{
            margin: 30px 0;
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.05);
        }}
        
        .chart-container h3 {{
            color: #333;
            margin-bottom: 15px;
            font-size: 1.3em;
        }}
        
        .viewer-container {{
            margin: 30px 0;
            text-align: center;
        }}
        
        footer {{
            background: #2c3e50;
            color: white;
            text-align: center;
            padding: 20px;
            font-size: 0.9em;
        }}
        
        .info-box {{
            background: #e3f2fd;
            border-left: 4px solid #2196f3;
            padding: 15px;
            margin: 20px 0;
            border-radius: 4px;
        }}
        
        .info-box strong {{
            color: #1976d2;
        }}
        
        @media print {{
            body {{
                background: white;
            }}
            .container {{
                box-shadow: none;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>Comprehensive van der Waals interaction analysis</h1>
        </header>
        
        <div class="summary">
            <div class="stat-card">
                <h3>Total Residues</h3>
                <div class="value">{len(all_energies)}</div>
            </div>
            <div class="stat-card">
                <h3>Total Energy</h3>
                <div class="value">{total_energy:.2f}</div>
                <div class="unit">kcal/mol</div>
            </div>
            <div class="stat-card">
                <h3>Mean Energy</h3>
                <div class="value">{mean_energy:.2f}</div>
                <div class="unit">kcal/mol</div>
            </div>
            <!--div class="stat-card">
                <h3>Median Energy</h3>
                <div class="value">{median_energy:.2f}</div>
                <div class="unit">kcal/mol</div>
            </div-->
            <div class="stat-card">
                <h3>Hotspots</h3>
                <div class="value">{len(hotspots)}</div>
                <div class="unit">E &lt; {energy_threshold:.1f} kcal/mol</div>
            </div>
        </div>
        
        <div class="content">
            {"" if not include_3d else f'''
            <div class="section">
                <h2>3D Molecular Viewer</h2>
                <div class="info-box">
                    <strong>Interactive Controls:</strong> 
                    <ul style="margin: 10px 0; padding-left: 20px;">
                        <li>Click and drag to rotate the structure</li>
                        <li>Scroll to zoom in/out</li>
                        <li><strong>Click</strong> on residues to show energy information</li>
                        <li>Click multiple residues to compare side-by-side</li>
                        <li>Click on a residue again to hide its label</li>
                        <li>Click empty space to clear all labels</li>
                    </ul>
                </div>
                <div class="viewer-container">
                    {viewer_html}
                </div>
            </div>
            '''}
            
            <div class="section">
                <h2>Energy Distribution</h2>
                <div class="chart-container">
                    {distribution_html}
                </div>
            </div>
            
            <div class="section">
                <h2>Top Energetic Hotspots</h2>
                <div class="info-box">
                    Residues with strongest favorable VDW interactions (most negative energies)
                </div>
                <div class="chart-container">
                    {hotspot_ranking_html}
                </div>
            </div>
            
            <div class="section">
                <h2>Per-Chain Energy Profiles</h2>
                {chain_profiles_html}
            </div>
            
            <!--div class="section">
                <h2>Chain Comparison</h2>
                <div class="chart-container">
                    {chain_comparison_html}
                </div>
            </div-->
            
        </div>
        
        <footer>
            <p>Generated VDW Energy Analysis Report | Cutoff: {cutoff:.1f} Å | Threshold: {energy_threshold:.1f} kcal/mol</p>
            <p>Blue = Favorable interactions | Red = Unfavorable interactions</p>
        </footer>
    </div>
</body>
</html>
"""
    
    # Write to file
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"\n3. HTML report generated: {output_file}")
    print("=" * 80)
    
    return os.path.abspath(output_file)


# ---------------------------
# Utility function to generate standalone HTML components
# ---------------------------
def generate_html_component(
    per_res_energy: Dict[str, float],
    component_type: str = "distribution",
    **kwargs
) -> str:
    """
    Generate a standalone HTML component for embedding.
    
    Args:
        per_res_energy: Dictionary of residue_id -> energy
        component_type: Type of component ("distribution", "chains", 
                       "hotspots", "profile", "matrix")
        **kwargs: Additional arguments passed to the component function
    
    Returns:
        HTML string with the component (includes necessary scripts)
    
    Example:
        >>> energies = per_residue_LJ_decomposition(cx, ["A"], ["B"])
        >>> html = generate_html_component(energies, "distribution")
        >>> # Embed in your web page
    """
    # Common header with required libraries
    header = """
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
"""
    
    if component_type == "distribution":
        content = create_plotly_distribution(per_res_energy, **kwargs)
    elif component_type == "chains":
        content = create_plotly_chain_comparison(per_res_energy, **kwargs)
    elif component_type == "hotspots":
        content = create_plotly_hotspot_ranking(per_res_energy, **kwargs)
    elif component_type == "profile":
        if 'chain_id' not in kwargs:
            return "<p>Error: chain_id required for profile component</p>"
        content = create_plotly_residue_profile(per_res_energy, **kwargs)
    elif component_type == "matrix":
        if 'cx' not in kwargs or 'chain_exprs' not in kwargs:
            return "<p>Error: cx and chain_exprs required for matrix component</p>"
        content = create_plotly_contact_matrix(**kwargs)
    else:
        return f"<p>Error: Unknown component type '{component_type}'</p>"
    
    return header + content