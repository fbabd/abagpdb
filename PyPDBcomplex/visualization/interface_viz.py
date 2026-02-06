from __future__ import annotations
from typing import Optional, List, Dict, Tuple
import json

# Try relative imports first (when used as part of abag2 package)
# Fall back to absolute imports (when used standalone)
try:
    from ..models import Complex
    from ..interface import InterfaceAnalysis, ResidueInterfaceData, compute_interface
except ImportError:
    # Standalone usage - expect types to be passed in
    # Type hints will be strings to avoid import errors
    Complex = 'Complex'
    InterfaceAnalysis = 'InterfaceAnalysis'
    ResidueInterfaceData = 'ResidueInterfaceData'
    compute_interface = None

def generate_interface_html(
    analysis: InterfaceAnalysis,
    cx: Complex,
    title: str = "Protein Interface Analysis",
    color_scheme: str = "default",
    box_size: int = 50 
) -> str:
    """
    Generate an interactive HTML visualization of interface contacts.
    
    Args:
        analysis: InterfaceAnalysis result
        cx: Complex structure
        title: Page title
        color_scheme: 'default', 'heatmap', or 'grouped'
    
    Returns:
        HTML string ready to save or display
    """
    
    # Color schemes
    schemes = {
        "default": {
            "interface": "#ff6b6b",
            "non_interface": "#e8e8e8",
            "group_A": "#4ecdc4",
            "group_B": "#ffe66d",
        },
        "heatmap": {
            # Will use gradient based on contact count
            "low": "#fff5f0",
            "high": "#7f0000",
            "non_interface": "#f0f0f0",
        },
        "grouped": {
            "group_A_interface": "#3498db",
            "group_B_interface": "#e74c3c",
            "both_groups": "#9b59b6",
            "non_interface": "#ecf0f1",
        }
    }
    
    colors = schemes.get(color_scheme, schemes["default"])
    
    # Build HTML sections
    html_chains = []
    
    for chain_id in sorted(analysis.chain_data.keys()):
        res_data_list = analysis.chain_data[chain_id]
        
        if not res_data_list:
            continue
        
        is_group_A = chain_id in analysis.group_A_chains
        is_group_B = chain_id in analysis.group_B_chains
        
        # Determine chain label
        if is_group_A and is_group_B:
            group_label = "Groups A & B"
        elif is_group_A:
            group_label = "Group A"
        elif is_group_B:
            group_label = "Group B"
        else:
            group_label = "Other"
        
        html_chains.append(
            _render_chain_table(
                chain_id, 
                res_data_list, 
                group_label,
                is_group_A,
                is_group_B,
                color_scheme,
                colors
            )
        )
    
    # Summary statistics
    summary_A = analysis.get_group_summary("A")
    summary_B = analysis.get_group_summary("B")
    
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
        .summary {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .summary-card {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 6px;
            border-left: 4px solid #3498db;
        }}
        .summary-card h3 {{
            margin: 0 0 10px 0;
            color: #2c3e50;
            font-size: 14px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        .summary-card .value {{
            font-size: 28px;
            font-weight: bold;
            color: #3498db;
        }}
        .summary-card .label {{
            font-size: 12px;
            color: #7f8c8d;
            margin-top: 5px;
        }}
        .chain-section {{
            margin: 30px 0;
        }}
        .chain-header {{
            background: #34495e;
            color: white;
            padding: 12px 15px;
            border-radius: 6px 6px 0 0;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}
        .chain-header h2 {{
            margin: 0;
            font-size: 18px;
        }}
        .chain-header .badge {{
            background: rgba(255,255,255,0.2);
            padding: 4px 12px;
            border-radius: 12px;
            font-size: 12px;
            font-weight: normal;
        }}
        .residue-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax({box_size}px, 1fr));
            gap: 4px;
            padding: 15px;
            background: #fafafa;
            border: 1px solid #ddd;
            border-top: none;
            border-radius: 0 0 6px 6px;
        }}
        .residue-cell {{
            aspect-ratio: 1;
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            border-radius: 4px;
            cursor: pointer;
            transition: all 0.2s;
            font-size: 11px;
            padding: 4px;
            position: relative;
        }}
        .residue-cell:hover {{
            transform: scale(1.15);
            z-index: 10;
            box-shadow: 0 4px 12px rgba(0,0,0,0.2);
        }}
        .residue-cell .resname {{
            font-weight: bold;
            font-size: 12px;
        }}
        .residue-cell .resseq {{
            font-size: 10px;
            opacity: 0.8;
        }}
        .residue-cell .contacts {{
            position: absolute;
            top: 2px;
            right: 2px;
            background: rgba(0,0,0,0.7);
            color: white;
            font-size: 9px;
            padding: 1px 4px;
            border-radius: 3px;
            font-weight: bold;
        }}
        .legend {{
            margin: 20px 0;
            padding: 15px;
            background: #f8f9fa;
            border-radius: 6px;
            display: flex;
            flex-wrap: wrap;
            gap: 20px;
            align-items: center;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 8px;
        }}
        .legend-box {{
            width: 30px;
            height: 30px;
            border-radius: 4px;
            border: 1px solid #ddd;
        }}
        .tooltip {{
            position: fixed;
            background: rgba(0,0,0,0.9);
            color: white;
            padding: 10px;
            border-radius: 6px;
            font-size: 12px;
            pointer-events: none;
            z-index: 1000;
            display: none;
            max-width: 300px;
        }}
        .tooltip.show {{
            display: block;
        }}
        .contact-list {{
            margin-top: 5px;
            padding-top: 5px;
            border-top: 1px solid rgba(255,255,255,0.3);
            font-size: 11px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>{title}</h1>
        
        <div class="summary">
            <div class="summary-card">
                <h3>Group A Chains</h3>
                <div class="value">{', '.join(summary_A['chains'])}</div>
                <div class="label">{summary_A['interface_residues']} interface residues ({summary_A['interface_percentage']:.1f}%)</div>
            </div>
            <div class="summary-card">
                <h3>Group B Chains</h3>
                <div class="value">{', '.join(summary_B['chains'])}</div>
                <div class="label">{summary_B['interface_residues']} interface residues ({summary_B['interface_percentage']:.1f}%)</div>
            </div>
            <div class="summary-card">
                <h3>Total Contacts</h3>
                <div class="value">{analysis.total_contacts}</div>
                <div class="label">Atom-atom contacts within {analysis.cutoff} Å</div>
            </div>
            {f'''<div class="summary-card">
                <h3>Buried Surface Area</h3>
                <div class="value">{analysis.bsa_total:.0f} Ų</div>
                <div class="label">{analysis.note}</div>
            </div>''' if analysis.bsa_total else ''}
        </div>
        
        <div class="legend">
            <strong>Legend:</strong>
            <div class="legend-item">
                <div class="legend-box" style="background: {colors.get('interface', colors.get('group_A_interface', '#ff6b6b'))}"></div>
                <span>Interface Residue</span>
            </div>
            <div class="legend-item">
                <div class="legend-box" style="background: {colors.get('non_interface', '#e8e8e8')}"></div>
                <span>Non-Interface</span>
            </div>
        </div>
        
        <div style="margin: 25px 0; padding: 20px; background: #f8f9fa; border-radius: 8px; border-left: 4px solid #3498db;">
            <div style="display: flex; align-items: center; gap: 20px; flex-wrap: wrap;">
                <label style="font-weight: bold; color: #2c3e50;">Box Size:</label>
                <input type="range" id="boxSizeSlider" min="30" max="100" value="{box_size}" 
                       style="flex: 1; min-width: 200px;">
                <span id="boxSizeValue" style="font-weight: bold; color: #3498db; min-width: 60px;">{box_size}px</span>
                <button onclick="resetBoxSize()" 
                        style="background: #3498db; color: white; border: none; padding: 8px 16px; 
                               border-radius: 6px; cursor: pointer; font-weight: 600;">
                    Reset
                </button>
            </div>
        </div>
        
        {''.join(html_chains)}
        
        <div class="tooltip" id="tooltip"></div>
    </div>
    
    <script>
        const tooltip = document.getElementById('tooltip');
        
        // Tooltip functionality
        document.querySelectorAll('.residue-cell').forEach(cell => {{
            cell.addEventListener('mouseenter', (e) => {{
                const info = JSON.parse(cell.dataset.info);
                let html = `
                    <div><strong>${{info.resname}} ${{info.chain}}:${{info.resseq}}</strong></div>
                    <div>Interface: ${{info.is_interface ? 'Yes' : 'No'}}</div>
                    <div>Contacts: ${{info.contact_count}}</div>
                `;
                
                if (info.partners && info.partners.length > 0) {{
                    html += `<div class="contact-list">Partners: ${{info.partners.join(', ')}}</div>`;
                }}
                
                tooltip.innerHTML = html;
                tooltip.classList.add('show');
            }});
            
            cell.addEventListener('mousemove', (e) => {{
                tooltip.style.left = (e.clientX + 15) + 'px';
                tooltip.style.top = (e.clientY + 15) + 'px';
            }});
            
            cell.addEventListener('mouseleave', () => {{
                tooltip.classList.remove('show');
            }});
        }});
        
        // Box size control functionality
        const boxSizeSlider = document.getElementById('boxSizeSlider');
        const boxSizeValue = document.getElementById('boxSizeValue');
        const residueGrids = document.querySelectorAll('.residue-grid');
        
        function updateBoxSize(size) {{
            residueGrids.forEach(grid => {{
                grid.style.gridTemplateColumns = `repeat(auto-fill, minmax(${{size}}px, 1fr))`;
            }});
            boxSizeValue.textContent = size + 'px';
        }}
        
        if (boxSizeSlider) {{
            boxSizeSlider.addEventListener('input', (e) => {{
                updateBoxSize(e.target.value);
            }});
        }}
        
        function resetBoxSize() {{
            const defaultSize = {box_size};
            if (boxSizeSlider) {{
                boxSizeSlider.value = defaultSize;
                updateBoxSize(defaultSize);
            }}
        }}
    </script>
</body>
</html>
"""
    
    return html


def _render_chain_table(
    chain_id: str,
    res_data_list: List[ResidueInterfaceData],
    group_label: str,
    is_group_A: bool,
    is_group_B: bool,
    color_scheme: str,
    colors: Dict
) -> str:
    """Render a single chain's residue grid."""
    
    # Calculate statistics
    total = len(res_data_list)
    interface_count = sum(1 for rd in res_data_list if rd.is_interface)
    total_contacts = sum(rd.contact_count for rd in res_data_list)
    
    # Build residue cells
    cells = []
    for rd in res_data_list:
        # Determine color
        if color_scheme == "heatmap":
            if rd.is_interface:
                # Gradient from low to high
                max_contacts = max((r.contact_count for r in res_data_list if r.is_interface), default=1)
                intensity = rd.contact_count / max_contacts if max_contacts > 0 else 0
                # Interpolate between low and high
                color = _interpolate_color(colors["low"], colors["high"], intensity)
            else:
                color = colors["non_interface"]
        elif color_scheme == "grouped":
            if rd.is_interface:
                if is_group_A and is_group_B:
                    color = colors["both_groups"]
                elif is_group_A:
                    color = colors["group_A_interface"]
                elif is_group_B:
                    color = colors["group_B_interface"]
                else:
                    color = "#95a5a6"
            else:
                color = colors["non_interface"]
        else:  # default
            color = colors["interface"] if rd.is_interface else colors["non_interface"]
        
        # Text color (dark on light backgrounds, light on dark)
        text_color = "#000" if _is_light_color(color) else "#fff"
        
        # Build cell data
        info = {
            "resname": rd.residue.resname,
            "chain": chain_id,
            "resseq": f"{rd.residue.resseq}{rd.residue.icode}".strip(),
            "is_interface": rd.is_interface,
            "contact_count": rd.contact_count,
            "partners": rd.partner_residues[:5]  # Limit to first 5 for tooltip
        }
        
        import json
        cell_html = f"""
        <div class="residue-cell" 
             style="background: {color}; color: {text_color};"
             data-info='{json.dumps(info)}'>
            <div class="resname">{rd.residue.resname}</div>
            <div class="resseq">{info['resseq']}</div>
            {f'<div class="contacts">{rd.contact_count}</div>' if rd.is_interface else ''}
        </div>
        """
        cells.append(cell_html)
    
    return f"""
    <div class="chain-section">
        <div class="chain-header">
            <h2>Chain {chain_id}</h2>
            <div>
                <span class="badge">{group_label}</span>
                <span class="badge">{interface_count}/{total} interface ({100*interface_count/total:.0f}%)</span>
                <span class="badge">{total_contacts} contacts</span>
            </div>
        </div>
        <div class="residue-grid">
            {''.join(cells)}
        </div>
    </div>
    """


def _interpolate_color(color1: str, color2: str, t: float) -> str:
    """Interpolate between two hex colors."""
    # Parse hex colors
    r1, g1, b1 = int(color1[1:3], 16), int(color1[3:5], 16), int(color1[5:7], 16)
    r2, g2, b2 = int(color2[1:3], 16), int(color2[3:5], 16), int(color2[5:7], 16)
    
    # Interpolate
    r = int(r1 + (r2 - r1) * t)
    g = int(g1 + (g2 - g1) * t)
    b = int(b1 + (b2 - b1) * t)
    
    return f"#{r:02x}{g:02x}{b:02x}"


def _is_light_color(hex_color: str) -> bool:
    """Determine if a color is light (for contrast)."""
    r, g, b = int(hex_color[1:3], 16), int(hex_color[3:5], 16), int(hex_color[5:7], 16)
    # Calculate perceived brightness
    brightness = (r * 299 + g * 587 + b * 114) / 1000
    return brightness > 128


# ============ Alternative: Simple Table Format ============

def generate_interface_table(
    analysis: InterfaceAnalysis,
    chains: Optional[List[str]] = None,
    show_non_interface: bool = False
) -> str:
    """
    Generate a simple text table of interface residues.
    
    Args:
        analysis: InterfaceAnalysis result
        chains: Specific chains to show (None = all)
        show_non_interface: Include non-interface residues
    
    Returns:
        Formatted text table
    """
    chains = chains or sorted(analysis.chain_data.keys())
    
    lines = []
    lines.append("=" * 90)
    lines.append(f"{'Chain':<8} {'Residue':<12} {'Type':<10} {'Contacts':<10} {'Partners':<40}")
    lines.append("=" * 90)
    
    for chain_id in chains:
        if chain_id not in analysis.chain_data:
            continue
        
        for rd in analysis.chain_data[chain_id]:
            if not show_non_interface and not rd.is_interface:
                continue
            
            status = "Interface" if rd.is_interface else "-"
            partners = ", ".join(rd.partner_residues[:3])
            if len(rd.partner_residues) > 3:
                partners += f" +{len(rd.partner_residues)-3} more"
            
            lines.append(
                f"{chain_id:<8} {rd.id_str:<12} {status:<10} {rd.contact_count:<10} {partners:<40}"
            )
    
    lines.append("=" * 90)
    return "\n".join(lines)



def generate_interface_html_with_3d(
    cx: Complex,
    analysis: InterfaceAnalysis,
    title: str = "Protein Interface Analysis with 3D Structure",
    color_scheme: str = "default",
    box_size: int = 50,
    show_3d: bool = True,
    cutoff_min: float = 3.0,
    cutoff_max: float = 8.0,
    cutoff_step: float = 0.5
) -> str:
    """
    Generate an interactive HTML visualization with both 2D grid and 3D structure views.
    
    Args:
        cx: Complex structure
        analysis: InterfaceAnalysis result
        title: Page title
        color_scheme: 'default', 'heatmap', or 'grouped'
        box_size: Initial size of residue boxes in 2D grid
        show_3d: Include 3D structure viewer
        cutoff_min: Minimum cutoff distance for slider
        cutoff_max: Maximum cutoff distance for slider
        cutoff_step: Step size for cutoff slider
    
    Returns:
        HTML string with interactive visualization
    """
    
    # Color schemes
    schemes = {
        "default": {
            "interface": "#ff6b6b",
            "non_interface": "#e8e8e8",
            "group_A": "#4ecdc4",
            "group_B": "#ffe66d",
            "group_A_3d": "0x4ecdc4",
            "group_B_3d": "0xffe66d",
            "interface_3d": "0xff6b6b",
        },
        "heatmap": {
            "low": "#fff5f0",
            "high": "#7f0000",
            "non_interface": "#f0f0f0",
            "interface_3d": "0xff0000",
            "interface": "#ff6b6b",
            "non_interface": "#e8e8e8",
        },
        "grouped": {
            "group_A_interface": "#3498db",
            "group_B_interface": "#e74c3c",
            "both_groups": "#9b59b6",
            "non_interface": "#ecf0f1",
            "group_A_3d": "0x3498db",
            "group_B_3d": "0xe74c3c",
            "interface": "#ff6b6b",
            "non_interface": "#e8e8e8",
        }
    }
    
    colors = schemes.get(color_scheme, schemes["default"])
    
    # Get PDB content for 3D viewer with fallback
    pdb_content = ""
    if hasattr(cx, 'to_pdb_string'):
        try:
            pdb_content = cx.to_pdb_string()
        except Exception as e:
            print(f"Warning: Could not generate PDB string: {e}")
    
    # If still no PDB content, try to generate minimal PDB from Complex object
    if not pdb_content:
        try:
            pdb_lines = []
            pdb_lines.append("HEADER    GENERATED STRUCTURE")
            atom_num = 1
            for chain_id, chain in cx.chains.items():
                for residue in chain.iter_residues():
                    for atom in residue.iter_atoms():
                        coord = atom.coord
                        pdb_lines.append(
                            f"ATOM  {atom_num:5d}  {atom.name:<4s}{residue.resname:>3s} "
                            f"{chain_id}{residue.resseq:4d}{residue.icode}   "
                            f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}  1.00 20.00           "
                            f"{atom.element:>2s}"
                        )
                        atom_num += 1
            pdb_lines.append("END")
            pdb_content = "\n".join(pdb_lines)
        except Exception as e:
            print(f"Warning: Could not generate minimal PDB: {e}")
            pdb_content = "HEADER    NO STRUCTURE DATA\nEND\n"
    
    # Build HTML sections for 2D grid
    html_chains = []
    for chain_id in sorted(analysis.chain_data.keys()):
        res_data_list = analysis.chain_data[chain_id]
        
        if not res_data_list:
            continue
        
        is_group_A = chain_id in analysis.group_A_chains
        is_group_B = chain_id in analysis.group_B_chains
        
        # Determine chain label
        if is_group_A and is_group_B:
            group_label = "Groups A & B"
        elif is_group_A:
            group_label = "Group A"
        elif is_group_B:
            group_label = "Group B"
        else:
            group_label = "Other"
        
        html_chains.append(
            _render_chain_table(
                chain_id, 
                res_data_list, 
                group_label,
                is_group_A,
                is_group_B,
                color_scheme,
                colors
            )
        )
    
    # Summary statistics
    summary_A = analysis.get_group_summary("A")
    summary_B = analysis.get_group_summary("B")
    
    # Prepare interface residue data for 3D highlighting
    interface_residues_json = _prepare_interface_residues_for_3d(analysis)
    
    # Prepare complex data for recomputation
    complex_data = {
        "group_A_chains": analysis.group_A_chains,
        "group_B_chains": analysis.group_B_chains,
        "current_cutoff": analysis.cutoff
    }
    
    # Escape backticks in PDB content for JavaScript template literal
    pdb_content_escaped = pdb_content.replace('`', '\\`')
    
    html = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>{title}</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.1/3Dmol-min.js"></script>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f5f5f5;
        }}
        .container {{
            max-width: 1800px;
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
        
        /* 3D Viewer Section */
        .viewer-section {{
            margin: 20px 0;
            background: #f8f9fa;
            border-radius: 8px;
            padding: 20px;
        }}
        .viewer-controls {{
            display: flex;
            gap: 20px;
            align-items: center;
            margin-bottom: 15px;
            flex-wrap: wrap;
        }}
        .control-group {{
            display: flex;
            flex-direction: column;
            gap: 5px;
        }}
        .control-group label {{
            font-size: 12px;
            font-weight: 600;
            color: #2c3e50;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        #viewer-3d {{
            width: 100%;
            height: 600px;
            background: #000;
            border-radius: 6px;
            position: relative;
        }}
        .slider-container {{
            display: flex;
            flex-direction: column;
            gap: 5px;
            min-width: 250px;
        }}
        .slider-row {{
            display: flex;
            align-items: center;
            gap: 10px;
        }}
        input[type="range"] {{
            flex: 1;
            min-width: 150px;
        }}
        .slider-value {{
            min-width: 60px;
            font-weight: bold;
            color: #3498db;
            font-size: 14px;
        }}
        button {{
            background: #3498db;
            color: white;
            border: none;
            padding: 8px 16px;
            border-radius: 6px;
            cursor: pointer;
            font-weight: 600;
            transition: background 0.2s;
        }}
        button:hover {{
            background: #2980b9;
        }}
        button:disabled {{
            background: #95a5a6;
            cursor: not-allowed;
        }}
        .spinner {{
            display: none;
            width: 20px;
            height: 20px;
            border: 3px solid #f3f3f3;
            border-top: 3px solid #3498db;
            border-radius: 50%;
            animation: spin 1s linear infinite;
            margin-left: 10px;
        }}
        @keyframes spin {{
            0% {{ transform: rotate(0deg); }}
            100% {{ transform: rotate(360deg); }}
        }}
        
        /* Summary Section */
        .summary {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .summary-card {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 6px;
            border-left: 4px solid #3498db;
        }}
        .summary-card h3 {{
            margin: 0 0 10px 0;
            color: #2c3e50;
            font-size: 14px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        .summary-card .value {{
            font-size: 28px;
            font-weight: bold;
            color: #3498db;
        }}
        .summary-card .label {{
            font-size: 12px;
            color: #7f8c8d;
            margin-top: 5px;
        }}
        
        /* Chain Grid Section */
        .chain-section {{
            margin: 30px 0;
        }}
        .chain-header {{
            background: #34495e;
            color: white;
            padding: 12px 15px;
            border-radius: 6px 6px 0 0;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}
        .chain-header h2 {{
            margin: 0;
            font-size: 18px;
        }}
        .chain-header .badge {{
            background: rgba(255,255,255,0.2);
            padding: 4px 12px;
            border-radius: 12px;
            font-size: 12px;
            font-weight: normal;
        }}
        .residue-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax({box_size}px, 1fr));
            gap: 4px;
            padding: 15px;
            background: #fafafa;
            border: 1px solid #ddd;
            border-top: none;
            border-radius: 0 0 6px 6px;
        }}
        .residue-cell {{
            aspect-ratio: 1;
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            border-radius: 4px;
            cursor: pointer;
            transition: all 0.2s;
            font-size: 11px;
            padding: 4px;
            position: relative;
        }}
        .residue-cell:hover {{
            transform: scale(1.15);
            z-index: 10;
            box-shadow: 0 4px 12px rgba(0,0,0,0.2);
        }}
        .residue-cell .resname {{
            font-weight: bold;
            font-size: 12px;
        }}
        .residue-cell .resseq {{
            font-size: 10px;
            opacity: 0.8;
        }}
        .residue-cell .contacts {{
            position: absolute;
            top: 2px;
            right: 2px;
            background: rgba(0,0,0,0.7);
            color: white;
            font-size: 9px;
            padding: 1px 4px;
            border-radius: 3px;
            font-weight: bold;
        }}
        .legend {{
            margin: 20px 0;
            padding: 15px;
            background: #f8f9fa;
            border-radius: 6px;
            display: flex;
            flex-wrap: wrap;
            gap: 20px;
            align-items: center;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 8px;
        }}
        .legend-box {{
            width: 30px;
            height: 30px;
            border-radius: 4px;
            border: 1px solid #ddd;
        }}
        .tooltip {{
            position: fixed;
            background: rgba(0,0,0,0.9);
            color: white;
            padding: 10px;
            border-radius: 6px;
            font-size: 12px;
            pointer-events: none;
            z-index: 1000;
            display: none;
            max-width: 300px;
        }}
        .tooltip.show {{
            display: block;
        }}
        .contact-list {{
            margin-top: 5px;
            padding-top: 5px;
            border-top: 1px solid rgba(255,255,255,0.3);
            font-size: 11px;
        }}
        
        .section-toggle {{
            background: #34495e;
            color: white;
            padding: 12px 15px;
            border-radius: 6px;
            cursor: pointer;
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-top: 20px;
        }}
        .section-toggle:hover {{
            background: #2c3e50;
        }}
        .section-content {{
            max-height: 2000px;
            overflow: hidden;
            transition: max-height 0.3s ease-out;
        }}
        .section-content.collapsed {{
            max-height: 0;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>{title}</h1>
        
        <div class="summary">
            <div class="summary-card">
                <h3>Total Contacts</h3>
                <div class="value" id="total-contacts">{analysis.total_contacts}</div>
                <div class="label">Atom-atom contacts at {analysis.cutoff}Å cutoff</div>
            </div>
            
            <div class="summary-card">
                <h3>Group A Interface</h3>
                <div class="value">{summary_A['interface_residues']}</div>
                <div class="label">
                    {summary_A['interface_percentage']:.1f}% of {summary_A['total_residues']} residues<br>
                    Chains: {', '.join(summary_A['chains'])}
                </div>
            </div>
            
            <div class="summary-card">
                <h3>Group B Interface</h3>
                <div class="value">{summary_B['interface_residues']}</div>
                <div class="label">
                    {summary_B['interface_percentage']:.1f}% of {summary_B['total_residues']} residues<br>
                    Chains: {', '.join(summary_B['chains'])}
                </div>
            </div>
            
            {f'''<div class="summary-card">
                <h3>Buried Surface Area</h3>
                <div class="value">{analysis.bsa_total:.0f}</div>
                <div class="label">Å² ({analysis.note})</div>
            </div>''' if analysis.bsa_total else ''}
        </div>
        
        {_render_3d_viewer_section(cutoff_min, cutoff_max, cutoff_step, analysis.cutoff) if show_3d else ''}
        
        <div class="section-toggle" onclick="toggleSection('grid-section')">
            <h2 style="margin: 0;">2D Residue Grid View</h2>
            <span id="grid-toggle-icon">▼</span>
        </div>
        
        <div id="grid-section" class="section-content">
            <div class="legend">
                <div style="font-weight: 600; margin-right: 10px;">Legend:</div>
                <div class="legend-item">
                    <div class="legend-box" style="background: {colors['interface']};"></div>
                    <span>Interface Residue</span>
                </div>
                <div class="legend-item">
                    <div class="legend-box" style="background: {colors['non_interface']};"></div>
                    <span>Non-Interface</span>
                </div>
            </div>
            
            <div style="margin: 20px 0;">
                <label style="font-weight: 600; margin-right: 10px;">Box Size:</label>
                <input type="range" id="boxSizeSlider" min="30" max="100" value="{box_size}" 
                       style="width: 200px; vertical-align: middle;">
                <span id="boxSizeValue" style="margin-left: 10px; font-weight: bold;">{box_size}px</span>
                <button onclick="resetBoxSize()" 
                        style="margin-left: 15px;">
                    Reset
                </button>
            </div>
            
            {''.join(html_chains)}
        </div>
        
        <div class="tooltip" id="tooltip"></div>
    </div>
    
    <script>
        // Global variables
        let viewer3d = null;
        let pdbData = `{pdb_content_escaped}`;
        let interfaceResidues = {interface_residues_json};
        let complexData = {json.dumps(complex_data)};
        let currentCutoff = {analysis.cutoff};
        
        // Initialize 3D viewer
        {_render_3d_viewer_javascript() if show_3d else ''}
        
        // Tooltip functionality for 2D grid
        const tooltip = document.getElementById('tooltip');
        
        document.querySelectorAll('.residue-cell').forEach(cell => {{
            cell.addEventListener('mouseenter', (e) => {{
                const info = JSON.parse(cell.dataset.info);
                let html = `
                    <div><strong>${{info.resname}} ${{info.chain}}:${{info.resseq}}</strong></div>
                    <div>Interface: ${{info.is_interface ? 'Yes' : 'No'}}</div>
                    <div>Contacts: ${{info.contact_count}}</div>
                `;
                
                if (info.partners && info.partners.length > 0) {{
                    html += `<div class="contact-list">Partners: ${{info.partners.join(', ')}}</div>`;
                }}
                
                tooltip.innerHTML = html;
                tooltip.classList.add('show');
            }});
            
            cell.addEventListener('mousemove', (e) => {{
                tooltip.style.left = (e.clientX + 15) + 'px';
                tooltip.style.top = (e.clientY + 15) + 'px';
            }});
            
            cell.addEventListener('mouseleave', () => {{
                tooltip.classList.remove('show');
            }});
        }});
        
        // Box size control for 2D grid
        const boxSizeSlider = document.getElementById('boxSizeSlider');
        const boxSizeValue = document.getElementById('boxSizeValue');
        const residueGrids = document.querySelectorAll('.residue-grid');
        
        function updateBoxSize(size) {{
            residueGrids.forEach(grid => {{
                grid.style.gridTemplateColumns = `repeat(auto-fill, minmax(${{size}}px, 1fr))`;
            }});
            boxSizeValue.textContent = size + 'px';
        }}
        
        if (boxSizeSlider) {{
            boxSizeSlider.addEventListener('input', (e) => {{
                updateBoxSize(e.target.value);
            }});
        }}
        
        function resetBoxSize() {{
            const defaultSize = {box_size};
            if (boxSizeSlider) {{
                boxSizeSlider.value = defaultSize;
                updateBoxSize(defaultSize);
            }}
        }}
        
        // Section toggle functionality
        function toggleSection(sectionId) {{
            const section = document.getElementById(sectionId);
            const icon = document.getElementById(sectionId.replace('-section', '-toggle-icon'));
            
            if (section.classList.contains('collapsed')) {{
                section.classList.remove('collapsed');
                icon.textContent = '▼';
            }} else {{
                section.classList.add('collapsed');
                icon.textContent = '▶';
            }}
        }}
    </script>
</body>
</html>
"""
    
    return html


def _render_3d_viewer_section(cutoff_min: float, cutoff_max: float, cutoff_step: float, current_cutoff: float) -> str:
    """Render the 3D viewer section HTML"""
    return f"""
        <div class="viewer-section">
            <h2 style="margin-top: 0;">3D Structure Viewer</h2>
            
            <div class="viewer-controls">
                
                <div class="control-group">
                    <label>View Style</label>
                    <div style="display: flex; gap: 10px;">
                        <button onclick="setStyle('cartoon')">Cartoon</button>
                        <button onclick="setStyle('stick')">Stick</button>
                        <button onclick="setStyle('sphere')">Sphere</button>
                    </div>
                </div>
                
                <div class="control-group">
                    <label>Actions</label>
                    <div style="display: flex; gap: 10px;">
                        <button onclick="resetView()">Reset View</button>
                        <button onclick="toggleSpin()">Toggle Spin</button>
                    </div>
                </div>
            </div>
            
            <div id="viewer-3d"></div>
        </div>
    """


def _render_3d_viewer_javascript() -> str:
    """Render the JavaScript code for 3D viewer functionality"""
    return """
        // Initialize 3D viewer when DOM is ready
        document.addEventListener('DOMContentLoaded', function() {
            initViewer3D();
        });
        
        function initViewer3D() {
            console.log('Initializing 3D viewer...');
            const element = document.getElementById('viewer-3d');
            
            if (!element) {
                console.error('3D viewer element not found!');
                return;
            }
            
            if (typeof $3Dmol === 'undefined') {
                console.error('3Dmol.js library not loaded!');
                alert('Error: 3Dmol.js library failed to load. Check your internet connection.');
                return;
            }
            
            if (!pdbData || pdbData.trim() === '' || pdbData === 'HEADER    NO STRUCTURE DATA') {
                console.error('No PDB data available!');
                element.innerHTML = '<div style="color: white; padding: 20px; text-align: center;">' +
                                   'No structure data available<br>PDB file may be empty or failed to load</div>';
                return;
            }
            
            console.log('PDB data length:', pdbData.length, 'characters');
            console.log('Interface residues:', Object.keys(interfaceResidues).length, 'chains');
            
            const config = { backgroundColor: 'white' };
            viewer3d = $3Dmol.createViewer(element, config);
            
            // Load PDB structure
            try {
                viewer3d.addModel(pdbData, 'pdb');
                console.log('PDB model loaded successfully');
            } catch (error) {
                console.error('Error loading PDB model:', error);
                element.innerHTML = '<div style="color: white; padding: 20px; text-align: center;">' +
                                   'Error loading structure<br>' + error.message + '</div>';
                return;
            }
            
            // Apply initial styling
            applyStructureStyle('cartoon');
            
            viewer3d.zoomTo();
            viewer3d.render();
            console.log('3D viewer initialized successfully');
        }
        
        function applyStructureStyle(style) {
            if (!viewer3d) return;
            
            viewer3d.setStyle({}, {});  // Clear all styles
            
            // Style for non-interface residues (gray)
            if (style === 'cartoon') {
                viewer3d.setStyle({}, {cartoon: {color: 'gray', opacity: 0.7}});
            } else if (style === 'stick') {
                viewer3d.setStyle({}, {stick: {color: 'gray', opacity: 0.7}});
            } else if (style === 'sphere') {
                viewer3d.setStyle({}, {sphere: {color: 'gray', opacity: 0.7, scale: 0.3}});
            }
            
            // Highlight interface residues
            highlightInterfaceResidues(style);
            
            viewer3d.render();
        }
        
        function highlightInterfaceResidues(style) {
            if (!viewer3d || !interfaceResidues) return;
            
            const styleMap = {
                'cartoon': (color) => ({cartoon: {color: color, opacity: 1.0}}),
                'stick': (color) => ({stick: {color: color, opacity: 1.0, radius: 0.25}}),
                'sphere': (color) => ({sphere: {color: color, opacity: 1.0, scale: 0.5}})
            };
            
            const getStyle = styleMap[style] || styleMap['cartoon'];
            
            // Highlight each interface residue
            Object.keys(interfaceResidues).forEach(chainId => {
                const residues = interfaceResidues[chainId];
                
                residues.forEach(res => {
                    const selection = {
                        chain: chainId,
                        resi: res.resseq
                    };
                    
                    // Use different colors for different groups
                    let color = 'red';  // Default interface color
                    if (complexData.group_A_chains.includes(chainId)) {
                        color = 'cyan';
                    } else if (complexData.group_B_chains.includes(chainId)) {
                        color = 'yellow';
                    }
                    
                    viewer3d.setStyle(selection, getStyle(color));
                });
            });
        }
        
        function setStyle(style) {
            applyStructureStyle(style);
        }
        
        function resetView() {
            if (!viewer3d) return;
            viewer3d.zoomTo();
            viewer3d.render();
        }
        
        let isSpinning = false;
        function toggleSpin() {
            isSpinning = !isSpinning;
            if (isSpinning) {
                viewer3d.spin(true);
            } else {
                viewer3d.spin(false);
            }
        }
        
        
    """


def _prepare_interface_residues_for_3d(analysis: InterfaceAnalysis) -> str:
    """Prepare interface residue data for 3D visualization as JSON"""
    interface_data = {}
    
    for chain_id, res_data_list in analysis.chain_data.items():
        interface_residues = [
            {
                "resseq": rd.residue.resseq,
                "resname": rd.residue.resname,
                "icode": rd.residue.icode,
                "contact_count": rd.contact_count
            }
            for rd in res_data_list if rd.is_interface
        ]
        
        if interface_residues:
            interface_data[chain_id] = interface_residues
    
    return json.dumps(interface_data)


def _render_chain_table(
    chain_id: str,
    res_data_list: List[ResidueInterfaceData],
    group_label: str,
    is_group_A: bool,
    is_group_B: bool,
    color_scheme: str,
    colors: Dict
) -> str:
    """Render a single chain's residue grid (same as original)."""
    
    # Calculate statistics
    total = len(res_data_list)
    interface_count = sum(1 for rd in res_data_list if rd.is_interface)
    total_contacts = sum(rd.contact_count for rd in res_data_list)
    
    # Build residue cells
    cells = []
    for rd in res_data_list:
        # Determine color
        if color_scheme == "heatmap":
            if rd.is_interface:
                max_contacts = max((r.contact_count for r in res_data_list if r.is_interface), default=1)
                intensity = rd.contact_count / max_contacts if max_contacts > 0 else 0
                color = _interpolate_color(colors["low"], colors["high"], intensity)
            else:
                color = colors["non_interface"]
        elif color_scheme == "grouped":
            if rd.is_interface:
                if is_group_A and is_group_B:
                    color = colors["both_groups"]
                elif is_group_A:
                    color = colors["group_A_interface"]
                elif is_group_B:
                    color = colors["group_B_interface"]
                else:
                    color = "#95a5a6"
            else:
                color = colors["non_interface"]
        else:  # default
            color = colors["interface"] if rd.is_interface else colors["non_interface"]
        
        # Text color (dark on light backgrounds, light on dark)
        text_color = "#000" if _is_light_color(color) else "#fff"
        
        # Build cell data
        info = {
            "resname": rd.residue.resname,
            "chain": chain_id,
            "resseq": f"{rd.residue.resseq}{rd.residue.icode}".strip(),
            "is_interface": rd.is_interface,
            "contact_count": rd.contact_count,
            "partners": rd.partner_residues[:5]
        }
        
        cell_html = f"""
        <div class="residue-cell" 
             style="background: {color}; color: {text_color};"
             data-info='{json.dumps(info)}'>
            <div class="resname">{rd.residue.resname}</div>
            <div class="resseq">{info['resseq']}</div>
            {f'<div class="contacts">{rd.contact_count}</div>' if rd.is_interface else ''}
        </div>
        """
        cells.append(cell_html)
    
    return f"""
    <div class="chain-section">
        <div class="chain-header">
            <h2>Chain {chain_id}</h2>
            <div>
                <span class="badge">{group_label}</span>
                <span class="badge">{interface_count}/{total} interface ({100*interface_count/total:.0f}%)</span>
                <span class="badge">{total_contacts} contacts</span>
            </div>
        </div>
        <div class="residue-grid">
            {''.join(cells)}
        </div>
    </div>
    """


def _interpolate_color(color1: str, color2: str, t: float) -> str:
    """Interpolate between two hex colors."""
    r1, g1, b1 = int(color1[1:3], 16), int(color1[3:5], 16), int(color1[5:7], 16)
    r2, g2, b2 = int(color2[1:3], 16), int(color2[3:5], 16), int(color2[5:7], 16)
    
    r = int(r1 + (r2 - r1) * t)
    g = int(g1 + (g2 - g1) * t)
    b = int(b1 + (b2 - b1) * t)
    
    return f"#{r:02x}{g:02x}{b:02x}"


def _is_light_color(hex_color: str) -> bool:
    """Determine if a color is light (for contrast)."""
    r, g, b = int(hex_color[1:3], 16), int(hex_color[3:5], 16), int(hex_color[5:7], 16)
    brightness = (r * 299 + g * 587 + b * 114) / 1000
    return brightness > 128

