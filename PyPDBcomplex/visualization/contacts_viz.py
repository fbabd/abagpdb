from __future__ import annotations
from typing import Optional, List, Dict, Tuple, Set
from collections import defaultdict
import json

from ..models import Complex
from ..contacts import ContactAnalysis, ContactType, Contact
from ..interface import InterfaceAnalysis


def generate_contact_html(
    analysis: ContactAnalysis,
    cx: Complex,
    interface_analysis: Optional[InterfaceAnalysis] = None,
    title: str = "Protein Contact Analysis",
    color_scheme: str = "by_type",
    box_size: int = 60
) -> str:
    """
    Generate interactive HTML visualization of molecular contacts.
    
    Args:
        analysis: ContactAnalysis result
        cx: Complex structure
        interface_analysis: Optional InterfaceAnalysis to highlight interface regions
        title: Page title
        color_scheme: 'by_type' (color by contact type) or 'by_count' (heatmap)
    
    Returns:
        HTML string ready to save or display
    """
    
    # Contact type colors
    contact_colors = {
        ContactType.HYDROGEN_BOND: "#3498db",      # Blue
        ContactType.SALT_BRIDGE: "#e74c3c",        # Red
        ContactType.DISULFIDE: "#f39c12",          # Orange
        ContactType.HYDROPHOBIC: "#2ecc71",        # Green
        ContactType.PI_STACKING: "#9b59b6",        # Purple
    }
    
    # Get summaries
    residue_summary = analysis.get_residue_summary()
    contact_counts = analysis.get_contact_counts()
    
    # Build chain sections
    html_chains = []
    
    for chain_id in sorted(cx.chains.keys()):
        chain = cx.chains[chain_id]
        residues = list(chain.iter_residues())
        
        if not residues:
            continue
        
        html_chains.append(
            _render_contact_chain(
                chain_id,
                residues,
                residue_summary,
                analysis,
                interface_analysis,
                color_scheme,
                contact_colors
            )
        )
    
    # Generate legend
    legend_html = _generate_legend(contact_colors, color_scheme) 
    
    size_control_html = f"""
        <div style="margin: 25px 0; padding: 20px; background: #f8f9fa; border-radius: 8px; border-left: 4px solid #667eea;">
            <div style="display: flex; align-items: center; gap: 20px; flex-wrap: wrap;">
                <label style="font-weight: bold; color: #2c3e50;">Box Size:</label>
                <input type="range" id="boxSizeSlider" min="40" max="120" value="{box_size}" 
                    style="flex: 1; min-width: 200px;">
                <span id="boxSizeValue" style="font-weight: bold; color: #667eea; min-width: 60px;">{box_size}px</span>
                <button onclick="resetBoxSize()" 
                        style="background: #667eea; color: white; border: none; padding: 8px 16px; 
                            border-radius: 6px; cursor: pointer; font-weight: 600;">
                    Reset
                </button>
            </div>
        </div>
        """
    
    # Summary statistics
    summary_html = _generate_summary_cards(contact_counts, analysis)
    
    # Main HTML
    html = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>{title}</title>
    <style>
        * {{
            box-sizing: border-box;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f5f7fa;
        }}
        .container {{
            max-width: 1600px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 12px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.07);
        }}
        h1 {{
            color: #2c3e50;
            margin-bottom: 10px;
            font-size: 28px;
        }}
        .subtitle {{
            color: #7f8c8d;
            margin-bottom: 25px;
        }}
        
        /* Summary Cards */
        .summary {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 15px;
            margin: 25px 0;
        }}
        .summary-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
            border-radius: 10px;
            color: white;
        }}
        .summary-card.type-hbond {{
            background: linear-gradient(135deg, #3498db 0%, #2980b9 100%);
        }}
        .summary-card.type-salt {{
            background: linear-gradient(135deg, #e74c3c 0%, #c0392b 100%);
        }}
        .summary-card.type-disulfide {{
            background: linear-gradient(135deg, #f39c12 0%, #d68910 100%);
        }}
        .summary-card.type-hydrophobic {{
            background: linear-gradient(135deg, #2ecc71 0%, #27ae60 100%);
        }}
        .summary-card.type-pi {{
            background: linear-gradient(135deg, #9b59b6 0%, #8e44ad 100%);
        }}
        .summary-card h3 {{
            margin: 0 0 8px 0;
            font-size: 13px;
            text-transform: uppercase;
            letter-spacing: 1px;
            opacity: 0.9;
        }}
        .summary-card .value {{
            font-size: 32px;
            font-weight: bold;
            margin: 0;
        }}
        
        /* Legend */
        .legend {{
            margin: 25px 0;
            padding: 20px;
            background: #f8f9fa;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }}
        .legend-title {{
            font-weight: bold;
            margin-bottom: 12px;
            color: #2c3e50;
        }}
        .legend-items {{
            display: flex;
            flex-wrap: wrap;
            gap: 20px;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 8px;
            font-size: 14px;
        }}
        .legend-box {{
            width: 24px;
            height: 24px;
            border-radius: 4px;
            border: 2px solid rgba(0,0,0,0.1);
        }}
        
        /* Chain Section */
        .chain-section {{
            margin: 35px 0;
        }}
        .chain-header {{
            background: linear-gradient(135deg, #34495e 0%, #2c3e50 100%);
            color: white;
            padding: 15px 20px;
            border-radius: 8px 8px 0 0;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}
        .chain-header h2 {{
            margin: 0;
            font-size: 20px;
            font-weight: 600;
        }}
        .chain-badges {{
            display: flex;
            gap: 10px;
        }}
        .badge {{
            background: rgba(255,255,255,0.2);
            padding: 5px 14px;
            border-radius: 20px;
            font-size: 13px;
            font-weight: 500;
        }}
        .badge.interface {{
            background: rgba(46, 204, 113, 0.3);
        }}
        
        /* Residue Grid */
        .residue-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax({box_size}px, 1fr));
            gap: 6px;
            padding: 20px;
            background: #fafbfc;
            border: 1px solid #e1e4e8;
            border-top: none;
            border-radius: 0 0 8px 8px;
        }}
        .residue-cell {{
            aspect-ratio: 1;
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            border-radius: 6px;
            cursor: pointer;
            transition: all 0.2s cubic-bezier(0.4, 0, 0.2, 1);
            font-size: 11px;
            padding: 6px;
            position: relative;
            border: 2px solid transparent;
            box-shadow: 0 1px 3px rgba(0,0,0,0.08);
        }}
        .residue-cell:hover {{
            transform: scale(1.2) translateY(-4px);
            z-index: 100;
            box-shadow: 0 8px 20px rgba(0,0,0,0.15);
            border-color: #667eea;
        }}
        .residue-cell.has-contacts {{
            border-color: rgba(0,0,0,0.15);
        }}
        .residue-cell.interface {{
            border: 3px solid #2ecc71;
            box-shadow: 0 2px 8px rgba(46, 204, 113, 0.3);
        }}
        .residue-cell .resname {{
            font-weight: bold;
            font-size: 13px;
            margin-bottom: 2px;
        }}
        .residue-cell .resseq {{
            font-size: 10px;
            opacity: 0.7;
        }}
        .residue-cell .contact-badge {{
            position: absolute;
            top: 3px;
            right: 3px;
            background: rgba(0,0,0,0.75);
            color: white;
            font-size: 10px;
            padding: 2px 5px;
            border-radius: 4px;
            font-weight: bold;
        }}
        .residue-cell .contact-indicators {{
            position: absolute;
            bottom: 3px;
            left: 3px;
            right: 3px;
            display: flex;
            gap: 2px;
            justify-content: center;
        }}
        .contact-dot {{
            width: 8px;
            height: 8px;
            border-radius: 50%;
            border: 1.5px solid rgba(255,255,255,0.9);
            box-shadow: 0 1px 2px rgba(0,0,0,0.3);
        }}
        
        /* Multi-color gradient background for cells with multiple contact types */
        .residue-cell.multi-contact {{
            background: linear-gradient(135deg, var(--color1) 0%, var(--color2) 50%, var(--color3) 100%) !important;
        }}
        .residue-cell.multi-contact-2 {{
            background: linear-gradient(135deg, var(--color1) 0%, var(--color2) 100%) !important;
        }}
        .residue-cell.multi-contact-3 {{
            background: linear-gradient(120deg, var(--color1) 0%, var(--color2) 40%, var(--color3) 100%) !important;
        }}
        .residue-cell.multi-contact-4 {{
            background: linear-gradient(135deg, var(--color1) 0%, var(--color2) 33%, var(--color3) 66%, var(--color4) 100%) !important;
        }}
        .residue-cell.multi-contact-5 {{
            background: conic-gradient(var(--color1) 0deg 72deg, var(--color2) 72deg 144deg, var(--color3) 144deg 216deg, var(--color4) 216deg 288deg, var(--color5) 288deg 360deg) !important;
        }}
        
        /* Tooltip */
        .tooltip {{
            position: fixed;
            background: rgba(0,0,0,0.92);
            color: white;
            padding: 14px;
            border-radius: 8px;
            font-size: 13px;
            pointer-events: none;
            z-index: 1000;
            display: none;
            max-width: 350px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.3);
        }}
        .tooltip.show {{
            display: block;
        }}
        .tooltip-header {{
            font-size: 15px;
            font-weight: bold;
            margin-bottom: 10px;
            padding-bottom: 8px;
            border-bottom: 1px solid rgba(255,255,255,0.2);
        }}
        .tooltip-section {{
            margin: 8px 0;
        }}
        .tooltip-section-title {{
            font-size: 11px;
            text-transform: uppercase;
            color: #aaa;
            margin-bottom: 5px;
            letter-spacing: 0.5px;
        }}
        .contact-item {{
            padding: 4px 0;
            display: flex;
            align-items: center;
            gap: 8px;
        }}
        .contact-type-badge {{
            padding: 2px 6px;
            border-radius: 3px;
            font-size: 10px;
            font-weight: bold;
        }}
        .partner-list {{
            margin-top: 4px;
        }}
        .partner-item {{
            padding: 3px 0;
            font-size: 12px;
        }}
        
        /* Contact Network View */
        .network-toggle {{
            margin: 20px 0;
            text-align: center;
        }}
        .network-toggle button {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border: none;
            padding: 12px 24px;
            border-radius: 6px;
            cursor: pointer;
            font-size: 14px;
            font-weight: 600;
            transition: all 0.2s;
        }}
        .network-toggle button:hover {{
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(102, 126, 234, 0.4);
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>{title}</h1>
        <div class="subtitle">Molecular contact analysis with {len(analysis.contacts)} total interactions</div>
        
        {summary_html}
        {legend_html}
        
        {size_control_html}
        
        {''.join(html_chains)}
        
        <div class="tooltip" id="tooltip"></div>
    </div>
    
    <script>
        const tooltip = document.getElementById('tooltip');
        const contactColors = {json.dumps({k.value: v for k, v in contact_colors.items()})};
        
        // Tooltip functionality
        document.querySelectorAll('.residue-cell').forEach(cell => {{
            cell.addEventListener('mouseenter', (e) => {{
                const info = JSON.parse(cell.dataset.info);
                
                let html = `<div class="tooltip-header">${{info.resname}} ${{info.chain}}:${{info.resseq}}</div>`;
                
                if (info.is_interface) {{
                    html += `<div class="tooltip-section">
                        <div class="tooltip-section-title">⭐ Interface Residue</div>
                    </div>`;
                }}
                
                if (info.contacts && info.contacts.length > 0) {{
                    html += `<div class="tooltip-section">
                        <div class="tooltip-section-title">Contacts (${{info.contacts.length}})</div>`;
                    
                    const contactsByType = {{}};
                    info.contacts.forEach(c => {{
                        if (!contactsByType[c.type]) contactsByType[c.type] = [];
                        contactsByType[c.type].push(c);
                    }});
                    
                    for (const [type, contacts] of Object.entries(contactsByType)) {{
                        const color = contactColors[type] || '#999';
                        html += `<div class="contact-item">
                            <span class="contact-type-badge" style="background: ${{color}}">
                                ${{type.replace('_', ' ').toUpperCase()}}
                            </span>
                            <span>${{contacts.length}}×</span>
                        </div>`;
                    }}
                    
                    html += `</div>`;
                    
                    // Show partners
                    if (info.partners && info.partners.length > 0) {{
                        html += `<div class="tooltip-section">
                            <div class="tooltip-section-title">Partner Residues</div>
                            <div class="partner-list">`;
                        
                        info.partners.slice(0, 5).forEach(p => {{
                            html += `<div class="partner-item">→ ${{p}}</div>`;
                        }});
                        
                        if (info.partners.length > 5) {{
                            html += `<div class="partner-item">... and ${{info.partners.length - 5}} more</div>`;
                        }}
                        
                        html += `</div></div>`;
                    }}
                }} else {{
                    html += `<div class="tooltip-section">No contacts detected</div>`;
                }}
                
                tooltip.innerHTML = html;
                tooltip.classList.add('show');
            }});
            
            cell.addEventListener('mousemove', (e) => {{
                const offset = 15;
                let x = e.clientX + offset;
                let y = e.clientY + offset;
                
                // Keep tooltip in viewport
                const rect = tooltip.getBoundingClientRect();
                if (x + rect.width > window.innerWidth) {{
                    x = e.clientX - rect.width - offset;
                }}
                if (y + rect.height > window.innerHeight) {{
                    y = e.clientY - rect.height - offset;
                }}
                
                tooltip.style.left = x + 'px';
                tooltip.style.top = y + 'px';
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


def _generate_summary_cards(contact_counts: Dict, analysis: ContactAnalysis) -> str:
    """Generate summary statistic cards."""
    cards = []
    
    type_map = {
        ContactType.HYDROGEN_BOND: ("H-Bonds", "type-hbond"),
        ContactType.SALT_BRIDGE: ("Salt Bridges", "type-salt"),
        ContactType.DISULFIDE: ("Disulfides", "type-disulfide"),
        ContactType.HYDROPHOBIC: ("Hydrophobic", "type-hydrophobic"),
        ContactType.PI_STACKING: ("π-π Stacking", "type-pi"),
    }
    
    for ctype, count in contact_counts.items():
        label, css_class = type_map.get(ctype, (ctype.value, ""))
        cards.append(f"""
        <div class="summary-card {css_class}">
            <h3>{label}</h3>
            <div class="value">{count}</div>
        </div>
        """)
    
    total = sum(contact_counts.values())
    cards.insert(0, f"""
    <div class="summary-card">
        <h3>Total Contacts</h3>
        <div class="value">{total}</div>
    </div>
    """)
    
    return f'<div class="summary">{"".join(cards)}</div>'


def _generate_legend(contact_colors: Dict, color_scheme: str) -> str:
    """Generate legend HTML."""
    items = []
    
    type_labels = {
        ContactType.HYDROGEN_BOND: "Hydrogen Bond",
        ContactType.SALT_BRIDGE: "Salt Bridge",
        ContactType.DISULFIDE: "Disulfide Bond",
        ContactType.HYDROPHOBIC: "Hydrophobic Contact",
        ContactType.PI_STACKING: "π-π Stacking",
    }
    
    for ctype, color in contact_colors.items():
        label = type_labels.get(ctype, ctype.value)
        items.append(f"""
        <div class="legend-item">
            <div class="legend-box" style="background: {color}"></div>
            <span>{label}</span>
        </div>
        """)
    
    if color_scheme == "by_type":
        # Show example of multi-contact gradient
        items.append(f"""
        <div class="legend-item">
            <div class="legend-box" style="background: linear-gradient(135deg, {contact_colors[ContactType.HYDROGEN_BOND]} 0%, {contact_colors[ContactType.SALT_BRIDGE]} 50%, {contact_colors[ContactType.HYDROPHOBIC]} 100%)"></div>
            <span>Multiple Contact Types (gradient)</span>
        </div>
        """)
    
    items.append("""
    <div class="legend-item">
        <div class="legend-box" style="background: #f5f5f5; border: 2px solid #ddd"></div>
        <span>No Contacts</span>
    </div>
    """)
    
  
    
    items.append("""
    <div class="legend-item" style="margin-top: 10px; width: 100%; padding-top: 10px; border-top: 1px solid #ddd;">
        <span style="font-size: 12px; color: #666;">
            💡 Residues with multiple contact types show gradient backgrounds. 
            Colored dots at bottom indicate all contact types present.
        </span>
    </div>
    """)
    
    return f"""
    <div class="legend">
        <div class="legend-title">Legend - Color Scheme: {color_scheme.replace('_', ' ').title()}</div>
        <div class="legend-items">
            {"".join(items)}
        </div>
    </div>
    """


def _render_contact_chain(
    chain_id: str,
    residues: List,
    residue_summary: Dict,
    analysis: ContactAnalysis,
    interface_analysis: Optional[InterfaceAnalysis],
    color_scheme: str,
    contact_colors: Dict
) -> str:
    """Render a single chain with contact visualization."""
    
    # Get interface residue IDs if available
    interface_resids = set()
    if interface_analysis and chain_id in interface_analysis.chain_data:
        interface_resids = {
            rd.residue.id_str 
            for rd in interface_analysis.chain_data[chain_id] 
            if rd.is_interface
        }
    
    # Debug: Check what residue IDs look like
    # print(f"Chain {chain_id} residue summary keys: {list(residue_summary.keys())[:5]}")
    # if residues:
    #     print(f"First residue id_str: {residues[0].id_str}")
    
    # Count contacts in this chain
    total_contacts = 0
    interface_contact_count = 0
    
    cells = []
    for residue in residues:
        resid = residue.id_str
        is_interface = resid in interface_resids
        
        # Get contact information
        res_contacts = residue_summary.get(resid, {})
        contact_total = res_contacts.get('total', 0)
        total_contacts += contact_total
        
        if is_interface:
            interface_contact_count += contact_total
        
        # Determine cell color and styling
        contact_types_with_colors = []
        for ctype in ContactType:
            if res_contacts.get(ctype.value, 0) > 0:
                contact_types_with_colors.append((ctype, contact_colors.get(ctype, "#999")))
        
        num_contact_types = len(contact_types_with_colors)
        
        if color_scheme == "by_type" and num_contact_types > 0:
            if num_contact_types == 1:
                # Single contact type - solid color
                color = contact_types_with_colors[0][1]
                multi_class = ""
                style_vars = ""
            else:
                # Multiple contact types - gradient background
                color = "#ffffff"  # Base color (won't be used due to gradient)
                multi_class = f"multi-contact-{min(num_contact_types, 5)}"
                
                # Create CSS variables for gradient colors
                style_vars = " ".join([
                    f"--color{i+1}: {color};" 
                    for i, (_, color) in enumerate(contact_types_with_colors[:5])
                ])
        elif color_scheme == "by_count" and contact_total > 0:
            # Heatmap based on total contacts
            max_contacts = max((s.get('total', 0) for s in residue_summary.values()), default=1)
            color = _get_heatmap_color(contact_total, max_contacts)
            multi_class = ""
            style_vars = ""
        else:
            # No contacts - use light gray
            color = "#f5f5f5"
            multi_class = ""
            style_vars = ""
        
        # Text color for contrast
        if num_contact_types > 1:
            # For gradients, use dark text with white shadow for readability
            text_color = "#000"
            text_shadow = "text-shadow: 0 0 3px white, 0 0 3px white, 0 0 3px white;"
        else:
            text_color = "#000" if _is_light_color(color) else "#fff"
            text_shadow = ""
        
        # Get contact details for tooltip
        contact_details = []
        partners = set()
        for contact in analysis.contacts:
            if contact.residue1_id == resid:
                contact_details.append({
                    "type": contact.type.value,
                    "partner": contact.residue2_id,
                    "distance": contact.distance
                })
                partners.add(contact.residue2_id)
            elif contact.residue2_id == resid:
                contact_details.append({
                    "type": contact.type.value,
                    "partner": contact.residue1_id,
                    "distance": contact.distance
                })
                partners.add(contact.residue1_id)
        
        # Build contact type indicators (colored dots)
        contact_dots = []
        for ctype in ContactType:
            if res_contacts.get(ctype.value, 0) > 0:
                color_dot = contact_colors.get(ctype, "#999")
                contact_dots.append(f'<div class="contact-dot" style="background: {color_dot}"></div>')
        
        # Cell info for tooltip
        info = {
            "resname": residue.resname,
            "chain": chain_id,
            "resseq": f"{residue.resseq}{residue.icode}".strip(),
            "is_interface": is_interface,
            "total_contacts": contact_total,
            "contacts": contact_details,
            "partners": sorted(list(partners))
        }
        
        interface_class = "interface" if is_interface else ""
        has_contacts_class = "has-contacts" if contact_total > 0 else ""
        
        cells.append(f"""
        <div class="residue-cell {interface_class} {has_contacts_class} {multi_class}"
             style="background: {color}; color: {text_color}; {text_shadow} {style_vars}"
             data-info='{json.dumps(info)}'>
            <div class="resname">{residue.resname}</div>
            <div class="resseq">{info['resseq']}</div>
            {f'<div class="contact-badge">{contact_total}</div>' if contact_total > 0 else ''}
            {f'<div class="contact-indicators">{"".join(contact_dots[:5])}</div>' if contact_dots else ''}
        </div>
        """)
    
    # Chain statistics
    interface_badge = ""
    if interface_resids:
        interface_badge = f'<span class="badge interface">{len(interface_resids)} interface residues</span>'
    
    return f"""
    <div class="chain-section">
        <div class="chain-header">
            <h2>Chain {chain_id}</h2>
            <div class="chain-badges">
                {interface_badge}
                <span class="badge">{len(residues)} residues</span>
                <span class="badge">{total_contacts} contacts</span>
            </div>
        </div>
        <div class="residue-grid">
            {"".join(cells)}
        </div>
    </div>
    """


def _get_dominant_contact_color(res_contacts: Dict, contact_colors: Dict, mixed_threshold: int = 2) -> str:
    """
    Get color for the dominant contact type.
    If multiple contact types are present, blend colors or show as mixed.
    """
    if not res_contacts:
        return "#e8e8e8"  # Gray for no contacts
    
    # Count how many contact types are present
    contact_types_present = []
    for ctype in ContactType:
        count = res_contacts.get(ctype.value, 0)
        if count > 0:
            contact_types_present.append((ctype, count))
    
    if not contact_types_present:
        return "#e8e8e8"
    
    # If only one type, return its color
    if len(contact_types_present) == 1:
        return contact_colors.get(contact_types_present[0][0], "#999")
    
    # If multiple types and above threshold, show as "mixed" (gradient)
    if len(contact_types_present) >= mixed_threshold:
        # Return a purple-ish color to indicate mixed interactions
        return "#8e44ad"
    
    # Otherwise return the dominant (highest count) type
    dominant_type = max(contact_types_present, key=lambda x: x[1])[0]
    return contact_colors.get(dominant_type, "#999")


def _get_heatmap_color(count: int, max_count: int) -> str:
    """Get heatmap color based on contact count."""
    if count == 0:
        return "#f0f0f0"
    
    # Gradient from light yellow to dark red
    intensity = min(count / max_count, 1.0) if max_count > 0 else 0
    
    # Color interpolation: light -> dark
    r = int(255)
    g = int(255 * (1 - intensity * 0.8))
    b = int(255 * (1 - intensity))
    
    return f"#{r:02x}{g:02x}{b:02x}"


def _is_light_color(hex_color: str) -> bool:
    """Determine if a color is light (for contrast)."""
    hex_color = hex_color.lstrip('#')
    r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
    brightness = (r * 299 + g * 587 + b * 114) / 1000
    return brightness > 128


# ==================== Example Usage ====================

def example_usage():
    """Demonstrate contact visualization."""
    from ..pdbparser import parse_pdb
    from ..contacts import analyze_contacts
    from ..interface import compute_interface
    
    # Load structure
    cx = parse_pdb("complex.pdb")
    
    # Analyze contacts
    contact_analysis = analyze_contacts(cx, "H,L", "A")
    
    # Optionally analyze interface too
    interface_analysis = compute_interface(cx, ["H", "L"], ["A"], cutoff=5.0)
    
    # Generate visualization
    html = generate_contact_html(
        contact_analysis,
        cx,
        interface_analysis=interface_analysis,
        title="Antibody-Antigen Contact Map",
        color_scheme="by_type"  # or "by_count"
    )
    
    # Save to file
    with open("contact_visualization.html", "w") as f:
        f.write(html)
    
    print("Visualization saved to contact_visualization.html")
    
    