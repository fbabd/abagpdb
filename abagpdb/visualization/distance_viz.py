from __future__ import annotations
from typing import List, Dict, Optional, Tuple, Union
import json
import math

try:
    import numpy as np
    _HAS_NUMPY = True
except:
    _HAS_NUMPY = False

from ..models import Complex
from ..distances import (
    atom_distance_matrix,
    residue_distance_matrix,
    pairwise_min_distance
)


def generate_distance_heatmap(
    distance_matrix: Union["np.ndarray", List[List[float]]],
    row_labels: List[str],
    col_labels: List[str],
    title: str = "Distance Matrix Heatmap",
    color_scheme: str = "blue_red",  # "blue_red", "viridis", "greyscale"
    cutoff: Optional[float] = None,
    show_values: bool = False,
    cluster: bool = False
) -> str:
    """
    Generate interactive heatmap visualization of distance matrix.
    
    Args:
        distance_matrix: 2D distance matrix
        row_labels: Labels for rows
        col_labels: Labels for columns
        title: Plot title
        color_scheme: Color gradient scheme
        cutoff: Highlight distances below cutoff
        show_values: Display numerical values on hover
        cluster: Perform hierarchical clustering
    
    Returns:
        HTML string with interactive heatmap
    """
    
    # Convert to numpy if needed
    if _HAS_NUMPY and not isinstance(distance_matrix, np.ndarray):
        distance_matrix = np.array(distance_matrix)
    
    # Calculate statistics
    if _HAS_NUMPY:
        flat = distance_matrix[~np.isnan(distance_matrix)].flatten()
        min_dist = float(np.min(flat)) if len(flat) > 0 else 0
        max_dist = float(np.max(flat)) if len(flat) > 0 else 10
        mean_dist = float(np.mean(flat)) if len(flat) > 0 else 5
    else:
        flat = [distance_matrix[i][j] for i in range(len(distance_matrix)) 
                for j in range(len(distance_matrix[0])) 
                if not math.isnan(distance_matrix[i][j])]
        min_dist = min(flat) if flat else 0
        max_dist = max(flat) if flat else 10
        mean_dist = sum(flat) / len(flat) if flat else 5
    
    # Color schemes
    color_schemes = {
        "blue_red": {
            "colors": ["#0d47a1", "#1976d2", "#42a5f5", "#90caf9", "#ffeb3b", "#ffa726", "#ff5722", "#b71c1c"],
            "description": "Blue (close) → Red (far)"
        },
        "viridis": {
            "colors": ["#440154", "#414487", "#2a788e", "#22a884", "#7ad151", "#fde724"],
            "description": "Purple (close) → Yellow (far)"
        },
        "greyscale": {
            "colors": ["#000000", "#404040", "#808080", "#b0b0b0", "#e0e0e0", "#ffffff"],
            "description": "Black (close) → White (far)"
        }
    }
    
    scheme = color_schemes.get(color_scheme, color_schemes["blue_red"])
    
    # Build data as JSON
    matrix_data = []
    for i in range(len(row_labels)):
        row = []
        for j in range(len(col_labels)):
            if _HAS_NUMPY:
                val = float(distance_matrix[i, j])
            else:
                val = float(distance_matrix[i][j])
            row.append(val if not math.isnan(val) else None)
        matrix_data.append(row)
    
    html = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>{title}</title>
    <style>
        * {{ box-sizing: border-box; }}
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
        }}
        .stats {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 15px;
            border-radius: 8px;
            color: white;
        }}
        .stat-label {{
            font-size: 12px;
            opacity: 0.9;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        .stat-value {{
            font-size: 24px;
            font-weight: bold;
            margin-top: 5px;
        }}
        .controls {{
            margin: 20px 0;
            padding: 15px;
            background: #f8f9fa;
            border-radius: 8px;
            display: flex;
            gap: 20px;
            flex-wrap: wrap;
            align-items: center;
        }}
        .control-group {{
            display: flex;
            align-items: center;
            gap: 10px;
        }}
        .control-group label {{
            font-size: 14px;
            font-weight: 500;
        }}
        .control-group input[type="range"] {{
            width: 200px;
        }}
        .control-group input[type="checkbox"] {{
            width: 18px;
            height: 18px;
        }}
        .heatmap-container {{
            overflow: auto;
            margin: 20px 0;
            border: 1px solid #ddd;
            border-radius: 8px;
            max-height: 800px;
        }}
        .heatmap {{
            display: grid;
            grid-template-columns: 100px repeat({len(col_labels)}, minmax(30px, 50px));
            gap: 1px;
            background: #ddd;
            padding: 1px;
            font-size: 10px;
        }}
        .heatmap-cell {{
            background: white;
            min-height: 30px;
            max-height: 50px;
            aspect-ratio: 1;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 9px;
            cursor: pointer;
            transition: all 0.2s;
            position: relative;
            overflow: hidden;
        }}
        .heatmap-cell:hover {{
            transform: scale(1.3);
            z-index: 10;
            box-shadow: 0 4px 12px rgba(0,0,0,0.2);
        }}
        .heatmap-header {{
            background: #34495e;
            color: white;
            font-weight: bold;
            position: sticky;
            top: 0;
            z-index: 5;
            padding: 3px;
            writing-mode: vertical-rl;
            text-align: center;
            font-size: 9px;
            min-height: 30px;
        }}
        .heatmap-row-label {{
            background: #ecf0f1;
            font-weight: bold;
            position: sticky;
            left: 0;
            z-index: 4;
            padding: 3px 5px;
            text-align: right;
            font-size: 9px;
            display: flex;
            align-items: center;
            justify-content: flex-end;
            min-height: 30px;
        }}
        .heatmap-corner {{
            background: #2c3e50;
            position: sticky;
            left: 0;
            top: 0;
            z-index: 10;
        }}
        
        /* Dynamic sizing controls */
        .size-controls {{
            margin: 15px 0;
            display: flex;
            gap: 15px;
            align-items: center;
            padding: 10px 15px;
            background: #e8f4f8;
            border-radius: 6px;
        }}
        .size-controls label {{
            font-size: 13px;
            font-weight: 500;
        }}
        .size-controls input[type="range"] {{
            width: 150px;
        }}
        .size-controls .size-value {{
            font-weight: bold;
            color: #2c3e50;
        }}
        .color-scale {{
            margin: 20px 0;
            padding: 20px;
            background: #f8f9fa;
            border-radius: 8px;
        }}
        .scale-gradient {{
            height: 30px;
            border-radius: 4px;
            margin: 10px 0;
            position: relative;
        }}
        .scale-labels {{
            display: flex;
            justify-content: space-between;
            font-size: 12px;
            color: #666;
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
            box-shadow: 0 4px 12px rgba(0,0,0,0.3);
        }}
        .tooltip.show {{ display: block; }}
        .highlight {{
            border: 3px solid #e74c3c !important;
            box-shadow: 0 0 10px rgba(231, 76, 60, 0.5) !important;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>{title}</h1>
        
        <div class="stats">
            <div class="stat-card">
                <div class="stat-label">Minimum Distance</div>
                <div class="stat-value">{min_dist:.2f} Å</div>
            </div>
            <div class="stat-card" style="background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);">
                <div class="stat-label">Mean Distance</div>
                <div class="stat-value">{mean_dist:.2f} Å</div>
            </div>
            <div class="stat-card" style="background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%);">
                <div class="stat-label">Maximum Distance</div>
                <div class="stat-value">{max_dist:.2f} Å</div>
            </div>
            <div class="stat-card" style="background: linear-gradient(135deg, #43e97b 0%, #38f9d7 100%);">
                <div class="stat-label">Matrix Size</div>
                <div class="stat-value">{len(row_labels)}×{len(col_labels)}</div>
            </div>
        </div>
        
        <div class="controls">
            <div class="control-group">
                <label for="cutoffSlider">Distance Cutoff:</label>
                <input type="range" id="cutoffSlider" min="{min_dist}" max="{max_dist}" 
                       value="{cutoff if cutoff else (min_dist + max_dist) / 2}" step="0.1">
                <span id="cutoffValue">{cutoff if cutoff else (min_dist + max_dist) / 2:.1f} Å</span>
            </div>
            <div class="control-group">
                <input type="checkbox" id="showValues" {"checked" if show_values else ""}>
                <label for="showValues">Show values</label>
            </div>
            <div class="control-group">
                <input type="checkbox" id="highlightCutoff" checked>
                <label for="highlightCutoff">Highlight below cutoff</label>
            </div>
        </div>
        
        <div class="size-controls">
            <label for="cellSizeSlider">Cell Size:</label>
            <input type="range" id="cellSizeSlider" min="10" max="80" value="40" step="5">
            <span class="size-value"><span id="cellSizeValue">40</span>px</span>
            <button class="control-btn" onclick="autoFitSize()" style="margin-left: 15px; padding: 6px 15px; font-size: 12px;">Auto Fit</button>
        </div>
        
        <div class="color-scale">
            <div style="font-weight: bold; margin-bottom: 10px;">Color Scale: {scheme['description']}</div>
            <div class="scale-gradient" id="colorGradient"></div>
            <div class="scale-labels">
                <span>{min_dist:.1f} Å</span>
                <span>{mean_dist:.1f} Å</span>
                <span>{max_dist:.1f} Å</span>
            </div>
        </div>
        
        <div class="heatmap-container">
            <div class="heatmap" id="heatmap"></div>
        </div>
        
        <div class="tooltip" id="tooltip"></div>
    </div>
    
    <script>
        const matrixData = {json.dumps(matrix_data)};
        const rowLabels = {json.dumps(row_labels)};
        const colLabels = {json.dumps(col_labels)};
        const minDist = {min_dist};
        const maxDist = {max_dist};
        const colors = {json.dumps(scheme['colors'])};
        
        let currentCutoff = {cutoff if cutoff else (min_dist + max_dist) / 2};
        let showValues = {json.dumps(show_values)};
        let highlightCutoff = true;
        let cellSize = 20;
        
        // Create gradient
        const gradient = document.getElementById('colorGradient');
        gradient.style.background = `linear-gradient(to right, ${{colors.join(', ')}})`;
        
        function getColor(distance) {{
            if (distance === null || isNaN(distance)) return '#f0f0f0';
            const normalized = (distance - minDist) / (maxDist - minDist);
            const idx = Math.floor(normalized * (colors.length - 1));
            return colors[Math.max(0, Math.min(colors.length - 1, idx))];
        }}
        
        function updateCellSize(size) {{
            cellSize = size;
            const heatmap = document.getElementById('heatmap');
            const labelWidth = Math.max(80, Math.min(150, cellSize * 4)); 
            heatmap.style.gridTemplateColumns = `${{labelWidth}}px repeat(${{colLabels.length}}, ${{cellSize}}px)`;
            
            // Update all cells
            document.querySelectorAll('.heatmap-cell').forEach(cell => {{
                cell.style.minHeight = cellSize + 'px';
                cell.style.maxHeight = cellSize + 'px';
                cell.style.fontSize = Math.max(6, Math.min(12, cellSize / 4)) + 'px';  
            }});
            
            document.querySelectorAll('.heatmap-header').forEach(cell => {{
                cell.style.minHeight = 2*cellSize + 'px';
                cell.style.fontSize = Math.max(8, Math.min(12, cellSize / 3)) + 'px';  
            }});
            
            document.querySelectorAll('.heatmap-row-label').forEach(cell => {{
                cell.style.minWidth = 7*cellSize + 'px';
                cell.style.fontSize = Math.max(8, Math.min(12, cellSize / 3)) + 'px'; 
            }}); 
            
        }}
        
        function autoFitSize() {{
            const container = document.querySelector('.heatmap-container');
            const availableWidth = container.clientWidth - 150; 
            const numCols = colLabels.length;
            const optimalSize = Math.floor(availableWidth / numCols);
            const size = Math.max(10, Math.min(10, optimalSize));  
            
            document.getElementById('cellSizeSlider').value = size;
            document.getElementById('cellSizeValue').textContent = size;
            updateCellSize(size);
        }}

        
        function renderHeatmap() {{
            const heatmap = document.getElementById('heatmap');
            heatmap.innerHTML = '';
            
            // Corner cell
            const corner = document.createElement('div');
            corner.className = 'heatmap-cell heatmap-corner';
            heatmap.appendChild(corner);
            
            // Column headers
            colLabels.forEach(label => {{
                const header = document.createElement('div');
                header.className = 'heatmap-cell heatmap-header';
                header.textContent = label;
                heatmap.appendChild(header);
            }});
            
            // Data rows
            matrixData.forEach((row, i) => {{
                // Row label
                const rowLabel = document.createElement('div');
                rowLabel.className = 'heatmap-cell heatmap-row-label';
                rowLabel.textContent = rowLabels[i];
                heatmap.appendChild(rowLabel);
                
                // Data cells
                row.forEach((distance, j) => {{
                    const cell = document.createElement('div');
                    cell.className = 'heatmap-cell';
                    
                    if (distance !== null) {{
                        cell.style.background = getColor(distance);
                        
                        if (showValues) {{
                            cell.textContent = distance.toFixed(1);
                        }}
                        
                        if (highlightCutoff && distance <= currentCutoff) {{
                            cell.classList.add('highlight');
                        }}
                        
                        cell.dataset.row = i;
                        cell.dataset.col = j;
                        cell.dataset.distance = distance;
                        
                        cell.addEventListener('mouseenter', showTooltip);
                        cell.addEventListener('mousemove', moveTooltip);
                        cell.addEventListener('mouseleave', hideTooltip);
                    }} else {{
                        cell.style.background = '#f0f0f0';
                        cell.textContent = '—';
                    }}
                    
                    heatmap.appendChild(cell);
                }});
            }});
        }}
        
        function showTooltip(e) {{
            const cell = e.target;
            const distance = cell.dataset.distance;
            if (!distance) return;
            
            const i = cell.dataset.row;
            const j = cell.dataset.col;
            
            const tooltip = document.getElementById('tooltip');
            tooltip.innerHTML = `
                <div><strong>${{rowLabels[i]}}</strong></div>
                <div><strong>${{colLabels[j]}}</strong></div>
                <div style="margin-top: 8px; padding-top: 8px; border-top: 1px solid rgba(255,255,255,0.3);">
                    Distance: <strong>${{parseFloat(distance).toFixed(2)}} Å</strong>
                </div>
                ${{parseFloat(distance) <= currentCutoff ? '<div style="color: #ff6b6b; margin-top: 5px;">✓ Within cutoff</div>' : ''}}
            `;
            tooltip.classList.add('show');
        }}
        
        function moveTooltip(e) {{
            const tooltip = document.getElementById('tooltip');
            tooltip.style.left = (e.clientX + 15) + 'px';
            tooltip.style.top = (e.clientY + 15) + 'px';
        }}
        
        function hideTooltip() {{
            document.getElementById('tooltip').classList.remove('show');
        }}
        
        // Controls
        document.getElementById('cutoffSlider').addEventListener('input', (e) => {{
            currentCutoff = parseFloat(e.target.value);
            document.getElementById('cutoffValue').textContent = currentCutoff.toFixed(1) + ' Å';
            renderHeatmap();
        }});
        
        document.getElementById('showValues').addEventListener('change', (e) => {{
            showValues = e.target.checked;
            renderHeatmap();
        }});
        
        document.getElementById('highlightCutoff').addEventListener('change', (e) => {{
            highlightCutoff = e.target.checked;
            renderHeatmap();
        }});
        
        document.getElementById('cellSizeSlider').addEventListener('input', (e) => {{
            const size = parseInt(e.target.value);
            document.getElementById('cellSizeValue').textContent = size;
            updateCellSize(size);
        }});
        
        // Initial render and auto-fit on load
        renderHeatmap();
        
        // Auto-fit on window resize
        let resizeTimer;
        window.addEventListener('resize', () => {{
            clearTimeout(resizeTimer);
            resizeTimer = setTimeout(() => {{
                if (document.getElementById('cellSizeSlider').value == 40) {{
                    autoFitSize();
                }}
            }}, 250);
        }});
    </script>
</body>
</html>
"""
    
    return html


def generate_contact_network(
    cx: Complex,
    selA: str,
    selB: str,
    cutoff: float = 5.0,
    level: str = "residue",
    title: str = "Contact Network"
) -> str:
    """
    Generate interactive network graph visualization of contacts.
    
    Args:
        cx: Complex structure
        selA: First selection (e.g., "H,L")
        selB: Second selection (e.g., "A")
        cutoff: Distance cutoff for contact
        level: "residue" or "atom"
        title: Plot title
    
    Returns:
        HTML with force-directed network graph
    """
    
    from ..distances import residue_distance_matrix, atom_distance_matrix
    
    # Get distance matrix
    if level == "residue":
        D, labelsA, labelsB = residue_distance_matrix(cx, selA, selB, mode="min")
    else:
        D, labelsA, labelsB = atom_distance_matrix(cx, selA, selB)
    
    # Build nodes and edges
    nodes = []
    edges = []
    
    # Nodes from group A
    for i, label in enumerate(labelsA):
        nodes.append({
            "id": f"A_{i}",
            "label": label,
            "group": "A",
            "index": i
        })
    
    # Nodes from group B
    for j, label in enumerate(labelsB):
        nodes.append({
            "id": f"B_{j}",
            "label": label,
            "group": "B",
            "index": j
        })
    
    # Edges (contacts within cutoff)
    contact_count = 0
    if _HAS_NUMPY:
        for i in range(len(labelsA)):
            for j in range(len(labelsB)):
                dist = float(D[i, j])
                if dist <= cutoff and not math.isnan(dist):
                    edges.append({
                        "source": f"A_{i}",
                        "target": f"B_{j}",
                        "distance": dist,
                        "weight": cutoff - dist  # Closer = stronger edge
                    })
                    contact_count += 1
    else:
        for i in range(len(labelsA)):
            for j in range(len(labelsB)):
                dist = D[i][j]
                if dist <= cutoff and not math.isnan(dist):
                    edges.append({
                        "source": f"A_{i}",
                        "target": f"B_{j}",
                        "distance": dist,
                        "weight": cutoff - dist
                    })
                    contact_count += 1
    
    html = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>{title}</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/d3/7.8.5/d3.min.js"></script>
    <style>
        * {{ box-sizing: border-box; }}
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
        }}
        .info {{
            background: #e8f4f8;
            padding: 15px;
            border-radius: 8px;
            margin: 20px 0;
        }}
        #network {{
            width: 100%;
            height: 700px;
            border: 1px solid #ddd;
            border-radius: 8px;
            background: #fafafa;
        }}
        .controls {{
            margin: 20px 0;
            display: flex;
            gap: 20px;
            flex-wrap: wrap;
        }}
        .control-btn {{
            padding: 10px 20px;
            border: none;
            border-radius: 6px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            cursor: pointer;
            font-weight: 500;
            transition: transform 0.2s;
        }}
        .control-btn:hover {{
            transform: translateY(-2px);
        }}
        .legend {{
            display: flex;
            gap: 30px;
            margin: 20px 0;
            padding: 15px;
            background: #f8f9fa;
            border-radius: 8px;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 10px;
        }}
        .legend-circle {{
            width: 20px;
            height: 20px;
            border-radius: 50%;
        }}
        .tooltip {{
            position: absolute;
            background: rgba(0,0,0,0.9);
            color: white;
            padding: 10px;
            border-radius: 6px;
            font-size: 12px;
            pointer-events: none;
            display: none;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>{title}</h1>
        
        <div class="info">
            <strong>Contact Network Analysis</strong><br>
            Found {contact_count} contacts within {cutoff} Å cutoff<br>
            Nodes: {len(nodes)} ({len(labelsA)} from group A, {len(labelsB)} from group B)<br>
            Click and drag nodes • Hover for details • Scroll to zoom
        </div>
        
        <div class="legend">
            <div class="legend-item">
                <div class="legend-circle" style="background: #3498db;"></div>
                <span>Group A ({selA})</span>
            </div>
            <div class="legend-item">
                <div class="legend-circle" style="background: #e74c3c;"></div>
                <span>Group B ({selB})</span>
            </div>
            <div class="legend-item">
                <div style="width: 30px; height: 2px; background: #95a5a6;"></div>
                <span>Contact (thickness = strength)</span>
            </div>
        </div>
        
        <div class="controls">
            <button class="control-btn" onclick="restartSimulation()">Restart Layout</button>
            <button class="control-btn" onclick="resetZoom()">Reset Zoom</button>
            <button class="control-btn" onclick="exportSVG()">Export SVG</button>
        </div>
        
        <div id="network"></div>
        <div class="tooltip" id="tooltip"></div>
    </div>
    
    <script>
        const nodes = {json.dumps(nodes)};
        const edges = {json.dumps(edges)};
        
        const width = document.getElementById('network').clientWidth;
        const height = 700;
        
        const svg = d3.select('#network')
            .append('svg')
            .attr('width', width)
            .attr('height', height);
        
        const g = svg.append('g');
        
        // Zoom behavior
        const zoom = d3.zoom()
            .scaleExtent([0.1, 10])
            .on('zoom', (event) => {{
                g.attr('transform', event.transform);
            }});
        
        svg.call(zoom);
        
        // Color scale
        const colorScale = {{
            'A': '#3498db',
            'B': '#e74c3c'
        }};
        
        // Force simulation
        const simulation = d3.forceSimulation(nodes)
            .force('link', d3.forceLink(edges).id(d => d.id).distance(d => d.distance * 10))
            .force('charge', d3.forceManyBody().strength(-300))
            .force('center', d3.forceCenter(width / 2, height / 2))
            .force('collision', d3.forceCollide().radius(25));
        
        // Draw edges
        const link = g.append('g')
            .selectAll('line')
            .data(edges)
            .join('line')
            .attr('stroke', '#95a5a6')
            .attr('stroke-width', d => Math.max(1, 5 * d.weight / {cutoff}))
            .attr('stroke-opacity', 0.6);
        
        // Draw nodes
        const node = g.append('g')
            .selectAll('g')
            .data(nodes)
            .join('g')
            .call(d3.drag()
                .on('start', dragStarted)
                .on('drag', dragged)
                .on('end', dragEnded));
        
        node.append('circle')
            .attr('r', 8)
            .attr('fill', d => colorScale[d.group])
            .attr('stroke', '#fff')
            .attr('stroke-width', 2);
        
        node.append('text')
            .text(d => d.label.split(' ')[0])  // Show just residue name
            .attr('x', 12)
            .attr('y', 4)
            .attr('font-size', '10px')
            .attr('fill', '#2c3e50');
        
        // Tooltip
        const tooltip = d3.select('#tooltip');
        
        node.on('mouseenter', (event, d) => {{
            const connections = edges.filter(e => e.source.id === d.id || e.target.id === d.id);
            tooltip.html(`
                <strong>${{d.label}}</strong><br>
                Group: ${{d.group}}<br>
                Connections: ${{connections.length}}
            `);
            tooltip.style('display', 'block');
        }})
        .on('mousemove', (event) => {{
            tooltip.style('left', (event.pageX + 10) + 'px')
                   .style('top', (event.pageY + 10) + 'px');
        }})
        .on('mouseleave', () => {{
            tooltip.style('display', 'none');
        }});
        
        link.on('mouseenter', (event, d) => {{
            tooltip.html(`
                <strong>${{d.source.label}}</strong> ↔ <strong>${{d.target.label}}</strong><br>
                Distance: ${{d.distance.toFixed(2)}} Å
            `);
            tooltip.style('display', 'block');
        }})
        .on('mousemove', (event) => {{
            tooltip.style('left', (event.pageX + 10) + 'px')
                   .style('top', (event.pageY + 10) + 'px');
        }})
        .on('mouseleave', () => {{
            tooltip.style('display', 'none');
        }});
        
        // Update positions
        simulation.on('tick', () => {{
            link
                .attr('x1', d => d.source.x)
                .attr('y1', d => d.source.y)
                .attr('x2', d => d.target.x)
                .attr('y2', d => d.target.y);
            
            node.attr('transform', d => `translate(${{d.x}},${{d.y}})`);
        }});
        
        function dragStarted(event, d) {{
            if (!event.active) simulation.alphaTarget(0.3).restart();
            d.fx = d.x;
            d.fy = d.y;
        }}
        
        function dragged(event, d) {{
            d.fx = event.x;
            d.fy = event.y;
        }}
        
        function dragEnded(event, d) {{
            if (!event.active) simulation.alphaTarget(0);
            d.fx = null;
            d.fy = null;
        }}
        
        function restartSimulation() {{
            simulation.alpha(1).restart();
        }}
        
        function resetZoom() {{
            svg.transition().duration(750).call(zoom.transform, d3.zoomIdentity);
        }}
        
        function exportSVG() {{
            const svgData = document.querySelector('#network svg').outerHTML;
            const blob = new Blob([svgData], {{type: 'image/svg+xml'}});
            const url = URL.createObjectURL(blob);
            const link = document.createElement('a');
            link.href = url;
            link.download = 'contact_network.svg';
            link.click();
        }}
    </script>
</body>
</html>
"""
    
    return html


def generate_distance_histogram(
    distance_matrix: Union["np.ndarray", List[List[float]]],
    title: str = "Distance Distribution",
    bins: int = 50,
    cutoffs: Optional[List[float]] = None
) -> str:
    """
    Generate interactive histogram of distance distributions.
    
    Args:
        distance_matrix: 2D distance matrix
        title: Plot title
        bins: Number of histogram bins
        cutoffs: List of distance cutoffs to highlight
    
    Returns:
        HTML with interactive histogram
    """
    cutoffs = cutoffs or [5.0, 8.0]
    
    # Flatten and remove NaNs
    if _HAS_NUMPY:
        flat = distance_matrix[~np.isnan(distance_matrix)].flatten()
        flat = flat[flat > 0]  # Remove self-distances
    else:
        flat = []
        for row in distance_matrix:
            for val in row:
                if not math.isnan(val) and val > 0:
                    flat.append(val)
    
    if len(flat) == 0:
        return "<html><body><h1>No valid distances to plot</h1></body></html>"
    
    # Calculate histogram
    min_val = min(flat)
    max_val = max(flat)
    bin_width = (max_val - min_val) / bins
    
    hist = [0] * bins
    for val in flat:
        bin_idx = min(int((val - min_val) / bin_width), bins - 1)
        hist[bin_idx] += 1
    
    # Bin edges
    bin_edges = [min_val + i * bin_width for i in range(bins + 1)]
    
    # Statistics
    mean_val = sum(flat) / len(flat)
    median_val = sorted(flat)[len(flat) // 2]
    
    # Count contacts within cutoffs
    cutoff_counts = {}
    for cutoff in cutoffs:
        count = sum(1 for v in flat if v <= cutoff)
        cutoff_counts[cutoff] = count
    
    html = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>{title}</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/4.4.0/chart.umd.min.js"></script>
    <style>
        * {{ box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f5f7fa;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 12px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.07);
        }}
        h1 {{
            color: #2c3e50;
            margin-bottom: 10px;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin: 25px 0;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
            border-radius: 10px;
            color: white;
        }}
        .stat-card:nth-child(2) {{
            background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
        }}
        .stat-card:nth-child(3) {{
            background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%);
        }}
        .stat-card h3 {{
            margin: 0 0 10px 0;
            font-size: 13px;
            opacity: 0.9;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        .stat-card .value {{
            font-size: 28px;
            font-weight: bold;
        }}
        .chart-container {{
            position: relative;
            height: 400px;
            margin: 30px 0;
        }}
        .cutoff-info {{
            margin: 25px 0;
            padding: 20px;
            background: #f8f9fa;
            border-radius: 8px;
        }}
        .cutoff-info h3 {{
            margin-top: 0;
            color: #2c3e50;
        }}
        .cutoff-list {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 15px;
            margin-top: 15px;
        }}
        .cutoff-item {{
            padding: 15px;
            background: white;
            border-left: 4px solid #3498db;
            border-radius: 4px;
        }}
        .cutoff-item .cutoff-value {{
            font-size: 20px;
            font-weight: bold;
            color: #3498db;
        }}
        .cutoff-item .cutoff-count {{
            font-size: 14px;
            color: #666;
            margin-top: 5px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>{title}</h1>
        
        <div class="stats-grid">
            <div class="stat-card">
                <h3>Total Distances</h3>
                <div class="value">{len(flat)}</div>
            </div>
            <div class="stat-card">
                <h3>Mean Distance</h3>
                <div class="value">{mean_val:.2f} Å</div>
            </div>
            <div class="stat-card">
                <h3>Median Distance</h3>
                <div class="value">{median_val:.2f} Å</div>
            </div>
        </div>
        
        <div class="chart-container">
            <canvas id="histogram"></canvas>
        </div>
        
        <div class="cutoff-info">
            <h3>Contacts Within Cutoff Distances</h3>
            <div class="cutoff-list">
                {"".join([f'''
                <div class="cutoff-item">
                    <div class="cutoff-value">≤ {cutoff} Å</div>
                    <div class="cutoff-count">{cutoff_counts[cutoff]} contacts ({100*cutoff_counts[cutoff]/len(flat):.1f}%)</div>
                </div>
                ''' for cutoff in cutoffs])}
            </div>
        </div>
    </div>
    
    <script>
        const ctx = document.getElementById('histogram').getContext('2d');
        
        const data = {{
            labels: {json.dumps([f"{bin_edges[i]:.1f}-{bin_edges[i+1]:.1f}" for i in range(bins)])},
            datasets: [{{
                label: 'Distance Count',
                data: {json.dumps(hist)},
                backgroundColor: 'rgba(52, 152, 219, 0.6)',
                borderColor: 'rgba(52, 152, 219, 1)',
                borderWidth: 1
            }}]
        }};
        
        const cutoffLines = {json.dumps(cutoffs)};
        
        const config = {{
            type: 'bar',
            data: data,
            options: {{
                responsive: true,
                maintainAspectRatio: false,
                plugins: {{
                    legend: {{
                        display: true,
                        position: 'top'
                    }},
                    title: {{
                        display: false
                    }},
                    tooltip: {{
                        callbacks: {{
                            label: function(context) {{
                                return 'Count: ' + context.parsed.y;
                            }}
                        }}
                    }},
                    annotation: {{
                        annotations: cutoffLines.reduce((acc, cutoff, idx) => {{
                            acc['cutoff' + idx] = {{
                                type: 'line',
                                xMin: cutoff,
                                xMax: cutoff,
                                borderColor: ['#e74c3c', '#f39c12', '#9b59b6'][idx % 3],
                                borderWidth: 2,
                                label: {{
                                    content: cutoff + ' Å cutoff',
                                    enabled: true,
                                    position: 'start'
                                }}
                            }};
                            return acc;
                        }}, {{}})
                    }}
                }},
                scales: {{
                    x: {{
                        title: {{
                            display: true,
                            text: 'Distance (Å)'
                        }},
                        ticks: {{
                            maxRotation: 45,
                            minRotation: 45
                        }}
                    }},
                    y: {{
                        title: {{
                            display: true,
                            text: 'Count'
                        }},
                        beginAtZero: true
                    }}
                }}
            }}
        }};
        
        new Chart(ctx, config);
    </script>
</body>
</html>
"""
    
    return html


# ===================================================================
# Unified Visualization Function
# ===================================================================

def visualize_distances(
    cx: Complex,
    selA: str,
    selB: Optional[str] = None,
    level: str = "residue",
    mode: str = "min",
    viz_type: str = "heatmap",
    output_file: str = "distance_viz.html",
    **kwargs
) -> str:
    """
    Unified function to visualize distances with multiple visualization types.
    
    Args:
        cx: Complex structure
        selA: First selection (e.g., "H,L")
        selB: Second selection (e.g., "A"), None for all-vs-all
        level: "residue" or "atom"
        mode: "min" or "CA" (for residue level)
        viz_type: Visualization type:
            - "heatmap": Interactive distance matrix heatmap
            - "network": Force-directed contact network graph
            - "histogram": Distance distribution histogram
            - "all": Generate all three and create index page
        output_file: Output HTML filename
        **kwargs: Additional parameters for specific visualizations
    
    Returns:
        Path to output HTML file
    
    Example:
        >>> # Heatmap of antibody-antigen distances
        >>> visualize_distances(cx, "H,L", "A", viz_type="heatmap", cutoff=5.0)
        
        >>> # Network graph of contacts
        >>> visualize_distances(cx, "H,L", "A", viz_type="network", cutoff=8.0)
        
        >>> # All visualizations
        >>> visualize_distances(cx, "H,L", "A", viz_type="all")
    """
    
    # Get distance matrix
    if level == "residue":
        D, labelsA, labelsB = residue_distance_matrix(
            cx, selA, selB, mode=mode
        )
    else:
        D, labelsA, labelsB = atom_distance_matrix(cx, selA, selB)
    
    # Generate visualization based on type
    if viz_type == "heatmap":
        html = generate_distance_heatmap(
            D, labelsA, labelsB,
            title=f"Distance Matrix: {selA} vs {selB or selA}",
            **kwargs
        )
        with open(output_file, "w") as f:
            f.write(html)
        return output_file
    
    elif viz_type == "network":
        if selB is None:
            raise ValueError("Network visualization requires two different selections")
        html = generate_contact_network(
            cx, selA, selB,
            cutoff=kwargs.get("cutoff", 5.0),
            level=level,
            title=f"Contact Network: {selA} ↔ {selB}"
        )
        with open(output_file, "w") as f:
            f.write(html)
        return output_file
    
    elif viz_type == "histogram":
        html = generate_distance_histogram(
            D,
            title=f"Distance Distribution: {selA} vs {selB or selA}",
            bins=kwargs.get("bins", 50),
            cutoffs=kwargs.get("cutoffs", [5.0, 8.0])
        )
        with open(output_file, "w") as f:
            f.write(html)
        return output_file
    
    elif viz_type == "all":
        # Generate all three visualizations
        base_name = output_file.replace(".html", "")
        
        # Heatmap
        heatmap_file = f"{base_name}_heatmap.html"
        heatmap_html = generate_distance_heatmap(
            D, labelsA, labelsB,
            title=f"Distance Matrix: {selA} vs {selB or selA}",
            **kwargs
        )
        with open(heatmap_file, "w") as f:
            f.write(heatmap_html)
        
        # Network (only if two selections)
        network_file = None
        if selB is not None:
            network_file = f"{base_name}_network.html"
            network_html = generate_contact_network(
                cx, selA, selB,
                cutoff=kwargs.get("cutoff", 5.0),
                level=level,
                title=f"Contact Network: {selA} ↔ {selB}"
            )
            with open(network_file, "w") as f:
                f.write(network_html)
        
        # Histogram
        histogram_file = f"{base_name}_histogram.html"
        histogram_html = generate_distance_histogram(
            D,
            title=f"Distance Distribution: {selA} vs {selB or selA}",
            bins=kwargs.get("bins", 50),
            cutoffs=kwargs.get("cutoffs", [5.0, 8.0])
        )
        with open(histogram_file, "w") as f:
            f.write(histogram_html)
        
        # Create index page
        index_html = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Distance Analysis Dashboard</title>
    <style>
        * {{ box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f5f7fa;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
        }}
        h1 {{
            color: #2c3e50;
            text-align: center;
            margin-bottom: 40px;
        }}
        .viz-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(350px, 1fr));
            gap: 25px;
        }}
        .viz-card {{
            background: white;
            padding: 30px;
            border-radius: 12px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.07);
            text-align: center;
            transition: transform 0.2s;
        }}
        .viz-card:hover {{
            transform: translateY(-5px);
            box-shadow: 0 8px 16px rgba(0,0,0,0.15);
        }}
        .viz-icon {{
            font-size: 60px;
            margin-bottom: 20px;
        }}
        .viz-card h2 {{
            color: #2c3e50;
            margin-bottom: 15px;
        }}
        .viz-card p {{
            color: #7f8c8d;
            margin-bottom: 25px;
            line-height: 1.6;
        }}
        .viz-btn {{
            display: inline-block;
            padding: 12px 30px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            text-decoration: none;
            border-radius: 6px;
            font-weight: 500;
            transition: transform 0.2s;
        }}
        .viz-btn:hover {{
            transform: scale(1.05);
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Distance Analysis Dashboard</h1>
        <p style="text-align: center; color: #7f8c8d; margin-bottom: 40px;">
            Analysis of {selA} vs {selB or selA} • Level: {level} • Mode: {mode}
        </p>
        
        <div class="viz-grid">
            <div class="viz-card">
                <div class="viz-icon">🔥</div>
                <h2>Distance Heatmap</h2>
                <p>Interactive color-coded matrix showing all pairwise distances. 
                   Hover for details, adjust cutoff, and identify close contacts.</p>
                <a href="{heatmap_file}" class="viz-btn">View Heatmap</a>
            </div>
            
            {f'''
            <div class="viz-card">
                <div class="viz-icon">🕸️</div>
                <h2>Contact Network</h2>
                <p>Force-directed graph visualization of contact network. 
                   Drag nodes, zoom, and explore interaction patterns.</p>
                <a href="{network_file}" class="viz-btn">View Network</a>
            </div>
            ''' if network_file else ''}
            
            <div class="viz-card">
                <div class="viz-icon">📊</div>
                <h2>Distance Histogram</h2>
                <p>Statistical distribution of distances showing frequency 
                   and contacts within common cutoff thresholds.</p>
                <a href="{histogram_file}" class="viz-btn">View Histogram</a>
            </div>
        </div>
    </div>
</body>
</html>
"""
        
        with open(output_file, "w") as f:
            f.write(index_html)
        
        return output_file
    
    else:
        raise ValueError(f"Unknown viz_type: {viz_type}. Use 'heatmap', 'network', 'histogram', or 'all'")


# ===================================================================
# Example Usage
# ===================================================================

def example_usage():
    """Demonstrate distance visualizations."""
    from ..pdbparser import parse_pdb
    
    cx = parse_pdb("complex.pdb")
    
    # 1. Interactive heatmap
    visualize_distances(
        cx, "H,L", "A",
        viz_type="heatmap",
        cutoff=5.0,
        color_scheme="blue_red",
        output_file="heatmap.html"
    )
    
    # 2. Contact network graph
    visualize_distances(
        cx, "H,L", "A",
        viz_type="network",
        cutoff=8.0,
        level="residue",
        output_file="network.html"
    )
    
    # 3. Distance histogram
    visualize_distances(
        cx, "H,L", "A",
        viz_type="histogram",
        bins=50,
        cutoffs=[4.0, 5.0, 8.0],
        output_file="histogram.html"
    )
    
    # 4. All visualizations with dashboard
    visualize_distances(
        cx, "H,L", "A",
        viz_type="all",
        cutoff=5.0,
        output_file="dashboard.html"
    )
    
