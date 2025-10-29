"""Visualization module for structures and knot analysis."""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Arc
import numpy as np
import config
from pathlib import Path


class StructureVisualizer:
    """Visualize secondary structures and knot analysis results."""
    
    def __init__(self, output_dir: str = config.OUTPUT_DIR):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.colors = config.COLOR_PALETTE
        plt.style.use('seaborn-v0_8-darkgrid')
    
    def draw_dot_bracket(self, sequence: str, structure: str, 
                         title: str = "Secondary Structure") -> plt.Figure:
        """Draw secondary structure as dot-bracket diagram."""
        fig, ax = plt.subplots(figsize=(12, 6))
        
        positions = self._calculate_arc_positions(structure, len(sequence))
        
        # Draw bases
        for i, (base, pos) in enumerate(zip(sequence, positions)):
            color = self._get_base_color(structure[i])
            ax.plot(pos[0], pos[1], 'o', markersize=8, color=color, zorder=3)
            ax.text(pos[0], pos[1] - 0.15, base, ha='center', fontsize=9, weight='bold')
        
        # Draw base pairs
        pairs = self._get_base_pair_arcs(structure)
        for i, j in pairs:
            pos_i = positions[i]
            pos_j = positions[j]
            ax.plot([pos_i[0], pos_j[0]], [pos_i[1], pos_j[1]], 
                   'k-', linewidth=1.5, alpha=0.6, zorder=1)
        
        ax.set_xlim(-1, len(sequence) + 1)
        ax.set_ylim(-1.5, 2)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title(title, fontsize=14, weight='bold', pad=20)
        
        return fig
    
    def _calculate_arc_positions(self, structure: str, length: int) -> np.ndarray:
        """Calculate 2D positions for base arc representation."""
        positions = np.zeros((length, 2))
        
        for i in range(length):
            angle = np.pi * i / (length - 1) if length > 1 else 0
            positions[i, 0] = np.cos(angle)
            positions[i, 1] = np.sin(angle)
        
        return positions
    
    def _get_base_pair_arcs(self, structure: str) -> list:
        """Extract base pair positions from structure."""
        pairs = []
        stack = []
        for pos, char in enumerate(structure):
            if char == '(':
                stack.append(pos)
            elif char == ')' and stack:
                pairs.append((stack.pop(), pos))
        return pairs
    
    def _get_base_color(self, char: str) -> str:
        """Get color based on structure character."""
        if char == '(':
            return self.colors['stem']
        elif char == ')':
            return self.colors['stem']
        elif char == '.':
            return self.colors['loop']
        else:
            return self.colors['bulge']
    
    def plot_knot_risk_landscape(self, window_results: list, 
                                 seq_id: str = "sequence") -> plt.Figure:
        """Plot knot risk across genomic positions."""
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8))
        
        positions = [r['window_start'] for r in window_results]
        complexity = [r['complexity_score'] for r in window_results]
        crossings = [r['crossing_count'] for r in window_results]
        
        # Complexity score
        colors = [self.colors[self._risk_to_color(r['risk_level'])] 
                 for r in window_results]
        ax1.bar(positions, complexity, color=colors, alpha=0.7, width=20)
        ax1.axhline(y=config.KNOT_COMPLEXITY_THRESHOLD, color='red', 
                   linestyle='--', label='Risk Threshold')
        ax1.set_ylabel('Knot Complexity Score', fontsize=11, weight='bold')
        ax1.set_title(f'Knot Risk Landscape: {seq_id}', fontsize=13, weight='bold')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Crossing number
        ax2.plot(positions, crossings, 'o-', color=self.colors['knot'], 
                linewidth=2, markersize=6)
        ax2.axhline(y=config.CROSSING_NUMBER_THRESHOLD, color='red', 
                   linestyle='--', label='Crossing Threshold')
        ax2.set_xlabel('Genomic Position (bp)', fontsize=11, weight='bold')
        ax2.set_ylabel('Base Pair Crossings', fontsize=11, weight='bold')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        return fig
    
    def plot_risk_distribution(self, summary: dict) -> plt.Figure:
        """Plot risk level distribution."""
        fig, ax = plt.subplots(figsize=(10, 6))
        
        risks = summary['risk_distribution']
        categories = ['LOW', 'MEDIUM', 'HIGH', 'CRITICAL']
        values = [risks[cat] for cat in categories]
        colors_list = [self.colors.get(c.lower(), '#999999') for c in categories]
        
        bars = ax.bar(categories, values, color=colors_list, alpha=0.7, edgecolor='black')
        
        for bar, val in zip(bars, values):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(val)}', ha='center', va='bottom', weight='bold')
        
        ax.set_ylabel('Number of Windows', fontsize=11, weight='bold')
        ax.set_title('Risk Level Distribution', fontsize=13, weight='bold')
        ax.grid(axis='y', alpha=0.3)
        
        return fig
    
    def _risk_to_color(self, risk_level: str) -> str:
        """Map risk level to color key."""
        mapping = {'LOW': 'stem', 'MEDIUM': 'loop', 'HIGH': 'bulge', 'CRITICAL': 'knot'}
        return mapping.get(risk_level, 'background')
    
    def save_figure(self, fig: plt.Figure, filename: str):
        """Save figure to output directory."""
        filepath = self.output_dir / filename
        fig.savefig(str(filepath), dpi=config.FIGURE_DPI, bbox_inches='tight')
        print(f"Saved: {filepath}")
        plt.close(fig)


def main():
    """Example usage."""
    viz = StructureVisualizer()


if __name__ == '__main__':
    main()
