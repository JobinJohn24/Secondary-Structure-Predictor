"""Bioinformatics Analysis - 5 Key Tests"""

import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, List
import seaborn as sns


class BioinformaticsAnalyzer:
    """5 key bioinformatics analyses for DNA structures."""
    
    def __init__(self, output_dir: str):
        self.output_dir = output_dir
        sns.set_style("whitegrid")
    
    def _metrics(self, kr: List[Dict]):
        """Calculate metrics with error handling."""
        gc, tm, hp, ent, cplx = [], [], [], [], []
        for r in kr:
            seq = r.get('sequence', '')
            st = r.get('structure', '')
            if seq and len(seq) > 0:
                gc_v = (seq.count('G') + seq.count('C')) / len(seq) * 100
                gc.append(gc_v)
                if len(seq) >= 14: 
                    tm.append(64 + 41 * (gc_v/100 - 0.5))
                mx, c = 1, 1
                for i in range(1, len(seq)):
                    c = c + 1 if seq[i] == seq[i-1] else 1
                    mx = max(mx, c)
                hp.append(mx)
            if st and len(st) > 0:
                ch = {}
                for x in st:
                    ch[x] = ch.get(x, 0) + 1
                e = sum(-c/len(st) * np.log2(c/len(st)) for c in ch.values())
                ent.append(e)
            cplx.append(r.get('complexity_score', 0))
        return gc, tm, hp, ent, cplx
    
    def analyze_all(self, kr: List[Dict]) -> Dict:
        """Generate all 5 tests."""
        print("\n[BioInfo] Running 5 bioinformatics tests...")
        gc, tm, hp, ent, cplx = self._metrics(kr)
        
        r = {'gc_m': np.mean(gc) if gc else 0, 'gc_s': np.std(gc) if gc else 0,
             'tm_m': np.mean(tm) if tm else 0, 'tm_s': np.std(tm) if tm else 0,
             'hp_m': np.mean(hp) if hp else 0, 'hp_p': (sum(1 for x in hp if x >= 5) / len(hp) * 100) if hp else 0,
             'e_m': np.mean(ent) if ent else 0, 'e_s': np.std(ent) if ent else 0}
        
        self._p1(gc, cplx)
        self._p2(tm)
        self._p3(hp, cplx)
        self._p4(ent)
        self._p5(kr)
        self._psummary(r)
        
        print(f"  ✓ Test 1: GC Content (Mean: {r['gc_m']:.2f}%)")
        print(f"  ✓ Test 2: Melting Temp (Mean: {r['tm_m']:.1f}°C)")
        print(f"  ✓ Test 3: Homopolymers ({r['hp_p']:.1f}% problematic)")
        print(f"  ✓ Test 4: Entropy (Mean: {r['e_m']:.3f})")
        print(f"  ✓ Test 5: Codon Bias")
        return r
    
    def _p1(self, gc: List, cplx: List):
        """Test 1: GC vs Complexity."""
        if len(gc) < 2: return
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.scatter(gc, cplx, alpha=0.6, s=100, c='#e74c3c', edgecolors='black')
        z = np.polyfit(gc, cplx, 1)
        ax.plot(sorted(gc), np.poly1d(z)(sorted(gc)), "r--", linewidth=2)
        ax.set_title('Test 1: GC Content vs Complexity', fontweight='bold', fontsize=12)
        ax.set_xlabel('GC Content (%)', fontweight='bold')
        ax.set_ylabel('Complexity Score', fontweight='bold')
        ax.grid(alpha=0.3)
        self._save(fig, "test1_gc_content.png")
    
    def _p2(self, tm: List):
        """Test 2: Melting Temp."""
        if not tm: return
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(tm, bins=12, color='#e67e22', edgecolor='black', alpha=0.7)
        ax.axvline(np.mean(tm), color='red', linestyle='--', linewidth=2, label=f'Mean: {np.mean(tm):.1f}°C')
        ax.set_title('Test 2: Melting Temperature Distribution', fontweight='bold', fontsize=12)
        ax.set_xlabel('Tm (°C)', fontweight='bold')
        ax.set_ylabel('Frequency', fontweight='bold')
        ax.legend()
        self._save(fig, "test2_tm_analysis.png")
    
    def _p3(self, hp: List, cplx: List):
        """Test 3: Homopolymer."""
        if not hp or not cplx: return
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.scatter(hp, cplx, alpha=0.6, s=100, c='#9b59b6', edgecolors='black')
        ax.axvline(5, color='red', linestyle='--', linewidth=2, label='Threshold (5bp)')
        ax.set_title('Test 3: Homopolymer vs Complexity', fontweight='bold', fontsize=12)
        ax.set_xlabel('Max Run Length (bp)', fontweight='bold')
        ax.set_ylabel('Complexity Score', fontweight='bold')
        ax.legend()
        ax.grid(alpha=0.3)
        self._save(fig, "test3_homopolymer.png")
    
    def _p4(self, ent: List):
        """Test 4: Entropy."""
        if not ent: return
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(ent, bins=15, color='#16a085', edgecolor='black', alpha=0.7)
        ax.axvline(np.mean(ent), color='red', linestyle='--', linewidth=2, label=f'Mean: {np.mean(ent):.3f}')
        ax.set_title('Test 4: Shannon Entropy Distribution', fontweight='bold', fontsize=12)
        ax.set_xlabel('Shannon Entropy', fontweight='bold')
        ax.set_ylabel('Frequency', fontweight='bold')
        ax.legend()
        self._save(fig, "test4_entropy.png")
    
    def _p5(self, kr: List):
        """Test 5: Codon Bias."""
        fig, ax = plt.subplots(figsize=(10, 6))
        cd = {}
        for r in kr:
            for i in range(0, len(r.get('sequence', '')) - 2, 3):
                c = r['sequence'][i:i+3]
                if len(c) == 3: cd[c] = cd.get(c, 0) + 1
        top = sorted(cd.items(), key=lambda x: x[1], reverse=True)[:12]
        if top:
            codons, counts = zip(*top)
            bars = ax.barh(codons, counts, color='#3498db', edgecolor='black')
            for b, cnt in zip(bars, counts):
                ax.text(cnt, b.get_y() + b.get_height()/2, f' {cnt}', va='center', fontweight='bold')
        ax.set_xlabel('Frequency', fontweight='bold')
        ax.set_title('Test 5: Top Codon Usage', fontweight='bold', fontsize=12)
        self._save(fig, "test5_codon_bias.png")
    
    def _psummary(self, r: Dict):
        """Summary report."""
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.axis('off')
        text = f"""BIOINFORMATICS SUMMARY
{'='*40}
Test 1 - GC Content:
  Mean: {r['gc_m']:.2f}% | Std: {r['gc_s']:.2f}%

Test 2 - Melting Temp:
  Mean: {r['tm_m']:.1f}°C | Std: {r['tm_s']:.1f}°C

Test 3 - Homopolymer:
  Mean: {r['hp_m']:.2f}bp | Problem: {r['hp_p']:.1f}%

Test 4 - Entropy:
  Mean: {r['e_m']:.3f} | Std: {r['e_s']:.3f}

Test 5 - Codon Bias: Complete"""
        ax.text(0.05, 0.95, text, fontfamily='monospace', fontsize=11, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        self._save(fig, "bioinformatics_summary.png")
    
    def _save(self, fig, name: str):
        filepath = f"{self.output_dir}/{name}"
        fig.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"    Saved: {name}")
        plt.close(fig)
