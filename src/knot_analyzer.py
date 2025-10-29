"""Knot detection and topological analysis module."""

import numpy as np
from typing import Dict, List, Tuple
import config


class KnotAnalyzer:
    """Identify and analyze knot-forming regions."""
    
    def __init__(self, writhe_threshold: float = config.WRITHE_THRESHOLD):
        self.writhe_threshold = writhe_threshold
        self.crossing_threshold = config.CROSSING_NUMBER_THRESHOLD
    
    def _get_base_pairs(self, structure: str) -> List[Tuple[int, int]]:
        """Extract base pair positions from dot-bracket notation."""
        pairs, stack = [], []
        for pos, char in enumerate(structure):
            if char == '(':
                stack.append(pos)
            elif char == ')' and stack:
                pairs.append((stack.pop(), pos))
        return sorted(pairs)
    
    def compute_writhe(self, structure: str) -> float:
        """Compute writhe as measure of topological twist."""
        writhe = 0.0
        stack = []
        for char in structure:
            if char == '(':
                stack.append(len(stack))
            elif char == ')' and stack:
                depth_diff = len(stack) - 1 - stack.pop()
                writhe += np.sign(depth_diff) if depth_diff != 0 else 0
        return writhe
    
    def identify_crossing_patterns(self, structure: str) -> List[Tuple[int, int]]:
        """Identify base pair crossings indicating knot-prone regions."""
        crossings, pairs = [], self._get_base_pairs(structure)
        for i, (p1_i, p1_j) in enumerate(pairs):
            for p2_i, p2_j in pairs[i + 1:]:
                if p1_i < p2_i < p1_j < p2_j or p2_i < p1_i < p2_j < p1_j:
                    crossings.append((min(p1_i, p2_i), max(p1_i, p2_i)))
        return crossings
    
    def _calc_nesting_complexity(self, structure: str) -> float:
        """Calculate complexity from nesting depth and transitions."""
        if not structure:
            return 0.0
        depth, max_depth, changes, prev = 0, 0, 0, 0
        for char in structure:
            depth += 1 if char == '(' else (-1 if char == ')' else 0)
            max_depth = max(max_depth, depth)
            changes += 1 if depth != prev else 0
            prev = depth
        nesting = min(max_depth / 10.0, 1.0)
        transitions = min(changes / len(structure), 1.0)
        return nesting * 0.6 + transitions * 0.4
    
    def detect_knots(self, structure: str) -> Dict:
        """Identify knot-prone regions with multiple criteria."""
        writhe = self.compute_writhe(structure)
        crossings = self.identify_crossing_patterns(structure)
        pairs = self._get_base_pairs(structure)
        nesting = max([0] + [max(0, max(0, *[1 if i < k < j else 0 for i, j in pairs])) for k in range(len(structure))]) if pairs else 0
        linking = nesting + writhe
        
        complexity = (
            min(abs(writhe) / 2.0, 1.0) * 0.3 +
            min(len(crossings) / 3.0, 1.0) * 0.3 +
            min(abs(linking) / 5.0, 1.0) * 0.2 +
            self._calc_nesting_complexity(structure) * 0.2
        )
        
        return {
            'writhe': writhe, 'crossing_count': len(crossings),
            'crossing_positions': crossings, 'linking_number': linking,
            'knot_prone': abs(writhe) >= self.writhe_threshold or len(crossings) >= self.crossing_threshold,
            'complexity_score': min(complexity, 1.0),
            'risk_level': ('LOW' if complexity < 0.1 else 'MEDIUM' if complexity < 0.2 else 'HIGH' if complexity < 0.3 else 'CRITICAL')
        }
    
    def analyze_window_knots(self, predictions: List[Dict]) -> List[Dict]:
        """Analyze knots across sliding window predictions."""
        results = []
        for pred in predictions:
            data = self.detect_knots(pred['structure'])
            data.update({'window_start': pred['window_start'], 'window_end': pred['window_end'], 'sequence_id': pred.get('sequence_id', 'unknown')})
            results.append(data)
        return results
    
    def summarize_knot_regions(self, window_results: List[Dict]) -> Dict:
        """Summarize high-risk knot regions across genome."""
        high_risk = [r for r in window_results if r['risk_level'] in ['HIGH', 'CRITICAL']]
        return {
            'total_windows_analyzed': len(window_results),
            'high_risk_windows': len(high_risk),
            'average_complexity': float(np.mean([r['complexity_score'] for r in window_results])),
            'max_complexity': float(max([r['complexity_score'] for r in window_results])),
            'high_risk_regions': high_risk,
            'risk_distribution': {
                'LOW': len([r for r in window_results if r['risk_level'] == 'LOW']),
                'MEDIUM': len([r for r in window_results if r['risk_level'] == 'MEDIUM']),
                'HIGH': len([r for r in window_results if r['risk_level'] == 'HIGH']),
                'CRITICAL': len([r for r in window_results if r['risk_level'] == 'CRITICAL'])
            }
        }
