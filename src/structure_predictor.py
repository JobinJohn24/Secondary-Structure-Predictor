"""Secondary structure prediction using thermodynamic models."""

import numpy as np
from typing import Dict, List, Tuple
import config


class StructurePredictor:
    """Predict DNA/RNA secondary structures using Nussinov algorithm."""
    
    def __init__(self, window_size: int = config.WINDOW_SIZE):
        self.window_size = window_size
        self.bp_complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A'}
        self.min_stem_length = config.MIN_STEM_LENGTH
    
    def can_pair(self, base1: str, base2: str) -> bool:
        """Check if two bases can form Watson-Crick base pair."""
        return base2 == self.bp_complement.get(base1, None)
    
    def nussinov_algorithm(self, sequence: str, min_loop: int = 1) -> np.ndarray:
        """
        Nussinov algorithm for secondary structure prediction.
        Returns score matrix where dp[i][j] = max base pairs in seq[i:j+1].
        """
        n = len(sequence)
        dp = np.zeros((n, n), dtype=float)
        
        for length in range(min_loop + 1, n + 1):
            for i in range(n - length + 1):
                j = i + length - 1
                
                # No pairing at position j
                dp[i][j] = dp[i][j - 1]
                
                # Try pairing i with k
                for k in range(i, j):
                    if self.can_pair(sequence[i], sequence[k]):
                        score = 1 + (dp[i + 1][k - 1] if i + 1 <= k - 1 else 0)
                        score += (dp[k + 1][j] if k + 1 <= j else 0)
                        dp[i][j] = max(dp[i][j], score)
        
        return dp
    
    def traceback(self, sequence: str, dp: np.ndarray, i: int, j: int) -> str:
        """Traceback to construct dot-bracket notation."""
        if i > j:
            return '.' * (j - i + 1)
        
        structure = ['.' for _ in range(len(sequence))]
        
        def traceback_recursive(i: int, j: int):
            if i > j:
                return
            
            if dp[i][j] == dp[i][j - 1]:
                traceback_recursive(i, j - 1)
            else:
                for k in range(i, j):
                    if self.can_pair(sequence[i], sequence[k]):
                        score = 1 + (dp[i + 1][k - 1] if i + 1 <= k - 1 else 0)
                        score += (dp[k + 1][j] if k + 1 <= j else 0)
                        
                        if dp[i][j] == score:
                            structure[i] = '('
                            structure[k] = ')'
                            traceback_recursive(i + 1, k - 1)
                            traceback_recursive(k + 1, j)
                            return
        
        traceback_recursive(i, j)
        return ''.join(structure)
    
    def predict_structure(self, sequence: str) -> Dict:
        """Predict secondary structure for a sequence."""
        dp = self.nussinov_algorithm(sequence)
        dot_bracket = self.traceback(sequence, dp, 0, len(sequence) - 1)
        
        return {
            'sequence': sequence,
            'structure': dot_bracket,
            'base_pairs': int(dp[0][len(sequence) - 1]),
            'score_matrix': dp,
            'length': len(sequence)
        }
    
    def sliding_window_prediction(self, sequence: str, stride: int = 10) -> List[Dict]:
        """Predict structures in sliding windows."""
        results = []
        
        for i in range(0, len(sequence) - self.window_size, stride):
            window = sequence[i:i + self.window_size]
            pred = self.predict_structure(window)
            pred['window_start'] = i
            pred['window_end'] = i + self.window_size
            results.append(pred)
        
        return results
    
    def extract_motifs(self, structure: str) -> Dict:
        """Extract secondary structure motifs."""
        motifs = {
            'stems': [],
            'loops': [],
            'bulges': [],
            'hairpins': []
        }
        
        i = 0
        while i < len(structure):
            if structure[i] == '(':
                stem_start = i
                while i < len(structure) and structure[i] == '(':
                    i += 1
                motifs['stems'].append((stem_start, i - 1))
            elif structure[i] == '.':
                loop_start = i
                while i < len(structure) and structure[i] == '.':
                    i += 1
                motifs['loops'].append((loop_start, i - 1))
            else:
                i += 1
        
        return motifs


def main():
    """Example usage."""
    predictor = StructurePredictor()
    # Example: predictor.predict_structure("GAATTCGC")


if __name__ == '__main__':
    main()
