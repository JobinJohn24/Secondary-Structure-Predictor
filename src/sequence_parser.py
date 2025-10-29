"""Sequence parsing and validation module for FASTA files."""

import re
from typing import Dict, List, Tuple
from pathlib import Path
import config


class SequenceParser:
    """Parse and validate DNA/RNA sequences from FASTA format."""
    
    def __init__(self, min_length: int = config.MIN_SEQUENCE_LENGTH):
        self.min_length = min_length
        self.valid_nucleotides = config.VALID_NUCLEOTIDES
        self.stats = {}
    
    def parse_fasta(self, filepath: str) -> Dict[str, str]:
        """Parse FASTA file and return sequence dictionary."""
        sequences = {}
        current_id = None
        current_seq = []
        
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('>'):
                    if current_id:
                        sequences[current_id] = ''.join(current_seq)
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line.upper())
            
            if current_id:
                sequences[current_id] = ''.join(current_seq)
        
        return sequences
    
    def validate_sequence(self, seq_id: str, sequence: str) -> Tuple[bool, str]:
        """Validate sequence format and content."""
        if len(sequence) < self.min_length:
            return False, f"Sequence too short: {len(sequence)} < {self.min_length}"
        
        invalid_chars = set(sequence) - self.valid_nucleotides
        if invalid_chars:
            return False, f"Invalid nucleotides: {invalid_chars}"
        
        if sequence.count('N') / len(sequence) > 0.1:
            return False, "Too many unknown nucleotides (N > 10%)"
        
        return True, "Valid"
    
    def preprocess_sequence(self, sequence: str) -> str:
        """Clean and normalize sequence."""
        seq = sequence.upper()
        seq = re.sub(r'[^ATGCUN]', '', seq)
        return seq
    
    def compute_stats(self, sequences: Dict[str, str]) -> Dict:
        """Compute nucleotide composition and GC content."""
        stats = {}
        
        for seq_id, seq in sequences.items():
            gc_count = seq.count('G') + seq.count('C')
            at_count = seq.count('A') + seq.count('T')
            n_count = seq.count('N')
            
            stats[seq_id] = {
                'length': len(seq),
                'gc_content': gc_count / len(seq) if seq else 0,
                'at_content': at_count / len(seq) if seq else 0,
                'n_content': n_count / len(seq) if seq else 0,
                'unique_length': len(seq) - n_count
            }
        
        return stats
    
    def filter_sequences(self, sequences: Dict[str, str]) -> Dict[str, str]:
        """Filter sequences based on validation rules."""
        filtered = {}
        
        for seq_id, sequence in sequences.items():
            valid, msg = self.validate_sequence(seq_id, sequence)
            if valid:
                filtered[seq_id] = self.preprocess_sequence(sequence)
        
        return filtered


def main():
    """Example usage."""
    parser = SequenceParser()
    # Example: parser.parse_fasta('sequences.fasta')


if __name__ == '__main__':
    main()
