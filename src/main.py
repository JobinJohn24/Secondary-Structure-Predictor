"""Main orchestration module for DNA secondary structure prediction pipeline."""

import json
from pathlib import Path
from sequence_parser import SequenceParser
from structure_predictor import StructurePredictor
from knot_analyzer import KnotAnalyzer
from visualization import StructureVisualizer
from bioinformatics_analyzer import BioinformaticsAnalyzer
import config


class DNAStructurePipeline:
    """Complete pipeline orchestrating all analysis steps."""
    
    def __init__(self, fasta_path: str = None, output_dir: str = config.OUTPUT_DIR):
        # Auto-detect FASTA file if not provided
        if fasta_path is None:
            fasta_path = self._find_fasta_file()
        
        self.fasta_path = fasta_path
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        self.parser = SequenceParser()
        self.predictor = StructurePredictor()
        self.analyzer = KnotAnalyzer()
        self.visualizer = StructureVisualizer(output_dir)
        self.bio_analyzer = BioinformaticsAnalyzer(output_dir)
        
        self.results = {}
    
    def _find_fasta_file(self):
        """Auto-detect FASTA file in common locations."""
        search_paths = [
            Path('.'),
            Path('data'),
            Path('./'),
            config.INPUT_DIR
        ]
        
        fasta_extensions = ['.fasta', '.fa', '.fna']
        
        for search_dir in search_paths:
            if not search_dir.exists():
                continue
            for ext in fasta_extensions:
                fasta_files = list(search_dir.glob(f'*{ext}'))
                if fasta_files:
                    print(f"Found FASTA file: {fasta_files[0]}")
                    return str(fasta_files[0])
        
        raise FileNotFoundError(
            "No FASTA file. Place a .fasta, .fa, or .fna file in the directory."
        )
    
    def run_analysis(self) -> dict:
        """Execute complete analysis pipeline."""
        print(f"Starting DNA Structure Analysis Pipeline...")
        print(f"Processing: {self.fasta_path}\n")
        
        try:
            # Step 1: Parse sequences
            print("[1/4] Parsing FASTA file...")
            sequences = self.parser.parse_fasta(self.fasta_path)
            sequences = self.parser.filter_sequences(sequences)
            stats = self.parser.compute_stats(sequences)
            
            if not sequences:
                print("  ERROR: No valid sequences found!")
                return {}
            
            print(f"  ✓ Loaded {len(sequences)} valid sequences")
            
            # Step 2: Predict structures
            print("\n[2/4] Predicting secondary structures...")
            structure_data = {}
            window_predictions = []
            
            for seq_id, sequence in sequences.items():
                pred = self.predictor.predict_structure(sequence)
                pred['sequence_id'] = seq_id
                structure_data[seq_id] = pred
                
                windows = self.predictor.sliding_window_prediction(
                    sequence, stride=15
                )
                for w in windows:
                    w['sequence_id'] = seq_id
                window_predictions.extend(windows)
            
            print(f"  ✓ Predicted structures for {len(sequences)} sequences")
            
            # Step 3: Analyze knots
            print("\n[3/4] Analyzing knot topology...")
            knot_results = self.analyzer.analyze_window_knots(window_predictions)
            summary = self.analyzer.summarize_knot_regions(knot_results)
            
            print(f"  ✓ Risk distribution: {summary['risk_distribution']}")
            
            # Step 3.5: Bioinformatics analysis (pass window_predictions which has sequence/structure)
            bio_report = self.bio_analyzer.analyze_all(window_predictions)
            
            # Step 4: Visualize results
            print("\n[4/4] Generating visualizations...")
            for seq_id, pred in structure_data.items():
                fig = self.visualizer.draw_dot_bracket(
                    pred['sequence'], 
                    pred['structure'],
                    title=f"Structure: {seq_id}"
                )
                self.visualizer.save_figure(fig, f"{seq_id}_structure.png")
            
            # Plot landscape
            fig = self.visualizer.plot_knot_risk_landscape(
                knot_results,
                seq_id=list(sequences.keys())[0] if sequences else "all"
            )
            self.visualizer.save_figure(fig, "knot_risk_landscape.png")
            
            fig = self.visualizer.plot_risk_distribution(summary)
            self.visualizer.save_figure(fig, "risk_distribution.png")
            
            # Step 5: Save results
            self._save_results(stats, structure_data, knot_results, summary)
            
            print(f"\n✓ Analysis complete! Results saved to {self.output_dir}\n")
            return {'stats': stats, 'structures': structure_data, 'knots': summary}
        
        except Exception as e:
            print(f"ERROR during analysis: {e}")
            raise


    def _save_results(self, stats, structures, knots, summary):
        """Save analysis results to JSON files."""
        output_json = {
            'sequence_statistics': stats,
            'knot_summary': summary,
            'high_risk_regions': summary['high_risk_regions'][:10]
        }
        
        with open(self.output_dir / 'results.json', 'w') as f:
            json.dump(output_json, f, indent=2)
        
        with open(self.output_dir / 'knot_details.json', 'w') as f:
            json.dump(knots, f, indent=2, default=str)


def main():
    """Main entry point - runs automatically without arguments."""
    try:
        # Specify example_sequences.fasta or let it auto-detect
        pipeline = DNAStructurePipeline(fasta_path='example_sequences.fasta')
        results = pipeline.run_analysis()
    except FileNotFoundError as e:
        print(f"Setup Error: {e}")
    except Exception as e:
        print(f"Fatal Error: {e}")


if __name__ == '__main__':
    main()
