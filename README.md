# Secondary Structure Predictor: DNA Risk Classification via Knot Theory

A bioinformatics pipeline integrating knot theory with genomic analysis to classify synthetic DNA sequences by structural and biological risk. Uses five complementary metrics for comprehensive risk assessment in genetic engineering applications.

## Overview

This project analyzes 20 diverse DNA sequences (low to medium-high risk profiles) leveraging knot invariants and biophysical properties to predict secondary structure formation. By combining topological analysis with traditional bioinformatics, it provides novel risk stratification for synthetic biology.

## Five Core Metrics

| # | Test | Metric | Finding | Application |
|---|------|--------|---------|-------------|
| 1 | Test 1 | **GC Content** | Mean 53.39% | Thermal stability; 50-60% optimal for synthetic biology |
| 2 | Test 2 | **Melting Temperature** | Mean 65.4°C | DNA denaturation prediction; critical for PCR design |
| 3 | Test 3 | **Homopolymer Complexity** | 16.1% problematic | Detects risky repetitive sequences |
| 4 | Test 4 | **Shannon Entropy** | Mean 1.163 | Lower entropy indicates higher secondary structure risk |
| 5 | Test 5 | **Codon Usage Bias** | GAT/GAC dominant | Expression-level risks in synthetic constructs |

## Results & Visualizations

### Test 1: GC Content Distribution
![Test 1: GC Content](./results/test1_gc_content.png)

### Test 2: Melting Temperature Distribution
![Test 2: Melting Temperature](./results/test2_tm_analysis.png)

### Test 3: Homopolymer vs Complexity
![Test 3: Homopolymer](./results/test3_homopolymer.png)

### Test 4: Shannon Entropy Distribution
![Test 4: Entropy](./results/test4_entropy.png)

### Test 5: Codon Usage Bias
![Test 5: Codon Bias](./results/test5_codon_bias.png)

### Integrated Risk Assessment
![Risk Landscape](./results/knot_risk_landscape.png)

![Risk Distribution](./results/risk_distribution.png)

## Dataset Profile

- **Sequences Analyzed**: 20 diverse DNA sequences
- **Risk Distribution**: 35 low-risk, 27 medium-risk, 0 high-risk, 0 critical-risk
- **Structural Variety**: Simple repeats, GC-rich, hairpin-prone, inverted structures, triplex-forming, palindromic, complex topologies

## Methodology

### Knot Theory Integration
Knot invariants (crossing numbers, complexity metrics) characterize DNA topology. Higher topological complexity correlates with secondary structure formation potential and functional compromise.

### Five-Metric Approach
Stratified risk classification using: thermodynamic stability (GC, Tm), structural propensity (homopolymers, entropy), and expression risk (codon bias).

### Integration Strategy
Combines topological analysis with biophysical metrics to identify sequences at risk of unwanted secondary structures, off-target interactions, and synthesis failures.

## Key Findings

**Effective Stratification**: Five-metric approach successfully classified sequences across risk spectrum  
**Strong Correlation**: Topological features show high predictive value for structural stability  
**GC Sweet Spot**: 50–60% range optimal for synthetic biology applications  
**Entropy Matters**: Lower entropy sequences demonstrate significantly higher secondary structure propensity  
**Codon Considerations**: GAT/GAC dominance indicates expression-level design considerations  

## Quick Start

```bash
python main.py
```

## Technical Stack

- **Language**: Python 3.x
- **Libraries**: NumPy, Pandas, Matplotlib, Seaborn
- **Methods**: Knot invariant computation, thermodynamic modeling, statistical analysis

## Project Structure

```
secondary_structure_predictor/
├── config.py                      # Configuration & thresholds
├── knot_analyzer.py               # Knot theory calculations
├── bioinformatics_analyzer.py      # Metric computation
├── main.py                        # Pipeline orchestration
├── sequence_parser.py             # FASTA parsing
├── structure_predictor.py         # Secondary structure prediction
├── visualization.py               # Visualization generation
├── example_sequences.fasta         # Test data
├── requirements.txt               # Dependencies
├── README.md                      # This file
├── QUICKSTART.md                  # Quick reference
└── results/
    ├── test1_gc_content.png
    ├── test2_tm_analysis.png
    ├── test3_homopolymer.png
    ├── test4_entropy.png
    ├── test5_codon_bias.png
    ├── knot_risk_landscape.png
    ├── risk_distribution.png
    └── sequence_*.png             # 20 individual analyses
```
