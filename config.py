"""Configuration and constants for DNA secondary structure prediction."""

from pathlib import Path
# Sequence parameters
MIN_SEQUENCE_LENGTH = 20
MAX_SEQUENCE_LENGTH = 5000
VALID_NUCLEOTIDES = {'A', 'T', 'G', 'C', 'U', 'N'}

# Structure prediction parameters
WINDOW_SIZE = 30
ENERGY_THRESHOLD = -5.0
MIN_STEM_LENGTH = 4
MAX_LOOP_SIZE = 30

# Knot detection parameters
KNOT_COMPLEXITY_THRESHOLD = 0.3
WRITHE_THRESHOLD = 8.0
CROSSING_NUMBER_THRESHOLD = 2

# Visualization parameters
FIGURE_DPI = 300
FIGURE_SIZE = (14, 10)
COLOR_PALETTE = {
    'stem': '#2E86AB',
    'loop': '#A23B72',
    'bulge': '#F18F01',
    'knot': '#C73E1D',
    'background': '#F6F6F6'
}

# File paths
OUTPUT_DIR = './results'
INPUT_DIR = Path('.')
CACHE_DIR = './.cache'

# Thermodynamic model
STRUCTURE_MODEL = 'RNA37'
