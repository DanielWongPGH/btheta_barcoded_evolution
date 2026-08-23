import sys
from pathlib import Path

project_root = Path(__file__).resolve().parents[1]
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from project_config import load_config

project_config = load_config()
local_paths = project_config['local']
remote_paths = project_config['remote']

project_root = local_paths['root']
notebooks_dir = local_paths['notebooks']
data_dir = local_paths['data']
pickled_dir = local_paths['pickles']
figures_dir = local_paths['figures']
figure_dir = figures_dir
tables_dir = local_paths['tables']
papers_dir = local_paths['papers']
provenance_dir = local_paths['provenance']
raw_processing_dir = local_paths['raw_processing']
plot_dir = local_paths['legacy_plots']
si_dir = local_paths['si_figures']

mpl_configs = {
    'lines.linewidth': 0.75,
    'font.family': 'Helvetica',
    'font.size': 8,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
    'savefig.dpi': 300,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'xtick.major.pad': 2,
    'xtick.minor.pad': 1.5,
    'ytick.major.pad': 2,
    'ytick.minor.pad': 1.5,
    'figure.facecolor': (1,1,1,1), #white
    'figure.edgecolor': (1,1,1,1), #white
    'figure.figsize': (6.5, 4),
    'axes.titlesize': 10,
    'legend.fontsize': 6,
    'axes.labelpad': 0}
