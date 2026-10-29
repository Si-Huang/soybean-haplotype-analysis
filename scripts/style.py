"""Shared figure semantics; sizes refer to the final standalone figure."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

POPULATIONS = ['Soja', 'Landrace', 'Improved']
POP_LABELS = {'Soja': 'Wild (G. soja)', 'Landrace': 'Landrace', 'Improved': 'Improved'}
POP_COLORS = {'Soja': '#3182A0', 'Landrace': '#DEA126', 'Improved': '#CC5656'}
BASE_COLORS = {'A': '#73A65B', 'C': '#E3A06F', 'G': '#7A9FC4', 'T': '#397CA0'}
HAP_COLORS = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00',
              '#56B4E9', '#7B61A8', '#8C6D31', '#B85450', '#5C8585', '#525252', '#A6A33A']
RESPONSE_COLORS = {'R': '#0072B2', 'ER': '#2266AD', 'N': '#888888', 'S': '#D95F02'}


def setup():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 12,
                         'axes.labelsize': 14, 'axes.titlesize': 16,
                         'xtick.labelsize': 12, 'ytick.labelsize': 12,
                         'legend.fontsize': 12, 'pdf.fonttype': 42,
                         'svg.fonttype': 'none', 'axes.spines.top': False,
                         'axes.spines.right': False, 'figure.facecolor': 'white'})


def save(fig, directory, name):
    directory.mkdir(parents=True, exist_ok=True)
    for extension in ['pdf', 'svg', 'png']:
        fig.savefig(directory / f'{name}.{extension}', dpi=300,
                    bbox_inches='tight', facecolor='white')
    plt.close(fig)
