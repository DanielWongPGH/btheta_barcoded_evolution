import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np; rnd = np.random.default_rng()
import pickle

KELLY_COLORS = ['#FFB300', '#803E75', '#FF6800', '#A6BDD7',
                '#C10020', '#CEA262', '#817066', '#007D34',
                '#F6768E', '#00538A', '#FF7A5C', '#53377A',
                '#FF8E00', '#B32851', '#F4C800', '#7F180D',
                '#93AA00', '#593315', '#F13A13', '#232C16']

def as_si(x, ndp):
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    return r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))


def hamming(str1, str2):
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))

def plot_diagonal(ax, **kwargs):
    if 'color' not in kwargs: kwargs['color'] = 'black'
    if 'linestyle' not in kwargs: kwargs['linestyle'] = 'dashed'
    if 'zorder' not in kwargs: kwargs['zorder'] = 0

    min_xy = max( [ax.get_xlim()[0], ax.get_ylim()[0]])
    max_xy = min( [ax.get_xlim()[1], ax.get_ylim()[1]])
    ax.axline((min_xy, min_xy), (max_xy, max_xy), **kwargs)

def turn_off_ax(ax):
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

def set_log(ax):
    ax.set_yscale('log')
    ax.set_xscale('log')

def make_outer(fig, nrows, ncols, xlabel=None, ylabel=None):
    outer = mpl.gridspec.GridSpec(nrows=nrows, ncols=ncols, figure=fig)
    outer_ax = fig.add_subplot(outer[:, :])
    outer_ax.set_xlabel(xlabel)
    outer_ax.set_ylabel(ylabel)
    turn_off_ax(outer_ax)
    return outer

def make_fig_and_gridspec(n_panels, n_cols=5, size_per_panel=(3,3)):
    n_rows = (n_panels // n_cols) + (n_panels % n_cols)
    xdim = size_per_panel[0] * n_cols
    ydim = size_per_panel[1] * n_rows

    fig = plt.figure(figsize=(xdim, ydim))
    outer = mpl.gridspec.GridSpec(nrows=n_rows, ncols=n_cols, figure=fig)

    return fig, outer
