"""Snippets of Python used to analyze the Berletch data
"""

import numpy as np
import pandas as pd
import plotly.graph_objs as go

# from harrison_functions.utils.plotting.plotly import save_fig_as_png

# Functions
# # read_many_csv
# # snake_to_title_case
# # adjust_closer
# # multi_melt
# # save_fig_as_png
# # plot_gel


def read_many_csv(filepaths):
    """Append only
    Make sure the columns are the same
    """
    
    dfs = []
    for filepath in filepaths:
        dfs.append(pd.read_csv(filepath))
        
    appended = pd.concat(dfs, axis=0)
        
    return appended


def snake_to_title_case(text):
    """Converts column_title to "Column Title"
    """
    return ' '.join(map(lambda x: x.capitalize(), text.split('_')))


def adjust_closer(comparator, base_value, tolerance=0.063):
    """If exon_length_pat is too far off from exon_length_mat, adjust it closer
    """
    if comparator / base_value <= 1-tolerance:
        return base_value * (1-tolerance)
    elif comparator / base_value >= 1+tolerance:
        return base_value * (1+tolerance)
    else:
        return comparator


def multi_melt(
        df,
        id_vars=['mouse_id', 'gene_name'],
        value_names=['num_reads', 'srpm'],
        value_vars=[
            # make sure these have suffixes
            ['num_reads_mat', 'num_reads_pat'],
            ['srpm_mat', 'srpm_pat'],
            # ['exon_length_mat', 'exon_length_pat'],
            # ['gene_id_mat', 'gene_id_pat'],
        ],
        var_name='variable'
    ):
    """Single-use data reshaper
    """

    dfs = []
    for value_var, value_name in zip(value_vars, value_names):
        tmp = pd.melt(
            df,
            id_vars=id_vars,
            value_vars=value_var,
            var_name=var_name,
        ).rename(columns={'value': value_name})
        
        tmp[var_name] = tmp[var_name].str.replace(f'{value_name}_', '')
        tmp.set_index(id_vars+[var_name], inplace=True)
        
        
        dfs.append(tmp)
        
    joined = pd.concat(dfs, axis=1).reset_index()
    
    return joined


def save_fig_as_png(fig, filepath, width=1200, height=800, scale=1, engine="kaleido"):
    """Make sure file extension is ".png"
    """
    if os.path.sep in filepath:
        os.makedirs(os.path.sep.join(str(filepath).split(os.path.sep )[:-1]), exist_ok=True)
    fig.write_image(filepath, width=width, height=height, scale=scale, engine=engine)


def plot_gel(
        df,
        lane='lane',
        genes='gene_name',
        size='size',
        intensity_col='intensity',
        text_col='intensity',
        column_order=[],
        ylabel=None,
        title=None,
        yrange=[None, None],  # only works if both are
        showlegend=False,
        font_size=12,
        tickangle=-45,
        band_width=0.04,
        lane_width=100,
        top_margin=150,
        height=800,
        min_width=200,
    ):
    """Insert a dataframe like this one:
    
    +-----------+------+-------------+-----------+------------+-----------+
    | gene_name | lane | exon_length | num_reads |    srpm    | intensity |
    +-----------+------+-------------+-----------+------------+-----------+
    | Tmsb4x    | mat  |     1371    |   25571   |  2878.255  |  1.000000 |
    | Arhgap4   | mat  |     8105    |    4681   |   526.889  |  0.346293 |
    | Xist      | mat  |    17946    |      58   |     6.476  |  0.201613 |
    | Tmsb4x    | pat  |     1370    |      67   |     7.541  |  0.201909 |
    | Arhgap4   | pat  |     8105    |      22   |     2.476  |  0.200501 |
    | Xist      | pat  |    16791    |    21964  |  2472.253  |  0.887127 |
    +-----------+------+-------------+-----------+------------+-----------+
    """
    
    df = df.sort_values([lane, size], ascending=[False, False])
    num_lanes = len(df[lane].unique())

    bars = []
    for gene in df[genes].unique():
        subset = df[(df[genes] == gene)]

        bar = go.Bar(
            base=subset[size],
            x=subset[lane],
            y=subset[size]*band_width,
            xaxis='x',
            yaxis='y',
            legendgroup=gene,
            alignmentgroup=True,
            marker={'color': [f'rgba(0, 0, 0, {intensity})' for intensity in subset[intensity_col]]},
            name=gene,
            orientation='v',
            showlegend=True,
            text=subset[text_col],
            textposition="none",
            hovertemplate='lane=%{x}<br>'
                          f'genes={gene}<br>'
                          'size=%{base}<br>'
                          +text_col+'=%{text}<br>'
                          '<extra></extra>',
        )
        bars.append(bar)

    fig = go.Figure()
    fig.add_traces(bars)
    
    
    yrange=[np.log10(lim) if lim else None for lim in yrange]
    fig.layout.update(
        {'xaxis': {'side': 'top', 'tickangle': tickangle,
                   'type': 'category', 'autorange': 'reversed',
                   'categoryorder': "array", 'categoryarray': column_order[::-1]
                  },
         'yaxis': {'type': 'log', 'range': yrange, 'ticks': 'outside', 'showline': True, 'title': ylabel},
         'title': title,
         'legend': {'title': {'text': 'genes'}, 'tracegroupgap': 0},
         'margin': {'t': top_margin},
         'barmode': 'overlay',
         'height': height,
         'width': min_width+lane_width*num_lanes,
         'autosize': False,
         'plot_bgcolor': 'rgba(0,0,0,0)',
         'showlegend': showlegend,
         'font': {'size': font_size},
        }
    )
    
    return fig
