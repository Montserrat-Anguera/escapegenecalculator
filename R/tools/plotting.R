import::here(plotly, 'plot_ly', 'add_trace', 'layout')
import::here(magrittr, '%>%')

## Functions
## plot_gel
## plot_multiscatter


#' Plot Gel
#' 
#' @description
#' Insert a dataframe like this one:
#' +-----------+------+-------------+-----------+------------+-----------+
#' | gene_name | lane | exon_length | num_reads |    srpm    | intensity |
#' +-----------+------+-------------+-----------+------------+-----------+
#' | Tmsb4x    | mat  |     1371    |   25571   |  2878.255  |  1.000000 |
#' | Arhgap4   | mat  |     8105    |    4681   |   526.889  |  0.346293 |
#' | Xist      | mat  |    17946    |      58   |     6.476  |  0.201613 |
#' | Tmsb4x    | pat  |     1370    |      67   |     7.541  |  0.201909 |
#' | Arhgap4   | pat  |     8105    |      22   |     2.476  |  0.200501 |
#' | Xist      | pat  |    16791    |    21964  |  2472.253  |  0.887127 |
#' +-----------+------+-------------+-----------+------------+-----------+
#' 
#' @export
plot_gel <- function(
    df,
    lane='lane',
    genes='gene_name',
    size='size',
    intensity_col='intensity',
    text_col='intensity',
    column_order=c(),
    ylabel=NULL,
    title=NULL,
    yrange=c(NA, NA),
    showlegend=FALSE,
    font_size=12,
    tickangle=-45,
    band_width=0.04,
    lane_width=100,
    top_margin=150,
    height=800,
    min_width=200
) {

    num_lanes = nrow(unique(df[lane]))
    yrange = log(yrange, base=10)

    fig <- plot_ly(
        height =  height,
        width = min_width+lane_width*num_lanes,
        type = 'bar'
    )
    for (gene in unique(df[[genes]])) {
        subset = df[(df[genes] == gene), ]
        fig <- fig %>% add_trace(
            base=subset[[size]],
            x=subset[[lane]],
            y=subset[[size]]*band_width,
            text=subset[[text_col]],
            legendgroup=gene,
            alignmentgroup=TRUE,
            marker=list(
                color= c(paste0('rgba(0, 0, 0, ', subset[[intensity_col]]))
            ),
            name=gene,
            orientation='v',
            showlegend=TRUE,
            textposition="none",
            hovertemplate=paste0(
                'lane=%{x}<br>',
                'gene=', gene, '<br>',
                'size=%{base}<br>',
                text_col, '=%{text}<br>',
                '<extra></extra>'
            ),
            type='bar'
        )
    }

    fig <- fig %>% layout(
        title = list(
            text = title,
            x = 0
        ),
        xaxis = list(
            side = 'top',
            tickangle = tickangle,
            type = 'category',
            # autorange = 'reversed',
            categoryorder = 'array',
            categoryarray = column_order
        ),
        yaxis = list(
            showgrid = FALSE,
            zeroline = FALSE,
            type = 'log',
            range = yrange,
            ticks = 'outside',
            title_text = ylabel
        ),
        legend = list(
            title = list(text = 'genes'),
            tracegroupgap = 0
        ),
        margin = list(
            t = top_margin
        ),
        barmode = 'overlay',
        autosize = FALSE,
        plot_bgcolor = 'rgba(0,0,0,0)',
        showlegend =showlegend,
        font = list(
            size = font_size
        ),
        hovermode = 'closest'
    )

    return(fig)
}


#' Plot Multiple Scatter
#' 
plot_multiscatter <- function(
    df, x, y, color, size=NULL,
    xlabel=NULL, ylabel=NULL, clabel=NULL, title=NULL,
    xmin=NULL, xmax=NULL,
    ymin=NULL, ymax=NULL,
    hover_data=c(),
    xaxis_type='linear',
    yaxis_type='linear',
    color_discrete_map=NULL
) {

    # Range params
    if (is.null(xmin)) {
        xmin <- min(df[[x]])
    }
    if (is.null(xmax)) {
        xmax <- max(df[[x]])
    }
    if (xaxis_type=='linear') {
        xrange = c(xmin - (xmax-xmin) * 0.05, xmax + (xmax-xmin) * 0.05)  # add padding
    }
    if (xaxis_type=='log') {
        xrange = sapply(c(xmin, xmax), function(x) if (x > 0) {log(x, base=10)} else {NULL} )
    }

    if (is.null(ymin)) {
        ymin <- min(df[[y]])
    }
    if (is.null(ymax)) {
        ymax <- max(df[[y]])
    }
    if (yaxis_type=='linear') {
        yrange = c(ymin - (ymax-ymin) * 0.1, ymax + (ymax-ymin) * 0.1)  # add padding
    }
    if (yaxis_type=='log') {
        yrange = sapply(c(ymin, ymax), function(y) if (y > 0) {log(y, base=10)} else {NULL} )
    }

    # hoverdata
    hovertext <- ''
    for (field in c(color, x, y, size, hover_data)) {
        if (!is.null(field)) {
            hovertext <- paste0(hovertext, field, "=", df[[field]], "<br>")
        }
    }

    fig <- plot_ly() %>% add_trace(
        data = df,
        x = df[[x]],
        y = df[[y]],
        color = df[[color]],
        colors = color_discrete_map,
        size = df[[size]],
        fill = '',
        type = 'scatter',
        mode = 'markers',
        hoverinfo = 'text',
        hovertext = hovertext
    ) %>% layout(
        title = list(
            text = title,
            x = 0
        ),
        xaxis = list(
            title_text = xlabel,
            showgrid = TRUE, zeroline = FALSE,
            range = xrange,
            type = xaxis_type
        ),
        yaxis = list(
            title_text = ylabel,
            showgrid = TRUE, gridcolor = '#E4EAF2', zeroline = FALSE,
            range = yrange,
            type = yaxis_type
        ),
        plot_bgcolor = 'rgba(0,0,0,0)',
        showlegend = TRUE,
        hovermode = 'closest'
    )

    return(fig)
}
