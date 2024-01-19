import::here(plotly, 'plot_ly', 'layout')

## Functions
## plot_multiscatter


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

    fig <- plot_ly(
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
        title = title,
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
