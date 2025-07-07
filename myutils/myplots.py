import plotly.graph_objects as go

def create_pie_chart(labels, values):
    pie = go.Figure(data=[
        go.Pie(
            labels=labels, values=values,
            hoverinfo="label+value",
            textinfo="label+percent",
            hole=0.4,
            rotation=270,
        )
    ])

    pie.update_layout(
        autosize=False,
        width=400,
        height=400,
        title=dict(
            text="Proportion of Atoms",
            font=dict(size=20),
            xanchor="center",
            yanchor="middle",
            x=0.5,
            pad=dict(t=0, b=0, r=0, l=0),
        ),
        legend=dict(
            orientation="h",
            yanchor="middle",
            y=-0.1,
            xanchor="center",
            x=0.5,
        ),
        margin=dict(t=45, b=15, r=20, l=20),
    )

    return pie

def create_heatmap(
    z_vals, x_labels, y_labels, z_labels, hover, bar_title, map_title,
    color_light="#9ad2ff", color_dark="#007bff"
):
    color_scheme = [[0, color_light], [1, color_dark]]
    heatmap = go.Figure(
        data=go.Heatmap(
        z=z_vals,
        zmin=0,
        zmax=1,
        x=x_labels,
        y=y_labels,
        text=z_labels,
        texttemplate="%{text}",
        hovertemplate=hover,
        hoverongaps=False,
        colorscale=color_scheme,
        colorbar=dict(
            title=dict(text=bar_title, side="right"),
            thickness=20,
            y=0.5
        )),
    )

    heatmap.update_layout(
        title=dict(
            text=map_title,
            font=dict(size=20),
            x=0.5,
            xanchor="center", # Center at x
        ),
        width=600,
        height=600,
        xaxis=dict(
            title="Molecules X", tickangle=-45, constrain="domain", domain=[0,1],
            tickfont=dict(size=10), tickvals=x_labels,
            ticktext=[f"{short[:10]}..." if len(short)>10 else short for short in x_labels],
            tickson="boundaries"
        ),
        yaxis=dict(
            title="Molecules Y", scaleanchor="x", scaleratio=1,
            tickangle=-45, constrain="domain", domain=[0,1],
            tickfont=dict(size=10), tickvals=y_labels,
            ticktext=[f"{short[:10]}..." if len(short)>10 else short for short in y_labels],
            tickson="boundaries"
        ),
    )
    return heatmap
