import dash_mantine_components as dmc
import dash
from dash import dcc, html
from dash_iconify import DashIconify

dash.register_page(__name__, path="/updates", name="Updates")

description = dmc.Container(
    [
        dmc.Title(
            "Website Updates", order=1, fw=700, fz=60, ta="center", my=15,
        ),
        dmc.Text(
            "Major changes/updates to this website are documented here.",
            fz=20, ta="center", my=15
        )
    ]
)

def timeline_item(date, text, size, color, icon=None):
    return dmc.TimelineItem(
        title=date,
        bullet=DashIconify(icon=icon) if icon else None,
        children=[
            dmc.Text(text, size=size, c=color)
        ]
    )

timeline = dmc.Timeline(
    active=0,
    bulletSize=30,
    lineWidth=4,
    children=[
        timeline_item(
            "Tue Jul 8 2025",
            ["Some changes..."],
            size="md",
            color="dimmed"
        ),
        timeline_item(
            "Mon Jul 7 2025",
            ["Website deployed."],
            size="md",
            color="dimmed",
            icon="octicon:rocket-16"
        ),      
    ],
    fz=35,
)

layout = dmc.Container(
    [
        description,
        dmc.Divider(mb=20, size="md"),
        dmc.Center(timeline)
    ],
    fluid=True
)