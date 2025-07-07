import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
import dash
from dash import dcc, html
from dash_iconify import DashIconify

dash.register_page(__name__, path="/lab", name="Lab")
description = dmc.Container(
    [
        dmc.Title(
            "Laboratory Help", order=1, fw=700, fz=60, ta="center", my=15,
        ),
        dmc.Text(
            "Help with laboratory procedures. " \
            "Read the descriptions for more information.", fz=20, ta="center", my=15
        )
    ]
)

feature_card_1 = dmc.Card(
    p="lg",
    children=[
        dmc.Title(
            "Generate Reagent Tables", order=3,
            fz=30, fw=700, mb=15,
        ),
        dmc.Image(
            src="/assets/images/aspirin reagents.png",
            alt="Aspirin Reagents",
            w=300,
            h="auto",
            style={"display": "block", "margin": "auto"},
            radius="lg"
        ),
        dmc.List(
            [
                dmc.ListItem("Draft a table of reagents prior to laboratory work.", mb=10,),
                dmc.ListItem("Properties include melting and boiling points, densities, and first-aid measures.", mb=10,),
            ],
            fz=16, m="md"
        ),

        dmc.Center(
            dmc.Anchor(
                dmc.Button(
                    "Generate Reagent Table",
                    fz=18,
                    w="fit-content",
                    leftSection=DashIconify(
                        icon="bi:table",
                        width=25,
                    ),
                    variant="light",
                    size="md",
                    radius="xl"
                ),
                href="/lab/reagent-table",
            ),
        ),

    ],
    shadow="md",
    radius="xl",
    withBorder=True,
    w=400
)

row_cards = dmc.Flex(
    justify="center", gap="xl", wrap="nowrap", align="start",
    children=[feature_card_1,],
)


layout = dmc.Container(
    [
        description,
        dmc.Divider(mb=10, size="md"),
        row_cards,
    ],
    fluid=True,
)