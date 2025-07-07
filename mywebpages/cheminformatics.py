import dash_mantine_components as dmc
import dash
from dash import dcc, html
from dash_iconify import DashIconify

dash.register_page(__name__, path="/cheminformatics", name="Cheminformatics")
description = dmc.Container(
    [
        dmc.Title(
            "Cheminformatics Tools", order=1, fw=700, fz=60, ta="center", my=15,
        ),
        dmc.Text(
            "Explore any of the available cheminformatics tools. " \
            "Read the descriptions for more information.", fz=20, ta="center", my=15
        )
    ]
)

feature_card_1 = dmc.Card(
    p="lg",
    children=[
        dmc.Title(
            "Tanimoto Similarity/Coefficient", order=3,
            fz=30, fw=700, mb=15,
        ),

        dmc.Image(
            src="/assets/images/tanimoto-amino_acids.png",
            alt="Tanimoto similarity matrix",
            w=300,
            h=300,
            style={"display": "block", "margin": "auto"},
            radius="lg"
        ),

        dmc.List(
            [
                dmc.ListItem("Generate an interactive Tanimoto similarity matrix/heatmap for requested molecules.", mb=10,),
                dmc.ListItem("Click on each tile to view highlighted similarities between two specific molecules.", mb=10,),
                dmc.ListItem("Discover what the Tanimoto coefficient is and its applications.", mb=10,),
            ],
            fz=16, m="md"
        ),

        dmc.Center(
            dmc.Anchor(
                dmc.Button(
                    "Run Tanimoto Similarity",
                    fz=18,
                    w="fit-content",
                    leftSection=DashIconify(
                        icon="dinkie-icons:display-dot-matrix-small",
                        width=25,
                    ),
                    variant="light",
                    size="md",
                    radius="xl"
                ),
                href="/cheminformatics/tanimoto",
            ),
        ),

    ],
    shadow="md",
    radius="xl",
    withBorder=True,
    w=400,
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
