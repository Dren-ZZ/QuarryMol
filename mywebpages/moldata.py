import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
import dash
from dash import dcc, html
from dash_iconify import DashIconify

dash.register_page(__name__, path="/moldata", name="Molecular Data")
description = dmc.Container(
    [
        dmc.Title(
            "Molecule Search", order=1, fw=700, fz=60, ta="center", my=15,
        ),
        dmc.Text(
            "Enter a molecule worth exploring into the search bar to view its properties.", fz=20, ta="center", my=15
        )
    ]
)

searchbar = dmc.Autocomplete(
    id={"type": "suggested-list", "index": "autocomplete"},
    variant="filled",
    placeholder='(e.g. "Tyrosine")',
    w="450px", h="60p",
    ta="center",
    styles={
        "input": {"fontSize": 16, "paddingTop": 15, "paddingBottom": 15,},
    },
)

search_button = dmc.ActionIcon(
    DashIconify(icon="material-symbols:search-rounded"),
    id="_search-internal",
    size="lg",
    loaderProps={"type": "dots"},
    variant="light"
)

search = dmc.Flex(
    [searchbar, search_button],
    gap="xs", wrap="nowrap", justify="center",
    align="self-end",
    mb=20,
)

moldata = dmc.Box(
    [
        dcc.Store(id="_search-output-moldata"),
        dmc.Container(id="_search-output", fluid=True),
    ]
)

layout = dmc.Container(
    [
        description,
        search,
        moldata,
    ],
    fluid=True,
)

