import dash
from dash import dcc
import dash_mantine_components as dmc
from dash_iconify import DashIconify
dash.register_page(__name__, path="/lab/reagent-table", name="Reagents")


reagent_instructions = dmc.Box(
    [
        dmc.Card(
            [
                dmc.Title("Instructions on Generating a Table of Reagents", order=3),
                dmc.List(
                    [
                        dmc.ListItem(dmc.Text("Search up a molecule.")),
                        dmc.ListItem(dmc.Text("Add the molecule-- molecule will appear in the space below.")),
                        dmc.ListItem(dmc.Text("Keep adding until desired amount reached.")),
                        dmc.ListItem(dmc.Text("Generate the table.")),
                        dmc.ListItem(dmc.Text("Copy into your lab notebook.")),
                    ],
                    type="ordered",
                    spacing="md"
                )
            ],
            withBorder=True,
            shadow="lg",
            p="md",
            radius="xl",
            mb=15,
            w="fit-content",
        )
    ],
    style={"display": "flex", "justifyContent": "center"},
    
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
    mb=15,
)

add_button = dmc.Button(
    "Add Molecule",
    leftSection=DashIconify(icon="ci:add-to-queue", width=25),
    id="_add-molecule",
    variant="light",
    w="fit-content",
    fz=18,
    size="md",
    radius="xl",
    loaderProps={"type": "dots"}
)

run_button = dmc.Button(
    "Generate Table",
    disabled=True,
    color="green",
    leftSection=DashIconify(icon="fluent-mdl2:generate", width=25),
    id={"type": "_gen-graphic", "index": "table"},
    variant="light",
    w="fit-content",
    fz=18,
    size="md",
    radius="xl",
    loaderProps={"type": "dots"}
)


buttons_search = dmc.Flex(
    [add_button, run_button],
    gap="xs", wrap="nowrap", justify="center",
    align="center",
)

search_group = dmc.Stack([buttons_search, searchbar], align="center")

mol_list = dmc.Box(
    children=[
        dmc.Text(id="_num-molecules", p="xs", fs="italic"),
        dmc.ScrollArea(
            h=400,
            id="_list-molecules",
        )
    ],
    # withBorder=True,
    # shadow="lg",
    # p="md",
    # radius="xl",
    # mb=15,
)


return_alert = dmc.Container(id="_feedback-alert", fluid=True)

reagent_data = dmc.Box(
    [   
        dcc.Store(id="_store-items", data=[]),
        return_alert,
        mol_list,
    ]
)

table = dmc.Container(
    id="_reagents-table-container",
    style={"display": "flex", "justifyContent": "center"},
    fluid=True
)


layout = dmc.Container(
    [
        dmc.Title(
            "Generate a Table of Reagents", order=1, fw=700, fz=60, ta="center", my=15,
        ),
        dmc.Text(
            "Draft a reagent table with necessary properties, including " \
            "melting point, boiling point, density, first aid, etc.", fz=20, ta="center", my=15
        ),
        dmc.Divider(mb=15, size="md"),
        reagent_instructions,
        dmc.Divider(mb=15, size="md"),
        search_group,
        dmc.Divider(mb=15, size="md"),
        dcc.Store(id="_reagents-table-data"),
        reagent_data,
        dmc.Divider(mb=15, size="md"),
        table
    ],
    fluid=True,
)

