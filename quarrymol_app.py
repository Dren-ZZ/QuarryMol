from dash import Dash
import dash
import dash_mantine_components as dmc
from dash_iconify import DashIconify

import myactions.callbacks as callbacks

dmc.add_figure_templates("mantine_light")
app = Dash(
    __name__, external_stylesheets=dmc.styles.ALL,
    use_pages=True, pages_folder="mywebpages",
    suppress_callback_exceptions=True,
    meta_tags=[
        {
            "name": "viewport",
            "content": "width=device-width, initial-scale=1.0, maximum-scale=1.2, minimum-scale=0.5,"
        }
    ],
)
server = app.server


FONT_SIZE_SM = 16

def nav_link(label, href,):
    return dmc.Anchor(label, href=href, fz=20, fw=700)

def menu_item(label, href,):
    return dmc.MenuItem(label, href=href, fz=20)

navbar_dropdown_menu = dmc.Menu(
    trigger="hover",
    openDelay=100,
    closeDelay=200,
    withArrow=True,
    position="bottom",
    zIndex=1000,
    children=[
        dmc.MenuTarget(
            dmc.Button(
                children="Tools", fz=20, fw=700,
                variant="light",
            )
        ),
        dmc.MenuDropdown([
            menu_item("Molecular Data", href="/moldata"),
            menu_item("Cheminformatics", href="/cheminformatics"),
            menu_item("Lab", href="/lab"),
        ])
    ],
)

navbar = dmc.Paper(
    p="md",
    shadow="sm",
    opacity=0.95,
    style={
        "position": "sticky",
        "top": 0,
        "zIndex": 1,
    },
    children=dmc.Flex(
        justify="space-between",  # left/right alignment
        children=[
            dmc.Anchor(
                dmc.Button(
                    dmc.Title(
                        "QuarryMol", order=1, 
                        fw=700,
                    ),
                    variant="subtle",
                    size="lg",
                ),
                href="/"
            ),
            dmc.Group(
                wrap="wrap",
                children=[
                    nav_link("Home", href="/",),
                    nav_link("Updates", href="/updates"),
                    navbar_dropdown_menu,
                ]
            )
        ]
    ),
)

footer = dmc.Paper(
    p="sm",
    style={
        "textAlign": "center",
        "marginTop": "auto",
    },
    children=[
        dmc.Text(
            "QuarryMol • 2025 • Darren Liu • Developed with Dash, Plotly, RDKit, & PubChem data.",
            fz="sm", ff="Helvetica",
        ),
        dmc.Flex(
            [
                DashIconify(icon="ic:outline-email"),
                dmc.Text(
                    "yuansheng.liu@gmail.com",
                    fz="sm",
                    ff="Helvetica"
                ),
            ],
            gap="xs",
            justify="center",
            align="center"
        )
    ]
)

app.layout = dmc.MantineProvider(
    defaultColorScheme="light",
    theme={
        "primaryColor": "blue",       # set default color
        "fontFamily": "Nunito, sans-serif",  # change font globally
    },
    children=dmc.Container(
        style={
            "display": "flex", "flexDirection": "column", "minHeight": "100vh",
        },
        children=[
            navbar,
            dash.page_container,
            footer,
        ],
        fluid=True,
    ),
)

callbacks.register_callbacks(app)

if __name__ == "__main__":
    app.run(debug=True)
