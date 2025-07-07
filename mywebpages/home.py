import dash_mantine_components as dmc
from dash_iconify import DashIconify
import dash
from dash import dcc, html


from myutils.mypubchem import get_compound

dash.register_page(__name__, path="/", name="Home")
compound_1 = get_compound("deoxyribonucleic acid")
compound_2 = get_compound("lornoxicam")

hero = dmc.Container(
    children=[
       dmc.Flex(
            [
                dmc.Stack(
                    [
                        dmc.Title(
                            "Welcome to QuarryMol", order=1,
                            fz=100, ta="center", fw=700,
                        ),
                        dmc.Text(
                            "Site for molecular data, cheminformatics, & lab facilitation.",
                            ta="center", my=15, fz=30,
                        ),
                    ],
                    gap="xs"
                ),
                dmc.Image(
                    src=compound_2.get("image"),
                    alt=compound_2.get("name", None),
                    w=250,
                    h=250,
                    radius="lg",
                )
            ],
            gap="xl",
            align="stretch",
            justify="center",
            p="xl"
        ),
        dmc.Flex(
            [
                dmc.Card(
                    [
                        dmc.Title(
                            "✦ Available Features ✦", order=1, fz=40, my=15, c="green"
                        ),
                        dmc.Text("Check them out in the cards below.", fz=20)
                    ],
                    w=500,
                    shadow="md",
                    radius="lg",
                    ta="center",
                    style={
                        "background": "linear-gradient(135deg, var(--mantine-color-green-1), var(--mantine-color-blue-1))"
                    }
                ),
                dmc.Divider(orientation="vertical", size="md", color="dark",),
                dmc.List(
                    [
                        dmc.ListItem(
                            ["Explore molecules."], c="blue"
                        ),
                        dmc.ListItem(
                            ["Analyze molecular data."], c="cyan"
                        ),
                        dmc.ListItem(
                            ["Assist with laboratory work."], c="teal"
                        )
                    ],
                    fz=25,
                    my=15,
                    spacing="md",
                    fw=700,
                    icon=DashIconify(icon="tabler:pin", width=30)
                )
            ],
            my=15,
            gap="md",
            justify="center",
            align="stretch"
        ),
        dmc.Center(
            html.A(
                dmc.ActionIcon(
                    DashIconify(
                        icon="tabler:chevron-down",
                        width=50,
                    ),
                    w=50,
                    h="auto",
                    radius="md",
                    variant="light"
                ),
                href="#features",
                style={"marginTop": 25, "textDecoration": "none"},
            )
        ),
    ],
    fluid=True,
    style={
        "height": "100vh", 
        "display": "flex", 
        "flexDirection": "column", 
        "justifyContent": "center",
    },
    bg="var(--mantine-color-gray-0)",
)

feature_card_1 = dmc.Card(
    p="md",
    children=[
        dmc.Title(
            "Molecular Search", order=3, fz=30, fw=700, mb=15, ta="center", c="blue"
        ),

        dmc.CardSection([
            dmc.Image(
                src=compound_1.get("image"),
                alt=compound_1.get("name", None),
                w=300,
                h=300,
                style={"display": "block", "margin": "auto"},
                radius="lg"
            ),
            dmc.Container(
                dmc.Anchor(
                    compound_1.get("name"), fs="italic", fz=10,
                    href=f"https://pubchem.ncbi.nlm.nih.gov/compound/{compound_1['PubChem CID']['value']}",
                    target="_blank",
                ),
                ta="right",
            )
        ]),

        dmc.Text(
            "Look up PubChem molecules by name and view their properties.",
            fz=16, m="md"
        ),
        dmc.Center(
            dmc.Anchor(
                dmc.Button(
                    "Search Molecule",
                    fz=18,
                    w="fit-content",
                    leftSection=DashIconify(
                        icon="line-md:file-search",
                        width=25,
                    ),
                    variant="light",
                    size="md",
                    radius="xl",
                    mb=15,
                    id="_moldata-btn"           
                ),
                href="/moldata",
            ),
        ),

    ],
    shadow="md",
    radius="xl",
    withBorder=True,
    w=350,
    # bg="var(--mantine-color-blue-1)",
)

feature_card_2 = dmc.Card(
    p="md",
    children=[
        dmc.Title(
            "Cheminformatics", order=3, fz=30, fw=700, ta="center", c="cyan"
        ),
        dmc.CardSection(
            dmc.Image(
                src="/assets/images/tanimoto-alkanes.png",
                alt="tanimoto",
                w=300,
                h=300,
                style={"display": "block", "margin": "auto"},
                radius="lg"
            )
        ),
        dmc.Text(
            "Visualize molecular similarities and run fingerprint comparisons.",
            fz=16, m="md"
        ),
        dmc.Center(
            dmc.Anchor(
                dmc.Button(
                    "Explore Tools",
                    fz=18, w="fit-content",
                    leftSection=DashIconify(
                        icon="lets-icons:molecule-light",
                        width=25,
                    ),
                    variant="light",
                    size="md",
                    radius="xl",
                    mb=15,
                ),
                href="/cheminformatics",
                id="_cheminformatics-btn"
            ),
        ),

    ],
    shadow="md",
    radius="xl",
    withBorder=True,
    w=350,
    # bg="var(--mantine-color-cyan-1)",
)

feature_card_3 = dmc.Card(
    p="md",
    children=[
        dmc.Title(
            "Lab Assistant", order=3, fz=30, fw=700, ta="center", c="teal"
        ),
        dmc.CardSection(
            dmc.Image(
                src="/assets/images/aspirin reagents.png",
                alt="Aspirin Reagents",
                w=300,
                h="auto",
                style={"display": "block", "margin": "auto"},
                radius="lg"
            )
        ),
        dmc.Text(
            "Generate reagent tables & other lab help.",
            fz=16, m="md"
        ),
        dmc.Center(
            dmc.Anchor(
                dmc.Button(
                    "Explore Lab Helpers",
                    fz=18, w="fit-content",
                    leftSection=DashIconify(
                        icon="solar:benzene-ring-outline",
                        width=25,
                    ),
                    variant="light",
                    size="md",
                    radius="xl",
                    mb=15,
                    id="_lab-btn"
                ),
                href="/lab"
            ),
        ),
    ],
    shadow="md",
    radius="xl",
    withBorder=True,
    w=350,
    # bg="var(--mantine-color-teal-1)",
)

row_cards = dmc.Stack(
    [
        dmc.Title("Features", order=1, fw=700, ta="center", fz=60),
        dmc.Flex(
            justify="center", gap="xl", align="stretch",
            children=[feature_card_1, feature_card_2, feature_card_3], mb=20
        ),
    ],
    id="features",
)



about = dmc.Flex(
    [
        dmc.Card(
            [
                dmc.Title("Mission", order=2, my=15, fz=30,),
                dmc.Text("QuarryMol offers chemistry and lab-related tools that streamline chemical education.", fz=20)
            ],
            p="md",
            w=400
        ),

        dmc.Card(
            dmc.CardSection(
                dmc.Image(
                    src="/assets/images/tryptophan-histidine_sims.png",
                    alt="similarity highlights",
                    w=400,
                    style={"display": "block", "margin": "auto"},
                    radius="lg"
                )
            ),
        )

    ],
    justify="center", gap="xs", align="center",
)

layout = dmc.Container(
    [
        hero,
        dmc.Divider(my=10, size="md"),
        row_cards,
        dmc.Divider(my=10, size="md"),
        about
    ],
    fluid=True,
)
