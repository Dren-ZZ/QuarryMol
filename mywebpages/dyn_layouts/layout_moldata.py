import dash_mantine_components as dmc
from dash import dcc, html
from dash_iconify import DashIconify
from myutils import mygraphics, myplots, mycheminformatics as mychinf
from myutils import mypubchem as mypcp
import pandas as pd
from rdkit import Chem


def create_component_properties(properties):
    properties_list = [
        dmc.ListItem(
            children=[dmc.Text(p["label"], fw=700), dmc.Code(p["value"])],
            icon=DashIconify(icon=p["icon"], width=30)
        )
        for p in properties
    ]
    
    format_list = dmc.List(
        properties_list,
        size="sm",
        spacing="md",
        styles={"itemLabel": {"whiteSpace": "nowrap"}},
        m=15,

    )

    return dmc.ScrollArea(format_list, type="hover", scrollHideDelay=0)


def create_pcp_properties_card(compound_name, properties):
    card = dmc.Card(
        [
            dmc.Title(compound_name, order=3, fz=30,),
            properties
        ],
        maw=800,
        h=400,
        shadow="md",
        radius="xl",
        withBorder=True,
    )

    return card


def create_image_card(src, alt, caption, idx):
    image = dmc.Image(
        src=src, alt=alt, w=200, h=200, radius="lg",
        style={"display": "block", "marginLeft": "auto", "marginRight": "auto"},
    )

    clickable_img = dmc.ActionIcon(
        image,
        id={"type": "card-image", "index": idx}, # Unique ID for callback.
        w=200, h=200,
        variant="subtle",
        
    )
    hidden_src = dmc.Box(
        src,
        id={"type": "card-image-src", "index": idx},
        style={"display": "none"}
    )
    
    cap = dmc.Text(caption, ff="Helvetica", fz=10, fs="italic", ta="end", mt=10)

    card = dmc.Card(
        [clickable_img, cap, hidden_src], p="md", shadow="md",
        radius="xl", withBorder=True,
    )

    return card


def create_pie_display(labels, values, caption):
    pie = dcc.Graph(
        figure=myplots.create_pie_chart(labels, values),
        style={
            "display": "block", "marginLeft": "auto", "marginRight": "auto",
        },
    )
    cap = dmc.Text(caption, ff="Helvetica", fz=10, fs="italic", ta="end")
    
    pad_pie = dmc.Card(
        [pie, cap],
        p="md", shadow="md", radius="xl",
        withBorder=True
    )
    
    return pad_pie


def create_modal_image(idx): # idx for gallery of images in case
    return dmc.Modal(
        id={"type": "image-modal", "index": idx},
        opened=False,
        withCloseButton=False,
        size="auto",
        radius="lg",
        zIndex=2,
        children=[
            dmc.Image(
                id={"type": "modal-image", "index": idx},
                alt="full"
            )
        ]
    )


def create_pcp_properties(compound, properties):
    properties_list = create_component_properties(properties)
    properties_card = create_pcp_properties_card(compound["name"], properties_list)
    image_card = create_image_card(compound["image"], alt="structure", caption="Made with RDKit", idx=0)
    modal_image = create_modal_image(idx=0)

    return dmc.Flex(
        [properties_card, image_card, modal_image,],
        gap="xs", wrap="nowrap", align="flex-start",
    )


def create_charts(compound):
    if not compound:
        return None
    
    title = dmc.Group(
        [
            dmc.Title("Charts & Graphs", order=3, fw=700, mb=15, fz=30,),
            DashIconify(icon="foundation:graph-pie", width=25),
        ],
        gap="xs"
    )


    smiles = compound["InChI"]["value"]
    mol = Chem.MolFromInchi(smiles)
    mol_hs = Chem.AddHs(mol)

    pie_data = mygraphics.create_atoms_pie_data(mol_hs)
    display_pie = create_pie_display(*pie_data, caption="Graphed with Plotly")
    
    charts = dmc.Flex(
        [
            display_pie,
        ],
        gap="xs", wrap="nowrap", align="flex-start",
    )

    pad_charts = dmc.Paper(
        [title, charts],
        withBorder=True,
        shadow="lg",
        p="md",
        radius="xl"
    )

    return pad_charts


def create_physdesc_carousel(descriptions):
    about_list = dmc.List(
        [
            dmc.ListItem(
                "Browse through each description from each source."
            ),

            dmc.ListItem(
                "Qualitative descriptions are all via PubChem, with their sources listed below.",
            ),

            dmc.ListItem(
                "These descriptions may seem brief. " \
                "Please visit the description's primary source for more information.",
            ),

            dmc.ListItem("Primary sources are labeled in the description's footer.")
        ], mb=10
    )

    title = dmc.Group(
        [
            dmc.Title("Physical/Qualitative Descriptions", order=3, fw=700, fz=30,),
            DashIconify(icon="flat-color-icons:about", width=25),
        ],
        gap="xs"
    )
    
    slides = [
        dmc.CarouselSlide(
            [
                dmc.Card(
                    [
                        dmc.ScrollArea(
                            dmc.Blockquote(
                                desc["values"], cite=desc["sources"],
                                ff="monospace", fz=12, ta="left",
                            ),
                            type="hover",
                            scrollHideDelay=0,
                            h=200
                        ),
                    ],
                    withBorder=True,
                    shadow="lg",
                    radius="xl",
                    mb=15,
                )
            ],
        )
        for desc in descriptions
    ]

    carousel = dmc.Carousel(
        slides,
        withIndicators=True,
        maw=700,
        emblaOptions = {"loop": True},
        classNames={"indicator": "dmc-indicator", "controls": "dmc-controls", "root": "dmc-root"}
    )

    pad_carousel = dmc.Card(
        [title, about_list, carousel],
        my="auto", withBorder=True,
        p="lg", shadow="lg", radius="xl", w="fit-content"
    )

    return pad_carousel


def create_physdesc(descriptions):
    if not descriptions:
        return []
    
    return create_physdesc_carousel(descriptions)


def create_exp_props_table(df_dict):
    if not df_dict:
        return None

    df = pd.DataFrame.from_dict(df_dict)
    title = dmc.Group(
        [
            dmc.Title("Experimental Properties", order=3, fw=700, fz=30,),
            DashIconify(icon="mdi:beaker-check-outline", width=25),
        ],
        gap="xs"
    )

    rows = []
    prop_groups = df.groupby("property")
    for prop, group in prop_groups:
        group_rows = group.reset_index() # Resets it to numbered default index

        for col, row in group_rows.iterrows():
            row_cells = []
            if col == 0: # First column
                row_cells.append(
                    dmc.TableTd(
                        prop,
                        tableProps={
                            "rowSpan": len(group_rows),
                            "width": 50
                        },
                        style={"verticalAlign": "top"},
                        fw=700,
                    )
                ) # Merges cells vertically

            row_cells.append(
                dmc.TableTd(
                    row["value"],
                    tableProps={"width": 300}, fz=12,
                    style={"verticalAlign": "top"},
                )
            )
            row_cells.append(
                dmc.TableTd(
                    row["source"], 
                    fs="italic", fz=10,
                    style={"verticalAlign": "top"},
                )
            )
            rows.append(dmc.TableTr(row_cells))

    head = dmc.TableThead(
        dmc.TableTr(
            [dmc.TableTh("Property"), dmc.TableTh("Value"), dmc.TableTh("Source")]
        )
    )
    caption = dmc.TableCaption("Table of selected experimental properties with their respective sources.")
    table = dmc.Table(
        [head, dmc.TableTbody(rows), caption],
        withRowBorders=True,
        highlightOnHover=True,
    )
    scroll_table = dmc.TableScrollContainer(
        table, maxHeight=800, minWidth=1000, type="scrollarea"
    )

    pad_table = dmc.Card(
        [title, scroll_table],
        withBorder=True, shadow="lg", p="md", radius="xl",
    )

    return pad_table
