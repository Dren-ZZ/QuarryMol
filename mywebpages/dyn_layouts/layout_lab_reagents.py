import dash_mantine_components as dmc
from dash import dcc, html
from dash_iconify import DashIconify
from myutils import mygraphics, myplots
from myutils import mypubchem as mypcp
import pandas as pd


def display_molecules_list(items):
    molecules = []

    for item in items:
        mol_image = dmc.Image(
            src=item.get("image"),
            alt="structure",
            w=200, h=200,
            fit="contain",
            radius="lg",
            mb=5,
        ) if item.get("image") else dmc.Avatar("?", size="lg")

        mol_badge = dmc.Badge(
            children=item.get("name"),
            variant="light",
            leftSection=dmc.ActionIcon(
                DashIconify(icon="ep:close-bold", width=10,),
                id={"type": "remove-button", "index": item.get("name", "")},
                radius="xl", size=10,
            ),
            mb=15,
            size="xl",
            radius="xl",
            fz=16,
        )

        molecules.append(dmc.Stack([mol_image, mol_badge], gap="xs", align="center"))

    count_molecules = dmc.Text(f"{len(molecules)} molecules added.", mb=15)

    mol_grid = dmc.SimpleGrid(
        cols=6,
        spacing="md",
        children=molecules
    )

    return [count_molecules, mol_grid]


def create_reagent_table_layout(df):
    if df.empty:
        return None

    title = dmc.Group(
        [
            dmc.Title("Reagents", order=3, fw=700, fz=30,),
            DashIconify(icon="streamline-ultimate:lab-tube-bold", width=25),
        ],
        gap="xs"
    )

    rows = []
    for _, row in df.iterrows():
        row_cells = []
        for val in row:
            if isinstance(val, list):
                cell = dmc.List([dmc.ListItem(item) for item in val], size="sm", spacing="xs")

            elif pd.isna(val):
                cell = dmc.Text("\u2014", size="sm") # em dash
            
            else:
                cell = dmc.Text(val, size="sm")
            
            row_cells.append(dmc.TableTd(cell, fz=12, style={"verticalAlign": "top"}))
        
        rows.append(dmc.TableTr(row_cells))
    
    head = dmc.TableThead(
        dmc.TableTr(
            [dmc.TableTh(col.title()) for col in df.columns]
        )
    )

    caption = dmc.TableCaption("Missing values? Refer to several online sources. This table is merely a guide.")
    table = dmc.Table(
        [head, dmc.TableTbody(rows), caption,],
        withRowBorders=True,
        highlightOnHover=True,
        striped=True
    )
    scroll_table = dmc.TableScrollContainer(
        table, maxHeight=800, minWidth=1000, type="scrollarea"
    )

    pad_table = dmc.Card(
        [title, scroll_table,],
        withBorder=True, shadow="lg", p="md", radius="xl",
    )

    file_formats = dmc.Group(
        [
            dmc.RadioGroup(
                id="_reagents-export-format",
                value="csv",
                label="Download Format",
                children=dmc.Group([
                    dmc.Radio(label="CSV (.csv)", value="csv"),
                    dmc.Radio(label="Excel (.xlsx)", value="excel"),
                    dmc.Radio(label="JSON (.json)", value="json"),
                    dmc.Radio(label="Markdown (.md)", value="markdown"),
                ]),
                mb=15,
                size="sm"
            ),

            dmc.ActionIcon(
                DashIconify(icon="mdi:file-download", width=25),
                id="_download-table-btn",
                variant="light",
                loaderProps={"type": "dots"}
            ),
            dcc.Download(id="_download-reagents")
        ],
        justify="center"
    )

    return dmc.Stack([pad_table, file_formats,])

