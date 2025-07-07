import dash_mantine_components as dmc
from dash import dcc, html
from dash_iconify import DashIconify
from myutils import mygraphics, myplots
from myutils import mypubchem as mypcp
import pandas as pd
from rdkit import Chem


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


def display_tanimoto_matrix(z_vals, x_labels, y_labels, z_labels, caption, title, blurb):
    figure = myplots.create_heatmap(
        z_vals=z_vals,
        x_labels=x_labels,
        y_labels=y_labels,
        z_labels=z_labels,
        hover="x: %{x}<br>y: %{y}<br> Tanimoto: %{z}<extra></extra>",
        bar_title="Tanimoto Coefficient",
        map_title="Tanimoto Similarity Heatmap"        
    )
    heatmap = dcc.Graph(
        id="_tanimoto-heatmap",
        figure=figure,
        style={
            "display": "block", 
            "marginLeft": "auto", 
            "marginRight": "auto",
        }
    )

    heading = dmc.Title(title, order=3, fw=700, mb=15, fz=30,)
    new_blurb = dmc.Text(blurb)
    cap = dmc.Text(caption, ff="Helvetica", fz=10, fs="italic", ta="end")
    
    return dmc.Card(
        [heading, new_blurb, heatmap, cap],
        withBorder=True,
        shadow="lg",
        p="md",
        radius="xl",
        mb=15,
    )



