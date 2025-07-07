import dash
from dash import Input, Output, State, dcc, clientside_callback
from dash import callback_context as ctx
import myutils.mypubchem as mypubchem
import dash_mantine_components as dmc
from dash_iconify import DashIconify
from rdkit import Chem
import pandas as pd
import numpy as np
import io
from tabulate import tabulate

from mywebpages.dyn_layouts import layout_moldata
from mywebpages.dyn_layouts import layout_cheminfo_tanimoto
from mywebpages.dyn_layouts import layout_lab_reagents

from myutils.mycheminformatics import generate_tanimoto_matrix, show_sims, draw_3D_molecule

def register_callbacks(app):
    @app.callback(
        Output({"type": "suggested-list", "index": "autocomplete"}, "data"),
        Input({"type": "suggested-list", "index": "autocomplete"}, "value"),
        prevent_initial_call=True,
    )
    def update_suggestions(query):
        if not query:
            return []
        
        results = mypubchem.get_suggested_compounds(query)
        if not results:
            return ["No results found"]
        
        return results
    
    @app.callback(
        [
            Output("_search-output-moldata", "data"),
        ],
        [
            Input("_search-internal", "n_clicks"),
            State({"type": "suggested-list", "index": "autocomplete"}, "value"),
        ],
        prevent_initial_call=True,
        running=[(Output("_search-internal", "loading"), True, False)],
    )
    def fetch_results_moldata(n_clicks, value):
        if not value:
            return [{"error": "Please enter a query."}]
        
        compound = mypubchem.get_compound(value)
        if not compound:
            return [{"error": f'"{value}" not found.'}]
        
        skip_keys = {"image", "name", "image3D"}
        descriptions = mypubchem.physical_descriptions(compound["PubChem CID"]["value"])
        exp_props = mypubchem.get_experimental_properties(compound["PubChem CID"]["value"])

        properties = [
            {
                "label": f"{k} ({v['units']})" if "units" in v else k,
                "value": v["value"],
                "icon": v["icon"]
            }
            for k, v in compound.items()
            if k not in skip_keys
        ]

        return [
            {
                "query": value,
                "compound": compound,
                "descriptions": descriptions,
                "image": compound["image"],
                "name": compound["name"],
                "pcp props": properties,
                "exp props": exp_props,
            }
        ]
    
    @app.callback(
        Output("_search-output", "children"),
        Input("_search-output-moldata", "data"),
        prevent_initial_call=True,
    )
    def render_results_moldata(data):
        if "error" in data:
            return dmc.Alert(data["error"], title="Error", color="red", variant="light", withCloseButton=True)
        
        value = data["query"]
        compound = data["compound"]
        descriptions = data["descriptions"]
        pcp_props = data["pcp props"]
        exp_props = data["exp props"]

        pcp_properties = layout_moldata.create_pcp_properties(compound, pcp_props)
        qual_props = layout_moldata.create_physdesc(descriptions)
        exp_table = layout_moldata.create_exp_props_table(exp_props)
        charts = layout_moldata.create_charts(compound)
        
        all_properties = [
            dmc.Paper(
                dmc.Title(
                    f'Results for "{value}"', order=2,
                    mb=20, ff="Quicksand", fz=40
                )
            ),
            dmc.Divider(size="md"),
            dmc.Box(pcp_properties, mb=10, p="md", w="fit-content"),
            dmc.Divider(size="md"),
            dmc.Box(charts, mb=10, p="md", w="fit-content"),
            dmc.Divider(size="md"),
            dmc.Box(qual_props, mb=10, p="md", w="fit-content"),
            dmc.Divider(size="md"),
            dmc.Box(exp_table, mb=10, p="md", w="fit-content"),
        ]

        return all_properties
    
    @app.callback(
        [
            Output({"type": "image-modal", "index": dash.MATCH}, "opened"),
            Output({"type": "modal-image", "index": dash.MATCH}, "src"),
        ],

        Input({"type": "card-image", "index": dash.MATCH}, "n_clicks"),
        State({"type": "card-image-src", "index": dash.MATCH}, "children"),
        prevent_initial_call=True,
    )
    def open_image_modal(n_clicks, src):
        return True, src

    @app.callback(
        Output({"type": "_gen-graphic", "index": dash.MATCH}, "disabled"),
        Input("_store-items", "data"),
        prevent_initial_call=True,
    )
    def toggle_button(items):
        # The `dash.MATCH` output corresponds to a single button
        triggered_output = ctx.outputs_list
        button_index = triggered_output["id"].get("index")

        if button_index == "heatmap":
            return len(items) < 2
        
        else:
            return len(items) < 1
    
    @app.callback(
        [
            Output("_store-items", "data"),
            Output("_feedback-alert", "children"),
        ],
        [
            Input("_add-molecule", "n_clicks"),
            Input({"type": "remove-button", "index": dash.ALL}, "n_clicks"),
        ],
        [
            State({"type": "suggested-list", "index": "autocomplete"}, "value"),
            State("_store-items", "data"),
        ],
        prevent_initial_call=True,
        running=[(Output("_add-molecule", "loading"), True, False)],
    )
    def update_mol_list(add_clicks, remove_clicks, query, items):
        triggered_id = ctx.triggered_id
        feedback = None

        if not query:
            feedback = dmc.Alert(
                f"Please enter a query.",
                withCloseButton=True,
                title="Error", 
                color="red",
                mb=15)

        else:
            compound = mypubchem.get_compound(query)
            
            # Handle add
            if triggered_id == "_add-molecule":
                if not compound:
                    feedback = dmc.Alert(
                        f'"{query} not found."',
                        withCloseButton=True,
                        title="Error",
                        color="red",
                        mb=15
                    )

                else:
                    name = compound.get("name")
                    if name and name.strip() not in [item["name"] for item in items]:
                        items.append(
                            {
                                "name": name.strip(),
                                "image": compound.get("image"),
                            }
                        )

                        feedback = dmc.Alert(
                            f'"{name}" successfully added.', 
                            withCloseButton=True,
                            title="Success", 
                            color="green",
                            mb=15
                        )
                    else:
                        feedback = dmc.Alert(
                            f'"{name}" already added.', 
                            withCloseButton=True,
                            title="Error", 
                            color="red", 
                            mb=15
                        )
            
            else:
                # Handle remove
                if triggered_id.get("type") == "remove-button":
                    target = triggered_id["index"]
                    items = [item for item in items if item["name"] != target]
                    feedback = dmc.Alert(
                        f'"{target}" successfully removed.', 
                        withCloseButton=True,
                        title="Success", 
                        color="green", 
                        mb=15
                    )


        return [items, feedback]
    

    @app.callback(
        [
            Output("_num-molecules", "children"),
            Output("_list-molecules", "children"),
        ],
        Input("_store-items", "data"),
        prevent_initial_call=True,
    )
    def render_list(items):
        mol_grid = layout_cheminfo_tanimoto.display_molecules_list(items)
        return mol_grid

    
    @app.callback(
        [
            Output("_tanimoto-heatmap-container", "children"),
            Output("_smiles-dict", "data")
        ],

        Input({"type": "_gen-graphic", "index": "heatmap"}, "n_clicks"),
        State("_store-items", "data"),
        prevent_initial_call=True,
        running=[(Output({"type": "_gen-graphic", "index": "heatmap"}, "loading"), True, False)],
    )
    def render_tanimoto(n_clicks, items):
        names = [item.get("name") for item in items]
        inchis = [mypubchem.get_compound(item.get("name"))["InChI"]["value"] for item in items]
        mol_rdkits = [Chem.MolFromInchi(inchi) for inchi in inchis]

        smiles_dict = {name: inchi for (name, inchi) in zip(names, inchis)} # Originally SMILES, now InChI-- may revert back once SMILES returns

        tanimoto = generate_tanimoto_matrix(mol_rdkits)
        tanimoto_matrix = tanimoto[1]
        tanimoto_text = np.where(np.isnan(tanimoto_matrix), "",
            np.round(tanimoto_matrix, 2).astype(str)
        )

        figure = layout_cheminfo_tanimoto.display_tanimoto_matrix(
            z_vals=tanimoto_matrix,
            x_labels=names,
            y_labels=names,
            z_labels=tanimoto_text,
            caption="Graphed with Plotly; Calculated with RDKit;",
            title="Generated Matrix",
            blurb="Click on any tile to visualize highlighted similarities."
        )
        
        return [figure, smiles_dict]
    

    @app.callback(
        [
            Output("_click-tanimoto-tile", "children"),
            Output("_mol-mapping", "children"),
        ],
        [
            Input("_tanimoto-heatmap", "clickData"),
            Input("_smiles-dict", "data")
        ],
        prevent_initial_call=True
    )
    def display_clicked_tile(clickData, load_mols_dict):
        if clickData is None:
            return ["", None]

        point = clickData["points"][0]
        mol_name_x = point["x"]
        mol_name_y = point["y"]
        mol_names = [mol_name_x, mol_name_y]
        tanimoto = point["z"]

        output_mols = [
            dmc.Text(f"Molecule X: {mol_name_x}"),
            dmc.Text(f"Molecule Y: {mol_name_y}"),
            dmc.Text(f"Tanimoto Coefficient: {tanimoto}")
        ]
            

        mol_x = Chem.MolFromInchi(load_mols_dict[mol_name_x])
        mol_y = Chem.MolFromInchi(load_mols_dict[mol_name_y])
        mol_img_tags = show_sims(mol_x, mol_y, highlight_color=(166/255, 219/255, 252/255, 1))

        mol_imgs = [
            dmc.Stack(
                [
                    dmc.Image(
                        src=tag,
                        radius="lg",
                        w=200,
                        h="auto"
                    ),
                    dmc.Text(name, ta="center", fz=10)
                ]
            )
            for (tag, name) in zip(mol_img_tags, mol_names)
        ]

        return [
            dmc.Card(
                children=output_mols,
                withBorder=True,
                shadow="lg",
                p="md",
                radius="xl",
                mb=15,
            ), 
            dmc.Card(
                dmc.Flex(children=mol_imgs, gap="sm",),
                withBorder=True,
                shadow="lg",
                p="md",
                radius="xl",
                mb=15,
            )
        ]
    
    
    @app.callback(
        [
            Output("_reagents-table-container", "children"),
            Output("_reagents-table-data", "data"),  # store data
        ],
        Input({"type": "_gen-graphic", "index": "table"}, "n_clicks"),
        State("_store-items", "data"),
        prevent_initial_call=True,
        running=[(Output({"type": "_gen-graphic", "index": "table"}, "loading"), True, False)],
    )
    def render_reagents(n_clicks, items):
        names = [item.get("name") for item in items]
        df = mypubchem.create_reagent_table(names)

        # store as list of dicts for JSON-safe transport
        json_data = df.to_dict(orient="records")

        table = layout_lab_reagents.create_reagent_table_layout(df)
        return [table, json_data]
    

    @app.callback(
        Output("_download-reagents", "data"),
        Input("_download-table-btn", "n_clicks"),
        [
            State("_reagents-table-data", "data"),
            State("_reagents-export-format", "value"),
        ],
        prevent_initial_call=True,
        running=[(Output("_download-table-btn", "loading"), True, False)],
    )
    def download_table(n_clicks, table_data, format_choice):
        if not table_data:
            return dash.no_update

        # Restore DataFrame
        df = pd.DataFrame(table_data)

        # Format logic
        if format_choice == "excel":
            df = df.applymap(lambda x: "\n".join(x) if isinstance(x, list) else x)
            buffer = io.BytesIO()
            df.to_excel(buffer, index=False, engine="openpyxl")
            buffer.seek(0)
            return dcc.send_bytes(buffer.getvalue(), filename="reagents.xlsx")

        elif format_choice == "csv":
            df = df.applymap(lambda x: "; ".join(x) if isinstance(x, list) else x)
            buffer = io.StringIO()
            df.to_csv(buffer, index=False)
            buffer.seek(0)
            return dcc.send_string(buffer.getvalue(), filename="reagents.csv")

        elif format_choice == "json":
            buffer = io.StringIO()
            df.to_json(buffer, orient="records", indent=2)
            buffer.seek(0)
            return dcc.send_string(buffer.getvalue(), filename="reagents.json")

        elif format_choice == "markdown":
            df = df.applymap(lambda x: "; ".join(x) if isinstance(x, list) else x)
            markdown_table = tabulate(df, headers="keys", tablefmt="github")
            return dcc.send_string(markdown_table, filename="reagents.md")

        return dash.no_update

