import requests
import re
import pandas as pd

import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import rdDepictor

from myutils.mygraphics import convert_to_png, rgb_rescaled, create_color_palette
import myutils.mycheminformatics as mychinfo

EXP_PROPS = [
    "boiling point",
    "melting point",
    "flash point",
    "density",
    "vapor pressure",
    "decomposition",
    "dissociation constants",
    "solubility",
    "vapor density",
]

REAGENT_PROPS = [
    "melting point",
    "boiling point",
    "density",
    "first aid measures"
]

def search_compound(query, namespace="name"):
    compounds = pcp.get_compounds(query, namespace=namespace)
    try:
        if compounds:
            return compounds[0] # For simplified purposes, just the first result is taken.
        
        else:
            quoted_query = f'"{query}"'
            params = {
                "db": "pccompound",
                "term": quoted_query,
                "retmode": "json"
            }
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            resp = requests.get(url, params=params).json()

            cid = resp["esearchresult"]["idlist"][0]
            return pcp.get_compounds(int(cid), namespace="cid")[0]
            
    except pcp.PubChemHTTPError:
        return "ERR_REQ"
    except:
        return None


def get_suggested_compounds(query, limit=10):
    params = {"limit": limit}
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/{query}/json"
    resp = requests.get(url, params=params)

    if resp.ok:
        return resp.json().get("dictionary_terms", {}).get("compound", [])
    
    return []


def get_compound(name, namespace="name"):
    compound = search_compound(name, namespace)

    if compound and (compound != "ERR_REQ"):
        mol = Chem.MolFromInchi(compound.inchi)
        filter_syns = [
            name for name in compound.synonyms
            if not re.fullmatch(r"\d{2,7}-\d{2}-\d", name)
            and not name.startswith("DTXSID")
            and not name.startswith("DTXCID")
            and not name.startswith("SCHEMBL")
            and not name.startswith("UNII")
            and not name.startswith("NSC")
            and not name.startswith("ChemIDplus")
            and not name.startswith("InChI=")
            and not re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z0-9]", name)
        ]
        
        return {
            "name": filter_syns[0] if filter_syns else compound.iupac_name,
            "image": convert_to_png(
                mychinfo.draw_compound(
                    mol, 700, 700, mychinfo.RAW_PALETTE_LIGHT
                )
            ),

            "PubChem CID": {
                "value": compound.cid if compound.cid else "Unavailable",
                "icon": "stash:data-numbers-solid"
            },
            "Formula": {
                "value": compound.molecular_formula if compound.molecular_formula else "Unavailable",
                "icon": "arcticons:chemistry"
            },
            "Molar Mass": {
                "value": compound.molecular_weight if compound.molecular_weight else "Unavailable",
                "icon": "streamline-ultimate:science-molecule-strucutre",
                "units": "g/mol"
            },
            "IUPAC": {
                "value": compound.iupac_name if compound.iupac_name else "Unavailable",
                "icon": "gridicons:nametag",
            },
            "SMILES": {
                "value": compound.isomeric_smiles if compound.isomeric_smiles else "Unavailable",
                "icon": "gravity-ui:molecule"
            },
            "InChIKey": {
                "value": compound.inchikey if compound.inchikey else "Unavailable",
                "icon": "gravity-ui:molecule"
            },
            "InChI": {
                "value": compound.inchi if compound.inchi else "Unavailable",
                "icon": "gravity-ui:molecule"
            },
            "Formal Charge": {
                "value": compound.charge if compound.charge else "Unavailable",
                "icon": "simple-icons:ionic"
            },
            "XLogP3": {
                "value": compound.xlogp if compound.xlogp else "Unavailable",
                "icon": "mdi:flask-outline"
            },
        }

    elif compound == "ERR_REQ":
        return "ERR_REQ"
    else:
        return None


def extract_properties_sections(sections, ref_lookup, desired_heading=None):
    if not sections:
        return []
    
    results = []
    for section in sections:
        # Only include this section if heading matches, or if no filter is set
        toc_heading = section.get("TOCHeading", "").strip().lower()
        if desired_heading and toc_heading != desired_heading.lower():
            # Still check child sections
            if "Section" in section:
                results.extend(extract_properties_sections(section["Section"], ref_lookup, desired_heading))
            continue
        
        # Extract content under "Information"
        contents = section.get("Information", [])
        for info in contents:
            values = []
            if ("Value" in info) and ("StringWithMarkup" in info["Value"]):
                values = [item["String"] for item in info["Value"]["StringWithMarkup"]]
            
            elif ("Value" in info) and ("Number" in info["Value"]):
                numbers = info["Value"]["Number"]
                unit = info["Value"].get("Unit", "")

                # Combine number(s) with unit
                values = [f"{num} {unit}".strip() for num in numbers]

            ref_num = info.get("ReferenceNumber", -1)
            source_name = ref_lookup.get(ref_num, None)
            citation = info.get("Reference", [source_name])
            
            results.append({
                "values": values,
                "sources": citation,
            })
        
        if "Section" in section:
            # Extend results for all sections
            results.extend(extract_properties_sections(section["Section"], ref_lookup, desired_heading))
    
    return results


def get_more_compound_properties(cid, property_):
    params = {
        "heading": property_,
    }
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
    resp = requests.get(url, params=params).json()
    record = resp.get("Record", {})
    sections = record.get("Section", [])

    references = record.get("Reference", [])
    ref_lookup = {
        ref["ReferenceNumber"]: ref.get("SourceName", "N/A")
        for ref in references
    }
    
    return extract_properties_sections(sections, ref_lookup, desired_heading=property_)


def physical_descriptions(cid):
    desc = get_more_compound_properties(cid, "Physical Description")
    return desc if desc else []


def get_experimental_properties(cid):
    records = []
    props_data = [
        get_more_compound_properties(cid, prop)
        for prop in EXP_PROPS
    ]

    for prop, entries in zip(EXP_PROPS, props_data):
        for entry in entries:
            values = entry.get("values", [])
            sources = entry.get("sources", [])
            
            for value, source in zip(values, sources):
                records.append(
                    {
                        "property": prop,
                        "value": value,
                        "source": source
                    }
                )
    
    return pd.DataFrame(records).to_dict()


def create_reagent_table(compounds):
    pubchem_compounds = [get_compound(compound) for compound in compounds]
    cids = [compound["PubChem CID"]["value"] for compound in pubchem_compounds]

    names = [compound["name"] for compound in pubchem_compounds]
    formulas = [compound["Formula"]["value"] for compound in pubchem_compounds]
    molar_masses = [compound["Molar Mass"]["value"] for compound in pubchem_compounds]


    properties = []
    for i, cid in enumerate(cids):
        molecule_data = {
            "name": names[i],
            "formula": formulas[i],
            "molar mass": molar_masses[i]
        }

        for prop in REAGENT_PROPS:
            item = get_more_compound_properties(cid, prop)

            first_item = ""
            if item and prop != "first aid measures":
                first_item = item[0]["values"][0]
            
            elif prop == "first aid measures":
                first_item = [fa["values"][0] for fa in item] if item else None

            else:
                first_item = None

            molecule_data[prop] = first_item
        
        properties.append(molecule_data)
    
    return pd.DataFrame(properties)
