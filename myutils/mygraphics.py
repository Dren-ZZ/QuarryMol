from PIL import Image
from io import BytesIO
import base64

from rdkit import Chem

PT = Chem.GetPeriodicTable()


def create_img_tag(img):
    buf = BytesIO()
    img.save(buf, format="PNG")
    img_b64 = base64.b64encode(buf.getvalue()).decode("utf-8")
    return f"data:image/png;base64,{img_b64}"


def convert_to_png(png_bytes):
    img = Image.open(BytesIO(png_bytes))
    return create_img_tag(img)


def rgb_rescaled(r, g, b, a=1):
    return (r/255, g/255, b/255, a)


def rgb_raw(r, g, b, a=1):
    return (r*255, g*255, b*255, a)

def create_color_palette(mol, my_palette, default):
    palette = {}
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        rgb = my_palette.get(atomic_num, default)
        palette[atomic_num] = rgb_rescaled(*rgb)
    
    return palette


def create_atoms_pie_data(mol):
    # colors = [rgb_raw(*color) for color in palette.values()]
    mol_hs = Chem.AddHs(mol)
    frequencies = {}

    all_atoms = mol_hs.GetAtoms()
    atom_symbols = [atom.GetSymbol() for atom in all_atoms]

    for atom in atom_symbols:
        if atom not in frequencies:
            frequencies[atom] = 1

        else:
            frequencies[atom] += 1
    
    
    labels = list(frequencies.keys())
    values = list(frequencies.values())

    return labels, values


