import dash
from dash import dcc
import dash_mantine_components as dmc
from dash_iconify import DashIconify
dash.register_page(__name__, path="/cheminformatics/tanimoto", name="Tanimoto")

TANIMOTO_EQUATION = r"$J(A,B) = \dfrac{A \cap B}{A \cup B} = \dfrac{A \cap B}{|A| + |B| - A \cap B}$"

tanimoto_description = dmc.Stack(
    [
        dcc.Markdown(
            """
            The Tanimoto Coefficient, or Jaccard Index, is a statistical 
            metric used to assess the similarity between two sample sets. 
            Tanimoto coefficients range from a value of 0.0 to 1.0 due to how the index is calculated. 
            The equation below shows the method for calculating the Tanimoto index:
            """,
        ),
        dcc.Markdown(
            TANIMOTO_EQUATION,
            mathjax=True,
            style={"textAlign": "center"}
        ),
        dcc.Markdown(
            f"""
            where $J(A,B)$ represents the Jaccard Index/Tanimoto Coefficient and A and B represent distinct sample sets. 
            Figures 1a and 1b visually represent $A \\cap B$ and $A \\cup B$, respectively.
            The *accepted* value for a high degree of similarity between two sets is that $J(A,B) > 0.85$.
            In a biochemical context, even if $J(A,B)>0.85$, a high similarity does _NOT NECESSARILY_ indicate
            similar bioactivities.
            """,
            mathjax=True,
            dangerously_allow_html=True
        )
    ],
    fz=14,
)  

tanimoto_imgs = dmc.Stack(
    [
        dmc.Stack(
            [
                dmc.Image(
                    src="/assets/images/intersection_of_sets.png",
                    alt="Intersection Image",
                    w=200, h="auto", fit="contain",
                    style={"display": "block", "margin": "auto"},
                ),
                dmc.Text(
                    "Figure 1a. Intersection between A & B.",
                    ta="center",
                    fs="italic",
                    ff="Helvetica",
                    fz=10
                )
            ]
        ),
        dmc.Stack(
            [
                dmc.Image(
                    src="/assets/images/union_of_sets.png", # make sure to include leading "/"
                    alt="Union Image",
                    w=200, h="auto", fit="contain",
                    style={"display": "block", "margin": "auto"},
                    p="sm"
                ),
                dmc.Text(
                    "Figure 1b. Union between A & B.",
                    ta="center",
                    fs="italic",
                    ff="Helvetica",
                    fz=10
                )
            ]
        ),
    ],
    gap="xs",
    p="sm"
)

description = dmc.Box(
    [
        dmc.Title(
            "Introduction — What is a Tanimoto Coefficient?", order=3,
        ),
        dmc.Flex(
            [
                tanimoto_description,
                tanimoto_imgs
            ],
            gap="sm"
        ),
    ],
    # withBorder=True,
    # shadow="lg",
    # p="lg",
    # radius="xl",
    # mb=15,
    # miw=650,
)

tanimoto_applications = dmc.Stack(
    [
        dcc.Markdown(
            """
            Calculating the Tanimoto Similarity Coefficient is done in the context of several real-world applications.
            For example, biogeographical research by Chung et al used the Jaccard/Tanimoto coefficient
            when analyzing the co-occurences of bird and fish species in certain areas.<sup>1</sup> In cheminformatics, it is used
            in quantifying the similarity between molecules of interest.
            """,
            dangerously_allow_html=True,
        ),
        dcc.Markdown(
            """
            **In fact, the purpose of this of this page**
            is inspired by the methods used in finding molecular structural similarities
            demonstrated by the work of Hadipour et al.<sup>2</sup> Specifically, Hadipour et al
            had several molecules evaluated and mapped on a heatmap, or matrix. The molecules were treated as sample sets
            in the form of molecular fingerprints, or computational representations of molecule structures.
            Thus, this page aims to **generate a simple Tanimoto similarity heatmap/matrix** 
            given the user's requested molecules.
            Additionally, the generated matrix has clickable tiles that highlight similar regions/substructures
            of the two sets of molecules. Figure 2 shows a generated Tanimoto similarity matrix for the 22 amino acids using the program below.
            Figure 3 shows the highlighted similarities mapped onto the amino acids histidine and tryptophan.
            """,
            dangerously_allow_html=True
        ),
        dcc.Markdown(
            """
            The sets used in molecular similarity for the heatmap are molecular fingerprints (often in a string of text).
            These fingerprints are not limited
            to one form, as there are different representations of fingerprinting a molecule. Binary fingerprints, for example,
            represent the presence (1) or absence (0) of structural features while circular fingerprints considers
            an atom of a molecule's nearby atoms within a specific radius. Some molecular fingerprinting methods go
            an extra step and encode information on the frequency of occurences of each substructure within a molecule (count-based fingerprints).
            This Tanimoto demonstration uses a **count-based Morgan fingerprinting (C-MF)** algorithm in determining molecular similarities through RDKit.
            More instructions in generating a simple Tanimoto heatmap can be found below.
            """,
            dangerously_allow_html=True
        ),
    ],
    fz=14
)

tanimoto_applications_imgs = dmc.Stack(
    [
        dmc.Stack(
            [
                dmc.Image(
                    src="/assets/images/tanimoto-amino_acids.png",
                    alt="Intersection Image",
                    w=300, h="auto", fit="contain",
                    style={"display": "block", "margin": "auto"},
                ),
                dmc.Text(
                    "Figure 2. Tanimoto Similarity Matrix for 22 Amino Acids",
                    ta="center",
                    fs="italic",
                    ff="Helvetica",
                    fz=10
                )
            ]
        ),
        dmc.Stack(
            [
                dmc.Image(
                    src="/assets/images/tryptophan-histidine_sims.png",
                    alt="Intersection Image",
                    w=150, fit="contain",
                    style={"display": "block", "margin": "auto"},
                ),
                dmc.Text(
                    "Figure 3. Highlighted similarities between histidine (left) and tryptophan (right)",
                    ta="center",
                    fs="italic",
                    ff="Helvetica",
                    fz=10
                )
            ]
        ),
    ],
    gap="xs",
    p="sm"
)

references = dmc.Box(
    [
        dmc.Title("References", order=5),
        dmc.List(
            [
                dmc.ListItem(
                    dcc.Markdown(
                        """
                        Chung, N. C.; Miasojedow, B.; Startek, M.; Gambin, A. 
                        Jaccard/Tanimoto Similarity Test and Estimation Methods for 
                        Biological Presence-Absence Data. *BMC Bioinformatics* **2019**, 20 (S15). 
                        https://doi.org/10.1186/s12859-019-3118-5.
                        """
                    )
                ),
                dmc.ListItem(
                    dcc.Markdown(
                        """
                        Hamid Hadipour; Liu, C.; Davis, R. L.; Cardona, S. T.; Hu, P. 
                        Deep Clustering of Small Molecules at Large-Scale via Variational 
                        Autoencoder Embedding and K-Means. **2022**, 23 (S4). 
                        https://doi.org/10.1186/s12859-022-04667-1.
                        """
                    )
                )
            ],
            fz=10,
            type="ordered",
            ff="Helvetica",
        )
    ],
)

full_tanimoto_applications = dmc.Box(
    [
        dmc.Title(
            "Applications of the Tanimoto Coefficient", order=3,
        ),
        dmc.Flex(
            [
                tanimoto_applications,
                tanimoto_applications_imgs,
            ],
            gap="sm"
        ),
    ],
    # withBorder=True,
    # shadow="lg",
    # p="lg",
    # radius="xl",
    # mb=15,
)

tanimoto_instructions = dmc.Box(
    [
        dmc.Title("Instructions on Generating a Tanimoto Similarity Matrix", order=3),
        dmc.List(
            [
                dmc.ListItem(dmc.Text("Search up a molecule.")),
                dmc.ListItem(dmc.Text("Add the molecule-- molecule will appear in the space below.")),
                dmc.ListItem(dmc.Text("Keep adding until desired amount reached.")),
                dmc.ListItem(dmc.Text("Generate the matrix. The number of molecules added reflects the width of the matrix.")),
                dmc.ListItem(
                    dmc.Text(
                        "Click on any tile to view highlighted similarities and morgan fingerprints between two molecules. " \
                        "The tile clicked will have a Tanimoto coefficient between the corresponding " \
                        "two molecules labeled on the horizontal and vertical axes."
                    )
                ),
            ],
            type="ordered",
            spacing="md"
        )
    ],
    # withBorder=True,
    # shadow="lg",
    # p="md",
    # radius="xl",
    mb=15,
)

full_description = dmc.Container(
    [
        dmc.Title(
            "Tanimoto Similarity Matrix", order=1, fw=700, fz=60, ta="center", my=15,
        ),
        dmc.Text(
            "Explore & visualize molecular similarities.", fz=20, ta="center", my=15
        ),
        dmc.Divider(mb=10, size="md"),
        dmc.Flex(
            [
                description, 
                # dmc.Divider(orientation="vertical"), 
                full_tanimoto_applications,
            ], gap="md"
        ),
        dmc.Divider(mb=10, size="xs", variant="dashed"),
        references,
        dmc.Divider(mb=10, size="md"),
        tanimoto_instructions,
    ],
    fluid=True
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
    "Generate Matrix",
    disabled=True,
    color="green",
    leftSection=DashIconify(icon="fluent-mdl2:generate", width=25),
    id={"type": "_gen-graphic", "index": "heatmap"},
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
tanimoto_heatmap = dmc.Box(
    id="_tanimoto-heatmap-container",
)

click_info = dmc.Box(
    id="_click-tanimoto-tile",

)

similarity_maps = dmc.Box(
    id="_mol-mapping",
)

tanimoto_data = dmc.Box(
    [   
        dcc.Store(id="_store-items", data=[]),
        return_alert,
        mol_list,
        dcc.Store(id="_smiles-dict")
    ]
)

tanimoto_mapping = dmc.Container(
    [
        dmc.Flex(
            [
                tanimoto_heatmap,
                dmc.Stack(
                    [click_info, similarity_maps]
                )
            ],
            gap="md",
            # justify="center",
        ),
        
    ],
    fluid=True,
)

layout = dmc.Container(
    [
        full_description,
        dmc.Divider(mb=15, size="md"),
        search_group,
        dmc.Divider(mb=15, size="md"),
        tanimoto_data,
        dmc.Divider(mb=15, size="md"),
        tanimoto_mapping,
    ],
    fluid=True,
)
