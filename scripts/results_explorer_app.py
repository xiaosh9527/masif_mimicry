import copy
import json
import os
import uuid
from io import StringIO

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import py3Dmol
from Bio.PDB import PDBIO, PDBParser
from dash import ALL, Dash, Input, Output, State, callback_context, dcc, html, no_update


# ------ Data paths ------
# Path to results file
results_path = "/scratch/ymeng/masif_seed/masif_mimicry/data/NUP98/search_results/021structure_C_mimicry_library_102aa_260522.csv"

# Path to target .pdb file
target_path = "/scratch/ymeng/masif_seed/masif_mimicry/data/NUP98/input/021structure_AB.pdb"

# Info about the ligand structure - this will be displayed as sticks
lig_name = "021"
lig_chain = "A"
lig_resnum = 1

# ------ Configuration ------
# Path to database of target structures
# domainome_db_dir = "/work/lpdi/users/shxiao/masif_seed/masif/data/masif_human_proteome_domains_merged/data_preparation/01-benchmark_pdbs"
# domainome_db_dir = "/work/lpdi/users/diazrovi/domaindome/20260221-AFDBv6_domaindome_DPAM_masif/dpam_domaindome_masif_db/data_preparation/01-benchmark_pdbs"
domainome_db_dir = "/scratch/ymeng/TED_domainome/output/data_preparation/01-benchmark_pdbs"

# Default column names
matched_protein_colname = "P1_id"             # column name of matched_protein in results file
flattened_transform_colname = "flattened_transform"     # column name of flattened_transform in results file
target_iface_resi_colname = "target_iface_resi"         # column name of target interface residue number in results file (optional, used to highlight target interface residues)
matched_iface_resi_colname = "matched_iface_resi"       # column name of matched interface residue number in results file (optional, used to highlight matched interface residues)


# Default filter (set to None to start unfiltered with no pre-loaded criteria)
#default_filter_path = None
default_filter_path = "/scratch/ymeng/masif_seed/masif_mimicry/data/NUP98/search_results/021structure_C_mimicry_postprocessed_filters_102aa_260522.json"

# Default plotting fields (set in notebook comments)
DEFAULT_X_AXIS = "MaSIF-score"
DEFAULT_Y_AXIS = "P1_source_TMscore"
DEFAULT_COLOR_BY = "cluster_size"
UNIPROT_ACCN_COL = "uniprot_accn"
VIS_ROWS_FILENAME = "visualized_rows.csv"
VIS_PML_FILENAME = "visualization.pml"
INSPECTION_COL = "inspection"

# ------ Functions ------
def apply_transform(structure, transform_str):
    """Apply a 4x4 affine transform to all atoms in a Bio.PDB structure."""
    transform_flat = np.fromstring(transform_str, sep=",")
    if transform_flat.size != 16:
        raise ValueError(f"Expected flattened transform of size 16 (4x4), got {transform_flat.size}")
    transform_matrix = transform_flat.reshape((4, 4))
    rotation = transform_matrix[:3, :3]
    translation = transform_matrix[:3, 3]
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coord = atom.get_coord()
                    new_coord = np.dot(rotation, coord) + translation
                    atom.set_coord(new_coord)


def structure_to_pdbstr(structure):
    """Convert Bio.PDB structure to a PDB string."""
    io = PDBIO()
    io.set_structure(structure)
    fileobj = StringIO()
    io.save(fileobj)
    return fileobj.getvalue()


def build_residue_bfactor_color_bins(pdb_str):
    """Map residues to pLDDT-style color bins based on B-factor."""
    residue_max_b = {}
    for line in pdb_str.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        try:
            bfactor = float(line[60:66].strip())
        except ValueError:
            bfactor = 0.0
        chain = line[21].strip() or ""
        resi = line[22:26].strip()
        if not resi:
            continue
        key = (chain, resi)
        prev = residue_max_b.get(key)
        residue_max_b[key] = bfactor if prev is None else max(prev, bfactor)

    bins = {"orange": {}, "yellow": {}, "cyan": {}, "blue": {}}
    for (chain, resi), bfactor in residue_max_b.items():
        if bfactor <= 50:
            color = "orange"
        elif bfactor <= 70:
            color = "yellow"
        elif bfactor < 90:
            color = "cyan"
        else:
            color = "blue"
        bins[color].setdefault(chain, []).append(resi)
    return bins


def apply_filtering_criteria(df, filtering_criteria):
    """Apply filtering criteria in the format {column: (min_val, max_val)}."""
    filtered = df
    for col, (min_val, max_val) in filtering_criteria.items():
        if col not in filtered.columns:
            continue
        if min_val is not None:
            filtered = filtered[filtered[col] >= min_val]
        if max_val is not None:
            filtered = filtered[filtered[col] <= max_val]
    return filtered


def build_filtering_criteria(filter_rows):
    """Normalize filter row state into filtering_criteria dict."""
    criteria = {}
    for row in filter_rows or []:
        col = row.get("column")
        if not col:
            continue
        min_val = row.get("min")
        max_val = row.get("max")
        criteria[col] = (min_val, max_val)
    return criteria


def normalize_imported_filters(imported_payload, numeric_columns_set):
    """Convert saved JSON payload into filter-store row format."""
    if not isinstance(imported_payload, dict):
        raise ValueError("Invalid filter JSON: expected an object mapping columns to {min,max}.")

    rows = []
    ignored_columns = []
    for col, bounds in imported_payload.items():
        if col not in numeric_columns_set:
            ignored_columns.append(col)
            continue
        if not isinstance(bounds, dict):
            raise ValueError(f"Invalid filter JSON entry for '{col}': expected an object.")
        rows.append(
            {
                "id": str(uuid.uuid4()),
                "column": col,
                "min": bounds.get("min"),
                "max": bounds.get("max"),
            }
        )
    return rows, ignored_columns


def load_default_filter_state(default_path, numeric_columns_set):
    """Load initial filter-store rows and dataset mode from default_filter_path."""
    if not default_path:
        return [], "unfiltered"
    path = default_path.strip() if isinstance(default_path, str) else ""
    if not path:
        return [], "unfiltered"
    if not os.path.isfile(path):
        print(f"Warning: default_filter_path not found: {path}")
        return [], "unfiltered"
    with open(path, "r", encoding="utf-8") as handle:
        imported_payload = json.load(handle)
    rows, ignored_columns = normalize_imported_filters(imported_payload, numeric_columns_set)
    if ignored_columns:
        print(f"Warning: default filter ignored columns: {', '.join(ignored_columns)}")
    if not rows:
        print(f"Warning: no valid filters loaded from {path}")
        return [], "unfiltered"
    return rows, "filtered"


def view_predicted_binder_html(
    row,
    matched_protein_colname,
    flattened_transform_colname,
    target_struct,
    domainome_db_dir,
    parser,
    lig_chain=None,
    lig_name=None,
    lig_resi=None,
    target_iface_resi_colname=None,
    matched_iface_resi_colname=None,
):
    """Return py3Dmol view HTML for the target-ligand-binder complex."""
    matched_color = "cyan"
    target_color = "lightgrey"
    ligand_color = "green"

    matched_id = row.get(matched_protein_colname)
    if not isinstance(matched_id, str) or not matched_id.strip():
        raise ValueError(f"Missing or invalid '{matched_protein_colname}' for selected row.")

    transform_value = row.get(flattened_transform_colname)
    if not isinstance(transform_value, str) or not transform_value.strip():
        raise ValueError(f"Missing or invalid '{flattened_transform_colname}' for selected row.")

    matched_protein_path = os.path.join(domainome_db_dir, f"{matched_id}.pdb")
    if not os.path.exists(matched_protein_path):
        raise FileNotFoundError(f"Matched protein PDB not found: {matched_protein_path}")

    matched_protein_struct = parser.get_structure("matched_protein", matched_protein_path)
    apply_transform(matched_protein_struct, transform_value)

    target_pdb_str = structure_to_pdbstr(target_struct)
    matched_pdb_str = structure_to_pdbstr(matched_protein_struct)
    matched_residue_color_bins = build_residue_bfactor_color_bins(matched_pdb_str)

    view = py3Dmol.view(width=700, height=600)
    view.addModel(target_pdb_str, "pdb")
    view.setStyle({"model": 0}, {"cartoon": {"color": target_color}, "line": {"colorscheme": f"{target_color}Carbon"}})
    view.addModel(matched_pdb_str, "pdb")
    view.setStyle({"model": 1}, {"cartoon": {"color": "orange"}, "line": {"colorscheme": "orangeCarbon"}})
    for color_name in ["yellow", "cyan", "blue"]:
        for chain, residues in matched_residue_color_bins[color_name].items():
            if not residues:
                continue
            selection = {"model": 1, "resi": residues}
            if chain:
                selection["chain"] = chain
            view.setStyle(
                selection,
                {"cartoon": {"color": color_name}, "line": {"colorscheme": f"{color_name}Carbon"}},
            )

    if (
        target_iface_resi_colname
        and target_iface_resi_colname in row
        and isinstance(row[target_iface_resi_colname], str)
        and row[target_iface_resi_colname].strip() != ""
    ):
        iface_resi_list = [res.strip() for res in row[target_iface_resi_colname].split(",") if res.strip()]
        for resi in iface_resi_list:
            view.setStyle(
                {"model": 0, "resi": resi},
                {
                    "stick": {"colorscheme": f"{target_color}Carbon", "radius": 0.25},
                    "cartoon": {"color": target_color},
                    "line": {"colorscheme": f"{target_color}Carbon"},
                },
            )

    if (
        matched_iface_resi_colname
        and matched_iface_resi_colname in row
        and isinstance(row[matched_iface_resi_colname], str)
        and row[matched_iface_resi_colname].strip() != ""
    ):
        matched_iface_resi_list = [res.strip() for res in row[matched_iface_resi_colname].split(",") if res.strip()]
        for resi in matched_iface_resi_list:
            iface_color = "orange"
            for color_name in ["yellow", "cyan", "blue"]:
                if any(resi in residues for residues in matched_residue_color_bins[color_name].values()):
                    iface_color = color_name
            view.setStyle(
                {"model": 1, "resi": resi},
                {
                    "cartoon": {"color": iface_color},
                    "line": {"colorscheme": f"{iface_color}Carbon"},
                    "stick": {"colorscheme": f"{iface_color}Carbon", "radius": 0.25},
                },
            )

    if lig_chain is not None and lig_name is not None and lig_resi is not None:
        ligand_style_selector = {"model": 0, "chain": lig_chain, "resn": lig_name, "resi": str(lig_resi)}
        view.setStyle(ligand_style_selector, {"stick": {"colorscheme": f"{ligand_color}Carbon"}})

    view.zoomTo()
    return view._make_html()


def build_filter_row(index, numeric_columns, row_state=None):
    """Create a single filter-row control block."""
    row_state = row_state or {}
    options = [{"label": col, "value": col} for col in numeric_columns]
    return html.Div(
        [
            dcc.Dropdown(
                id={"type": "filter-column", "index": index},
                options=options,
                value=row_state.get("column"),
                placeholder="Column",
                clearable=True,
                style={"minWidth": "220px"},
            ),
            dcc.Input(
                id={"type": "filter-min", "index": index},
                type="number",
                value=row_state.get("min"),
                placeholder="Min",
                style={"width": "120px"},
            ),
            dcc.Input(
                id={"type": "filter-max", "index": index},
                type="number",
                value=row_state.get("max"),
                placeholder="Max",
                style={"width": "120px"},
            ),
            html.Button("Remove", id={"type": "filter-remove", "index": index}, n_clicks=0),
        ],
        style={"display": "flex", "gap": "8px", "alignItems": "center", "marginBottom": "8px"},
    )


def build_uniprot_url(accession):
    """Return UniProt entry URL for a valid accession string."""
    if not isinstance(accession, str):
        return None
    accn = accession.strip()
    if not accn:
        return None
    return f"https://www.uniprot.org/uniprotkb/{accn}/entry"


def ensure_output_dir(dir_path):
    """Create output directory when missing and return normalized path."""
    normalized = os.path.abspath(os.path.expanduser((dir_path or "").strip()))
    if not normalized:
        raise ValueError("Please provide a destination directory path.")
    os.makedirs(normalized, exist_ok=True)
    if not os.path.isdir(normalized):
        raise OSError(f"Failed to create or access destination directory: {normalized}")
    return normalized


def normalize_inspection_value(value):
    """Normalize inspection labels to -1/0/1 integer values."""
    try:
        parsed = int(float(value))
    except (TypeError, ValueError):
        return 0
    if parsed > 0:
        return 1
    if parsed < 0:
        return -1
    return 0


def inspection_label(value):
    """Return a friendly inspection label."""
    normalized = normalize_inspection_value(value)
    if normalized == 1:
        return "Positive"
    if normalized == -1:
        return "Negative"
    return "Pending"


def _row_to_comparable_series(row_series, columns):
    """Convert row values into a stable comparable representation."""
    aligned = row_series.reindex(columns)
    return aligned.map(lambda val: "__NaN__" if pd.isna(val) else str(val))


def append_unique_row_to_csv(row_series, csv_path, columns):
    """
    Append selected row if not present already.

    Returns tuple: (updated_dataframe, was_appended).
    """
    row_df = pd.DataFrame([row_series.reindex(columns)])
    row_df = row_df.reindex(columns=columns)

    if not os.path.exists(csv_path):
        row_df.to_csv(csv_path, index=False)
        return row_df, True

    existing_df = pd.read_csv(csv_path)
    existing_df = existing_df.reindex(columns=columns)

    candidate = _row_to_comparable_series(row_series, columns)
    if not existing_df.empty:
        existing_cmp = existing_df.apply(lambda row: _row_to_comparable_series(row, columns), axis=1)
        duplicate_mask = existing_cmp.eq(candidate, axis=1).all(axis=1)
        if bool(duplicate_mask.any()):
            return existing_df, False

    updated_df = pd.concat([existing_df, row_df], ignore_index=True)
    updated_df.to_csv(csv_path, index=False)
    return updated_df, True


def render_pml_from_csv(df_saved_rows, pml_path):
    """Rewrite visualization.pml from all currently saved rows."""
    lines = [f"remote_loadpdb {target_path}"]
    for _, row in df_saved_rows.iterrows():
        matched_protein = str(row.get(matched_protein_colname, "")).strip()
        flattened_transform = str(row.get(flattened_transform_colname, "")).strip()
        if not matched_protein or not flattened_transform:
            continue
        escaped_transform = flattened_transform.replace('"', '\\"')
        lines.append(f"remote_loadAF2 {matched_protein}")
        lines.append(f'apply_transform {matched_protein}, "{escaped_transform}"')
        lines.append(f"coloraf {matched_protein}")
        matched_iface_resi = row.get(matched_iface_resi_colname, "")
        if pd.notna(matched_iface_resi):
            iface_resi_list = [res.strip() for res in str(matched_iface_resi).split(",") if res.strip()]
            if iface_resi_list:
                resi_selector = "+".join(iface_resi_list)
                lines.append(f'cmd.show("sticks", "{matched_protein} and resi {resi_selector}")')
    lines.append('cmd.show("lines"     ,"all")')
    lines.append('cmd.remove("(all) and hydro")')
    lines.append('util.cnc("all",_self=cmd)')
    with open(pml_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


# ------ Main App ------
df_results = pd.read_csv(results_path).reset_index(drop=True)
if INSPECTION_COL not in df_results.columns:
    df_results[INSPECTION_COL] = 0
df_results[INSPECTION_COL] = df_results[INSPECTION_COL].apply(normalize_inspection_value).astype(int)
all_columns = df_results.columns.tolist()
numeric_columns = [col for col in all_columns if pd.api.types.is_numeric_dtype(df_results[col])]

if not numeric_columns:
    raise ValueError("No numeric columns available for plotting/filtering.")

default_x = DEFAULT_X_AXIS if DEFAULT_X_AXIS in all_columns else numeric_columns[0]
default_y = DEFAULT_Y_AXIS if DEFAULT_Y_AXIS in all_columns else numeric_columns[min(1, len(numeric_columns) - 1)]
default_color = DEFAULT_COLOR_BY if DEFAULT_COLOR_BY in all_columns else numeric_columns[0]

default_filter_rows, default_dataset_mode = load_default_filter_state(
    default_filter_path, set(numeric_columns)
)
default_filters_json_path = (
    default_filter_path.strip()
    if default_filter_path
    else results_path.replace(".csv", "_filters.json")
)

parser = PDBParser()
target_struct = parser.get_structure("target", target_path)

app = Dash(__name__)

app.layout = html.Div(
    [
        dcc.Store(id="filter-store", data=default_filter_rows),
        dcc.Store(id="selected-row-store", data=None),
        dcc.Store(id="displayed-row-order-store", data=[]),
        html.H2("MaSIF-Neosurf Results Explorer (Dash)"),
        html.Div(
            [
                html.Div(
                    [
                        html.H4("Plot"),
                        dcc.Graph(
                            id="scatter-plot",
                            clear_on_unhover=True,
                            style={"height": "600px", "width": "100%"},
                            config={"responsive": True},
                        ),
                                                html.Div(
                            [
                                html.Label("Select point by row"),
                                dcc.Dropdown(
                                    id="row-selector-dropdown",
                                    options=[],
                                    value=None,
                                    placeholder="Select a row...",
                                    clearable=True,
                                ),
                            ],
                            style={"marginBottom": "12px"},
                        ),
                    ],
                    style={"padding": "12px", "border": "1px solid #ddd", "borderRadius": "8px"},
                ),
                html.Div(
                    [
                        html.H4("3D Viewer"),
                        html.Div(id="viewer-container", children="Click a point to render structure."),
                        html.Div(
                            [
                                html.Div(id="inspection-progress", children="No point selected."),
                                html.Div(
                                    [
                                        html.Button("Previous", id="inspection-prev-btn", n_clicks=0, accessKey="a", title="Shortcut: Alt+Shift+A"),
                                        html.Button("Next", id="inspection-next-btn", n_clicks=0, accessKey="d", title="Shortcut: Alt+Shift+D"),
                                    ],
                                    style={"display": "flex", "gap": "8px", "marginTop": "8px"},
                                ),
                                html.Div(
                                    [
                                        html.Button("Positive", id="inspection-positive-btn", n_clicks=0, accessKey="1", title="Shortcut: Alt+Shift+1"),
                                        html.Button("Pending", id="inspection-pending-btn", n_clicks=0, accessKey="2", title="Shortcut: Alt+Shift+2"),
                                        html.Button("Negative", id="inspection-negative-btn", n_clicks=0, accessKey="3", title="Shortcut: Alt+Shift+3"),
                                    ],
                                    style={"display": "flex", "gap": "8px", "marginTop": "8px"},
                                ),
                                html.Div(
                                    "Keyboard shortcuts use browser access keys (often Alt+Shift+key on Linux).",
                                    style={"marginTop": "6px", "fontSize": "12px", "color": "#555"},
                                ),
                            ],
                            style={"marginTop": "8px", "marginBottom": "12px"},
                        ),
                        html.H4("Save Visualization", style={"marginTop": "12px"}),
                        html.Label("Visualization output directory"),
                        dcc.Input(
                            id="viewer-save-dir",
                            type="text",
                            value=results_path.replace(".csv", "_visualization"),
                            style={"width": "100%", "marginBottom": "8px"},
                        ),
                        html.Button("Save visualization", id="save-visualization-btn", n_clicks=0),
                        html.Div(id="viewer-save-status-msg", style={"marginTop": "8px"}),
                    ],
                    style={"padding": "12px", "border": "1px solid #ddd", "borderRadius": "8px", "minHeight": "520px"},
                ),
                html.Div(
                    [
                        html.H4("Controls & Filters"),
                        html.Div(
                            [
                                html.Label("X axis"),
                                dcc.Dropdown(
                                    id="x-axis-dropdown",
                                    options=[{"label": col, "value": col} for col in all_columns],
                                    value=default_x,
                                    clearable=False,
                                ),
                            ],
                            style={"marginBottom": "8px"},
                        ),
                        html.Div(
                            [
                                html.Label("Y axis"),
                                dcc.Dropdown(
                                    id="y-axis-dropdown",
                                    options=[{"label": col, "value": col} for col in all_columns],
                                    value=default_y,
                                    clearable=False,
                                ),
                            ],
                            style={"marginBottom": "8px"},
                        ),
                        html.Div(
                            [
                                html.Label("Color by"),
                                dcc.Dropdown(
                                    id="color-by-dropdown",
                                    options=[{"label": col, "value": col} for col in all_columns],
                                    value=default_color,
                                    clearable=False,
                                ),
                            ],
                            style={"marginBottom": "12px"},
                        ),
                        html.H4("Filter Builder"),
                        html.Button("Add filter", id="add-filter-btn", n_clicks=0),
                        html.Div(id="filter-rows-container", style={"marginTop": "10px"}),
                        html.Label("Dataset mode", style={"marginTop": "8px", "display": "block"}),
                        dcc.RadioItems(
                            id="dataset-mode",
                            options=[
                                {"label": "Unfiltered", "value": "unfiltered"},
                                {"label": "Filtered", "value": "filtered"},
                            ],
                            value=default_dataset_mode,
                            inline=True,
                        ),
                        html.Div(id="dataset-shape-info", style={"marginTop": "6px"}),
                        html.H4("Save / Export", style={"marginTop": "14px"}),
                        html.Label("Filtering criteria path (.json)"),
                        dcc.Input(
                            id="filters-save-path",
                            type="text",
                            value=default_filters_json_path,
                            style={"width": "100%", "marginBottom": "8px"},
                        ),
                        html.Div(
                            [
                                html.Button("Save filters", id="save-filters-btn", n_clicks=0),
                            ],
                            style={"display": "flex", "gap": "8px", "marginBottom": "8px"},
                        ),
                        html.Label("Import filtering criteria path (.json)"),
                        dcc.Input(
                            id="filters-import-path",
                            type="text",
                            value=default_filters_json_path,
                            style={"width": "100%", "marginBottom": "8px"},
                        ),
                        html.Div(
                            [
                                html.Button("Import filters", id="import-filters-btn", n_clicks=0),
                            ],
                            style={"display": "flex", "gap": "8px", "marginBottom": "8px"},
                        ),
                        html.Label("Filtered dataframe path (.csv)"),
                        dcc.Input(
                            id="filtered-csv-save-path",
                            type="text",
                            value=results_path.replace(".csv", "_filtered.csv"),
                            style={"width": "100%", "marginBottom": "8px"},
                        ),
                        html.Div(
                            [
                                html.Button("Save filtered CSV", id="save-filtered-csv-btn", n_clicks=0),
                            ],
                            style={"display": "flex", "gap": "8px", "marginBottom": "8px"},
                        ),

                        html.Div(id="save-status-msg"),
                    ],
                    style={"padding": "12px", "border": "1px solid #ddd", "borderRadius": "8px", "minHeight": "300px"},
                ),
                html.Div(
                    [html.H4("UniProt Viewer"), html.Div(id="uniprot-container", children="No UniProt accession selected.")],
                    style={"padding": "12px", "border": "1px solid #ddd", "borderRadius": "8px", "minHeight": "520px"},
                ),
                html.Div(
                    [html.H4("Selected Row"), html.Div(id="row-summary", children="No point selected.")],
                    style={"padding": "12px", "border": "1px solid #ddd", "borderRadius": "8px", "minHeight": "300px"},
                ),
            ],
            style={"display": "grid", "gridTemplateColumns": "minmax(0, 1fr) minmax(0, 1fr)", "gap": "16px", "alignItems": "stretch"},
        ),

    ],
    style={"margin": "0 auto", "padding": "16px"},
)


@app.callback(
    Output("filter-store", "data"),
    Input("add-filter-btn", "n_clicks"),
    Input({"type": "filter-remove", "index": ALL}, "n_clicks"),
    Input({"type": "filter-column", "index": ALL}, "value"),
    Input({"type": "filter-min", "index": ALL}, "value"),
    Input({"type": "filter-max", "index": ALL}, "value"),
    State("filter-store", "data"),
    prevent_initial_call=True,
)
def update_filter_store(add_clicks, remove_clicks, columns, mins, maxs, filter_store):
    """Maintain normalized filter row state."""
    current_store = copy.deepcopy(filter_store or [])
    ctx = callback_context
    if not ctx.triggered:
        return current_store

    trigger = ctx.triggered[0]["prop_id"].split(".")[0]
    if trigger == "add-filter-btn":
        current_store.append({"id": str(uuid.uuid4()), "column": None, "min": None, "max": None})
        return current_store

    if trigger.startswith("{"):
        trigger_id = json.loads(trigger)
        if trigger_id.get("type") == "filter-remove":
            target_idx = trigger_id.get("index")
            return [row for i, row in enumerate(current_store) if i != target_idx]

    for i, row in enumerate(current_store):
        row["column"] = columns[i] if i < len(columns) else row.get("column")
        row["min"] = mins[i] if i < len(mins) else row.get("min")
        row["max"] = maxs[i] if i < len(maxs) else row.get("max")
    return current_store


@app.callback(
    Output("filter-rows-container", "children"),
    Input("filter-store", "data"),
)
def render_filter_rows(filter_store):
    """Render dynamic filter controls from store data."""
    rows = filter_store or []
    if not rows:
        return html.Div("No filters yet. Click 'Add filter'.")
    return [build_filter_row(i, numeric_columns, row_state=row) for i, row in enumerate(rows)]


@app.callback(
    Output("scatter-plot", "figure"),
    Output("dataset-shape-info", "children"),
    Input("x-axis-dropdown", "value"),
    Input("y-axis-dropdown", "value"),
    Input("color-by-dropdown", "value"),
    Input("dataset-mode", "value"),
    Input("filter-store", "data"),
    Input("selected-row-store", "data"),
)
def update_scatter_plot(x_axis, y_axis, color_by, dataset_mode, filter_store, selected_row_idx):
    """Render scatter plot from either full or filtered dataset."""
    filtering_criteria = build_filtering_criteria(filter_store)
    df_filtered = apply_filtering_criteria(df_results, filtering_criteria)
    df_display = (df_results if dataset_mode == "unfiltered" else df_filtered).copy()
    df_display["_row_idx"] = df_display.index

    fig = px.scatter(
        df_display,
        x=x_axis,
        y=y_axis,
        color=color_by,
        hover_data=[matched_protein_colname] if matched_protein_colname in df_display.columns else None,
        custom_data=["_row_idx"],
        color_continuous_scale="viridis",
        height=600,
    )

    if selected_row_idx is not None and selected_row_idx in df_display.index:
        selected_row = df_display.loc[selected_row_idx]
        fig.add_trace(
            go.Scatter(
                x=[selected_row[x_axis]],
                y=[selected_row[y_axis]],
                mode="markers",
                marker={
                    "symbol": "x",
                    "size": 16,
                    "color": "red",
                    "line": {"width": 2, "color": "red"},
                },
                name="Selected",
                hoverinfo="skip",
                showlegend=False,
            )
        )

    info = f"Rows shown: {df_display.shape[0]} / {df_results.shape[0]}"
    return fig, info


@app.callback(
    Output("row-selector-dropdown", "options"),
    Input("x-axis-dropdown", "value"),
    Input("y-axis-dropdown", "value"),
    Input("color-by-dropdown", "value"),
    Input("dataset-mode", "value"),
    Input("filter-store", "data"),
)
def update_row_selector_options(x_axis, y_axis, color_by, dataset_mode, filter_store):
    """Build dropdown options from displayed rows and keep current selection valid."""
    filtering_criteria = build_filtering_criteria(filter_store)
    df_filtered = apply_filtering_criteria(df_results, filtering_criteria)
    df_display = df_results if dataset_mode == "unfiltered" else df_filtered

    def safe_value(val):
        if pd.isna(val):
            return "-"
        return val

    options = []
    for row_idx, row in df_display.iterrows():
        matched_value = safe_value(row.get(matched_protein_colname, "-"))
        x_val = safe_value(row.get(x_axis, "-"))
        y_val = safe_value(row.get(y_axis, "-"))
        c_val = safe_value(row.get(color_by, "-"))
        i_val = inspection_label(row.get(INSPECTION_COL, 0))
        label = (
            f"idx={row_idx} | matched_protein={matched_value} | "
            f"{x_axis}={x_val} | {y_axis}={y_val} | {color_by}={c_val} | inspection={i_val}"
        )
        options.append({"label": label, "value": int(row_idx)})

    return options


@app.callback(
    Output("displayed-row-order-store", "data"),
    Input("dataset-mode", "value"),
    Input("filter-store", "data"),
)
def update_displayed_row_order(dataset_mode, filter_store):
    """Track row indices currently displayed in the scatter plot."""
    filtering_criteria = build_filtering_criteria(filter_store)
    df_filtered = apply_filtering_criteria(df_results, filtering_criteria)
    df_display = df_results if dataset_mode == "unfiltered" else df_filtered
    return [int(idx) for idx in df_display.index.tolist()]


@app.callback(
    Output("viewer-container", "children"),
    Output("row-summary", "children"),
    Output("uniprot-container", "children"),
    Output("selected-row-store", "data"),
    Output("inspection-progress", "children"),
    Input("scatter-plot", "clickData"),
    Input("row-selector-dropdown", "value"),
    Input("inspection-prev-btn", "n_clicks"),
    Input("inspection-next-btn", "n_clicks"),
    Input("inspection-positive-btn", "n_clicks"),
    Input("inspection-pending-btn", "n_clicks"),
    Input("inspection-negative-btn", "n_clicks"),
    Input("x-axis-dropdown", "value"),
    Input("y-axis-dropdown", "value"),
    Input("color-by-dropdown", "value"),
    State("displayed-row-order-store", "data"),
    State("selected-row-store", "data"),
)
def update_selection(
    click_data,
    dropdown_row_idx,
    prev_clicks,
    next_clicks,
    positive_clicks,
    pending_clicks,
    negative_clicks,
    x_axis,
    y_axis,
    color_by,
    displayed_row_order,
    selected_row_store,
):
    """Update 3D viewer + row summary from clicked point or dropdown selection."""
    ctx = callback_context
    triggered = ctx.triggered[0]["prop_id"].split(".")[0] if ctx.triggered else ""
    displayed_row_order = displayed_row_order or []

    row_idx = None
    if triggered == "scatter-plot" and click_data and click_data.get("points"):
        customdata = click_data["points"][0].get("customdata")
        row_idx = customdata[0] if isinstance(customdata, (list, tuple)) and customdata else customdata
    elif triggered == "row-selector-dropdown":
        row_idx = dropdown_row_idx
    elif triggered in {"inspection-prev-btn", "inspection-next-btn"}:
        row_idx = selected_row_store
        if displayed_row_order:
            try:
                current_pos = displayed_row_order.index(int(row_idx))
            except (TypeError, ValueError):
                current_pos = 0
            except Exception:
                current_pos = 0
            if triggered == "inspection-prev-btn":
                next_pos = max(0, current_pos - 1)
            else:
                next_pos = min(len(displayed_row_order) - 1, current_pos + 1)
            row_idx = displayed_row_order[next_pos]
    elif triggered in {"inspection-positive-btn", "inspection-pending-btn", "inspection-negative-btn"}:
        row_idx = selected_row_store
        if row_idx is None and displayed_row_order:
            row_idx = displayed_row_order[0]
        if row_idx is not None:
            try:
                row_idx = int(row_idx)
                if row_idx in df_results.index:
                    if triggered == "inspection-positive-btn":
                        df_results.at[row_idx, INSPECTION_COL] = 1
                    elif triggered == "inspection-negative-btn":
                        df_results.at[row_idx, INSPECTION_COL] = -1
                    else:
                        df_results.at[row_idx, INSPECTION_COL] = 0
                    if displayed_row_order:
                        try:
                            current_pos = displayed_row_order.index(row_idx)
                        except ValueError:
                            current_pos = -1
                        if current_pos >= 0:
                            row_idx = displayed_row_order[min(len(displayed_row_order) - 1, current_pos + 1)]
            except (TypeError, ValueError):
                row_idx = None
    else:
        row_idx = selected_row_store

    try:
        row_idx = int(row_idx)
    except (TypeError, ValueError):
        row_idx = None

    if row_idx is None:
        return (
            "Click a point or select a row from dropdown to render structure.",
            "No point selected.",
            "No UniProt accession selected.",
            None,
            "No point selected.",
        )

    if row_idx not in df_results.index:
        return (
            "Selected point has no valid row index.",
            "No row summary available.",
            "No UniProt accession selected.",
            None,
            "Selected row index is not valid for current dataframe.",
        )

    row = df_results.loc[row_idx]
    first_cols = list(df_results.columns[:20])
    special = [col for col in [x_axis, y_axis, color_by, INSPECTION_COL] if col not in first_cols]
    #hidden_summary_cols = {"target_path", "matched_protein_path"}
    hidden_summary_cols = {}
    selected_cols = [
        col for col in first_cols + special
        if col in df_results.columns and col not in hidden_summary_cols
    ]

    try:
        viewer_html = view_predicted_binder_html(
            row,
            matched_protein_colname,
            flattened_transform_colname,
            target_struct,
            domainome_db_dir,
            parser,
            lig_chain=lig_chain,
            lig_name=lig_name,
            lig_resi=lig_resnum,
            target_iface_resi_colname=target_iface_resi_colname,
            matched_iface_resi_colname=matched_iface_resi_colname,
        )
        viewer = html.Iframe(srcDoc=viewer_html, style={"width": "100%", "height": "600px", "border": "none"})
    except Exception as exc:
        viewer = html.Pre(f"Unable to render structure: {exc}")

    uniprot_url = build_uniprot_url(row.get(UNIPROT_ACCN_COL))
    if uniprot_url is None:
        uniprot_view = html.Div(f"No valid `{UNIPROT_ACCN_COL}` value for selected row.")
    else:
        uniprot_view = html.Div(
            [
                html.Div(
                    [
                        "Open in new tab: ",
                        html.A(uniprot_url, href=uniprot_url, target="_blank", rel="noopener noreferrer"),
                    ],
                    style={"marginBottom": "8px"},
                ),
                html.Iframe(src=uniprot_url, style={"width": "100%", "height": "600px", "border": "none"}),
            ]
        )

    summary_rows = [html.Tr([html.Th("index"), html.Td(str(row_idx))])]
    summary_rows.extend([html.Tr([html.Th(col), html.Td(str(row[col]))]) for col in selected_cols])
    summary_table = html.Table(summary_rows, style={"width": "100%", "borderCollapse": "collapse"})
    if displayed_row_order and row_idx in displayed_row_order:
        row_pos = displayed_row_order.index(row_idx) + 1
        total = len(displayed_row_order)
    else:
        row_pos = 1
        total = max(1, len(displayed_row_order))
    matched_value = row.get(matched_protein_colname, "-")
    progress_text = (
        f"Row {row_pos} / {total} (global idx={row_idx}) | "
        f"inspection={inspection_label(row.get(INSPECTION_COL, 0))} | "
        f"{matched_protein_colname}={matched_value}"
    )
    return viewer, summary_table, uniprot_view, row_idx, progress_text


@app.callback(
    Output("viewer-save-status-msg", "children"),
    Input("save-visualization-btn", "n_clicks"),
    State("viewer-save-dir", "value"),
    State("selected-row-store", "data"),
    State("row-selector-dropdown", "value"),
    prevent_initial_call=True,
)
def save_visualization_assets(n_clicks, viewer_save_dir, selected_row_idx, dropdown_row_idx):
    """Persist selected row to visualized_rows.csv and regenerate visualization.pml."""
    if not n_clicks:
        return ""
    row_candidate = selected_row_idx if selected_row_idx is not None else dropdown_row_idx
    if row_candidate is None:
        return "Error: No row selected. Click a point or choose a row first."

    try:
        row_idx = int(row_candidate)
    except (TypeError, ValueError):
        return "Error: Selected row index is invalid."

    if row_idx not in df_results.index:
        return f"Error: Selected row index {row_idx} not found in dataset."

    row = df_results.loc[row_idx]
    columns = list(df_results.columns)

    try:
        output_dir = ensure_output_dir(viewer_save_dir)
    except (ValueError, OSError) as exc:
        return f"Error: {exc}"

    csv_path = os.path.join(output_dir, VIS_ROWS_FILENAME)
    pml_path = os.path.join(output_dir, VIS_PML_FILENAME)

    try:
        saved_df, was_appended = append_unique_row_to_csv(row, csv_path, columns)
        render_pml_from_csv(saved_df, pml_path)
    except Exception as exc:  # noqa: BLE001
        return f"Error while saving visualization assets: {exc}"

    action_msg = "appended new row" if was_appended else "row already exists (no duplicate appended)"
    return (
        f"[click {n_clicks}] Saved visualization assets in: {output_dir} | "
        f"{action_msg} | rows in {VIS_ROWS_FILENAME}: {saved_df.shape[0]} | "
        f"rewrote {VIS_PML_FILENAME}"
    )


@app.callback(
    Output("save-status-msg", "children"),
    Output("filter-store", "data", allow_duplicate=True),
    Input("save-filters-btn", "n_clicks"),
    Input("import-filters-btn", "n_clicks"),
    Input("save-filtered-csv-btn", "n_clicks"),
    State("filter-store", "data"),
    State("filters-save-path", "value"),
    State("filters-import-path", "value"),
    State("filtered-csv-save-path", "value"),
    prevent_initial_call=True,
)
def save_filters_and_csv(
    save_filters_clicks,
    import_filters_clicks,
    save_csv_clicks,
    filter_store,
    filters_save_path,
    filters_import_path,
    filtered_csv_save_path,
):
    """Save filter criteria JSON or filtered CSV to user-provided file paths."""
    ctx = callback_context
    if not ctx.triggered:
        return "", no_update

    trigger = ctx.triggered[0]["prop_id"].split(".")[0]
    filtering_criteria = build_filtering_criteria(filter_store)

    if trigger == "save-filters-btn":
        output_path = (filters_save_path or "").strip()
        if not output_path:
            return "Error: Please provide a path for filtering criteria JSON.", no_update
        parent_dir = os.path.dirname(output_path) or "."
        if not os.path.isdir(parent_dir):
            return f"Error: Directory does not exist: {parent_dir}", no_update

        serializable_criteria = {
            col: {"min": min_val, "max": max_val}
            for col, (min_val, max_val) in filtering_criteria.items()
        }
        try:
            with open(output_path, "w", encoding="utf-8") as handle:
                json.dump(serializable_criteria, handle, indent=2)
            return f"Saved filtering criteria to: {output_path}", no_update
        except OSError as exc:
            return f"Error saving filtering criteria: {exc}", no_update

    if trigger == "import-filters-btn":
        import_path = (filters_import_path or "").strip()
        if not import_path:
            return "Error: Please provide a path for filtering criteria JSON to import.", no_update
        if not os.path.isfile(import_path):
            return f"Error: File does not exist: {import_path}", no_update
        try:
            with open(import_path, "r", encoding="utf-8") as handle:
                imported_payload = json.load(handle)
            imported_rows, ignored_columns = normalize_imported_filters(imported_payload, set(numeric_columns))
        except (OSError, json.JSONDecodeError, ValueError) as exc:
            return f"Error importing filtering criteria: {exc}", no_update

        if not imported_rows:
            ignored_msg = f" Ignored columns: {', '.join(ignored_columns)}." if ignored_columns else ""
            return f"No valid numeric filter columns found in imported file.{ignored_msg}", []

        ignored_msg = f" Ignored columns: {', '.join(ignored_columns)}." if ignored_columns else ""
        return (
            f"Imported {len(imported_rows)} filters from: {import_path}.{ignored_msg}",
            imported_rows,
        )

    if trigger == "save-filtered-csv-btn":
        output_path = (filtered_csv_save_path or "").strip()
        if not output_path:
            return "Error: Please provide a path for filtered CSV.", no_update
        parent_dir = os.path.dirname(output_path) or "."
        if not os.path.isdir(parent_dir):
            return f"Error: Directory does not exist: {parent_dir}", no_update

        filtered_df = apply_filtering_criteria(df_results, filtering_criteria)
        try:
            filtered_df.to_csv(output_path, index=False)
            return f"Saved filtered CSV to: {output_path} (rows: {filtered_df.shape[0]})", no_update
        except OSError as exc:
            return f"Error saving filtered CSV: {exc}", no_update

    return "", no_update


if __name__ == "__main__":
    app.run(debug=True)
