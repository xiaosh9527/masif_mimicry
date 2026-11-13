from __future__ import annotations

import csv
from pathlib import Path

import streamlit as st
from pydantic import ValidationError

from masif_mimicry.config import JobConfig, dump_job
from masif_mimicry.pipeline import prepare_run
from masif_mimicry.results import find_result_csvs

st.set_page_config(page_title="MaSIF-mimicry", page_icon="🧬", layout="wide")
st.title("MaSIF-mimicry")
st.caption("Configure, validate, and inspect a protein-surface mimicry screen.")

setup_tab, results_tab = st.tabs(["New screen", "Results"])

with setup_tab:
    with st.form("job"):
        st.subheader("Experiment")
        col1, col2 = st.columns(2)
        name = col1.text_input("Name", "znf692-vs-mtor")
        project_dir = col2.text_input("Experiment directory", str(Path.cwd() / "experiment"))

        st.subheader("Target")
        target_path = st.text_input("Target PDB path")
        col1, col2, col3 = st.columns(3)
        query_chain = col1.text_input("Query chain", "C")
        receptor_chain = col2.text_input("Receptor chain", "B")
        ligand = col3.text_input("Ligand", "B:Y70", help="CHAIN:RESNAME")
        col1, col2 = st.columns(2)
        residue = col1.number_input("Focus residue", min_value=1, value=423)
        atom = col2.text_input("Focus atom", "CA")

        st.subheader("Search")
        proteins = st.text_area(
            "Processed database proteins (one per line)", "P42345-F1-dom-01_A"
        )
        col1, col2 = st.columns(2)
        preset = col1.selectbox("Preset", ["balanced", "fast", "sensitive"])
        count_clashes = col2.checkbox("Count clashes", True)

        st.subheader("Runtime")
        backend = st.selectbox("Backend", ["slurm", "apptainer", "local"])
        suite_root = st.text_input("MaSIF suite root")
        image = st.text_input(
            "Apptainer image", disabled=backend == "local", help="Required for Slurm and Apptainer"
        )
        prepare = st.form_submit_button("Validate and prepare", type="primary")

    if prepare:
        try:
            directory = Path(project_dir).expanduser().resolve()
            raw = {
                "name": name,
                "target": {
                    "structure": Path(target_path).expanduser().resolve(),
                    "query_chain": query_chain,
                    "receptor_chain": receptor_chain or None,
                    "focus": {"residue": residue, "atom": atom},
                    "ligand": ligand or None,
                },
                "database": {
                    "name": "masif_human_proteome",
                    "proteins": [value.strip() for value in proteins.splitlines() if value.strip()],
                },
                "search": {"preset": preset, "count_clashes": count_clashes},
                "postprocessing": {"enabled": bool(ligand), "truncation_length": 86},
                "runtime": {
                    "backend": backend,
                    "suite_root": Path(suite_root).expanduser().resolve(),
                    "image": Path(image).expanduser().resolve() if image else None,
                    "submit": False,
                },
                "output_root": directory / "runs",
            }
            job = JobConfig.model_validate(raw)
            directory.mkdir(parents=True, exist_ok=True)
            job_path = directory / "job.yaml"
            dump_job(job, job_path)
            plan = prepare_run(job)
            st.success(f"Prepared {plan.run_dir}")
            st.code(plan.script_path.read_text(encoding="utf-8"), language="bash")
        except (OSError, ValueError, ValidationError) as exc:
            st.error(str(exc))

with results_tab:
    run_path_value = st.text_input("Run directory", key="results_path")
    if run_path_value:
        run_path = Path(run_path_value).expanduser().resolve()
        if not run_path.is_dir():
            st.warning("Run directory does not exist.")
        else:
            csv_paths = find_result_csvs(run_path)
            st.metric("Result files", len(csv_paths))
            for csv_path in csv_paths:
                st.subheader(str(csv_path.relative_to(run_path)))
                with csv_path.open(newline="", encoding="utf-8") as handle:
                    rows = list(csv.DictReader(handle))
                st.dataframe(rows, use_container_width=True)
