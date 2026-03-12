import mrich
import typer
from .bulkdock import BulkDock
from typing import Optional
from typing_extensions import Annotated
from .config import VARIABLES

HELP = """
💪 BulkDock: Manage batches of Fragmenstein restrained protein-ligand docking jobs
"""

app = typer.Typer(help=HELP, no_args_is_help=True)
engine = BulkDock()


@app.command()
def status(
    user: str = None,
):
    """Show the status of running BulkDock jobs"""

    from .status import status

    status(user=user)


@app.command()
def export(
    target: str,
    tag: str,
    best_by_compound: bool = False,
    metadata: bool = False,
    subsites: bool = True,
    tags: bool = True,
    generate_pdbs: bool = False,
    submitter_name: str = None,
    submitter_institution: str = None,
    submitter_email: str = None,
    ref_url: str = None,
    max_energy_score: float | None = 0.0,
    max_distance_score: float | None = 2.0,
    require_outcome: str | None = "acceptable",
    posebusters: bool = True,
    output: str | None = None,
    debug: bool = False,
):
    """Export poses from a successful output into a Fragalysis-ready format"""

    pose_filter_methods = []
    if posebusters:
        pose_filter_methods.append("posebusters")

    if require_outcome in ["None", "False"]:
        require_outcome = None

    engine.export(
        target=target,
        tag=tag,
        best_by_compound=best_by_compound,
        generate_pdbs=generate_pdbs,
        submitter_name=submitter_name,
        submitter_institution=submitter_institution,
        submitter_email=submitter_email,
        ref_url=ref_url,
        max_energy_score=max_energy_score,
        max_distance_score=max_distance_score,
        require_outcome=require_outcome,
        pose_filter_methods=pose_filter_methods,
        output=output,
        metadata=metadata,
        tags=tags,
        subsites=subsites,
        debug=debug,
    )


@app.command()
def place(
    target: str,
    file: str,
    split: int = 1_000,
    stagger: int = 0.25,
    dependency: int = 0,
    reference: Annotated[
        str,
        typer.Option(
            help="Name of reference pose, if none is specified will ensemble dock against inspirations"
        ),
    ] = "",
):
    """Start a placement job.

    Input file must be a CSV with the first column containing the SMILES and all subsequent columns containing observation shortcodes for inspiration hits

    """

    # :param target: to place against
    # :param name: or path to input file (must be in configured INPUTS directory)
    # :param split: split the input file into batches of this size

    engine.submit_placement_jobs(
        target,
        file,
        split=split,
        stagger=stagger,
        dependency=dependency,
        reference=reference,
    )


@app.command()
def configure(
    variable: Annotated[
        str, typer.Argument(help=f"variable to configure, options are: {VARIABLES}")
    ],
    value: str,
):
    """Configure"""

    mrich.h2("BulkDock.configure")
    mrich.var("variable", variable)
    mrich.var("value", value)
    assert variable in VARIABLES
    engine.set_config_value(variable, value)

@app.command()
def requeue(
    job_id: list[int],
    stagger: float = 0,
):
    """Requeue a bulkdock placement job with given ID(s)"""

    import time

    for i in job_id:
        if stagger and i > 0:
            with mrich.clock("Staggering job submission..."):
                time.sleep(stagger)
        engine.requeue_placement_job(job_id=i)

@app.command()
def create_directories():
    """Create directories as configured"""
    engine.create_directories()


@app.command()
def extract(target: str):
    """Extract target zip file"""
    engine.extract_target(target)


@app.command()
def setup(target: str):
    """Setup HIPPO Database for a target"""
    engine.setup_hippo(target)


def main():
    app()


if __name__ == "__main__":
    main()
