import mrich
import typer
from typing_extensions import Annotated
from .bulkdock import BulkDock

HELP = """
💪 BulkDock: INTERNAL CLI ONLY!
"""

app = typer.Typer()
engine = BulkDock()

"""Inner CLI for batch jobs"""


@app.command()
def place(
    target: str,
    file: str,
    reference: Annotated[
        str,
        typer.Option(
            help="Name of reference pose, if none is specified will ensemble dock against inspirations"
        ),
    ] = "",
    debug: bool = False,
):
    """Run Bulkdock.place"""
    mrich.h3("bulkdock.batch.place")
    mrich.var("target", target)
    mrich.var("file", file)
    mrich.var("reference", reference)
    mrich.var("debug", debug)
    engine.place(target, file, reference=reference, debug=debug)


@app.command()
def combine(
    target: str, 
    csv_file: str,
    batch_size: Annotated[
        int,
        typer.Option(
            help="Select one of multiple batch sizes"
        ),
    ] = 0,
    skip_names: bool = False,
):
    """Combine split SDF outputs from placement jobs"""

    import subprocess
    from pandas import DataFrame
    from pathlib import Path
    from rich.table import Table
    from math import ceil
    from molparse.rdkit import sdf_combine

    mrich.h3("bulkdock.batch.combine")
    mrich.var("csv_file", csv_file)

    batch_size = batch_size or None

    csv_path = engine.get_infile_path(csv_file)
    mrich.var("csv_path", csv_path)

    commands = [f"grep -v smiles {str(csv_path.resolve())} | wc -l"]
    x = subprocess.run(commands, shell=True, stdout=subprocess.PIPE)
    num_compounds = int(x.stdout.decode())

    mrich.var("num_compounds", num_compounds)

    # name of the input file
    key = Path(csv_file).name

    assert key.endswith(".csv") or key.endswith(".sdf")
    key = key.removesuffix(".csv").removesuffix(".sdf")

    mrich.var("key", key)

    pattern = f"{key}*.sdf"

    files = list(engine.output_dir.glob(pattern))

    if not files:
        mrich.error(f"Did not find any files in {engine.output_dir}/{pattern}")
        raise FileNotFoundError(
            f"Did not find any files in {engine.output_dir}/{pattern}"
        )

    if len(files) > 1:

        df = []
        for i,file in enumerate(files):

            d = dict(key=key, file=file)

            file_name = file.name

            detail = file_name.removeprefix(key).removesuffix(".sdf")

            fields = [s for s in detail.split("_") if s]

            if not len(fields) == 3:
                mrich.warning(f"Ignoring {file_name=}")
                continue

            d["batch_size"] = int(fields[0].removeprefix("split"))
            d["batch_index"] = int(fields[1].removeprefix("batch"))
            d["job_id"] = int(fields[2])

            df.append(d)

        df = DataFrame(df)

        df = df.sort_values(by="batch_index")

        mrich.print(df.drop(columns=["file"]))

        if batch_size:
            df = df[df["batch_size"]==batch_size]

        batch_sizes = df["batch_size"].to_numpy()
        if not (batch_sizes[0] == batch_sizes).all():
            raise NotImplementedError("Not currently supporting multiple batch sizes. Re-run with batch-size option")

        files = []

        batch_size = batch_sizes[0]

        mrich.var("batch_size", batch_size)

        expected_batch_count = ceil(num_compounds / batch_size)

        mrich.var("expected_batch_count", expected_batch_count)

        if len(df) > expected_batch_count:
            mrich.warning("Too many batches")

        elif len(df) < expected_batch_count:
            mrich.warning("Missing batches")

        for i in range(expected_batch_count):

            subdf = df[df["batch_index"] == i]

            if len(subdf) == 0:
                mrich.error(f"Missing batch {i}")
                continue

            elif len(subdf) > 1:
                mrich.warning(f"Multiple batches w/ {i=}: {subdf}")
                row = subdf.iloc[0]
            else:
                row = subdf.iloc[0]

            files.append(row["file"])

    mrich.var("files", files)
    
    animal = engine.get_animal(target)

    for i,file in enumerate(files):

        mrich.h1(f"{i}/{len(files)}")

        mrich.bold(file)
        try:
            animal.load_sdf(
                target=target, 
                path=file, 
                inspiration_col="inspiration_ids", 
                reference_col="reference_id",
                compound_tags = [key],
                pose_tags = ["fragmenstein_placed", key],
                skip_names = False,
            )
        except AssertionError:
            mrich.error("Something wrong with", file)
            continue

    mrich.success("Done")

@app.command()
def collate(
    outname: str, target: str, json_path: str, tag: str = "Fragmenstein placed"
):

    import json
    from pathlib import Path

    mrich.var("outname", outname)
    mrich.var("target", target)
    mrich.var("json_path", json_path)
    mrich.var("tag", tag)

    animal = engine.get_animal(target)

    subdir = engine.get_scratch_subdir(f"{target}_inputs")

    json_path = subdir / Path(json_path).name

    assert json_path.exists(), "JSON Input does not exist"

    job_ids = json.load(open(json_path, "rt"))
    job_ids = set(job_ids)

    assert job_ids, "Null Job IDs"

    mrich.var("job_ids", job_ids)

    tagged_poses = animal.poses(tag=tag)

    mrich.var("tagged_poses", tagged_poses)

    pose_ids = set()
    for i, pose in mrich.track(enumerate(tagged_poses), prefix="Getting poses"):

        mrich.set_progress_field("progress", f"{i+1}/{len(tagged_poses)}")
        mrich.set_progress_field("found", len(pose_ids))

        job_id = int(Path(pose.path).parent.parent.name)

        if job_id not in job_ids:
            continue

        pose_ids.add(pose.id)

    poses = animal.poses[pose_ids]

    mrich.var("#poses", len(poses))

    outfile = engine.get_outfile_path(outname)
    poses.write_sdf(outfile, name_col="id")

    return outfile


def main():
    app()


if __name__ == "__main__":
    main()
