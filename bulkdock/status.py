import mrich
import subprocess
from rich.table import Table
from rich.panel import Panel
from richqueue.slurm import combined_df, get_user
from richqueue.table import color_by_state, COLUMNS
from datetime import timedelta
from richqueue.tools import human_timedelta  # , human_timedelta_to_seconds


def status(user: str = None):

    user = user or get_user()

    df = combined_df(user=user)

    # filter dataframe
    df = df[df["name"].str.startswith("BulkDock")]
    df = df[df["job_state"].isin(["RUNNING", "PENDING"])]
    df = df[["name", "job_id", "run_time", "standard_output", "job_state"]]
    df = df.sort_values(by="job_id")

    # create the table
    table = Table(title="Active BulkDock Jobs", box=None, header_style="")

    table.add_column(**COLUMNS["job_id"])
    table.add_column("[cyan3 underline]Command", style="cyan3")
    table.add_column("[cyan3 underline]Target", style="cyan3")
    table.add_column(**COLUMNS["name"])
    table.add_column(**COLUMNS["run_time"])
    table.add_column(**COLUMNS["job_state"])
    table.add_column(
        "[cornflower_blue underline]Attempted", justify="right", style="cornflower_blue"
    )
    table.add_column(
        "[bold underline]Progress", justify="right", style="bold cornflower_blue"
    )
    table.add_column(
        "[bold underline]Performance", justify="right", style="bold cornflower_blue"
    )
    table.add_column(
        "[bold underline]Locked", justify="right", style="bold cornflower_blue"
    )
    table.add_column(
        "[cornflower_blue underline]Remaining", justify="right", style="cornflower_blue"
    )

    for i, row in df.iterrows():

        values = []

        name = row["name"]

        command, target, name = name.split(":")
        command = command.removeprefix("BulkDock.")

        # calculate placement progress

        grep = [f'grep "Placement task " {row.standard_output} | tail -n 1']

        x = subprocess.run(
            grep, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        progress = x.stdout.decode()

        grep = [
            f'grep "SQLite Database was locked, retrying..." {row.standard_output} | wc -l'
        ]

        x = subprocess.run(
            grep, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        locked = x.stdout.decode()

        run_seconds = human_timedelta_to_seconds(row.run_time)

        if command == "place":
            try:
                i, n = progress.split("Placement task ")[-1].split(" ")[0].split("/")

                i = int(i)
                n = int(n)

                fraction = i / n

                progress = color_by_fraction(fraction)

                # calculate performance

                performance = color_by_performance(run_seconds / i)

                # calculate remaining estimate

                remaining = (n - i) * (run_seconds / i)
                remaining = human_timedelta(timedelta(seconds=remaining))

            except ValueError:
                progress = color_by_fraction(0)
                performance = ""
                remaining = ""
                i = 0

        else:

            progress = ""
            performance = ""
            remaining = ""
            i = ""

        # calculate lock fraction
        if run_seconds:
            locked = f"{int(locked) / run_seconds*100:.1f} %"
        else:
            locked = ""

        values.append(str(row.job_id))
        values.append(command)
        values.append(target)
        values.append(name)
        values.append(str(row.run_time))
        values.append(color_by_state(row.job_state))
        values.append(str(i))
        values.append(progress)
        values.append(performance)
        values.append(locked)
        values.append(str(remaining))

        table.add_row(*values)

    mrich.print(Panel(table, expand=False))


def human_timedelta_to_seconds(s):

    values = s.split()

    seconds = 0

    for value in values:
        if value.endswith("s"):
            seconds += int(value[:-1])
        elif value.endswith("m"):
            seconds += 60 * int(value[:-1])
        elif value.endswith("h"):
            seconds += 60 * 60 * int(value[:-1])
        elif value.endswith("d"):
            seconds += 24 * 60 * 60 * int(value[:-1])

    return seconds


def color_by_fraction(fraction):

    if fraction > 0.75:
        color = "bright_green"
    elif fraction > 0.50:
        color = "bright_yellow"
    elif fraction > 0.25:
        color = "dark_orange"
    else:
        color = "red"

    return f"[{color}]{fraction*100:.1f} %"


def color_by_fraction_inverse(fraction):

    if fraction > 0.75:
        color = "dark_orange"
    elif fraction > 0.50:
        color = "bright_yellow"
    elif fraction > 0.25:
        color = "bright_green"
    else:
        color = "red"

    return f"[{color}]{fraction*100:.1f} %"


def color_by_performance(performance):

    if performance < 10:
        color = "bright_green"
    elif performance < 20:
        color = "bright_yellow"
    elif performance < 60:
        color = "dark_orange"
    else:
        color = "red"

    return f"[{color}]{performance:.2f} s/it"
