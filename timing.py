import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@click.command()
@click.argument("ref_file")
@click.argument("cur_file")
def cli(ref_file, cur_file):
    ref = pd.read_csv(
        ref_file,
        delim_whitespace=True,
        names=["routine", "time", "ncalls"],
        index_col=0,
    )
    cur = pd.read_csv(
        cur_file,
        delim_whitespace=True,
        names=["routine", "time", "ncalls"],
        index_col=0,
    )

    cur = cur.sort_values(by=["time"], ascending=False)
    ref = ref.sort_values(by=["time"], ascending=False)

    cur = cur.rename(index=lambda s: s.replace("_lr", ""))
    cur = cur.rename(index=lambda s: s.replace("_acc", ""))

    selection = set(
        [
            "radiation_two_stream:calc_two_stream_gammas_lw",
            "radiation_two_stream:calc_reflectance_transmittance_lw",
            "radiation_adding_ica_lw:adding_ica_lw",
            "radiation_two_stream:calc_no_scattering_transmittance_lw",
            "radiation_adding_ica_lw:calc_fluxes_no_scattering_lw",
            "radiation_cloud_generator:cloud_generator",
            "radiation_adding_ica_lw:fast_adding_ica_lw",
            "radiation_lw_derivatives:modify_lw_derivatives_ica",
            "radiation_lw_derivatives:calc_lw_derivatives_ica",
            "radiation_mcica_lw:solver_mcica_lw",
        ]
    )

    selection = selection.intersection(set(cur.index)).intersection(set(ref.index))

    a = cur.loc[selection, "time"]
    b = ref.loc[selection, "time"]

    print("\ncurrent time")
    print("-------------------------------------")
    print(a)

    print("\nreference time")
    print("-------------------------------------")
    print(b)


if __name__ == "__main__":
    cli()
