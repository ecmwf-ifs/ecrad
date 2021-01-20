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
        ]
    )

    selection = selection.intersection(set(cur.index)).intersection(set(ref.index))

    print(selection)

    a = cur.loc[selection, "time"].values
    b = ref.loc[selection, "time"].values

    diff = np.abs((a - b) / a) * 100

    fig, ax = plt.subplots(ncols=1, nrows=1)

    ax.barh(range(len(selection)), diff)

    plt.show()


if __name__ == "__main__":
    cli()
