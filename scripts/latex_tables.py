import numpy as np
import os
import utils
import pandas as pd


def create_observation_comparison_table(name):
    modelled = utils.get_salinity_and_flux_at_measurement_cells(name)
    data = pd.read_csv("/home/superuser/objective_3/data/observations/observations.csv", index_col=0,
                       dtype={'point': str, 'concentration_average': float, 'discharge_min': float, 'discharge_max': float})
    table = pd.DataFrame(index=modelled.index, columns=["C [PSU] (modelled)", r"$q_z$ [m/day] (modelled)",
                                                        "C [PSU] (Hoffman)", "$q_z$ [m/day] (Harding)"])
    data.index = modelled.index
    for point in modelled.index:
        table.loc[point, "C [PSU] (modelled)"] =  f"{modelled.loc[point, 'conc_min']:.1f} - {modelled.loc[point, 'conc_max']:.1f}"
        table.loc[point, r"$q_z$ [m/day] (modelled)"] = f"{modelled.loc[point, 'discharge_min']:.1e}}}$ - {modelled.loc[point, 'discharge_max']:.1e}}}$"
        if point[0] == "m":
            table.loc[point, "C [PSU] (Hoffman)"] = f"{data.loc[point, "concentration_average"]:.1f}"
            table.loc[point, "$q_z$ [m/day] (Harding)"] = "--"
        else:
            table.loc[point, "C [PSU] (Hoffman)"] = "--"
            table.loc[point, "$q_z$ [m/day] (Harding)"] = f"{data.loc[point, 'discharge_min']:.1e}}}$ - {data.loc[point, 'discharge_max']:.1e}}}$"
    string = table.to_latex()
    string = str.replace(string, "e+0", r"$\times 10^{")
    string = str.replace(string, "e-0", r"$\times 10^{-")
    with open(f"/home/superuser/objective_3/outputs/{name}/table.txt", "w") as text_file:
        text_file.write(string)

def pockmark_discharges_by_scenario(set, names):
    table = pd.DataFrame(columns=["discharge", "conc"]*len(names), index=[i for i in range (1, 7)]+['total'])
    for i, name in enumerate(names):
        data = utils.get_final_pockmark_results(name)
        table.iloc[:-1, 2*i] = data['discharge']
        table.iloc[:-1, 2*i+1] = data['concentration']
        table.loc['total'].iloc[2*i] = data['discharge'].sum()
        table.loc['total'].iloc[2*i+1] = np.sum([discharge*conc for discharge, conc in zip(data['discharge'], data['concentration'])])/np.sum(data['discharge'])

    string = table.to_latex(float_format="%.2f")
    if not os.path.exists(f"/home/superuser/objective_3/outputs/comparisons/scenarios_{set}"):
        os.makedirs(f"/home/superuser/objective_3/outputs/comparisons/scenarios_{set}")
    with open(f"/home/superuser/objective_3/outputs/comparisons/scenarios_{set}/pockmark_discharges_by_scenario.txt", "w") as text_file:
        text_file.write(string)

def create_multi_hypothesis_scenario_table(names, type):

    modelled = [utils.get_salinity_and_flux_at_measurement_cells(name) for name in names]
    data = pd.read_csv("/home/superuser/objective_3/data/observations/observations.csv", index_col=0,
                       dtype={'point': str, 'concentration_average': float, 'discharge_min': float,
                              'discharge_max': float})

    if type == "conc":
        table = pd.DataFrame(index=['mc2-1', 'mc2-2', 'mc2-3'], columns=[f"C [PSU] {name}" for name in names]+["C [PSU] (Hoffman et al. 2023)"])

        for point in table.index:
            for i, name in enumerate(names):
                point_temp = point.replace("-", "")
                table.loc[point, f"C [PSU] {name}"] = f"{modelled[i].loc[point_temp, 'conc_min']:.1f} -- {modelled[i].loc[point_temp, 'conc_max']:.1f}"
            table.loc[point, "C [PSU] (Hoffman et al. 2023)"] = f"{data.loc[point, 'concentration_average']:.1f}"

    string = table.to_latex(float_format="%.2f")
    if not os.path.exists(f"/home/superuser/objective_3/outputs/comparisons"):
        os.makedirs(f"/home/superuser/objective_3/outputs/comparisons")
    with open(f"/home/superuser/objective_3/outputs/comparisons/{type}_{"_".join(names)}.txt", "w") as text_file:
        text_file.write(string)


if __name__=="__main__":
    # create_observation_comparison_table("1average")
    pockmark_discharges_by_scenario('base', ['1average', '3slr2', 'hk1ave_new', "hkslr2_new", 'c1ave_new', "cslr2_new"])
    # create_observation_comparison_table("hk1summer")
    # create_multi_hypothesis_scenario_table(["1summer", "hk1summer"], "conc")
