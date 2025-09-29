import pandas as pd
from project_base import project_dir

def find_medium_summer_values(df):
    df['date'] = pd.to_datetime(df['date'])
    recent = df[df['date'] >= '2016-01-01']
    recent['date'] = recent['date'].dt.month
    summer = recent[(recent['date'] >= 1) | (recent['date'] <= 3)]
    mcew_median = summer['McEw_Sh'].median()
    somes_median = summer['Somes'].median()
    return mcew_median, somes_median

def find_summer_values(print=False):
    current = pd.read_csv(project_dir.joinpath('data', 'regional_modelling', 'current_outputs.csv'))
    SLR1 = pd.read_csv(project_dir.joinpath('data', 'regional_modelling', 'SLR1_outputs.csv'))
    SLR2 = pd.read_csv(project_dir.joinpath('data', 'regional_modelling', 'SLR2_outputs.csv'))

    names = ['current', 'SLR1', 'SLR2']
    dfs = [current, SLR1, SLR2]

    values = {}
    for name, df in zip(names, dfs):
        values[name] = {}
        mcew_median, somes_median = find_medium_summer_values(df)
        if print:
            print(f"{name} McEwan minimum = {df['McEw_Sh'].min()}")
            print(f"{name} Somes minimum = {df['Somes'].min()}")
            print(f"{name} McEwan median = {mcew_median}")
            print(f"{name} Somes median = {somes_median}")
        values[name]['mcew_median'] = mcew_median
        values[name]['somes_median'] = somes_median
        values[name]['mcew_min'] = df['McEw_Sh'].min()
        values[name]['somes_min'] = df['Somes'].min()

    if not print:
        return values

if __name__=="__main__":
    find_summer_values()

