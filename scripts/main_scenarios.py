import modelling
import plots


def new_grid_scenarios():

    # modelling.create_run_model("1average", case='base', time='current',
    #                            scenario='summer_average', init='fresh', tscenario='steady')
    # plots.plot_salinity_at_pockmarks("1average")
    # modelling.create_run_model("2slr1", case='base', time='SLR1',
    #                            scenario='summer_average', init='1average', tscenario='slr1')
    # modelling.create_run_model("3slr2", case='base', time='SLR2',
    #                            scenario='summer_average', init='1average', tscenario='slr2')
    #
    # modelling.create_run_model("4minimum", case='base', time='current',
    #                            scenario='summer_minimum', init='1average', tscenario='transient')
    # modelling.create_run_model("5slr1_min", case='base', time='SLR1',
    #                            scenario='summer_minimum', init='2slr1', tscenario='transient')
    # modelling.create_run_model("6slr2_min", case='base', time='SLR2',
    #                            scenario='summer_minimum', init='3slr2', tscenario='transient')
    #
    # modelling.create_run_model("c7swi_b", case='base', time='current',
    #                            scenario='swi_b', init='1average', tscenario='transient', hypothesis="conduits")
    # modelling.create_run_model("8slr1_swi", case='base', time='SLR1',
    #                            scenario='swi', init='2slr1', tscenario='transient')
    # modelling.create_run_model("9slr2_swi", case='base', time='SLR2',
    #                            scenario='swi', init='3slr2', tscenario='transient')

    names = ['1average'] # '3slr2', '4minimum', '5slr1_min', '6slr2_min', '7swi', '8slr1_swi', '9slr2_swi']
    for name in names:
        plots.plot_salinity_at_pockmarks(name)

def new_grid_high_k_scenarios():
    # modelling.create_run_model("hk1average", case='base', time='current',
    #                             scenario='summer_average', init='fresh', tscenario='steady', hypothesis="high_pockmark_k")
    # modelling.create_run_model("hk2slr1", case='base', time='SLR1',
    #                            scenario='summer_average', init='hk1average', tscenario='slr1', hypothesis="high_pockmark_k")
    # modelling.create_run_model("hk3slr2", case='base', time='SLR2',
    #                            scenario='summer_average', init='hk1average', tscenario='slr2', hypothesis="high_pockmark_k")
    #
    # modelling.create_run_model("hk4minimum", case='base', time='current',
    #                            scenario='summer_minimum', init='hk1average', tscenario='transient', hypothesis="high_pockmark_k")
    # modelling.create_run_model("hk5slr1_min", case='base', time='SLR1',
    #                            scenario='summer_minimum', init='hk2slr1', tscenario='transient', hypothesis="high_pockmark_k")
    # modelling.create_run_model("hk6slr2_min", case='base', time='SLR2',
    #                            scenario='summer_minimum', init='hk3slr2', tscenario='transient', hypothesis="high_pockmark_k")
    #
    # modelling.create_run_model("hk7swi", case='base', time='current',
    #                            scenario='swi', init='hk1average', tscenario='transient', hypothesis="high_pockmark_k")
    # modelling.create_run_model("hk8slr1_swi", case='base', time='SLR1',
    #                            scenario='swi', init='hk2slr1', tscenario='transient', hypothesis="high_pockmark_k")
    # modelling.create_run_model("hk9slr2_swi", case='base', time='SLR2',
    #                            scenario='swi', init='hk3slr2', tscenario='transient', hypothesis="high_pockmark_k")

    names = ['hk4minimum', 'hk5slr1_min', 'hk6slr2_min', 'hk7swi', 'hk8slr1_swi', 'hk9slr2_swi']
    for name in names:
        plots.plot_salinity_at_pockmarks(name)

def new_grid_conduit_scenarios():
    # modelling.create_run_model("c1average", case='base', time='current',
    #                             scenario='summer_average', init='fresh', tscenario='steady', hypothesis="conduits")
    # modelling.create_run_model("c2slr1", case='base', time='SLR1',
    #                            scenario='summer_average', init='c1average', tscenario='slr1', hypothesis="conduits")
    # modelling.create_run_model("c3slr2", case='base', time='SLR2',
    #                            scenario='summer_average', init='c1average', tscenario='slr2', hypothesis="conduits")
    #
    # modelling.create_run_model("c4minimum", case='base', time='current',
    #                            scenario='summer_minimum', init='c1average', tscenario='transient', hypothesis="conduits")
    # modelling.create_run_model("c5slr1_min", case='base', time='SLR1',
    #                            scenario='summer_minimum', init='c2slr1', tscenario='transient', hypothesis="conduits")
    # modelling.create_run_model("c6slr2_min", case='base', time='SLR2',
    #                            scenario='summer_minimum', init='c3slr2', tscenario='transient', hypothesis="conduits")
    #
    # modelling.create_run_model("c7swi", case='base', time='current',
    #                            scenario='swi', init='c1average', tscenario='transient', hypothesis="conduits")
    # modelling.create_run_model("c8slr1_swi", case='base', time='SLR1',
    #                            scenario='swi', init='c2slr1', tscenario='transient', hypothesis="conduits")
    # modelling.create_run_model("c9slr2_swi", case='base', time='SLR2',
    #                            scenario='swi', init='c3slr2', tscenario='transient', hypothesis="conduits")

    names = ['c1average', 'c2slr1', 'c3slr2', 'c4minimum', 'c5slr1_min', 'c6slr2_min', 'c7swi', 'c8slr1_swi', 'c9slr2_swi']
    for name in names:
        plots.plot_salinity_at_pockmarks(name)

def low_swi_scenarios():
    # modelling.create_run_model("big1_swi", case='base', time='current',
    #                            scenario='swi_b', init='1average', tscenario='transient')
    modelling.create_run_model("hkbig1_swi_n", case='base', time='current',
                               scenario='swi_b', init='hk1ave_new', tscenario='transient', hypothesis="high_pockmark_k")
    modelling.create_run_model("cbig1_swi_n", case='base', time='current',
                               scenario='swi_b', init='c1ave_new', tscenario='transient', hypothesis="conduits")
    names = ['hkbig1_swi_n', 'cbig1_swi_n']
    for name in names:
        plots.plot_salinity_at_pockmarks(name)

def low_long_swi_scenarios():
    # modelling.create_run_model("bl1_swi", case='base', time='current',
    #                            scenario='swi_b', init='1average', tscenario='transient_l')
    # modelling.create_run_model("hkbl1_swi_n", case='base', time='current',
    #                            scenario='swi_b', init='hk1ave_new', tscenario='transient_l', hypothesis="high_pockmark_k")
    # modelling.create_run_model("cbl1_swi_n", case='base', time='current',
    #                            scenario='swi_b', init='c1ave_new', tscenario='transient_l', hypothesis="conduits")
    modelling.create_run_model("hhhkbl1_swi_n", case='base', time='current',
                               scenario='swi_b', init="hhhk_average", tscenario='transient_l', hypothesis="high_pockmark_k", patch_multiplier=100)
    # names = ["bl1_swi", 'hkbl1_swi_n', 'cbl1_swi_n']
    # for name in names:
    #     plots.plot_salinity_at_pockmarks(name)

def new_conduits_and_patches_steady():
    # modelling.create_run_model("hk1ave_new", case='base', time='current',
    #                              scenario='summer_average', init='fresh', tscenario='steady', hypothesis="high_pockmark_k")
    # modelling.create_run_model("c1ave_new", case='base', time='current',
    #                              scenario='summer_average', init='fresh', tscenario='steady', hypothesis="conduits")
    # plots.plot_alternate_hypotheses_pockmarks_and_model_vs_observed(['1average', 'hk1ave_new', 'c1ave_new'], 1)
    # modelling.create_run_model("hkslr2_new", case='base', time='SLR2',
    #                            scenario='summer_average', init="hk1ave_new", tscenario='slr2', hypothesis="high_pockmark_k"
    #                            )
    # modelling.create_run_model("cslr2_new", case='base', time='SLR2',
    #                            scenario='summer_average', init="c1ave_new", tscenario='slr2', hypothesis="conduits")
    names = ["hk1ave_new", "c1ave_new"]
    for name in names:
       plots.plot_salinity_at_pockmarks(name)


def dropping_head_runs():
    # modelling.create_run_model("dropping", case='base', time='current', scenario='dropping_head',
    #                         init='', tscenario='dropping_head')
    # modelling.create_run_model("hkdropping", case='base', time='current', scenario='dropping_head',
    #                            init='fresh', tscenario='dropping_head', hypothesis="high_pockmark_k")
    modelling.create_run_model("cdropping", case='base', time='current', scenario='dropping_head',
                               init='c1ave_new', tscenario='dropping_head', hypothesis="conduits")
    # plots.plot_salinity_at_pockmarks("dropping")

def patches_higher_k():
    modelling.create_run_model("hhk_average", case='base', time='current',
                               scenario='summer_average', init='fresh', tscenario='steady', hypothesis="high_pockmark_k", patch_multiplier=20)
    modelling.create_run_model("hhk_slr2", case='base', time='SLR2',
                               scenario='summer_average', init="hhk_average", tscenario='slr2', hypothesis="high_pockmark_k", patch_multiplier=20)
    modelling.create_run_model('hhk_swib', case='base', time='current',
                               scenario='swi_b', init='hhk_average', tscenario='transient', hypothesis="high_pockmark_k", patch_multiplier=20)
    names = ['hhk_average', 'hhk_slr2', 'hhk_swib']
    for name in names:
        plots.plot_salinity_at_pockmarks(name)

def patches_even_higher_k():
    modelling.create_run_model("hhhk_average", case='base', time='current',
                               scenario='summer_average', init='fresh', tscenario='steady', hypothesis="high_pockmark_k", patch_multiplier=100)
    modelling.create_run_model("hhhk_slr2", case='base', time='SLR2',
                               scenario='summer_average', init="hhk_average", tscenario='slr2', hypothesis="high_pockmark_k", patch_multiplier=100)
    modelling.create_run_model('hhhk_swib', case='base', time='current',
                               scenario='swi_b', init='hhk_average', tscenario='transient', hypothesis="high_pockmark_k", patch_multiplier=100)
    names = ['hhhk_average', 'hhhk_slr2', 'hhhk_swib']
    for name in names:
        plots.plot_salinity_at_pockmarks(name)

def conduits_higher_k():
    modelling.create_run_model("ch_average", case='base', time='current',
                               scenario='summer_average', init='fresh', tscenario='steady', hypothesis="conduits", conduit_k=100)
    modelling.create_run_model("ch_slr2", case='base', time='SLR2',
                               scenario='summer_average', init="ch_average", tscenario='slr2', hypothesis="conduits", conduit_k=100)
    modelling.create_run_model('ch_swib', case='base', time='current',
                               scenario='swi_b', init='ch_average', tscenario='transient', hypothesis="conduits", conduit_k=100)
    names = ['ch_average', 'ch_slr2', 'ch_swib']
    for name in names:
        plots.plot_salinity_at_pockmarks(name)

def results():
    # names = ['t1summer', 't2slr', 't3low', 't4low_slr', 't5swi', 't6swi_slr']
    # names = ['c1summer', 'c2slr', 'c3low', 'c4low_slr', 'c5swi', 'c6swi_slr']
    # latex_tables.pockmark_discharges_by_scenario("first_scenarios", names)
    # plots.plot_changes_and_examples('conduits', names)
    pass


if __name__ == "__main__":
    # new_grid_high_k_scenarios()
    # new_grid_conduit_scenarios()
    # results()
    new_conduits_and_patches_steady()
    # low_swi_scenarios()
    # dropping_head_runs()
    # patches_higher_k()
    # conduits_higher_k()
    # low_long_swi_scenarios()
    # new_grid_scenarios()
    # new_grid_scenarios()