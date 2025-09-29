import plots
import utils
# modelling.create_run_steady_state('test_short', tscenario='short')
# plots.plot_boundary_conditions(False)
# plots.plot_salinity_at_pockmarks('disu_test')
#modelling.create_run_steady_state('medium', tscenario='short')
# utils.get_final_results("disu_test")
# names = ['t1summer', 't2slr', 't3low', 't4low_slr', 't5swi', 't6swi_slr']
# names = ['1summer', '2slr', '3low', '4low_slr', '5swi', '6swi_slr']
# plots.plot_changes_and_examples('first', names)

#plots.plot_salinity_at_pockmarks('1average')
utils.calc_pockmark_and_model_areas(buffer=50)
# main_scenarios.new_conduits_and_patches_steady()