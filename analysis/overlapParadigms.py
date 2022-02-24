##############################################################
### Definitions of assembly neurons in different paradigms ###
###           of overlap between cell assemblies           ###
##############################################################

### Copyright 2020-2022 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: jannik.lubo[at]gmx.de

import numpy as np

# coreDefinitions
# Returns the neuron numbers belonging to each of three assemblies in a given paradigm
# paradigm: name of the paradigm to consider
# core_size: the size of one assembly
# return: the three assemblies
def coreDefinitions(paradigm, core_size = 600):

	# full_overlap
	# Returns the neuron number belonging to each of three assemblies with equal overlaps
	# overlap: overlap between each two cell assemblies (0.1 equals "OVERLAP10" and so on)
	# return: the three assemblies
	def full_overlap(overlap):
		tot_wo_overlap = 1-overlap
		half_overlap = overlap / 2
		tot2_wo_overlap_oh = 2-overlap-half_overlap
		core1 = np.arange(core_size)
		core2 = np.arange(int(np.round(tot_wo_overlap*core_size)), int(np.round(tot_wo_overlap*core_size))+core_size)
		core3 = np.concatenate(( np.arange(int(np.round(tot2_wo_overlap_oh*core_size)), int(np.round(tot2_wo_overlap_oh*core_size))+int(np.round(tot_wo_overlap*core_size))), \
			                 np.arange(int(np.round(half_overlap*core_size))), \
			                 np.arange(int(np.round(tot_wo_overlap*core_size)), int(np.round(tot_wo_overlap*core_size))+int(np.round(half_overlap*core_size))) ))
		
		return (core1, core2, core3)

	# no_ABC_overlap
	# Returns the neuron number belonging to each of three assemblies with equal overlaps but no common overlap
	# overlap: overlap between each two cell assemblies (0.1 equals "OVERLAP10 no ABC" and so on)
	# return: the three assemblies
	def no_ABC_overlap(overlap):
		tot_wo_overlap = 1-overlap
		tot2_wo_overlap_oh = 2 - 2*overlap
		core1 = np.arange(core_size)
		core2 = np.arange(int(np.round(tot_wo_overlap*core_size)), int(np.round(tot_wo_overlap*core_size))+core_size)
		core3 = np.concatenate(( np.arange(int(np.round(tot2_wo_overlap_oh*core_size)), int(np.round(tot2_wo_overlap_oh*core_size))+int(np.round(tot_wo_overlap*core_size))), \
			                 np.arange(int(np.round(overlap*core_size))) ))
		return (core1, core2, core3)

	# no_AC_no_ABC_overlap
	# Returns the neuron number belonging to each of three assemblies with "no AC, no ABC" overlap
	# overlap: overlap between each two cell assemblies (0.1 equals "OVERLAP10 no AC, no ABC" and so on)
	# return: the three assemblies
	def no_AC_no_ABC_overlap(overlap):
		tot_wo_overlap = 1-overlap
		tot2_wo_overlap_oh = 2 - 2*overlap
		core1 = np.arange(core_size)
		core2 = np.arange(int(np.round(tot_wo_overlap*core_size)), int(np.round(tot_wo_overlap*core_size))+core_size)
		core3 = np.arange(int(np.round(tot2_wo_overlap_oh*core_size)), int(np.round(tot2_wo_overlap_oh*core_size))+core_size)
		return (core1, core2, core3)

	# no_BC_no_ABC_overlap
	# Returns the neuron number belonging to each of three assemblies with "no BC, no ABC" overlap
	# overlap: overlap between each two cell assemblies (0.1 equals "OVERLAP10 no BC, no ABC" and so on)
	# return: the three assemblies
	def no_BC_no_ABC_overlap(overlap):
		tot_wo_overlap = 1-overlap
		core1 = np.arange(core_size)
		core2 = np.arange(int(np.round(tot_wo_overlap*core_size)), int(np.round(tot_wo_overlap*core_size))+core_size)
		core3 = np.concatenate(( np.arange(int(np.round(tot_wo_overlap*core_size))+core_size, 2*int(np.round(tot_wo_overlap*core_size))+core_size), \
			                 np.arange(int(np.round(overlap*core_size))) ))
		return (core1, core2, core3)

	# handling the different overlap paradigms:
	if paradigm == "NOOVERLAP":
		core1 = np.arange(core_size)
		core2 = np.arange(core_size, 2*core_size)
		core3 = np.arange(2*core_size, 3*core_size)
	elif paradigm == "OVERLAP10":
		core1, core2, core3 = full_overlap(0.1)
	elif paradigm == "OVERLAP10 no ABC":
		core1, core2, core3 = no_ABC_overlap(0.1)
	elif paradigm == "OVERLAP10 no AC, no ABC":
		core1, core2, core3 = no_AC_no_ABC_overlap(0.1)
	elif paradigm == "OVERLAP10 no BC, no ABC":
		core1, core2, core3 = no_BC_no_ABC_overlap(0.1)
	elif paradigm == "OVERLAP15":
		core1, core2, core3 = full_overlap(0.15)
	elif paradigm == "OVERLAP15 no ABC":
		core1, core2, core3 = no_ABC_overlap(0.15)
	elif paradigm == "OVERLAP15 no AC, no ABC":
		core1, core2, core3 = no_AC_no_ABC_overlap(0.15)
	elif paradigm == "OVERLAP15 no BC, no ABC":
		core1, core2, core3 = no_BC_no_ABC_overlap(0.15)
	elif paradigm == "OVERLAP20":
		core1, core2, core3 = full_overlap(0.2)
	elif paradigm == "OVERLAP20 no ABC":
		core1, core2, core3 = no_ABC_overlap(0.2)
	elif paradigm == "OVERLAP20 no AC, no ABC":
		core1, core2, core3 = no_AC_no_ABC_overlap(0.2)
	elif paradigm == "OVERLAP20 no BC, no ABC":
		core1, core2, core3 = no_BC_no_ABC_overlap(0.2)
	else:
		raise ValueError("Unknown paradigm: " + paradigm)
		
	return (core1, core2, core3)
