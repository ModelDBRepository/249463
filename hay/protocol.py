from neuron import h


condSuffixes = ['Ca_HVA','Ca_LVAst','Ih','Im','K_Pst','K_Tst','NaTa_t','Nap_Et2','SK_E2','SKv3_1']
condNames = ['gCa_HVAbar_Ca_HVA', 'gCa_LVAstbar_Ca_LVAst', 'gIhbar_Ih', 'gImbar_Im', 'gK_Pstbar_K_Pst', 'gK_Tstbar_K_Tst', 'gNaTa_tbar_NaTa_t', 'gNap_Et2bar_Nap_Et2', 'gSK_E2bar_SK_E2', 'gSKv3_1bar_SKv3_1']

defValsSomatic = {
  'S_gIhbar_Ih': 0.0002,
  'S_g_pas': 0.0000338,
  'S_decay_CaDynamics_E2': 460.0,
  'S_gamma_CaDynamics_E2': 0.000501,
  'S_gCa_LVAstbar_Ca_LVAst': 0.00343,
  'S_gCa_HVAbar_Ca_HVA': 0.000992,
  'S_gSKv3_1bar_SKv3_1': 0.693,
  'S_gSK_E2bar_SK_E2': 0.0441,
  'S_gK_Tstbar_K_Tst': 0.0812,
  'S_gK_Pstbar_K_Pst': 0.00223,
  'S_gNap_Et2bar_Nap_Et2': 0.00172,
  'S_gNaTa_tbar_NaTa_t': 2.04
}
defValsApical = {
  'A_decay_CaDynamics_E2': 122,
  'A_gamma_CaDynamics_E2': 0.000509,
  'A_gSK_E2bar_SK_E2': 0.0012,
  'A_gSKv3_1bar_SKv3_1': 0.000261,
  'A_gNaTa_tbar_NaTa_t': 0.0213,
  'A_gImbar_Im': 0.0000675,
  'A_g_pas': 0.0000589
}
defValsBasal = {
  'B_gIhbar_Ih': 0.0002,
  'B_g_pas': 0.0000467
}
defValsApical = {
  'A_decay_CaDynamics_E2': 122,
  'A_gamma_CaDynamics_E2': 0.000509,
  'A_gSK_E2bar_SK_E2': 0.0012,
  'A_gSKv3_1bar_SKv3_1': 0.000261,
  'A_gNaTa_tbar_NaTa_t': 0.0213,
  'A_gImbar_Im': 0.0000675,
  'A_g_pas': 0.0000589
}
defValsBasal = {
  'B_gIhbar_Ih': 0.0002,
  'B_g_pas': 0.0000467
}
defValsDists = {'D_aIh': -0.8696,
                'D_bIh': 2.0870,
                'D_aCa_LVAst': 0.0187,
                'D_aCa_HVA': 0.000555
}

variables = [['S_gIhbar_Ih', 0.5, 2.0],
             ['S_g_pas', 0.5, 2.0],
             ['S_decay_CaDynamics_E2', 0.8, 1.25],
             ['S_gamma_CaDynamics_E2', 0.8, 1.25],
             ['S_gCa_LVAstbar_Ca_LVAst', 0.5, 2.0],
             ['S_gCa_HVAbar_Ca_HVA', 0.2, 1.25], #observed to induce bursting if (radically) decreased
             ['S_gSKv3_1bar_SKv3_1', 0.5, 2.0],
             ['S_gSK_E2bar_SK_E2', 0.5, 2.0],
             ['S_gK_Tstbar_K_Tst', 0.5, 2.0],
             ['S_gK_Pstbar_K_Pst', 0.5, 2.0],
             ['S_gNap_Et2bar_Nap_Et2', 0.5, 2.0],
             ['S_gNaTa_tbar_NaTa_t', 0.5, 2.0],
             ['A_decay_CaDynamics_E2', 0.8, 1.25],
             ['A_gamma_CaDynamics_E2', 0.8, 1.25],
             ['A_gSK_E2bar_SK_E2', 0.5, 2.0],
             ['A_gSKv3_1bar_SKv3_1', 0.5, 2.0],
             ['A_gNaTa_tbar_NaTa_t', 0.8, 4.0], #observed to induce bursting if (radically) increased
             ['A_gImbar_Im', 0.5, 2.0],
             ['A_g_pas', 0.5, 2.0],
             ['B_gIhbar_Ih', 0.5, 2.0],
             ['B_g_pas', 0.5, 2.0],
             ['D_aIh',  0.5, 2.0],
             ['D_bIh', 0.5, 2.0],
             ['D_aCa_LVAst', 0.5, 2.0],
             ['D_aCa_HVA', 0.5, 2.0]]

defVals = defValsSomatic.copy()
defVals.update(defValsApical)
defVals.update(defValsBasal)
defVals.update(defValsDists)


#                   NAME    TYPE      WHERE    AT A FIXED POS?   AMPLITUDE NAME
#                                              OR AT A DISTANCE                
#                                              FROM SOMA?                      
stimulus_types = [ ["st1",  "IClamp", "soma", ["fixed", 0.5],    "amp" ],
                   ["syn1", "epsp",   "apic", ["distance", 620], "imax" ] ]

#                   WHAT TYPE OF OUTPUT? VOLTAGE/[CA] AT A
#                   FIXED TIME (LAST TVEC BEFORE GIVEN TIME)
#                   OR MAXIMUM/NSPIKES DURING A GIVEN INTERVAL?
data_storage_types = [ ["fixed", 13000],
                       ["max", [10000,10200] ],
                       ["trace", [9950+5*x for x in range(0,51)] ],
                       ["highrestrace", [9950+x for x in range(0,251)] ],
                       ["highrestraceandspikes", [9950+x for x in range(0,251)] ],
                       ["nspikes", [12000, 15000] ],
                       ["nspikesandothers", [12000, 15000] ] ]


#            STIMULUS_TYPE   AMPLITUDE NAME  PARAMETERS        
#                                                              
#                                                              
stimuli = [ [ 0, [ ["del", 10000], ["dur", 3000] ] ],
            [ 0, [ ["del", 10000], ["dur", 100] ] ],
            [ 1, [ ["onset", 10000], ["tau0", 0.5], ["tau1", 5.0] ] ],
            [ 0, [ ["del", 10000], ["dur", 5000] ] ],
            [ 0, [ ["del", 10000], ["dur", 5] ] ] ]

#                NAME     WHAT    WHERE
#                                   WHICH    AT WHICH
#                                   BRANCH   LOCATIONS
recordings = [ [ "vsoma",  "v",   [ ["soma", [0.5] ] ] ],
               [ "vdend",  "v",   [ ["apic", [0.05*x for x in range(0,21)] ], 
                                    ["dend", [0.05*x for x in range(0,21)] ] ] ],
               [ "cadend", "cai", [ ["apic", [0.05*x for x in range(0,21)] ] ] ] ]

#                STIMULUS INDEX
#                        AMPLITUDES
#                                                 RECORDING INDEX
#                                                         DATA STORAGE INDEX
setup =      [ [ [1],   [0.25, 0.5],              0,      4 ],            # MEMBRANE POTENTIAL TRACE RESPONSE TO A 100 ms DC,
               [ [4],   [1.9],                    0,      4 ],            # MEMBRANE POTENTIAL TRACE RESPONSE TO A SHORT SOMATIC STIMULUS
               [ [2],   [1.5],                    0,      4 ],            # MEMBRANE POTENTIAL TRACE RESPONSE TO A SHORT SOMATIC STIMULUS
               [ [4,2], [[1.9, 0.5]],             0,      4 ],            # MEMBRANE POTENTIAL TRACE RESPONSE TO COMBINATION OF SHORT SOMATIC AND APICAL STIMULI
               [ [3],   [0.78, 1.0, 1.9],         0,      6 ] ]           # NUMBERS OF SPIKES, CVISI AND NSPIKESPERBURST AS A RESPONSE TO A LONG DC

stop_if_no_spikes = [ [0, 1], [0], [0], [1], [1, 1, 1] ]
stop_if_more_spikes_or_as_many_as = [ [2, 5], [5], [5], [5], [100, 200, 300] ]  # when nSpikes >= x)

# Return the full matrix of the step-wise fitting variables
def get_variables():
  global variables
  return variables
def get_defvals():
  global defVals
  return defVals

# Return the stimulus array
def get_stimuli():
  global stimuli
  return stimuli

# Return the stimulus types array
def get_stimulus_types():
  global stimulus_types
  return stimulus_types

# Return the data storage array
def get_data_storage_types():
  global data_storage_types
  return data_storage_types

# Return the recordings array
def get_recordings():
  global recordings
  return recordings

# Return the stimulus-recording setup array
def get_setup():
  global setup
  return setup

def get_nspike_restrictions():
  global stop_if_no_spikes, stop_if_more_spikes_or_as_many_as
  return [stop_if_no_spikes, stop_if_more_spikes_or_as_many_as]




