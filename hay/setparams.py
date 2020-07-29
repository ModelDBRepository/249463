from neuron import h
import mutation_stuff
import protocol

defValsMut = mutation_stuff.getdefvals()
defVals = protocol.get_defvals()

def setparams(params):
  params_copy = params.copy()
  keys = params.keys()
  for ikey in range(0,len(keys)):
    params_copy[keys[ikey]] = params[keys[ikey]]*defVals[keys[ikey]]

  #Apply the new parameter values:                                                                                                                                                                        
  IhChanged = False
  CaHVAChanged = False
  CaLVAChanged = False
  aIh = defVals['D_aIh']
  bIh = defVals['D_bIh']
  print "setparams: params_copy = "+str(params_copy)
  for ikey in range(0,len(keys)):
    key = keys[ikey]
    section = []
    if key[0:2] == "A_" or key[0:2] == "B_" or key[0:2] == "S_" or key[0:2] == "*_":
      if key[0:2] == "A_":
        section = "apical"
      if key[0:2] == "B_":
        section = "basal"
      if key[0:2] == "S_":
        section = "somatic"
      if key[0:2] == "*_":
        section = "all"
      if len(section) > 0:
        print("forsec L5PC."+section+" "+key[2:]+" = "+str(params_copy[key]))
        h("forsec L5PC."+section+" "+key[2:]+" = "+str(params_copy[key]))
    elif key[0:2] == "D_":
      if key[2:] == "aIh":
        IhChanged = True
        aIh = params_copy[key]
      if key[2:] == "bIh":
        IhChanged = True
        bIh = params_copy[key]
      if key == "aCa_HVA":
        CaHVAChanged = True
      if key[2:] == "aCa_LVAst":
        CaLVAChanged = True
  if IhChanged:
    print("""L5PC.distribute_channels(\"apic\",\"gIhbar_Ih\",2,"""+str(aIh)+""",3.6161,0.0,"""+str(bIh)+""",0.00020000000)""")
    h("""L5PC.distribute_channels(\"apic\",\"gIhbar_Ih\",2,"""+str(aIh)+""",3.6161,0.0,"""+str(bIh)+""",0.00020000000)""")
  if CaLVAChanged:
    print("""L5PC.distribute_channels(\"apic\",\"gCa_LVAstbar_Ca_LVAst\",3,1.000000,0.01,685.000000,885.000000,"""+str(params_copy['D_aCa_LVAst'])+""")""")
    h("""L5PC.distribute_channels(\"apic\",\"gCa_LVAstbar_Ca_LVAst\",3,1.000000,0.01,685.000000,885.000000,"""+str(params_copy['D_aCa_LVAst'])+""")""")
  if CaHVAChanged:
    print("""L5PC.distribute_channels(\"apic\",\"gCa_HVAbar_Ca_HVA\",3,1.000000,0.100000,685.000000,885.000000,"""+str(params_copy['D_aCa_HVA'])+""")""")
    h("""L5PC.distribute_channels(\"apic\",\"gCa_HVAbar_Ca_HVA\",3,1.000000,0.100000,685.000000,885.000000,"""+str(params_copy['D_aCa_HVA'])+""")""")



