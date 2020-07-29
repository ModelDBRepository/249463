def getMT():
 #       List of genes
 #         List of mutations
 #           List of groups of variables
 #             Pair (variable list + range)
 #               List of variables
 MT = [] 
 #CACNA1C:
 MT.append([ [ [ ['offm_cah', 'offmt_cah'], -25.9 ],                       #http://www.ncbi.nlm.nih.gov/pubmed/19265197
               [ ['offh_cah', 'offht_cah'], -27.0 ] ],
             [ [ ['offm_cah', 'offht_cah'], -37.3 ],                       #http://www.ncbi.nlm.nih.gov/pubmed/19265197
               [ ['offh_cah', 'offht_cah'], -30.0 ] ],
             [ [ ['offm_cah', 'offmt_cah'], [-31.4, +7.0] ],               #http://www.ncbi.nlm.nih.gov/pubmed/21685391                                                        
               [ ['slom_cah', 'slomt_cah'], [0.85, 1.45] ],
               [ ['offh_cah', 'offht_cah'], [-28.5, +16.3] ],
               [ ['sloh_cah', 'sloht_cah'], [0.72, 1.38] ] ],
             [ [ ['offm_cah', 'offmt_cah'], [-38.5, +12.9] ],              #http://www.ncbi.nlm.nih.gov/pubmed/16157588                                                        
               [ ['slom_cah', 'slomt_cah'], [0.46, 1.56] ] ],
             [ [ ['offm_cah', 'offmt_cah'], [-27.8, +8.7] ],               #http://www.ncbi.nlm.nih.gov/pubmed/18836301                                                        
               [ ['slom_cah', 'slomt_cah'], [0.89, 1.14] ],
               [ ['offh_cah', 'offht_cah'], [-19.1, +4.7] ] ],
             [ [ ['offm_cah', 'offmt_cah'], [-11.2, +1.0] ],               #http://www.ncbi.nlm.nih.gov/pubmed/15299022                                                        
               [ ['offh_cah', 'offht_cah'], -3.1 ],
               [ ['sloh_cah', 'sloht_cah'], 1.24 ] ] ])

 #CACNB2:
 MT.append([ [ [ ['offh_cah', 'offht_cah'], -5.2 ],                        #http://www.ncbi.nlm.nih.gov/pubmed/19358333
               [ ['sloh_cah', 'sloht_cah'], 0.69 ] ], 
             [ [ 'tauhmax_cah', 1.7 ] ],                                   #http://www.ncbi.nlm.nih.gov/pubmed/7723731
             [ [ ['offm_cah', 'offmt_cah'], [-4.9, 4.9] ],                 #http://www.ncbi.nlm.nih.gov/pubmed/19723630
               [ ['offh_cah', 'offht_cah'], [-5.1, 5.1] ],
               [ 'taummax_cah', [0.6, 1.68] ],
               [ 'tauhmax_cah', [0.6, 1.66] ] ],
             [ [ 'taummax_cah', 1.26] ] ])                                 #http://www.ncbi.nlm.nih.gov/pubmed/20025708
 #CACNA1D:
 MT.append([ [ [ ['offm_cah', 'offmt_cah'], -10.9 ],                       #http://www.ncbi.nlm.nih.gov/pubmed/21998309 and
               [ ['slom_cah', 'slomt_cah'], 0.73 ],                        #http://www.ncbi.nlm.nih.gov/pubmed/21998310
               [ ['offh_cah', 'offht_cah'], [-3.0, 3.5] ],                 #(42A)
               [ ['sloh_cah', 'sloht_cah'], 0.81 ],
               [ 'tauhmax_cah', 1.25 ] ],
             [ [ ['offm_cah', 'offmt_cah'], [-10.6, 3.4] ],                #http://www.ncbi.nlm.nih.gov/pubmed/21998309 and
               [ ['slom_cah', 'slomt_cah'], [0.8, 1.12] ],                 #http://www.ncbi.nlm.nih.gov/pubmed/21998310
               [ ['offh_cah', 'offht_cah'], [-5.3, 1.2] ],                 #(43S)
               [ ['sloh_cah', 'sloht_cah'], 0.66 ],
               [ 'tauhmax_cah', 0.72 ] ],
             [ [ ['offm_cah', 'offmt_cah'], 6.6 ],                         #http://www.ncbi.nlm.nih.gov/pubmed/20951705 and
               [ ['slom_cah', 'slomt_cah'], [0.75, 1.19] ],                #http://www.ncbi.nlm.nih.gov/pubmed/21054386
               [ 'tauhmax_cah', [0.5, 1.12] ] ],                           #(CaV1.3 KO)
             [ [ ['offm_cah', 'offmt_cah'], -9.8 ],                        #http://www.ncbi.nlm.nih.gov/pubmed/25620733                                                        
               [ ['slom_cah', 'slomt_cah'], 0.8 ],
               [ ['offh_cah', 'offht_cah'], -15.4 ],
               [ ['sloh_cah', 'sloht_cah'], 1.05 ] ],
             [ [ ['offm_cah', 'offmt_cah'], [-24.2, +6.1] ],               #http://www.ncbi.nlm.nih.gov/pubmed/23913004                                                        
               [ ['slom_cah', 'slomt_cah'], [0.7, 1.24] ],
               [ ['offh_cah', 'offht_cah'], -14.5 ],
               [ ['sloh_cah', 'sloht_cah'], [0.72, 1.28] ],
               [ 'taummax_cah', 3.52 ] ],
             [ [ ['offm_cah', 'offmt_cah'], -17.8 ],                       #http://www.ncbi.nlm.nih.gov/pubmed/22760075                                                        
               [ ['slom_cah', 'slomt_cah'], 0.81 ],
               [ 'taummax_cah', [0.77, 1.31] ] ] ])

 #CACNA1I:
 MT.append([ [ [ ['offm_car'], 1.3 ],                                      #http://www.ncbi.nlm.nih.gov/pubmed/15254077
               [ ['offh_car'], 1.6 ],
               [ 'taummax_car', [0.87, 1.45] ],
               [ 'tauhmax_car', 0.8 ] ],
             [ [ 'offm_car', -4.3 ],                                       #http://www.ncbi.nlm.nih.gov/pubmed/12080115                                                        
               [ 'slom_car', 1.14 ],
               [ 'offh_car', -4.4 ],
               [ 'sloh_car', [0.89, 1.04] ],
               [ 'taummax_car', 0.53 ],
               [ 'tauhmax_car', 0.46 ] ] ])
 #CACNA1S:
 MT.append([ [ [ 'taummax_cah', 0.67 ] ],                                  #http://www.ncbi.nlm.nih.gov/pubmed/20861472
             [ [ ['offm_cah', 'offmt_cah'], -30.02 ],                      #http://www.ncbi.nlm.nih.gov/pubmed/19134469
               [ ['slom_cah', 'slomt_cah'], 0.62 ],
               [ 'taummax_cah', 0.49] ],
             [ [ ['offm_cah', 'offmt_cah'], -4.4 ],                        #http://www.ncbi.nlm.nih.gov/pubmed/24240197                                                        
               [ ['slom_cah', 'slomt_cah'], 0.95 ],
               [ ['offh_cah', 'offht_cah'], 20.6 ],
               [ ['sloh_cah', 'sloht_cah'], 0.97] ] ])

 #ATP2A:
 MT.append([ [ [ 'gamma_cad', 0.6 ] ] ])                                   #http://www.ncbi.nlm.nih.gov/pubmed/10970890
                
 #ATP2B:
 MT.append([ [ [ 'taur_cad', 1.97 ] ],                                     #http://www.ncbi.nlm.nih.gov/pubmed/22789621 
             [ [ 'taur_cad', 1.5 ],                                        #http://www.ncbi.nlm.nih.gov/pubmed/21232211
               [ 'cainf_cad' , 1.4 ] ],
             [ [ 'taur_cad', 4.45 ] ],                                     #http://www.ncbi.nlm.nih.gov/pubmed/17234811
             [ [ 'cainf_cad', 1.1 ] ] ])                                   #http://www.ncbi.nlm.nih.gov/pubmed/22047666
 #NRGN:
 MT.append([ [ [ 'gamma_cad', 0.4 ] ] ])                                   #http://www.ncbi.nlm.nih.gov/pubmed/15564582
 #SCN1A:
 MT.append([ [ [ ['offm_na', 'offmt_na'], -0.3 ],                          #http://www.ncbi.nlm.nih.gov/pubmed/18632931
               [ ['offh_na', 'offht_na'], 5 ],
               [ ['slom_na', 'slomt_na'], 1.15 ],
               [ ['sloh_na', 'sloht_na'], 1.23 ] ],
             [ [ ['offm_na', 'offmt_na'], 2.8 ],                           #http://www.ncbi.nlm.nih.gov/pubmed/18632931
               [ ['offh_na', 'offht_na'], 9.6 ],
               [ ['slom_na', 'slomt_na'], 0.984 ],
               [ ['sloh_na', 'sloht_na'], 1.042 ] ],
             [ [ ['offm_na', 'offmt_na'], -4.0 ],                          #http://www.ncbi.nlm.nih.gov/pubmed/21864321                                                        
               [ ['offh_na', 'offht_na'], -5.8 ],
               [ ['slom_na', 'slomt_na'], 0.92 ],
               [ ['sloh_na', 'sloht_na'], 1.13 ],
               [ ['i1_na', 'i2_na'], 1.47 ] ],
             [ [ ['offm_na', 'offmt_na'], -8.1 ],                          #http://www.ncbi.nlm.nih.gov/pubmed/21864321                                                        
               [ ['offh_na', 'offht_na'], 2.2 ],
               [ ['slom_na', 'slomt_na'], 0.97 ],
               [ ['sloh_na', 'sloht_na'], 0.97 ],
               [ ['i1_na', 'i2_na'], 1.59 ] ],
             [ [ ['offm_na', 'offmt_na'], 6.0 ],                           #http://www.ncbi.nlm.nih.gov/pubmed/23398611                                                        
               [ ['slom_na', 'slomt_na'], 1.16 ],
               [ ['i1_na', 'i2_na'], 1.29 ] ],
             [ [ ['offm_na', 'offmt_na'], 10.0 ],                          #http://www.ncbi.nlm.nih.gov/pubmed/16326807                                                        
               [ ['offh_na', 'offht_na'], -0.6 ],
               [ ['slom_na', 'slomt_na'], 1.15 ],
               [ ['sloh_na', 'sloht_na'], 1.14 ] ] ])

 #SCN9A: (Probably no way to include the effects in this model..!)
 MT.append([ [ [ ['offm_na', 'offmt_na'], 0 ] ],                           #http://www.ncbi.nlm.nih.gov/pubmed/22136189
             [ [ ['offm_na', 'offmt_na'], 0 ] ],                           #http://www.ncbi.nlm.nih.gov/pubmed/18945915
             [ [ ['offm_na', 'offmt_na'], 0 ] ],                           #http://www.ncbi.nlm.nih.gov/pubmed/16392115
             [ [ ['offm_na', 'offmt_na'], 0 ] ] ])                         #http://www.ncbi.nlm.nih.gov/pubmed/15958509
               
 #KCNS3:
 MT.append([ [ [ ['a0_kslow', 'a1_kslow', 'a2_kslow'], 2.0 ],              #http://www.ncbi.nlm.nih.gov/pubmed/10484328
               [ ['b0_kslow', 'b11_kslow', 'b2_kslow', 'bb0_kslow', 'bb1_kslow', 'bb2_kslow'], 2.5 ],
               [ 'sloh_kslow', 0.5 ] ] ])
 #KCNN3:
 MT.append([ [ [ 'offc_sk', 0.86 ],                                        #http://www.ncbi.nlm.nih.gov/pubmed/14978258
               [ 'sloc_sk', 1.24 ] ],
             [ [ 'sloc_sk', 0.84 ] ] ])                                    #http://www.ncbi.nlm.nih.gov/pubmed/17167222

 #HCN1:
 MT.append([ [ [ ['off_iH', 'offt1_iH', 'offt2_iH'], -26.5 ],              #http://www.ncbi.nlm.nih.gov/pubmed/17185333
               [ ['slo_iH', 'slot1_iH', 'slot2_iH'], 0.64 ] ],
             [ [ ['off_iH', 'offt1_iH', 'offt2_iH'], [-25.9, 17.7] ],      #http://www.ncbi.nlm.nih.gov/pubmed/12668666                                                        
               [ ['slo_iH', 'slot1_iH', 'slot2_iH'], 0.6 ] ] ])
 #KCNB1:
 MT.append([ [ [ ['offma_kslow', 'offmb_kslow'], 5 ],                      #http://www.ncbi.nlm.nih.gov/pubmed/21455829 (T203K)
               [ 'offh_kslow', 3 ],
               [ ['sloma_kslow', 'slomb_kslow'], 1.11 ],
               [ 'sloh_kslow', 0.86 ],
               [ ['a0_kslow', 'a1_kslow', 'a2_kslow'], 0.5 ],
               [ ['b0_kslow', 'b11_kslow', 'b2_kslow', 'bb0_kslow', 'bb1_kslow', 'bb2_kslow'], 0.53 ] ],
             [ [ ['offma_kslow', 'offmb_kslow'], 1 ],                      #http://www.ncbi.nlm.nih.gov/pubmed/21455829 (T203D)
               [ 'offh_kslow', -6 ],
               [ ['sloma_kslow', 'slomb_kslow'], 1.22 ],
               [ 'sloh_kslow', 1.0 ],
               [ ['a0_kslow', 'a1_kslow', 'a2_kslow'], 0.89 ],
               [ ['b0_kslow', 'b11_kslow', 'b2_kslow', 'bb0_kslow', 'bb1_kslow', 'bb2_kslow'], 1.13 ] ],
             [ [ ['offma_kslow', 'offmb_kslow'], 6 ],                      #http://www.ncbi.nlm.nih.gov/pubmed/21455829 (S347K)
               [ 'offh_kslow', -8 ],
               [ ['sloma_kslow', 'slomb_kslow'], 1.33 ],
               [ 'sloh_kslow', 1.0 ],
               [ ['a0_kslow', 'a1_kslow', 'a2_kslow'], 0.5 ],
               [ ['b0_kslow', 'b11_kslow', 'b2_kslow', 'bb0_kslow', 'bb1_kslow', 'bb2_kslow'], 0.87 ] ],
             [ [ ['offma_kslow', 'offmb_kslow'], -28 ],                    #http://www.ncbi.nlm.nih.gov/pubmed/21455829 (S347D)
               [ 'offh_kslow', -27 ],
               [ ['sloma_kslow', 'slomb_kslow'], 1.11 ],
               [ 'sloh_kslow', 0.71 ],
               [ ['a0_kslow', 'a1_kslow', 'a2_kslow'], 1.13 ],
               [ ['b0_kslow', 'b11_kslow', 'b2_kslow', 'bb0_kslow', 'bb1_kslow', 'bb2_kslow'], 2.27 ] ],
             [ [ ['offma_kslow', 'offmb_kslow'], 14 ],                     #http://www.ncbi.nlm.nih.gov/pubmed/21455829 (T203W)
               [ 'offh_kslow', -21 ],
               [ ['sloma_kslow', 'slomb_kslow'], 2.0 ],
               [ 'sloh_kslow', 1.0 ],
               [ ['a0_kslow', 'a1_kslow', 'a2_kslow'], 0.39 ],
               [ ['b0_kslow', 'b11_kslow', 'b2_kslow', 'bb0_kslow', 'bb1_kslow', 'bb2_kslow'], 1.2 ] ],
             [ [ ['offma_kslow', 'offmb_kslow'], -13 ],                    #http://www.ncbi.nlm.nih.gov/pubmed/21455829 (S347W)
               [ 'offh_kslow', -13 ],
               [ ['sloma_kslow', 'slomb_kslow'], 1.33 ],
               [ 'sloh_kslow', 0.71 ],
               [ ['a0_kslow', 'a1_kslow', 'a2_kslow'], 0.95 ],
               [ ['b0_kslow', 'b11_kslow', 'b2_kslow', 'bb0_kslow', 'bb1_kslow', 'bb2_kslow'], 5.13 ] ] ])
 #CACNB2 reprise:
 MT.append([ [ [ ['offm_cah', 'offmt_cah'], 3.00 ],                        #http://www.ncbi.nlm.nih.gov/pubmed/19723630 (N1 vs N4)
               [ ['offh_cah', 'offht_cah'], 3.48 ],
               [ 'taummax_cah', 1.01 ],                                    #(tau_h_slow considered)                               
               [ 'tauhmax_cah', 0.89 ] ],
             [ [ ['offm_cah', 'offmt_cah'], -1.11 ],                       #http://www.ncbi.nlm.nih.gov/pubmed/19723630 (N3 vs N4)
               [ ['offh_cah', 'offht_cah'], 5.14 ],
               [ 'taummax_cah', 0.6 ],                                     #(tau_h_slow considered)                               
               [ 'tauhmax_cah', 1.48 ] ],
             [ [ ['offm_cah', 'offmt_cah'], -1.88 ],                       #http://www.ncbi.nlm.nih.gov/pubmed/19723630 (N5 vs N4)
               [ ['offh_cah', 'offht_cah'], 2.69 ],
               [ 'taummax_cah', 0.6 ],                                     #(tau_h_slow considered)                               
               [ 'tauhmax_cah', 1.35 ] ] ])
 #KCNMA1:
 MT.append([ [ [ 'offm_bk', [-13.8, 57.0] ],                               #http://www.ncbi.nlm.nih.gov/pubmed/16100257 (E912A, D916A, N918A, Q920A, D923A)
               [ ['ctm_bk', 'ctmmax_bk'], [0.51, 2.75] ] ],
             [ [ 'offm_bk', 49.2 ] ],                                      #http://www.ncbi.nlm.nih.gov/pubmed/24067659 (Slo1C-KvT, Slo1C-Kv-minT)
             [ [ 'offm_bk', 16 ],                                          #http://www.ncbi.nlm.nih.gov/pubmed/11880513 (hslo)
               [ ['ctm_bk', 'ctmmax_bk'], 5.91 ] ],
             [ [ 'offm_bk', [-5, 15] ],                                    #http://www.ncbi.nlm.nih.gov/pubmed/18719396 (e9alt, e9+e9alt) (voltage midpoint for [Ca]=5nM, tau for [Ca]=0.1uM)
               [ ['ctm_bk', 'ctmmax_bk'], 0.13 ] ],
             [ [ 'offm_bk', -58 ] ] ])                                     #http://www.ncbi.nlm.nih.gov/pubmed/11880513 (e20,e21,e22 vs ZERO)



             

 MTgenes = ['CACNB2', 'CACNA1D', 'CACNA1I', 'CACNA1S', 'ATP2A2', 'ATP2B2', 'NRGN', 'SCN1A', 'SCN9A', 'KCNS3', 'HCN1', 'KCNB1', 'CACNB2', 'KCNMA1'] 

 print "MT loaded successfully:"
 print MT
 return MT


def getdefvals():
 defVals = {'gamma_cad': 1.00,
           'taur_cad': 80.0,
           'depth_cad': 0.1,
           'cainf_cad': 1e-4,
           'offm_cah': -14.17,
           'offmt_cah': -26.31,
           'offh_cah': -22.63,
           'offht_cah': 19.73,
           'slom_cah': 9.76,
           'slomt_cah': 31.25,
           'sloh_cah': 6.6,
           'sloht_cah': 1/0.047,
           'taummax_cah': 0.97,
           'tauhmax_cah': 70,
           'offm_car': -23.0,
           'offmt_car': -23.0,
           'offh_car': -79.0,
           'offht_car': -79.0,
           'slom_car': 7.4,
           'slomt_car': 31.25,
           'sloh_car': 7.8,
           'sloht_car': 1/0.047,
           'taummax_car': 5.5,
           'tauhmax_car': 771.0,
           'eh_iH': -33.0,
           'off_iH': -91.0,
           'slo_iH': 6.0,
           'offt1_iH': -0.0,
           'offt2_iH': -0.0,
           't0_iH': 1.0/0.0003933,
           't1_iH': 1.0/0.0877,
           'slot1_iH': 1.0/0.0249,
           'slot2_iH': 1.0/0.062,
           'offh_kslow': -58,
           'sloh_kslow': 11.0,
           'a0_kslow': 1.0/0.0052,
           'a1_kslow': 1.0/0.01938,
           'a2_kslow': 1.0/0.0053,
           'offma_kslow': 11.1,
           'sloma_kslow': 13.1,
           'offmb_kslow': -1.27,
           'slomb_kslow': 71,
           'b0_kslow': 360.0,
           'b11_kslow': 1010.0,
           'b2_kslow': 23.7,
           'offht1_kslow': -54.0,
           'offht2_kslow': -75.0,
           'sloht2_kslow': 48.0,
           'bb0_kslow': 2350.0,
           'bb1_kslow': 1380.0,
           'bb2_kslow': -210.0,
           'offht3_kslow': 0.0,
           'offht4_kslow': 0.0,
           'sloht3_kslow': 1.0/0.01118,
           'sloht4_kslow': 1.0/0.0306,
           'offn_iA': -47.0,
           'slon_iA': 29.0,
           'offl_iA': -66.0,
           'slol_iA': 10,
           'offmt_iA': -71.0,
           'slomt_iA': 59.0,
           'taummin_iA': 0.34,
           'taumdiff_iA': 0.92,
           'offht_iA': -73.0,
           'sloht_iA': 23.0,
           'tauhmin_iA': 8.0,
           'tauhdiff_iA': 49.0,
           'offm_na': -38.0,
           'slom_na': 10.0,
           'offh_na': -66.0,
           'sloh_na': 6.0,
           'a1_na': 0.058,
           'a2_na': 0.114,
           'offmt_na': -36.0,
           'slomt_na': 28.0,
           'i1_na': 0.28,
           'i2_na': 16.7,
           'offht_na': -60.0,
           'sloht_na': 25.0,
           'b0inv_sk': 1.0/0.06,
           'offc_sk': (0.06/1.3e4)**0.25,
           'sloc_sk': 4.0,
           'zhalf_bk': 0.01,
           'offm_bk': -28.9,
           'slom_bk': 6.2,
           'ctm_bk': 0.000505,
           'ctmmax_bk': 1.0,
           'offmt1_bk': -86.4,
           'slomt1_bk': 10.1,
           'offmt2_bk': 33.3,
           'slomt2_bk': 10.0,
           'ctauz_bk': 1.0,
           'ch_bk': 0.085,
           'offh_bk': -32.0,
           'sloh_bk': 5.8,
           'cth_bk': 0.0019,
           'cthmax_bk': 1.0,
           'offht1_bk': -48.5,
           'sloht1_bk': 5.2,
           'offht2_bk': 54.2,
           'sloht2_bk': 12.9}
 print "defVals loaded successfully:"
 print defVals
 return defVals

def getsuffixes():
 return ['cad','cah','car','iH','kslow','iA','na','sk','bk']
def getsuffixexceptions():
 return [[],[],[],[['eh_iH','eh']],[],[],[],[],[]]

def getgenenames():
    return ['CACNA1C', 'CACNB2', 'CACNA1D', 'CACNA1I', 'CACNA1S', 'ATP2A2', 'ATP2B2', 'NRGN', 'SCN1A', 'SCN9A', 'KCNS3', 'KCNN3', 'HCN1', 'KCNB1', 'CACNB2', 'KCNMA1'] 
 #return ['CACNB2b', 'CACNA1D', 'CACNA1I', 'CACNA1S', 'ATP2A2', 'ATP2B2', 'NRGN', 'SCN1A', 'SCN9A', 'KCNS3', 'HCN1'] 

           
def getArticles():
 #       List of genes
 #         List of mutations
 Articles = [] 
 #CACNA1C:
 Articles.append([ 'http://www.ncbi.nlm.nih.gov/pubmed/19265197', 'http://www.ncbi.nlm.nih.gov/pubmed/19265197' ])
 #CACNB2:
 Articles.append([ 'http://www.ncbi.nlm.nih.gov/pubmed/19358333', 'http://www.ncbi.nlm.nih.gov/pubmed/7723731', 'http://www.ncbi.nlm.nih.gov/pubmed/19723630'])
 #CACNA1D:
 Articles.append([ 'http://www.ncbi.nlm.nih.gov/pubmed/21998309 and http://www.ncbi.nlm.nih.gov/pubmed/21998310 (42A)',
                   'http://www.ncbi.nlm.nih.gov/pubmed/21998309 and http://www.ncbi.nlm.nih.gov/pubmed/21998310 (43S)',
                   'http://www.ncbi.nlm.nih.gov/pubmed/20951705 and http://www.ncbi.nlm.nih.gov/pubmed/21054386 (CaV1.3 KO)'])
 #CACNA1I:
 Articles.append([ 'http://www.ncbi.nlm.nih.gov/pubmed/15254077'])
 #CACNA1S:
 Articles.append([ 'http://www.ncbi.nlm.nih.gov/pubmed/20861472', 'http://www.ncbi.nlm.nih.gov/pubmed/19134469'])
 #ATP2A:
 Articles.append([ 'http://www.ncbi.nlm.nih.gov/pubmed/10970890'])
 #ATP2B:
 Articles.append([ 'http://www.ncbi.nlm.nih.gov/pubmed/22789621', 'http://www.ncbi.nlm.nih.gov/pubmed/21232211', 'http://www.ncbi.nlm.nih.gov/pubmed/17234811'])
 #NRGN:
 Articles.append([ 'http://www.ncbi.nlm.nih.gov/pubmed/15564582'])
 #SCN1A:
 Articles.append([ 'http://www.ncbi.nlm.nih.gov/pubmed/18632931', 'http://www.ncbi.nlm.nih.gov/pubmed/18632931'])
 #SCN9A: (Probably no way to include the effects in this model..!)
 Articles.append([ 'http://www.ncbi.nlm.nih.gov/pubmed/22136189', 'http://www.ncbi.nlm.nih.gov/pubmed/18945915', 'http://www.ncbi.nlm.nih.gov/pubmed/16392115', 'http://www.ncbi.nlm.nih.gov/pubmed/15958509'])
 #KCNS3:
 Articles.append([ 'http://www.ncbi.nlm.nih.gov/pubmed/10484328'])
 #KCNN3:
 Articles.append([ 'http://www.ncbi.nlm.nih.gov/pubmed/14978258'])
 #HCN1:
 Articles.append([ 'http://www.ncbi.nlm.nih.gov/pubmed/17185333'])
 #KCNB1:
 Articles.append([ 'http://www.ncbi.nlm.nih.gov/pubmed/21455829 (T203K)', 'http://www.ncbi.nlm.nih.gov/pubmed/21455829 (T203D)', 'http://www.ncbi.nlm.nih.gov/pubmed/21455829 (S347K)',
                   'http://www.ncbi.nlm.nih.gov/pubmed/21455829 (S347D)', 'http://www.ncbi.nlm.nih.gov/pubmed/21455829 (T203W)', 'http://www.ncbi.nlm.nih.gov/pubmed/21455829 (S347W)'])
 #CACNB2 reprise:
 Articles.append([ 'http://www.ncbi.nlm.nih.gov/pubmed/19723630 (N1 vs N4)', 'http://www.ncbi.nlm.nih.gov/pubmed/19723630 (N3 vs N4)', 'http://www.ncbi.nlm.nih.gov/pubmed/19723630 (N5 vs N4)'])
 #KCNMA1:
 Articles.append([ 'http://www.ncbi.nlm.nih.gov/pubmed/16100257 (E912A, D916A, N918A, Q920A, D923A)', 'http://www.ncbi.nlm.nih.gov/pubmed/24067659 (Slo1C-KvT, Slo1C-Kv-minT)',
                   'http://www.ncbi.nlm.nih.gov/pubmed/11880513 (hslo)', 'http://www.ncbi.nlm.nih.gov/pubmed/2921853 (e9alt, e9+e9alt)', 'http://www.ncbi.nlm.nih.gov/pubmed/11880513 (e20,e21,e22 vs ZERO)'])
 return Articles
