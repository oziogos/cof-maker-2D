# Each fragment has a 'src' string containing the relative path to the
# molecular fragment xyz file, a 'sites' list that includes the 
# atomic ids of the docking atoms, and an optional 'centers' list with
# the atomic ids used to calculate the geometric center of the molecule.
# 
# The ids for planar molecules should be entered in a clockwise fashion.
#
fragments={
    
    # --- CORES ----
    
    # porphyrins and phthalocyanines
    'porphyrin':{
        'src':'data/cores/porphyrin.xyz',
        'sites':[6,12,24,14]
    },
    'TIP-2H':{
        'src':'data/cores/TIP-2H.xyz',
        'sites':[23,32,1,12]
    },
    'Pc':{
        'src':'data/cores/phthalocyanine_BO2H.xyz',
        'sites':[53,61,58,56]
    },

    # --- LINKERS ---
    
    # rod-like
    'Ac2':{
        'src':'data/linkers/ethynylene.xyz',
        'sites':[4,1],
        'centers':[3,2]
    },
    'Ph2Ac2':{
        'src':'data/linkers/Ph2Ac2.xyz',
        'sites':[1,14],
        'centers':[3,5,1]
    },
    
    # PAHs
    'phenyl':{
        'src':'data/linkers/phenyl.xyz',
        'sites':[7,10],
        'centers':[10,12,8]
    },
    'phenanthrene':{
        'src':'data/linkers/phenanthrene.xyz',
        'sites':[1,13],
        'centers':[1,3,5]
    },
    'dibenzoanthracene':{
        'src':'data/linkers/dibenzoanthracene.xyz',
        'sites':[21,16],
        'centers':[9,4,7]
    },
    'pyrene':{
        'src':'data/linkers/pyrene.xyz',
        'sites':[1,13],
        'centers':[10,14,12]
    },
    'picene':{
        'src':'data/linkers/picene.xyz',
        'sites':[1,21],
        'centers':[12,14,10]
    },
    
    # articulated
    'biphenyl':{
        'src':'data/linkers/biphenyl.xyz',
        'sites':[4,10],
        # using the first phenyl ring for alignment reference
        'centers':[3,1,5]
    },
    'TBAA':{
        'src':'data/linkers/TBAA.xyz',
        'sites':[2,21],
        'centers':[16,37,34,35,18,15]
    },
    
    # PAHs with boronic acid
    'benzene_BO2H':{
        'src':'data/linkers/benzene_BO2H.xyz',
        'sites':[11,12],
        'centers':[3,5,1]
    },
    'anthracene_BO2H':{
        'src':'data/linkers/anthracene_BO2H.xyz',
        'sites':[17,20],
        'centers':[9,7,5]
    },
    
    # PAHs with sulfur
    'thienothiophene':{
        'src':'data/linkers/thienothiophene.xyz',
        'sites':[1,7],
        'centers':[5,2,8,6]
    },
    'BDT':{
        'src':'data/linkers/benzodithiophene.xyz',
        'sites':[8,12],
        'centers':[1,3,5]
    },
    'ADT':{
        'src':'data/linkers/ADT.xyz',
        'sites':[16,19],
        'centers':[9,7,5]
    },
    
}


# fragments={

    




#     'TIP':{
#         'src':'cores/aligned_preview_TIP_b3lyp_6-311gs_sym-011.xyz',
#         'sites':[32,1,12,23]
#     },
#     'TIP_v2':{
#         'src':'cores/aligned_preview_TIP_b3lyp_6-311gs_sym-011.xyz',
#         'sites':[31,2,13,22]
#     },
#     'TIP-2H':{
#         'src':'cores/preview_TIP-2H_b3lyp_6-311gs_sym-020.xyz',
#         'sites':[32,1,12,23]
#     },
#     'TIP-2H_v2':{
#         'src':'cores/preview_TIP-2H_b3lyp_6-311gs_sym-020.xyz',
#         'sites':[31,2,13,22]
#     },
#     'ethylene':{
#         'src':'ethylene.xyz',
#         'sites':[1,2],
#         'centers':[6,5,3,4]
#     },
#     'pyrene':{
#         'src':'linkers/aligned_pyrene-009.xyz',
#         'sites':[6,2,14,12],
#         #'sites':[24,7,8,23],
#         #'sites':[6,1,14,13],
#         'centers':(3,9,11,5)
#     },
#     'Ph2Ac2':{
#         'src':'linkers/aligned_Ph2Ac2-013.xyz',
#         'sites':[14,1],
#         'centers':[12,5,3,16]
#     },
#     'X':{
#         'src':'X.xyz',
#         'sites':[64,27,21,70]
#     },
#     'T':{
#         'src':'triangle.xyz',
#         'sites':[39,38,37]
#     },
#     'L':{
#         'src':'biphenyl.xyz',
#         'sites':[10,1],
#         'centers':[3,5,1]
#     },
#     'H':{
#         'src':'H.xyz',
#         'sites':[9,19,30,39]
#     },
#     'OTPB':{
#         'src':'cores/OTPB-010.xyz',
#         'sites':[21,13,1]
#     },
#     'anthraB':{
#         'src':'linkers/anthraB.xyz',
#         'sites':[17,20],
#         'centers':[9,7,5]
#     },
#     'TPA':{
#         'src':'cores/TPA.xyz',
#         'sites':[5,11,17]
#     },
#     'TPB':{
#         'src':'cores/TPB.xyz',
#         'sites':[10,16,22]
#     },
#     'TPBor':{
#         'src':'cores/TPBor.xyz',
#         'sites':[5,11,17]
#     },
#     'boroxine':{
#         'src':'cores/boroxine.xyz',
#         'sites':[3,5,1]
#     },
#     'fluorene':{
#         'src':'fluorene.xyz',
#         'sites':[2,11],
#         'centers':[4,8,9]
#     },
#     'PhPy':{
#         'src':'PhPy.xyz',
#         'sites':[21,6],
#         'centers':[10,12,8]
#     },
#     'cyclopent':{
#         'src':'cyclopent.xyz',
#         'sites':[3,5],
#         'centers':[3,5,1]
#     },
#     'benzB':{
#         'src':'linkers/benzB.xyz',
#         'sites':[11,12],
#         'centers':[3,5,1]
#     },
#     'PhNCPh':{
#         'src':'PhNCPh.xyz',
#         'sites':[7,1],
#         'centers':[13,16,9]
#     },
#     'pyrene-link':{
#         'src':'linkers/aligned_pyrene-009.xyz',
#         'sites':[1,13],
#         'centers':(3,9,11,5)
#     },
#     'TIP_non_planar':{
#         'src':'cores/TIP_non_planar.xyz',
#         'sites':[32,1,12,23]
#     },
#     'CTPB':{
#         'src':'cores/CTPB_6-311gs-009.xyz',
#         'sites':[17,3,11]
#     },
#     'TIP_non_planar_v2':{
#         'src':'cores/TIP_non_planar.xyz',
#         'sites':[31,2,13,22]
#     },
#     'Zn-TIP':{
#         'src':'cores/Zn-TIP.xyz',
#         'sites':[32,1,12,23]
#     },
#     'Zn-TIP_v2':{
#         'src':'cores/Zn-TIP.xyz',
#         'sites':[31,2,13,22]
#     },
# }