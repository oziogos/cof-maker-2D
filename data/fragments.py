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
