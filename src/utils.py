# USER FUNCTIONS
import numpy as np
import math
import matplotlib.path as mplPath
# preview function
def xyz_preview(filename,header,atoms,species,x,y,z,writemode,*extra_points):
    """Generates xyz preview files. Arguments:
    - filename of the xyz file
    - second line header comment
    - number of atoms
    - species and coordinates arrays
    - write mode: 'w' first time; 'a' all the other times to append
    - optional: extra points to print as 3 elements lists, separated by commas"""
    with open(filename,mode=writemode) as fp:
        print(atoms+len(extra_points),file=fp)
        print(header,file=fp)
        for i in range(atoms):
            print(f'{species[i]}\t{x[i]:.8f}\t{y[i]:.8f}\t{z[i]:.8f}',file=fp)
        for i in extra_points:
            print(f'Xx\t{i[0]:.8f}\t{i[1]:.8f}\t{i[2]:.8f}',file=fp)
    fp.close()
# read xyz file function
def read_xyz(filename):
    """Read xyz file - typical format"""
    with open(filename,mode='r') as fp: # read file
        data=fp.read()
    fp.close() # done reading
    xyz0=data.split('\n') # tokenize with respect to new lines
    atoms=int(xyz0[0]) # get the number of atoms from first line
    xyz0[0:2]=[] # get rid of first two lines
    # serialize
    xyz=''.join([i+' ' for i in xyz0])
    # save species and coordinates
    species=xyz.split()[0:4*atoms:4]
    x=np.array([float(i) for i in (xyz.split()[1:4*atoms:4])]) # use list comprehensions for casting
    y=np.array([float(i) for i in (xyz.split()[2:4*atoms:4])])
    z=np.array([float(i) for i in (xyz.split()[3:4*atoms:4])])
    return [atoms,species,x,y,z]
# arithmetic center
def geo_center(x,y,z,*atoms):
    """Calculate the geometric center in 3D.
    Optional: define atomic indices as extra arguments in order to calculata the center using only these atoms."""
    collection=[int(i)-1 for i in atoms]
    if len(collection)!=0:
        return [x[collection].mean(),y[collection].mean(),z[collection].mean()]
    else:
        return [x.mean(),y.mean(),z.mean()]
# align to point
def align_to_point(x,y,z,X,Y,Z,*atoms):
    """Align molecule to given coordinates (X,Y,Z). Optional input the atomic ids to calculate the geometric center."""
    X_internal,Y_internal,Z_internal=geo_center(x,y,z,*atoms) # calculate center using unpacked atoms list
    x-=X_internal # center molecule
    y-=Y_internal
    z-=Z_internal
    x+=X # align
    y+=Y
    z+=Z
    X_internal,Y_internal,Z_internal=geo_center(x,y,z,*atoms) # recalculate center
    return [X_internal,Y_internal,Z_internal]
# rotation function
def rotate(ux,uy,uz,X,Y,Z,theta,x,y,z):
    """Rotation function
    https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle"""
    if isinstance(x,np.ndarray)==False: # treat the case when the input is a single point - not ndarray
        x=np.array([x])
        y=np.array([y])
        z=np.array([z])
    x_internal=np.zeros(x.size)
    y_internal=np.zeros(x.size)
    z_internal=np.zeros(x.size)
    u=np.array([ux,uy,uz],dtype=float) # define vector
    u/=np.linalg.norm(u) # normalize
    theta=theta*math.pi/180.0
    x_internal=(u[0]**2*(1.0-math.cos(theta))+math.cos(theta))*(x-X)+(u[0]*u[1]*(1.0-math.cos(theta))-u[2]*math.sin(theta))*(y-Y)+(u[0]*u[2]*(1.0-math.cos(theta))+u[1]*math.sin(theta))*(z-Z)+X
    y_internal=(u[0]*u[1]*(1.0-math.cos(theta))+u[2]*math.sin(theta))*(x-X)+(u[1]**2*(1.0-math.cos(theta))+math.cos(theta))*(y-Y)+(u[1]*u[2]*(1.0-math.cos(theta))-u[0]*math.sin(theta))*(z-Z)+Y
    z_internal=(u[0]*u[2]*(1.0-math.cos(theta))-u[1]*math.sin(theta))*(x-X)+(u[1]*u[2]*(1.0-math.cos(theta))+u[0]*math.sin(theta))*(y-Y)+(u[2]**2*(1.0-math.cos(theta))+math.cos(theta))*(z-Z)+Z
    return [x_internal,y_internal,z_internal]
# identify H atoms and crop
def crop_H_from_anchors(species,x,y,z,sites):
    """Crop H atoms connected to anchor sites and rescale the sites array"""
    H_thres=1.4 # bond length threshold
    atoms=len(species) # get number of atoms
    print_list=[i for i in range(atoms)] # indices list; non-zero will be mapped to H atoms for exclusion
    H_list=[i for i in range(atoms) if species[i]=='H'] # find the indices of H atoms in the molecule
    for i in sites: # nested loop: iterate over anchors and all H atoms
        for j in H_list:
            r=np.linalg.norm([x[j]-x[i],y[j]-y[i],z[j]-z[i]])
            if r<H_thres:
                print_list[j]=-1 # flag H atoms linked to an anchor
    map_sites=[-1 for i in range(atoms)] # since H atoms are deleted, we'll need to rescale atomic ids...
    j=-1
    for i in range(atoms):
        if print_list[i]!=-1:
            j+=1
            map_sites[i]=j
    map_dict={i:map_sites[i] for i in range(atoms)} # create a dictionary between old and new indices
    print_list=[print_list[i] for i in range(atoms) if print_list[i]!=-1] # filter out excluded H atoms
    return [atoms-len(sites),
            [species[i] for i in print_list],x[print_list],y[print_list],z[print_list],
            [map_dict[sites[i]] for i in range(len(sites))]]
# write PBC
def write_pbc(p1,p2,p3,species,x,y,z):
    """Read atomic info and three points to define the 2D cell - without the (0,0,0) origin - and apply PBC"""
    delta=1e-4 # tolerance
    cell=mplPath.Path(np.array([[0,0],p1,p2,p3,[0,0]])) # define area
    inside=cell.contains_points(tuple((x[i]+delta,y[i]+delta) for i in range(x.size)))
    inside_indices=[i for i in range(x.size) if inside[i]==True]
    # return atoms, species, x, y, z
    return [len(inside_indices),[species[i] for i in inside_indices],x[inside_indices],y[inside_indices],z[inside_indices]]
def replicate(primitive,A,B,A_rep,B_rep):
    repeat=A_rep*B_rep
    atoms=primitive.atoms*repeat
    species=[]
    x=np.zeros(atoms)
    y=np.zeros(atoms)
    z=np.zeros(atoms)
    for i in range(repeat): species+=primitive.species
    l=-1
    for i in range(A_rep):
        for j in range(B_rep):
            for k in range(primitive.atoms):
                l+=1
                x[l]=primitive.x[k]+i*A[0]+j*B[0]
                y[l]=primitive.y[k]+i*A[1]+j*B[1]
                z[l]=primitive.z[k]+i*A[2]+j*B[2]
    return [atoms, species, x, y, z]
def tri_wrap(a,b,c,x,y,z):
    
    def proj(u,v):
        return np.dot(u,v)/np.dot(u,u)
    
    a_hat=a/np.linalg.norm(a)
    b_hat=b/np.linalg.norm(b)
    c_hat=c/np.linalg.norm(c)
    
    projaa0=proj(a_hat,[1,0,0]);
    projab0=proj(a_hat,[0,1,0]);
    projac0=proj(a_hat,[0,0,1]);
    projba0=proj(b_hat,[1,0,0]);
    projbb0=proj(b_hat,[0,1,0]);
    projbc0=proj(b_hat,[0,0,1]);
    projca0=proj(c_hat,[1,0,0]);
    projcb0=proj(c_hat,[0,1,0]);
    projcc0=proj(c_hat,[0,0,1]);
    
    denom=-projac0*projbb0*projca0+projab0*projbc0*projca0+projac0*projba0*projcb0-projaa0*projbc0*projcb0-projab0*projba0*projcc0+projaa0*projbb0*projcc0
    
    a_norm_REF=math.sqrt(a[0]**2)
    b_norm_REF=math.sqrt(b[0]**2+b[1]**2)
    c_norm_REF=math.sqrt(c[0]**2+c[1]**2+c[2]**2)
        
    nx=np.zeros(x.size)
    ny=np.zeros(x.size)
    nz=np.zeros(x.size)
    
    for i in range(x.size):
        # projections
        nom1=-1.0*(projbc0*projcb0*x[i]-projbb0*projcc0*x[i]-projbc0*projca0*y[i]+projba0*projcc0*y[i]+projbb0*projca0*z[i]-projba0*projcb0*z[i])
        nom2=(projac0*projcb0*x[i]-projab0*projcc0*x[i]-projac0*projca0*y[i]+projaa0*projcc0*y[i]+projab0*projca0*z[i]-projaa0*projcb0*z[i])
        nom3=-1.0*(projac0*projbb0*x[i]-projab0*projbc0*x[i]-projac0*projba0*y[i]+projaa0*projbc0*y[i]+projab0*projba0*z[i]-projaa0*projbb0*z[i])
        
        # coordinates with respect to oblique system
        alpha=nom1/denom
        beta=nom2/denom
        gamma=nom3/denom
        
        # +/- direction
        if (alpha*a_hat[0]-0.0)*(a[0]-0.0)+(alpha*a_hat[1]-0.0)*(0.0-0.0)+(alpha*a_hat[2]-0.0)*(0.0-0.0) < 0.0 :
            a_sign=-1.0
        else:
            a_sign=1.0
        if(beta*b_hat[0]-0.0)*(b[0]-0.0)+(beta*b_hat[1]-0.0)*(b[1]-0.0)+(beta*b_hat[2]-0.0)*(0.0-0.0) < 0.0 :
            b_sign=-1.0
        else:
            b_sign=1.0
        if(gamma*c_hat[0]-0.0)*(c[0]-0.0)+(gamma*c_hat[1]-0.0)*(c[1]-0.0)+(gamma*c_hat[2]-0.0)*(c[2]-0.0) < 0.0 :
            c_sign=-1.0
        else:
            c_sign=1.0
    
        # image calculation
        a_norm=(alpha*a_hat[0])**2 + (alpha*a_hat[1])**2 + (alpha*a_hat[2])**2
        a_norm=math.sqrt(a_norm)
        b_norm=(beta*b_hat[0])**2 + (beta*b_hat[1])**2 + (beta*b_hat[2])**2
        b_norm=math.sqrt(b_norm)
        c_norm=(gamma*c_hat[0])**2 + (gamma*c_hat[1])**2 + (gamma*c_hat[2])**2
        c_norm=math.sqrt(c_norm)
        
        nx[i]=math.floor(a_norm/a_norm_REF)
        if a_sign < 0.0:
            nx[i]=-nx[i]-1
        ny[i]=math.floor(b_norm/b_norm_REF)
        if b_sign < 0.0:
            ny[i]=-ny[i]-1
        nz[i]=math.floor(c_norm/c_norm_REF)
        if c_sign < 0.0:
            nz[i]=-nz[i]-1
            
    return nx,ny,nz
def QE_scf(primitive,symmetry,system='test'):
    #system='test'
    masses={
        'H':1.00794,
        'B': 10.811, 
        'C': 12.0107, 
        'N': 14.0067, 
        'O': 15.999,
        'S': 32.065, 
        'F': 18.9984, 
        'Cl': 35.453, 
        'Br': 79.904, 
        'I': 126.90447,
        'Zn': 65.38,
    }
    # initialize: use if no initialization dict is givem
    pseudo='../pseudo/'
    ecutwfc=40
    ecutrho=480
    mixing_beta=0.7
    conv_thr='1d-8'    
    PBE_rrkjus_pseudos={
        'H':'H.pbe-rrkjus_psl.1.0.0.UPF',
        'B':'B.pbe-n-rrkjus_psl.1.0.0.UPF',
        'C':'C.pbe-n-rrkjus_psl.1.0.0.UPF',
        'N':'N.pbe-n-rrkjus_psl.1.0.0.UPF',
        'O':'O.pbe-n-rrkjus_psl.1.0.0.UPF',
        'S':'S.pbe-n-rrkjus_psl.1.0.0.UPF',
        'Zn':'Zn.pbe-dnl-rrkjus_psl.1.0.0.UPF'
    }
    # lattice types
    if symmetry is 'tetragonal':
        # tetragonal
        ibrav=6
        lattice_A=primitive.a[0]
        lattice_B=primitive.b[1]
        lattice_C=primitive.c[2]
        cosAB='0.0'
        cosAC='0.0'
        cosBC='0.0'
    # setup CONTROL
    control_str=f'&CONTROL\n calculation=\'scf\',\n restart_mode=\'from_scratch\',\n prefix=\'{system}\',\n pseudo_dir=\'{pseudo}\',\n outdir=\'./\',\n/'
    # species information
    unique_species=list(set(primitive.species))
    atom_types=len(unique_species)
    # setup SYSTEM
    system_str=f'&SYSTEM\n ibrav={ibrav},\n A={lattice_A:.6f},\n B={lattice_B:.6f},\n C={lattice_C:.6f},\n cosAB={cosAB},\n cosAC={cosAC},\n cosBC={cosBC},\n nat={primitive.atoms},\n ntyp={atom_types},\n ecutwfc={ecutwfc},\n ecutrho={ecutrho},\n/'
    # setup ELECTRONS
    electrons_str=f'&ELECTRONS\n mixing_beta={mixing_beta},\n conv_thr={conv_thr},\n/'
    # setup ATOMIC_SPECIES
    species_record=''
    for i in unique_species:
        species_record+=' '+i+' '+str(masses[i])+' '+PBE_rrkjus_pseudos[i]+'\n'
    species_record=species_record.rstrip('\n')
    atomic_str=f'ATOMIC_SPECIES\n{species_record}'
    # setup ATOMIC_POSITIONS
    positions_str='ATOMIC_POSITIONS (crystal)\n'
    if ibrav==6:
        for i in range(primitive.atoms):
            positions_str+=f'{primitive.species[i]}\t\
            {primitive.x[i]/lattice_A:.6f}\t{primitive.y[i]/lattice_B:.6f}\t{primitive.z[i]/lattice_C:.6f}\n'
    positions_str=positions_str.rstrip('\n')
    # setup KPOINTS
    kpoints_str='K_POINTS gamma'
#     print(control_str)
#     print(system_str)
#     print(electrons_str)
#     print(atomic_str)
#     print(positions_str)
#     print(kpoints_str)
    # write to file
    with open('scf.in',mode='w') as fp:
        print(control_str,file=fp)
        print(system_str,file=fp)
        print(electrons_str,file=fp)
        print(atomic_str,file=fp)
        print(positions_str,file=fp)
        print(kpoints_str,file=fp)
    fp.close()
