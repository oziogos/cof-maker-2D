from src.utils import *
import logging
import numpy as np
import math
import copy
import os
import shutil

class Common:
    
    def preview(self,mode='w',filename='preview.dat',header=' ',*points):
        """Generate an xyz preview file using the xyz_preview() function
        from the utilities.
        
        Arguments:
         - mode: 'w' (default) or 'a' for overwrite or append
         - filename: the default is 'preview.dat'
         - header: comment string for the second line of the xyz file
         - optional arguments: lists of cartesian points to be added; [x,y,z] format
           e.g. .preview('w','preview.dat',' ',[0,0,0],[0.5,0,0])
        
        If no arguments are given, the default action is to overwrite the
        'preview.dat' file.
        """
        xyz_preview(filename,header,self.atoms,self.species,self.x,self.y,self.z,mode,*points)
        logging.info(f'writing xyz file: {filename}; mode {mode}')
    
class Brick(Common):
    """This is a molecular fragment (aka Brick) for the creation for complex molecules.
    
    Arguments:
    - source: an appropriate dictionary for reading from file OR a Molecule object
    - sites: optional 'sites' list when reading Molecule object
    - centers: optional 'centers' list when reading Molecule object
    
    Basic atomic attributes:
     - atoms (int)
     - species (str list)
     - x, y, z (float ndarray): atomic coordinates
     
    Basic geometric attributes:
     - centers (int list): a list with atomic ids to define the geometric center
     - X, Y, Z (float): coordinates of the geometric center
     - normal (float ndarray): perpendicular vector for planar arrangements
       * automatically calculated using the 'centers' list as the cross product
         of the vectors starting from the geometric center and pointing to the
         second and first atoms defined in the 'centers' list
       * normalized vector
    """
    def __init__(self,source,sites=None,centers=None):
        # log
        logging.info('BRICK INITIALIZATION')
        # read dictionary entry
        if isinstance(source,dict)==True:
            # save the dict entry as a brick attribute
            self.source=source
            # read the docking sites
            self.sites=source['sites']
            # read the atomic ids for the geometric center or use 'sites' if
            # 'centers' isn't present in the dictionary
            if 'centers' in source: 
                self.centers=source['centers']
            else:
                self.centers=self.sites
            # read atomic info
            self.atoms,self.species,self.x,self.y,self.z=read_xyz(source['src'])
            # log
            logging.info(f'creating brick object from file: {self.source}')
        # read Molecule entry
        if isinstance(source,Molecule)==True:
            # create dict entry from Molecule object and save as attribute
            self.source={'src':f'internal Molecule object'}
            # read sites and centers
            if sites is None:
                self.sites=source.molecule_sites
            else:
                self.sites=sites
            self.source['sites']=self.sites
            if centers is None:
                self.centers=self.sites
            else:
                self.centers=centers
                self.source['centers']=centers
            # read atomic info
            self.atoms=source.atoms
            self.species=source.species
            self.x=source.x
            self.y=source.y
            self.z=source.z
            # log
            logging.info(f'creating brick object from molecule object: {self.source}')
        # rescale from atomic to python indices for 'sites'
        logging.info(f'sites in atomic indices: {self.sites}')
        self.sites=[int(i)-1 for i in self.sites]
        logging.info(f'sites in python indices: {self.sites}')
        # create the dictionary with H atoms linked to anchoring sites
        self.find_H()
        # initialize the docking axes and normal vector
        self.resolve_axes()
        # initialize the docking record 
        # this is an inter-brick connectivity record
        self.registry=[]
        # this is a dictionary listing the availability of each docking axis per docking site
        self.avail_dock={i:{j:True for j in range(len(self.axes[i]))} for i in range(len(self.sites))}
        # log
        logging.info(f'available docking axes per site: {self.avail_dock}')
        logging.info('END OF BRICK INITIALIZATION')
    
    def init_planar(self,site=None,direction=None,axis=0):
        """Optional initialization for planar molecules.
        
        Arguments:
        - site: selected site index
        - direction: alignment vector
        - axis: axis index belonging to the selected site
        
        When called without arguments, the brick's geometric center is aligned
        to (0,0,0) and an appropriate rotation aligns the normal vector with +z,
        thus aligning the planar molecule with the xy plane.
        
        If arguments are provided, the method requires a site index, a reference
        direction, and an axis index (default is 0). Then the defined axis of the
        selected site will be aligned with the given direction.
        
        This method is suited for aligning planar molecules, provided that the 
        docking axes are perpendicular to the normal vector.
        """
        # move molecular center to to (0,0,0)
        self.displace(0,0,0)
        # align with xy plane; skip if already aligned
        if abs(self.normal[2])<1.0:
            # this is the angle of the normal vector with +z
            theta=math.acos(self.normal[2])*180/math.pi
            # resolve rotation direction for +z alignment
            p=np.cross([0,0,1],self.normal)
            # align
            self.x,self.y,self.z=rotate(p[0],p[1],p[2],0,0,0,-theta,self.x,self.y,self.z)
            # update vectors
            self.resolve_axes()
        # conditional alignment of 'site' axis with given 'direction'
        if site is not None and direction is not None:
            # store vectors
            a=np.array(direction,dtype=float)
            b=self.axes[site][axis]
            # normalize input
            a/=np.linalg.norm(a)
            # products
            # TODO: handle the collinear and opposite vector cases
            dot=np.dot(a,b)
            c=np.cross(a,b)
            # align
            self.rotate(*c,-math.acos(dot)*180/math.pi,*self.resolve_site(site))
            # recenter to (0,0,0)
            self.displace(0,0,0)
            # update vectors
            self.resolve_axes()
        
    def find_H(self):
        """Resolve the H atoms linked to each docking site and store 
        atomic indices in the 'H_record' dictionary using python indices.
        """
        # bond length threshold
        H_thres=1.3
        # find the indices of H atoms in the molecule
        H_list=[i for i in range(self.atoms) if self.species[i]=='H']
        # initialize dictionary
        self.H_record={}
        # nested loop: iterate over anchors and all H atoms
        for i in self.sites:
            # initialize as a list
            occurrences=[]
            for j in H_list:
                # calc distance
                r=np.linalg.norm([self.x[j]-self.x[i],self.y[j]-self.y[i],self.z[j]-self.z[i]])
                # check distance
                if r<H_thres:
                    # augment the list
                    occurrences.append(j)
                    # overwrite entry in dictionary
                    self.H_record[i]=occurrences
        # log
        logging.info(f'H atoms linked to docking sites (python indices): {self.H_record}')
        
    def resolve_axes(self):
        """Resolve the docking axes per site using the information on the covalently linked
        H atoms. Store the normalized ndarray vectors in the 'axes' dictionary (sequential
        indexing starting from zero).
        
        For non-linear molecules (i.e. 'centers' list has more than two entries), calculate
        the normalized perpendicular vector 'normal' (ndarray) from the cross product of the
        vectors starting from the geometric center and pointing to the second and first atoms
        listed inside 'centers'.
        """
        # initialize dictionary
        self.axes={}
        for key,value in self.H_record.items():                        
            self.axes[self.sites.index(key)]=[np.array(
                [self.x[i]-self.x[key],
                 self.y[i]-self.y[key],
                 self.z[i]-self.z[key]]) for i in value]
            # normalize
            self.axes[self.sites.index(key)]=[i/np.linalg.norm(i) for i in self.axes[self.sites.index(key)]]
        # log
        logging.info(f'resolved docking axes per site: {self.axes}')
        # calculate center
        self.X,self.Y,self.Z=geo_center(self.x,self.y,self.z,*list(self.centers))
        # calculate normal vector
        # initialize
        self.normal=None
        # skip linear molecules!
        if len(self.centers)>2:
            self.normal=np.cross(
                [self.x[self.centers[1]-1]-self.X,self.y[self.centers[1]-1]-self.Y,self.z[self.centers[1]-1]-self.Z],
                [self.x[self.centers[0]-1]-self.X,self.y[self.centers[0]-1]-self.Y,self.z[self.centers[0]-1]-self.Z])
            # normalize
            self.normal/=np.linalg.norm(self.normal)
        # log
        logging.info(f'normal vector: {self.normal}')
    
    def displace(self,X,Y,Z,centers=None):
        """Move a brick.
        
        Arguments:
        - X, Y, Z : the point in space where the brick will be placed
        - centers: an optional list of atomic ids
        
        The method uses the align_to_point() function from the utilities. 
        If 'centers' is not defined, the method uses the brick's 'centers'
        list in order to calculate the geometric center and align it with
        (X,Y,Z). If a single atomic id is defined in 'centers', then the
        method will align the molecule by a rigid translation in order to
        bring the selected atom to (X,Y,Z).
        """
        if centers is None: centers=self.centers
        self.X,self.Y,self.Z=align_to_point(self.x,self.y,self.z,X,Y,Z,*centers)
        
    def rotate(self,ux,uy,uz,theta,X=None,Y=None,Z=None):
        """Rotate a brick.
        
        Arguments:
        - ux, uy, uz: the rotation direction
        - theta: the rotation angle abour (ux,uy,uz)
        - X, Y, Z: optional rotation center
        
        The method uses the rotate() function from the utilities in order to
        carry out a rigid body rotation. If X, Y, Z are not defined, the rotation
        takes place with respect to the brick's geometric center.
        
        After every rotation, the resolve_axes() method is called in order to
        recalculate all vectors AND the geometric center.
        """
        if X is None: X=self.X
        if Y is None: Y=self.Y
        if Z is None: Z=self.Z
        self.x,self.y,self.z=rotate(ux,uy,uz,X,Y,Z,theta,self.x,self.y,self.z)
        self.resolve_axes()
    
    def resolve_site(self,site):
        """Return the cartesian coordinates of the selected docking site."""
        return [self.x[self.sites[site]],self.y[self.sites[site]],self.z[self.sites[site]]]
    
    def __repr__(self):
        axes=''.join(
            f' {self.species[self.sites[key]]} {self.sites[key]+1}:\n {[list(j) for i,j in enumerate(value)]}\n' 
            for key,value in self.axes.items())
        avail_dock=''.join([f'  (site index {site}, axis {check}) is available for docking\n' 
                            for site,axis in self.avail_dock.items() for check in axis if axis[check]==True])
        if avail_dock=='': avail_dock='  saturated brick\n'
        if self.registry==[]:
            state='undocked brick\n'
        else:
            state=' part of molecule\n'
            # unfold registry
            for record in self.registry:
                state+=f'  linked via (site index {record[0]}, axis {record[0]} {self.axes[record[0]][record[1]]}) to element {record[2]}\n'
        return(f'Brick object\n'+
               f' src: {self.source["src"]}\n'+
               f' atoms: {self.atoms}\n'+
               f' sites: {[int(i)+1 for i in self.sites]}\n'+
               f' centers: {self.centers}\n'+
               f' axes:\n'+axes+
               f' normal: {self.normal}\n'+
               f' geo center: {[self.X,self.Y,self.Z]}\n'+
               f' state: {state}'+
               f'{avail_dock}')

class xyz(Common):
    
    def __init__(self,atoms,species,x,y,z,a=None,b=None,c=None):
        self.atoms=atoms
        self.species=species
        self.x=x
        self.y=y
        self.z=z
        self.a=a
        self.b=b
        self.c=c
    
class Molecule(Common):
    
    def __init__(self):
        # the molecule is a list of brick objects; initialize here
        self.element=[]      
        
    def add(self,new_brick):
        # adding a new brick to the molecule
        self.element.append(new_brick)

    def dock(self,parent,child,psite,csite,roll=0,paxis=0,caxis=0):
        # update the connectivity registry
        self.element[parent].registry.append([psite,paxis,child,csite,caxis])
        self.element[child].registry.append([csite,caxis,parent,psite,paxis])
        # update the docking record of parent and child bricks
        self.element[parent].avail_dock[psite][paxis]=False
        self.element[child].avail_dock[csite][caxis]=False
        # retrieve the coords of the parent molecule linking site and store in 'point'
        point=self.element[parent].resolve_site(psite)
        # align the linking site of the child molecule to 'point'
        self.element[child].displace(*point,[self.element[child].sites[csite]+1])
        # retrieve the linking axes
        a=self.element[parent].axes[psite][paxis]
        b=self.element[child].axes[csite][caxis]
        # check alignment state
        state='general'
        if abs(np.dot(a,b)-1.0) <= 1.0e-6:
            state='collinear'
        elif abs(np.dot(a,b)+1.0) <= 1.0e-6:
            state='opposite'
        # force opposite alignment 
        if state=='general':
            c=list(np.cross(a,b))
            self.element[child].rotate(*c,-math.acos(np.dot(a,b))*180/math.pi+180,*point)
        elif state=='collinear':
            c=list(self.element[parent].normal)
            self.element[child].rotate(*c,180,*point)
        # slide child molecule
        d=1.5
        direction=[a[0],a[1],a[2]]
        self.element[child].x+=d*direction[0]
        self.element[child].y+=d*direction[1]
        self.element[child].z+=d*direction[2]
        # roll! the normal vector is used for alignment purposes...
        if self.element[parent].normal is not None and self.element[child].normal is not None:
            a=self.element[parent].normal
            b=self.element[child].normal
            c=np.cross(a,b)
            dot=np.dot(c,np.array(direction))
            self.element[child].rotate(*direction,-np.sign(dot)*math.acos(np.dot(a,b))*180/math.pi+roll,*point)
        elif self.element[parent].normal is None and self.element[child].normal is not None:
            a=np.array([0,0,1])
            b=self.element[child].normal
            c=np.cross(a,b)
            dot=np.dot(c,np.array(direction))
            self.element[child].rotate(*direction,-np.sign(dot)*math.acos(np.dot(a,b))*180/math.pi+roll,*point)
        # recalculate center
        self.element[child].X,self.element[child].Y,self.element[child].Z=geo_center(
            self.element[child].x,self.element[child].y,self.element[child].z,*list(self.element[child].centers))
                
    def generate_atomic(self):
        self.atoms=0
        self.species=[]
        
        self.atomscount=[]
        for i in self.element:
            self.atomscount.append(self.atoms)
            self.atoms+=i.atoms
            self.species+=i.species
        self.x=np.hstack([i.x for i in self.element])
        self.y=np.hstack([i.y for i in self.element])
        self.z=np.hstack([i.z for i in self.element])
        
        erase_list=[]
        for elem in self.element:
            local=[]
            for site_index,dock_dict in elem.avail_dock.items():
                for H,state in dock_dict.items():
                    if state==False:
                        H_site=elem.sites[site_index]
                        local.append(elem.H_record[H_site][H])
            erase_list.append(local)        
                
        exclude=[j+self.atomscount[i] for i in range(len(self.atomscount)) for j in erase_list[i]]
        
        j=-1
        self.rescale=[]
        for i in range(self.atoms):
            if i not in exclude:
                j+=1
            self.rescale.append(j)
        
        self.atoms-=len(exclude)
        self.species=[self.species[i] for i in range(len(self.species)) if i not in exclude]
        self.x=np.delete(self.x,exclude)
        self.y=np.delete(self.y,exclude)
        self.z=np.delete(self.z,exclude)
        
        self.molecule_sites=[]
        i=-1
        for elem in self.element:
            i+=1
            for site,axis in elem.avail_dock.items():
                for check in axis:
                        if axis[check]==True:
                            self.molecule_sites.append(self.rescale[elem.sites[site]+1+self.atomscount[i]-1]+1)
        

    def preview_centers(self,filename='preview.dat'):
        self.generate_atomic()
        self.preview('w',filename)
        X=np.array([self.element[i].X for i in range(len(self.element))])
        Y=np.array([self.element[i].Y for i in range(len(self.element))])
        Z=np.array([self.element[i].Z for i in range(len(self.element))])
        xyz_preview(filename,' ',len(self.element),[''.join('Xx') for i in range(len(self.element))],X,Y,Z,'a')
        
def COF_4fold(fragments,core,linker1,linker2=None,roll1=0,roll2=0,repeat1=3,repeat2=3,z_vacuum=15,symmetry=None,save=False,outdir='output'):
    
    # treat string or brick input for the LINKERS
    if linker2 is None and isinstance(linker1,str)==True:
        linker2=linker1
    if linker2 is None and isinstance(linker1,Brick)==True:
        linker2=copy.deepcopy(linker1)
    
    if os.path.exists(outdir)==False:
        os.mkdir(outdir)
    
    directory=f'{outdir}/COF_4fold_{core}_{linker1}_{linker2}_{roll1}_{roll2}_{repeat1}_{repeat2}_{z_vacuum}_{symmetry}_{save}'
    if os.path.exists(directory)==False:
        os.mkdir(directory)
    
    # for debugging purposes...
    logging.getLogger().setLevel(logging.INFO)
    log = logging.getLogger()  # root logger
    for hdlr in log.handlers[:]:  # remove all old handlers
        log.removeHandler(hdlr)
    logging.basicConfig(filename=f'{directory}/log.dat', filemode='w', format='%(asctime)s - %(message)s')
    # initialize
    framework=Molecule()
    # place 4-fold core
    framework.add(Brick(fragments[core]))
    # align core to xy; align first site axis with +y
    framework.element[0].init_planar(0,[0,1,0])
    # dock first linker
    parent,child,psite,csite,roll=(0,1,0,0,roll1)
    if isinstance(linker1,str)==True:
        framework.add(Brick(fragments[linker1]));
    elif isinstance(linker1,Brick):
        framework.add(linker1);
    framework.dock(parent,child,psite,csite,roll)
    # calculate the A lattice vector
    A_lattice=framework.element[child].resolve_site(1)+framework.element[child].axes[1][0]*1.5-framework.element[parent].resolve_site(2)
    # dock second linker
    parent,child,psite,csite,roll=(0,2,1,0,roll2)
    if isinstance(linker2,str)==True:
        framework.add(Brick(fragments[linker2]));
    elif isinstance(linker2,Brick):
        framework.add(linker2);
    framework.dock(parent,child,psite,csite,roll)
    # calculate the B lattice vector
    B_lattice=framework.element[child].resolve_site(1)+framework.element[child].axes[1][0]*1.5-framework.element[parent].resolve_site(3)    
    # clear docking sites for periodic repetition
    framework.element[0].avail_dock[2][0]=False
    framework.element[0].avail_dock[3][0]=False
    framework.element[1].avail_dock[1][0]=False
    framework.element[2].avail_dock[1][0]=False
    framework.preview_centers(f'{directory}/preview.dat')
    # resolve the C cell constant using the vacuum value
    C_lattice=np.cross(B_lattice,A_lattice)
    C_lattice*=z_vacuum/np.linalg.norm(C_lattice)
    # alignments
    # align B_lattice with +x
    B_norm=B_lattice/np.linalg.norm(B_lattice)
    if abs(np.dot([1,0,0],B_norm))<1.0:
        u=np.cross([1,0,0],B_norm)
        theta=math.acos(np.dot([1,0,0],B_norm))*180/math.pi
        [framework.x,framework.y,framework.z]=rotate(u[0],u[1],u[2],0,0,0,-theta,framework.x,framework.y,framework.z)
        A_lattice=np.array([i[0] for i in rotate(u[0],u[1],u[2],0,0,0,-theta,*A_lattice)])
        B_lattice=np.array([i[0] for i in rotate(u[0],u[1],u[2],0,0,0,-theta,*B_lattice)])
        C_lattice=np.array([i[0] for i in rotate(u[0],u[1],u[2],0,0,0,-theta,*C_lattice)])
    # align A_lattice with xy
    if abs(np.dot(B_lattice/np.linalg.norm(B_lattice),A_lattice/np.linalg.norm(A_lattice)))<1.0:
        u=np.cross(B_lattice,A_lattice)
        u/=np.linalg.norm(u)
        theta=math.acos(np.dot([0,0,1],u))*180/math.pi
        if abs(np.dot([0,0,1],u))<1.0:
            p=np.cross([0,0,1],u)
            [framework.x,framework.y,framework.z]=rotate(p[0],p[1],p[2],0,0,0,-theta,framework.x,framework.y,framework.z)
            A_lattice=np.array([i[0] for i in rotate(p[0],p[1],p[2],0,0,0,-theta,*A_lattice)])
            B_lattice=np.array([i[0] for i in rotate(p[0],p[1],p[2],0,0,0,-theta,*B_lattice)])
            C_lattice=np.array([i[0] for i in rotate(p[0],p[1],p[2],0,0,0,-theta,*C_lattice)])
    # shift for vacuum
    framework.x+=C_lattice[0]/2
    framework.y+=C_lattice[1]/2
    framework.z+=C_lattice[2]/2
    # shift for core centering
    framework.x-=A_lattice[0]/2+B_lattice[0]/2+C_lattice[0]/2
    framework.y-=A_lattice[1]/2+B_lattice[1]/2+C_lattice[1]/2
    # repeat
    [atoms,species,x,y,z]=replicate(framework,A_lattice,B_lattice,repeat1,repeat2)
    xyz_preview(f'{directory}/repeat.dat',' ',atoms,species,x,y,z,'w',
                [0,0,0],list(A_lattice),list(B_lattice),list(C_lattice),list(A_lattice+B_lattice))
    # periodic check
    [nx,ny,nz]=tri_wrap(B_lattice,A_lattice,C_lattice,x,y,z)
    # create unit cell
    atoms_PBC=0
    for i in range(nx.size):
        if nx[i]==0 and ny[i]==0 and nz[i]==0:
            atoms_PBC+=1
    species_PBC=[]
    x_PBC=np.zeros(atoms_PBC)
    y_PBC=np.zeros(atoms_PBC)
    z_PBC=np.zeros(atoms_PBC)
    j=-1
    for i in range(nx.size):
        if nx[i]==0 and ny[i]==0 and nz[i]==0:
            j+=1
            x_PBC[j]=x[i]
            y_PBC[j]=y[i]
            z_PBC[j]=z[i]
            species_PBC.append(species[i])
    xyz_preview(f'{directory}/unit_cell.dat',' ',atoms_PBC,species_PBC,x_PBC,y_PBC,z_PBC,'w',
                [0,0,0],list(A_lattice),list(B_lattice),list(C_lattice),list(A_lattice+B_lattice))
    unit=xyz(atoms_PBC,species_PBC,x_PBC,y_PBC,z_PBC)
    # repeat for preview
    [atoms,species,x,y,z]=replicate(unit,A_lattice,B_lattice,repeat1,repeat2)
    xyz_preview(f'{directory}/unit_cell_repeat.dat',' ',atoms,species,x,y,z,'w',
                [0,0,0],list(A_lattice),list(B_lattice),list(C_lattice),list(A_lattice+B_lattice))
    # symmetrize - conditional
    # unit cell matrix
    a=np.array([[B_lattice[0],A_lattice[0],C_lattice[0]],
                [B_lattice[1],A_lattice[1],C_lattice[1]],
                [B_lattice[2],A_lattice[2],C_lattice[2]]])
    # fractional coordinates
    x_fraq=np.zeros(atoms_PBC)
    y_fraq=np.zeros(atoms_PBC)
    z_fraq=np.zeros(atoms_PBC)
    for i in range(atoms_PBC):
        # cartesian coords
        b=np.array([x_PBC[i],y_PBC[i],z_PBC[i]])
        [x_fraq[i],y_fraq[i],z_fraq[i]]=np.linalg.solve(a,b)
    # check symmetry input
    if symmetry is 'tetragonal_xy':
        # get mean tetragonal length
        L_xy=0.5*(np.linalg.norm(B_lattice)+np.linalg.norm(A_lattice))
        L_z=np.linalg.norm(C_lattice)
        # new coords
        x_PBC=x_fraq*L_xy
        y_PBC=y_fraq*L_xy
        z_PBC=z_fraq*L_z
        # new cell
        B_lattice=np.array([L_xy,0,0])
        A_lattice=np.array([0,L_xy,0])
        C_lattice=np.array([0,0,L_z])
        # overwrite data and generate preview files
        xyz_preview(f'{directory}/unit_cell_sym.dat',' ',atoms_PBC,species_PBC,x_PBC,y_PBC,z_PBC,'w',
                    [0,0,0],list(A_lattice),list(B_lattice),list(C_lattice),list(A_lattice+B_lattice))
        unit=xyz(atoms_PBC,species_PBC,x_PBC,y_PBC,z_PBC,B_lattice,A_lattice,C_lattice)
        # write QE
        QE_scf(unit,'tetragonal',f'{core}_{linker1}_{linker2}_{roll1}_{roll2}')
        # repeat
        [atoms,species,x,y,z]=replicate(unit,A_lattice,B_lattice,repeat1,repeat2)
        xyz_preview(f'{directory}/unit_cell_repeat_sym.dat',' ',atoms,species,x,y,z,'w',
                    [0,0,0],list(A_lattice),list(B_lattice),list(C_lattice),list(A_lattice+B_lattice))
    # backup files
    if linker1!=linker2:
        system_string=f'COF_4fold_{core}_{linker1}_{linker2}_roll1_{roll1}_roll2_{roll2}_{symmetry}'
    else:
        system_string=f'COF_4fold_{core}_{linker1}_roll1_{roll1}_roll2_{roll2}_{symmetry}'
    if save is True:
        if symmetry is not None:
            os.system(f'mv scf.in {directory}/scf_{system_string}.in')
            os.system(f'cp {directory}/unit_cell_sym.dat {directory}/unit_cell_sym_{system_string}.dat')
            os.system(f'cp {directory}/unit_cell_repeat_sym.dat {directory}/unit_cell_repeat_sym_{system_string}.dat')

    return directory

def generate_preview(out,outdir):
    preview=[]
    output=[i for i in os.listdir(outdir) if not i.startswith('.') and i!=out.split('/')[-1]]
    for i in output:
        fp=open(f'{outdir}/{i}/unit_cell_repeat_sym.dat',mode='r')
        instance=fp.readlines()
        atoms=instance[0].strip()
        preview.append(atoms)
        preview.append(i)
        for j in instance[2::]:
            preview.append(j.strip())
        fp.close()
    with open(out,mode='w') as fp:
        for i in preview:
            print(i,file=fp)

class COF_test():
    def __init__(self):
        self.output=[]
    def run_test(self,fragments,test_dict,target):
        for core in test_dict[target]['cores']:
            for linker,linker_value in test_dict[target]['linkers'].items():
                eval(test_dict[target]['function'])(fragments,core,linker,roll1=linker_value["roll1"],roll2=linker_value["roll2"],symmetry=linker_value["symmetry"],save=test_dict[target]["save"],outdir=test_dict[target]["outdir"])
        generate_preview(f'{test_dict[target]["outdir"]}/preview.dat',test_dict[target]["outdir"])
        self.output.append(f'{test_dict[target]["outdir"]}/preview.dat')
    def cleanup(self,test_dict,target):
        for i in self.output:
            if os.path.isdir(i) == True:
                shutil.rmtree(i)
            if os.path.isfile(i) == True:
                os.remove(i)
        shutil.rmtree(test_dict[target]["outdir"])