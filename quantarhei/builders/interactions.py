# -*- coding: utf-8 -*-

import scipy.constants as const
from ..core.units import eps0_int
import numpy as np


def dipole_dipole_interaction(r1, r2, d1, d2, epsr):
    """ Calculates interaction between two dipoles
    
    
    Parameters
    ----------
    
    r1 : array like
        position of molecule 1
        
    r2 : array like 
        position of molecule 2
        
    d1 : array like
        dipole moment of the first molecule
        
    d2 : array like
        dipole moment of the second molecule
        
    epsr : float
        relative permitivity of the environment
        
    
    """
    
    R = r1 - r2
    RR = np.sqrt(np.dot(R,R))
    
    prf = 1.0/(4.0*const.pi*eps0_int)
    
    cc = (np.dot(d1,d2)/(RR**3)
        - 3.0*np.dot(d1,R)*np.dot(d2,R)/(RR**5))
    
    return prf*cc/epsr    
    

def dipole_dipole(center1,dipole1,center2,dipole2,*args):
    ''' Calculates interaction between two dipoles
    
    Pro dipoly v AU = naboj*Bohr a polohy b Bohrech vypocita
    interakcni energii v cm-1 
    
    '''
    r12=np.zeros(3)    
    r12=center2-center1
    R=np.sqrt(np.dot(r12,r12))
    
    d12=np.dot(dipole1,dipole2)
    dr1=np.dot(dipole1,r12)/R
    dr2=np.dot(dipole2,r12)/R    
    
    Edip_dip_Ha=(d12 - 3*dr1*dr2)/(R**3)   # interaction energy in hartree
    nma = 1.0
    Edip_dip_cm1=Edip_dip_Ha*nma.HaToInvcm

    if args:
        for a in args:
            if a=='Hartree':
                return Edip_dip_Ha
            elif a=='cm-1':
                return Edip_dip_cm1
                
    return Edip_dip_cm1
    
def Oscilator3D(rr1, bond1, AtType1, NMN1, TotDip1, rr2, bond2, AtType2,
                NMN2, TotDip2, *args):
    '''Interaction energy in inverse centimeters
    
    For position and dipole in atomic units (position in bohr radius and
    dipole in charge*bohr) calculates interaction energy in inverse centimeters
    
    Parameters
    ----------
    
    rr : float array, shape=(Natoms,3) 
        Atomic position in Bohr (AtomicUnits) (or atomic orbital position) 
        
    bond : int list
        Indexes of atoms with any kind of bond in between e.g bond
        between 1-6 and 3-5 bond=[[1,6],[3,5]] 
        
    AtType :
        Atom types (same length as rr) for example AtType=['C','N','C','C',...]
        
    NMN1 :
        Oscilator normal mode which will be used for calculation
        of interaction energy
        
    TotDip :
        Total transition dipole of molecule in charge*bohr (AtomicUnits)
        obtained from quantum chemistry or from experiment
    *args :
        Units of result = 'Hartree' or 'cm-1'. If not defined interaction 
        energy will be in inverse centimeters (cm-1) 
        
    '''
    
    scale_by_overlap=False    
#    use_model_carotenoid=False
#    
#    if use_model_carotenoid:
#        center=np.sum(rr1,0)
#        center=center/len(rr1)
#        VecX=rr1[len(rr1)//2+2,:]-rr1[len(rr1)//2,:]
#        VecY=rr1[len(rr1)//2+1,:]-rr1[len(rr1)//2,:]
#        rr1=pos.prepare_alkene(len(rr1),Position=center,vec_x=VecX,vec_y=VecY)
#        
#        center=np.sum(rr2,0)
#        center=center/len(rr2)
#        VecX=rr2[len(rr2)//2+2,:]-rr2[len(rr2)//2,:]
#        VecY=rr2[len(rr2)//2+1,:]-rr2[len(rr2)//2,:]
#        rr2=pos.prepare_alkene(len(rr2),Position=center,vec_x=VecX,vec_y=VecY)
#
    
    def molecule_osc_3D(rr,bond,factor,NMN,TotDip,*args):
        only_next_neighbour=False        
        TotDip_norm=np.sqrt(np.dot(TotDip,TotDip))        
        
        Ndip=len(bond)
        ro=np.zeros((Ndip,3),dtype='f8')
        do=np.zeros((Ndip,3),dtype='f8')
        
        # Place unity dipole moments in centers of all bonds
        for ii in range(Ndip):
            do[ii,:]=rr[bond[ii,1],:]-rr[bond[ii,0],:]
            norm=np.sqrt(np.dot(do[ii,:],do[ii,:]))
            # scaling dipole by defined value
            do[ii,:]=do[ii,:]/norm*factor[ii]           
            ro[ii,:]=(rr[bond[ii,1],:]+rr[bond[ii,0],:])/2
        
        # Calculate dipole-dipole interaction energy between all dipoles
        hh=np.zeros((Ndip,Ndip),dtype='f8')
        for ii in range(Ndip):
            for jj in range(ii+1,Ndip):
                hh[ii,jj]=dipole_dipole(ro[ii,:],do[ii,:],ro[jj,:],do[jj,:])
                hh[jj,ii]=hh[ii,jj]

        if only_next_neighbour:
            for ii in range(Ndip):
                for jj in range(ii+2,Ndip):
                    hh[ii,jj]=0.0
                    hh[jj,ii]=hh[ii,jj]       
        
        # Calculate normal modes (eigenvectors and eigenvalues of hh)
        # val= vector with eigenvalues and vec= matrix with 
        # eigenvectors in columns
        val,vec=np.linalg.eigh(hh) 
#        for ii in range(Ndip):
#            print('Eigenvalue: ',val[ii])
#            fig = plt.figure()
#            plt.plot(range(Ndip),vec[:,ii])
#            plt.show
        Dip=np.dot(vec[:,NMN],do)
        norm=np.sqrt(np.dot(Dip,Dip))
        
        for ii in range(Ndip):
            do[ii,:]=do[ii,:]*vec[ii,NMN]*TotDip_norm/norm
        
        return ro,do
            
            
    Ndip1=len(bond1)
    Ndip2=len(bond2)
            
    is_units=False
    if args:
        for a in args:
            is_units=True
            if a=='Hartree':
                units='Hartree'
            elif a=='cm-1':
                units='cm-1'    
    
    # Elementary dipole size is dependent on atom types (dipole in center of N-C bond should be 
    # different than dipole in center of C-C and also in centrer N-N). So far all dipoles are the 
    # same for all atom types. This should be changed and be dependent on AtType
    if scale_by_overlap:
        pass
#        factor1=np.zeros(Ndip1,dtype='f8')
#        factor2=np.zeros(Ndip2,dtype='f8')
#        coef=np.array([0.06899907,0.31642396,0.74430829])   # only for C-C bonds and for pz orbitals 6-31G basis
#        exp=np.array([7.86827235,1.88128854,0.54424926])    # only for C-C bonds and for pz orbitals 6-31G basis
#        r1=np.array([0.0,0.0,0.0])
#        for ii in range(Ndip1):
#            R=np.sqrt(np.dot(rr1[bond1[ii,1],:]-rr1[bond1[ii,0],:],rr1[bond1[ii,1],:]-rr1[bond1[ii,0],:]))
#            r2=np.array([R,0.0,0.0])            
#            Dip=qch.dipole_STO(r1,r2,coef,coef,exp,exp,[0,0,1],[0,0,1]) # only dependent on radial distance and not actual positions - it could be made more accurate this is only simplest approximation
#            factor1[ii]=np.sqrt(np.dot(Dip,Dip)) 
#        for ii in range(Ndip2):
#            R=np.sqrt(np.dot(rr2[bond2[ii,1],:]-rr2[bond2[ii,0],:],rr2[bond2[ii,1],:]-rr2[bond2[ii,0],:]))
#            r2=np.array([R,0.0,0.0])            
#            Dip=qch.dipole_STO(r1,r2,coef,coef,exp,exp,[0,0,1],[0,0,1]) # only dependent on radial distance and not actual positions - it could be made more accurate this is only simplest approximation            
#            factor2[ii]=np.sqrt(np.dot(Dip,Dip)) 
    else:
        # Possible scaling according to atom types (here everything is 1)
        factor1=np.ones(Ndip1,dtype='f8')
        factor2=np.ones(Ndip2,dtype='f8')
    
    ro1,do1=molecule_osc_3D(rr1,bond1,factor1,NMN1,TotDip1)
    ro2,do2=molecule_osc_3D(rr2,bond2,factor2,NMN2,TotDip2)
    
    #print('TotalDip1:',np.sum(do1,0))
    #print('TotalDip2:',np.sum(do2,0))    
    
    res=0.0
    if is_units:
        for ii in range(Ndip1):
            for jj in range(Ndip2):
                res += dipole_dipole(ro1[ii,:],do1[ii,:],ro2[jj,:],do2[jj,:],units)
    else:
        for ii in range(Ndip1):
            for jj in range(Ndip2):
                res += dipole_dipole(ro1[ii,:],do1[ii,:],ro2[jj,:],do2[jj,:])
            
    return res
    
def GuessBonds(rr, bond_length=4.0, **kwargs):
    ''' Function guesses pairs of atoms between which bond might occure.
    
    
    Parameters
    ----------
    
    rr : float array, shape=(Natoms,3)
        Atomic positions in Bohr (AtomicUnits) (or atomic orbital position)
        
    bond_length : float
        Bond length offset for pair of atoms in Bohr (AtomicUnits). 
        Should be dependent on atomic types (H-... bond reuire smaller offset
        than f.e. C-C bond)
        
        
    **kwargs :
        for example: AtType={['C','N','C','C',...]} (dictionary type) 
    
        example of **kwargs: 
            dictionary={'AtType': numpy.array([1,2,3,4,5])},
            dictionary['AtType'] -> array([1, 2, 3, 4, 5]),  
            list(dictionary.keys()) -> ['AtType']
            
    '''    
    
    Nat=len(rr)
    
    is_AtType=False
    for key in list(kwargs.keys()):
        if key=='AtType':
            AtType=kwargs['AtType']
            is_AtType=True
    
    Bonds=[]
    if is_AtType:
        for ii in range(Nat):
            for jj in range(ii+1,Nat):
                dr=rr[jj,:]-rr[ii,:]
                if AtType[ii]=='H' or AtType[jj]=='H':
                    if np.sqrt(np.dot(dr,dr))<bond_length/1.8:
                        Bonds.append([ii,jj])
                else:
                    if np.sqrt(np.dot(dr,dr))<bond_length:
                        Bonds.append([ii,jj])
    else:
        for ii in range(Nat):
            for jj in range(ii+1,Nat):
                dr=rr[jj,:]-rr[ii,:]
                if np.sqrt(np.dot(dr,dr))<bond_length:
                    Bonds.append([ii,jj])
    
    return np.array(Bonds,dtype='i8')
