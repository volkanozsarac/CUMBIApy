import numpy as np

def circ_un(Ec,fpc,eco,espall,dels):
    """
    Details
    -------
    Mander model for unconfined normal weight concrete
    
    Parameters
    ----------
    Ec : float
        Young's modulus of concrete. (MPa)
    fpc : float
        Concrete compressive strength (MPa).
    eco : float
        concrete strain at peak concrete strength.
    espall : float
        spalling strain.
    dels : float
        strain increment.

    Returns
    -------
    ec : numpy.darray
        concrete strain vector.
    fc : numpy.darray
        concrete strength vector (MPa).

    """
    
    ec = np.arange(0,espall,dels)
    fc = ec * 0
    Esecu = fpc/eco
    ru = Ec/(Ec-Esecu)
    xu = ec/eco
    
    for i in range(len(ec)):
        if ec[i] < 2*eco:
            fc[i] = fpc*xu[i]*ru/(ru - 1 + xu[i]**ru)
        
        elif ec[i] >= 2*eco and ec[i] <= espall:
            fc[i] = fpc*(2*ru/(ru - 1 +2**ru))*(1-(ec[i] - 2*eco)/(espall - 2*eco))
         
        elif ec[i] >= espall:
            break
        
    return ec, fc

def circ_unlw(Ec,fpc,eco,espall,dels):
    """
    Details
    -------
    Mander model for unconfined lightweight concrete
    
    Parameters
    ----------
    Ec : float
        Young's modulus of concrete. (MPa)
    fpc : float
        Concrete compressive strength (MPa).
    eco : float
        concrete strain at peak concrete strength.
    espall : float
        spalling strain.
    dels : float
        strain increment.

    Returns
    -------
    ec : numpy.darray
        concrete strain vector.
    fc : numpy.darray
        concrete strength vector (MPa).

    """
    
    ec = np.arange(0,espall,dels)
    fc = ec * 0
    Esecu = fpc/eco
    ru  = Ec/(Ec-Esecu)
    xu  = ec/eco
    ru2 = Ec/(Ec-1.8*fpc/eco)
    
    for i in range(len(ec)):
        if ec[i] < eco:
            fc[i] = fpc*xu[i]*ru/(ru - 1 + xu[i]**ru)
        
        if ec[i] >= eco and ec[i] < 1.3*eco:
            fc[i] = fpc*xu[i]*ru2/(ru2-1+xu[i]**ru2)
            
        if ec[i] >= 1.3*eco and ec[i] <= espall:
            fc[i] = fpc*(1.3*ru2/(ru2-1+1.3**ru2))*(1-(ec[i]-1.3*eco)/(espall-1.3*eco))
            
        elif ec[i] >= espall:
            break
        
    return ec, fc

def circ_conf(Ec,Ast,Dh,clb,s,fpc,fy,eco,esm,D,dels,stype):
    """
    Details
    -------
    Mander model for confined normal weight concrete, 
    solid circular sections
    
    Parameters
    ----------
    Ec : float
        Young's modulus of concrete. (MPa)
    Ast : float
        total long. steel area (mm**2).
    Dh : float
        diameter of transverse reinf. (mm).
    clb : float
        cover to longitudinal bars (mm).
    s : float
        spacing of transverse steel (mm).
    fpc : float
        Concrete compressive strength (MPa).
    fy : float
        long steel yielding stress (MPa).
    eco : float
        concrete strain at peak concrete strength.
    esm : float
        max transv. steel strain (usually ~0.10-0.15).
    D : float
        section diameter (mm).
    dels : float
        strain increment.
    stype : str
        transv. reinf. type: 'spirals' or 'hoops'.

    Returns
    -------
    ec : numpy.darray
        concrete strain vector.
    fc : numpy.darray
        concrete strength vector (MPa).
    """
    sp  = s - Dh
    Ash = 0.25*np.pi*(Dh**2)
    ds   = D - 2*clb + Dh  # core diameter
    ros  = 4*Ash/(ds*s)    # transv. steel area ratio
    Ac   = 0.25*np.pi*(ds**2) # core area
    rocc = Ast/Ac          # long. steel area ratio
    if stype == 'spirals':
        ke   = (1-sp/(2*ds))/(1-rocc)
    elif stype == 'hoops':
        ke   = ((1-sp/(2*ds))/(1-rocc))**2
        
    fpl  = 0.5*ke*ros*fy

    fpcc = (-1.254 + 2.254*np.sqrt(1 + 7.94*fpl/fpc) - 2*fpl/fpc)*fpc
    # fpcc = (1+fpl/(2*fpc))*fpc
    
    ecc  = eco*(1 + 5*(fpcc/fpc-1))
    Esec = fpcc/ecc
    r    = Ec/(Ec-Esec)
    ecu  = 1.5*(0.004 + 1.4*ros*fy*esm/fpcc)
    
    ec = np.arange(0,ecu,dels)
    x  = (1/ecc)*ec;
    fc = fpcc*x*r/(r-1+x**r);
        
    return ec, fc

def circ_conflw(Ec,Ast,Dh,clb,s,fpc,fy,eco,esm,D,dels,stype):
    """
    Details
    -------
    Mander model for confined light weight concrete, 
    solid circular sections
    
    Parameters
    ----------
    Ec : float
        Young's modulus of concrete. (MPa)
    Ast : float
        total long. steel area (mm**2).
    Dh : float
        diameter of transverse reinf. (mm).
    clb : float
        cover to longitudinal bars (mm).
    s : float
        spacing of transverse steel (mm).
    fpc : float
        Concrete compressive strength (MPa).
    fy : float
        long steel yielding stress (MPa).
    eco : float
        concrete strain at peak concrete strength.
    esm : float
        max transv. steel strain (usually ~0.10-0.15).
    D : float
        section diameter (mm).
    dels : float
        strain increment.
    stype : str
        transv. reinf. type: 'spirals' or 'hoops'.

    Returns
    -------
    ec : numpy.darray
        concrete strain vector.
    fc : numpy.darray
        concrete strength vector (MPa).
    """
    sp  = s - Dh
    Ash = 0.25*np.pi*(Dh**2)
    ds   = D - 2*clb + Dh  # core diameter
    ros  = 4*Ash/(ds*s)    # transv. steel area ratio
    Ac   = 0.25*np.pi*(ds**2) # core area
    rocc = Ast/Ac          # long. steel area ratio
    if stype == 'spirals':
        ke   = (1-sp/(2*ds))/(1-rocc)
    elif stype == 'hoops':
        ke   = ((1-sp/(2*ds))/(1-rocc))**2
        
    fpl  = 0.5*ke*ros*fy

    fpcc = (1+fpl/(2*fpc))*fpc
    
    ecc  = eco*(1 + 5*(fpcc/fpc-1))
    Esec = fpcc/ecc
    r    = Ec/(Ec-Esec)
    ecu  = 1.5*(0.004 + 1.4*ros*fy*esm/fpcc)
    
    ec = np.arange(0,ecu,dels)
    x  = (1/ecc)*ec;
    fc = fpcc*x*r/(r-1+x**r);
        
    return ec, fc