# =============================================================================
# Start of User Input Data
# =============================================================================

# input data:

name = 'Outputs'      # identifies actual work, the output file will be name.xls

# section properties:
D       = 1067-22*4         # section diameter (mm)
clb     = 50-22             # cover to longitudinal bars (mm)

# member properties

L             = 1651        # member clear length (mm)
bending       = 'single'    # single or double
ductilitymode = 'uniaxial'  # biaxial or uniaxial

# reinforcement details:

nbl     = 30             # number of longitudinal bars
Dbl     = 32             # long. bar diameter (mm)   
Dh      = 16             # diameter of transverse reinf. (mm)
stype   = 'spirals'      # 'spirals' or 'hoops'*
s       = 100            # spacing of transverse steel (mm)*

# aplieed loads:

P      =  1325           # axial load kN (-) tension (+)compression

# material models (input the 'name' of the file with the stress-strain relationship
# to use the default models: Mander model for confined or unconfined  concrete type 'mc' or 'mu'.
# For lightweight confined concrete type 'mclw' 
# King model for the steel 'ks', Raynor model for steel 'ra':

confined   = 'mc'
unconfined = 'mu'
rebar      = 'ks'

# material properties 

fpc     = 1.3*28  # concrete compressive strength (MPa)
Ec      = 0       # concrete modulus of elasticity (MPa) or
                  # input 0 for automatic calculation using
                  # 5000(fpc)^0.5
eco     = 0.002   # unconfined strain (usually 0.002 for normal weight or 0.004 for lightweight)*
esm     = 0.10    # max transv. steel strain (usually ~0.10-0.15)*
espall  = 0.0064  # max uncon. conc. strain (usually 0.0064)

fy      = 1.1*414 # long steel yielding stress (MPa)
fyh     = 1.1*414 # transverse steel yielding stress (MPa)
Es      = 200000  # steel modulus of elasticity
fsu     = 1.4*fy  # long steel max stress (MPa)*
esh     = 0.008   # long steel strain for strain hardening (usually 0.008)*
esu     = 0.10    # long. steel maximum strain (usually ~0.10-0.15)*

Ey     =  350     # slope of the yield plateau (MPa)
C1     =  3.5     # defines strain hardening curve in the Raynor model [2-6]

# this information is used only if the default material models are selected

# strain limits for yield surface (interaction diagram);

csid = 0.004  # concrete
ssid = 0.015  # steel

# Deformation Limit States:

ecser = 0.004; esser = 0.015   # concrete (ecser) and steel (esser) serviceability strain
ecdam = 0.018; esdam = 0.060   # concrete (ecser) and steel (esser) damage control strain
                               # (to use the 2/3 of the ultimate concrete strain just tipe 'twth'
# temperature information (in case of freezing conditions)
temp  = 40           # temperature of the specimen in celsius
kLsp  = 0.022        # constant to calculate Lsp = kLsp*fy*Dbl
                     # (usually 0.022 at ambient temp. or 0.011 at -40C)

# analysis control parameters:
itermax    = 1000       # max number of iterations
ncl        = 40         # # of concrete layers
tolerance  = 0.001      # x fpc x Ag
dels       = 0.0001     # delta strain for default material models
                     
# =============================================================================
# End of User Input Data
# =============================================================================

import os
import sys
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import mander
import steel

models_path = 'models' # user specified models

# concrete modulus of elaticity
if Ec == 0:
    Ec = 5000*(fpc**0.5)                  

# tensile strength, considering temperature effect
if temp < 0:
    fct = (1-0.0105*temp)*0.56*(fpc**0.5)
elif temp >= 0:
    fct = 0.56*(fpc**0.5)

eccr = fct/Ec                   # concrete strain for cracking

Dsp = D - 2 * clb + Dh          # core diameter
dcore = clb - Dh * 0.5          # distance to the core
P      = P * 1000               # axial load in Newtons
Ast    = nbl * 0.25 * np.pi * (Dbl**2) # total long. steel area mm2

tcl = D / ncl                   # thickness of concrete layers
yl  = tcl * np.arange(1,ncl+1)    # border distance conc. layer

esser = -esser
esdam = -esdam

ecun, fcun = mander.circ_un(Ec,fpc,eco,espall,dels)
ec, fc = mander.circ_conf(Ec,Ast,Dh,clb,s,fpc,fy,eco,esm,D,dels,stype)

if unconfined == 'mu':
    ecun, fcun = mander.circ_un(Ec,fpc,eco,espall,dels)
elif unconfined == 'mc':
    ecun, fcun = mander.circ_conf(Ec,Ast,Dh,clb,s,fpc,fy,eco,esm,D,dels,stype)
elif unconfined == 'mclw':
    ecun, fcun = mander.circ_conflw(Ec,Ast,Dh,clb,s,fpc,fy,eco,esm,D,dels,stype)
else:
    AUX = np.loadxtxt(os.path.join(models_path,unconfined+'txt'))
    ecun = AUX[:,0]
    fcun = AUX[:,1]
    
if confined == 'mu':
    ec, fc = mander.circ_un(Ec,fpc,eco,espall,dels)
elif confined == 'mc':
    ec, fc = mander.circ_conf(Ec,Ast,Dh,clb,s,fpc,fy,eco,esm,D,dels,stype)
elif confined == 'mclw':
    ec, fc = mander.circ_conflw(Ec,Ast,Dh,clb,s,fpc,fy,eco,esm,D,dels,stype)
else:
    AUX = np.loadxtxt(os.path.join(models_path,confined+'txt'))
    ec = AUX[:,0]
    fc = AUX[:,1]    

if rebar == 'ks':
    es, fs = steel.king(Es,fy,fsu,esh,esu,dels)
elif rebar == 'ra':
    es, fs = steel.raynor(Es,fy,fsu,esh,esu,dels,C1,Ey)
else:
    AUX = np.loadxtxt(os.path.join(models_path,rebar+'txt'))
    es = AUX[:,0]
    fs = AUX[:,1]    
  
ecu = ec[-1]                         # maximum strain confined concrete
ecumander = ecu/1.5                  # ultimate strain predicted by the original mander model

if ecdam == 'twth':
    ecdam = ecumander

# vector with strains of confined concrete
ec = np.append(-1e10, ec)
ec = np.append(ec,ec[-1]+dels)
ec = np.append(ec,1e10)
# vector with stresses of confined concrete  
fc = np.append(0, fc)
fc = np.append(fc, 0)
fc = np.append(fc, 0)

# vector with strains of unconfined concrete
ecun = np.append(-1e10, ecun)
ecun = np.append(ecun,ecun[-1]+dels)
ecun = np.append(ecun,1e10)
# vector with stresses of unconfined concrete
fcun = np.append(0, fcun)
fcun = np.append(fcun, 0)
fcun = np.append(fcun, 0)

# maximum strain steel
esu = es[-1]
# vector with strains of the steel
es = np.append(es, es[-1]+dels)
es = np.append(es, 1e10)
# vector with stresses of the steel
fs = np.append(fs, 0)
fs = np.append(fs, 0)

esaux = np.zeros(len(es))
fsaux = 0*esaux

for i in range(len(es)):
    esaux[i] = es[len(es)-i-1]
    fsaux[i] = fs[len(fs)-i-1]

# vector with strains of the steel
es = np.append(-esaux, es[3:])
# vector with stresses of the steel
fs = np.append(-fsaux, fs[3:])

fig, ax = plt.subplots(figsize = (8, 6))
ax.fill_between(ec, 0, fc, edgecolor = 'black', facecolor='green', interpolate=True, alpha = 0.5, label = 'Confined Concrete')
ax.fill_between(ecun, 0, fcun, edgecolor = 'black', facecolor='blue', interpolate=True, alpha = 0.5, label = 'Unconfined Concrete')
ax.set_ylabel('Stress [MPa]')
ax.set_xlabel('Strain')
ax.set_xlim([0,1.05*ec[-3]])
ax.set_ylim([0,1.05*np.max(fc)])
ax.legend(loc = 'upper right')
ax.grid(True)
ax.set_title('Stress-Strain Relation for Concrete')

fig, ax = plt.subplots(figsize = (8, 6))
ax.fill_between(es, 0, fs, edgecolor = 'black', facecolor='Red', interpolate=True, alpha = 0.5, label = 'Reinforcing Steel')
ax.set_ylabel('Stress [MPa]')
ax.set_xlabel('Strain')
ax.set_xlim([1.05*es[3],1.05*es[-3]])
ax.set_ylim([1.05*fs[3],1.05*np.max(fs)])
ax.legend(loc = 'upper left')
ax.grid(True)
ax.set_title('Stress-Strain Relation for Reinforcing Steel')

# ============================== CONCRETE LAYERS ============================

# add layers to consider unconfined concrete
yl = np.append(np.append(yl, dcore), D-dcore)
yl.sort()
# confined concrete layers
yc = yl - dcore;
yc = yc[np.where((0 < yc)*(yc <= Dsp))[0]] 

# total area of each layer
Atemp = ((D/2)**2)*np.arccos(1-2*yl/D)-(D/2-yl)*((D*yl-yl**2)**0.5)
Atc = Atemp-np.append(0,Atemp[:-1])

# total area of each conf. layer
Atemp = ((Dsp/2)**2)*np.arccos(1-2*yc/Dsp)-(Dsp/2-yc)*((Dsp*yc-yc**2)**0.5)
Atcc = Atemp-np.append(0,Atemp[:-1])

conclay = []
k = 0
for i in range(len(yl)):
    if yl[i]<=dcore or yl[i]>D-dcore:
        conclay.append([Atc[i], 0])
        
    if yl[i]>dcore and yl[i]<=D-dcore:
        conclay.append([Atc[i]-Atcc[k], Atcc[k]])
        k += 1

conclay = np.asarray(conclay)
yl.shape = (len(yl),1)
ycenter = np.append(yl[0]/2, 0.5*(yl[:-1]+yl[1:])); ycenter.shape = (len(ycenter),1)

# [center_layer|A_uncon|A_conf|d_top_layer]
conclay = np.concatenate((ycenter,conclay,yl),axis = 1)

# ================================    REBARS     =====================================

Asb   = 0.25*np.pi*(Dbl**2)
r     = 0.5*(D-2*clb-Dbl)
theta = (2*np.pi/nbl)*np.arange(0,nbl)
distld = (0.5*(D-2*r)+r*np.sin(theta)*np.tan(0.5*theta))
distld.sort()         # y coordinate of each bar

# =============================== CORRECTED AREAS ======================================

# Substract the steel area
for i in range(nbl):
    aux = np.where(yl > distld[i])[0][0]
    conclay[aux,2] = conclay[aux,2] - Asb
    if conclay[aux,2] < 0:
        print('decrease # of layers')
        sys.exit()
        
# ============  Define vector (def) with the deformations in the top concrete ==================

df = np.arange(1e-4,20*ecu,1e-4)
    
if ecu > 0.0018:
    df = np.append(df[df<=16e-4],np.arange(18e-4,20*ecu,2e-4))
    
if ecu > 0.0025:
    df = np.append(df[df<=20e-4],np.arange(25e-4,20*ecu,5e-4))
    
if ecu > 0.006:
    df = np.append(df[df<=50e-4],np.arange(60e-4,20*ecu,10e-4))
    
if ecu > 0.012: 
    df = np.append(df[df<=100e-4],np.arange(120e-4,20*ecu,20e-4))

npts = len(df)

if P > 0:
    for k in range(npts):
        f1 = interp1d(ecun,fcun); temp1 = np.sum(f1(df[0]*np.ones(len(yl)))*conclay[:,1])
        f2 = interp1d(ec,fc); temp2 = np.sum(f2(df[0]*np.ones(len(yl)))*conclay[:,2])
        f3 = interp1d(es,fs); temp3 = np.sum(Asb*f3(df[0]*np.ones(len(distld))))
        compch = temp1 + temp2 + temp3
        if compch < P:
            df = df[1:]

npts = len(df)

# ===============ITERATIVE PROCESS TO FIND THE MOMENT - CURVATURE RELATION: ==============================

msg = 0                                 # stop conditions

curv   = [0]                            # curvatures
mom    = [0]                            # moments
ejen   = [0]                            # neutral axis
DF     = [0]                            # force eqilibrium
vniter = [0]                            # iterations
coverstrain = [0]                          
corestrain  = [0]                          
steelstrain = [0]

tol = tolerance*0.25*np.pi*(D**2)*fpc        # tolerance allowed
x = [D/2]                                    # location of N.A.
for k in range(npts):
    lostmomcontrol = max(mom)
    if mom[k]<(0.8*lostmomcontrol):
        msg = 4
        break
    
    F = 10*tol
    niter = -1                                  
    while abs(F)>tol:
        niter = niter + 1
        eec = (df[k]/x[niter])*(conclay[:,0]-(D-x[niter]))    # vector with the strains in the concrete
        ees = (df[k]/x[niter])*(distld-(D-x[niter]))          # vector with the strains in the steel

        fcunconf = interp1d(ecun,fcun)(eec)        # vector with stresses in the unconfined concr.           
        fcconf = interp1d(ec,fc)(eec)              # vector with stresses in the confinded concr.    
        fsteel = interp1d(es,fs)(ees)              # vector with stresses in the steel
        FUNCON = fcunconf*conclay[:,1]
        FCONF  = fcconf*conclay[:,2]
        FST    = Asb*fsteel;
        F      = np.sum(FUNCON) + np.sum(FCONF) + np.sum(FST) - P
        if F>0:
            x.append(x[niter] - 0.05*x[niter])
        
        elif F<0:
            x.append(x[niter] + 0.05*x[niter])

        if niter>itermax:
            msg = 3
            break

    cores = (df[k]/x[niter])*abs(x[niter]-dcore)
    TF = confined == unconfined
    if not TF:
        if cores >= ecu:
            msg = 1
            break

    elif TF:
        if df[k] >= ecu:
            msg = 1
            break
        
    if abs(ees[0]) > esu:
        msg = 2
        break

    ejen.append(x[niter])
    DF.append(x[niter])
    vniter.append(niter)
    temp = (sum(FUNCON*conclay[:,0]) + sum(FCONF*conclay[:,0]) + sum(FST*distld) - P*(D/2))/(10**6)
    mom.append(temp)

    if mom[k+1] < 0:
        mom[k+1] = -0.01*mom[k+1]

    curv.append(1000*df[k]/x[niter])
    coverstrain.append(df[k])
    corestrain.append(cores)
    steelstrain.append(ees[0])
    x[0] = x[niter]
    del x[1:]
    if msg!=0:
        break
    
Agross = 0.25*np.pi*(D**2)
AsLong = nbl*Asb
LongSteelRatio  = nbl*Asb/Agross
TransvSteelRatio = np.pi*Dh*Dh/(s*Dsp)
AxialRatio  = P/(fpc*Agross)

Mn     = interp1d(coverstrain,mom)(0.004)
esaux    = interp1d(mom,steelstrain)(Mn)
if esaux <- 0.015:
    Mn    = interp1d(steelstrain,mom)(-0.015)
    
cMn = interp1d(mom,ejen)(Mn)

fycurv = interp1d(steelstrain,curv)(-fy/Es)          # curvature for first yield
fyM    = interp1d(curv,mom)(fycurv)                  # moment for first yield

eqcurv = max((Mn/fyM)*fycurv,fycurv)

curvbilin = np.append(np.append(0,eqcurv), curv[-1])
mombilin  = np.append(np.append(0,Mn), mom[-1])

SectionCurvatureDuctility = curv[-1]/eqcurv

plt.figure()
plt.plot(curvbilin,mombilin,c='red')
plt.plot(curv,mom,c='blue',linestyle='--')
plt.grid(True)
plt.xlabel('Curvature [m$^{-1}$]',fontsize=16)
plt.ylabel('Moment [kN.m]',fontsize=16)
plt.title('Moment - Curvature Relation',fontsize=16)

# ===============FIND THE FORCE - DISPLACEMENT RELATION: ==============================

Lsp = np.zeros(len(steelstrain))
for j in range(len(steelstrain)):
    ffss = (-steelstrain[j]*Es)
    if ffss > fy:
        ffss = fy
    Lsp[j] = kLsp*ffss*Dbl     # Strain penetration length      
    
kkk = min(0.2*(fsu/fy-1),0.08)

if bending == 'single':
    Lp = max(kkk*L + kLsp*fy*Dbl,2*kLsp*fy*Dbl)           # Plastic hinge length
    LBE = L
elif bending == 'double':
    Lp = max(kkk*L + kLsp*fy*Dbl,2*kLsp*fy*Dbl)           # Plastic hinge length
    LBE = L/2
else:
    print('bending should be specified as single or double')
    sys.exit()


# Lp = 467.5;
# Lsp(:) = 0.5*Lp;

# Moyer - Kowalsky Buckling model
bucritMK = 0
CuDu   = curv/eqcurv

steelstrain = np.asarray(steelstrain)
if SectionCurvatureDuctility > 4:
    
    esgr4  = -0.5*interp1d(CuDu,steelstrain)(4)   # resdidual growth strain at ductility 4
    escc   = 3*((s/Dbl)**(-2.5))                  # allowable steel compression strain
    
    esgr = np.zeros(len(steelstrain))
    for i in range(len(steelstrain)):
        
        if CuDu[i] < 1:
            esgr[i] = 0
            
        elif CuDu[i] < 4 and CuDu[i] > 1:
            esgr[i] = (esgr4/4)*CuDu[i]

        elif CuDu[i] > 4:
            esgr[i] = -0.5*steelstrain[i]

    esfl = escc-esgr
    
    plt.figure()
    if -steelstrain[-1] >= esfl[-1]:
        bucritMK = 1
        fail     = esfl + steelstrain
        failCuDuMK = interp1d(fail,CuDu)(0)
        failesfl = interp1d(fail,esfl)(0)
        failss   = -interp1d(fail,steelstrain)(0)
        plt.plot(failCuDuMK,failss,'m', marker = '.', MarkerEdgeColor = 'k', MarkerFaceColor = 'g', MarkerSize = 12, label = None)
        plt.xlabel('Curvature Ductility',fontsize=16)
        plt.ylabel('Steel Tension Strain',fontsize=16)
        plt.legend()
        plt.title('Moyer - Kowalsky Buckling Model',fontsize=16)
        
    plt.plot(CuDu,-steelstrain, c = 'red', label = 'Column strain ductility behavior')
    plt.plot(CuDu,esfl,c = 'blue', linestyle='--', label = 'Flexural Tension Strain')
    plt.grid(True)
    plt.xlabel('Curvature Ductility',fontsize=16)
    plt.ylabel('Steel Tension Strain',fontsize=16)
    plt.legend()
    plt.title('Moyer - Kowalsky Buckling Model',fontsize=16)
        
# Berry - Eberhard Buckling model

bucritBE = 0

if AxialRatio >= 0.30:
    C0=0.006; C1=7.190; C2=3.129; C3=0.651; C4=0.227;        # model constants
else:
    C0=0.0010; C1=7.30; C2=1.30; C3=1.30; C4=3.00;           # model constants

# effective confinement ratio
roeff = TransvSteelRatio*fyh/fpc                             

# plastic rotation at the onset of bar buckling
rotb  = C0*(1+C1*roeff)*((1+C2*P/(Agross*fpc))**(-1))*(1+C3*LBE/D+C4*Dbl*fy/D)    
plrot = (curv-fycurv)*(Lp)/1000

plt.figure()
if max(plrot) > rotb:
    bucritBE = 1
    failBE = plrot - rotb
    failplrot  = interp1d(failBE,plrot)(0)
    failCuDuBE = interp1d(failBE,CuDu)(0)
    plt.plot(failCuDuBE,failplrot, linestyle='', marker = '.', MarkerEdgeColor = 'k', MarkerFaceColor = 'g', MarkerSize = 12, label = 'Buckling')

plt.plot(CuDu,rotb*np.ones(len(CuDu)), c = 'red', label = 'Plastic Rotation for Buckling')
plt.plot(CuDu,plrot,c = 'blue', linestyle='--', label = 'Plastic Rotation')    
plt.xlabel('Curvature Ductility',fontsize=16)
plt.ylabel('Plastic Rotation',fontsize=16)
plt.legend()
plt.title('Berry - Eberhard Buckling Model',fontsize=16)
plt.grid(True)
    
# Flexure deflection:
displf = np.zeros(len(curv))
mom = np.asarray(mom)
if bending == 'single':
    for i in range(len(curv)):
        if coverstrain[i] < eccr:
            displf[i] = curv[i]*((L/1000)**2)/3
            
        if coverstrain[i] > eccr and curv[i] < fycurv:
            displf[i] = curv[i] * (((L+Lsp[i])/1000)**2)/3

        if curv[i] >= fycurv:
            displf[i] = (curv[i]-fycurv*(mom[i]/fyM))*(Lp/1000)*((L+Lsp[i]-0.5*Lp)/1000) + \
                        (fycurv*(((L+Lsp[i])/1000)**2)/3)*(mom[i]/fyM)

    Force = mom/(L/1000)
    
elif bending == 'double':
    for i in range(len(curv)):
        if coverstrain[i] < eccr:
            displf[i] = curv[i]*((L/1000)**2)/6

        if coverstrain[i] > eccr and curv[i] < fycurv:
            displf[i] = curv[i] * (((L+2*Lsp[i])/1000)**2)/6

        if curv[i] >= fycurv:
            displf[i] = (curv[i]-fycurv*(mom[i]/fyM))*(Lp/1000)*((L+2*(Lsp[i]-0.5*Lp))/1000) + \
                        (fycurv*(((L+2*Lsp[i])/1000)^2)/6)*(mom[i]/fyM);

    Force = 2*mom/(L/1000)

else:
    print('bending should be specified as single or double')
    sys.exit()
    
# Shear deflection:

G     = 0.43*Ec
As    = 0.9*Agross
Ig    = np.pi*(D**4)/64
Ieff  = (Mn*1000/(Ec*(10**6)*eqcurv))*(10**12)

beta  = min(0.5+20*LongSteelRatio,1)

if bending == 'single':
    alpha = min(max(1,3-L/D),1.5)
    
elif bending == 'double':
    alpha = min(max(1,3-L/(2*D)),1.5)

Vc1   = 0.29*alpha*beta*0.8*(fpc**(1/2))*Agross/1000

kscr  = ((0.39*TransvSteelRatio)*0.25*Es*((0.8*D/1000)**2)/(0.25+10*(0.39*TransvSteelRatio)))*1000

if bending == 'single':
    ksg   = (G*As/L)/1000
    kscr  = (kscr/L)
    forcebilin = mombilin/(L/1000)
    
elif bending == 'double':    
    ksg   = (G*As/(L/2))/1000
    kscr  = (kscr/(L/2))
    forcebilin = 2*mombilin/(L/1000)
    
kseff = ksg*(Ieff/Ig)
aux = (Vc1/kseff)/1000
aux2 = 0
momaux = mom*1
displsh = np.zeros(len(curv))
for i in range(len(curv)):
    if momaux[i] <= Mn and Force[i] < Vc1:
        displsh[i] = (Force[i]/kseff)/1000

    if momaux[i] <= Mn and Force[i] >= Vc1:
        displsh[i] = ((Force[i]-Vc1)/kscr)/1000+aux

    if momaux[i] > Mn:
        momaux = 4*momaux
        aux3 = i - aux2
        aux2 = aux2 + 1
        displsh[i] = (displf[i]/displf[i-1])*displsh[i-1]

displ = displsh + displf

# bilinear approx:
dy1 = interp1d(curv,displ)(fycurv)
dy  = (Mn/fyM)*dy1
du  = displ[-1]
displbilin  = np.append(np.append(0,dy),du)
Dduct = displ/dy
DisplDuct = max(Dduct)
dy1f = interp1d(curv,displf)(fycurv)
dyf  = (Mn/fyM)*dy1f

# Shear Strength:
Vs    = (0.5*np.pi*(0.25*np.pi*(Dh**2))*fyh/np.tan(np.pi/6)*(D-clb+0.5*Dh-cMn)/s)/1000
Vsd   = (0.5*np.pi*(0.25*np.pi*(Dh**2))*fyh/np.tan((35/180)*np.pi)*(D-clb+0.5*Dh-cMn)/s)/1000
beta  = min(0.5+20*LongSteelRatio,1)
Dductf = displ/dyf

if bending == 'single':
    alpha = min(max(1,3-L/D),1.5)
    if P > 0:
        Vp = (P*(D-cMn)/(2*L))/1000
    else:
        Vp = 0

elif bending == 'double':
    alpha = min(max(1,3-L/(2*D)),1.5)
    if P > 0:
        Vp = (P*(D-cMn)/(L))/100
    else:
        Vp = 0

Vc = np.zeros(len(Dductf))
if ductilitymode == 'uniaxial':
    for i in range(len(Dductf)):
        Vc[i] = alpha*beta*min(max(0.05,0.37-0.04*Dductf[i]),0.29)*0.8*(fpc**(1/2))*Agross/1000

elif ductilitymode == 'biaxial':        
    for i in range(len(Dductf)):
        Vc[i] = alpha*beta*min(max(0.05,0.33-0.04*Dductf[i]),0.29)*0.8*(fpc**(1/2))*Agross/1000
        
Vcd = 0.862*Vc
Vpd = 0.85*Vp
V   = Vc + Vs + Vp
Vd  = 0.85 * (Vcd + Vsd + Vpd)
criteria = 1
if V[-1] < Force[-1]:
    failure   = V-Force
    faildispl = interp1d(failure,displ)(0)
    failforce = interp1d(displ,Force)(faildispl)
    failduct  = interp1d(displ,Dduct)(faildispl)
    failmom   = interp1d(displ,mom)(faildispl)
    failcurv  = interp1d(displ,curv)(faildispl)
    failCuDu  = interp1d(displ,CuDu)(faildispl)
    
    if bending == 'single':
        if faildispl <= 2*dy:
            criteria = 2
        elif faildispl < 8*dy:
            criteria = 3
        else:
            criteria = 4
    elif bending == 'double':
        if faildispl <= 1*dy:
            criteria = 2
        elif faildispl < 7*dy:
            criteria = 3
        else:
            criteria = 4
            
Ieq = (Mn/(eqcurv*Ec))/1000 # equivalent I for NLTHA
Bi = 1/(((mombilin[1])/(curvbilin[1]))/((mombilin[2]-mombilin[1])/(curvbilin[2]-curvbilin[1]))) # Bilinear factor

# Limit States:
displdam = 0; displser = 0; Dductdam = 0; Dductser = 0;
curvdam = 0;   curvser  = 0; CuDudam  = 0; CuDuser = 0;
coverstraindam = 0; coverstrainser = 0;
steelstraindam = 0; steelstrainser = 0;
momdam = 0; momser = 0; Forcedam = 0; Forceser = 0; 

if max(coverstrain) > ecser or max(abs(steelstrain)) > abs(esser):
    
    if max(coverstrain) > ecdam or max(abs(steelstrain)) > abs(esdam):

        displdamc = interp1d(coverstrain,displ,fill_value="extrapolate")(ecdam)
        displdams = interp1d(steelstrain,displ,fill_value="extrapolate")(esdam)
        displdam  = min (displdamc,displdams)
        Dductdam  = interp1d(displ,Dduct,fill_value="extrapolate")(displdam)
        curvdam   = interp1d(displ,curv,fill_value="extrapolate")(displdam)
        CuDudam   = interp1d(displ,CuDu,fill_value="extrapolate")(displdam)
        coverstraindam = interp1d(displ,coverstrain,fill_value="extrapolate")(displdam)
        steelstraindam = interp1d(displ,steelstrain,fill_value="extrapolate")(displdam)
        momdam   = interp1d(displ,mom,fill_value="extrapolate")(displdam)
        Forcedam = interp1d(displ,Force,fill_value="extrapolate")(displdam)

    displserc = interp1d(coverstrain,displ,fill_value="extrapolate")(ecser)
    displsers = interp1d(steelstrain,displ,fill_value="extrapolate")(esser)
    displser  = min(displserc,displsers)
    Dductser  = interp1d(displ,Dduct,fill_value="extrapolate")(displser)
    curvser   = interp1d(displ,curv,fill_value="extrapolate")(displser)
    CuDuser   = interp1d(displ,CuDu,fill_value="extrapolate")(displser)
    coverstrainser = interp1d(displ,coverstrain,fill_value="extrapolate")(displser)
    steelstrainser = interp1d(displ,steelstrain,fill_value="extrapolate")(displser)
    momser   = interp1d(displ,mom,fill_value="extrapolate")(displser)
    Forceser = interp1d(displ,Force,fill_value="extrapolate")(displser)
    
outputlimit = [coverstrainser,steelstrainser,momser,Forceser,curvser,CuDuser,displser,Dductser, \
                coverstraindam,steelstraindam, momdam, Forcedam, curvdam, CuDudam, displdam, Dductdam, \
                max(coverstrain), min(steelstrain), mom[-1], Force[-1], max(curv), max(CuDu), max(displ), max(Dduct)]
for i in range(len(outputlimit)):
    outputlimit[i] = float(outputlimit[i])

# Plot the capacity curves
fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
if criteria != 1 and (bucritMK == 1 and bucritBE==0): 
    buckldispl = interp1d(CuDu,displ)(failCuDuMK)
    bucklforce = interp1d(CuDu,Force)(failCuDuMK)
    ax.plot(faildispl,failforce, c = 'magenta', linestyle = '', marker = '.', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'shear failure')
    ax.plot(buckldispl,bucklforce, c = 'red', linestyle = '', marker = '*', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'buckling (M & K)') 

elif criteria != 1 and (bucritMK == 1 and bucritBE==1): 
    buckldispl = interp1d(CuDu,displ)(failCuDuMK)
    bucklforce = interp1d(CuDu,Force)(failCuDuMK)
    buckldisplBE = interp1d(CuDu,displ)(failCuDuBE)
    bucklforceBE = interp1d(CuDu,Force)(failCuDuBE)
    ax.plot(faildispl,failforce, c = 'magenta', linestyle = '', marker = '.', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'shear failure')    
    ax.plot(buckldispl,bucklforce, c = 'red', linestyle = '', marker = '*', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'buckling (M & K)')    
    ax.plot(buckldisplBE,bucklforceBE, c = 'red', linestyle = '', marker = 's', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'buckling (B & E)')    
    
elif criteria !=1 and (bucritMK == 0 and bucritBE==1): 
    buckldisplBE = interp1d(CuDu,displ)(failCuDuBE)
    ax.plot(faildispl,failforce, c = 'magenta', linestyle = '', marker = '.', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'shear failure')    
    ax.plot(buckldisplBE,bucklforceBE, c = 'red', linestyle = '', marker = 's', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'buckling (B & E)')    

elif criteria !=1 and (bucritMK == 0 and bucritBE==0):
    ax.plot(faildispl,failforce, c = 'magenta', linestyle = '', marker = '.', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'shear failure')    

elif criteria == 1 and (bucritMK == 1 and bucritBE==0):
    buckldispl = interp1d(CuDu,displ)(failCuDuMK)
    bucklforce = interp1d(CuDu,Force)(failCuDuMK)
    ax.plot(buckldispl,bucklforce, c = 'red', linestyle = '', marker = '*', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'buckling (M & K)')    

elif criteria == 1 and (bucritMK == 1 and bucritBE==1): 
    buckldispl = interp1d(CuDu,displ)(failCuDuMK)
    bucklforce = interp1d(CuDu,Force)(failCuDuMK)
    buckldisplBE = interp1d(CuDu,displ)(failCuDuBE)
    bucklforceBE = interp1d(CuDu,Force)(failCuDuBE)
    ax.plot(buckldispl,bucklforce, c = 'red', linestyle = '', marker = '*', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'buckling (M & K)')    
    ax.plot(buckldisplBE,bucklforceBE, c = 'red', linestyle = '', marker = 's', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'buckling (B & E)')    

elif criteria == 1 and (bucritMK == 0 and bucritBE==1): 
    buckldisplBE = interp1d(CuDu,displ)(failCuDuBE)
    bucklforceBE = interp1d(CuDu,Force)(failCuDuBE)
    ax.plot(buckldisplBE,bucklforceBE, c = 'red', linestyle = '', marker = 's', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=15, label = 'buckling (B & E)')   

elif criteria == 1 and (bucritMK == 0 and bucritBE==0): 
    pass

pointsdam = np.where(displ<=displdam)[0]
pointsser = np.where(displ<=displser)[0]
ax.plot(displ,Force, c = 'black', linestyle = '-', label = 'total response')
ax.plot(displbilin,forcebilin, c = 'blue', linestyle = '--', label = 'bilinear approximation')
ax.plot(displ,V, c = 'red', linestyle = ':', label = 'shear capacity (assessment)')
ax.plot(displ,Vd, c = 'magenta', linestyle = ':', label = 'shear capacity (design)')
ax.fill_between(displ,0,Force, edgecolor = 'black', facecolor=(0.4, 0.8, 0.6), interpolate=True, alpha = 0.4, label = 'ultimate zone')
ax.fill_between(displ[pointsdam],0,Force[pointsdam], edgecolor = 'black', facecolor=(0.4, 0.6, 0.6), interpolate=True, alpha = 0.4, label = 'damage control zone')
ax.fill_between(displ[pointsser],0,Force[pointsser], edgecolor = 'black', facecolor=(0.4, 0.4, 0.6), interpolate=True, alpha = 0.4, label = 'serviceability zone')
ax.grid(True)
ax.set_xlabel('Displacement [m]', fontsize=16)
ax.set_ylabel('Force [kN]', fontsize=16)
ax.set_title('Force - Displacement Relation', fontsize=16)
plt.legend(bbox_to_anchor=(1.02, 0.5), loc='center left', frameon=False)
ylims = ax.get_ylim(); ax.set_ylim((0,ylims[1]))
xlims = ax.get_xlim(); ax.set_xlim((0,xlims[1]))

# ==========================================================================

output = np.asarray([np.asarray(coverstrain), np.asarray(corestrain), np.asarray(ejen), steelstrain, \
          mom, np.asarray(curv), Force, displsh, displf, displ, V, Vd]).T
    
outputbilin = np.asarray([curvbilin, mombilin, displbilin, forcebilin]).T

Acore = 0.25*np.pi*(Dsp**2)

# compression force for yield surface
PCid = interp1d(ec,fc)(csid)*(Acore-AsLong)+interp1d(ecun,fcun)(csid)*(Agross-Acore)+AsLong*interp1d(es,fs)(csid)
# tensile force for yield surface
PTid = AsLong*interp1d(es,fs)(ssid)

fid = open(name +'.xlsx','w')

fid.write('Circular Section\n\n')
if confined == 'mclw':
    fid.write('Lightweight concrete\n')
else:
    fid.write('Normalweight concrete\n')

fid.write('Diameter:\t%5.1f\tmm\n' % D)
fid.write('Cover to longitudinal bars:\t%.1f\tmm\n' % clb)
fid.write('Number of longitudinal bars:\t%.0f\t\n' % nbl)
fid.write('Diameter of longitudinal bars:\t%.1f\tmm\n' % Dbl)
fid.write('Diameter of transverse steel:\t%.1f\tmm\n' % Dh)
fid.write('Spacing of transverse steel:\t%.1f\tmm\n' % s)

if stype == 'spirals':
    fid.write('Type of tranverse reinforcement: Spirals\n')
elif stype == 'hoops':
    fid.write('Type of tranverse reinforcement: Hoops\n')

fid.write('Axial load:\t%8.2f\tkN\n' % float(P/1000))
fid.write('Concrete compressive strength:\t%3.2f\tMPa\n' % fpc)
fid.write('Long steel yielding stress:\t%4.2f\tMPa\n' % fy)
fid.write('Long steel max. stress:\t%4.2f\tMPa\n' % max(fs))
fid.write('Transverse steel yielding stress:\t%4.2f\tMPa\n' % fyh)
fid.write('Member Length:\t%5.1f\tmm\n' % L)

if bending == 'single':
    fid.write('Single Bending\n')
elif bending == 'double':
    fid.write('Double Bending\n')

if ductilitymode == 'uniaxial':
    fid.write('Uniaxial Bending\n')    
elif ductilitymode == 'biaxial':
    fid.write('Biaxial Bending\n')        

fid.write('Longitudinal Steel Ratio:\t%1.3f\n' % LongSteelRatio)
fid.write('Transverse Steel Ratio:\t%1.3f\n' % TransvSteelRatio)
fid.write('Axial Load Ratio:\t%1.3f\n\n' % AxialRatio)

fid.write('Cover\tCore\tN.A\tSteel\tMoment\tCurvature\tForce\tSh displ.\tFl displ.\tTotal displ.\tShear(assess.)\tShear(design)\n')
fid.write('Strain\tStrain\t[mm]\tStrain\t[kN.m]\t[1/m]\t[kN]\t[m]\t[m]\t[m]\t[kN]\t[kN]\n')
for i in range(output.shape[0]):
    fid.write('%1.5f\t%1.5f\t%4.2f\t%1.5f\t%8.2f\t%1.5f\t%8.2f\t%1.5f\t%1.5f\t%1.5f\t%8.2f\t%8.2f\n' % tuple(output[i,:]))
fid.write('\n')
fid.write('Bilinear Approximation:\n')
fid.write('Curvature\tMoment\tDispl.\tForce\n')
fid.write('[1/m]\t[kN.m]\t[m]\t[kN]\n')
for i in range(outputbilin.shape[0]):
    fid.write('%1.5f\t%8.2f\t%1.5f\t%8.2f\n' % tuple(outputbilin[i,:]))
fid.write('\n')

if msg == 1:
   fid.write('*** Concrete strain exceeds maximum ***')
elif msg == 2:
   fid.write('*** Steel strain exceeds maximum ***')
elif msg == 3:
   fid.write('*** Number of iteration exceeds maximum ***')
elif msg == 4:
   fid.write('*** Excessive loss of strength ***')

fid.write('\nMoment for First Yielding:\t%8.2f\tkN.m\n' % fyM)
fid.write('Curvature for First Yielding:\t%1.5f\t1/m\n' % fycurv)
fid.write('Potential Section Nominal Moment:\t%8.2f\tkN.m\n' % Mn)
fid.write('Equivalent Curvature:\t%1.5f\t1/m\n' % eqcurv)
fid.write('Potential Section Curvature Ductility:\t%3.2f\n' % SectionCurvatureDuctility)
fid.write('Potential Displacement Ductility:\t%3.2f\n' % DisplDuct)
fid.write('\n')

if criteria == 1:
    fid.write('*** flexural failure ***')
elif criteria == 2:
    fid.write('*** brittle shear failure ***')
    fid.write('\nDisplacement for Shear Failure:\t%1.5f\tm\n' % faildispl)
    fid.write('Displacement Ductility at Shear Failure:\t%8.2f\n' % failduct)
    fid.write('Force for Shear Failure:\t%8.2f\tkN\n' % failforce)
    fid.write('Curvature for Shear Failure:\t%1.5f\t1/m\n' % failcurv)
    fid.write('Curvature Ductility at Shear Failure:\t%8.2f\n' % failCuDu)
    fid.write('Moment for Shear Failure:\t%8.2f\tkN.m\n' % failmom)
        
elif criteria == 3:
    fid.write('*** shear failure at some ductility ***')
    fid.write('\nDisplacement for Shear Failure:\t%1.5f\tm\n' % faildispl);
    fid.write('Displacement Ductility at Shear Failure:\t%8.2f\n' % failduct)
    fid.write('Force for Shear Failure:\t%8.2f\tkN\n' % failforce)
    fid.write('Curvature for Shear Failure:\t%1.5f\t1/m\n' % failcurv)
    fid.write('Curvature Ductility at Shear Failure:\t%8.2f\n' % failCuDu)
    fid.write('Moment for Shear Failure:\t%8.2f\tkN.m\n' % failmom)
        
elif criteria == 4:
    fid.write('*** ductil shear failure ***')
    fid.write('\nDisplacement for Shear Failure:\t%1.5f\tm\n' % faildispl)
    fid.write('Displacement Ductility at Shear Failure:\t%8.2f\n' % failduct)
    fid.write('Force for Shear Failure:\t%8.2f\tkN\n' % failforce)
    fid.write('Curvature for Shear Failure:\t%1.5f\t1/m\n' % failcurv)
    fid.write('Curvature Ductility at Shear Failure:\t%8.2f\n' % failCuDu)
    fid.write('Moment for Shear Failure:\t%8.2f\tkN.m\n' % failmom)
        
if bucritMK == 1:
    bucklDd =   interp1d(CuDu,Dduct)(failCuDuMK)
    bucklcurv = interp1d(CuDu,curv)(failCuDuMK)
    bucklmom  = interp1d(CuDu,mom)(failCuDuMK)
    fid.write('Moyer - Kowalsky buckling model:\n')
    fid.write('\nCurvature Ductility for Buckling:\t%8.2f\n' % failCuDuMK)
    fid.write('Curvature at Buckling:\t%3.5f\tm\n' % bucklcurv)
    fid.write('Displacement Ductility at Buckling:\t%8.2f\n' % bucklDd)
    fid.write('Displacement at Buckling:\t%3.5f\tm\n' % buckldispl)
    fid.write('Force for Buckling:\t%8.2f\tkN\n' % bucklforce)
    fid.write('Moment for Buckling:\t%8.2f\tkN\n' % bucklmom)
    
if bucritBE == 1:
    bucklDdBE = interp1d(CuDu,Dduct)(failCuDuBE)
    bucklcurvBE = interp1d(CuDu,curv)(failCuDuBE)
    bucklmomBE  = interp1d(CuDu,mom)(failCuDuBE)
    fid.write('Berry - Eberhard buckling model:\n')
    fid.write('\nCurvature Ductility for Buckling:\t%8.2f\n' % failCuDuBE)
    fid.write('Curvature at Buckling:\t%3.5f\tm\n' % bucklcurvBE) 
    fid.write('Displacement Ductility at Buckling:\t%8.2f\n' % bucklDdBE)
    fid.write('Displacement at Buckling:\t%3.5f\tm\n' % buckldisplBE)
    fid.write('Force for Buckling:\t%8.2f\tkN\n' % bucklforceBE)
    fid.write('Moment for Buckling:\t%8.2f\tkN\n' % bucklmomBE)
    
fid.write('\n')
fid.write('*** Potential Deformation Limit States (serviceability/damage control/ultimate) ***\n')
fid.write('Cover\tSteel\tMoment\tForce\tCurvature\tCurvature\tDisplacement\tDisplacement\n')
fid.write('Strain\tStrain\t[kN.m]\t[kN]\t[1/m]\tDuctility\t[m]\tDuctility\n')
for i in range(1,4):
    fid.write('%1.5f\t%1.5f\t%8.2f\t%8.2f\t%1.5f\t%3.2f\t%2.5f\t%3.2f\n' % tuple(outputlimit[(i-1)*8:(8*i)]))
fid.write('\nDeformation Limit States Citeria :\n')
fid.write('serviceability concrete strain:\t%1.4f\n' % ecser)
fid.write('serviceability steel strain:\t%1.4f\n' % esser)
fid.write('damage control concrete strain:\t%1.4f\n' % ecdam)
fid.write('damage control steel strain:\t%1.4f\n' % esdam)

if confined == 'mc':
    fid.write('\nOriginal Mander Model Ultimate Concrete Strain:\t%1.4f\n' % ecumander)

fid.write('\nfor NLTHA:\n')
fid.write('E: \t%10.2f \tPa\n' % float(Ec*(10**6)))
fid.write('G: \t%10.2f \tPa\n' % float(G*(10**6)))
fid.write('A: \t%10.4f \tm2\n' % float(Agross/(10**6)))
fid.write('I: \t%10.6f \tm4\n' % float(Ieq))
fid.write('Bi-Factor: \t%1.3f\n' % float(Bi))
fid.write('Hinge Length: \t%1.3f \tm\n' % float(Lp/1000))
fid.write('Tension Yield: \t%10.2f \tN\n' % float(PTid))
fid.write('Compression Yield: \t%10.2f \tN\n' % float(PCid))
fid.write('Moment Yield: \t%10.2f \tN-m\n' % float(Mn*1000))

# vector with axial loads for interaction diagram
temp1 = np.append(np.arange(-0.90*PTid,0,0.30*PTid),0).tolist()
temp2 = np.append(np.arange(0.05*fpc*Agross,0.6*PCid,0.05*fpc*Agross),0.6*PCid).tolist()
PP = temp1 + temp2 +[0.7*PCid] + [0.8*PCid] + [0.9*PCid]
nPP = len(PP)
Mni = []
eci = []
esi = []
Msgs = []

for i in range(nPP):
    
    df = np.arange(1e-4,20*ecu,1e-4)
        
    if ecu > 0.0018:
        df = np.append(df[df<=16e-4],np.arange(18e-4,20*ecu,2e-4))
        
    if ecu > 0.0025:
        df = np.append(df[df<=20e-4],np.arange(25e-4,20*ecu,5e-4))
        
    if ecu > 0.006:
        df = np.append(df[df<=50e-4],np.arange(60e-4,20*ecu,10e-4))
        
    if ecu > 0.012: 
        df = np.append(df[df<=100e-4],np.arange(120e-4,20*ecu,20e-4))
    
    npts = len(df)
    
    if PP[i] > 0:
        for k in range(npts):
            f1 = interp1d(ecun,fcun); temp1 = np.sum(f1(df[0]*np.ones(len(yl)))*conclay[:,1])
            f2 = interp1d(ec,fc); temp2 = np.sum(f2(df[0]*np.ones(len(yl)))*conclay[:,2])
            f3 = interp1d(es,fs); temp3 = np.sum(Asb*f3(df[0]*np.ones(len(distld))))
            compch = temp1 + temp2 + temp3
            if compch < PP[i]:
                df = df[1:]
    
    npts = len(df)

    msg = 0                                 # stop conditions
    
    curv   = [0]                            # curvatures
    mom    = [0]                            # moments
    ejen   = [0]                            # neutral axis
    DF     = [0]                            # force eqilibrium
    vniter = [0]                            # iterations
    coverstrain = [0]                          
    corestrain  = [0]                          
    steelstrain = [0]
    
    x = [D/2]
    for k in range(npts):
        lostmomcontrol = max(mom)
        if mom[k]<(0.8*lostmomcontrol):
            msg = 4
            break
        
        F = 10*tol
        niter = -1                                  
        while abs(F)>tol:
            niter = niter + 1
            eec = (df[k]/x[niter])*(conclay[:,0]-(D-x[niter]))    # vector with the strains in the concrete
            ees = (df[k]/x[niter])*(distld-(D-x[niter]))          # vector with the strains in the steel
    
            fcunconf = interp1d(ecun,fcun,fill_value="extrapolate")(eec)        # vector with stresses in the unconfined concr.           
            fcconf = interp1d(ec,fc,fill_value="extrapolate")(eec)              # vector with stresses in the confinded concr.    
            fsteel = interp1d(es,fs,fill_value="extrapolate")(ees)              # vector with stresses in the steel
            FUNCON = fcunconf*conclay[:,1]
            FCONF  = fcconf*conclay[:,2]
            FST    = Asb*fsteel
            F      = np.sum(FUNCON) + np.sum(FCONF) + np.sum(FST) - PP[i]
            
            if F>0:
                x.append(x[niter] - 0.05*x[niter])
            
            elif F<0:
                x.append(x[niter] + 0.05*x[niter])
    
            if niter>itermax:
                msg = 3
                break
    
        cores = (df[k]/x[niter])*abs(x[niter]-dcore)
        TF = confined == unconfined
        if not TF:
            if cores >= ecu:
                msg = 1
                break
    
        elif TF:
            if df[k] >= ecu:
                msg = 1
                break
            
        if abs(ees[0]) > esu:
            message = 2
            break
    
        ejen.append(x[niter])
        DF.append(x[niter])
        vniter.append(niter)
        temp = (sum(FUNCON*conclay[:,0]) + sum(FCONF*conclay[:,0]) + sum(FST*distld) - PP[i]*(D/2))/(10**6)
        mom.append(temp)
    
        if mom[k+1] < 0:
            mom[k+1] = -0.01*mom[k+1]
    
        curv.append(1000*df[k]/x[niter])
        coverstrain.append(df[k])
        corestrain.append(cores)
        steelstrain.append(ees[0])
        x[0] = x[niter]
        del x[1:]
        if msg!=0:
            break
    
    Mni.append(interp1d(coverstrain,mom)(csid))
    esaux    = interp1d(coverstrain,steelstrain)(csid)
    cr = 0 # concrete control
    if abs(esaux) > abs(ssid) or np.isnan(Mni[i]):
        cr = 1 # steel concrete
        Mni[i]    = interp1d(steelstrain,mom)(-ssid)
    
    if cr == 0:
        eci.append(csid)
        esi.append(esaux)
    if cr == 1:
        esi.append(-ssid)
        eci.append(float(interp1d(steelstrain,coverstrain)(-ssid)))    
        
    Msgs.append(msg)

Mni = np.asarray([0]+Mni+[0])
PPn  = np.asarray([-PTid] + PP + [PCid])

MB = np.max(Mni)
PB = PPn[np.where(Mni==MB)[0]][0]

PB13 = (1/3)*PB
MB13 = float(interp1d(PPn,Mni)(PB13))

PB23 = (2/3)*PB
MB23 = float(interp1d(PPn,Mni)(PB23))

MB0 = float(interp1d(PPn,Mni)(0))

PPL = np.asarray([-PTid, 0, PB13, PB23, PB, PCid])
MnL = np.asarray([0, MB0, MB13, MB23, MB, 0])

plt.figure(figsize=(8,6))
plt.plot(Mni,PPn/1000, c = 'red', linestyle = '-', marker = 'o', MarkerEdgeColor = 'black', MarkerFaceColor = 'red', markersize=6, label = 'Interaction Diagram')
plt.plot(MnL,PPL/1000, c = 'blue', linestyle = '--', marker = 's', MarkerEdgeColor = 'black', MarkerFaceColor = 'blue', markersize=6, label = 'Approximation for NLTHA')
plt.xlabel('Moment [kN.m]', fontsize = 16)
plt.ylabel('Axial Load [kN]', fontsize = 16)
plt.title('Interaction Diagram', fontsize=16)
plt.grid(True)
        
fid.write('\n\n*** Interaction Surface ***\n')
fid.write('Concrete limit strain:\t%1.4f\n' % csid)
fid.write('Steel limit strain:\t%1.4f\n' % ssid)
fid.write('\nMoment\tAxial Load\n') 
fid.write('[kN-m]\t[kN]\n')
for i in range(len(Mni)):
    fid.write('%8.2f\t%8.2f\n' % (Mni[i], PPn[i]))
    
fid.write('\n')
fid.write('NLTHA Approximation:\n\n')
fid.write('PT:\t%6.1f\tkN\n' % float(-PTid/1000))
fid.write('PC:\t%6.1f\tkN\n' % float(PCid/1000))
fid.write('PB:\t%6.1f\tkN\tMB:\t%6.1f\tkN-m\n' % (float(PB/1000),MB))
fid.write('(1/3)PB:\t%6.1f\tkN\t(1/3)MB:\t%6.1f\tkN.m\n' % (float(PB13/1000),MB13))
fid.write('(2/3)PB:\t%6.1f\tkN\t(2/3)MB:\t%6.1f\tkN.m\n' % (float(PB23/1000),MB23))

plt.figure(figsize=(8,6))
plt.plot(PPn[1:-1]/1000, eci, c = 'blue', linestyle = '--', marker = 's', MarkerEdgeColor = 'black', MarkerFaceColor = 'blue', markersize=6)
plt.xlabel('Axial Force [kN]', fontsize = 16)
plt.ylabel('Concrete Strain', fontsize = 16)
plt.grid(True)

plt.figure(figsize=(8,6))
plt.plot(PPn[1:-1]/1000, esi, c = 'blue', linestyle = '--', marker = 's', MarkerEdgeColor = 'black', MarkerFaceColor = 'blue', markersize=6)
plt.xlabel('Axial Force [kN]', fontsize = 16)
plt.ylabel('Steel Strain', fontsize = 16)
plt.grid(True)

fid.close()