import matplotlib.patches as patches
import numpy as np
import awkward as ak
def f_1m(x):
    """return cross section [pb] of e/mu HNL of mass (x)[GeV] at ct=1m"""
    x0 = np.array([1,2,4,5,7,10])
    #y0 = np.array([13.57,0.4238,0.01289,0.004156,0.0007452,0.000121])    ## old xsec 
    y0 = np.array([8.492,0.2653,0.00809,2.612E-03,4.721E-04,7.751E-05])    
    return np.exp(np.poly1d(np.polyfit(x0,np.log(y0),5))(x))

def f_xsec(m):
    """ return a function of f(ctau[mm]) = pb for mass m[GeV]"""
    def xsec_m(x):
        return f_1m(m)/(x/1000.)
    return xsec_m
def f_v2(m):
    """return a function of f(ctau[mm]) = V2  for mass m[GeV]"""
    def ctau_m(x):
        return f_xsec(m)(x)*4.8759E-05 ## const. between xsec and v^2
    return ctau_m

def f_1m_tau(x):
    """ function for tau"""
    x0 = np.array([1,2,4,7,10])
    y0 = np.array([30700,963,17.9,0.646,0.0912])/1000.
    return np.exp(np.poly1d(np.polyfit(x0,np.log(y0),4))(x))

def f_xsec_tau(m):
    """ return a function of f(ctau[mm]) = pb for mass m[GeV]"""
    def xsec_m_tau(x):
        return f_1m_tau(m)/(x/1000.)
    return xsec_m_tau
def f_v2_tau(m):
    """ return a function of f(ctau[mm]) = V2  for mass m[GeV]"""
    def ctau_m_tau(x):
        return f_xsec_tau(m)(x)*4.8759E-05 ## const. between xsec and v^2
    return ctau_m_tau

from coffea.nanoevents.methods import vector
def pack(events,obj_str):
    """ Pack branches that starts with obj_str from events into an ak Array obj.
        e.g. muons  = util.pack(events,"Muon_")
    """
    obj  = ak.zip(
                {k.replace(obj_str,""):getattr(events,k) for k in events.fields if k.startswith(obj_str)}
                ,with_name="PtEtaPhiMLorentzVector",
                behavior=vector.behavior
               )
    return obj

def delta_r_pairs(pair_dict):
    """ Make delta_r pairs from a dict:
        e.g. pair_dict = {"cls":cluster,"bjets":bjets}
    """
    pairs = ak.cartesian(pair_dict,axis=1,nested=True)
    main_obj = getattr(pairs,pairs.fields[0])
    sec_obj = getattr(pairs,pairs.fields[1])
    dr_pairs = main_obj.delta_r(sec_obj)
    return dr_pairs,pairs

def delta_phi_pairs(pair_dict):
    """ Make delta_phi pairs from a dict:
        e.g. pair_dict = {"cls":cluster,"bjets":bjets}
    """
    pairs = ak.cartesian(pair_dict,axis=1,nested=True)
    main_obj = getattr(pairs,pairs.fields[0])
    sec_obj = getattr(pairs,pairs.fields[1])
    dphi_pairs = main_obj.delta_phi(sec_obj)
    return dphi_pairs,pairs

def dr_min(obj1,obj2,with_name):
    """ Compute the min. dr of obj2 w.r.t. obj1, return as a field of obj1    """

    dr,pairs = delta_r_pairs({"obj1":obj1,"obj2":obj2})
    #min dR pair
    dr_min = ak.min(dr,axis=2)
    return ak.with_field(obj1,dr_min,with_name)

def match(obj1,obj2,with_name,dR_cut=0.4):
    """
     Bool of the closest obj2 < dR_cut w.r.t. obj1, return as a field of obj1
    """
    dr,pairs = delta_r_pairs({"obj1":obj1,"obj2":obj2})
    #min dR pair
    dr_min = ak.min(dr,axis=2)
    #boolean of matched or not
    matched = dr_min<dR_cut
    
    return ak.with_field(obj1,matched,with_name)

def dphi_min(obj1,obj2,with_name):
    """
     Compute the min. dphi of obj2 w.r.t. obj1, return as a field of obj1
    """
    dphi,pairs = delta_phi_pairs({"obj1":obj1,"obj2":obj2})
    #min dphi pair
    dphi_min = ak.min(dphi,axis=2)
    return ak.with_field(obj1,dphi_min,with_name)

#Drawing utils
def drawCSCsteel(ax,hORv='v'): 
    y_max = ax.get_ylim()[1]
    if hORv=='v':
        ax.axvline(632);ax.axvline(632+39)
        ax.axvline(724);ax.axvline(724+64)
        ax.axvline(849);ax.axvline(849+62)
        ax.axvline(970);ax.axvline(970+32)
    else:
        ax.axhline(632);ax.axhline(632+39)
        ax.axhline(724);ax.axhline(724+64)
        ax.axhline(849);ax.axhline(849+62)
        ax.axhline(970);ax.axhline(970+32)   
        
    ax.text(570, y_max*1.02, 'ME1/1', fontsize=12)
    ax.text(670, y_max*1.02, 'ME1/2-3', fontsize=12)
    ax.text(800, y_max*1.02, 'ME2', fontsize=12)
    ax.text(920, y_max*1.02, 'ME3', fontsize=12)
    ax.text(1015, y_max*1.02,'ME4', fontsize=12)        
    return 

def drawCSCr(ax):
    y_max = ax.get_ylim()[1]
    ax.axvline(350,linestyle="--",color='grey')
    ax.text(350-110,y_max*0.05, "Inner ring", fontsize=15)
    ax.text(350+15 ,y_max*0.05, "Outer ring", fontsize=15)
    return ax

def drawCSCz(ax,text_loc=0.7):    
    ax.set_xlim(550,1075)
    (xmin,xmax) = ax.get_xlim()

    y_max = ax.get_ylim()[1]

    preME11 = patches.Rectangle((xmin, 0), 568-xmin, 2,color='grey',alpha=0.3)
    ME11_12 = patches.Rectangle((632, 0), 39, 2,color='grey',alpha=0.3)
    ME12_2  = patches.Rectangle((724, 0), 65, 2,color='grey',alpha=0.3)
    ME2_3   = patches.Rectangle((849, 0), 62, 2,color='grey',alpha=0.3)
    ME3_4   = patches.Rectangle((970, 0), 32, 2,color='grey',alpha=0.3)
    beyond  = patches.Rectangle((1050, 0),50, 2,color='grey',alpha=0.3)

    ax.text(570, y_max*1.02, 'ME1/1', fontsize=12)
    ax.text(670, y_max*1.02, 'ME1/2-3', fontsize=12)
    ax.text(800, y_max*1.02, 'ME2', fontsize=12)
    ax.text(920, y_max*1.02, 'ME3', fontsize=12)
    ax.text(1015, y_max*1.02,'ME4', fontsize=12)
    ax.text(xmin+5 ,y_max*0.15, "Steel", fontsize=15,rotation=90)
    ax.text(xmax-20,y_max*0.15, "Beyond CMS", fontsize=15,rotation=90)

    ax.add_patch(preME11)
    ax.add_patch(ME11_12)
    ax.add_patch(ME12_2)
    ax.add_patch(ME2_3)
    ax.add_patch(ME3_4)
    ax.add_patch(beyond)
    return ax


def drawDTz(ax):
    y_max = ax.get_ylim()[1]
    ax.axvline(140,linestyle="--",color='grey')
    ax.axvline(400,linestyle="--",color='grey')    

    ax.text(140-110,y_max*0.05, "Wheel 0", fontsize=15)
    ax.text(220 ,y_max*0.05, "Wheel 1", fontsize=15)
    ax.text(500 ,y_max*0.05, "Wheel 2", fontsize=15)    
    return ax

def drawDTr(ax,text_loc=0.7):    
    ax.set_xlim(350,760)
    (xmin,xmax) = ax.get_xlim()

    y_max = ax.get_ylim()[1]

    preMB1 = patches.Rectangle((xmin, 0), 405-xmin, 2,color='grey',alpha=0.3)
    MB1_2  = patches.Rectangle((465, 0) , 20, 2,color='grey',alpha=0.3)
    MB2_3  = patches.Rectangle((540, 0) , 55, 2,color='grey',alpha=0.3)
    MB3_4  = patches.Rectangle((640, 0) , 60, 2,color='grey',alpha=0.3)
    

    ax.text(420, y_max*1.02, 'MB1', fontsize=12)
    ax.text(500, y_max*1.02, 'MB2', fontsize=12)
    ax.text(600, y_max*1.02, 'MB3', fontsize=12)
    ax.text(710, y_max*1.02, 'MB4', fontsize=12)
    ax.axvline(740,linestyle="--",color='grey')    
    
    ax.text(xmin+5 ,y_max*0.15, "Steel/Solenoid", fontsize=15,rotation=90)
    ax.text(xmax-10,y_max*0.15, "Beyond CMS", fontsize=15,rotation=90)
    

    ax.add_patch(preMB1)
    ax.add_patch(MB1_2)
    ax.add_patch(MB2_3)
    ax.add_patch(MB3_4)
    return ax

def drawRZ(ax):
    ax.set_ylim(0,750)
    ax.set_xlim(0,1100)
    (xmin,xmax) = ax.get_xlim()

    y_max = ax.get_ylim()[1]
    
    MB1 = patches.Rectangle((MB_lim,402),661-MB_lim,449-402,color='grey',alpha=0.3)
    MB2 = patches.Rectangle((MB_lim,490),661-MB_lim,533-490,color='grey',alpha=0.3)
    MB3 = patches.Rectangle((MB_lim,597),661-MB_lim,636-597,color='grey',alpha=0.3)
    MB4 = patches.Rectangle((MB_lim,700),661-MB_lim,738-700,color='grey',alpha=0.3)

    solenoid= patches.Rectangle((MB_lim,295),661-MB_lim,380-295,color='grey',alpha=0.3)

    ME22 = patches.Rectangle((791,357),850-791,700-357,color='grey',alpha=0.3)
    ME32 = patches.Rectangle((911,357),970-911,700-357,color='grey',alpha=0.3)
    ME42 = patches.Rectangle((1002,357),1063-1002,700-357,color='grey',alpha=0.3)

    ME21 = patches.Rectangle((789,139),850-789,345-139,color='grey',alpha=0.3)
    ME31 = patches.Rectangle((915,160),970-915,345-160,color='grey',alpha=0.3)
    ME41 = patches.Rectangle((1002,178),1063-1002,345-178,color='grey',alpha=0.3)

    ME11 = patches.Rectangle((580,100),632-580,275-100,color='grey',alpha=0.3)
    ME12 = patches.Rectangle((668,275),724-668,465-275,color='grey',alpha=0.3)
    ME13 = patches.Rectangle((686,505),724-686,700-505,color='grey',alpha=0.3)
    
    list_of_boxes = [ME11,ME12,ME13,ME21,ME22,ME31,ME32,ME41,ME42,MB1,MB2,MB3,MB4,solenoid]
    
    for box in list_of_boxes:
        ax.add_patch(box)
    
    ax.text(410, 415, "MB1",fontsize=12)
    ax.text(410, 500, "MB2",fontsize=12)
    ax.text(410, 605, "MB3",fontsize=12)
    ax.text(410, 708, "MB4",fontsize=12)
    
    ax.text(665, 705, "ME1/3",fontsize=12)
    ax.text(780, 705, "ME2/2",fontsize=12)
    ax.text(900, 705, "ME3/2",fontsize=12)
    ax.text(990, 705, "ME4/2",fontsize=12)


    #ax.text(680, 140, "Steel",fontsize=12)
    #ax.text(410,660, "Steel",fontsize=12)
    ax.text(410, 330, "Solenoid",fontsize=12)
    #ax.text(430, 150, "HCAL",fontsize=12)

    ax.text(615, 110, "ME1/1",rotation='vertical',fontsize=12)
    ax.text(695, 331, "ME1/2",rotation='vertical',fontsize=12)
    ax.text(830, 145, "ME2/1",rotation='vertical',fontsize=12)
    ax.text(955, 170, "ME3/1",rotation='vertical',fontsize=12)
    ax.text(1040, 190, "ME4/1",rotation='vertical',fontsize=12)
    return ax


def plotEff(h_list,ax,axis='e'):
    for h in h_list:
        h_pass=h.integrate('selection',slice(1,None))
        h_pass.label='Cluster Efficiency'
        
        hist.plotratio(
            ax=ax,
            num  =h_pass.project(axis),
            denom=h.project(axis),
            error_opts={'marker': '.'},
            unc='num',
            label=h.identifiers('dataset')[0].label,
            clear=False
        )

    ax.legend(loc='best')
    return ax

