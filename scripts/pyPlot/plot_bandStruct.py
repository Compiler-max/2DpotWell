import matplotlib.pyplot as plt
import numpy as np
import re
from util_2dPW_Interf import read_AbIn_energies
from util_2dPW_Interf import read_w90_energies
from util_2dPW_Interf import getData


#adjust colors etc.
marker_abi  = '+-'
markerS_abi = 1.2
line_abi    = 0.6
color_abi   = 'black'

marker_w90  = '*-'
markerS_w90 = 1.0
line_w90    = 0.5
color_w90   = 'red'

show_Bfield_box = True
show_rashba_box = True


def plotBands(w90_dir=".",out_dir=".", minEn=0, maxEn=0, show_Bfield_box=False, show_rashba_box=False):
    #GET ABINITIO ENERGIES
    do_ABiN, nQ, nSolve, qpts, en_abi = read_AbIn_energies(out_dir+'/enABiN.txt')
  
    #GET INTERPOLATED ENERGIES
    do_w90, nK, nWfs, kpts, en_W90, velo_w90 = read_w90_energies(w90_dir+'/wf1_geninterp.dat')

    #PLOT
    fig, ax = plt.subplots(1,1)

    #PLOT ABINITIO
    if do_ABiN:
        qpath, qplot, qticks = get_BZ_path(qpts)
        #
        #
        for n in range(0,nSolve):
            for path in range(0,len(qpath)):
                enPlot = []
                for q in range(0,len(qpath[path])):
                    q_idx = qpath[path][q][0]
                    enPlot.append(en_abi[q_idx,n])
                #
                xPlot = qplot[path][0:len(qplot[0])]
                yPlot = enPlot[0:len(enPlot)]
                ax.plot(xPlot, yPlot,  marker_abi, markersize=markerS_abi,linewidth=line_abi, color=color_abi)
        print('...finished abinitio plot')


    #PLOT W90 INTERPOLATION
    if do_w90:
        kpath, kplot, kticks = get_BZ_path(kpts)   
        #
        for n in range(0,nWfs):
            for path in range(0,len(kpath)):
                enPlot = []
                for k in range(0,len(kpath[path])):
                    k_idx = kpath[path][k][0]
                    enPlot.append(en_W90[k_idx,n])
                #
                xPlot = kplot[path][0:len(kplot[0])]
                yPlot = enPlot[0:len(enPlot)]
                ax.plot(xPlot, yPlot,  marker_w90, markersize=markerS_w90,linewidth=line_w90, color=color_w90)
        print('...finished w90 interpolation plot')
       


    #X TICK-LABELS
    xtickLabel = np.array(['M','X',r'$\Gamma$','M','Y',r'$\Gamma$'])
    if do_ABiN:
        if len(qticks) != len(xtickLabel)   :
            print('ERROR: xtickLabel has wrong size')
        ax.set_xticks(qticks)
        ax.set_xticklabels(xtickLabel,fontsize=14)
        ax.set_xlim(qticks[0],qticks[-1])
    #if no abinit found try to use the w90 labels
    elif do_w90:
        if len(kticks) != len(xtickLabel)   :
            print('ERROR: xtickLabel has wrong size')
        ax.set_xticks(kticks)
        ax.set_xticklabels(xtickLabel,fontsize=14)
        ax.set_xlim(kticks[0],kticks[-1])

    #Y-LABELS 
    if maxEn-minEn > 0: 
        ax.set_ylim(minEn,maxEn)
    #
    plt.ylabel('E [eV]',fontsize=16)

    #GET PARAMETERS
    alpha_rashba    = getData('alpha_rashba'    ,out_dir+'/polOutput.txt')
    B_ext           = getData('magnetic_field'  ,out_dir+'/polOutput.txt')


    #PLOT PARAMETER BOX
    if show_Bfield_box:
        ax.text(0.8, 0.9, 'B_z='+str(B_ext[2])+' T', 
                horizontalalignment='center', 
                verticalalignment='center',
                transform = ax.transAxes, 
                fontsize=16, bbox={'facecolor':'white', 'alpha':0.8, 'pad':10})
    if show_rashba_box:
        ax.text(0.7, 0.9, r'$\alpha_{Rashba}=$'+str(alpha_rashba)+r'eV $\AA$', 
                horizontalalignment='center', 
                verticalalignment='center',
                transform = ax.transAxes, 
                fontsize=16, bbox={'facecolor':'white', 'alpha':0.8, 'pad':10})
    #
    #SAVE FILE
    file_name = 'en_B'+str(B_ext[2])+'_aRashb'+str(alpha_rashba)+'bands.pdf'
    file_path = out_dir+'/'+file_name
    plt.savefig(file_path,bbox_inches='tight')
    print('saved plot to file',file_path)
    #    
    #
    return ax

  




def get_BZ_path(qpts):
    #get q point path
    tol         = 1e-8
    qxmin       = min(qpts[:,0])
    qymin       = min(qpts[:,1]) 

    M_vec       = np.array([qxmin,qymin,0.0])
    M_norm      = np.linalg.norm(M_vec)
    eM_vec      = M_vec / M_norm

    qpathI   = []#np.empty([0,0,0],dtype=int)
    qpathII  = []
    qpathIII = []
    qpathIV  = []
    qpathV   = []


    for ind,q_vec in enumerate(qpts):
        if q_vec[0]<= 0.0 and q_vec[1]<=0.0:
            q_norm  = np.linalg.norm(q_vec)
            eQ_vec  = np.array([0.0,0.0,0.0])
            if( abs(q_norm) > 1e-8):
                eQ_vec  = np.divide(q_vec, q_norm)

            #I  :    M->X
            if abs(q_vec[0]-qxmin) < tol:
                qpathI.append( ([ind,q_vec[0],q_vec[1]])        )
                
            #II :    X->Gamma
            if abs(q_vec[1]) < tol:
                qpathII.append( ([ind,q_vec[0],q_vec[1]])        )


            #III:   Gamma->M
            if np.linalg.norm( eQ_vec - eM_vec) < tol:
                if q_norm <= M_norm:
                    qpathIII.append( ([ind,q_vec[0],q_vec[1]])        )
                else:
                    print('WARNING: found q-vector of of bounds:',q_vec)

            #IV :    M->Y
            if abs(q_vec[1]-qymin) < tol:
                qpathIV.append( ([ind,q_vec[0],q_vec[1]])      )

            #V  :   Y->Gamma
            if abs(q_vec[0]) < tol:
                qpathV.append( ([ind,q_vec[0],q_vec[1]])        )


    #add gamma point to start of pathIII
    qpathTmp = ([qpathII[-1]])
    for qpt in qpathIII:
        qpathTmp.append((qpt))
    qpathIII = qpathTmp

    qpathI  = sorted(qpathI,    key= lambda elem: elem[2]               )
    qpathII = sorted(qpathII,   key= lambda elem: elem[1]               )
    qpathIII= sorted(qpathIII,  key= lambda elem: elem[1], reverse=True )
    qpathIV = sorted(qpathIV,   key= lambda elem: elem[1]               )
    qpathV  = sorted(qpathV,    key= lambda elem: elem[2]               )


    print('start point:',qpathI[0][:])

       #check if q indices of start and end points match
    if qpathI[-1][0] is not qpathII[0][0]:
        print('WARNING: Q-path not smoth at I->II')
    if qpathII[-1][0] is not  qpathIII[0][0]:
        print('WARNING: Q-path not smoth at II->III')
    if qpathIII[-1][0] is not qpathIV[0][0]:
        print('WARNING: Q-path not smoth at III->IV')
    if qpathIV[-1][0] is not  qpathV[0][0]:
        print('WARNING: Q-path not smoth at IV->V')
  


    qpath   = []
    qpath.append(   qpathI     )
    qpath.append(   qpathII    )
    qpath.append(   qpathIII   )
    qpath.append(   qpathIV    )
    qpath.append(   qpathV     )



    qplot   = []
    offset  = 0.0
    qpVect  = []
    xticks  = [0.0]
    for path in range(0,len(qpath)):
        #get length
        start  = np.array([  qpath[path][ 0][1], qpath[path][ 0][2]         ])
        end    = np.array([  qpath[path][-1][1], qpath[path][-1][2]         ])
        qpVect.append(end-start)
        length  =np.linalg.norm(qpVect[path])

        #make linspace
        tmp = np.linspace(offset,offset+length,len(qpath[path]))  
        qplot.append( tmp ) 

        offset = offset + length
        xticks.append(offset)
    qplot = np.array(qplot)

    return qpath, qplot, xticks








#TEST
#axTest = plotBands()
#plt.show()

