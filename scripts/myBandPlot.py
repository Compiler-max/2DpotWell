import matplotlib.pyplot as plt
import numpy as np
import re


#dir to search
dim			= 2 
#w90_dir     = 'w90files'
#out_dir     = 'output'

def plotBands(w90_dir,out_dir, minEn, maxEn):
    print('Hello from the plot function')
    #GET ABINITIO ENERGIES
    print('*')
    print('*')
    print('*')
    print('********read abinitio energies:********************************')

    #Get external field
    with open(out_dir+'/'+'enABiN.txt','r') as f:
        parser = False
        for count,line in enumerate(f.readlines()):
            if 'B_ext' in line:
                #strip text around values
                sub     = line.partition('=')
                sub     = sub[2].partition('T')[0]
                   
                #find all floats in the string "sub"
                B_ext   = re.findall(r"[-+]?\d*\.\d+|\d+",sub)
                #convert to float point
                B_ext   = np.array(B_ext).astype(np.float)   

    if( len(B_ext) != 3):
        print('WARNING: problems reading the external magnetic field (should be 3d vector)')
    print('detected B_ext=',B_ext)




    outFile     = open(out_dir+'/'+'enABiN.txt','r')
    data        = np.genfromtxt(outFile, dtype=[('int'),('float64'),('float64'),('float64'),('float64')])
    #
    qpts = []
    en_abi= []
    for counter, line in enumerate(data):
        en_abi.append( line[4] )
        if line[0] > len(qpts):
            qpts.append(    ([line[1],line[2],line[3]])   )   
    #
    qpts    = np.array(qpts)
    nQ      = len(qpts)
    nSolve  = int(  len(en_abi)/ nQ )
    en_abi= np.reshape(en_abi,(nQ,nSolve))
    #
    print('detected nQ  =',nQ   )
    print('detected nSolve=',nSolve )



    #GET INTERPOLATED ENERGIES
    print('*')
    print('*')
    print('*')
    print('********read w90 energies:********************************')
    w90interp   = open(w90_dir+'/'+'wf1_geninterp.dat','r')
    data        = np.genfromtxt(w90interp, dtype=[('int'),('float64'),('float64'),('float64'),('float64'),('float64'),('float64'),('float64')])
    #
    kpts = []
    en_W90= []
    for counter, line in enumerate(data):
        en_W90.append( line[4] )
        if line[0] > len(kpts):
            kpts.append(    ([line[1],line[2],line[3]])   )   
    #
    kpts    = np.array(kpts)
    nK      = len(kpts)
    nWfs    = int(  len(en_W90)/ nK )
    en_W90   = np.reshape(en_W90,(nK,nWfs))
    #
    print('detected nK  =',nK   )
    print('detected nWfs=',nWfs )




    #SEARCH K-SPACE PATHS
    print('*')
    print('*')
    print('*')
    print('********get k-space pathes********************************')
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
            eQ_vec  = q_vec / q_norm

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

    #print('sorted q-paths (abi):')
    #print('I  :',qpath[0]   )
    #print('II :',qpath[1]   )
    #print('III:',qpath[2]   )
    #print('IV :',qpath[3]   )
    #print('V  :',qpath[4]   )



    #get k point path
    tol     = 1e-8
    kxmin   = min(kpts[:,0])
    kymin   = min(kpts[:,1]) 

    M_vec       = np.array([kxmin,kymin,0.0])
    M_norm      = np.linalg.norm(M_vec)
    eM_vec      = M_vec / M_norm


    kpathI   = []#np.empty([0,0,0],dtype=int)
    kpathII  = []
    kpathIII = []
    kpathIV  = []
    kpathV   = []


    
    for ind,k_vec in enumerate(kpts):
        if k_vec[0]<= 0.0 and k_vec[1]<=0.0:
            k_norm  = np.linalg.norm(k_vec)
            eK_vec  = k_vec / k_norm

            #I  :    M->X
            if abs(k_vec[0]-kxmin) < tol:
                kpathI.append( ([ind,k_vec[0],k_vec[1]])        )
                
            #II :    X->Gamma
            if abs(k_vec[1]) < tol:
                kpathII.append( ([ind,k_vec[0],k_vec[1]])        )


            #III:   Gamma->M
            if np.linalg.norm( eK_vec - eM_vec) < tol:
                if k_norm <= M_norm:
                    kpathIII.append( ([ind,k_vec[0],k_vec[1]])        )
                else:
                    print('WARNING: found k-vector of of bounds:',k_vec)

            #IV :    M->Y
            if abs(k_vec[1]-kymin) < tol:
                kpathIV.append( ([ind,k_vec[0],k_vec[1]])      )

            #V  :   Y->Gamma
            if abs(k_vec[0]) < tol:
                kpathV.append( ([ind,k_vec[0],k_vec[1]])        )


    #add gamma point to start of pathIII
    kpathTmp = ([kpathII[-1]])
    for kpt in kpathIII:
        kpathTmp.append((kpt))
    kpathIII = kpathTmp

    kpathI  = sorted(kpathI,    key= lambda elem: elem[2]               )
    kpathII = sorted(kpathII,   key= lambda elem: elem[1]               )
    kpathIII= sorted(kpathIII,  key= lambda elem: elem[1], reverse=True )
    kpathIV = sorted(kpathIV,   key= lambda elem: elem[1]               )
    kpathV  = sorted(kpathV,    key= lambda elem: elem[2]               )

    #check smoothness of path
    if kpathI[-1][0] is not kpathII[0][0]:
        print('WARNING: K-path not smoth at I->II')
    if kpathII[-1][0] is not kpathIII[0][0]:
        print('WARNING: K-path not smoth at II->III')
    if kpathIII[-1][0] is not kpathIV[0][0]:
        print('WARNING: K-path not smoth at III->IV')
    if kpathIV[-1][0] is not kpathV[0][0]:
        print('WARNING: K-path not smoth at IV->V')
  

    kpath   = []
    kpath.append(kpathI)
    kpath.append(kpathII)
    kpath.append(kpathIII)
    kpath.append(kpathIV)
    kpath.append(kpathV)


    #print('sorted k-paths (w90):')
    #print('I  :',kpath[0])
    #print('II :',kpath[1])
    #print('III:',kpath[2])
    #print('IV :',kpath[3])
    #print('V  :',kpath[4])




    #PLOT:
    print('*')
    print('*')
    print('**********prepare Plots now:******************')


    if( abs(qxmin-kxmin) > 1e-8 ):
        print('ERROR: wannier kmesh seems to live on different unit cell (x-coord) ')
    if( abs(qymin-kymin) > 1e-8 ):
        print('ERROR: wannier kmesh seems to live on different unit cell (y-coord) ')


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


    kplot   = []
    offset  = 0.0
    kpVect  = []
    for path in range(0,len(kpath)):
        start  = np.array([  kpath[path][ 0][1], kpath[path][ 0][2]         ])
        end    = np.array([  kpath[path][-1][1], kpath[path][-1][2]         ])
        kpVect.append(end-start)
        length  =np.linalg.norm(kpVect[path])

        tmp = np.linspace(offset,offset+length,len(kpath[path]))  
        kplot.append( tmp ) 

        offset = offset + length
    kplot = np.array(kplot)




    #TEST IF WANNIER AND ABINITIO GIVE SAME PATHS
    tol = 1e-8
    for path in range(0,len(qpVect)):
        if np.linalg.norm(qpVect[path]-kpVect[path]) > tol:
            print('Error: path #',path,' differs from abinit to w90')
            print('qVect=',qpVect[path])
            print('kVect=',kpVect[path])
        else:
            print('paths #',path,' match')











    fig, ax = plt.subplots(1,1)

    #PLOT ABINITIO
    marker_abi  = '+-'
    markerS_abi = 1.2
    line_abi    = 0.6
    color_abi   = 'black'
    for n in range(0,nSolve):

        for path in range(0,len(qpath)):
            enPlot = []
            for q in range(0,len(qpath[path])):
                q_idx = qpath[path][q][0]
                enPlot.append(en_abi[q_idx,n])
            
            xPlot = qplot[path][0:len(qplot[0])]
            yPlot = enPlot[0:len(enPlot)]
            ax.plot(xPlot, yPlot,  marker_abi, markersize=markerS_abi,linewidth=line_abi, color=color_abi)





    #PLOT W90 INTERPOLATION
    marker_w90  = '*-'
    markerS_w90 = 1.0
    line_w90    = 0.5
    color_w90   = 'red'

    for n in range(0,nWfs):

        for path in range(0,len(kpath)):
            enPlot = []
            for k in range(0,len(kpath[path])):
                k_idx = kpath[path][k][0]
                enPlot.append(en_W90[k_idx,n])
            
            xPlot = kplot[path][0:len(kplot[0])]
            yPlot = enPlot[0:len(enPlot)]
            ax.plot(xPlot, yPlot,  marker_w90, markersize=markerS_w90,linewidth=line_w90, color=color_w90)
       


    #LABELS
    xtickLabel = np.array(['M','X',r'$\Gamma$','M','Y',r'$\Gamma$'])
    if len(xticks) != len(xtickLabel)   :
        print('ERROR: xtickLabel has wrong size')


    ax.set_xticks(xticks)

    ax.set_xticklabels(xtickLabel,fontsize=14)
    ax.set_xlim(xticks[0],xticks[-1])
    ax.set_ylim(minEn,maxEn)
    #
    plt.ylabel('E [eV]',fontsize=16)



    ax.text(0.8, 0.9, 'B_z='+str(B_ext[2])+' T', horizontalalignment='center', verticalalignment='center',transform = ax.transAxes, fontsize=16, 
            bbox={'facecolor':'white', 'alpha':0.8, 'pad':10})


    return ax

    plt.savefig('enB'+str(B_ext[2])+'bands.pdf')
    #plt.show()

