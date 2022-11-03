import numpy as np
import ctypes
#import ReadBinary as RB
from scipy.optimize import fmin
import time 
import multiprocessing as mp
import emcee
#from mpi4py import MPI
from PyAstronomy import pyasl
Chi2Lib = ctypes.cdll.LoadLibrary('./chi2_calculator.so')

map_set_name = 'ob161195_test'

def chi2(tempparms,data,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,box_size,use_mcmc):
    initial_type = ctypes.c_float*4
    initial = initial_type()
    initial[0] = tempparms[0]
    initial[1] = tempparms[1]
    initial[2] = tempparms[2]
    initial[3] = tempparms[3]

    chi2 = 0.
    ichi2type = ctypes.c_float
    ichi2 = ichi2type()
    ip = ctypes.byref(ichi2)

    box_size_type = ctypes.c_float
    box_size_final = box_size_type()
    box_size_final = box_size

    Chi2Lib.wrapgetchi2.argtypes = [ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float), ctypes.c_int, ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.c_float,ctypes.POINTER(ctypes.c_float)]
    Chi2Lib.wrapprintlc.argtypes = [ctypes.POINTER(ctypes.c_float),ctypes.c_int,ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.c_float]

    for idata in data:
        nhjd = len(idata)
        temptype = ctypes.c_float*nhjd
        hjds = temptype()
        iflux = temptype()
        iferr = temptype()
        for i in range(nhjd):
            hjds[i]=idata[i,0]
            iflux[i]=idata[i,1]
            iferr[i]=idata[i,2]
        Chi2Lib.wrapgetchi2(hjds,iflux,iferr,nhjd,initial,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,box_size_final,ip)
        #if len(idata)==887:
        #    Chi2Lib.wrapprintlc(hjds,nhjd,initial,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,box_size_final)

        if ichi2.value <0 :
            chi2 = np.inf
            break
        chi2 += ichi2.value
    #print chi2
    if use_mcmc == True:
        return -0.5*chi2
    else:
        return chi2


def grid(n):
    print (n)
    global use_mcmc
    #mapdir = "./map_set_kb220371_version2/"#test/"
    
    map_content = np.load(mapdir+"%s.npz"%n,allow_pickle=True)
    all_layer_corner_mag_raw = map_content['all_layer_corner_mag']
    all_layer_whether_densed_raw = map_content['all_layer_whether_densed']
    all_layer_sequence_number_in_next_layer_file_raw = map_content['all_layer_sequence_number_in_next_layer_file']
    layer_length_raw = map_content['layer_length']
    box_size = map_content['box_size'][0]
    
    nlayer = len(layer_length_raw)

    layer_length_accumulated_in_front_of_each_layer_raw = []
    accumulation = 0
    for i in range(nlayer):
        layer_length_accumulated_in_front_of_each_layer_raw.append(accumulation)
        accumulation += layer_length_raw[i] 

    layer_length_raw = list(map(lambda x:int(x+0.5),layer_length_raw))
    sum_layer_length = sum(layer_length_raw)
    layer_length_accumulated_in_front_of_each_layer_raw = list(map(lambda x:int(x+0.5),layer_length_accumulated_in_front_of_each_layer_raw))

    layer_length_accumulated_type = ctypes.c_int*nlayer
    layer_length_accumulated_in_front_of_each_layer = layer_length_accumulated_type()
    for i in range(nlayer):
        layer_length_accumulated_in_front_of_each_layer[i] = layer_length_accumulated_in_front_of_each_layer_raw[i]

    all_layer_corner_mag_type = ctypes.c_float*(4*sum_layer_length)
    all_layer_corner_mag = all_layer_corner_mag_type()
    for i in range(4*sum_layer_length):
        all_layer_corner_mag[i] = all_layer_corner_mag_raw[i]
    
    all_layer_whether_densed_type = ctypes.c_bool*(sum_layer_length)
    all_layer_whether_densed = all_layer_whether_densed_type()
    for i in range(sum_layer_length):
        all_layer_whether_densed[i] = all_layer_whether_densed_raw[i]

    #all_layer_sequence_number_in_next_layer_file_raw = list(map(lambda x:int(x+0.5),all_layer_sequence_number_in_next_layer_file_raw))
    
    all_layer_sequence_number_in_next_layer_file_type = ctypes.c_int*(sum_layer_length)
    all_layer_sequence_number_in_next_layer_file = all_layer_sequence_number_in_next_layer_file_type()
    for i in range(sum_layer_length):
        try:
            all_layer_sequence_number_in_next_layer_file[i] = int(all_layer_sequence_number_in_next_layer_file_raw[i]+0.5)
            
        except: 
            all_layer_sequence_number_in_next_layer_file[i] = 65535
            continue
    """
    magmap,magmap_sparse = RB.ReadBin(mapdir+'%d.dat'%n,filelength=sum(filelength))
   
    datatype = ctypes.c_float*filelength
    data = datatype()
    for n in range(filelength):
        data[n] = data_raw[n]
    datatype1 = ctypes.c_float*filelength1
    data1 = datatype1()
    for n in range(filelength1):
        data1[n] = data_raw[filelength+n]
    """
    
    #shifttype = ctypes.c_float*2
    #shift_c = shifttype()
    #shift_c[0] = shift_x
    #shift_c[1] = shift_y
    #nalpha = 1
    #alphalist = [315.]
    nalpha = 16
    alphalist = np.linspace(0,360.0-360.0/nalpha,nalpha)
    allparm = []
    allchi2 = []
    
    """jiyuan"""
    #q = parms[n][1]
    """end"""

    for alpha in alphalist:
        #parmbest = [t0,u0,te,alpha]
        tempparms = [t0,u0,te,alpha]
        parmbest,chi2min,iter,funcalls,warnflag,allevcs = fmin(chi2,tempparms,args=(data,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,box_size,False),full_output=True,retall=True,disp=0,maxiter=500,maxfun=1000)
        if use_mcmc==True:
            nburn_in = 500
            nsample = 1000
            ndim = len(parmbest)
            nwalkers = 2*ndim
            pos = [parmbest+1e-4*np.random.randn(ndim) for i in range(nwalkers)]
            
            #pos = []
            #for i in range(nwalkers):
            #    pos.append(np.array([parmbest[0]+parmbest[2]*0.01*np.random.randn(1)[0],parmbest[1]+parmbest[1]*0.01*np.random.randn(1)[0],parmbest[2]+parmbest[2]*0.05*np.random.randn(1)[0],parmbest[3]+parmbest[3]*0.01*np.random.randn(1)[0]]))

            sampler = emcee.EnsembleSampler(nwalkers,ndim,chi2,args=(data,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,box_size,True),threads=1)
            ## run EMCEE ##
            pos,lnprob,rstate = sampler.run_mcmc(pos,nburn_in)  ## burn-in
            sampler.reset()
            sampler.run_mcmc(pos,nsample)
            ## save EMCEE results ##
            chain = sampler.chain.reshape((-1,ndim),order='F')
            #fsfbs = np.array(sampler.blobs).reshape((-1,2*len(data)),order='F')
            chi2s = -2*sampler.lnprobability.reshape(-1,order='F')
            chi2min = min(chi2s)
            parmbest = chain[chi2s==chi2min][0]

        allparm.append(parmbest)
        allchi2.append(chi2min)
    allchi2 = np.array(allchi2)
    allparm = np.array(allparm)
    arg = np.argmin(allchi2)
    parmopt = allparm[arg]
    chi2opt = allchi2[arg]


    return chi2opt,parmopt[0],parmopt[1],parmopt[2],parmopt[3],parms[n][0],parms[n][1],parms[n][2]

if __name__ == '__main__':
    eventname = "ob161195"
    source_alpha = 269. 
    source_delta = -30.
    datadir = './data/%s/'%eventname
    datanames = ['KMTC01_I.pysis.dat','KMTC41_I.pysis.dat','KMTC42_I.pysis.dat','KMTA01_I.pysis.dat','KMTA41_I.pysis.dat','KMTA42_I.pysis.dat','KMTS01_I.pysis.dat','KMTS41_I.pysis.dat','KMTS42_I.pysis.dat']#'KMTC03_I_bin.pysis.dat','KMTA03_I_bin.pysis.dat','KMTS03_I_bin.pysis.dat']
    fluxnames = []#'KMTC03_I_bin.pysis.dat','KMTA03_I_bin.pysis.dat','KMTS03_I_bin.pysis.dat']
    #'KMTC01_I.pysis.dat','KMTC41_I.pysis.dat','KMTA01_I.pysis.dat','KMTA41_I.pysis.dat','KMTS01_I.pysis.dat','KMTS41_I.pysis.dat']
    errfac = [1.599,1.554,1.595,1.545,1.795,1.757,1.439,1.840,1.563]
    errsys = [0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003]
    inverse = [False,False,False,False,False,False,False,False,False]
    jd_to_hjd = [False,False,False,False,False,False,False,False,False]

    t0 = 7568.769
    u0 = 0.05317
    te = 9.96
    tbegin = 7540
    tend = 7600
    
    test = False

    use_mcmc =True 



    parms = np.load('parms_%s.npy'%map_set_name)
    #parms = np.load('parms_kb220371_all_q_close.npy')
#    args = range(len(parms))
    ##############################################
    ''' hongjing: pick out a fraction of parms '''
    args = []

    def select(ilogs,ilogq,ilogrho):
        ls,lq,lrho = round(ilogs,3),round(ilogq,3),round(ilogrho,3)
        #logs_range = [-0.034-0.025,-0.034+0.025]
        #logq_range = [-3.773-0.25,-3.773+0.25]
        #logrho_range = [-2.487-0.25,-2.487+0.25]
        
        logs_range   = [-0.005,0.014]
        logq_range   = [-4.7,-4.3]
        logrho_range = [-2.5,-2.5]
        logs_step = 0.001
        logq_step = 0.05
        logrho_step = 0.3

        #logs_range = [-0.25,0.]
        #logq_range = [-5.,-3.5]
        #logrho_range = [-4.,-2.5]
        #logs_step = 0.002
        #logq_step = 0.05
        #logrho_step = 0.3

        logss = np.round(np.arange(logs_range[0],logs_range[1]+0.5*logs_step,logs_step),5)
        logqs = np.round(np.arange(logq_range[0],logq_range[1]+0.5*logq_step,logq_step),5)
        logrhos = np.round(np.arange(logrho_range[0],logrho_range[1]+0.5*logrho_step,logrho_step),5)
        if   (ls not in logss):
            return False
        elif (lq not in logqs):
            return False
        elif (lrho not in logrhos):
            return False

        ### exclude a sub region ###
        #elif lq <= (-3.5/0.8*(ls+0.2)-5.5):
        #    return False
        #elif lq >= (-2.0/0.7*(ls-0.0)-2.5):
        #    return False
        ### end sub region ###

        else:
            return True
    """
    for i,pi in enumerate(parms):
        ilogs,ilogq,ilogrho = pi
        if select(ilogs,ilogq,ilogrho): # == True:
            args.append(i)
        else:
            continue
    """

    mapdir = "./map_set_%s/"%map_set_name
    
    for i,pi in enumerate(parms):
        ilogs,ilogq,ilogrho = pi
        try:
            f = np.load(mapdir+"%s.npz"%i,allow_pickle=True)
            #f.close()
            if select(ilogs,ilogq,ilogrho): # == True:
                args.append(i)
                #parm_select.append([ilogs,ilogq,ilogrho])
                #p.append([i,ilogs,ilogq,ilogrho])
            else:
                continue

        except:
            print('map: %s do not exist, so skips'%pi)
            continue
    
    '''hongjing: pick out a fraction of parms '''
    '''               - END -                 '''
    #############################################
    print('Number of grid: %d'%(len(args)))

    initial_type = ctypes.c_float*4

    data = []
    for i in range(len(datanames)):
        iname = datanames[i]
        tempdata = np.loadtxt(datadir+iname,usecols=(0,1,2))
        if tempdata[0,0]>2450000:
            tempdata[:,0]-=2450000
        if tempdata[0,0]>50000:
            tempdata[:,0]-=50000
        if jd_to_hjd[i] == True:
            hjd = []
            for t in tempdata[:,0]:
                hjd.append(pyasl.helio_jd(t+50000,source_alpha,source_delta)-50000)
            tempdata[:,0] = hjd
        arg = (tempdata[:,0]>tbegin)*(tempdata[:,0]<tend)
        tempdata = tempdata[arg]
        if len(tempdata)==0:
            print ("No data in %s satisfy the time domain"%iname)
            continue
        if inverse[i] == True:
            tempdata[:,1] *= -1.
        if iname not in fluxnames:
            tempdata[:,2] = errfac[i]*np.sqrt(tempdata[:,2]**2+errsys[i]**2)
            tempdata[:,1] = 10.**(0.4*(18.-tempdata[:,1]))
            tempdata[:,2] = tempdata[:,2]*tempdata[:,1]*np.log(10.)/2.5
        else:
            tempdata[:,2] = tempdata[:,2]*errfac[i]
        data.append(tempdata)
        print ("%s has %d data points" %(iname,len(tempdata)))

    data = np.array(data)
    print ("Start Grid Search")

    #pool = mp.Pool(processes=300)


    #start = time.clock() 
    if test == True:
        print (grid(3))
    else:
        pool = mp.Pool(processes=192)
        results = pool.map(grid,args)
        results = np.array(results,dtype='float')
        print (results)
        np.save("result/%s_adaptive_map_test.npy"%(eventname),results)
    #print ("Total time: %f"%(time.clock()-start))
    exit(0)
