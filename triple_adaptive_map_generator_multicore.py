import numpy as np
from MulensModel.binarylens import BinaryLens
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpathes

import multiprocessing as mp
import time
import sys
import ctypes

import TripleLensing
TRIL = TripleLensing.TripleLensing()

#s = 0.93#0.8
#q = 1.2e-4#0.1
#rho = 1e-2
box_size = 2.0#0.4#10.0#0.8

#bilens = BinaryLens(mass_1=1/(1.0+q), mass_2=q/(1.0 + q), separation=s)


def densing_a_square( current_layer , current_serial_number , current_corner_mag , mlens3,zlens3,zlensv2 , rho , shift_x , shift_y , threshold_coefficient ):
    
    current_nmesh = 2**current_layer
    current_step_length = 2.*box_size / current_nmesh
    
    mag_interpolation_center = 0.25 * ( current_corner_mag[0] + current_corner_mag[1] + current_corner_mag[2] + current_corner_mag[3] )
    mag_interpolation_up     =  0.5 * ( current_corner_mag[2] + current_corner_mag[3] )
    mag_interpolation_down   =  0.5 * ( current_corner_mag[0] + current_corner_mag[1] )
    mag_interpolation_right  =  0.5 * ( current_corner_mag[1] + current_corner_mag[3] )
    mag_interpolation_left   =  0.5 * ( current_corner_mag[0] + current_corner_mag[2] )
    
    if abs(mag_interpolation_center) > 10000 or abs(mag_interpolation_up) > 10000 or abs(mag_interpolation_down) > 10000 or abs(mag_interpolation_right) > 10000 or abs(mag_interpolation_left) > 10000 :
        return [] , [] , False , False 
    #the second False means this map is wrong due to the mistake of VBBL
    
    arg_y = current_serial_number // current_nmesh
    arg_x = current_serial_number % current_nmesh

    y = ( arg_y + 0.5 ) * current_step_length - box_size + shift_y
    x = ( arg_x + 0.5 ) * current_step_length - box_size + shift_x
 
    next_nmesh = current_nmesh * 2
    next_step_length = 0.5 * current_step_length

    #mag_vbbl_center = bilens.vbbl_magnification( x , y , rho , u_limb_darkening=None , accuracy=1e-3 )
    #mag_vbbl_up     = bilens.vbbl_magnification( x , y + next_step_length , rho , u_limb_darkening=None , accuracy=1e-3 )
    #mag_vbbl_down   = bilens.vbbl_magnification( x , y - next_step_length , rho , u_limb_darkening=None , accuracy=1e-3 )
    #mag_vbbl_right  = bilens.vbbl_magnification( x + next_step_length , y , rho , u_limb_darkening=None , accuracy=1e-3 )
    #mag_vbbl_left   = bilens.vbbl_magnification( x - next_step_length , y , rho , u_limb_darkening=None , accuracy=1e-3 )

    mag_vbbl_center = TRIL.TriLightCurve(mlens3, zlensv2, [x - zlens3[4]] , [y - zlens3[5]] , rho, 45, 2, 1e-6, 1e-3)
    mag_vbbl_up     = TRIL.TriLightCurve(mlens3, zlensv2, [x - zlens3[4]] , [y + next_step_length - zlens3[5]] , rho, 45, 2, 1e-6, 1e-3)
    mag_vbbl_down   = TRIL.TriLightCurve(mlens3, zlensv2, [x - zlens3[4]] , [y - next_step_length - zlens3[5]] , rho, 45, 2, 1e-6, 1e-3)
    mag_vbbl_right  = TRIL.TriLightCurve(mlens3, zlensv2, [x + next_step_length - zlens3[4]] , [y - zlens3[5]] , rho, 45, 2, 1e-6, 1e-3)
    mag_vbbl_left   = TRIL.TriLightCurve(mlens3, zlensv2, [x - next_step_length - zlens3[4]] , [y - zlens3[5]] , rho, 45, 2, 1e-6, 1e-3)

    #if np.abs( mag_interpolation - mag_vbbl_center ) > threshold_coefficient * pow( mag_vbbl_center , 0.5 ) :
    if np.abs( mag_interpolation_center - mag_vbbl_center[0] ) > threshold_coefficient * pow( mag_vbbl_center[0] , 0.5 ) or np.abs( mag_interpolation_up - mag_vbbl_up[0] ) > threshold_coefficient * pow( mag_vbbl_up[0] , 0.5 ) or np.abs( mag_interpolation_down - mag_vbbl_down[0] ) > threshold_coefficient * pow( mag_vbbl_down[0] , 0.5 ) or np.abs( mag_interpolation_right - mag_vbbl_right[0] ) > threshold_coefficient * pow( mag_vbbl_right[0] , 0.5 ) or np.abs( mag_interpolation_left - mag_vbbl_left[0] ) > threshold_coefficient * pow( mag_vbbl_left[0] , 0.5 ) or current_layer <= 4:
        
        next_serial_number_0 = 2 * arg_y * next_nmesh + 2 * arg_x
        next_serial_number_1 = 2 * arg_y * next_nmesh + ( 2 * arg_x + 1 )
        next_serial_number_2 = ( 2 * arg_y + 1 ) * next_nmesh + 2 * arg_x
        next_serial_number_3 = ( 2 * arg_y + 1 ) * next_nmesh + ( 2 * arg_x + 1 )

        #mag_vbbl_up    = bilens.vbbl_magnification( x + shift_x , y + next_step_length + shift_y , rho, u_limb_darkening=None, accuracy=1e-3)
        #mag_vbbl_down  = bilens.vbbl_magnification( x + shift_x , y - next_step_length + shift_y , rho, u_limb_darkening=None, accuracy=1e-3)
        #mag_vbbl_right = bilens.vbbl_magnification( x + next_step_length + shift_x , y + shift_y , rho, u_limb_darkening=None, accuracy=1e-3)
        #mag_vbbl_left  = bilens.vbbl_magnification( x - next_step_length + shift_x , y + shift_y , rho, u_limb_darkening=None, accuracy=1e-3)
    
        next_corner_mag_0 = [ current_corner_mag[0] , mag_vbbl_down[0] , mag_vbbl_left[0] , mag_vbbl_center[0] ]
        next_corner_mag_1 = [ mag_vbbl_down[0] , current_corner_mag[1] , mag_vbbl_center[0] , mag_vbbl_right[0] ]
        next_corner_mag_2 = [ mag_vbbl_left[0] , mag_vbbl_center[0] , current_corner_mag[2] , mag_vbbl_up[0] ]
        next_corner_mag_3 = [ mag_vbbl_center[0] , mag_vbbl_right[0] , mag_vbbl_up[0] , current_corner_mag[3] ]

        return [ next_serial_number_0 , next_serial_number_1 , next_serial_number_2 , next_serial_number_3 ] , next_corner_mag_0 + next_corner_mag_1 + next_corner_mag_2 + next_corner_mag_3 , True , True
    
    else:

        return [] , [] , False , True
'''
def once(i):
    current_serial_number = current_serial_number_sum[i]
    current_corner_mag = [ current_corner_mag_sum[4*i] , current_corner_mag_sum[4*i+1] , current_corner_mag_sum[4*i+2] , current_corner_mag_sum[4*i+3] ]
    return densing_a_square( current_layer , current_serial_number , current_corner_mag , mlens3,zlens3,zlensv2 , rho , shift_x , shift_y , threshold_coefficient = 1e-2 )
'''

from multiprocessing import Pool


def generating_next_layer( current_layer , current_serial_number_sum , current_corner_mag_sum , mlens3,zlens3,zlensv2 , rho , shift_x , shift_y ):
    
    next_serial_number_sum = []
    next_corner_mag_sum = []

    current_whether_densed_sum = []
    current_sequence_number_in_next_layer_file_sum = []

    current_sequence_number_in_next_layer_file = 0

    global worker

    def worker(i):
        current_serial_number = current_serial_number_sum[i]
        current_corner_mag = [ current_corner_mag_sum[4*i] , current_corner_mag_sum[4*i+1] , current_corner_mag_sum[4*i+2] , current_corner_mag_sum[4*i+3] ]
        res = densing_a_square( current_layer , current_serial_number , current_corner_mag , mlens3,zlens3,zlensv2 , rho , shift_x , shift_y , threshold_coefficient = 1e-2 )
        print("worker", i, res[1])
        return res

    p = Pool(processes= min(len(current_serial_number_sum), 128) )
    mags_3L1S = p.map(worker, list(range(len(current_serial_number_sum))))
    p.close()


    for i in range( len(current_serial_number_sum) ):
        
        current_serial_number = current_serial_number_sum[i]
        current_corner_mag = [ current_corner_mag_sum[4*i] , current_corner_mag_sum[4*i+1] , current_corner_mag_sum[4*i+2] , current_corner_mag_sum[4*i+3] ]

        #next_serial_number , next_corner_mag , current_whether_densed , whether_correct_map = densing_a_square( current_layer , current_serial_number , current_corner_mag , mlens3,zlens3,zlensv2 , rho , shift_x , shift_y , threshold_coefficient = 1e-2 )
        next_serial_number , next_corner_mag , current_whether_densed , whether_correct_map = mags_3L1S[i]
        if whether_correct_map == False :
            return next_serial_number_sum , next_corner_mag_sum , current_whether_densed_sum , current_sequence_number_in_next_layer_file_sum , False
        
        next_serial_number_sum.extend( next_serial_number )
        next_corner_mag_sum.extend( next_corner_mag )

        current_whether_densed_sum.append( current_whether_densed )
        if current_whether_densed :
            current_sequence_number_in_next_layer_file_sum.append( current_sequence_number_in_next_layer_file )
            current_sequence_number_in_next_layer_file += 1
        else :
            current_sequence_number_in_next_layer_file_sum.append( None )
        
    return next_serial_number_sum , next_corner_mag_sum , current_whether_densed_sum , current_sequence_number_in_next_layer_file_sum , True

### the following is to draw the caustics###
def getCaustic(separation,mass_ratio,npt=1000):
    masses = np.array([1.,mass_ratio])
    totalMass = sum(masses)
    masses /= totalMass
    nlens = len(masses)
    offset = masses[1]*separation
    zlens1 = np.complex(-offset,0.)
    zlens2 = np.complex(separation-offset,0.)
    zlenses = np.array([zlens1,zlens2])
    ######
    f0 = np.zeros(2*nlens+1)*1j
    gc = np.zeros([nlens,2*nlens])*1j
    fc = np.zeros([nlens,2*nlens])*1j
    phis = np.linspace(0,2*np.pi,npt+1)[:npt]
    zcauList,zcriList = [],[]
    for phi in phis:
        f0[0] = zlenses[0]**2
        f0[1] = -2.*zlenses[0]
        f0[2] = 1.
        k = 1
        for ilens in range(1,nlens):
            k += 2
            f0[k+1] = f0[k-1]
            f0[k] = f0[k-2]-2.*f0[k-1]*zlenses[ilens]
            for j in range(k-2,0,-1):
                f0[j+1] = f0[j-1]-2.*f0[j]*zlenses[ilens]+f0[j+1]*zlenses[ilens]**2
            f0[1] = -2.*f0[0]*zlenses[ilens]+f0[1]*zlenses[ilens]**2
            f0[0] = f0[0]*zlenses[ilens]**2
        for ilens in range(nlens):
            gc[ilens,2*nlens-1] = f0[2*nlens]
            for j in range(2*nlens,1,-1):
                gc[ilens,j-2] = gc[ilens,j-1]*zlenses[ilens]+f0[j-1]
            fc[ilens,2*nlens-2] = gc[ilens,2*nlens-1]
            for j in range(2*nlens-1,1,-1):
                fc[ilens,j-2] = fc[ilens,j-1]*zlenses[ilens]+gc[ilens,j-1]
        hc = np.zeros(2*nlens+1)*1j
        eiphi = np.exp(1j*phi)
        hc[2*nlens] = f0[2*nlens]*eiphi
        hc[2*nlens-1] = f0[2*nlens-1]*eiphi
        for order in range(2*nlens-1,0,-1):
            hc[order-1] = f0[order-1]*eiphi
            secondTerm = 0.
            for ilens in range(nlens):
                secondTerm += masses[ilens]*fc[ilens,order-1]
            hc[order-1] -= secondTerm
        orders = 2*nlens+1
        coeffs = np.zeros(orders)*1j
        for ith in range(orders):
            coeffs[ith] = hc[orders-ith-1]
        zcri = np.roots(coeffs)
        zcau = np.zeros_like(zcri)*1j
        for eachRoot in range(len(zcri)):
            zcau[eachRoot] = zcri[eachRoot]
            for ilens in range(nlens):
                zcau[eachRoot] -= masses[ilens]/(np.conj(zcri[eachRoot])-np.conj(zlenses[ilens]))
        zcauList.extend(zcau)
        zcriList.extend(zcri)
    return np.array(zcauList),np.array(zcriList)
"""
if __name__ == '__main__':

    fig = plt.figure(dpi=180, figsize=(8,8))
    ax = fig.add_subplot(111)
    
    current_layer = 0
    current_serial_number_sum = [0]
     
    current_corner_mag_0 = bilens.vbbl_magnification( -box_size , -box_size , rho, u_limb_darkening=None, accuracy=1e-3)
    current_corner_mag_1 = bilens.vbbl_magnification( box_size , -box_size , rho, u_limb_darkening=None, accuracy=1e-3)
    current_corner_mag_2 = bilens.vbbl_magnification( -box_size , box_size , rho, u_limb_darkening=None, accuracy=1e-3)
    current_corner_mag_3 = bilens.vbbl_magnification( box_size , box_size , rho, u_limb_darkening=None, accuracy=1e-3)
    
    current_corner_mag_sum = [ current_corner_mag_0 , current_corner_mag_1 , current_corner_mag_2 , current_corner_mag_3 ]
    
    while(1):
        
        print('the following is layer = %s'%current_layer)
        #print(current_serial_number_sum)
        print( 'number of total square       = %s'%( (2**current_layer) * (2**current_layer) ) )
        print( 'number of real exists square = %s'%len(current_serial_number_sum) )
        #print(current_corner_mag_sum)

        current_nmesh = 2**current_layer
        current_step_length = 2*box_size / current_nmesh

        for i in current_serial_number_sum :
            arg_y = i // current_nmesh
            arg_x = i % current_nmesh
            
            lower_left_y = ( arg_y ) * current_step_length - box_size
            lower_left_x = ( arg_x ) * current_step_length - box_size

            lower_left_quarter = np.array( [ lower_left_x , lower_left_y ] )
            rect = mpathes.Rectangle(lower_left_quarter, current_step_length , current_step_length , color='r',linewidth = 0.5, fill = 0)
            ax.add_patch(rect)
        
        if current_layer == 10 :
            break

        next_serial_number_sum , next_corner_mag_sum , current_whether_densed_sum , current_sequence_number_in_next_layer_file_sum = generating_next_layer( current_layer , current_serial_number_sum , current_corner_mag_sum )

        print( 'whether densed = %s'%current_whether_densed_sum )
        print( 'sequence number in next layer file = %s\n'%current_sequence_number_in_next_layer_file_sum )

        current_layer = current_layer + 1
        current_serial_number_sum = next_serial_number_sum
        current_corner_mag_sum = next_corner_mag_sum

    zcauList,zcriList = getCaustic( s , q )
    #print(zcauList)
    ax.scatter( np.real(zcauList) , np.imag(zcauList) , s = 0.1 , marker = '.' ,color = 'blue')
    ### end ###

    plt.axis('equal')

    plt.show()


"""

if __name__ == '__main__':

    #fig = plt.figure(dpi=180, figsize=(8,8))
    #ax = fig.add_subplot(111)

    #print(arg)
    #print(parm[arg])
    #logs1,logq1,logrho = parm[arg]
    
    prms = np.load('lens_1_param.npz')
    mlens3=prms['mlens3'];
    zlens3=prms['zlens3'];

    zlensv2 = [0] * 6
    zlensv2[0] = zlens3[0] - zlens3[4]
    zlensv2[2] = zlens3[2] - zlens3[4]
    zlensv2[4] = zlens3[4] - zlens3[4]
    zlensv2[1] = zlens3[1] - zlens3[5]
    zlensv2[3] = zlens3[3] - zlens3[5]
    zlensv2[5] = zlens3[5] - zlens3[5]

    #q = 10.**logq1
    #s = 10.**logs1
    rho = 1e-3
    
    #shift_x = (-1.)*(s-1./s)*q/(1.+q)
    shift_x = 0.
    shift_y = 0.
    
    #bilens = BinaryLens(mass_1=1/(1.0+q), mass_2=q/(1.0 + q), separation=s)
    
    all_layer_serial_number = []
    all_layer_corner_mag = []
    all_layer_whether_densed = []
    all_layer_sequence_number_in_next_layer_file = []
    layer_length = []

    current_layer = 0
    current_serial_number_sum = [0]
     
    #current_corner_mag_0 = bilens.vbbl_magnification( -box_size + shift_x , -box_size + shift_y , rho, u_limb_darkening=None, accuracy=1e-3)
    #current_corner_mag_1 = bilens.vbbl_magnification( box_size + shift_x , -box_size + shift_y , rho, u_limb_darkening=None, accuracy=1e-3)
    #current_corner_mag_2 = bilens.vbbl_magnification( -box_size + shift_x , box_size + shift_y , rho, u_limb_darkening=None, accuracy=1e-3)
    #current_corner_mag_3 = bilens.vbbl_magnification( box_size + shift_x , box_size + shift_y , rho, u_limb_darkening=None, accuracy=1e-3)
   

    current_corner_mag_0 = TRIL.TriLightCurve(mlens3, zlensv2, [-box_size + shift_x - zlens3[4]] , [-box_size + shift_y - zlens3[5]] , rho, 45, 2, 1e-6, 1e-3)
    current_corner_mag_1 = TRIL.TriLightCurve(mlens3, zlensv2, [box_size + shift_x - zlens3[4]] , [-box_size + shift_y - zlens3[5]] , rho, 45, 2, 1e-6, 1e-3)
    current_corner_mag_2 = TRIL.TriLightCurve(mlens3, zlensv2, [-box_size + shift_x - zlens3[4]] , [box_size + shift_y - zlens3[5]] , rho, 45, 2, 1e-6, 1e-3)
    current_corner_mag_3 = TRIL.TriLightCurve(mlens3, zlensv2, [box_size + shift_x - zlens3[4]] , [box_size + shift_y - zlens3[5]] , rho, 45, 2, 1e-6, 1e-3)

    current_corner_mag_sum = [ current_corner_mag_0[0] , current_corner_mag_1[0] , current_corner_mag_2[0] , current_corner_mag_3[0] ]
    
    all_layer_serial_number.extend(current_serial_number_sum)
    all_layer_corner_mag.extend(current_corner_mag_sum)
    layer_length.append(len(current_serial_number_sum))

    while(1):
        
        print('the following is layer = %s'%current_layer)
        #print(current_serial_number_sum)
        print( 'number of total square       = %s'%( (2**current_layer) * (2**current_layer) ) )
        print( 'number of real exists square = %s'%len(current_serial_number_sum) )
        print(max(current_corner_mag_sum))
        
        current_nmesh = 2**current_layer
        current_step_length = 2*box_size / current_nmesh
        """
        for i in current_serial_number_sum :
            arg_y = i // current_nmesh
            arg_x = i % current_nmesh
            
            lower_left_y = ( arg_y ) * current_step_length - box_size
            lower_left_x = ( arg_x ) * current_step_length - box_size

            lower_left_quarter = np.array( [ lower_left_x , lower_left_y ] )
            rect = mpathes.Rectangle(lower_left_quarter, current_step_length , current_step_length , color='r',linewidth = 0.5, fill = 0)
            ax.add_patch(rect)
        """
        if current_layer == 14 :
            current_whether_densed_sum_final = []
            current_sequence_number_in_next_layer_file_sum_final = []
            for i in range(len(current_serial_number_sum)):
                current_whether_densed_sum_final.append(False)
                current_sequence_number_in_next_layer_file_sum_final.append(None)
            all_layer_whether_densed.extend(current_whether_densed_sum_final)
            all_layer_sequence_number_in_next_layer_file.extend(current_sequence_number_in_next_layer_file_sum_final)

            break

        next_serial_number_sum , next_corner_mag_sum , current_whether_densed_sum , current_sequence_number_in_next_layer_file_sum , whether_correct_map_2 = generating_next_layer( current_layer , current_serial_number_sum , current_corner_mag_sum ,mlens3,zlens3,zlensv2,rho, shift_x, shift_y )
        
        #if whether_correct_map_2 == False :
        #    print('map %s is wrong'%parm[arg])
        #    return

        all_layer_serial_number.extend(next_serial_number_sum)
        all_layer_corner_mag.extend(next_corner_mag_sum)
        all_layer_whether_densed.extend(current_whether_densed_sum)
        all_layer_sequence_number_in_next_layer_file.extend(current_sequence_number_in_next_layer_file_sum)
        layer_length.append(len(next_serial_number_sum))

        #print( 'whether densed = %s'%current_whether_densed_sum )
        #print( 'sequence number in next layer file = %s\n'%current_sequence_number_in_next_layer_file_sum )
        
        current_layer = current_layer + 1
        current_serial_number_sum = next_serial_number_sum
        current_corner_mag_sum = next_corner_mag_sum
    """
    zcauList,zcriList = getCaustic( s , q )
    #print(zcauList)
    ax.scatter( np.real(zcauList) , np.imag(zcauList) , s = 0.1 , marker = '.' ,color = 'blue')
    ### end ###
    
    plt.axis('equal')

    plt.show()
    """
    #print(all_layer_serial_number)
    #print(all_layer_corner_mag)
    #print(all_layer_whether_densed)
    #print(all_layer_sequence_number_in_next_layer_file)
    box_size_array = np.array([box_size])
    all_layer_serial_number_array = np.array(all_layer_serial_number)
    all_layer_corner_mag_array = np.array(all_layer_corner_mag)
    all_layer_whether_densed_array = np.array(all_layer_whether_densed)
    all_layer_sequence_number_in_next_layer_file_array = np.array(all_layer_sequence_number_in_next_layer_file)
    layer_length_array = np.array(layer_length)

    np.savez('./map_set_triple/triple_map',all_layer_serial_number=all_layer_serial_number_array , all_layer_corner_mag=all_layer_corner_mag_array , all_layer_whether_densed=all_layer_whether_densed_array , all_layer_sequence_number_in_next_layer_file=all_layer_sequence_number_in_next_layer_file_array , layer_length=layer_length_array , box_size=box_size_array)
    
"""
def saveparm():
    logsstep = 0.05
    #logsstart,logsend = -0.033,-0.033
    logsstart,logsend = 0.,1.5
    #logsstart,logsend = -0.035,-0.03

    logqstep = 0.1
    #logqstart,logqend = -3.842,-3.842
    logqstart,logqend = -5.5,1.
    #logqstart,logqend = -4.2,-3.6
    #logqstart,logqend = -3.8,-3.75

    logrhostep = 0.5
    #logrhostart,logrhoend = -2.558,-2.558
    logrhostart,logrhoend = -3.5,-2.5
    #logrhostart,logrhoend = -2.6,-2.6

    # alphas = np.linspace(0,180,19)
    allparm = []
    logs1 = logsstart
    n=0
    while(logs1<=logsend):
        logq1 = logqstart
        while(logq1<=logqend):
            logrho = logrhostart
            while(logrho<=logrhoend):
                allparm.append([logs1,logq1,logrho])
                n += 1
                logrho = round(logrho+logrhostep,5)
            logq1 = round(logq1+logqstep,5)
        logs1 = round(logs1+logsstep,5)
    allparm = np.array(allparm)
    #np.save('parms_for_lc_contrast.npy',allparm)
    np.save('parms_kb190505_all_q_large_s.npy',allparm)
    return n
def loadparm(parm,n):
    #tempparm = np.load('parms_for_lc_contrast.npy')
    tempparm = np.load('parms_kb190505_all_q_large_s.npy')
    for i in range(n):
        for j in range(3):
            parm[i][j] = float(tempparm[i,j])

if __name__ == '__main__':
    n = saveparm()
    print ('number of map = %d'%n)

    parm = mp.Array(ctypes.c_float, n*3)
    parm = np.ctypeslib.as_array(parm.get_obj())
    parm = parm.reshape(n,3)
    loadparm(parm,n)

#    start = time.clock()
    pool = mp.Pool(processes=132)
    args = range(len(parm))
    pool.map(generating_a_map,args)
    pool.close()
    pool.join()



#    print ("Total time: %f"%(time.clock()-start))
"""
