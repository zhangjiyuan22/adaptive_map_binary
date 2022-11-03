import numpy as np
import ctypes
#import ReadBinary as RB
from scipy.optimize import fmin
import time
import multiprocessing as mp
import emcee
#from mpi4py import MPI
from MulensModel.binarylens import BinaryLens
from PyAstronomy import pyasl
Chi2Lib = ctypes.cdll.LoadLibrary('./chi2_calculator.so')

x = -0.0112
y = 0.

if __name__ == '__main__':

    mapdir = "ob161195_test"
    mapname = "0"
    
    map_index = int(mapname)
    #print(map_index)

    print('(x,y) = (%s,%s)'%(x,y))
    map_parms = np.load('parms_'+mapdir+'.npy') 
    #print(map_parms[0])
    
    print('log s = %s'%map_parms[map_index][0])
    print('log q = %s'%map_parms[map_index][1])
    print('log rho = %s'%map_parms[map_index][2])
    s = 10**(map_parms[map_index][0])
    q = 10**(map_parms[map_index][1])
    rho = 10**(map_parms[map_index][2])

    map_content = np.load('map_set_'+mapdir+"/%s.npz"%mapname,allow_pickle=True)
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
    
    A = 0.
    iAtype = ctypes.c_float
    iA = iAtype()
    ip = ctypes.byref(iA)

    box_size_type = ctypes.c_float
    box_size_final = box_size_type()
    box_size_final = box_size

    Chi2Lib.wrapinterpolating_to_get_magnification.argtypes = [ctypes.c_float,ctypes.c_float,ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_bool),ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.c_float,ctypes.POINTER(ctypes.c_float)]
    Chi2Lib.wrapinterpolating_to_get_magnification(x,y,all_layer_corner_mag,all_layer_whether_densed,all_layer_sequence_number_in_next_layer_file,layer_length_accumulated_in_front_of_each_layer,box_size_final,ip)

    print('A_interpolation = %s'%iA.value)

    bilens = BinaryLens(mass_1=1/(1.0+q), mass_2=q/(1.0 + q), separation=s)
    A_vbbl = bilens.vbbl_magnification( x , y , rho, u_limb_darkening=None, accuracy=1e-3)
    print('A_vbbl = %s'%A_vbbl)

    print('absolute error = %s'%(np.abs(iA.value-A_vbbl)))
    print('absolute error threshold = 0.01*sqrt(A) = %s'%(0.01*(A_vbbl**0.5)))
















