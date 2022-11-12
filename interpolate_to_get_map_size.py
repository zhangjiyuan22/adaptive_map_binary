%matplotlib notebook
    
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interpn


### the following is to set the max layer of maps ###

# This set of map' box_size = 3.5, which leads to the minimum subbox's side-length = 2*3.5/(2**max_layer)

# This set of map is generated at max_layer = 16, thus the value of max_layer can be set from 0 to 16 for test

# Under current box_size, max_layer =  15 leads to the minimum subbox's side-length = 2.1*10**(-4), 
# which is sufficient for most events. 

max_layer = 16

###                      end                      ###


all_map_layer_length_raw_list_new = np.loadtxt('all_map_layer_length_raw_list_new.txt', 
                                           dtype=int)
print(all_map_layer_length_raw_list_new)

sum_layer_length_list_new = []
for i_layer_length_raw_list in all_map_layer_length_raw_list_new :
    sum_layer_length_list_new.append( sum(i_layer_length_raw_list[0:(max_layer+1):1]) )

print('len(sum_layer_length_list_new) = %s'%len(sum_layer_length_list_new))


storage_size_from_calculate_list = []

for i_sum_layer_length in sum_layer_length_list_new :
    #storage_size_from_calculate_list.append( i_sum_layer_length * (64+4*64+1+16) / 8 / 1024 )  # k bytes
    storage_size_from_calculate_list.append( i_sum_layer_length * (0+4*32+1+32/4) / 8 / 1024 )  # k bytes
    

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111)

ax.plot(list(range(17019)), storage_size_from_calculate_list, linewidth = 1, label = 'from calculate', c = 'C1', marker = 'o')



sum_storage_size_from_calculate = sum(storage_size_from_calculate_list)/1024/1024   # giga bytes

ax.set_xlabel('map index',fontsize = 15)
ax.set_ylabel('storage size per map / kilo bytes',fontsize = 15)
ax.legend(fontsize = 15)

print('sum_storage_size_from_calculate = %s giga bytes'%sum_storage_size_from_calculate)


plt.show()



# the following is the 3D interpolation part

arr  = np.zeros( (31,61,9) )
for i in range(31):
    for j in range(61):
        for k in range(9):
            arr[i][j][k] = storage_size_from_calculate_list[ i*61*9 + j*9 + k ]
#print(arr)

logs_array = np.array(np.linspace(-1.5,0,31))
logq_array = np.array(np.linspace(-6,0,61))
logrho_array = np.array(np.linspace(-4,-1.6,9))
points = (logs_array, logq_array, logrho_array)
#print(len(logs_array)*len(logq_array)*len(logrho_array))




### the following is to set the (logs,logq,logrho) pairs at which interpolation is needed ###
# The format like [(logs1,logq1,logrho1), (logs2,logq2,logrho2), ...] is accepted.

xi = [ (-1.01, -3.01, -3.15), (-0.01, -0.01, -3) ]

###                                      end                                              ###

print('(logs,logq,logrho) pairs = %s'%xi)
result = interpn(points, arr, xi, method='linear', bounds_error=True)
print('interpolation result = %s kilo bytes'%result)
print('')





# the following is to calculate real existing cell number of each layer for a map by interpolation


### the following is to set the (logs,logq,logrho) pairs at which interpolation is needed ###
# The format like (logs,logq,logrho) is accepted.
yi = (-0.025, -0.05, -3.15)
###                                      end                                              ###


real_existing_cell_number_of_each_layer = []

for i_layer in range(max_layer+1) :
    arr_two  = np.zeros( (31,61,9) )
    for i in range(31):
        for j in range(61):
            for k in range(9):
                arr_two[i][j][k] = all_map_layer_length_raw_list_new[ i*61*9 + j*9 + k ][i_layer]
    result = interpn(points, arr_two, yi, method='linear', bounds_error=True)
    real_existing_cell_number_of_each_layer.append(result[0])

print('(logs,logq,logrho) pairs = %s'%[yi])
print('real existing cell number of each layer from interpolation = \n%s'%real_existing_cell_number_of_each_layer)

