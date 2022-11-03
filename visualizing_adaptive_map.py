%matplotlib notebook

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpathes

def getCaustic(separation,mass_ratio,npt=1000):
    masses = np.array([1.,mass_ratio])
    totalMass = sum(masses)
    masses /= totalMass
    print(masses)
    nlens = len(masses)
    offset = masses[1]*separation
    zlens1 = np.complex(-offset,0.)
    zlens2 = np.complex(separation-offset,0.)
    zlenses = np.array([zlens1,zlens2])
    
    phis = np.linspace(0,2*np.pi,npt+1)[:npt]
    zcauList,zcriList = [],[]
    for phi in phis:
        eiphi = np.exp(1j*phi)
        coeffs = np.zeros(5)*1j
        
        coeffs[0] = eiphi
        coeffs[1] = -2 * eiphi * (zlens1 + zlens2)
        coeffs[2] = ( (zlens1+zlens2)**2 + 2*zlens1*zlens2 )*eiphi - (masses[0]+masses[1])
        coeffs[3] = 2*(masses[0]*zlens2+masses[1]*zlens1) - 2*zlens1*zlens2*(zlens1 + zlens2)*eiphi
        coeffs[4] = eiphi*(zlens1**2)*(zlens2**2) - (masses[0]*(zlens2**2) - masses[1]*(zlens1**2))
        
        zcri = np.roots(coeffs) # four solutions
        zcau = np.zeros_like(zcri)*1j
        for eachRoot in range(len(zcri)):
            zcau[eachRoot] = zcri[eachRoot]
            for ilens in range(nlens):
                zcau[eachRoot] -= masses[ilens]/(np.conj(zcri[eachRoot])-np.conj(zlenses[ilens]))
        zcauList.extend(zcau)
        zcriList.extend(zcri)
    return np.array(zcauList),np.array(zcriList)



if __name__ == '__main__':

    path = "/work/zhangjiyuan/adaptive_map_binary_github/"
    mapdir = "ob161195_test"
    
    mapname = "0"
    
    
    map_index = int(mapname)
    map_parms = np.load(path+'parms_'+mapdir+'.npy')
    
    logs = map_parms[map_index][0]
    logq = map_parms[map_index][1]
    logrho = map_parms[map_index][2]
    print('log s = %s'%logs)
    print('log q = %s'%logq)
    print('log rho = %s'%logrho)

    
    map_content = np.load(path+'map_set_'+mapdir+"/%s.npz"%mapname,allow_pickle=True)
    all_layer_serial_number_raw = map_content['all_layer_serial_number']
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

    fig = plt.figure(dpi=180, figsize=(8,8))
    ax = fig.add_subplot(111)
    
    for current_layer in range(nlayer):

        current_nmesh = 2**current_layer
        current_step_length = 2*box_size / current_nmesh
        
        current_serial_number_sum = all_layer_serial_number_raw[ layer_length_accumulated_in_front_of_each_layer_raw[current_layer] : 
                                                                 layer_length_accumulated_in_front_of_each_layer_raw[current_layer] 
                                                               + layer_length_raw[current_layer] : 
                                                                 1 ]
        
        
        print('the following is layer = %s'%current_layer)
        
        print( 'number of total square       = %s'%( (2**current_layer) * (2**current_layer) ) )
        print( 'number of real exists square = %s'%len(current_serial_number_sum) )
        print('')
        
        for i in current_serial_number_sum :
            arg_y = i // current_nmesh
            arg_x = i % current_nmesh

            lower_left_y = ( arg_y ) * current_step_length - box_size
            lower_left_x = ( arg_x ) * current_step_length - box_size

            lower_left_quarter = np.array( [ lower_left_x , lower_left_y ] )
            rect = mpathes.Rectangle(lower_left_quarter, current_step_length , current_step_length , color='r',linewidth = 0.5, fill = 0)
            ax.add_patch(rect)
    
    npt = 5000
    
    zcauList,zcriList = getCaustic(10.**logs,10.**logq, npt)
    
    zcauList_rearrange = []
    zcauList_rearrange.append([zcauList[0],zcauList[1],zcauList[2],zcauList[3]])

    for i in range(1,npt):
        position1 = np.argmin(np.abs(zcauList_rearrange[i-1]-zcauList[4*i+0]))
        position2 = np.argmin(np.abs(zcauList_rearrange[i-1]-zcauList[4*i+1]))
        position3 = np.argmin(np.abs(zcauList_rearrange[i-1]-zcauList[4*i+2]))
        position4 = np.argmin(np.abs(zcauList_rearrange[i-1]-zcauList[4*i+3]))

        zcauList_i_temp = [0,0,0,0]
        zcauList_i_temp[position1] = zcauList[4*i+0]
        zcauList_i_temp[position2] = zcauList[4*i+1]
        zcauList_i_temp[position3] = zcauList[4*i+2]
        zcauList_i_temp[position4] = zcauList[4*i+3]

        zcauList_rearrange.append(zcauList_i_temp)

    print(len(zcauList_rearrange))

    zcauList_part1 = []
    zcauList_part2 = []
    zcauList_part3 = []
    zcauList_part4 = []
    for i in range(len(zcauList_rearrange)):
        zcauList_part1.append(zcauList_rearrange[i][0])
        zcauList_part2.append(zcauList_rearrange[i][1])
        zcauList_part3.append(zcauList_rearrange[i][2])
        zcauList_part4.append(zcauList_rearrange[i][3])
        
    #ax.scatter( np.real(zcauList) , np.imag(zcauList) , s = 0.1 , marker = '.' ,color = 'blue')
    ax.plot( np.real(zcauList_part1) , np.imag(zcauList_part1)  ,color = 'blue')
    ax.plot( np.real(zcauList_part2) , np.imag(zcauList_part2)  ,color = 'blue')
    ax.plot( np.real(zcauList_part3) , np.imag(zcauList_part3)  ,color = 'blue')
    ax.plot( np.real(zcauList_part4) , np.imag(zcauList_part4)  ,color = 'blue')

    plt.axis('equal')
    plt.title('logs=%s, logq=%s, logrho=%s'%(logs,logq,logrho))
    plt.show()
        
        

    

