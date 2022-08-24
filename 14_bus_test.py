import numpy as np
import pandas as pd
from dlpf import DLPF
from fdlpf import FDLPF
import time
import matplotlib.pyplot as plt
import pandapower as pp
import pandapower.networks as pn

net=pn.case14()
pp.runpp(net)
print(net['res_bus'])


branch_data=pd.read_excel('14_bus_test.xlsx',sheet_name=0)
bus_data = pd.read_excel('14_bus_test.xlsx', sheet_name=1)
print(branch_data)
print(bus_data)
start=time.time()
test_obj=DLPF(branch_data=branch_data,bus_data=bus_data)
#test_obj=FDLPF(branch_data=branch_data,bus_data=bus_data)#FDLPF结果存在nan
test_obj.rundlpf()
end=time.time()
test_obj.show_result()
print(test_obj.bus_result)
print(test_obj.pf_result)

print('DLPF运算时间',end-start)


# start=time.time()
# test_obj=FDLPF(branch_data=branch_data,bus_data=bus_data)
# #test_obj=FDLPF(branch_data=branch_data,bus_data=bus_data)#FDLPF结果存在nan
# test_obj.rundlpf()
# end=time.time()
# test_obj.show_result()
# print(test_obj.bus_result)
# print(test_obj.pf_result)
#
# print('FDLPF运算时间',end-start)

plt.figure()
plt.title('v comparison')
x=test_obj.bus_result['num'];y=test_obj.bus_result['v']
plt.plot(x,y)
x=[];y=[]
for i in range(len(net['res_bus'])):
    x.append(i+1)
    temp=net['res_bus'].loc[i,'vm_pu']
    y.append(temp)
plt.plot(x,y)
plt.legend(['dlpf','pandapower'])

plt.show()


