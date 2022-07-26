import pandapower as pp
import pandapower.networks as pn
import numpy as np
import time
#net=pp.networks.case9()
net=pn.case14()
start=time.time()
pp.runpp(net)
end=time.time()
print('运算时间',end-start)
print(net['bus'])
print(net['line'])
print(net['res_bus'])
pp.to_excel(net,filename='ieee14_case.xlsx')