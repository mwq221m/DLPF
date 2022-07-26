import numpy as np
import pandas as pd
from dlpf import DLPF
import time


branch_data=pd.read_excel('6_bus_test.xlsx',sheet_name=0)
bus_data = pd.read_excel('6_bus_test.xlsx', sheet_name=1)
start=time.time()
test_obj=DLPF(branch_data=branch_data,bus_data=bus_data)
test_obj.rundlpf()
test_obj.show_result()
print(test_obj.bus_result)
print(test_obj.pf_result)
end=time.time()
print('运算时间',end-start)