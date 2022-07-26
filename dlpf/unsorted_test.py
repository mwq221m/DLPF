import numpy as np
from dlpf import DLPF
import pandas as pd
'''
不按照文献中rsl顺序排列节点测试
将文件中节点顺序逆序排列 并将1和2节点替换
'''

branch_data=pd.read_excel('unsorted_test.xlsx',sheet_name=0)
bus_data = pd.read_excel('unsorted_test.xlsx', sheet_name=1)
test_obj=DLPF(branch_data=branch_data,bus_data=bus_data)
test_obj.rundlpf()
test_obj.show_result()
print(test_obj.bus_result)
print(test_obj.pf_result)