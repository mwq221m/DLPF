import numpy as np
import pandas as pd
import time

class DLPF():
    def __init__(self,branch_data,bus_data):
        self.version='1.0.2'
        self.branch_data=branch_data
        self.bus_data=bus_data
        self.bus_num=len(self.bus_data)
        self.branch_num=len(self.branch_data)
        self.Y = np.zeros((self.bus_num, self.bus_num),dtype='complex')
        self.Y_without_shunt = np.zeros((self.bus_num, self.bus_num),dtype='complex')
        self.G=np.zeros((self.bus_num,self.bus_num),dtype='complex')
        self.B = np.zeros((self.bus_num, self.bus_num),dtype='complex')
        self.B_without_shunt = np.zeros((self.bus_num, self.bus_num),dtype='complex')
        self.R_bus_list=[]
        self.S_bus_list = []
        self.L_bus_list = []
        for i in range(self.bus_num):
            temp=self.bus_data.iloc[i]
            type_str=temp['type']
            num_str=temp['num']
            if type_str=='R':
                self.R_bus_list.append(int(num_str))
            if type_str=='S':
                self.S_bus_list.append(int(num_str))
            if type_str=='L':
                self.L_bus_list.append(int(num_str))
        self.R_bus_num=len(self.R_bus_list)
        self.S_bus_num = len(self.S_bus_list)
        self.L_bus_num = len(self.L_bus_list)
        self.R_bus_index_list=np.array(self.R_bus_list)-1
        self.S_bus_index_list = np.array(self.S_bus_list) - 1
        self.L_bus_index_list = np.array(self.L_bus_list) - 1
        self.p_tilde=np.zeros(self.S_bus_num+self.L_bus_num)
        self.q_tilde=np.zeros(self.L_bus_num)
        self.H=np.zeros((self.S_bus_num+self.L_bus_num,self.S_bus_num+self.L_bus_num),dtype='complex')
        self.N=np.zeros((self.S_bus_num+self.L_bus_num,self.L_bus_num),dtype='complex')
        self.M=np.zeros((self.L_bus_num,self.S_bus_num+self.L_bus_num),dtype='complex')
        self.L=np.zeros((self.L_bus_num,self.L_bus_num),dtype='complex')

        self.p_s=[]
        self.p_l=[]
        self.q_l=[]
        self.theta_r=[]
        self.v_r=[]
        self.v_s=[]
        for i in range(self.bus_num):
            temp=self.bus_data.iloc[i]
            type=temp['type']
            if type=='R':
                self.v_r.append(temp['v'])
                self.theta_r.append(temp['theta']*np.pi/180)
            if type=='S':
                self.v_s.append(temp['v'])
                self.p_s.append(temp['p'])
            if type=='L':
                self.p_l.append(temp['p'])
                self.q_l.append(temp['q'])

        self.SL_list=self.S_bus_list+self.L_bus_list




        for i in range(self.branch_num):
            temp=self.branch_data.iloc[i]
            if temp['status']==1:
                #print('*****')
                start_idx=int(temp['start'])-1
                end_idx=int(temp['end'])-1
                r=float(temp['r'])
                x=float(temp['x'])
                z=r+x*1j
                y=1/z
                half_b=temp['b/2']
                transformer=temp['transformer']
                if transformer==0:
                    self.Y_without_shunt[start_idx,start_idx]+=y
                    self.Y_without_shunt[start_idx,end_idx]+=-y
                    self.Y_without_shunt[end_idx,start_idx]+=-y
                    self.Y_without_shunt[end_idx,end_idx]+=y
                    self.Y[start_idx, start_idx] += y+half_b*(1j)
                    self.Y[start_idx, end_idx] += -y
                    self.Y[end_idx, start_idx] += -y
                    self.Y[end_idx, end_idx] += y+half_b*(1j)
                else:
                    '''导纳矩阵生成时考虑到变压器的影响'''
                    t=temp['t']
                    theta=temp['theta']
                    if np.isnan(half_b):
                        half_b=0#拿14节点测试时发现变压器支路忘记将该变量置零 逻辑上也容易忽视 所以在这里自动修正
                    self.Y_without_shunt[start_idx, start_idx] += y*1/t**2
                    self.Y_without_shunt[start_idx, end_idx] += -y*(1/(t*np.exp(-1j*theta)))
                    self.Y_without_shunt[end_idx, start_idx] += -y*(1/(t*np.exp(1j*theta)))
                    self.Y_without_shunt[end_idx, end_idx] += y
                    self.Y[start_idx, start_idx] +=( y + half_b*(1j))*1/t**2
                    self.Y[start_idx, end_idx] += -y*(1/(t*np.exp(-1j*theta)))
                    self.Y[end_idx, start_idx] += -y*(1/(t*np.exp(1j*theta)))
                    self.Y[end_idx, end_idx] += y + half_b*(1j)
                #self.G=np.real(self.Y_without_shunt)
        #self.B_without_shunt=np.imag(self.Y_without_shunt)
                #self.Y=self.Y_without_shunt.copy()

        self.B_without_shunt = np.imag(self.Y_without_shunt)
        for i in range(self.bus_num):
            temp=self.bus_data.iloc[i]
            shunt=temp['shunt']
            '''shunt不包括支路电纳，已在支路信息里增加b/2'''
            if shunt==1:
                g=temp['g']
                b=temp['b']
                y=g+b*1j
                self.Y[i,i]+=y
        self.G=np.real(self.Y)
        self.B=np.imag(self.Y)

        for i in range(self.S_bus_num+self.L_bus_num):#H矩阵形成
            for j in range(self.S_bus_num+self.L_bus_num):
                i_idx=int(self.SL_list[i]-1)
                j_idx=int(self.SL_list[j]-1)
                self.H[i,j]=-self.B_without_shunt[i_idx,j_idx]
                #print(-self.B_without_shunt[i_idx,j_idx])

        for i in range(self.S_bus_num+self.L_bus_num):#N矩阵
            for j in range(self.L_bus_num):
                i_idx = int(self.SL_list[i] - 1)
                j_idx=int(self.L_bus_list[j]-1)
                self.N[i,j]=self.G[i_idx,j_idx]

        for i in range(self.L_bus_num):#M
            for j in range(self.S_bus_num+self.L_bus_num):
                i_idx=int(self.L_bus_list[i]-1)
                j_idx = int(self.SL_list[j] - 1)
                self.M[i,j]=-self.G[i_idx,j_idx]

        for i in range(self.L_bus_num):#L
            for j in range(self.L_bus_num):
                i_idx = int(self.L_bus_list[i] - 1)
                j_idx = int(self.L_bus_list[j] - 1)
                self.L[i,j]=-self.B[i_idx,j_idx]

    def rundlpf(self):
        pq_temp=self.p_s+self.p_l+self.q_l
        pq_temp=np.array(pq_temp).reshape(-1,1)
        v_theta_temp=self.theta_r+self.v_r+self.v_s
        v_theta_temp=np.array(v_theta_temp).reshape(-1,1)
        array_temp=np.zeros((self.S_bus_num+2*self.L_bus_num,2*self.R_bus_num+self.S_bus_num),dtype='complex')
        '''公式12中的矩阵可分为四块'''
        for i in range(self.S_bus_num+self.L_bus_num):
            for j in range(self.R_bus_num):
                i_idx=int(self.SL_list[i]-1)
                j_idx=int(self.R_bus_list[j]-1)
                array_temp[i,j]=self.B_without_shunt[i_idx,j_idx]

        self.RS_list=self.R_bus_list+self.S_bus_list
        for i in range(self.S_bus_num+self.L_bus_num):
            for j in range(self.R_bus_num+self.S_bus_num):
                i_idx = int(self.SL_list[i] - 1)
                j_idx=int(self.RS_list[j]-1)
                array_temp[i,j+self.R_bus_num]=-self.G[i_idx,j_idx]

        for i in range(self.L_bus_num):
            for j in range(self.R_bus_num):
                i_idx=int(self.L_bus_list[i]-1)
                j_idx=int(self.R_bus_list[j]-1)
                array_temp[i+self.S_bus_num+self.L_bus_num,j]=self.G[i_idx,j_idx]

        for i in range(self.L_bus_num):
            for j in range(self.R_bus_num+self.S_bus_num):
                i_idx=int(self.L_bus_list[i]-1)
                j_idx=int(self.RS_list[j]-1)
                array_temp[i+self.S_bus_num+self.L_bus_num,j+self.R_bus_num]=self.B[i_idx,j_idx]

        pq_tilde=pq_temp+array_temp@v_theta_temp
        #print(array_temp)
        #print(pq_temp)
        #print(v_theta_temp)
        #print(pq_tilde)
        self.p_tilde=pq_tilde[:self.S_bus_num+self.L_bus_num,0]
        self.q_tilde=pq_tilde[self.S_bus_num+self.L_bus_num:,0]

        self.H_tilde=self.H-self.N@np.linalg.inv(self.L)@self.M
        self.L_tilde=self.L-self.M@np.linalg.inv(self.H)@self.N

        self.theta_tilde=np.linalg.inv(self.H_tilde)@self.p_tilde-np.linalg.inv(self.H_tilde)@self.N@np.linalg.inv(self.L)@self.q_tilde
        self.theta_tilde*=180/np.pi
        self.theta_tilde=np.real(self.theta_tilde)

        self.v_tilde=np.linalg.inv(self.L_tilde)@self.q_tilde-np.linalg.inv(self.L_tilde)@self.M@np.linalg.inv(self.H)@self.p_tilde
        self.v_tilde=np.real(self.v_tilde)

    def show_result(self):
        temp_list=[]
        #v_list=[]#使用列表添加只能按照rsl顺序编号节点 不然节点电压没问题但潮流功率有问题
        #theta_list=[]
        v_dic={}
        theta_dic={}
        for i in range(self.bus_num):
            temp=self.bus_data.iloc[i]
            type=temp['type']
            num=temp['num']
            if type=='R':
                temp_list.append({'num':num,'type':type,'v':temp['v'],'theta':temp['theta']})
                #v_list.append(temp['v'])
                #theta_list.append(temp['theta'])
                v_dic[num]=temp['v']
                theta_dic[num]=temp['theta']
            if type=='S':
                idx=np.where(np.array(self.S_bus_list)==num)[0][0]
                temp_list.append({'num':num,'type':type,'v':temp['v'],'theta':self.theta_tilde[idx]})
                #v_list.append(temp['v'])
                #theta_list.append(self.theta_tilde[idx])
                v_dic[num]=temp['v']
                theta_dic[num]=self.theta_tilde[idx]
            if type=='L':
                idx = np.where(np.array(self.L_bus_list) == num)[0][0]
                temp_list.append({'num':num,'type':type,'v':self.v_tilde[idx],'theta':self.theta_tilde[idx+self.S_bus_num]})
                #v_list.append(self.v_tilde[idx])
                #theta_list.append(self.theta_tilde[idx+self.S_bus_num])
                v_dic[num]=self.v_tilde[idx]
                theta_dic[num]=self.theta_tilde[idx+self.S_bus_num]
        self.bus_result=pd.DataFrame(temp_list)#生成dataframe数据也耗费时间

        temp_list=[]

        for i in range(self.branch_num):
            temp=self.branch_data.iloc[i]
            status=temp['status']
            if status==1:
                start=int(temp['start'])
                end=int(temp['end'])
                transformer=temp['transformer']
                start_idx=int(start-1)
                end_idx=int(end-1)
                r=temp['r']
                x=temp['x']
                z=r+x*1j
                y=1/z
                g=np.real(y)
                b=np.imag(y)
                #p=g*(v_list[start_idx]-v_list[end_idx])-b*(theta_list[start_idx]-theta_list[end_idx])*np.pi/180
                #q=-g*(theta_list[start_idx]-theta_list[end_idx])*np.pi/180-b*(v_list[start_idx]-v_list[end_idx])
                if transformer==0:
                    p=g*(v_dic[start]-v_dic[end])-b*(theta_dic[start]-theta_dic[end])*np.pi/180
                    q=-g*(theta_dic[start]-theta_dic[end])*np.pi/180-b*(v_dic[start]-v_dic[end])
                else:
                    t=temp['t']
                    theta=temp['theta']
                    p = g/t * (v_dic[start]/t - v_dic[end]) - b/t * (theta_dic[start] - theta_dic[end]-theta) * np.pi / 180
                    q = -g/t * (theta_dic[start] - theta_dic[end]-theta) * np.pi / 180 - b/t * (v_dic[start]/t - v_dic[end])
                temp_list.append({'start':start,'end':end,'p':p,'q':q})
        self.pf_result=pd.DataFrame(temp_list)





















if __name__=='__main__':
    #test_data=pd.read_excel('test.xlsx')
    #print(test_data)
    branch_data=pd.read_excel('test.xlsx',sheet_name=0)
    bus_data = pd.read_excel('test.xlsx', sheet_name=1)
    print(branch_data)
    print(bus_data)
    start = time.time()#大部分时间来自于文件读取
    test_obj=DLPF(branch_data=branch_data,bus_data=bus_data)
    #print(test_obj.Y)
    #print(test_obj.B)
    #print(test_obj.R_bus_list)
    #print(test_obj.R_bus_num)
    #print(test_obj.L)
    #print(test_obj.N)
    #print(test_obj.M)
    #print(test_obj.H)
    #print(test_obj.p_s)
    #print(test_obj.v_s)
    test_obj.rundlpf()
    #print(test_obj.p_tilde)
    #print(test_obj.q_tilde)
    #print('theta',test_obj.theta_tilde)
    #print('mag',test_obj.v_tilde)
    test_obj.show_result()
    print(test_obj.bus_result)
    print(test_obj.pf_result)
    end=time.time()
    print('运行时间',end-start)