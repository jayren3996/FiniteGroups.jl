include("$(@__DIR__)/../src/FiniteGroups.jl")
using .FiniteGroups, LinearAlgebra

py"""
from mpl_toolkits import mplot3d
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import time
from sympy.combinatorics import Permutation
from sympy.interactive import init_printing
import copy
import scipy
from scipy.stats import ortho_group

##——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————##
## def fundemental var and fun
sigma_0 = np.identity(2, dtype=complex);
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex);
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex);
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex);

sigma_1 = sigma_x
sigma_2 = sigma_y
sigma_3 = sigma_z

def direct_sum(rep1, rep2):
    _new_matrix = np.zeros(np.add(rep1.shape,rep2.shape), dtype = complex)
    _new_matrix[:rep1.shape[0],:rep1.shape[1]] = rep1
    _new_matrix[rep1.shape[0]:,rep1.shape[1]:] = rep2
    return _new_matrix 

Sigma = np.array([sigma_0, sigma_1, sigma_2, sigma_3])

#因为一些特殊的原因，这个地方的ij指标反了，由这样的函数生成出来的mass是和S = sigma_z sigma_0 sigma_0反对易的
def paulimatrix3(j,i,k):
    _matrix1 = np.kron(Sigma[i], Sigma[j]) 
    return np.kron(_matrix1, Sigma[k])

def paulimatrix2(i,j):
    return np.kron(Sigma[i], Sigma[j])

def chop(expr, delta=10**-5):
    return np.ma.masked_inside(expr, -delta, delta).filled(0)

##——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————##
Mass_term = np.array([paulimatrix3(0,1,0),
                      paulimatrix3(0,1,1),
                      paulimatrix3(0,1,3),
                      paulimatrix3(0,2,2),
                      
                      paulimatrix3(1,1,0),
                      paulimatrix3(1,1,1),
                      paulimatrix3(1,1,3),
                      paulimatrix3(1,2,2),
                      
                      paulimatrix3(2,1,2),
                      paulimatrix3(2,2,0),
                      paulimatrix3(2,2,1),
                      paulimatrix3(2,2,3),
                      
                      paulimatrix3(3,1,0),
                      paulimatrix3(3,1,1),
                      paulimatrix3(3,1,3),
                      paulimatrix3(3,2,2),
                      ])



def four_dim_rep_pointgroup(any_point_group_rep):
    if any_point_group_rep.shape[0] >=4:
        print('矩阵维度大于四，点群的线性表示的维度最高是3，所以出错')
        final_matrix = any_point_group_rep
    elif any_point_group_rep.shape[0] ==3:
        final_matrix = direct_sum(any_point_group_rep, np.array([[1]]))
    elif any_point_group_rep.shape[0] ==2:
        final_matrix = direct_sum(any_point_group_rep, sigma_0)
    elif any_point_group_rep.shape[0] ==1:
         final_matrix = direct_sum(any_point_group_rep, sigma_0)
         final_matrix = direct_sum(final_matrix, np.array([[1]]))
    else:
        final_matrix = np.kron(sigma_0, sigma_0)
    return final_matrix
        
        
def flavor_real_copro_rep(rep, dia_or_offdia):
    if dia_or_offdia == 1:
        _fla_rep = np.kron(sigma_0, rep)
    elif dia_or_offdia == -1:
        _fla_rep = np.kron(sigma_x, rep)
    else:
        print('错误的和S的对易关系，只能是1或者-1')
        _fla_rep = np.zeros((8,8))
    return _fla_rep

def mass_rep_fun(_PG_rep, dia_or_offdia):
    _four_dim_rep = four_dim_rep_pointgroup(_PG_rep)
    _eight_dim_copro_rep = flavor_real_copro_rep(_four_dim_rep, dia_or_offdia)
    _mass_rep = np.zeros((Mass_term.shape[0],Mass_term.shape[0]), dtype = complex)
    for j in range(Mass_term.shape[0]):
        _mass_after_transfor = np.dot(_eight_dim_copro_rep,np.dot(Mass_term[j],_eight_dim_copro_rep.conjugate().T))
        for i in range(Mass_term.shape[0]):
            _mass_rep[i,j]= np.trace(np.dot(Mass_term[i],_mass_after_transfor))
            
    return _mass_rep/8
  

def judge_comm_anticomm(m1,m2):
    if np.dot(m1,m2)-np.dot(m2,m1) == np.zeros(m1.shape):
        print
##——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————##
##transform the matrix form to "the array form" of permutation

def tran_Mat_to_Permu(matrix):
    _permu_array_form = list(range(matrix.shape[0]))
    for i in range(matrix.shape[0]):
        vec = np.zeros(matrix.shape[0])
        vec[i] = 1
        _permu_array_form[i]=np.nonzero(np.dot(matrix, vec))[0][0]
    return _permu_array_form

def index_to_conv(list_p):
    _list = copy.deepcopy(list_p)
    for i in range(len(_list)):
        for j in range(len(_list[i])):
            _list[i][j] +=1
    return _list



##——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————##

T_three_dim_rep = np.array([[0.5,-np.sqrt(2)/2,0.5],
                            [-np.sqrt(2)/2,0,np.sqrt(2)/2],
                            [0.5,np.sqrt(2)/2,0.5]])
# T_three_dim_rep = np.array([[0,0,1],
#                             [-1,0,0],
#                             [0,1,0]])
rep = T_three_dim_rep


##——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————##

# opera = np.load("Oh_7.npy")
# pairingsymm = np.load("Oh_pairing_4.npy")
# mass_rep_char =np.zeros(len(opera),dtype=complex)
# mass_rep_eigval=np.zeros(len(opera),dtype=object)
# for i in range(len(opera)): 
#     massrep = mass_rep_fun(opera[i],int(pairingsymm[i][0][0]))
#     mass_rep_char[i] = np.trace(massrep)
#     mass_rep_eigval[i],eigvec = np.linalg.eig(massrep)
#     print(chop(np.prod(mass_rep_eigval[i])))
#     print(mass_rep_eigval[i])
# np.save("Oh_7irrep_4Pairingsymm.npy",chop(mass_rep_char))





##——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————



##当mass在操作下会变成另一个mass而非一组mass的线性组合时，可以化成置换操作向轮换操作的分解问题。

# def mass_invar_space1():
#     _mass_rep = chop(np.real(mass_rep_fun(rep,-1)))

#     init_printing(perm_cyclic=True,pretty_print=False)
#     p = Permutation(tran_Mat_to_Permu(_mass_rep)).full_cyclic_form
#     p_conv = index_to_conv(p)
#     print('mass的不变子空间', p_conv)
#     print('———————————————————————————————————————————————————————————')

#     for i in range(len(p)):
#         if len(p[i])>4:
#             print('mass 的最小子空间维数大于4,')
#         for k in range(len(p[i])):
#             p[i][(k+1) %len(p[i])]
#             print('第%d个不变子空间的基的变换,第%d个mass-->第%d个mass的系数是%d+%dj'%(i+1, p_conv[i][k], p_conv[i][(k+1) %len(p[i])], _mass_rep[p[i][(k+1) %len(p[i])],p[i][k]].real, _mass_rep[p[i][(k+1) %len(p[i])],p[i][k]].imag))
#         for l in range(len(p[i])-1):
#             for j in range(len(p[i])-l-1):     
#                 if (np.dot(Mass_term[p[i][j]],Mass_term[p[i][j+l+1]])-np.dot(Mass_term[p[i][j+l+1]],Mass_term[p[i][j]])==np.zeros(Mass_term[p[i][j]].shape)).all():
#                     print('第%d个不变子空间里面的第%d个基和第%d个基是对易的'%(i+1,j+1,j+l+2))
#                 elif (np.dot(Mass_term[p[i][j]],Mass_term[p[i][j+l+1]])+np.dot(Mass_term[p[i][j+l+1]],Mass_term[p[i][j]])==np.zeros(Mass_term[p[i][j]].shape)).all():
#                     print('第%d个不变子空间里面的第%d个基和第%d个基是反对易的'%(i+1,j+1,j+l+2))
#                 else:
#                     print('第%d个不变子空间里面的第%d个基和第%d个基即不对易也不反对易'%(i+1,j+1,j+l+2))
#         print('———————————————————————————————————————————————————————————')
        
##——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
##当mass在操作下会变成另一个mass而非一组mass的线性组合时，可以化成置换操作向轮换操作的分解问题。


        
    
mass_rep = chop((mass_rep_fun(rep,-1)))
inv_space=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
for i in range(mass_rep.shape[0]):
    inv_space[i].insert(len(inv_space[i]), i)
    for j in range(mass_rep.shape[1]):
        if mass_rep[i][j] != 0:
            inv_space[i].insert(len(inv_space[i]), j)

invar_space = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
invar_space_nonoverlap=[]
invar_space_nonoverlap_before = []
invar_space_nonoverlap_after =inv_space
for overindex in range(100):   
    if invar_space_nonoverlap_before == invar_space_nonoverlap_after:
        break 
    else:
        invar_space_nonoverlap_before = invar_space_nonoverlap_after
        for i2 in range(len(invar_space_nonoverlap_before)):
            set_union = set(invar_space_nonoverlap_before[i2])  
            for j2 in range(len(invar_space_nonoverlap_before)):
                if set_union & set(invar_space_nonoverlap_before[j2])!=set():
                    set_union = set_union | set(invar_space_nonoverlap_before[j2])

            invar_space[i2]=list(np.sort(list(set_union)))
            invar_space_nonoverlap_after = invar_space
            
invar_space_nonoverlap = []  
for l in range(len(invar_space_nonoverlap_after)):
    if invar_space_nonoverlap_after[l] not in invar_space_nonoverlap:
        invar_space_nonoverlap.append(invar_space_nonoverlap_after[l])
            
        
            




    
p = invar_space_nonoverlap
p_conv = index_to_conv(p)
print('mass的不变子空间', p_conv)
print('———————————————————————————————————————————————————————————')
for i in range(len(p)):
    mass_invar_space_rep=np.zeros((len(p[i]),len(p[i])),dtype=float)
    for k in range(len(p[i])):
        for kk in range(len(p[i])):
            mass_invar_space_rep[k,kk] = np.real(mass_rep[p[i][k],p[i][kk]])
    print('第%d个不变子空间的基矢是mass'%i,p_conv[i],'mass不变子空间的变换矩阵是\n',mass_invar_space_rep)

    for l in range(len(p[i])-1):
        for j in range(len(p[i])-l-1):     
            if (np.dot(Mass_term[p[i][j]],Mass_term[p[i][j+l+1]])-np.dot(Mass_term[p[i][j+l+1]],Mass_term[p[i][j]])==np.zeros(Mass_term[p[i][j]].shape)).all():
                print('第%d个不变子空间里面的第%d个基和第%d个基是对易的'%(i+1,j+1,j+l+2))
            elif (np.dot(Mass_term[p[i][j]],Mass_term[p[i][j+l+1]])+np.dot(Mass_term[p[i][j+l+1]],Mass_term[p[i][j]])==np.zeros(Mass_term[p[i][j]].shape)).all():
                print('第%d个不变子空间里面的第%d个基和第%d个基是反对易的'%(i+1,j+1,j+l+2))
            else:
                print('第%d个不变子空间里面的第%d个基和第%d个基即不对易也不反对易'%(i+1,j+1,j+l+2))
        print('-----------------------------------------------------')
    print('———————————————————————————————————————————————————————————')
    
"""