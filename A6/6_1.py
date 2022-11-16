import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
import math
import scipy
import pandas as pd
def alpha1(x,y):
    dy_dx= y[1]
    dy1_dx=-2*(-(x**2)/(2)+((em)))*y[0]
    return [dy_dx,dy1_dx]
    
def sub_plot(ax,a,b,d,title):
    ax.scatter(a,b,label="Numerical Value")
    ax.plot(a,d,label="Inbuilt Solution")
    ax.set_title(title)
    ax.set_xlabel("\u03BE")
    ax.set_ylabel("u(\u03BE)")
    ax.legend()
    ax.grid(True)

def numerov(xi,xf,N,e,n):
    N = N + 1
    x=np.linspace(xi,xf,N)
    h = (xf-xi)/(N)
    alpha= 2*(-(x**2)/(2) + e)
    C = []
    global em 
    em =e 
    for i in range(N):
        C_i = 1+(h**2/12)*alpha[i]
        C.append(C_i)
    U = np.zeros((N))
    if n%2==0:
        u_0=1
        u_1= (6 - 5*C[0])/C[1]
        u_prime=0
    else:
        u_0=0
        u_1= h
        u_prime=1
    U[0] = u_0
    U[1] = u_1
    for i in range(1,N-1):
        U[i+1] = (1/C[i+1])*((12-10*C[i])*U[i]-C[i-1]*U[i-1])
    a=integrate.simps(U**2,x)
    N=U/(np.sqrt(a)) 
    sol=integrate.solve_ivp(alpha1,[xi,xf],[u_0,u_prime] , dense_output=True)   
    b=sol.sol(x)
    c=integrate.simps(np.array(b[0]**2),x)
    N1=b[0]/(np.sqrt(c))
    return N,x ,alpha , N1

def e_range(u,n_node,E_min,E_max) :
    I = []
    E = (E_min+E_max)/2
    for i in range(len(u)):
        if (u[i-1]*u[i]) < 0:
           I.append(i)
    N_node = len(I)
    if N_node > int((n_node)/2):
       E_max = E
    else:
       E_min = E
    return len(I),E_min,E_max

def E(xi,xf,N,n_node,E_min,E_max,n,tol):
    for i in range(10):
        U ,x,alpha ,analitic  = numerov(xi,xf,N,(E_max+E_min)/2,n)
        I ,E_min_new,E_max_new = e_range(U,n_node,E_min,E_max)
        if abs(E_max_new - E_min_new)<tol:
            break
        else:
            E_min = E_min_new
            E_max = E_max_new
    if n/2 == 1: E_min_new =E_min + 4
    return E_min_new,E_max_new,U,alpha ,analitic
eigen_value =[]
analytical=[]
n=[0,1,2,3,4,5]
for i in n:
    E_min,E_max,U,alpha , analytic= E(0,4,100,i,0,8,i,10**(-10))
    anal=(i+1/2)
    analytical.append(anal)
    eigen_value.append((E_min+E_max)/2)
data={"Eigen Value Numerically" : eigen_value , "Analytically": analytical}
print(pd.DataFrame(data))
plt.scatter(n,eigen_value,label='n vs E_n')
plt.grid()
plt.legend()
plt.xlabel(" N")
plt.xlabel(" Eigen Value")
plt.title("E_n as function of n")
plt.show()
slope, intercept, r, p, se=scipy.stats.linregress(n,eigen_value)
print("...............................")
print( " The Slope of the curve fitted graph is : ",slope)
print("...............................")
print( " The Intercept of the curve fitted graph is : ",intercept)

plt.plot(np.power(n,2),eigen_value )
plt.xlabel(" N_square")
plt.ylabel(" Eigen Value")
plt.grid()
plt.title("E_n as function of n")
plt.show()
slope, intercept, r, p, se=scipy.stats.linregress(np.power(n,2),eigen_value)
print("...............................")
print( " The Slope of the curve fitted graph is : ",slope)
print("...............................")
print( " The Intercept of the curve fitted graph is : ",intercept)

def parity(n,U, b):
    if n%2==0:
        U1=U[::-1]
        U2=U1[1:]
        b1=b[::-1]
        b2=b1[1:]
    elif n%2!= 0:
        U1=-U[::-1]
        U2=U1[1:]
        b1=-b[::-1]
        b2=b1[1:]
    total=np.concatenate((U2,U))
    tot_inbuillt=np.concatenate((b2,b))
    return total ,tot_inbuillt


fig2=plt.figure(figsize=(10,10))
fig2, ((axx1, axx2), (axx3, axx4), (axx5, axx6)) = plt.subplots(3,2)
fig2.suptitle('u(x) vs x for N intervals')
dict = {'0': axx1,'1': axx2,'2': axx3,'3': axx4,'4': axx5,'5': axx6}
fig2.delaxes(axx6)
for j in range(5):
    (E_min,E_max,U,alpha , analytic),a= E(0,5,100,j,0,8,j,10**(-10)),1
    if j==2 or j==3:
        sub_plot(dict[str(j)],np.linspace(-5,5,201),-parity(j,U,analytic)[0],-parity(j,U,analytic)[1],"For n =" +str(j))
    else:
        if j==1 : a=-1
        sub_plot(dict[str(j)],np.linspace(-5,5,201),parity(j,U,analytic)[0] ,a*parity(j,U,analytic)[1],"For n =" +str(j))
plt.tight_layout()
plt.show()

for j in range(5):
    (E_min,E_max,U,alpha , analytic),a= E(0,4,100,j,0,8,j,10**(-10)),1
    if j==2 or j==3:
        plt.scatter(np.linspace(-4,4,201),(-parity(j,U,analytic)[0])**2 , label='Numerical ')

        plt.plot(np.linspace(-4,4,201),(-parity(j,U,analytic)[1])**2 , label="Inbuilt") 
    else:
        if j==1 : a=-1 
        plt.scatter(np.linspace(-4,4,201),(parity(j,U,analytic)[0] )**2, label='Numerical')
        plt.plot(np.linspace(-4,4,201),(a*parity(j,U,analytic)[1])**2 , label="Inbuilt")
    plt.title("Probabilty Density for n =" +str(j))
    plt.grid()
    plt.xlabel("\u03BE")
    plt.ylabel("u**2(\u03BE)")
    plt.legend()
    plt.show()

#For x_max=10
fig2=plt.figure(figsize=(10,10))
fig2, ((axx1, axx2), (axx3, axx4), (axx5, axx6)) = plt.subplots(3,2)
fig2.suptitle('u(x) vs x for N intervals')
dict = {'0': axx1,'1': axx2,'2': axx3,'3': axx4,'4': axx5,'5': axx6}
fig2.delaxes(axx6)
for j in range(5):
    (E_min,E_max,U,alpha , analytic),a= E(0,10,100,j,0,8,j,10**(-10)),1
    if j==2 or j==3:
        sub_plot(dict[str(j)],np.linspace(-10,10,201),-parity(j,U,analytic)[0],-parity(j,U,analytic)[1],"For n =" +str(j))
    else:
        if j==1 : a=-1
        sub_plot(dict[str(j)],np.linspace(-10,10,201),parity(j,U,analytic)[0] ,a*parity(j,U,analytic)[1],"For n =" +str(j))
plt.tight_layout()
plt.show()

for j in range(5):
    (E_min,E_max,U,alpha , analytic),a= E(0,10,100,j,0,8,j,10**(-10)),1
    if j==2 or j==3:
        plt.scatter(np.linspace(-10,10,201),(-parity(j,U,analytic)[0])**2 , label='Numerical ')
        plt.plot(np.linspace(-10,10,201),(-parity(j,U,analytic)[1])**2 , label="Inbuilt") 
    else:
        if j==1 : a=-1 
        plt.scatter(np.linspace(-10,10,201),(parity(j,U,analytic)[0] )**2, label='Numerical')
        plt.plot(np.linspace(-10,10,201),(a*parity(j,U,analytic)[1])**2 , label="Inbuilt")
    plt.title("Probabilty Density for n =" +str(j))
    plt.grid()
    plt.legend()
    plt.show()

#For x_max=50

fig2=plt.figure(figsize=(10,10))
fig2, ((axx1, axx2), (axx3, axx4), (axx5, axx6)) = plt.subplots(3,2)
fig2.suptitle('u(x) vs x for N intervals')
dict = {'0': axx1,'1': axx2,'2': axx3,'3': axx4,'4': axx5,'5': axx6}
fig2.delaxes(axx6)
for j in range(5):
    (E_min,E_max,U,alpha , analytic),a= E(0,50,100,j,0,8,j,10**(-10)),1
    if j==2 or j==3:
        sub_plot(dict[str(j)],np.linspace(-50,50,201),-parity(j,U,analytic)[0],-parity(j,U,analytic)[1],"For n =" +str(j))
    else:
        if j==1 : a=-1
        sub_plot(dict[str(j)],np.linspace(-50,50,201),parity(j,U,analytic)[0] ,a*parity(j,U,analytic)[1],"For n =" +str(j))
plt.tight_layout()
plt.show()

for j in range(5):
    (E_min,E_max,U,alpha , analytic),a= E(0,50,100,j,0,8,j,10**(-10)),1
    if j==2 or j==3:
        plt.scatter(np.linspace(-50,50,201),(-parity(j,U,analytic)[0])**2 , label='Numerical ')
        plt.plot(np.linspace(-50,50,201),(-parity(j,U,analytic)[1])**2 , label="Inbuilt") 
    else:
        if j==1 : a=-1 
        plt.scatter(np.linspace(-50,50,201),(parity(j,U,analytic)[0] )**2, label='Numerical')
        plt.plot(np.linspace(-50,50,201),(a*parity(j,U,analytic)[1])**2 , label="Inbuilt")
    plt.title("Probabilty Density for n =" +str(j))
    plt.grid()
    plt.legend()
    plt.xlabel("\u03BE")
    plt.ylabel("u**2(\u03BE)")
    plt.show()
