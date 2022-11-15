from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.special import assoc_laguerre
import pandas as pd

def Numerov(IC,a,b,N,alp,l):
    alp = np.vectorize(alp)
    x = np.linspace(a,b,int(N)+1)
    h = x[1] - x[0]
    u = np.zeros(len(x))    
    c = (1 + ((h**2/12)*alp(x)))
    if IC[0] == 1: 
       IC[1] = (6-5*c[0])/c[1]
    else:
        IC[1] = h**(l+1)
    u[:2] = IC
    for i in range(2,len(x)):    
        u[i] = (1/c[i])*(((12-(10*c[i-1]))*u[i-1]) - (c[i-2]*u[i-2]))   
    return x,u
def solve_eigen(x_min,x_max,e,l):    
    def slope(x):
        k = (e + (2/x) - (l*(l+1)/(x**2)))
        return k
    x,si_o = Numerov([0,0],x_min,x_max,1000,slope,l)
    a = integrate.simps(np.power(si_o,2),x)
    nom_sk = si_o/np.sqrt(abs(a))   
    return x,nom_sk
def eigen_val(n,l):
    x_min = 1e-14 ; x_max = 5*n**2
    ep = np.linspace(-1,0,100)
    cmp = np.zeros(len(ep))
    ind = []
    for i in range(len(ep)):
        x,Si = solve_eigen(x_min,x_max,ep[i],l)
        cmp[i] = Si[-1]
          
    for i in range(len(cmp)-1):
        if cmp[i]*cmp[i+1] <= 0:
           ind.append(i)    
    p = n-l-1       
    e0 = ep[ind[p]] ; e1 = ep[ind[p]+1]
    return e0,e1

def plotting1(x,y,title,key,key2,key3,key4,s,func = None,):
    #key4 = 1 -- n  and key4 = 0 -- l
    dic2 = {1:'n',0:'l'}
    m = np.shape(y)[0]
    if key == 'p':
       fig,ax = plt.subplots()
       ax.set_title(title)
    for i in range(m):
        dic = {0:f'for {dic2[key4]} = {i+key4}',1:'Numerical'}
        if key == 'w':
           fig,ax = plt.subplots()
           ax.set_title(f'{title} for {dic2[abs(key4-1)]} = {s} and {dic2[key4]} = {i+key4} State')
        ax.plot(x[i],y[i],'o',label=f'Wave Function {dic[key3]}')
        if key2 == 'y':
            if key4 == 1:
               ax.plot(x[i],func(x[i],i+1,s),label = 'Analytic')
            else:
                ax.plot(x[i],func(x[i],s,i),label = 'Analytic')
        ax.set_xlabel('x')
        ax.set_ylabel('U')
        ax.legend()
    plt.show() 
def Analytic(r,n,l):
    return ((2*r/n)**(l)*assoc_laguerre(2*r/n,n-l-1,2*l+1))/np.exp(r/n)*r

prob_an = lambda r,n,l : Analytic(r,n,l)**2
area = lambda r,n,l : integrate.simps(prob_an(r,n,l),r)
new_anl = lambda r,n,l : Analytic(r,n,l)/np.sqrt(area(r,n,l))
new_pro = lambda r,n,l : new_anl(r,n,l)**2 

def eigenshoot(n,l,tol):
    x_min = 1e-7 ; x_max = 5*n**2
    def phi(ep):
        res1,res2 = solve_eigen(x_min, x_max, ep,l)
        y_c = res2
        return y_c[-1],res1,res2
    e_min,e_max = eigen_val(n,l)

    for i in range(10000):
        e = (e_min + e_max)/2
        la,x,u = phi(e)
        
        if la < 0 :
           if (n-l)%2 ==1 :
              e_max = e
           else:
              e_min = e 
        else:
           if (n-l)%2 == 1:
              e_min = e
           else:
              e_max = e 
        if abs(la) < tol :
           return e,x,u
    return 0

eigenshoot = np.vectorize(eigenshoot,otypes = [float,np.ndarray,np.ndarray])

#table for energy
N = np.arange(1,6,1)

e,x,K = eigenshoot(N,0,1e-10)
e1,x1,u1 = eigenshoot(N[1:],1,1e-10)
e2,x2,u2 = eigenshoot(N[2:],2,1e-10)

e_z1 = np.concatenate(([0],e1))
e_z2 = np.concatenate(([0,0],e2)) 

data = {'state':N,'Energy l = 0':e,'Energy L = 1':e_z1,'Energy l = 2':e_z2}
df = pd.DataFrame(data)
print(df)
df.to_csv('table.csv')

#plotting of wv function
plotting1(x,K ,'wavefunction','w','y',1,1,0,new_anl)

#plotting of probablity density
def plotpdf(n):
    l = []
    for i in range(n):
        l.append(i)
    e3,x3,u_p = eigenshoot(n,l,1e-10)
    u_21 = np.power(u_p,2)
    plotting1(x3,u_21,f'Probability Density n = {n}','p','y',0,0,n,new_pro)

plotpdf(1)
plotpdf(2)
plotpdf(3)

#probablty
X = x[1]
U = K[1]
p = np.arange(0.5,20.5,0.5)
prob = np.zeros(len(p))
for i in range(len(p)):
    i1 = np.where(abs(X - p[i]) < 0.001)[0]    
    prob[i] = integrate.simps(np.power(U[:int(i1)+1],2),X[:int(i1)+1])
fig,ax = plt.subplots()
ax.plot(p,prob,'+')
ax.set_title('Probability vs Radius ')
ax.set_xlabel('r')
ax.set_ylabel('Probability')
plt.plot()