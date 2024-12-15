import numpy as np
import matplotlib.pyplot as plt

N, M = 25, 25    #количество шагов по x,t
X, T = 1., 1.   #границы расчета
t0 = 0.5        #точки, в которых проверяется сходимость
x0 = 0.5
def solve(N,M,X,T):
    h = X/(N-1)    #шаг по  x,t
    tau = T/(M-1)
    eps=0.0001     #невязка
    
    u = np.zeros((N,M),dtype=float)  #сетка
    x = np.linspace(0,X,N)
    t = np.linspace(0,T,M)
    
    u[:,0] = 4/np.pi * np.arctan(x-2) + 2         #начальные и граничные условия
    u[0,:] = (2 - 4/np.pi * np.arctan(2))*np.exp(-t)
    
    def f(x11,x10,x01,x00):   #разностная схема
        return (x11-x10)/2/tau + (x10**2-x00**2)/4/h + (x01-x00)/2/tau + (x11**2-x01**2)/4/h
        
    def df(x11):
        return 1/2/tau + x11/2/h

    def step(x10,x01,x00):    #схема бегущего счета
        y1 = x00
        d = eps + 1
        while (d > eps):
            y0=y1
            y1 = y0 - f(y0,x10,x01,x00)/df(y0)
            d = abs(y1-y0)
        return y1

    for i in range(1,N):
        for j in range (1,M):
            u[i,j] = step(u[i,j-1],u[i-1,j],u[i-1,j-1])
    return u,x,t,h,tau

u1,x1,t1,h1,tau1 = solve(N,M,X,T)
u2,x2,t2,h2,tau2 = solve(2*N,2*M,X,T)
u3,x3,t3,h3,tau3 = solve(4*N,4*M,X,T)
u4,x4,t4,h4,tau4 = solve(8*N,8*M,X,T)

fig = plt.figure()           #график для N,M
ax = fig.gca(projection='3d')
T1, X1 = np.meshgrid(t1, x1)

surf = ax.plot_surface(X1, T1, u1)
#surf = ax.plot_surface(T1, X1, u1)

fig_ = plt.figure()         #график для 8N.8M
ax_ = fig_.gca(projection='3d')
T4, X4 = np.meshgrid(t4, x4)

surf_ = ax_.plot_surface(T4, X4, u4)

u1_2 = [u1[i,int(t0/tau1)] for i in range(0,N)]   #сходимость
u2_2 = [u2[i,int(t0/tau2)] for i in range(0,2*N)]
u3_2 = [u3[i,int(t0/tau3)] for i in range(0,4*N)]
u4_2 = [u4[i,int(t0/tau4)] for i in range(0,8*N)]


fig1 = plt.figure()
ax1 = fig1.gca()
ax2 = fig1.gca()
ax3 = fig1.gca()
ax4 = fig1.gca()

surf1 = ax1.plot(x1, u1_2)
surf2 = ax2.plot(x2, u2_2)
surf3 = ax3.plot(x3, u3_2)
surf4 = ax4.plot(x4, u4_2)
plt.show()

u1_0 = [u1[int(x0/h1),i] for i in range(0,M)]    # сходимость
u2_0 = [u2[int(x0/h2),i] for i in range(0,2*M)]
u3_0 = [u3[int(x0/h3),i] for i in range(0,4*M)]
u4_0 = [u4[int(x0/h4),i] for i in range(0,8*M)]


fig2 = plt.figure()
ax11 = fig2.gca()
ax21 = fig2.gca()
ax31 = fig2.gca()
ax41 = fig2.gca()

surf11 = ax11.plot(t1, u1_0)
surf21 = ax21.plot(t2, u2_0)
surf31 = ax31.plot(t3, u3_0)
surf41 = ax41.plot(t4, u4_0)
plt.show()


u4_b = (4/np.pi) * np.arctan(x4-2) + 2         #проверка ну и гу
u4_g = (2 - 4/np.pi * np.arctan(2))*np.exp(-t4)

print(u4[:,0] - u4_b)
print(u4[0,:] - u4_g)
