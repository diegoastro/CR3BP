#%% Sistema de ecuaciones y condiciones iniciales-------------------------------------------------------

import numpy as np
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from timeit import default_timer as timer
from matplotlib.animation import FuncAnimation

def mov(x, y, vx, vy, t):
    
  x1, x2 = rs[0], rsj[0]
  epsilon = 1e-6
  r1 = np.maximum(epsilon, (x1 - x) ** 2 + y ** 2) ** (3 / 2)
  r2 = np.maximum(epsilon, (x2 - x) ** 2 + y ** 2) ** (3 / 2)
  dvx = 2 * vy + x - (1 - m) * (x - x1) / r1 - m * (x - x2) / r2
  dvy = y - ((1 - m) / r1 + m / r2) * y - 2 * vx
  return np.array([vx, vy, dvx, dvy], float)


# Condiciones iniciales

m =  ((30 * 1.89e27) / (2e30 + 30*1.89e27))

l1 = [1-m - (m/3)**(1/3) + 0.02, 0.001] #l1 dentro de la esfera de Hill

l2 = [1-m + (m/3)**(1/3) + 0.01, 0.001] #l2 dentro de la esfera de Hill

l3 = [-(m+1) + 7/12*m - 0.01, 0.001] #l3

l4 = [0.5*(1-2*m) + 0.001, np.sqrt(3)/2] #l4

l5 = [0.5*(1-2*m) + 0.001, -np.sqrt(3)/2] #l5

vp = [0,0.001]

rs = [-m, 0]
rsj = [1 - m, 0]

#Numero de pasos para la integracion

t = 0
t_final = 5 * 2 * np.pi
dt = 0.00001
num_steps = int((t_final - t) / dt)

#%% RK4

x, y, vx, vy = l3[0], l3[1], vp[0], vp[1] #aqui se cambian las condiciones iniciales para ver las distintas trayectorias

# Inicializamos listas para posición y velocidad
xp = []
yp = []
vxp = []
vyp = []

start = timer()
for i in range(num_steps):
    
    xp.append(x)
    yp.append(y)
    
    vxp.append(vx)
    vyp.append(vy)
    
    k1 = dt * mov(x, y, vx, vy, t)
    k2 = dt * mov(x + 0.5 * k1[2], y + 0.5 * k1[3], vx + 0.5 * k1[0], vy + 0.5 * k1[1], t + 0.5 * dt)
    k3 = dt * mov(x + 0.5 * k2[2], y + 0.5 * k2[3], vx + 0.5 * k2[0], vy + 0.5 * k2[1], t + 0.5 * dt)
    k4 = dt * mov(x + k3[2], y + k3[3], vx + k3[0], vy + k3[1], t + dt)

    x += (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6
    y += (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6
    vx += (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) / 6
    vy += (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]) / 6

    
    t += dt

end = timer()
print('Tiempo de ejecución RK4', end - start)

#%% Verlet

x, y, vx, vy = l4[0], l4[1], vp[0], vp[1] #aqui se cambian las condiciones iniciales para ver las distintas trayectorias

xp_v = []
yp_v = []
vxp_v = []
vyp_v = []

start = timer()

for j in range(num_steps):
    
    xp_v.append(x)
    yp_v.append(y)
    
    vxp_v.append(vx)
    vyp_v.append(vy)
    
    result = mov(x, y, vx, vy, t)
    vx, vy, dvx, dvy = result[0], result[1], result[2], result[3]
    
    vx = vx + dvx * dt
    vy = vy + dvy * dt
    
    x = x + vx * dt
    y = y + vy * dt
    
    t += dt

end = timer()
print('Tiempo de ejecución Verlet:', end - start)

#%% Euler

x, y, vx, vy = l3[0], l3[1], vp[0], vp[1] #aqui se cambian las condiciones iniciales para ver las distintas trayectorias

xp_euler = []
yp_euler = []

start = timer()

for k in range(num_steps):
    
    xp_euler.append(x)
    yp_euler.append(y)
    
    dx, dy, dvx, dvy = mov(x, y, vx, vy, t)
    
    x += dt * dx
    y += dt * dy
    vx += dt * dvx
    vy += dt * dvy
    
    t += dt

end = timer()

print('Tiempo de ejecución Euler(RK1)', end - start)


#%%-------------------RK4-------------------------------------------------------------------------------
plt.figure()
plt.plot(xp, yp, color = 'blue', linewidth = 0.5 )
plt.plot(rsj[0], rsj[1], 'o', color='orange', label=r'Planeta ($30M_J$)', markersize = 10)
plt.plot(rs[0], rs[1], 'o', color='yellow', markeredgecolor="#FD7813", label='Sol', markersize = 15)
plt.plot(l3[0], l3[1], '+', color = 'black', label = r'$\mathcal{L}_3$')#Cambiar aqui segun que condicion inicial se este usando
plt.xlabel(r'$X$')
plt.ylabel(r'$Y$')
#plt.legend(loc = 'best')
plt.tight_layout()
plt.savefig('C:/Users/usuario/Desktop/Facultad/8vo Semestre/Fisica Comp/PROYECTO FINAL/L3/rk4_l3_zoom.pdf', bbox_inches='tight')

#%%-------------------Euler-----------------------------------------------------------------------------
plt.figure()
plt.plot(xp_euler,yp_euler, color = 'green', linewidth=0.5 )
plt.plot(rsj[0], rsj[1], 'o', color='orange', label=r'Planeta ($30M_J$)', markersize = 10)
plt.plot(rs[0], rs[1], 'o', color='yellow', markeredgecolor="#FD7813", label='Sol', markersize = 15)
plt.plot(l3[0], l3[1], '+', color = 'black', label = r'$\mathcal{L}_3$')#Cambiar aqui segun que condicion inicial se este usando
plt.xlabel(r'$X$')
plt.ylabel(r'$Y$')
plt.legend(loc = 'best')
plt.tight_layout()
plt.savefig('C:/Users/usuario/Desktop/Facultad/8vo Semestre/Fisica Comp/PROYECTO FINAL/L3/euler_l3.pdf', bbox_inches='tight')

#%%-------------------Verlet----------------------------------------------------------------------------
plt.figure()
plt.plot(xp_v, yp_v, color='purple' ,linewidth=0.5 )
plt.plot(rsj[0], rsj[1], 'o', color='orange', label=r'Planeta ($30M_J$)', markersize = 10)
plt.plot(rs[0], rs[1], 'o', color='yellow', markeredgecolor="#FD7813", label='Sol', markersize = 15)
plt.plot(l5[0], l5[1], '+', color = 'black', label = r'$\mathcal{L}_5$')#Cambiar aqui segun que condicion inicial se este usando
plt.xlabel(r'$X$')
plt.ylabel(r'$Y$')
plt.legend(loc = 'best')
plt.tight_layout()
plt.savefig('C:/Users/usuario/Desktop/Facultad/8vo Semestre/Fisica Comp/PROYECTO FINAL/L5/verlet_l5.pdf', bbox_inches='tight')

#%%Animacion

fig, ax = plt.subplots()
posx = []
posy = []

# Plot los puntos fijos fuera de la función animate
ax.plot(rsj[0], rsj[1], 'o', color='orange', label=r'Planeta ($30M_J$)', markersize = 10)
ax.plot(rs[0], rs[1], 'o', color='yellow', markeredgecolor="#FD7813", label='Sol', markersize = 15)
ax.plot(l4[0], l4[1], '+', color = 'black', label = r'$\mathcal{L}_4$')#Cambiar aqui segun que condicion inicial se este usando
ax.set_xlabel(r'$X$')
ax.set_ylabel(r'$Y$')
ax.legend(loc='best')
ax.set_title(r'Trayectoria de una Partícula lanzada desde $\mathcal{L}_4$')#Cambiar aqui segun que condicion inicial se este usando
# Inicializar la línea de trayectoria (vacía)
line, = ax.plot([], [], color = 'purple' ,label='Trayectoria de la partícula')

def animate(i):
    index = 2000 * i  # Aumenta el índice en 2 en cada paso, para mayor resolucion bajamos el index, mas alto el valor, mas rapida la animacion
    if index < len(xp_v):
        posx.append(xp_v[index])
        posy.append(yp_v[index])

        # Actualizar los datos de la línea de trayectoria
        line.set_data(posx, posy)
        ax.relim()  # Recalcular los límites del gráfico
        ax.autoscale_view()  # Autoscalar el gráfico
        ax.legend(loc='best')

ani = FuncAnimation(fig, animate, frames=len(xp_v)//2, interval=1, repeat=False)

plt.show()

