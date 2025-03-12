import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation 
from matplotlib.animation import PillowWriter 

def main():
    class Wave:
        def __init__(self, L=10, Nx=100, c=1, dt=0.01, Nt=200, damping_width=20):
            self.L = L  
            self.Nx = Nx  
            self.c = c  
            self.dx = L / Nx  
            self.dt = dt  
            self.Nt = Nt  
            self.s = self.c * self.dt / self.dx  

            if self.s > 1:
                raise ValueError('No se cumple la condición CFL')  
            else:
                print('Estabilidad garantizada')

            # Condición inicial: Onda gaussiana
            self.x = np.linspace(0, L, Nx)  
            self.set_initial_conditions()

            #Funcion de profundidad del oceano

            self.depth = self.compute_depth(self.x)


            self.damping_zone = np.ones(Nx)
            d = damping_width
            for i in range(Nx-d,Nx):
                self.damping_zone[i] = np.exp(-3*((i - (Nx - d))/d)**2)

            # Gráfico
            self.fig, self.ax = plt.subplots()
            self.ax.set_xlim(0, L)
            self.ax.set_ylim(-1, 1)
            self.line, = self.ax.plot(self.x, self.u, lw=2)

        def set_initial_conditions(self):
            x0 = self.L/3
            sigma = self.L/10
            self.u = np.exp(-((self.x - x0)**2)/(2*sigma**2))
            self.u_prev = np.copy(self.u)
            self.u_prev[1:-1] = self.u[1:-1] - 0.5 * self.s**2 * (self.u[2:] - 2*self.u[1:-1] + self.u[:-2])
            self.u_next = np.copy(self.u)

        def compute_depth(self, x):
            max_depth = 100
            min_depth = 10
            coast_start = np.random.uniform(self.L *0.6, self.L * 0.8)
            depth = np.full_like(x,min_depth)

            deep_region = x < coast_start
            depth[deep_region] = max_depth

            #costa 
            slope_region = x >= coast_start
            depth[slope_region] = depth[slope_region] = min_depth + (max_depth - min_depth) * ((x[slope_region] - coast_start) / (self.L - coast_start))

            return depth

        def equation(self, frame):
            """ Actualiza la onda usando diferencias finitas """
            for i in range(1, self.Nx - 1):
                #Como la velocidad depende de la profundidad (c = sqrt(g*h))
                c_local = np.sqrt(9.81 * self.depth[i])
                s_local = min(c_local * self.dt/ self.dx, 1.0)

                #Verificacion de valores NaN, inf, etc.

                #if not np.isfinite(s_local):
                    #print(f' s_local invalido en i={i}: {s_local}')
                    #print(f' dt={self.dt}, dx={self.dx}, depth={self.depth[i]}')

                #print(f's_local={s_local}, CFL={s_local**2}')

                self.u_next[i] = 2 * self.u[i] - self.u_prev[i] + (s_local ** 2) * (self.u[i + 1] - 2 * self.u[i] + self.u[i - 1])
            # Barrera absorbente en la parte derecha 

            self.u_next *= self.damping_zone
            self.u_next = np.maximum(self.u_next,0)

            # Condición de frontera reflejante (rebote)
            #self.u_next[0] = 0  
            #self.u_next[-1] = 0  

            # Actualizar estados
            self.u_prev[:] = self.u  
            self.u[:] = self.u_next  

            if np.max(self.u) < 0.01:
                self.set_initial_conditions()

            # Actualizar la animación
            self.line.set_ydata(self.u)
            return self.line,

        def simulation(self, speed=20):
            ani = animation.FuncAnimation(self.fig, self.equation, frames=100, interval=15, blit=True)
            #ani.save('wave_simulation.gif', writer=PillowWriter(fps=30))
            plt.show()

    wave_simulation = Wave()
    wave_simulation.simulation(speed=20)

if __name__ == '__main__':
    main()
