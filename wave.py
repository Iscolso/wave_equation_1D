import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
def main():
    class Wave():
        def __init__(self, dt=0.01, Nx=100, c=1, L=10, Nt=200, damping_width=50):
            self.dt = dt
            self.Nx = Nx
            self.c = c
            self.dx = L/Nx
            self.L = L
            self.Nt = Nt 
            self.s = self.c * self.dt / self.dx

            if self.s > 1:
                raise ValueError ('No se cumple la estabilidad.')
            else:
                print('Estabilidad garantizada!')

            #initial condition: gaussian wave

            x = np.linspace(0,L,Nx)
            x0 = L/3
            sigma = L/10
            self.u = np.exp(-((x-x0)**2)/(2*sigma**2))

            #Initialization 

            self.u_prev = np.copy(self.u)
            self.u_next = np.copy(self.u)

            #ocean depth

            self.depth = self.compute_depth(x)

            self.damping_zone = np.ones(Nx)
            d = damping_width
            for i in range(Nx-d, Nx):
                self.damping_zone[i] = np.exp(-3*((i - (Nx - d))/2)**2)

            #Graphic

            self.fig, self.ax = plt.subplots()
            self.ax.set_xlim(0,L)
            self.ax.set_ylim(-1,1)
            self.line = self.ax.plot(x,self.u, lw=2)

        

        def compute_depth(self, x):
            max_depth = 100
            min_depth = 10
            depth = max_depth - ( max_depth - min_depth) * np.sin(np.pi * x / self.L)
            return depth 

        def equation(self, frame):
            for i in range (1, self.Nx-1):
                c_local = np.sqrt(9.81*self.depth[i])
                s_local = min(c_local * self.dt/self.dx, 1.0)


                self.u_next[i] = 2 * self.u[i] - self.u_prev[i] + (s_local**2) * (self.u[i + 1] - 2 * self.u[i] + self.u[i - 1])

            self.u_next *= self.damping_zone
            self.u_next = np.maximum(self.u_next, 0)


            self.u_prev[:] = self.u
            self.u[:] = self.u_next

            self.line.set_ydata(self.u)

            return self.line,

        def simulation(self, speed=120):
            ani = animation.FuncAnimation(self.fig, self.equation, frames=100, interval= 15, blit=True)
            plt.show()

    wave_sim = Wave()
    wave_sim.simulation(speed=120)


if __name__ == '__main__':
    main()