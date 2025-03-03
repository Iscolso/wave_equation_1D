import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation 

def main():
    class Wave:
        def __init__(self, L=10, Nx=100, c=1, dt=0.01, Nt=200):
            self.L = L  
            self.Nx = Nx  
            self.c = c  
            self.dx = L / Nx  
            self.dt = dt  
            self.Nt = Nt  
            self.s = self.c * self.dt / self.dx  

            if self.s > 1:
                raise ValueError('❌ No se cumple la condición CFL')  
            else:
                print('✅ Estabilidad garantizada')

            # Condición inicial: Onda gaussiana
            x = np.linspace(0, L, Nx)  
            x0 = L / 3  
            sigma = L / 10  
            self.u = np.exp(-((x - x0) ** 2) / (2 * sigma ** 2))  

            # Inicialización correcta de `u_prev` para que la onda se mueva
            self.u_prev = np.copy(self.u)
            self.u_prev[1:-1] = self.u[1:-1] - 0.5 * self.s**2 * (self.u[2:] - 2*self.u[1:-1] + self.u[:-2])

            self.u_next = np.copy(self.u)

            # Gráfico
            self.fig, self.ax = plt.subplots()
            self.ax.set_xlim(0, L)
            self.ax.set_ylim(-1, 1)
            self.line, = self.ax.plot(x, self.u, lw=2)

        def equation(self, frame):
            """ Actualiza la onda usando diferencias finitas """
            for i in range(1, self.Nx - 1):
                self.u_next[i] = 2 * self.u[i] - self.u_prev[i] + (self.s ** 2) * (self.u[i + 1] - 2 * self.u[i] + self.u[i - 1])

            # Condición de frontera reflejante (rebote)
            self.u_next[0] = 0  
            self.u_next[-1] = 0  

            # Actualizar estados
            self.u_prev[:] = self.u  
            self.u[:] = self.u_next  

            # Actualizar la animación
            self.line.set_ydata(self.u)
            return self.line,

        def simulation(self, speed=50):
            ani = animation.FuncAnimation(self.fig, self.equation, frames=self.Nt, interval=speed, blit=True)
            plt.show()

    wave_simulation = Wave()
    wave_simulation.simulation(speed=50)

if __name__ == '__main__':
    main()
