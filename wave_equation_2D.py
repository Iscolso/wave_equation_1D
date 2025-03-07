import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation

def main():
    class Wave2D:
        def __init__ (self, Lx=10,Ly=10,Nx=100, Ny=100, c=1, dt=0.01, Nt=200, damping_width=10):
            self.Lx = Lx
            self.Ly = Ly
            self.Nx = Nx
            self.Ny = Ny
            self.c = c
            self.dt = dt
            self.dx = Lx/Nx
            self.dy = Ly/Ny
            self.Nt = Nt

            #Condicion CFL

            self.sx = self.c * self.dt/self.dx
            self.sy = self.c * self.dt/self.dy

            if self.sx > 1 or self.sy > 1:
                raise ValueError('No se cumple la condicion CFL')
            else:
                print('Estabilidad garantizada')

            self.u = np.zeros((Nx,Ny))
            self.u_prev = np.zeros((Nx,Ny))
            self.u_next = np.zeros((Nx,Ny))

            #Condicion inicial 
            x = np.linspace(0,Lx,Nx)
            y = np.linspace(0,Ly,Ny)
            X,Y = np.meshgrid(x,y)
            x0, y0 = Lx/2, Ly/2
            sigma = Lx/10
            self.u = np.exp(-((X-x0)**2 + (Y-y0)**2)/(2*sigma**2))

            #Iniciacion correcta de u

            self.u_prev[:,:] = self.u[:,:] 

            #barrera absorbente 

            self.damping_zone = np.ones((Nx,Ny))
            d = damping_width

            for i in range(Nx):
                for j in range(Ny):
                    if i > Nx - d or j > Ny - d or i < d or j < d:
                        self.damping_zone[i, j] = np.exp(-((min(i, j, Nx-i, Ny-j)) / d) ** 2)

            # Configuración de la figura
            self.fig, self.ax = plt.subplots()
            self.im = self.ax.imshow(self.u, cmap='plasma', vmin=-1, vmax=1, extent=[0, Lx, 0, Ly])
            self.ax.set_xlabel("X")
            self.ax.set_ylabel("Y")

        def equation(self, frame):
            """Actualiza la onda usando diferencias finitas en 2D"""
            for i in range(1, self.Nx - 1):
                for j in range(1, self.Ny - 1):
                    self.u_next[i, j] = (2 * self.u[i, j] - self.u_prev[i, j] + self.sx**2 * (self.u[i+1, j] + self.u[i-1, j] - 2 * self.u[i, j]) + self.sy**2 * (self.u[i, j+1] + self.u[i, j-1] - 2 * self.u[i, j]))

             # Aplicar la barrera absorbente en los bordes
            self.u_next *= self.damping_zone  

            # Actualizar estados
            self.u_prev[:, :] = self.u[:, :]
            self.u[:, :] = self.u_next[:, :]

            # Actualizar la animación
            self.im.set_array(self.u)
            return self.im,

        def simulation(self, speed=50):
            ani = animation.FuncAnimation(self.fig, self.equation, frames=self.Nt, interval=speed, blit=True)
            plt.show()

# Ejecutar la simulación
    wave_simulation = Wave2D()
    wave_simulation.simulation(speed=50)



if __name__ == '__main__':
    main()