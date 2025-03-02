import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation 

def main():
    class Wave:
        def __init__(self, L=10 , Nx=100,c=1,dt=0.01, Nt=200):
            self.L = L #Longitud del dominio
            self.Nx = Nx #Numero de pasos espaciales 
            self.c = c #constante de propagacion
            self.dx = L/Nx #Espacios de separacion
            self.dt = dt #espaciado 
            self.Nt = Nt # numero de pasos de tiempo 
            self.s = self.c * self.dt / self.dx
            if self.s > 1:
                raise ValueError('No se cumple la condicion CFL') 
            else:
                print('Estabilidad garantizada')

            self.u = np.zeros(Nx)

            #Condiciones iniciales con onda gaussiana 

            x = np.linspace(0,L, Nx) #puntos del dominio 
            x0 = L/3 #Inicio de la onda
            sigma = L/10  #Longitud de la onda 
            self.u = np.exp(-((x-x0)**2) / (2*sigma**2))

            
            self.u_prev = np.copy(self.u)  # Copia inicial
            #self.u_prev[1:-1] = self.u[1:-1] - 0.5 * self.s**2 * (self.u[2:] - 2*self.u[1:-1] + self.u[:-2])

            self.u_next = np.copy(self.u)

            #Grafico

            self.fig, self.ax = plt.subplots()
            self.ax.set_xlim(0, L)
            self.ax.set_ylim(-1, 1)
            self.line, = self.ax.plot(np.linspace(0, L, Nx), self.u)
            self.scatter = self.ax.scatter(x,self.u,c=self.u, cmap='plasma', vmin=-1, vmax=1)
        
        def equation(self, frame):
            #Empleamos la ecuacion de ondas

            for i in range(1, self.Nx-1):
                self.u_next[i] = ((2*self.u[i]) - self.u_prev[i] + ((self.s)**2) * ((self.u[i+1]) - (2*self.u[i]) + self.u[i-1]))
                
            #Condiciones de frontera

            #self.u_next[0] = self.u[0] * 0.5  # Se disipa la mitad
            self.u_next[-1] = self.u[-1] * 0.5  # Se disipa al final también


            self.u_prev[:] = self.u
            self.u[:] = self.u_next
            self.scatter.set_array(self.u)
            return self.scatter,
        
        def simulation(self, speed= 100):
            ani = animation.FuncAnimation(self.fig, self.equation, frames=self.Nt, interval=50)
            plt.show()

    wave_simulation = Wave()
    wave_simulation.simulation(speed=100)
    

if __name__ == '__main__':
    main()