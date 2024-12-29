import matplotlib.pyplot as plt
import numpy as np

class TuringPattern(object):
    def __init__(self,sizex,sizey,dx,dt,Du,Dv,feed_rate,kill_rate,epoch): 
        self.sizex = sizex
        self.sizey = sizey
        self.dx = dx
        self.dt = dt
        self.Du = Du
        self.Dv = Dv
        self.feed_rate = feed_rate
        self.kill_rate = kill_rate
        self.epoch = epoch

        self.U = np.ones((sizex, sizey), dtype=float)
        self.V = np.zeros((sizex, sizey), dtype=float)
        self.V[50:60, 50:70] = 1
        self.V[60:80, 70:80] = 1

    # 拉普拉斯算子
    def laplacian(self,in_array):
        center = -in_array
        direct_neighbors = 0.20 * (
            np.roll(in_array, 1, axis=0)
            + np.roll(in_array, -1, axis=0)
            + np.roll(in_array, 1, axis=1)
            + np.roll(in_array, -1, axis=1)
        )
        diagonal_neighbors = 0.05 * (
            np.roll(np.roll(in_array, 1, axis=0), 1, axis=1)
            + np.roll(np.roll(in_array, -1, axis=0), 1, axis=1)
            + np.roll(np.roll(in_array, -1, axis=0), -1, axis=1)
            + np.roll(np.roll(in_array, 1, axis=0), -1, axis=1)
        )

        out_array = center + direct_neighbors + diagonal_neighbors
        return out_array
    
    def f_function(self):
        return - self.U * self.V**2 + self.feed_rate * (1 - self.U)
    
    def g_function(self):
        return self.U * self.V**2 - (self.kill_rate + self.feed_rate) * self.V
    
    def run_epoch(self):
        laplacian_u = self.laplacian(self.U)
        laplacian_v = self.laplacian(self.V)

        deltaU = self.Du * laplacian_u + self.f_function()
        deltaV = self.Dv * laplacian_v + self.g_function()

        self.U += self.dt * deltaU
        self.V += self.dt * deltaV

    def run(self):
        for run_time in range(epoch):
            self.run_epoch()

    def show(self):
        plt.imshow(self.V, cmap="plasma", interpolation="nearest")
        plt.colorbar()
        s="final TuringPattern"
        plt.title(s)
        plt.show()




dx = 1
dt = 0.25  
Du = 1  
Dv = 0.5  
feed_rate = 0.055 
kill_rate = 0.062  
sizex = 128 
sizey = 128
epoch = 100000 
T=TuringPattern(sizex,sizey,dx,dt,Du,Dv,feed_rate,kill_rate,epoch)

T.run()
T.show()

























