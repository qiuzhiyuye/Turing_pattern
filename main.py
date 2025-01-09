import matplotlib.pyplot as plt
import numpy as np
import argparse


class TuringPattern(object):
    def __init__(self, sizex, sizey, dx, dt, Du, Dv, feed_rate, kill_rate, epoch, equation,name):
        self.sizex = sizex
        self.sizey = sizey
        self.dx = dx
        self.dt = dt
        self.Du = Du
        self.Dv = Dv
        self.feed_rate = feed_rate
        self.kill_rate = kill_rate
        self.epoch = epoch
        self.equation = equation
        self.name = name

        self.U = np.ones((sizex, sizey), dtype=float)
        self.V = np.zeros((sizex, sizey), dtype=float)
        self.V[50:60, 50:70] = 1
        self.V[60:80, 70:80] = 1

    # 拉普拉斯算子
    # 通过中心点周围8个点对其下一时刻值进行更新
    def laplacian(self, in_array):
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
        if(self.equation == "Gray-Scott"):
            return - self.U * self.V ** 2 + self.feed_rate * (1 - self.U)

    def g_function(self):
        if(self.equation == "Gray-Scott"):
            return self.U * self.V ** 2 - (self.kill_rate + self.feed_rate) * self.V

    def difference_method(self, in_array):
        # mu = self.dt / self.dx ** 2
        sum_direction_array = 0.25 * (np.roll(in_array, 1, axis=0)
                                      + np.roll(in_array, -1, axis=0)
                                      + np.roll(in_array, 1, axis=1)
                                      + np.roll(in_array, -1, axis=1))
        # out_array = -mu * (sum_direction_array - 4 * self.U) / 4
        out_array = sum_direction_array - in_array
        return out_array

    def run_epoch_difference(self):
        self.U += self.dt * (self.Du * self.difference_method(self.U) + self.f_function())
        self.V += self.dt * (self.Dv * self.difference_method(self.V) + self.g_function())

    # 有限差分迭代计算
    # def run_difference(self):
    #     for run_time in range(self.epoch):
    #         self.run_epoch_difference()

    def run_difference(self):
        for run_time in range(self.epoch):
            self.run_epoch_difference()
            if run_time % 100 == 0:
                self.show_process(run_time)  # 保存每轮迭代后的图像

    def show_process(self, run_time):
        plt.imshow(self.V, cmap="plasma", interpolation="nearest")
        plt.colorbar()
        s = f"{self.name}_final_TuringPattern_epoch_{run_time}.png"
        plt.title(s)
        # plt.show()
        plt.pause(0.01)
        plt.clf()

    def run_epoch(self):
        laplacian_u = self.laplacian(self.U)
        laplacian_v = self.laplacian(self.V)

        deltaU = self.Du * laplacian_u + self.f_function()
        deltaV = self.Dv * laplacian_v + self.g_function()

        self.U += self.dt * deltaU
        self.V += self.dt * deltaV

    def run(self):
        for run_time in range(self.epoch):
            self.run_epoch()

    def show(self):
        plt.imshow(self.V, cmap="plasma", interpolation="nearest")
        plt.colorbar()
        s = self.name+"final TuringPattern"
        plt.title(s)
        plt.show()

# dx = 1
# dt = 0.25
# Du = 1
# Dv = 0.5
# feed_rate = 0.055
# kill_rate = 0.062
# feed_rate = 0.039
# kill_rate = 0.058
# sizex = 128
# sizey = 128
# epoch = 100000
# T = TuringPattern(sizex, sizey, dx, dt, Du, Dv, feed_rate, kill_rate, epoch)

T = TuringPattern(128, 128, 1, 0.25, 1, 0.5, 0.039, 0.058, 100000,"Gray-Scott","leopard print")

T.run_difference()
# T.run()
T.show()
