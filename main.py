import matplotlib.pyplot as plt
import numpy as np
import argparse
import matplotlib.animation as animation


class TuringPattern(object):
    def __init__(self, sizex, sizey, dx, dt, Du, Dv, feed_rate, kill_rate, epoch, equation, name, process):
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
        self.process = process

        self.U = np.ones((sizex, sizey), dtype=float)
        self.V = np.zeros((sizex, sizey), dtype=float)
        # self.V[50:60, 50:70] = 1
        # self.V[60:80, 70:80] = 1
        self.Initial()

    def Initial(self):
        if (self.name == 'zebra'):
            # center_x, center_y = self.sizex // 2, self.sizey // 2
            center_x, center_y = 0, 0
            x, y = np.ogrid[:self.sizex, :self.sizey]
            radius = self.sizex // 5
            mask = (x - center_x) ** 2 + (y - center_y) ** 2 <= radius ** 2
            self.V[mask] = 1
            radius = radius - 5
            mask = (x - center_x) ** 2 + (y - center_y) ** 2 <= radius ** 2
            self.V[mask] = 0
            radius = radius - 5
            mask = (x - center_x) ** 2 + (y - center_y) ** 2 <= radius ** 2
            self.V[mask] = 1
            radius = radius - 5
            mask = (x - center_x) ** 2 + (y - center_y) ** 2 <= radius ** 2
            self.V[mask] = 0
            radius = radius - 5
            mask = (x - center_x) ** 2 + (y - center_y) ** 2 <= radius ** 2
            self.V[mask] = 1
            radius = radius - 5
            mask = (x - center_x) ** 2 + (y - center_y) ** 2 <= radius ** 2
            self.V[mask] = 0
            radius = radius - 5
            mask = (x - center_x) ** 2 + (y - center_y) ** 2 <= radius ** 2
            self.V[mask] = 1
        else:
            lenth = self.sizex // 10
            self.V[self.sizex // 2 - lenth:self.sizex // 2, self.sizex // 2 - lenth:self.sizex // 2 + lenth] = 1
            self.V[self.sizex // 2:self.sizex // 2 + lenth * 2, self.sizex // 2 + lenth:self.sizex // 2 + lenth * 2] = 1

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
        if self.equation == "Gray-Scott":
            return - self.U * self.V ** 2 + self.feed_rate * (1 - self.U)

    def g_function(self):
        if self.equation == "Gray-Scott":
            return self.U * self.V ** 2 - (self.kill_rate + self.feed_rate) * self.V

    def difference_method(self, in_array):
        temp_array = np.zeros((self.sizex, self.sizey), dtype=float)
        for i in range(in_array.shape[0]):
            if i == 0:
                # 外向递推 Un1 = 2Un-Un_1得到-1
                # temp_array[i, :] += (in_array[i, :] * 2 - in_array[i + 1, :]) + in_array[i + 1, :]
                temp_array[i, :] += 2*in_array[i, :]
            elif i == in_array.shape[0] - 1:
                # temp_array[i, :] += (in_array[i, :] * 2 - in_array[i - 1, :]) + in_array[i - 1, :]
                temp_array[i, :] += 2*in_array[i, :]
            else:
                temp_array[i, :] += (in_array[i - 1, :] + in_array[i + 1, :])
        for i in range(in_array.shape[1]):
            if i == 0:
                # 外向递推 Un1 = 2Un-Un_1得到-1
                # temp_array[:, i] += (in_array[:, i] * 2 - in_array[:, i + 1]) + in_array[:, i + 1]
                temp_array[:, i] += 2*in_array[:, i]
            elif i == in_array.shape[1] - 1:
                # temp_array[:, i] += (in_array[:, i] * 2 - in_array[:, i - 1]) + in_array[:, i - 1]
                temp_array[:, i] += 2*in_array[:, i]
            else:
                temp_array[:, i] += (in_array[:, i - 1] + in_array[:, i + 1])
        sum_direction_array = 0.25 * (np.roll(in_array, 1, axis=0)
                                      + np.roll(in_array, -1, axis=0)
                                      + np.roll(in_array, 1, axis=1)
                                      + np.roll(in_array, -1, axis=1))
        out_array = temp_array*0.25 - in_array
        return out_array

    def run_epoch_difference(self):
        self.U += self.dt * (self.Du * self.difference_method(self.U) + self.f_function())
        self.V += self.dt * (self.Dv * self.difference_method(self.V) + self.g_function())

    # 有限差分迭代计算
    def run_difference(self):
        # for run_time in range(self.epoch):
        #     self.run_epoch_difference()
        #     if self.process and (run_time % 100 == 0):
        #         self.show_process(run_time)  # 保存每轮迭代后的图像
        # break
        ani = self.create_animation()

    def show_process(self, run_time):
        plt.imshow(self.V, cmap="plasma", interpolation="nearest")
        plt.colorbar()
        s = f"{self.name}_final_TuringPattern_epoch_{run_time}.png"
        plt.title(s)
        # plt.show()
        plt.pause(0.001)
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

    def show_result(self):
        plt.imshow(self.V, cmap="plasma", interpolation="nearest")
        plt.colorbar()
        s = self.name + " final TuringPattern"
        plt.title(s)
        plt.show()

    def animate(self, i):
        for step in range(200):
            self.run_epoch_difference()
        self.im.set_array(self.V)
        return [self.im]

    def create_animation(self):
        fig, ax = plt.subplots()
        self.im = ax.imshow(self.V, cmap="plasma", interpolation="nearest")
        plt.colorbar(self.im)
        s = self.name + " TuringPattern Animation"
        plt.title(s)
        ani = animation.FuncAnimation(fig, self.animate, frames=self.epoch, interval=0.01, blit=True)
        plt.show()
        return ani


T = TuringPattern(128, 128, 1, 0.25, 1, 0.5, 0.039, 0.058, 30000, "Gray-Scott", "zebra", True)

# T = TuringPattern(256, 256, 1, 0.25, 1, 0.5, 0.039, 0.058, 30000, "Gray-Scott", "leopard print",True)

# T = TuringPattern(256, 256, 1, 0.25, 1, 0.5, 0.026, 0.061, 30000, "Gray-Scott", "leopard print",True)

# T = TuringPattern(256, 256, 1, 0.25, 1, 0.5, 0.0343, 0.0618, 50000, "Gray-Scott", "苏眉鱼纹路1",True)

# T = TuringPattern(128, 128, 1, 0.25, 1, 0.5, 0.0517, 0.0628, 50000, "Gray-Scott", "苏眉鱼纹路2",True)  #效果不错

# T = TuringPattern(256, 256, 1, 0.25, 1, 0.5, 0.098, 0.0555, 20000, "Gray-Scott", "巨蜥",True)

# T = TuringPattern(256, 256, 1, 0.25, 1, 0.5, 0.0517, 0.0628, 50000, "Gray-Scott", "苏眉鱼纹路2",True)  


# T.run()

T.run_difference()
# T.show_result()

# ani = T.create_animation()
