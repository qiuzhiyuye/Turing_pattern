### Turing Pattern

原有代码的迭代计算过程为(通过中心点以及其周围八个点对$t=n+1$时刻数据值进行更新，其中对角线的影响更小)

$$
U_{i,j}^{n+1} = U_{i,j}^n + dt(Du*laplacian+f)
$$

$$
laplacian = -U_{i,j}^n+\frac{1}{5}(U_{i,j+1}^n+U_{i,j-1}^n+U_{i-1,j}^n+U_{i+1,j}^n)+\frac{1}{20}(U_{i+1,j+1}^n+U_{i-1,j-1}^n+U_{i-1,j+1}^n+U_{i+1,j-1}^n)
$$


采用中心差分进行替换

$$
U_{i,j}^{n+1} = U_{i,j}^n+dt*(Du*difference+f)
$$

$$
difference = -U_{i,j}^n+\frac{1}{4}(U_{i,j+1}^n+U_{i,j-1}^n+U_{i-1,j}^n+U_{i+1,j}^n)
$$
