# 计算流体力学 期末大作业

## 1. 问题描述

题目要求对一个无粘可压流体流过一个折角的流动进行数值模拟. 我们预期的结果是一个标准的斜激波. 我们需要计算的是流场的速度、压力、密度等参数的分布. 

题给条件为:   入流边界速度为$u=686.47 \mathrm{m/s}$, $v=0$, 压强为$99719 \mathrm{Pa}$, 温度为$293.15 \mathrm{K}$, $\gamma = 1.4$, 折角大小$\theta=15^\circ$. 无量纲化后的流动参数为$u=2.0$, $v = 0$, $p = 0.714497$, $\rho = 1.0$.

斜激波有理论解, 我们可以通过计算结果与理论解的比较来验证我们的数值模拟的正确性. 对于题给条件, 我们可以计算得到激波倾角约为$45^\circ$, 激波前后压强比约为$2.195$, 激波后的马赫数约为$1.446$

### 1.1 控制方程

无粘可压流体对应的是Euler方程, 即

$$
\frac{\partial U}{\partial t} + \frac{\partial F}{\partial x} + \frac{\partial G}{\partial y} = 0
$$

其中

$$
U = \begin{bmatrix} \rho \\ \rho u \\ \rho v \\ \rho E \end{bmatrix}, \quad 
F = \begin{bmatrix} \rho u \\ \rho u^2 + p \\ \rho uv \\ \rho u H \end{bmatrix}, \quad 
G = \begin{bmatrix} \rho v \\ \rho uv \\ \rho v^2 + p \\ \rho v H \end{bmatrix}, 
$$

这里, $\rho$ 是密度, $u, v$ 是速度, $p$ 是压强, $E = \dfrac{p}{\gamma - 1} + \dfrac{1}{2} \rho (u^2 + v^2)$.

### 1.2 边界条件

根据帮助文档给出的边界条件, 左边界是入流边界, 右边界和上边界是出流边界, 下边界比较难以处理. 

起初我们尝试给定下边界速度为0, 密度和压强为入流边界的值, 但是这样会导致下边界的流动速度不符合物理规律, 基本无法收敛. 我们改为使用反射边界条件, 即下边界的法向速度为0, 不限制切向速度. 

### 1.3 网格划分

在帮助文档给出的网格划分的基础上, 我们对x和y方向的网格数进行了调整, 此外, 我们调整了y方向网格的划分方式. 划分后的网格如图. 

![[Result/Mesh.png|480]]

### 1.4 离散格式

空间上我们使用有限体积法离散控制方程, 原方程在单个控制体上的半离散格式为

$$
\frac{dU}{dt} = -\sum_i \frac{\left|l_i\right|}{\left|\Omega\right|} \vec{F} \cdot \vec{n_i} \ ,
$$

其中$\left|l_i\right|$是边界的长度, $\left|\Omega\right|$是控制体的面积, $\vec{n_i}$是边界的外法向, $\vec{F} = [F, G]$是界面上的通量.

时间上我们使用向前Euler法进行时间离散, 虽然简单但在我们的问题上表现足够好. 

蔡庆东老师讲义上比较高级的说法是，这种求解方法的计算过程包含三个步骤:  

1. 重构 Reconstruction:  利用$t_n$时刻的网格均值构造分片常数的空间分布; 
2. 演化 Evolution:  在网格交界面处，利用Riemann 问题精确解给出$t_{n+1} - t_{n}$时间间隔内的通量, 然后将所有的Riemann解综合起来, 给出$t_{n+1}$时刻空间上的整体解$U(x, t{n+1})$；
3. 投影 Projection:  将$t_{n+1}$时刻的空间整体解在每个网格内做平均，给出$t_{n+1}$时刻网格上的平均值. 在实际计算中, 由于Riemann解的复杂性，我们并不是真的去对函数$U(x, t_{n+1})$ 做平均, 而是利用演化步的Riemann解计算出界面处的数值通量, 然后直接计算网格内的均值. 

### 1.5 重构与激波捕捉

在有限体积法中, 我们存储的是每个控制体的平均值, 为了计算界面上的对流通量$\vec{F}$, 我们需要对这些平均值进行重构. 界面上通量的计算总体上有中心型格式和迎风型格式两大类, 在有激波的情况下, 为了避免数值震荡, 我们通常使用 (某种意义上的) 迎风型格式, 有时需要在间断处降低阶数. 这里我们都采用NND格式. 对物理量$\phi$的重构公式为

$$
\phi_{i+1/2} = \phi_i + \frac{1}{2} \mathrm{minmod}(\phi_{i+1}-\phi_i,\phi_i-\phi_{i-1}), 
$$

其中

$$
\mathrm{minmod}(a, b) = \left\{
\begin{aligned}
&\min(a, b), \quad \text{if } a > 0 \text{ and } b > 0, \\
&\max(a, b), \quad \text{if } a < 0 \text{ and } b < 0, \\
&0, \quad \text{otherwise}.
\end{aligned}
\right.
$$

在下面计算数值通量时, 我们会根据需要对流场的原始变量或是特征变量进行重构. 

### 1.6 数值通量

重构得到界面上的物理量后, 我们需要计算界面上的通量. 单元界面上的通量计算常常可以被看作Riemann问题. 通量计算有FVS和FDS两大类方法, 我们这里分别采用FVS中的Roe格式和FDS中的Van Leer格式. 特别地, 我们还尝试了KFVS格式.

#### 1.6.1 Roe格式

Roe 格式的基本理念是使用界面左右两侧的状态来构造一个平均值(Roe平均), 进而构造矩阵$\tilde A$替代原始的Jacobi矩阵$A = \partial F / \partial U$, 将非线性问题转化为线性问题. 

Roe 平均涉及使用密度加权平均从状态 $U_0$ 和 $U_r$ 计算平均状态 $\bar{U}$:   
$$
\begin{aligned}
\bar{\rho} &= \sqrt{\rho_0 \rho_r}, \\ 
\bar{u} &= \frac{\sqrt{\rho_0} u_0 + \sqrt{\rho_r} u_r}{\sqrt{\rho_0} + \sqrt{\rho_r}},\\
\bar{h} &= \frac{\sqrt\rho_0 H_0 + \sqrt\rho_r H_r}{\sqrt\rho_0 + \sqrt\rho_r}.
\end{aligned}
$$
然后使用这些平均构造通量的Jacobi矩阵$\tilde A$

Roe格式的通量计算公式为

$$
\tilde{F}(U_0, U_r) = \frac{1}{2} [ F(U_0) + F(U_r) - |\tilde{A}|(U_0, U_r) (U_r - U_0) ],
$$
其中, 
$$
A = L^{-1} \Lambda L, \quad
\left| \tilde A\right| = L^{-1} \left|\Lambda\right| L, 
$$
$\Lambda$是对角的特征值矩阵. 

#### Van Leer 格式

Van Leer 通量分裂方法基于特征分解，通过计算界面处的法向马赫数 $M_n = \frac{V}{c}$，其中 $V$ 是界面法向速度，$c$ 是声速，来处理流动界面上的通量问题。此方法通过对流动的原始变量进行插值到界面处，以形成左右两个状态，从而实现通量的计算。

**马赫数的计算与分裂**

首先计算左右侧的法向马赫数 $M_0$ 和 $M_r$，然后根据以下公式计算正部分 $M_0^+$ 和负部分 $M_r^-$: 

- **$M_0^+$ 定义为**: 
$$
  M_0^+ = 
  \begin{cases} 
  M_0 & \text{if } M_0 \geq 1, \\
  \frac{1}{4} (M_0 + 1)^2 & \text{if } |M_0| < 1, \\
  0 & \text{if } M_0 \leq -1.
  \end{cases}
$$

- **$M_r^-$ 定义为**: 
$$
  M_r^- = 
  \begin{cases} 
  0 & \text{if } M_r \geq 1, \\
  -\frac{1}{4} (M_r - 1)^2 & \text{if } |M_r| < 1, \\
  M_r & \text{if } M_r \leq -1.
  \end{cases}
$$

最终的界面马赫数由 $M_n = M_0^+ + M_r^-$ 给出。

**超声速情况 ($|M_n| \geq 1$)**

- 如果 $M_n \geq +1$，则界面处只有正通量 $F_0^+ = F_0$，负通量 $F_0^- = 0$。
- 如果 $M_n \leq -1$，则界面处只有负通量 $F_0^+ = 0$，正通量 $F_0^- = F_r$。

**亚声速情况 ($|M_n| < 1$)**

界面正负通量由以下公式给出: 
$$
F_0^\pm = \left\{
\begin{array}{l}
f_\text{mass}^\pm \\
f_\text{mass}^\pm \left[ n_x \left(-V \pm 2c\right)/\gamma + u \right] \\
f_\text{mass}^\pm \left[ n_y \left(-V \pm 2c\right)/\gamma + v \right] \\
f_\text{energy}^\pm
\end{array}
\right.
$$

$$
F_0 = F_0^+ + F_0^-
$$

其中的质量通量和能量通量为: 
$$
f_\text{mass}^+ = \rho_0 c_0 \frac{(M_0 + 1)^2}{4}, \quad f_\text{mass}^- = -\rho_r c_r \frac{(M_r - 1)^2}{4},
$$
$$
f_\text{energy}^\pm = f_\text{mass}^\pm \left[ \frac{(\gamma - 1)V \pm 2c)^2}{2(\gamma^2 - 1)} + \frac{u^2 + v^2 + w^2 - V2}{2} \right],
$$

$V$是法向速度. 

#### 1.6.3 边界处理

在边界上, 二阶重构会导致流场难以收敛, 造成如下图所示的数值振荡. 

![[Result/Log/HighOrder.png|480]]

我们在边界上采用一阶特征重构的Roe格式. 

#### 1.6.4 KFVS格式

KFVS(动理学矢通量分裂, kinetic flux vector splitting)格式是一种基于无碰撞的Boltzmann方程的矢通量分裂格式, 物理背景比较清晰. 在宏观尺度上, KFVS会过渡到Euler方程. 对于可压缩流体, KFVS格式表现良好. 在形式上, KFVS格式也是一种FVS格式, 且具有一定的迎风性质. 

这种格式的原理是先由宏观参量得到界面两侧气体粒子的分布函数, 再由分布函数计算单元界面上守恒量的通量, 而这和我们前述两种Riemann Solver做的事情本质上是一样的. 于是格式可以被直接嵌入到一般的有限体积法程序中. 

具体推导过程不再赘述, 直接给出x方向上数值通量的表达式: 

$$
\begin{align*} 
F_1^\pm &= \rho u X^\pm \pm p Y, \\ 
F_2^\pm &= (p + \rho u^2) X^\pm \pm \rho u Y, \\ 
F_3^\pm &= \rho u v X^\pm \pm \rho v Y, \\ 
F_4^\pm &= \left( \frac{\gamma}{\gamma - 1} p + \frac{1}{2} \rho (u^2 + v^2) \right) u X^\pm \\
&\pm \left( \frac{\gamma + 1}{2(\gamma - 1)} p + \frac{1}{2} \rho (u^2 + v^2) \right) Y,
\end{align*}
$$
$$
F = \begin{bmatrix}
                  F_1^+ + F_1^- \\
                  F_2^+ + F_2^- \\
                  F_3^+ + F_3^- \\
                  F_4^+ + F_4^- \\
\end{bmatrix}
$$

其中

$$
X = \frac{1 + \text{erf}(S)}{2}, \quad Y_i = \frac{\exp(-S^2)}{\sqrt{\pi \beta}}, \quad S_i = u_i \sqrt{\beta}, \quad \beta = \frac{1}{2RT}
$$

具体实现过程中, 我们对流场的原始变量使用NND方法重构, 得到边界值, 利用KFVS计算通量. 

气体动理学格式(GKS格式)是由香港科技大学的徐昆等人提出的一种新的数值通量计算格式, 格式基于对Boltzmann方程中的碰撞项进行BGK近似. KFVS可以看作是GKS的一种退化情形.

参考文献:   

Mandal, J. C., & Deshpande, S. M. (1994). Kinetic flux vector splitting for euler equations. Computers & Fluids, 23(2), 447–478. doi:  10.1016/0045-7930(94)90050-7 

李启兵, 徐昆.  气体动理学格式研究进展[J]. 力学进展, 2012, 42(5):   522-537. doi:   10.6052/1000-0992-11-149

## 2. 计算结果

### 2.1 Roe格式

无量纲压强:  

![[Result/Roe/Pressure.png|480]]

马赫数:  

![[Result/Roe/Mach.png|480]]

压强比和马赫数:  

```bash
Pressure Ratio before and after shockwave:    2.205398473437763
Mach Number after shockwave:    1.4334660698992898
```

### 2.2 Van Leer格式

无量纲压强:  

![[Result/V-L/Pressure.png|480]]

马赫数:  

![[Result/V-L/Mach.png|480]]

压强比和马赫数:  

```bash
Pressure Ratio before and after shockwave:    2.1647915516432534
Mach Number after shockwave:    1.4502738096147407
```

### 2.3 KFVS格式

无量纲压强:  

![[Result/KFVS/Pressure.png|480]]

马赫数:  

![[Result/KFVS/Mach.png|480]]

压强比和马赫数:  

```bash
Pressure Ratio before and after shockwave:   2.175568442117899
Mach Number after shockwave:   1.4480725961547125
```

### 2.4 结果分析

可以发现, 两种格式的计算结果都与斜激波理论值(2.195, 1.446)较为接近, 且比较好地捕捉到了激波的位置. 同时, 两种格式没有给出明显的数值震荡, 也没有引入过多的数值耗散. 总的来说, 我们的数值模拮是成功的.

其中数值误差的来源除了来自计算过程中的近似, 还与选取的激波前后数据点的位置有关. 由于激波的位置是通过观察计算结果得到的, 所以适当调整激波前后数据点的位置可能会得到更好的结果.

图片上来看, KFVS计算结果和理论值更加接近, 激波更加尖锐, 效果更好, 但这有待我们在更多问题中进一步验证. 

## 3. 结论与展望

### 3.1 结论

本次大作业中, 我们成功地对一个无粘可压流体流过一个折角的流动进行了数值实验, 并得到了较为满意的结果. 我们的数值模拮是成功的, 但仍有一些改进的空间.

我们原本打算使用有限差分方法, 但是在给出的网格不正交, 有限差分方法, 特别是基于通量分裂的有限差分方法, 在处理通量上面有困难:  在非正交坐标系下, 通量Jacobi矩阵的形式不确定, 通量分裂的方法不再适用. 我们转而使用有限体积法. 

起初, 我们给定的边界条件不符合物理规律, 计算难以收敛. 

![[Result/Log/Legacy/Wrong BC.png|480]]

更改边界条件后有改善, 但是仍然会发生数值振荡. 

![[Result/Log/HighOrder.png|480]]

发生震荡的位置可以看成一个准一维问题, 我们尝试对其建模, 采用了不同的数值通量处理方式, 但是效果仍然不理想. 

我们尝试了使用Julia语言, 用有限体积法对一维激波在边界处的反射问题进行了数值模拟. 

```julia
using Plots

function main(max_iter: : Int)
    n = 20

    x = collect(0: 1/n: 1)
    p = ones(n+4-1) * 0.714497
    u = ones(n+4-1) * 0.9659258262890683
    v = ones(n+4-1) * 0.25881904510252074
    rho = ones(n+4-1)

    p_new = ones(n+4-1) * 0.714497
    u_new = ones(n+4-1) * 0.9659258262890683
    v_new = ones(n+4-1) * 0.25881904510252074
    rho_new = ones(n+4-1)


    gamma = 1.4

    dt = 0.005


    for cnt = 1: max_iter
        u[1] = u[3]
        u[2] = u[3]
        v[1] = v[3]
        v[2] = v[3]
        rho[1] = rho[3]
        rho[2] = rho[3]
        p[1] = p[3]
        p[2] = p[3]

        u[end] = u[end-2]
        u[end-1] = u[end-2]

        v[end] = -v[end-2]
        v[end-1] = -v[end-2]

        rho[end] = rho[end-2]
        rho[end-1] = rho[end-2]
        p[end] = p[end-2]
        p[end-1] = p[end-2]

        for i = 3: n+1
            
            E = p[i] / (gamma - 1) + 0.5 * rho[i] * (u[i]^2 + v[i]^2)
            E_up = p[i+1] / (gamma - 1) + 0.5 * rho[i+1] * (u[i+1]^2 + v[i+1]^2)
            E_down = p[i-1] / (gamma - 1) + 0.5 * rho[i-1] * (u[i-1]^2 + v[i-1]^2)
            
            h = (E + p[i]) / rho[i]
            h_up = (E_up + p[i+1]) / rho[i+1]
            h_down = (E_down + p[i-1]) / rho[i-1]

            # println(u)

            u_up = u[i+1]
            u_down = u[i-1]
            v_up = v[i+1]
            v_down = v[i-1]
            rho_up = rho[i+1]
            rho_down = rho[i-1]
            p_up = p[i+1]
            p_down = p[i-1]

            # ============================
            # 我们在这里重构并处理数值通量
            # ============================

            # Apply Finite Volume Method
            rho_new[i] = rho[i] - (Flux_up[1] - Flux_down[1]) * dt / (1.0/n)
            u_new[i] = u[i] - (Flux_up[2] - Flux_down[2]) * dt / (1.0/n) / rho[i]
            v_new[i] = v[i] - (Flux_up[3] - Flux_down[3]) * dt / (1.0/n) / rho[i]
            E = E - (Flux_up[4] - Flux_down[4]) * dt / (1.0/n)
            p_new[i] = (gamma - 1) * (E - 0.5 * rho[i] * (u[i]^2 + v[i]^2))
            # if p_new[i] < 0
            #     p_new[i] = 0
            # end

        end
        rho = copy(rho_new)
        u = copy(u_new)
        v = copy(v_new)
        p = copy(p_new)

    end

    # visualize p
    plot(x[2: end-1], p[3: end-2], label="p", xlabel="x", ylabel="p", title="1D Shockwave", legend=: topleft)

end

main(30)
```

![[Result/Log/1D.png|480]]

最后我们试着改变边界条件, 在边界上给定速度, 密度和压强的值, 数值振荡的现象得到了很好的抑制. 

定位了问题在于边界条件上后, 我们回到二维代码, 边界条件形式不能修改, 但是我们尝试在边界上降阶, 问题终于得到了缓解. 

### 3.2 展望

我们可以尝试使用更高阶的重构格式, 如MUSCL格式, WENO格式等, 以提高数值精度, 进一步减小数值耗散.

我们可以尝试使用其他的数值通量处理方法, 如GKS格式, 以更好捕捉物理现象.

此外, 我们可以尝试使用更复杂的网格划分, 如非结构网格, 在激波附近使用更细的网格, 以提高数值精度. 

