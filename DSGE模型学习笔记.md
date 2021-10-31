## DSGE模型学习笔记

（目前还有几个地方没弄懂：式(6.33)到(6.34)的推导——`$1-\sigma$`是怎么出来的）

>**参考文献和学习资料**：
>
>连享会课件
>
>RBC模型：
>
>King, R. G., M. W. Watson, 1998, The solution of singular linear difference systems under rational expectations, International Economic Review, 39 (4): 1015-1026.
>
>Blanchard, O. J., C. M. Kahn, 1980, The solution of linear difference models under rational expectations, Econometrica, 48 (5): 1305-1311.
>
>McCandless, G.(2008), “The ABCs of RBCs”, Harvard University Press.
>
>Dejong, D. and Dave, C. (2007), “Structural Macroeconometrics”, Princeton University Press.
>
>NK模型：
>
>The Basic New Keynesian Model, Drago Bergholt.
>
>Christiano, Eichenbaum, and Evans(2005), Nominal Rigidities and the Dynamic Effects of a Shock to Monetary Policy, Journal of Politic Economics.
>
>An Estimated Dynamic Stochastic General Equilibrium Mmdel of The Euro Area（贝叶斯估计）
>
>A Baseline DSGE Model, Fernández-Villaverde and Rubio-Ramírez.（中等规模DSGE模型）

### ==0. 滤波方法==

* DSGE模型的研究对象是经济波动。因此滤波的目的是将经济变量分解为趋势项和周期项（假设两者不相关）。

* 常用的方法：

1. **对时间回归获得残差**，用于确定趋势的去除(detrend)

2. **自回归获得残差**

3. **HP滤波**，基于傅里叶变换，需要选择平滑系数`$\lambda$`

4. **BP滤波**，本质是多点滑动平均，使用过去和未来的数据进行平滑

### ==1. RBC模型（不考虑货币）==

* #### 不同的资源配置方式（决定了不同的考察角度）：

1. 中央计划者问题

2. 竞争均衡问题

* 根据福利经济学第二定理，两者等价。我们先从中央计划者的角度考察问题：

#### ==基本模型==

* #### 偏好

* 人们的偏好取决于消费`$C$`、闲暇`$L$`和折现因子`$\beta$`

```math
U=\sum_{t=0}^{\infty}{\beta}^t u(C_t,L_t) \tag{1.1}
```

* 其中偏好函数`$u()$`需满足一阶偏导大于零、凸性和**稻田条件**：

```math
\lim_{C \to 0}{u_C(C,L)=\infty}

\lim_{C \to \infty}{u_C(C,L)=0}

\lim_{L \to 0}{u_L(C,L)=\infty}

\lim_{L \to \infty}{u_L(C,L)=0}
```

* 折现因子`$\beta$`满足`$0<\beta<1$`

* #### 生产函数

```math
Y_t=F(K_t,N_t) \tag{1.2}
```

* 产出`$Y$`取决于资本`$K$`和劳动`$N$`。

* 生产函数的特点：**规模报酬不变**、资本存量是**前定(predetermined)** 的（即可以在上一期决定，状态变量），且同样要满足稻田条件。

* #### 资本积累方程

```math
K_{t+1}=I_t+(1-\delta)K_t \quad 0<\delta<1 \tag{1.3}
```

* 类似索洛模型，当期资本存量等于当期投资`$I_t$`加上上一期资本存量减去折旧（折旧率为`$\delta$`）。

* #### 资源约束

* 不考虑人口增长。没有（国外的）借贷，所以每期的产出只能消费或者储蓄；同时每个人只有有限的时间在工作和闲暇中做出选择：

```math
L_t+N_t \leq 1 \tag{1.4}

C_t+I_t \leq Y_t \tag{1.5}
```

* 所有变量在每一期都满足非负条件 。

* #### 社会计划者的问题

* 社会计划者实际上是要给出所有家庭和生产部门一条消费、投资的“路径”，使得社会的福利最大化。

* 以(1.2)到(1.5)所有的约束条件为前提，最终需要求解这样一个最大化问题：

```math
\max_{(C_t,L_t,Y_t,I_t,K_{t+1},N_t)^{\infty}_{t=0}}{\sum_{t=0}^{\infty}}\beta^t u(C_t,L_t) \quad  for \quad t=0,1,2 \cdots ,\infty \tag{1.6}
```

* 将约束条件代入(1.6)式，可得如下形式：

```math
\max_{(C_t,L_t,Y_t,I_t,K_{t+1},N_t)^{\infty}_{t=0}}{\sum_{t=0}^{\infty}}\beta^t u(F(K_t,N_t)+(1-\delta)K_t-K_{t+1},1-N_t)  \tag{1.7}
```

* #### 欧拉方程

* 对`$K_{t+1}$`求导：

```math
-u_C(C_t,1-N_t)+\beta \bigg\{u_C(C_{t+1},1-N_{t+1}) \times \bigg[F_K(K_{t+1},N_{t+1})+(1-\delta)\bigg]\bigg\}=0 \tag{1.8}
```

* **注意**：第一项是对当期（t期）中的`$K_{t+1}$`求导，而第二项是对下一期（t+1期）中的`$K_{t+1}$`求导（形式上表现为对下一期的`$K_{t}$`求导），因此第二项前有一个`$\beta$`。

* 对`$N_{t}$`求导：

```math
u_C(C_t,1-N_t)\times F_N(K_t,N_t)+u_L(C_t,1-N_t)\times (-1) \tag{1.9}
```

* 同时满足条件：

```math
C_t=F(K_t,N_t)+(1-\delta)K_t-K_{t+1}
```

* 第一个方程的含义：**消费的跨期安排**，当期放弃一单位的消费在下一期可以带来相同的效用。


* 第二个方程的含义：**劳动力的供给**，放弃一单位的闲暇可以获得恰好足以弥补效用损失的商品。

* 构成一个二阶（非线性）差分方程组。为什么是二阶？因为`$C_{t+1}$`中有`$K_{t+2}$`。

* #### 横截条件

* 为了保证方程组唯一解的存在性，需要引入横截条件，其经济含义是当趋于无穷期时资本存量的边际效用趋于零：

```math
\lim_{T \to 0} \beta^T u_C(C_T,1-N_T)K_{T+1}=0 \tag{1.10}
```

* #### 动态规划解法
* 可以将给定`$K_t$`时，中央计划者面临的决策问题写成以下递归形式（称为**贝尔曼方程**）：

```math
v(K_t)=\max_{0\leq K_{t+1},0<N_t<1} u\bigg(F(K_t,N_t)+(1-\delta)K_t-K_{t+1},1-N_t\bigg)+\beta v(K_{t+1}) \tag{1.11}
```

* 此时`$K_t$`给定，需要对`$K_{t+1}$`和`$N_t$`做出决策，解是一组决策的 **“政策函数”**。


* 最终得到的结果和欧拉方程结果一致（课件上命名为1.12-1.14，此处略）。

* #### 竞争均衡问题中的厂商

* 接下来，我们将目光转向**竞争均衡问题**的考察：

* 企业的价值：（共有`$J$`个厂商，厂商`$j$`在时间`$t$`的价值为`$Q_t^j$`，`$\pi_t^j$`是每股分红或利润，`$q_t^j$`是股价，`$S_t^j$`是发行股数，在不失一般性的情况下可以将其设定为1）

```math
Q_t^j=(q_t^j+\pi_t^j)S_t^j \tag{1.15}
```

* 企业利润的决定：（`$Y_t^j$`为对应产出，`$r_t^j$`为租金，`$K_t^j$`为使用的资本，`$w_t^j$`为工资，`$N_t^j$`为使用的劳动）

```math
\pi_t^j = Y_t^j-r_tK_t^j-w_tN_t^j \tag{1.16}
```

* 产出仍然遵循生产函数：

```math
Y_t \leq F(K_t,N_t) \tag{1.17}
```

* #### 竞争均衡问题中的家庭部门

* 类似于式(1.1)，家庭部门的偏好仍有：（同质）

```math
U_t^i=\sum_{t=0}^{\infty}{\beta}^t u(C_t^i,1-N_t^i) \tag{1.18}
```

* 面临以下的预算约束：

```math
C_t^i+K_{t+1}^i-(1-\delta)K_t^i+\sum_{j}q_t^js_{t+1}^{ij} \leq w_tN_t^i+r_tK_t^i+\sum_j(q_t^j+\pi_t^j)s_t^{ij} \tag{1.19}
```
* 类似地，还有所有变量的非负性约束。

* 其中当期家庭的资本存量`$K_t^i$`和投资比重`$s_t^{i1} \cdots s_t^{1j}$`是状态变量，而消费`$C_i$`和劳动`$N^i$`则是控制变量。

* #### 竞争均衡

* 最终的均衡表现为一系列的资源分配和价格使得各个微观主体效用最大化和市场出清：

```math
\sum_JY_t^i=\sum_I[C_t^i+K_{t+1}^i-(1-\delta)K_t^i] \tag{1.20}

\sum_JN_t^j=\sum_IN_t^i \tag{1.21}

\sum_JK_t^j=\sum_IK_t^i \tag{1.22}

\sum_Js_t^{ij}=1 \quad ,j=1,2,\cdots,J \tag{1.23}
```

#### ==局部近似==

* 二阶非线性差分方程组没有解析解，需要通过线性近似的方式求得近似解。

* 什么时候不适用局部近似：存在多重均衡（有多个稳态）、发生大的冲击（如金融危机）时（因为此时偏离稳态太远）

* 基本思想：使用**泰勒展开**的线性部分对某个点附近的函数值进行近似（局部近似，也存在全局近似方法，比较前沿）

* 稳态的定义——资本存量`$K$`、劳动力投入`$N$`和消费`$C$`均不再发生变化（写作`$\overline{K}$`、`$\overline{N}$`和`$\overline{C}$`）

* 考虑欧拉方程和约束条件：

```math
-u_C(C_t,1-N_t)+\beta \bigg\{u_C(C_{t+1},1-N_{t+1}) \times \bigg[F_K(K_{t+1},N_{t+1})+(1-\delta)\bigg]\bigg\}=0 \tag{1.24}

u_C(C_t,1-N_t)\times F_N(K_t,N_t)+u_L(C_t,1-N_t)\times (-1) \tag{1.25}

C_t=F(K_t,N_t)+(1-\delta)K_t-K_{t+1} \tag{1.26}
```

* 可以推导出在稳态下`$\overline{K}$`、`$\overline{N}$`和`$\overline{C}$`满足以下函数形式：

```math
F_K(\overline{K},\overline{N})=\frac{1}{\beta}-1+\delta \tag{1.27}

F_N(\overline{K},\overline{N}) = \frac{u_L(\overline{C},1-\overline{N})}{u_C(\overline{C},1-\overline{N})} \tag{1.28}

\overline{C}+\delta\overline{K}=F(\overline{K},\overline{N}) \tag{1.29}
```

* 在实际应用中往往不是对水平值进行线性化，而是**先对变量执行对数化**：（`$\hat{x}_t$`可以理解为偏离稳态的百分比）

```math
\hat{x}_t\equiv \log X_t-\log\overline{X} \tag{1.30}

\hat{x}_t\equiv \log X_t-\log\overline{X} = \log(\frac{X_t}{\overline{X}})=\log(1+\% change) \approx \% change \tag{1.31}
```

* #### 对数线性化的三块“砖头”：

* **第一块“砖头”——对`$X_t$`本身对数线性化**：

```math
X_t = \overline{X}(\frac{X_t}{X}) = \overline{X}\exp^{(\log X_t/\overline{X})} = \overline{X}\exp^{\hat{x}_t} \tag{1.32}
```

* 对上式进行一阶泰勒近似（围绕`$\hat{x}_t=0$`展开，代表在稳态附近展开）：

```math
X_t = \overline{X}\exp^{\hat{x}_t} \approx \overline{X}\exp^0+\overline{X}\exp^0(\hat{x}_t-0)=\overline{X}(1+\hat{x}_t) \tag{1.33}
```

* **第二块“砖头”——`$X_tY_t$`的对数线性化**：（`$\hat{x}_t\hat{y}_t$`是高阶项，舍去）

```math
\begin{aligned}
X_tY_t &\approx X(1+\hat{x}_t)Y(1+\hat{y}_t) \\
&= \overline{X} \overline{Y}(1+\hat{x}_t+\hat{y}_t+\hat{x}_t\hat{y}_t) \\
&\approx \overline{X} \overline{Y}(1+\hat{x}_t+\hat{y}_t) \tag{1.34}
\end{aligned}
```

* **第三块“砖头”——`$f(X_t)$`的对数线性化**：（其中`$\xi\equiv  \frac{f'(\overline{X})\overline{X}}{f(\overline{X})}$`，为弹性）

```math
\begin{aligned}
f(X_t) &\approx f(\overline{X})+f'(\overline{X})(X_t-\overline{X}) \\
&= f(\overline{X})+f'(\overline{X})\overline{X}(X_t/\overline{X}-1) \\
&= f(\overline{X})+f(\overline{X})\frac{f'(\overline{X})\overline{X}}{f(\overline{X})}(1+\hat{x}_t-1) \\
&= f(\overline{X})+f(\overline{X})\xi\hat{x}_t \\
&= f(\overline{X})(1+\xi\hat{x}_t) \tag{1.35} \\ 
\end{aligned}
```

* #### 线性方程组的求解：
* 在使用对数线性化方法对欧拉方程和约束条件（即式1.24-1.26）进行处理后，可将结果写为以下变量对稳态的偏离（即`$\hat{x_t}$`）形式：

```math
\xi_{cc}\hat{c}_t-\xi_{cl} \frac{\overline{N}}{1-\overline{N}} \hat{n}_t = \xi_{cc}\hat{c}_{t+1} + \bigg(1-\beta(1-\delta)\bigg)(\eta_{kn}\hat{n}_{t+1}+\eta_{kk}\hat{k}_{t+1}) \tag{1.36}

\eta_{nn}\hat{n}_t+\eta_{nk}\hat{k}_t = (\xi_{lc}-\xi_{cc})\hat{c}_t+(\xi_{cl}-\xi_{ll}) \frac{\overline{N}}{1-\overline{N}}\hat{n}_t \tag{1.37}

s_c\hat{c}_t+\frac{s_i}{\delta}\hat{k}_{t+1} = ((1-\alpha)+s_i\frac{1-\delta}{\delta})\hat{k}_t+\alpha\hat{n}_t \tag{1.38}
```
1. 其中`$\xi_{ab}$`代表在稳态时a的边际效用对b的弹性，即`$\xi_{ab} = \frac{\partial u_a(a,b)/u_a(a,b)}{\partial b/b}$`
2. 其中`$\eta_{ab}$`代表在稳态时a的边际产出对b的弹性，即`$\eta_{ab} = \frac{\partial F_a(a,b)/F_a(a,b)}{\partial b/b}$`
3. `$s_c$`和`$s_i$`分别代表消费和投资在稳态时占产出的份额
4. `$\alpha$`代表劳动在产出中的贡献，即劳动的产出弹性：`$\alpha=\frac{F_N(\overline{K},\overline{N})\overline{N}}{F(\overline{K},\overline{N})}$`

* 此时可以通过代入的方法消去变量`$\hat{n}_t$`，将剩下的部分整理后写成矩阵形式：

```math
\begin{bmatrix}
\hat{k}_{t+1} \\
\hat{c}_{t+1}
\end{bmatrix} = W \begin{bmatrix}
\hat{k}_{t} \\
\hat{c}_{t}
\end{bmatrix} \tag{1.39}
```

* #### 待定系数法求解差分方程（Uhlig 1999）

* 因为这是一个线性系统，其解一定是一个线性函数。我们可以将解（政策函数）写成下列形式：

```math
\hat{k}_{t+1} = \phi_1 \hat{k}_t \tag{1.40}

\hat{c}_t = \phi_2 \hat{k}_t \tag{1.41}
```

* 再将(1.40)和(1.41)代入(1.39)可得：（`$w_{ij}$`是`$W$`矩阵的第i行第j列元素）

```math
\begin{bmatrix}
\phi_1 \\
\phi_2\phi_1
\end{bmatrix} \hat{k}_t = \begin{bmatrix}
w_{11}+w_{12}\phi_2 \\
w_{21}+w_{22}\phi_2
\end{bmatrix} \hat{k}_t \tag{1.42}
```

* 可以解得：

```math
\phi_1 = w_{11}+w_{12}\phi_2 \tag{1.43}

\phi_2\phi_1 = w_{21}+w_{22}\phi_2 \tag{1.44}
```

* 当对模型中的参数赋值后，`$W$`内元素的值是已知的，因此可以据此求出政策函数(1.40)和(1.41)。

* 对比现实中美国宏观经济的数据，模型中还缺少两个要素：**变量的共同增长趋势、经济波动（随机项）**。

* 具体地，如果假定效用函数满足如下 **“可加可分”** 的形式，生产函数满足**柯布-道格拉斯函数**的形式，即有：（`$\theta_l$`是使N在稳态时达到0.2的常量，A是使Y在稳态时达到1的常量，这两个值不会对决策产生影响）

```math
u(C_t,L_t) = \log(C_t)+\theta_l\log(L_t) \tag{1.45}

Y_t  =AK_t^{1-\alpha}N_t^\alpha \tag{1.46}
```

* 进行对数线性化后得到：

```math
-\hat{c}_t = -\hat{c}_{t+1}+\bigg(1-\beta(1-\delta)\bigg)(\alpha\hat{n}_{t+1}-\alpha\hat{k}_{t+1}) \tag{1.47}

-(1-\alpha)\hat{n}_t+(1-\alpha)\hat{k}_t = \hat{c}_{t}+\frac{\overline{N}}{1-\overline{N}}\hat{n}_t \tag{1.48}

s_c\hat{c}_t+\frac{s_i}{\delta}\hat{k}_{k+1} = (1-\alpha+s_i \frac{1-\delta}{\delta})\hat{k}_t+\alpha\hat{n}_t \tag{1.49}
```

#### ==平衡增长路径==

* 为了实现变量的共同增长，在生产函数中添加劳动效率`$X_t$`：

```math
Y_t  =AK_t^{1-\alpha}(N_tX_t)^\alpha \tag{1.50}
```

* 此时`$X_tN_t$`称为有效劳动，技术的永久增长被限制在劳动效率`$X_t$`上。

* 在这样的假定之下，任意一个经济变量（包括投资`$I_t$`，消费`$C_t$`，资本存量`$K_t$`，产出`$Y_t$`，注意，不包括劳动`$N_t$`和闲暇`$L_t$`）`$Z_t$`的增长率都是恒定且相等的：

```math
Z_{t+1}/Z_t = \gamma_z \tag{1.51}
```

* 接下来还需对效用函数施加约束，要求效用函数是 **“可分的”**（具体证明见King(1998)的附录）：

```math
u(C_t,L_t) = \frac{C_t^{1-\sigma}}{1-\sigma}v(L_t) 
\quad if \quad \sigma\in(0,1)\bigcup(1,\infty) \tag{1.52}

u(C_t,L_t)=\log(C_t)+v(L_t) \quad if \quad \sigma=1 \tag{1.53}
```

* #### 平稳化转换

* 在引入劳动效率`$X_t$`后，由于变量具有了增长趋势，不再平稳，因此需要进行平稳化转换。具体地，定义`$i_t=I_t/X_t$`，`$y_t=Y_t/X_t$`，`$k_{t+1}=K_{t+1}/X_{t+1}$`，此时欧拉方程变为：（`$\gamma_x$`为`$X_t$`的增长率，即`$\gamma_x=\frac{X_{t+1}}{X_t}$`）

```math
-u_C(c_t,1-N_t)+\beta\gamma_x^{-\sigma}\times\{u_C(c_{t+1},1-N_{t+1})[(1-\alpha)A(\frac{k_{t+1}}{N_{t+1}})^{-\alpha}]+1-\delta\} = 0 \tag{1.54}

-u_I(c_t,1-N_t)+u_C(c_t,1-N_t)\alpha A(\frac{k_t}{N_t})^{1-\alpha} = 0 \tag{1.55}
```

* 还应该包括约束条件（平稳化转换后的形式）：

```math
c_t+\gamma_xk_{t+1}-A_tk_t^{1-\alpha}(N_t)^\alpha+(1-\delta)k_t = 0
```

* 同样采用可加可分形式的效用函数，对数线性化后得到：**(注意与1.47-1.49，即没有平衡增长路径的结果的区别)**

```math
-\hat{c}_t = -\hat{c}_{t+1}+(1-\frac{\beta}{\gamma_x}(1-\delta))(\alpha\hat{n}_{t+1}-\alpha\hat{k}_{t+1}) \tag{1.56}

-(1-\alpha)\hat{n}_t+(1-\alpha)\hat{k}_t = \hat{c}_{t}+\frac{\overline{N}}{1-\overline{N}}\hat{n}_t \tag{1.57}

s_c\hat{c}_t+\frac{s_i\gamma_x}{\delta+\gamma_x-1}\hat{k}_{k+1} = (1-\alpha+s_i \frac{1-\delta}{\delta+\gamma_x-1})\hat{k}_t+\alpha\hat{n}_t \tag{1.58}
```

#### ==引入不确定性==

* 最重要的一点——**将`$A$`替换成`$A_t$`**：（注意这里也引入了平衡增长路径）

```math
Y_t  =A_tK_t^{1-\alpha}(N_tX_t)^\alpha \tag{1.59}
```

* `$A_t$`是一个**外生随机过程**，服从 **AR(1)** 过程的形式：

```math
\begin{aligned}
A_t &= \overline{A}e^{a_t}\\
a_t &= \rho a_{t-1}+\epsilon_t \tag{1.60}
\end{aligned}
```

* 其他方程保持和之前一致：

* 效用函数（一般形式）：

```math
u(C_t,1-N_t) \quad (0<\beta<1) \tag{1.61}
```

* 资本积累：

```math
K_{t+1}=I_t+(1-\delta)K_t, \quad (0<\delta<1) \tag{1.62}
```

* 资源约束：

```math
C_t+I_t=F(K_t,N_t) \tag{1.63}
```

* 此时中央计划者需要最大化一个**期望效用**：（`$I_0$`是0时期的信息集，给定初始值`$K_0$`和`$A_o$`，控制变量满足非负条件）

```math
\max_{\{K_{t+1},N_t\}_{t=0}^{\infty}}E\bigg[\sum_{t=0}^{\infty}\beta^tu(F(K_t,N_t)+(1-\delta)K_t-K_{t+1},1-N_t|I_0\bigg] \tag{1.64}
```

* 求得的一阶条件如下（比之前多了个期望符号）：

```math
-u_C(C_t,1-N_t)+\beta^{-\sigma}\times E_t\bigg\{u_C(C_{t+1},1-N_{t+1})\bigg[(1-\alpha)A_{t+1}(\frac{K_{t+1}}{N_{t+1}})^{-\alpha}\bigg]+1-\delta\bigg\} = 0 \tag{1.65}

-u_I(C_t,1-N_t)+u_C(C_t,1-N_t)\alpha A_t(\frac{K_t}{N_t})^{1-\alpha} = 0 \tag{1.66}

\lim_{T \to 0} \beta^T u_C(C_T,1-N_T)K_{T+1}=0 \tag{1.67}
```

* 同时需要满足约束条件：

```math
C_t=F(K_t,N_t)+(1-\delta)K_t-K_{t+1}
```

* 接下来进行平稳化转换（比之前多了个期望符号）：

```math
-u_C(c_t,1-N_t)+\beta\gamma_x^{-\sigma}\times E_t\{u_C(c_{t+1},1-N_{t+1})[(1-\alpha)A(\frac{k_{t+1}}{N_{t+1}})^{-\alpha}]+1-\delta\} = 0 \tag{1.68}

-u_I(c_t,1-N_t)+u_C(c_t,1-N_t)\alpha A(\frac{k_t}{N_t})^{1-\alpha} = 0 \tag{1.69}

c_t+\gamma_xk_{t+1}-A_tk_t^{1-\alpha}(N_t)^\alpha+(1-\delta)k_t = 0 \tag{1.70}
```

* 需要寻找的政策函数：

```math
\begin{aligned}
k_{t+1} &= k(k_t,a_t) \\
c_{t} &= c(k_t,a_t) \\
N_{t} &= n(k_t,a_t)
\end{aligned}
```

* 使用同样的方法，进行对数线性化：

```math
-\xi_{cc}\hat{c}_t-\xi_{cl} \frac{\overline{N}}{1-\overline{N}} \hat{n}_t =E[ \xi_{cc}\hat{c}_{t+1} -\xi_{cl} \frac{\overline{N}}{1-\overline{N}} \hat{n}_{t+1}+ (1-\beta\gamma_x^{-\sigma}(1-\sigma))(a_{t+1}+\alpha\hat{n}_{t+1}-\alpha\hat{k}_{t+1}) \tag{1.71}]

a_{t}+(1-\alpha)\hat{n}_{t}+(1-\alpha)\hat{k}_{t} = (\xi_{lc}-\xi_{cc})\hat{c}_t+(\xi_{cl}-\xi_{ll}) \frac{\overline{N}}{1-\overline{N}}\hat{n}_t \tag{1.72}

s_c\hat{c}_t+\frac{s_i\gamma_x}{\delta+\gamma_x-1}\hat{k}_{k+1} = a_t+(1-\alpha+s_i \frac{1-\delta}{\delta+\gamma_x-1})\hat{k}_t+\alpha\hat{n}_t \tag{1.73}
```

* 然后写成矩阵形式：

```math
E\begin{bmatrix}
\hat{k}_{t+1} \\
a_{t+1} \\
\hat{c}_{t+1}
\end{bmatrix} = W\begin{bmatrix}
\hat{k}_{t} \\
a_t \\
\hat{c}_t
\end{bmatrix} \tag{1.74}
```

* 其中：

```math
W = \begin{bmatrix}
w_{11} &w_{12} &w_{13} \\
0 &\rho &0 \\
w_{31} &w_{32} &w_{33}
\end{bmatrix}
```

* 再次使用**待定系数法**，预设线性函数形式的解：

```math
\hat{k}_{t+1} = \phi_{11} \hat{k}_t +\phi_{12} a_t \tag{1.75}

\hat{c}_t = \phi_{21} \hat{k}_t +\phi_{22} a_t \tag{1.76}
```

* 整理为以下形式并对应求解：（选择一个满足横截条件的解）

```math
\begin{bmatrix}
\phi_{11} &\phi_{12}\\
0 &\rho \\
\phi_{21}\phi_{11} &\phi_{21}\phi_{12}+\phi_{22}\rho
\end{bmatrix}
\begin{bmatrix}
\hat{k}_t \\
a_t
\end{bmatrix} = 
\begin{bmatrix}
w_{11}+w_{12}\phi_{21} &w_{12}+w_{13}\phi_{22}\\
0 &\rho \\
w_{31}+w_{33}\phi_{21} &w_{32}+w_{33}\phi_{22}
\end{bmatrix}
\begin{bmatrix}
\hat{k}_t \\
a_t
\end{bmatrix} \tag{1.77} 

\phi_{11} = w_{11}+w_{12}\phi_{21} \tag{1.78}

\phi_{12} = w_{12}+w_{13}\phi_{22} \tag{1.79}

\phi_{21}\phi_{11} = w_{31}+w_{33}\phi_{21} \tag{1.80}

\phi_{21}\phi_{12}+\phi_{22}\rho = w_{32}+w_{33}\phi_{22} \tag{1.81}
```

#### ==状态空间方程==

* **任何DSGE模型都可以整理成为状态空间方程**

* 状态空间方程可以写成以下形式：（其中第一个方程称为**状态方程**，第二个方程称为**观测方程**，也就是政策函数）

```math
s_{t+1}  =Gs_t+Fe_{t+1} \tag{1.82}

z_t = Hs_t \tag{1.83}
```

* 其中`$s_t$`是状态变量(`$\hat{k}_t$`和`$a_t$`)，`$z_t$`包含了除状态变量以外的其他变量，`$e_t$`是技术冲击`$\epsilon_t$`。

#### ==总结==

* RBC模型的求解思路（待定系数法）：

> 1. 写出模型和最优化形式
> 2. 求出一阶条件（均衡条件，非线性差分方程组）
> 3. 对数线性化（泰勒展开），得到线性差分方程组
> 4. 使用待定系数法，求出政策函数（控制变量和状态变量的函数关系）

* RBC模型的拓展：

> 1. 平衡增长路径（加入劳动效率）
> 2. 加入随机波动（加入技术冲击即`$A_t$`）


### ==2.求解方法==

* 待定系数法只适用于比较简单的，状态变量较少的情况，当状态变量增多时，需要使用其他方法。

#### ==Blanchard-Kahn方法==

* 应用该方法需要满足**Blanchard-Kahn条件**（唯一解条件）：系数矩阵特征值分解后的**不稳定特征值**（即绝对值大于等于1的特征值）数量恰好等同于**控制变量**个数。

* Blanchard-Kahn方法是用于求解带期望的线性差分方程的方法，将方程写为如下形式（`$\epsilon_t$`是一个`$((m+n)\times1)$`向量，下同）：

```math
A
\begin{bmatrix}
x_{t+1} \\
E_t y_{t+1}
\end{bmatrix} = B
\begin{bmatrix}
x_{t} \\
y_{t}
\end{bmatrix} + C\epsilon_t \tag{2.1}
```
* 其中`$x_t$`是状态变量（n个），`$y_t$`是控制变量（m个）。在本模型中前者包括`$k_t,a_t$`，后者包括`$c_t,i_t,n_t,y_t$`。

* 接下来进行变换（**要求矩阵`$A$`可逆**）：

```math
\begin{bmatrix}
x_{t+1} \\
E_t y_{t+1}
\end{bmatrix} =A^{-1}B
\begin{bmatrix}
x_{t}\\ y_{t}
\end{bmatrix} + A^{-1}C\epsilon_t=F
\begin{bmatrix}
x_{t} \\
y_{t}
\end{bmatrix}+G\epsilon_t \tag{2.2}
```

* 使用特征值分解处理矩阵`$F$`：

```math
F=H \Lambda H^{-1} \tag{2.3}
```

* 可以通过选择特征根及特征向量的顺序，使矩阵`$\Lambda$`内的特征值从小到大排列：

```math
|\lambda_1|<|\lambda_2|<\cdots<|\lambda_{n+m}| \tag{2.4}
```

* 将(2.3)式和(2.4)式合并可得：

```math
\begin{bmatrix}
x_{t+1} \\
E_t y_{t+1}
\end{bmatrix} = H \Lambda H^{-1}+
\begin{bmatrix}
x_{t} \\
y_{t}
\end{bmatrix}+G\epsilon_t \tag{2.5}

H^{-1}
\begin{bmatrix}
x_{t+1} \\
E_t y_{t+1}
\end{bmatrix} = \Lambda H^{-1}+
\begin{bmatrix} x_{t} \\
y_{t}
\end{bmatrix}+H^{-1}G\epsilon_t \tag{2.6}
```

* 使用分块矩阵定义：

```math
H^{-1}  =
\begin{bmatrix}
H_{11} & H_{12} \\ 
H_{21} & H_{22}
\end{bmatrix} \tag{2.7}
```

* 代入式(2.6)可得： 

```math
\begin{bmatrix}
H_{11} & H_{12} \\ 
H_{21} & H_{22}
\end{bmatrix}
\begin{bmatrix}
x_{t+1}\\ E_t y_{t+1}
\end{bmatrix} = 
\begin{bmatrix}
\Lambda_1 & 0 \\ 
0 & \Lambda_2
\end{bmatrix}
\begin{bmatrix}
H_{11} & H_{12} \\ 
H_{21} & H_{22}
\end{bmatrix}
\begin{bmatrix}
x_{t} \\
y_{t}
\end{bmatrix} + 
\begin{bmatrix}
H_{11} & H_{12} \\ 
H_{21} & H_{22}
\end{bmatrix}
\begin{bmatrix}
G_{1} \\
G_{2}
\end{bmatrix}
\epsilon_t \tag{2.8}
```

* 注意式(2.8)中的一个事实——`$\Lambda_1$`的行数等同于`$F$`中稳定特征值的个数，`$\Lambda_2$`的行数等同于`$F$`中不稳定特征值的个数，因此可以将方程分为两部分求解——稳定部分和不稳定部分，此时Blanchard-Kahn条件的要求就显得十分直观了。

* 将上式再作合并简化：

```math
\begin{bmatrix}
\tilde{x}_{t} \\
\tilde{y}_{t}
\end{bmatrix} =
\begin{bmatrix}
H_{11} & H_{12} \\ 
H_{21} & H_{22}
\end{bmatrix}
\begin{bmatrix}
x_{t} \\
y_{t}
\end{bmatrix} \tag{2.9}

\begin{bmatrix}
\tilde{G}_{1} \\
\tilde{G}_{2}
\end{bmatrix} =
\begin{bmatrix}
H_{11} & H_{12} \\ 
H_{21} & H_{22}
\end{bmatrix}
\begin{bmatrix}
G_{1} \\
G_{2}
\end{bmatrix} \tag{2.10}
```

* 代入式(2.8)可得：

```math
\begin{bmatrix}
\tilde{x}_{t+1} \\
E\tilde{y}_{t+1}
\end{bmatrix} =
\begin{bmatrix}
\Lambda_1 & 0 \\ 
0 & \Lambda_2
\end{bmatrix}
\begin{bmatrix}
\tilde{x}_{t} \\
\tilde{y}_{t}
\end{bmatrix} +
\begin{bmatrix}
\tilde{G}_{1} \\
\tilde{G}_{2}
\end{bmatrix} \epsilon_t \tag{2.11}
```

* 方程的稳定部分可用一般方法求解，以下考虑不稳定部分，有：

```math
E_t\tilde{y}_{t+1} = \Lambda_2\tilde{y}_t+\tilde{G}_2\epsilon_t \tag{2.12}
```

* 通过移项和向前一期可得：

```math
\tilde{y}_t = \Lambda_2^{-1}E_t\tilde{y}_{t+1}=\Lambda_2^{-1}\tilde{G}_2\epsilon_t \tag{2.13}

\tilde{y}_{t+1} = \Lambda_2^{-1}E_{t+1}\tilde{y}_{t+2}=\Lambda_2^{-1}\tilde{G}_2\epsilon_{t+1} \tag{2.14}
```

* 接下来将式(2.14)代入式(2.13)，使用了**期望迭代法则**（tower性质：保留较小的信息集）：

```math
\begin{aligned}
\tilde{y}_t &= \Lambda_2^{-1}E_{t+1}(\tilde{y}_{t+2}-\Lambda_2^{-1}\tilde{G}_2\epsilon_{t+1})-\tilde{G}_2\epsilon_{t} \\
&= \Lambda_2^{-2}E_{t}\tilde{y}_{t+2}-\Lambda_2^{-2}E_t(\tilde{G}_2\epsilon_{t+1})-\Lambda_2^{-1}\tilde{G}_2\epsilon_{t} \\
&\cdots \\
&= -\Lambda_2^{-1}\tilde{G}_2\epsilon_{t}-\Lambda_2^{-2}E_t(\tilde{G}_2\epsilon_{t+1})-\Lambda_2^{-3}E_t(\tilde{G}_2\epsilon_{t+2})\cdots \tag{2.15-16}
\end{aligned}
```

* 可以发现`$E_t\epsilon_{t+1}=E_t\epsilon_{t+2}=\cdots=0$`，上式可化为：

```math
\tilde{y}_t = -\Lambda_2^{-1}\tilde{G}_2\epsilon_{t} \tag{2.17}
```

* 接下来从式(2.9)出发，代入式(2.17)：**（得到的式(2.18)可以认为就是状态空间方程中的观测方程）**

```math
\begin{aligned}
\tilde{y}_t &= H_{21}x_t+H_{22}y_t \\
\Rightarrow y_t &= H_{22}^{-1}\tilde{y}_t-H_{22}^{-1}H_{21}x_t \\
\Rightarrow y_t &= -H_{22}^{-1}\Lambda_2^{-1}\tilde{G}_2\epsilon_{t}-H_{22}^{-1}H_{21}x_t \tag{2.18}
\end{aligned}
```

* 再从式(2.8)出发，代入式(2.18)：**（得到的式(2.19)可以认为是状态空间方程中的状态方程）**

```math
\begin{aligned}
x_{t+1} &= F_{11}x_t+F_{22}y_t+G_1\epsilon_t \\
\Rightarrow x_{t+1} &= F_{11}x_t+F_{22}(-H_{22}^{-1}\Lambda_2^{-1}\tilde{G}_2\epsilon_{t}-H_{22}^{-1}H_{21}x_t)+G_1\epsilon_t \\
\Rightarrow x_{t+1} &= (F_{11}-H_{22}^{-1}H_{21})x_t+(G_1-F_{22}-H_{22}^{-1}\Lambda_2^{-1})\epsilon_t \tag{2.19}
\end{aligned}
```

#### ==广义Schur法（QZ算法）==

* 先考察Hansen(1985)的一个模型：

```math
\overline{K}\hat{K}_{t+1} = \overline{Y}\hat{Y}_t-\overline{C}\hat{C}_t+(1-\delta)\overline{K}\hat{K}_t \tag{2.20}

\hat{\lambda}_t = \gamma\hat{\lambda}_{t-1}+\epsilon_t \tag{2.21}

0 = \hat{\lambda}_{t}-\theta\hat{Y}_t+\theta\hat{K}_t-(1-\theta)\hat{C}_t \tag{2.22}

0 = \hat{K}_t+\hat{r}_t-\hat{Y}_t \tag{2.23}

\hat{C}_t = E_t\hat{C}_{t+1}-\beta\overline{r}E_t\hat{r}_{t+1} \tag{2.24}
```

* 整理成状态空间方程的矩阵形式：

```math
\begin{bmatrix}
x_{t+1} \\
E_ty_{t+1}
\end{bmatrix} = 
\begin{bmatrix}
\hat{K}_{t+1} \\
\hat{\lambda}_t \\
\hat{Y}_t \\
E_t\hat{C}_{t+1} \\
E_t\hat{r}_{t+1} \tag{2.25}
\end{bmatrix}

A
\begin{bmatrix}
x_{t+1} \\
E_ty_{t+1}
\end{bmatrix} = B
\begin{bmatrix}
x_{t} \\
y_{t}
\end{bmatrix} + G\epsilon_t \tag{2.26}
```

* 其中：（可以看出矩阵`$A$`是不可逆的）


```math
A =
\begin{bmatrix}
\overline{K} &0 &-\overline{Y} &0 &0 \\
0 &1 &0 &0 &0 \\
0 &-1 &\theta &0 &0 \\
0 &0 &1 &0 &0 \\
0 &0 &0 &1 &-\overline{r}\beta \tag{2.27}
\end{bmatrix}
```

* 对于这类问题，我们需要使用广义Schur法处理，其原理类似于Blanchard-Kahn方法，可以同时拆分矩阵`$A$`和`$B$`。

* 首先，使用**QZ分解**处理矩阵`$A$`和`$B$`：

```math
A = QTZ',\quad B=QSZ'
```

* 其中矩阵`$Q$`和`$Z$`是正交矩阵（其转置等于其逆矩阵），矩阵`$S$`和`$T$`是上三角矩阵：

```math
QQ'=Q^TQ=I=ZZ'=Z'Z \tag{2.28}
```

* 定义系统的特征根为：（其中`$s_{ii}$`和`$t_{ii}$`是矩阵`$S$`和`$T$`的对角线元素，同样按照绝对值从小到大排序）

```math
\lambda_{ii}=\frac{s_{ii}}{t_{ii}} \tag{2.29}
```

* 模型的确定部分（不包含`$\epsilon_t$`）有：

```math
QTZ'
\begin{bmatrix}
x_{t+1} \\
E_t y_{t+1}
\end{bmatrix} = QSZ'
\begin{bmatrix}
x_{t} \\
y_{t}
\end{bmatrix} \tag{2.30}
```

* 两边同乘`$Q'$`，并将矩阵根据状态变量和控制变量分块可得：

```math
\begin{bmatrix}
T_{11} & T_{12} \\ 
0 & T_{22}
\end{bmatrix}
\begin{bmatrix}
Z'_{11} & Z'_{12} \\ 
Z'_{21} & Z'_{22}
\end{bmatrix}
\begin{bmatrix}
x_{t+1}\\ E_t y_{t+1}
\end{bmatrix} = 
\begin{bmatrix}
S_{11} & S_{12} \\ 
0 & S_{22}
\end{bmatrix}
\begin{bmatrix}
Z'_{11} & Z'_{12} \\ 
Z'_{21} & Z'_{22}
\end{bmatrix}
\begin{bmatrix}
x_{t} \\
y_{t}
\end{bmatrix} + 
\begin{bmatrix}
Q'_{11} & Q'_{12} \\ 
Q'_{21} & Q'_{22}
\end{bmatrix}
\begin{bmatrix}
G_{1} \\
G_{2}
\end{bmatrix}
\epsilon_t \tag{2.31}
```

* 与Blanchard-Kahn方法类似，我们将不稳定特征值放在下方，因此式(2.31)的不稳定部分可写作：

```math
0 = S_{22}Z'_{21}x_t+S_{22}Z'_{22}y_t+[Q'_{21}G_1+Q'_{22}G_2]\epsilon_t \tag{2.32}

\begin{aligned}
\Rightarrow y_t &= -[S_{22}Z'_{22}]^{-1}S_{22}Z'_{21}x_t-[S_{22}Z'_{22}]^{-1}[Q'_{21}G_1+Q'_{22}G_2]\epsilon_t \\
&= -(Z'_{22})^{-1}Z'_{21}x_t-(Z'_{22})^{-1}S_{22}^{-1}[Q'_{21}G_1+Q'_{22}G_2]\epsilon_t \tag{2.33}
\end{aligned}
```

* 做记号：

```math
N = (Z'_{22})^{-1}Z'_{21} \tag{2.34}

L = (Z'_{22})^{-1}S_{22}^{-1}[Q'_{21}G_1+Q'_{22}G_2] \tag{2.35}
```

* 此时方程的控制变量部分可写作 **（观测方程）**：

```math
y_t = -N_tx_t-L\epsilon_t \tag{2.36}
```

* 对上式求期望可得：

```math
E_ty_{t+1} = -(Z'_{22})^{-1}Z'_{21}x_{t+1} = -Nx_{t+1} \tag{2.37}
```

* 将式(2.37)代入式(2.26)：

```math
\begin{bmatrix}
A_{11} & A_{12} \\ 
A_{21} & A_{22}
\end{bmatrix}
\begin{bmatrix}
I \\ 
-N
\end{bmatrix} x_{t+1} =
\begin{bmatrix}
B_{11} & B_{12} \\ 
B_{21} & B_{22}
\end{bmatrix}
\begin{bmatrix}
x_{t} \\
y_{t}
\end{bmatrix} +
\begin{bmatrix}
G_{1} \\
G_{2}
\end{bmatrix}
\epsilon_t \tag{2.38}
```

* 再将式(2.36)代入式(2.38)：

```math
\begin{bmatrix}
A_{11} & A_{12} \\ 
A_{21} & A_{22}
\end{bmatrix}
\begin{bmatrix}
I \\ 
-N
\end{bmatrix} x_{t+1} =
\begin{bmatrix}
B_{11} & B_{12} \\ 
B_{21} & B_{22}
\end{bmatrix}
\begin{bmatrix}
I \\ 
-N
\end{bmatrix} x_{t} +
\begin{bmatrix}
G_{1}-B_{12}L \\
G_{2}-B_{22}L
\end{bmatrix}
\epsilon_t \tag{2.39}
```

* 展开稳定的（上半部分）方程可得：

```math
[A_{11}-A_{12}N]x_{t+1} = [B_{11}-B_{12}N]x_t+[G_1-A_{12}L]\epsilon_t \tag{2.40}

\Rightarrow x_{t+1} = [A_{11}-A_{12}N]^{-1}[B_{11}-B_{12}N]x_t+[A_{11}-A_{12}N]^{-1}[G_1-A_{12}L]\epsilon_t \tag{2.41}
```

* 这就是稳定部分的解 **（状态方程）**。

### ==3.校准==

* 校准更加像是一种通过选取参数使理论模型“接近”真实的数据生成过程、“匹配”数据特征（典型事实）的方法，不同于计量中的参数估计。

* 宏观研究中的四类问题结构：
1. How much of fact X can be explained with impulses of type Y?
2. Is it possible to generate features F by using theory T?
3. Can we reduce the discrepancy D of the theory from the data by using feature F?
4. How much do endogenous variables change if the process for the exogenous variables is altered?

* **校准的常用方法**：

1. 利用稳态条件，让稳态的经济变量和实际经济运行匹配
2. 利用微观计量的相关研究，如调查数据
3. 使用估计的方法，即运用计量方法，如GMM、ML等

* **匹配矩条件**：模型校准后，可以通过状态空间模型表示出来，也就是一个时间序列模型。可以通过考察模型模拟数据的矩和真实数据的矩的差异，对模型准确性进行评估。

### ==4.MIU（Money-in-Utility）模型和CIA（Cash-in-Advance）模型==

#### ==MIU模型==

* **在效用函数中添加货币因素**：（其中`$m_t$`是实际货币余额——名义货币的购买力：`$m_t=\frac{M_t}{P_t}$`）

```math
U = E_0\sum_{t=0}^{\infty}\beta^tu(C_t,m_t) \quad 0<\beta<1 \tag{4.1}
```

* 效用函数仍然满足凸性、二阶可微可导和稻田条件。同时注意闲暇被实际货币余额取代了。

* 资本积累方程没有改动：

```math
K_{t+1}=I_t+(1-\delta)K_t \quad 0<\delta<1 \tag{4.2}
```

* 家庭部门面临的预算约束：（`$i_t$`表示**名义利率**，`$t_t$`表示一次付清的转移支付）

```math
C_t+I_t+\frac{B_t}{P_t}+\frac{M_t}{P_t} \leq w_t+r_tK_t+(1+i_{t-1})\frac{B_{t-1}}{P_t}+\frac{M_{t-1}}{P_t}+\frac{T_t}{P_t} \tag{4.3}
```

* 将价格写成通货膨胀`$\pi$`的形式，上式可表示为：（要求`$b_t=\frac{B_t}{P_t}>-\overline{b}$`，避免无限借债出现庞氏骗局）

```math
C_t+I_t+b_t+m_t \leq w_t+r_tK_t+\frac{1+i_{t-1}}{1+\pi_t}b_{t-1}+\frac{m_{t-1}}{1+\pi_t}+t_t \tag{4.4}
```

* 家庭需要对`$\{C_t,b_t,m_t,K_{t+1}\}_{t=0}^{\infty}$`进行决策，每期的`$\{\pi_t,i_t,w_t,r_t,t_t\}_{t=0}^{\infty}$`以及这些变量的初值是已知的。

* 厂商的设定仍然沿用柯布-道格拉斯生产函数和AR(1)形式的全要素生产率，即：

```math
Y_t  =A_tK_t^{1-\alpha}(N_t)^\alpha \tag{4.5}

A_t = \overline{A}e^{a_t} \tag{4.6}

a_t = \rho^\alpha a_{t-1}+\epsilon_t^a \tag{4.7}
```

* 厂商需要通过选择`$N_t,K_t$`最大化利润`$Y_t-r_tK_t^{d}-w_tN_t^d$`：

* 货币的发行者——政府的预算约束：（左侧为政府的收入项，右侧是政府的支出项）

```math
\frac{M_t^s}{P_t}-\frac{M_t^s}{P_t}+\frac{B_{t}^s}{P_t} = \frac{T_t}{P_t} + (1+i_{t-1})\frac{B_{t-1}^s}{P_t} \tag{4.8}

m_t^s-\frac{m_{t-1}^s}{1+\pi_t}+b_t^s- \frac{1+i_{t-1}}{1+\pi_t}b_{t-1}^s=t_t \tag{4.9}
```

* 如果政府不发行债务，则可以得到一个更简单的方程：

```math
t_t = m_t^s-\frac{m_{t-1}^s}{1+\pi_t} \tag{4.10}
```

* 政府发行货币的数量规则：（其中`$\theta_t=\frac{M_t}{M_{t-1}}-\mu-1$`）

```math
\theta_t = \rho^\theta \theta_{t-1}+\epsilon_t^\theta \tag{4.11}
```

* 结合约束条件，可以得到：（作为均衡条件还需要加上预算约束`$C_t=A_tK_t^{1-\alpha}+(1-\delta)K_t-K_{t+1}$`）

```math
-u_C(C_t,m_t)+\beta E_t\bigg[u_C(C_{t+1},M_{t+1}((1-\alpha)A_{t+1}K_{t+1}^{-\alpha}+1-\delta)\bigg] = 0 \tag{4.12}

-u_C(C_t,m_t)+\beta E_t\bigg[u_C(C_{t+1},m_{t+1})\frac{1+i_t}{1+\pi_{t+1}}\bigg] = 0 \tag{4.13}

u_m(C_t,m_t)+\beta E_t\bigg[u_c(C_{t+1},m_{t+1})\frac{1}{1+\pi_{t+1}}\bigg]-u_C(C_t,m_t)=0 \tag{4.14}
```

* 其中式(4.12)为消费-投资替代的欧拉方程，式(4.13)为消费-债券替代的欧拉方程，式(4.14)为货币需求函数（消费-货币替代方程）。

* 同时需要满足横截条件（保证唯一解）：

```math
\lim_{T \to \infty} \beta^TE_t\bigg[u_C(C_T,m_T)K_{T+1}\bigg]=0 

\lim_{T \to \infty} \beta^TE_t\bigg[u_C(C_T,m_T)b_T\bigg]=0

\lim_{T \to \infty} \beta^TE_t\bigg[u_C(C_T,m_T)m_T\bigg]=0
```

* 将式(4.13)代入式(4.14)可得：

```math
\begin{aligned}
u_m(C_t,m_t) &= -\beta E_t\bigg[u_C(C_{t+1},m_{t+1})\frac{1}{1+\pi_{t+1}}\bigg] \\
\frac{u_m(C_t,m_t)}{u_C(C_t,m_t)} &= 1-\beta E_t\bigg[\frac{u_C(C_{t+1},m_{t+1})}{u_C(C_{t},m_{t})}\frac{1}{1+\pi_{t+1}}\bigg] \\
\frac{u_m(C_t,m_t)}{u_C(C_t,m_t)} &= 1-\frac{1}{1+i_t} \\
\frac{u_m(C_t,m_t)}{u_C(C_t,m_t)} &= \frac{i_t}{1+i_t} \tag{4.15}
\end{aligned}
```

* 式(4.15)刻画了持有货币的机会成本，也可以认为是**货币需求函数**（实际货币需求量与名义利率的关系）。

* 从式(4.12)可以推出一个稳态条件：

```math
\beta((1-\alpha)\overline{A}\overline{K}^{-\alpha}+1-\delta) = 1 \tag{4.16}
```

* 同时根据另外几个稳态条件：`$\overline{I}=\delta K$`，`$\overline{Y}=\overline{A}\overline{K}^{1-\alpha}$`，`$\overline{C}+\overline{I}=\overline{Y}$`，可以解出上述所有稳态变量的值，可以发现此时货币变量`$m_t$`对稳态无影响，表现出货币中性。

* 此时考虑`$m_t$`的稳态，如果存在稳态，则必须存在一个正值`$\overline{m}$`满足：

```math
u_m(\overline{C},\overline{m}) = (1-\frac{\beta}{1+\mu})u_C(\overline{C},\overline{m})
```

* 为了保证解的存在性，我们还需要对效用函数施加可加可分的条件：

```math
u(C,m)=u^1(C)+u^2(m) \to u_m^2(\overline{m}) = (1-\frac{\beta}{1+\mu})u_c^1(\overline{C})>0
```

* 货币中性：货币发行数量（水平值）不会影响实际变量，这种情况的模型也称为**弹性价格模型**。

* 货币超中性：不仅货币发行数量（水平值）不会影响实际变量，货币发行速度（增长率）也不影响实际变量。（MIU模型是一个**短期中性**，**长期超中性**的模型，后者是因为在稳态（长期）中，货币增长速度`$\mu$`并不影响实际变量）

* 如果选择以下的效用函数形式（可加可分）：

```math
u(C_t,m_t) = \frac{C_t^{1-\alpha}-1}{1-\sigma}+\phi\frac{m_t^{1-\chi}-1}{1-\chi}
```

* 可以推出如下的均衡条件：

```math
-C_t^{-\sigma}+\beta E_t\bigg[C_{t+1}^{-\sigma}(1-\alpha)A_{t+1}K_{t+1}^{-\alpha}+1-\delta\bigg] = 0

-C_t^{-\sigma}+\beta E_t\bigg[C_{t+1}^{-\sigma}\frac{1+i_t}{1+\pi_{t+1}}\bigg] = 0

\phi m_t^{-\chi}-\frac{i_t}{1+i_t}C_t^{-\sigma} = 0

A_tK_t^{1-\alpha}+(1-\delta)K_t-K_{t+1} = C_t
```

* 其中：

```math
\frac{M_t}{M_{t-1}} = 1+\mu+\theta

\theta_t = \rho^{\theta}\theta_{t-1}+\epsilon_t^{\theta}

\log(A_t) = (1-\rho^{\alpha})\log(\overline{A})+\rho^{\alpha}\log A_{t-1}+\epsilon_t^\alpha
```

* 对数线性化后有：

```math
-\sigma \hat{c}_t = -\sigma E_t \hat{c}_{t+1}+(1-\beta(1-\delta))E_t\bigg[a_{t+1}-\alpha\hat{k}_{t+1}\bigg]

-\sigma \hat{c}_t = -\sigma E_t \hat{c}_{t+1}+\hat{i}_t-E_t\hat{\pi}_{t+1}

-\chi\hat{m}_t = -\sigma\hat{c}_t+(\frac{\beta}{1+\mu-\beta})\hat{i}_t

s_x\hat{c}_t+\frac{s_i}{\delta}\hat{k}_{t+1} = a_t+\bigg((1-\alpha)+s_i\frac{1-\delta}{\delta}\bigg)\hat{k}_t

\hat{m}_t = \hat{m}_{t-1}-\hat{\pi}_t+\theta_t

\theta_t = \rho^{\theta}\theta_{t-1}+\epsilon_t^{\theta}

a_t = \rho^\alpha a_{t-1}+\epsilon_t^\alpha
```

#### ==CIA模型==

* 另外一种将货币引入RBC模型的方法是**增加约束条件（式4.18）**，此时效用函数不改变：（需要满足可加可分的形式）

```math
U=E_0\sum_{t=0}^{\infty}{\beta}^t u(C_t,L_t) \tag{4.17}
```

* 家庭的消费必须在事先持有货币：

```math
P_tC_t \leq M_{t-1}+T_t \tag{4.18}
```

* 同样地，给出家庭的预算约束和资本积累方程：（和MIU模型的不同之处在于由于CIA模型保留了劳动和闲暇，因此此时的`$w_t$`代表工资率，需要乘上`$N_t$`）

```math
C_t+I_t+b_t+m_t \leq w_tN_t+r_tK_t+\frac{1+i_{t-1}}{1+\pi_t}b_{t-1}+\frac{m_{t-1}}{1+\pi_t}+t_t \tag{4.19}

K_{t+1}=I_t+(1-\delta)K_t \quad 0<\delta<1 \tag{4.20}
```

* 家庭需要对`$\{C_t,b_t,m_t,N_t,K_{t+1}\}_{t=0}^{\infty}$`进行决策，每期的`$\{\pi_t,i_t,w_t,r_t,t_t\}_{t=0}^{\infty}$`以及这些变量的初值是已知的。

* 可以得到以下的均衡条件：（模型中有两个约束条件，因此存在两个拉格朗日乘子，其中`$\lambda_t$`代表预算约束或者财富的边际效用，`$\Xi_t$`代表现金约束或者流动性的价值）

```math
\lambda-u_C(C_t,1-N_t)+\Xi_t = 0 \tag{4.21}

-\lambda+\beta E_t\bigg[\lambda_{t+1}(1-\alpha)A_{t+1}(\frac{K_{t+1}}{N_{t+1}})^{-\alpha}+1-\delta\bigg] = 0 \tag{4.22}

-u_I(C_t,1-N_t)+\lambda_t\alpha A_t(\frac{K_{t}}{N_{t}})^{1-\alpha} = 0 \tag{4.23}

-\lambda_t+\beta E_t\bigg[\lambda_{t+1}\frac{1+i_t}{1+\pi_{t+1}}\bigg] = 0 \tag{4.24}

-\lambda_t+\beta E_t\bigg[u_C(C_{t+1},1-N_{t+1})\frac{1}{1+\pi_{t+1}}\bigg] = 0 \tag{4.25}

m_t = C_t \tag{4.26}
```

* 其中式(4.21)到(4.25)是由拉格朗日函数分别求导获得，式(4.26)是由式(4.10)和式(4.18)获得（同时满足`$C_t=A_t K_t^{1-\alpha}N_t^\alpha +  (1-\delta)K_t-K_{t+1}$`的约束条件。）。

* 货币需求函数可以表示为：

```math
\beta E_t\bigg[\frac{\Xi}{1+\pi_{t+1}}\bigg] = \frac{i_t}{1+i_t}\lambda_t
```

### ==5.Dynare实操==
> [Dynare手册](https://www.dynare.org/manual/)

* 一个Dynare（校准）程序的**四个要素**：
1. 变量和参数：var/varexo
2. 模型：model; end; **注意内生变量个数需要等于方程个数**
3. 初始值：initval;
4. 冲击：shocks; end;

* 一个需要注意的地方：在描述资本存量（以及其他状态变量）`$K_t$`时，我们通常认为指的是t期期初，但在Dynare中要写成`$K_{t-1}$`代表t-1期期末。

* 一个简单RBC模型的Dynare实现：

```matlab
var y i k a c;%指定变量名，Dynare会自动判断变量是控制变量还是状态变量
varexo e;%外生变量（当作冲击项）

parameters alpha beta delta1 rho sigma sigmae;%指定参数并赋值

alpha = 0.33;
beta = 0.99;
delta1 = 0.025;
rho = 0.95;
sigma = 1;
sigmae = 0.01;

model;
% 由于Dynare会自动执行线性化，而我们需要的结果是对数线性化，因此在实际运行程序时变量的值均为对数值；对应地，满足的模型条件需要将变量取指数。
exp(c)^(-sigma) = beta*(exp(c(+1))^(-sigma)*(alpha*exp(a(+1))*exp(k)^(alpha-1)+(1-delta1)));
exp(y) = exp(a)*exp(k(-1))^(alpha);
exp(k) = exp(i)+(1-delta1)*exp(k(-1));
exp(y) = exp(c)+exp(i);
a = rho*a(-1)+e;
end;

initval;%指定变量初值
%不指定初值也是可以的，但在一些规模较大的模型中，指定初值可以提高运算速度
k = log(29);
y = log(3);
a = 0;
c = log(2.5);
i = log(1.5);
end;

shocks;
var e = sigmae^2;
end;

steady;

stoch_simul(hp_filter=1600, order=1, irf=40); 
%模拟求解，其中指定H-P滤波的lambda值为1600，一阶近似（即线性化），脉冲响应图像绘制40期
```

* 模型程序写好后，还需在MATLAB的命令行中执行：（filename表示文件名）

```matlab
dynare filename.mod
```

#### ==Dynare运行结果的存储和调用==

* 变量x对变量e的脉冲响应函数储存在x_e中

* 写出状态空间方程形式的模型：

```math
\begin{cases}
\begin{aligned}
s_{t} &= As_{t-1}+Be_{t+1} \\
x_t &= \Phi s_t \\
\end{aligned}
\end{cases}
\\
\begin{aligned}
\Rightarrow x_t &= \Phi As_{t-1}+\Phi Be_{t+1} \\
\Rightarrow \begin{bmatrix} s_t \\
x_t
\end{bmatrix} &= \begin{bmatrix} A \\
\Phi A
\end{bmatrix}s_{t-1} +
\begin{bmatrix} B \\
\Phi B
\end{bmatrix} \epsilon_t
\end{aligned}
```

* 矩阵`$\Psi = \begin{bmatrix} A \\ \Phi A \end{bmatrix}$`储存在**oo_.dr.ghx**中，矩阵`$\Omega = \begin{bmatrix} B \\ \Phi B \end{bmatrix}$`储存在**oo_.dr.ghu**中，稳态值存储在**oo_.dr.ys**中。

* 一个基于上述参数的Matlab模拟程序（画图部分省略）：

```Matlab
m = 2; % number of states
n = 3; % number of controls
psi = oo_.dr.ghx;
omega = oo_.dr.ghu;
ss = oo_.dr.ys;

T = 200; % number of periods to simulate

es = sigmae*randn(T,1);
Xsim = zeros(n+m,T);

Xsim(:,1) = omega*es(1,1);

for j = 2:T
   Xsim(:,j) = psi*Xsim(3:4,j-1) + omega*es(j,1);
end
```

### ==6.NK（新凯恩斯主义）模型==

* RBC模型和VAR证据不符合，因此我们需要引入新的机制和假设使模型体现出货币非中性。一种方法是保留弹性价格，引入不完全信息、有限参与等机制；另一种更常用的方法是直接引入凯恩斯主义的**名义刚性**。

* NK模型一个重要的特点在于引入了垄断竞争——**赋予企业定价权**。

#### ==模型基本设定==

* #### 效用函数

* 在MIU模型的基础上对效用函数进行拓展：

```math
U = E_0\sum_{t=0}^{\infty}\beta^tu(C_t,1-N_t,m_t) \quad 0<\beta<1 \tag{6.1}
```

* #### Dixit-Stiglitz加总

* 商品是一个0到1的**连续统**（continuum），且每个商品和每个厂商一一对应。`$\theta$`刻画了每个商品之间的的**替代弹性**：

```math
C_t = \bigg[\int_0^1C_t(i)^{\frac{\theta-1}{\theta}}di\bigg]^{\frac{\theta}{\theta-1}} \quad \theta>1 \tag{6.2}
```

* 其中`$C_t(i)$`代表商品i在t期的消费数量，`$C_t$`相当于一个加总后的“消费包”（大于1说明能提供大于一个“消费包”的效用）。我们定义购买一个消费包的价格为`$P_t$`，每个商品的价格为`$P_t(i)$`，那么对消费包内部具有以下优化条件：(此时`$P_t$`正是这个优化问题的拉格朗日乘子)

```math
\min_{c_t(i)} \int_0^1P_t(i)C_t(i)di \quad s.t.  \quad C_t \geq 1 \tag{6.3}

\begin{aligned}
\Rightarrow \mathcal{L} &= \int_0^1P_t(i)C_t(i)di - P_t(C_t-1) \\
&= \int_0^1P_t(i)C_t(i)di - P_t(\bigg[\int_0^1C_t(i)^{\frac{\theta-1}{\theta}}di\bigg]^{\frac{\theta}{\theta-1}}) \\
\Rightarrow \frac{\partial L}{\partial C_t(i)} &= P_t(i)-P_t\times(\frac{\theta}{\theta-1})\times\bigg[\int_0^1C_t(i)^{\frac{\theta-1}{\theta}}di\bigg]^{\frac{\theta}{\theta-1}-1}\times \frac{\theta-1}{\theta}C_t(i)^{\frac{\theta-1}{\theta}-1} = 0
\end{aligned}

\Rightarrow P_t(i) = P_t(\frac{C_t(i)}{C_t})^{-\frac{1}{\theta}} \tag{6.4}

\Leftrightarrow C_t(i) = C_t(\frac{P_t(i)}{P_t})^{-\theta} \tag{6.5}
```
* 式(6.4)和式(6.5)以`$C_t=1$`为前提。可以看出，相对价格`$\frac{P_t(i)}{P_t}$`是单个商品购买数量的决定因素 **（式(6.5)也是每个企业的生产函数）**。

* 将式(6.4)代入式(6.3)可得：

```math
\begin{aligned}
\int_0^1P_t(i)C_t(i)di &= P_i(\int_0^1C_t(i)^{\frac{\theta-1}{\theta}}di)C_t^{\frac{1}{\theta}} \\
&= P_tC_t^{\frac{\theta-1}{\theta}}C_t^{\frac{1}{\theta}} \\
&= P_tC_t \\
&= P_t
\end{aligned}
```

* 将式(6.5)代入式(6.3)：（这被称为**Dixit-Stiglitz价格指数**，代表整个社会的物价指数`$P_t$`是如何由每个商品的价格组成的）

```math
P_t = (\int_0^1P_t(i)^{1-\theta}di)^{\frac{1}{1-\theta}}
```

* #### 预算约束

* 家庭部门的预算约束写作：（符号沿用之前的设定，`$W_t$`是名义工资，`$T_t$`是一次付清的转移支付，`$Q_t$`是从公司获得的名义利润）

```math
\int_0^1P_t(i)C_t(i)di+B_t+M_t \leq W_iN_t+M_{t-1}+(1+i_{t-1})B_{t-1}+T_t+Q_t \tag{6.6}
```

* 由于我们已经将消费品组合并优化成为消费包，所以式(6.6)可以写作：

```math
P_tC_t+B_t+M_t \leq W_iN_t+M_{t-1}+(1+i_{t-1})B_{t-1}+T_t+Q_t \tag{6.7}
```

* #### 均衡条件

* 此时消费者的跨期决策问题转化为通过选择`$\{C_t,m_t,N_t\}_{t=0}^{\infty}$`使其效用最大化，一阶条件写作：（其中`$w_t=\frac{W_t}{P_t}$`表示实际工资）

```math
-u_C(C_t,m_t,1-N_t) = \beta E_t\bigg[u_C(C_{t+1},M_{t+1},1-N_{t+1})\frac{1+i_t}{1+\pi_{t+1}}\bigg] = 0 \tag{6.8}

\frac{u_m(C_t,m_t,1-N_t)}{u_C(C_t,m_t,1-N_t)} = \frac{i_t}{1+i_t}\tag{6.9}

u_C(C_t,m_t,1-N_t)w_t = u_l(C_t,m_t,1-N_t)\tag{6.10}
```

* #### 厂商部门

* 厂商也连续分布在一个`$[0,1]$`的连续统上。同时为了方便讨论，厂商只使用劳动进行生产：

```math
Y_t(i) = A_tN_t^\alpha(i) \quad 0<\alpha<1,A>0 \tag{6.11}
```

* 每个厂商的需求函数就是前面的式(6.5)：

```math
Y_t(i) = C_t(i) = (\frac{P_t(i)}{P_t})^{-\theta}C_t \tag{6.12}
```

* 假设不存在跨期考虑，厂商只需要最大化其当期（实际）利润：（第二步需要代入式(6.11)和(6.12)）

```math
\begin{aligned}
\frac{Q_t(i)}{P_t} &= \frac{P_t(i)}{P_t}Y_t(i)-w_tN_t(i) \\
&= (\frac{P_t(i)}{P_t})^{1-\theta}C_t-\frac{w_t}{A_t^{\frac{1}{\alpha}}}(\frac{P_t(i)}{P_t})^{-\frac{\theta}{\alpha}}C_t^{\frac{1}{\alpha}} \tag{6.13}
\end{aligned}
```

* 由于上式中企业只能调整自己的价格`$P_t(i)$`，可以得到下面的一阶条件：

```math
(1-\theta)(\frac{P_t(i)}{P_t})^{-\theta}\frac{C_t}{P_t}+\frac{w_t}{A_t^{\frac{1}{\alpha}}}\frac{\theta}{\alpha}(\frac{P_t(i)}{P_t})^{-\frac{\theta}{\alpha}-1}\frac{C_t^{\frac{1}{\alpha}}}{P_t} = 0

\Leftrightarrow (1-\theta)Y_t(i)+\frac{\theta}{\alpha}w_tN_t(i)\frac{P_t}{P_t(i)} = 0

\Leftrightarrow (1-\theta)P_t(i)Y_t(i)+\frac{\theta}{\alpha}W_tN_t(i) = 0

\Leftrightarrow P_t(i) = \frac{\theta}{\theta-1}\frac{W_tN_t(i)}{\alpha Y_t(i)} \tag{6.14}
```

* 可以看出边际成本为`$MC_t(i)=\frac{W_tN_t(i)}{\alpha Y_t(i)}$`，且在垄断竞争下厂商设定的价格是高于其边际成本的，其比例`$\frac{\theta}{\theta-1}$`称为**价格加成**(markup)：（注意当`$\theta$`趋于无穷大时，商品趋于完美替代，此时价格加成趋于1）

```math
P_t(i) = \omega MC_t(i) \tag{6.15}
```

* 边际成本`$MC_t(i)$`就是下述最优化问题的拉格朗日乘子：

```math
\min_{N_t(i)}W_tN_t(i) \quad s.t. \quad Y_t(i) \geq \overline{Y} \tag{6.16}
```

* 其一阶条件为：

```math
W_t = MC_t(i)\alpha A_t(N_t(i))^{\alpha-1} \tag{6.17}

MC_t(i) = \frac{W_tN_t(i)}{\alpha Y_t(i)} \tag{6.18}
```

* #### 政府部门

* 此时政府部门货币发行仍然沿用MIU模型中的设定：

```math
\frac{M_t}{M_{t-1}}-\mu-1 = \theta_t

\theta_t = \rho^\theta \theta_{t-1}+\epsilon_t^\theta 
```

* #### 对称均衡

* 假定厂商是**同质的**，因此每个产品的价格和劳动力投入都是相等的，即`$P_t(i)=P_t,N_t(i)=N_t$`。

* 作为均衡的一般条件，还需要满足家庭部门的效用最大化和市场出清。

* #### 模型具体设定

* 我们给定一个具体的模型，其中效用函数是可加可分的：

```math
\begin{aligned}
u(C_t,1-N_t,m_t) &= \tilde{u}(C_t,N_t,m_t) \\
&= \frac{C_t^{1-\alpha}}{1-\sigma}+\phi_m\frac{m_t^{1-\chi}}{1-\chi}+\phi \frac{N_t^{1-\xi}}{1-\xi} \quad \sigma,\chi,\xi, \phi_m,\phi_n>0 \tag{6.19}
\end{aligned}
```

* 可以推出以下均衡条件：

```math
1 = \beta E_t\bigg[(\frac{C_{t+1}}{C_t})^{-\sigma}\frac{1+i_t}{1+\pi_{t+1}}\bigg] \tag{6.20}

Y_t = C_t \tag{6.21}

\phi_m m_t^{-\chi} = \frac{i_t}{1+i_t}C_t^{-\sigma} \tag{6.22}

\frac{m_t}{m_{t-1}}(1+\pi_t) = 1+\mu+\theta_t \tag{6.23}

\phi_n N_t^\xi = \frac{a}{\omega}\frac{Y_t}{N_t}C_t^{-\alpha} \tag{6.24}

Y_t = A_tN_t^\alpha \tag{6.25}
```

* 如果将式(6.20)（家庭最大化效用对商品的需求）和式(6.21)（市场可以提供的商品）合并，我们可以得到式(6.26)**（产品市场均衡）**，也就是**IS曲线**（`$Y_t$`和`$i_t$`具有负相关关系）。

```math
Y_t^{-\sigma} = \beta E_t\bigg[Y_{t+1}^{-\sigma}\frac{1+i_t}{1+\pi_{t+1}}\bigg] \tag{6.26}
```

* 同样地，将式(6.21)到(6.23)组合起来可以获得**货币市场的均衡**，也就是**LM曲线**（`$Y_t$`和`$i_t$`具有正相关关系）：

```math
\frac{i_t}{1+i_t} = \phi_m(\frac{1+\mu+\theta}{1+\pi_t}m_{t-1})^{-\chi}Y_t^\sigma \tag{6.27}
```

* 将式(6.24)和(6.25)组合起来可以得到**总供给曲线（AS曲线）**：

```math
Y_t = A_t^{\frac{1+\xi}{1+\xi+\alpha(\sigma-1)}}(\frac{\alpha}{\phi_n \omega})^{\frac{\alpha}{1+\xi+\alpha(\sigma-1)}} \tag{6.28}
```

* 将以上三个方程做对数线性化可得：（内生变量只有`$\hat{y}_t,\hat{i}_t,\hat{\pi}_t$`）

```math
\hat{y}_t = E_t\hat{y}_{t+1}-\frac{1}{\sigma}(\hat{i}_t-E_t\hat{\pi}_{t+1}) \tag{6.29}

\sigma \hat{y}_t = \chi(\theta_t+\hat{m}_{t-1}-\hat{\pi}_t)+\frac{\beta}{1+\mu-\beta}\hat{i}_t \tag{6.30}

\hat{y}_t = \frac{1+\xi}{1+\xi+\alpha(\sigma-1)}a_t \tag{6.31}
```

* 以上是为了在模型中引入垄断竞争而做的模型假设，但根据脉冲响应函数的结果，模型仍然表现出货币中性（产出只和`$A_t$`有关）。因此我们需要对模型修改以表现出**价格黏性**。

#### ==提前一期定价（过渡版本）==

* 假定企业在提前一期定价，同时承诺在下一期可以任意满足需求。因此企业需要最大化其对折现后下一期实际利润的期望：（`$\lambda_t$`是预算约束条件的拉格朗日乘子，也就是`$u_C(C_t,1-N_t,m_t)=(C_{t})^{-\sigma}$`，结合式(6.20)可知`$\beta\frac{\lambda_t}{\lambda_{t-1}}=\frac{1}{\frac{1+i_{t-1}}{1+\pi_t}}$`（实际利率的倒数），也就是说`$\beta\frac{\lambda_t}{\lambda_{t-1}}$`是 **“随机折现因子”**。

* 考虑这个最优化问题的目标函数，和式(6.13)过程相同：

```math
\begin{aligned}
E_{t-1}\bigg[\beta\frac{\lambda_t}{\lambda_{t-1}}\frac{Q_t(i)}{P_t}\bigg] &= E_{t-1}\bigg[\beta \frac{\lambda_t}{\lambda_{t-1}}\bigg(\frac{P_t(i)}{P_t}Y_t(i)-w_t N_t(i)\bigg)\bigg] \\
&= E_{t-1}\bigg[\beta\frac{\lambda_t}{\lambda_{t-1}}\bigg((\frac{P_t(i)}{P_t})^{1-\theta}C_t-\frac{w_t}{A_t^{\frac{1}{\alpha}}}(\frac{P_t(i)}{P_t})^{-\frac{\theta}{\alpha}}C_t^{\frac{1}{\alpha}}\bigg)\bigg] \tag{6.32}
\end{aligned}
```

* 求解这个最优化问题，同式(6.14)：（`$\beta$`和`$\lambda_{t-1}$`在t-1期是定值，因此在此可直接约去）

```math
E_{t-1}\bigg[\lambda_t\bigg((1-\theta)(\frac{P_t(i)}{P_t})^{-\theta}\frac{C_t}{P_t}+\frac{w_t}{A_t^{\frac{1}{\alpha}}}\frac{\theta}{\alpha}(\frac{P_t(i)}{P_t})^{-\frac{\theta}{\alpha}-1}\frac{C_t^{\frac{1}{\alpha}}}{P_t}\bigg)\bigg] = 0

\Leftrightarrow E_{t-1}\bigg[\lambda_t\bigg((1-\theta)(\frac{P_t(i)}{P_t})Y_t(i)+\frac{\theta}{\alpha}w_tN_t(i)\bigg)\bigg] = 0

\Leftrightarrow E_{t-1}\bigg[\lambda_t\bigg((\frac{P_t(i)}{P_t})Y_t(i)+\frac{\omega}{\alpha}W_tN_t(i)\bigg)\bigg] = 0

\Leftrightarrow E_{t-1}\bigg[\lambda_tY_t(i)\bigg(\frac{P_t(i)}{P_t}-\omega\frac{MC_t(i)}{P_t(i)}\bigg)\bigg] = 0 \tag{6.33}
```

* 在一个对称均衡中，我们可以得到市场总体产出和价格的关系：

```math
E_{t-1}\bigg[Y_t^{1-\sigma}\bigg(1-\omega\frac{MC_t}{P_t}\bigg)\bigg] = 0

\Leftrightarrow E_{t-1}\bigg[\bigg(Y_t^{1-\sigma}-\frac{\omega\phi_n}{\alpha}N_t^{1+\xi}\bigg)\bigg] = 0

\Leftrightarrow E_{t-1}\bigg[Y_t^{1-\sigma}-\frac{\omega\phi_n}{\alpha}A_t^{-\frac{1+\xi}{\alpha}}Y_t^{\frac{1+\xi}{\alpha}}\bigg] = 0

\Leftrightarrow E_{t-1}\bigg[Y_t^{1-\sigma}\bigg(1-A_t^{-\frac{1+\xi}{\alpha}}Y_t^{\frac{1+\xi+a(\sigma-1)}{\alpha}}\bigg)\bigg] = 0

\Leftrightarrow E_{t-1}\bigg[Y_t^{1-\sigma}\bigg(1-(\frac{Y_t}{Y_t^n})^{\frac{1+\xi+a(\sigma-1)}{\alpha}}\bigg)\bigg] = 0 \tag{6.34}
```

* 其中`$Y_t^n$`定义为自然产出 **（没有黏性价格的产出）**，也就是式(6.28)给出的产出：

```math
Y_t = A_t^{\frac{1+\xi}{1+\xi+\alpha(\sigma-1)}}(\frac{\alpha}{\phi_n \omega})^{\frac{\alpha}{1+\xi+\alpha(\sigma-1)}} \tag{6.35}
```

* 两者的关系经对数线性化后可写作：（可以看出式(6.37)
和式(6.31)是一样的，`$\hat{y}_t-\hat{y}_t^n$`就是**产出缺口**）

```math
E_{t-1}\bigg[\hat{y}_t-\hat{y}_t^n\bigg] = 0 \tag{6.36}

\hat{y}_t^n = \frac{1+\xi}{1+\xi+\alpha(\sigma-1)}a_t \tag{6.37}
```

* 在提前一期定价的模型设定下，货币供给增加确实会导致产出增加。但当货币供给持久性增加时，产出只会增加一期。因此我们需要引入提前多期定价来解决这个问题。

#### ==Calvo定价==

> 参考资料：
> Woodford(2003) CH.3 ; Walsh(2003) CH.5.4.

* 假设在每一期有`$\gamma$`的比例的厂商不改变其价格（沿用上一期的价格`$P_{t-1}$`），`$1-\gamma$`的厂商改变其商品价格（记新选择的价格为`$P_t^*$`）。

* 在引入Calvo定价之后，需要对Dixit-Stiglitz价格指数重新考察：

```math
\begin{aligned}
P_t &= \bigg((1-\gamma)(P_t^*)^{1-\theta}+\gamma\int_0^1P_{t-1}(i)^{1-\theta}di\bigg)^{\frac{1}{1-\theta}} \\
&= \bigg((1-\gamma)(P_t^*)^{1-\theta}+\gamma P_{t-1}^{1-\theta}\bigg)^{\frac{1}{1-\theta}} \tag{6.38}
\end{aligned}
```

* 此时企业要最大化其未来无穷期（不能改变价格的情况折现后）的利润，目标函数写作：（对照式(6.32)）

```math
E_{t}\bigg[\sum_{s=t}^\infty(\beta\gamma)^{s-t} \frac{\lambda_s}{\lambda_{t}}\bigg(\frac{P_t(i)}{P_s}Y_s(i)-w_s N_s(i)\bigg)\bigg] \\
= E_{t}\bigg[\sum_{s=t}^\infty(\beta\gamma)^{s-t} \frac{\lambda_s}{\lambda_{t}}\bigg((\frac{P_t(i)}{P_s})^{1-\theta}C_s-\frac{w_s}{A_s^{\frac{1}{\alpha}}}(\frac{P_t(i)}{P_s})^{-\frac{\theta}{\alpha}}C_s^{\frac{1}{\alpha}}\bigg)\bigg] \tag{6.39}
```

* 一阶条件写作：（对照式(6.33)）

```math
E_{t}\bigg[\sum_{s=t}^\infty(\beta\gamma)^{s-t} \frac{\lambda_s}{\lambda_{t}}\bigg((1-\theta)(\frac{P_t(i)}{P_s})^{-\theta}Y_s(i)+\frac{\theta}{\alpha}w_sN_s(i)\bigg)\bigg] = 0

\Leftrightarrow E_{t}\bigg[\sum_{s=t}^\infty(\beta\gamma)^{s-t} \frac{\lambda_s}{\lambda_{t}}Y_s(i)\bigg(\frac{P_t(i)}{P_s}-\omega\frac{MC_s(i)}{P_s(i)}\bigg)\bigg] = 0 \tag{6.40}
```

* 在对称均衡的情况下，有`$P_t(i)=P_t^*$`，上式写作：（`$\lambda_t$`被约掉了）

```math
E_{t}\bigg[\sum_{s=t}^\infty(\beta\gamma)^{s-t} \lambda_sY_s^t\bigg(\frac{P_t^*}{P_s}-\omega\frac{MC_s^t}{P_s}\bigg)\bigg] = 0 \tag{6.41}
```

* 其中`$Y_s^t=(\frac{P_t^*}{P_s})^{-\theta}Y_s$`，代表在t期最后一次设定价格的企业在s期的产出，且：

```math
Y_s = \bigg(\int_0^1Y_t(i)^{\frac{\theta-1}{\theta}}di\bigg)^{\frac{\theta}{\theta-1}} \tag{6.42}
```

* 同理，`$MC_s^t$`代表在t期最后一次设定价格的企业在s期的边际成本：

```math
MC_s^t = \frac{W_s}{\alpha}A_s^{-\frac{1}{\alpha}}(Y_s^t)^{\frac{1-\alpha}{\alpha}} \tag{6.43}
```

* 将`$Y_s^t=(\frac{P_t^*}{P_s})^{-\theta}Y_s$`代入式(6.41)，可得：

```math
E_{t}\bigg[\sum_{s=t}^\infty(\beta\gamma)^{s-t} \lambda_s(\frac{P_t^*}{P_s})^{-\theta}Y_s^{1-\sigma}\bigg(\frac{P_t^*}{P_s}-\omega\frac{MC_s^t}{P_s}\bigg)\bigg] = 0

\Leftrightarrow E_{t}\bigg[\sum_{s=t}^\infty(\beta\gamma)^{s-t} \lambda_s(\frac{P_t^*}{P_s})^{1-\theta}Y_s^{1-\sigma}\bigg] = E_{t}\bigg[\sum_{s=t}^\infty(\beta\gamma)^{s-t} \lambda_s(\frac{P_t^*}{P_s})^{-\theta}Y_s^{1-\sigma}\omega\frac{MC_s^t}{P_s}\bigg]

\Leftrightarrow E_{t}\bigg[\sum_{s=t}^\infty(\beta\gamma)^{s-t} \lambda_s(\frac{P_t}{P_s})^{1-\theta}Y_s^{1-\sigma}\bigg] = \frac{P_t}{P_t^*}\omega E_{t}\bigg[\sum_{s=t}^\infty(\beta\gamma)^{s-t} \lambda_s(\frac{P_t}{P_s})^{-\theta}Y_s^{1-\sigma}\omega\frac{MC_s^t}{P_s}\bigg] \tag{6.44}
```

* 可以解出每次重新定价的价格选择函数：

```math
\frac{P_t^*}{P_t} = 
\omega\frac{E_{t}\bigg[\sum_{s=t}^\infty(\beta\gamma)^{s-t} \lambda_s(\frac{P_s}{P_t})^{\theta}Y_s^{1-\sigma}\omega\frac{MC_s^t}{P_s}\bigg]}{E_{t}\bigg[\sum_{s=t}^\infty(\beta\gamma)^{s-t} \lambda_s(\frac{P_s}{P_t})^{1-\theta}Y_s^{1-\sigma}\bigg]} \tag{6.45}
```

* #### 菲利普斯曲线

* 定义`$p_t^*=\frac{P_t^*}{P_t}$`，`$mc_s^t=\frac{MC_s^t}{P_s}$`，此时可以将式(6.45)写作：(注意`$\frac{P_s}{P_t}=\frac{P_s}{P_{s-1}}\frac{P_{s-1}}{P_{s-2}}\cdots\frac{P_{t+1}}{P_t}=(1+\pi_s)(1+\pi_{s-1})\cdots(1+\pi_t)$`）

```math
p_t^* = \frac{\omega}{1+\pi_t}\frac{E_{t}\bigg[\sum_{s=t}^\infty(\beta\gamma)^{s-t} \bigg(\prod_{j=0}^{s-j}(1+\pi_{t+j})\bigg)^{\theta}Y_s^{1-\sigma}mc_s^t\bigg]}{E_{t}\bigg[\sum_{s=t}^\infty(\beta\gamma)^{s-t} \bigg(\prod_{j=0}^{s-j}(1+\pi_{t+j})\bigg)^{\theta-1}Y_s^{1-\sigma}\bigg]} \tag{6.46}
```

* 对分子执行对数线性化：

```math
(1-\beta \gamma) E_{t}\left[\sum_{s=t}^{\infty}(\beta \gamma)^{s-t}\left((1-\sigma) \hat{y}_{s}+\hat{m} c_{s}^{t}+\sum_{j=0}^{s-t} \theta \hat{\pi}_{t+j}\right)\right]
```

* 对分母执行对数线性化：

```math
(1-\beta \gamma) E_{t}\left[\sum_{s=t}^{\infty}(\beta \gamma)^{s-t}\left((1-\sigma) \hat{y}_{s}+\sum_{j=0}^{s-t}(\theta-1) \hat{\pi}_{t+j}\right)\right]
```

* 综上，式(6.45)的对数线性化为：

```math
\begin{aligned}
\hat{p}_{t}^{*}+\hat{\pi}_{t}=&(1-\beta \gamma) E_{t}\left[\sum_{s=t}^{\infty}(\beta \gamma)^{s-t}\left(\hat{mc}_{s}^{t}+\sum_{j=0}^{s-t} \hat{\pi}_{t+j}\right)\right] \\
=&(1-\beta \gamma)\left(\hat{m} c_{s}^{t}+\hat{\pi}_{t}\right)+\\
& \beta \gamma\left(\hat{\pi}_{t}+(1-\beta \gamma) E_{t}\left[\sum_{s=t+1}^{\infty}(\beta \gamma)^{s-t-1}\left(\hat{mc} _{s}^{t}+\sum_{j=0}^{s-t-1} \hat{\pi}_{t+1+j}\right)\right]\right)
\end{aligned} \tag{6.47}
```

* 可以据此求出`$\hat{mc}_{s}^{t}$`的表达式：（其中`$mc_s=\frac{w_s}{\alpha}A_s^{-\frac{1}{\alpha}}Y_s^{\frac{1-\alpha}{\alpha}}$`）

```math
\begin{aligned}
\hat{mc}_{s}^{t} &=\frac{w_{s}}{\alpha} A_{s}^{-\frac{1}{\alpha}} Y_{s}^{t} \frac{1-\alpha}{\alpha}=\frac{w_{s}}{\alpha} A_{s}^{-\frac{1}{\alpha}}\left(\frac{P_{t}^{*}}{P_{s}}\right)^{-\theta \frac{1-\alpha}{\alpha}} Y_{s}^{t \frac{1-\alpha}{\alpha}} \\
&=m c_{s}\left(p_{t}^{*}\right)^{-\theta \frac{1-\alpha}{\alpha}}\left(\frac{\pi_{j=0}^{s-t}\left(1+\pi_{t+j}\right)}{1+\pi_{t}}\right)^{\theta \frac{1-\alpha}{\alpha}}
\end{aligned} \tag{6.48}
```

* 再将式(6.49)代入式(6.47)可得：

```math
\begin{aligned}
\hat{p}_{t}^{*}+\hat{\pi}_{t}=&(1-\beta \gamma)\left(\hat{mc}_{t}-\theta \frac{1-\alpha}{\alpha} \hat{p}_{t}^{*}+\hat{\pi}\right)+\beta \gamma \hat{\pi}_{t} \\
+& \beta \gamma(1-\beta \gamma) E_{t}\left[\sum_{s=t+1}^{\infty}(\beta \gamma)^{s-t-1}\right.\\
&\left.\left(\hat{mc}_{s}-\theta \frac{1-\alpha}{\alpha} \hat{p}_{t}^{*}+\left(1+\theta \frac{1-\alpha}{\alpha}\right)^{s-t-1} \hat{\pi}_{t+1+j}\right)\right]
\end{aligned}

\begin{aligned}
\left(1+\theta \frac{1-\alpha}{\alpha}\right) \hat{p}_{t}^{*}=&(1-\beta \gamma) \hat{m}_{t}+\\
& \beta \gamma(1-\beta \gamma) E_{t}\left[\sum_{s=t+1}^{\infty}(\beta \gamma)^{s-t-1}\left(\hat{mc}_{s}+\left(1+\theta \frac{1-\alpha}{\alpha}\right)^{s-t-1} \hat{\pi}_{t+1+j}\right)\right]
\end{aligned}

\left(1+\theta \frac{1-\alpha}{\alpha}\right) \hat{p}_{t}^{*}=(1-\beta \gamma) \hat{mc}_{t}+\beta \gamma\left(1+\theta \frac{1-\alpha}{\alpha}\right) E_{t}\left[\hat{p}_{t+1}^{*}-\hat{\pi}_{t+1}\right] \tag{6.50}
```

* 最后一步的推导：

```math
\hat{m} c_{s}^{t+1} = \hat{mc}_{s}-\theta \frac{1-\alpha}{\alpha} \hat{p}_{t+1}^{*} \\
+\theta \frac{1-\alpha}{\alpha}\left(\sum_{j=0}^{s-t-1} \hat{\pi}_{t+1+j}\right)-\theta \frac{1-\alpha}{\alpha} \hat{\pi}_{t+1}

\begin{aligned}
\Leftrightarrow \hat{mc}_{s}^{t+1}+\sum_{j=0}^{s-t-1} \hat{\pi}_{t+1+j}=& \hat{mc}_{s}-\theta \frac{1-\alpha}{\alpha} \hat{p}_{t+1}^{*} \\
&+\left(1+\theta \frac{1-\alpha}{\alpha}\right)\left(\sum_{j=0}^{s-t-1} \hat{\pi}_{t+1+j}\right)-\theta \frac{1-\alpha}{\alpha} \hat{\pi}_{t+1}
\end{aligned}

\Leftrightarrow \hat{mc}_{S}+\left(1+\theta \frac{1-\alpha}{\alpha}\right) \sum_{j=0}^{s-t-1} \hat{\pi}_{t+1+j}=\hat{m c_{s}}^{t+1}+\sum_{j=0}^{s-t-1} \hat{\pi}_{t+1+j}+\theta \frac{1-\alpha}{\alpha}\left(\hat{p}_{t+1}^{*}+\hat{\pi}_{t+1}\right)

\begin{aligned}
&\beta \gamma(1-\beta \gamma) E_{t}\left[\sum_{s=t+1}^{\infty}(\beta \gamma)^{s-t-1}\left(\hat{m} c_{s}+\left(1+\theta \frac{1-\alpha}{\alpha}\right)^{s-t-1} \hat{\pi}_{t+1+j}\right)\right]\\
&=\beta \gamma(1-\beta \gamma) E_{t}\left[\sum_{s=t+1}^{\infty}(\beta \gamma)^{s-t-1}\left(\hat{m} c_{s}^{t+1}+\sum_{j=0}^{s-t-1} \hat{\pi}_{t+1+j}+\theta \frac{1-\alpha}{\alpha}\left(\hat{p}_{t+1}^{*}+\hat{\pi}_{t+1}\right)\right)\right] \\
&=\beta \gamma\left(1+\theta \frac{1-\alpha}{\alpha}\right) E_{t}\left[\hat{p}_{t+1}^{*}-\hat{\pi}_{t+1}\right]
\end{aligned}
```

* 再考虑式(6.38)

```math
P_{t} =\left((1-\gamma)\left(P_{t}^{*}\right)^{1-\theta}+\gamma P_{t-1}^{1-\theta}\right)^{\frac{1}{1-\theta}} \tag{6.51}

\Leftrightarrow 1 = (1-\gamma)\left(p_{t}^{*}\right)^{1-\theta}+\gamma\left(1+\pi_{t}\right)^{\theta-1} \tag{6.52}
```

* 对式(6.52)执行对数线性化：

```math
\hat{p}_{t}^{*}=\frac{\gamma}{1-\gamma} \hat{\pi}_{t} \tag{6.53}
```

* 再将式(6.53)代入(6.50)，得到NK模型下的**菲利普斯曲线**：（通货膨胀来自于社会边际成本及其预期的变化）

```math
\hat{\pi}_{t}=\kappa \hat{m} c_{t}+\beta E_{t}\left[\hat{\pi}_{t+1}\right] \tag{6.54}

\kappa=\frac{(1-\beta \gamma)(1-\gamma)}{\gamma\left(1+\theta \frac{1-\alpha}{\alpha}\right)}>0, \quad 0<\beta<1
```

* 可以认为当期通货膨胀是未来无穷期通货膨胀的一个折现：

```math
\hat{\pi}_{t}=\kappa E_{t} \sum_{s=t}^{\infty} \beta^{s-t} \hat{mc}_{s} \tag{6.55}
```

* 从定义（式(6.24)）上看，平均边际成本可以写作：

```math
mc_{t}=\frac{\phi_{n} N_{t}^{\xi}}{Y_{t}^{-\sigma}} \frac{N_{t}}{\alpha Y_{t}} \tag{6.56}
```

* `$Y_t$`写作：

```math
\begin{aligned}
Y_{t} &=\left(\int_{0}^{1} Y_{t}(i)^{\frac{\theta-1}{\theta}} d i\right)^{\frac{\theta}{\theta-1}} \\
&=A_{t}\left(\int_{0}^{1} N_{t}(i)^{\alpha \frac{\theta-1}{\theta}} d i\right)^{\frac{\theta}{\theta-1}}
\end{aligned} \tag{6.57}
```

* 将两者对数线性化后得：

```math
\hat{y}_t=a_t+\alpha\hat{n}_t \tag{6.58}

\begin{aligned}
\hat{mc}_{t} &=(1+\xi) \hat{n}_{t}+(\sigma-1) \hat{y}_{t} \\
&=-\frac{1+\xi}{\alpha} \hat{a}_{t}+\left(\frac{1+\xi}{\alpha}+\sigma-1\right) \hat{y}_{t} \\
&=\frac{1+\xi+\alpha(\sigma-1)}{\alpha}\left(\hat{y}_{t}-\hat{y}_{t}^{n}\right)
\end{aligned} \tag{6.59}
```

* 再将式(6.54)中的边际成本替换为式(6.59)中的表达形式，可以得到标准形式的菲利普斯曲线（通货膨胀和产出缺口的关系）或总供给曲线：

```math
\hat{\pi}_{t}=\kappa \frac{1+\xi+\alpha(\sigma-1)}{\alpha}\left(\hat{y}_{t}-\hat{y}_{t}^{n}\right)+\beta E_{t}\left[\hat{\pi}_{t+1}\right] \tag{6.60}
```

#### ==泰勒规则==

* 为了在新凯恩斯主义模型中研究货币政策，我们引入泰勒规则：（`$\pi^*$`是货币政策目标，`$\epsilon^i$`是外部冲击）

```math
\frac{1+i_{t}}{1+\bar{i}}=\Psi\left(\frac{1+i_{t-1}}{1+\bar{i}}, \frac{1+\pi_{t}}{1+\pi^{*}}, \frac{Y_{t}}{Y_{t}^{n}}, e_{t}^{\epsilon^{i}}\right) \tag{6.61}
```

* 泰勒规则是央行对货币政策调整时参考的一种准则，综合考虑当前利率、通货膨胀、产出和外部冲击。

* 对数线性化形式的泰勒规则可以写作：（其中`$0\leq\rho_i\leq1$`刻画利率惯性，`$\rho_\pi,\rho_y\geq0$`衡量央行对通胀和产出的反应程度）

```math
\hat{i}_{t}=\rho_{i} \hat{i}_{t-1}+\left(1-\rho_{i}\right)\left(\rho_{\pi} \hat{\pi}_{t}+\rho_{y} \hat{x}_{t}\right)+\epsilon_{t}^{i} \tag{6.62}
```

* **上式(式(6.62)）、IS曲线（式(6.29)）和Calvo定价下的总供给（AS）曲线（式(6.60)）构成了一个完备的NK模型，称为新凯恩斯主义的三个方程：**（其中`$\lambda=\frac{1+\xi+\alpha(\sigma-1)}{\alpha}$`，`$\hat{x}_t=\hat{y}_{t}-\hat{y}_{t}^{n}$`代表产出缺口，`$u_t=E_ty_{t+1}^n-y_t^n$`只依赖于供给冲击）

```math
\hat{x}_{t} = E_{t} \hat{x}_{t+1}-\frac{1}{\sigma}\left(\hat{i}_{t}-E_{t} \hat{\pi}_{t+1}\right)+u_{t} \tag{6.63}

\hat{\pi}_{t} =\kappa \lambda \hat{x}_{t}+\beta E_{t}\left[\hat{\pi}_{t+1}\right] \tag{6.64}
```

* 要使得这个非线性差分方程组有解，需要满足泰勒原则(Taylor Principle)，对系数有一定约束条件。

### ==7.贝叶斯估计==

* 近年来，研究者往往使用贝叶斯估计而非极大似然估计和校准

#### ==卡尔曼滤波和似然函数==

>参考：Technology Shocks in the New Keynesian Model(2003), Peter N. Ireland.（把fminu改成fminunc或者fminsearch）

* 将DSGE模型写作状态空间方程的形式：

* 状态变量：

```math
x_t = F(\theta)x_{t-1}+e_t \tag{7.1}

e_t = G(\theta)\epsilon_t \tag{7.2}

E(e_te'_t) = G(\theta)E(\epsilon_t\epsilon'_t)G(\theta)' = Q(\theta) \tag{7.3}
```

* 部分变量是无法观测的，因此有：

```math
X_t = H(\theta)'x_t \tag{7.4}

X_t = H(\theta)'x_t+u_t, \quad E(u_tu'_t)=\Sigma_u \tag{7.5}
```

* 卡尔曼滤波实际上是将计算联合概率的问题转换成了计算条件概率的问题。

* 假设所有变量都是可以观测的，即`$X_t \equiv x_t$`

* 预测方程如下：

```math
\hat{x}_t = F(\theta)x_{t-1}
```

* 在已知下一期真实值后，预测误差写作:

```math
\hat{e}_t = x_t-F(\theta)x_{t-1}
```

* （这一期的）条件似然函数满足：（`$X^{t-1}$`表示t-1期已知的值）

```math
L(X_t|X^{t-1}) = p_e(\hat{e}_t)
```

* 再将上述条件似然函数连乘，即可得到似然函数：

```math
L(X) = \prod_{t=1}^TL(X_t|X^{t-1})
```

* 对于无法完全观测的情况，需要先预测无法观测的值。

#### ==贝叶斯估计==

* 历史上难以应用贝叶斯估计的原因主要是因为难以求解后验分布的条件期望。数值积分方法由于参数维度过高难以模拟，因此一般使用MCMC方法。

* 先验分布的设定——类似于校准。对于影响稳态参数的参数，其先验分布往往使用样本信息；对于影响内部机制的参数，往往使用微观证据（如Calvo定价中的`$\gamma$`），比较难估计的参数可以直接校准。还需要考虑参数的值域。

* 寻找替代分布——**重要性采样**(importance sampling)，MCMC方法其实就是一种重要性采样的算法，意在寻找替代分布。

* 选择一个可以抽取的替代分布`$H(\Theta)$`和重要性权重`$\omega(\Theta)=\frac{P\left(\Theta \mid Y^{T}\right)}{h(\Theta)}$`：

```math
\begin{aligned}
E[g(\Theta)] &=\frac{\int g(\Theta) \frac{P\left(\Theta \mid Y^{T}\right)}{h(\Theta)} h(\Theta) d \Theta}{\int \frac{P\left(\Theta \mid Y^{T}\right)}{h(\Theta)} h(\Theta) d \Theta} \\
&=\frac{\int g(\Theta) \omega(\Theta) h(\Theta) d \Theta}{\int \omega(\Theta) h(\Theta) d \Theta}
\end{aligned}
```

* 据此可以实现近似积分符号：

```math
E[g(\Theta)] \approx \frac{\sum_{m=1}^{M} \omega\left(\Theta^{(m)}\right) g\left(\Theta^{(m)}\right)}{\sum_{m=1}^{M} \omega\left(\Theta^{(m)}\right)}
```

* 如何寻找替代分布？（通常选择厚尾的分布）一般来说Dynare会自动寻找替代分布。同时不能独立抽取（方差会越来越大，而且计算量很大）。

* #### M-H算法

* 核心概念：接受-拒绝法则，通过比较抽取的结果的概率，如果增大则保留，否则拒绝。

![M-H](https://note.youdao.com/yws/res/6875/14D95EF05B5D41269316CD32C0F0FCC3)

* 其中：

```math
q\left(\Theta_{i+1} \mid \Theta_{i}\right)=\min \left[1, \frac{P\left(\Theta_{i+1}^{*} \mid Y^{T}\right)}{P\left(\Theta_{i}^{*} \mid Y^{T}\right)} \frac{h\left(\Theta_{i}^{*} ; \psi\right)}{h\left(\Theta_{i+1}^{*} ; \psi\right)}\right]
```

* 更详细的技术细节见课件。

#### ==贝叶斯估计的Dynare实操==

```Matlab
var y c k i l y_l w r z;%内生变量
varexo e;%外生变量
parameters beta psi delta alpha rho epsilon;

model;
(1/c) = beta*(1/c(+1))*(1+r(+1)-delta);
psi*c/(1-l) = w;
c+i = y;
y = (k(-1)^alpha)*(exp(z)*l)^(1-alpha);
w = y*((epsilon-1)/epsilon)*(1-alpha)/l;
r = y*((epsilon-1)/epsilon)*alpha/k(-1);
i = k-(1-delta)*k(-1);
y_l = y/l;
z = rho*z(-1)+e;
end;

varobs y;%使用产出的数据来估算参数

initval;
  k = 9;
  c = 0.76;
  l = 0.3;
  w = 2.07;
  r = 0.03;
  z = 0;
  e = 0;
end;

estimated_params;%设定的先验分布形式，均值和方差
alpha, beta_pdf, 0.35, 0.02;
beta, beta_pdf, 0.99, 0.002;
delta, beta_pdf, 0.025, 0.003;
psi, gamma_pdf, 1.75, 0.1;
rho, beta_pdf, 0.95, 0.05;
epsilon, gamma_pdf, 10, 0.5;
stderr e, inv_gamma_pdf, 0.01, inf;
end;

estimation(datafile=simuldataRBC,nobs=200,first_obs=500,mh_replic=2000,mh_nblocks=2,mh_drop=0.45,mh_jscale=0.8,mode_compute=4);
%选项从左到右：数据来源，使用的观测样本数目，第一个使用的观测值序号（前面的舍去），（MH算法的）链长度（抽取多少次），（MH算法的）链数目，（MH算法的）接受拒绝率，（MH算法的）c值（见课件39页），计算众数的方法
```

* 返回结果：灰色的线是先验分布，黑色的线是后验分布，绿线代表后验分布众数（概率密度最大的点），还会返回先验均值和后验均值的对照表。