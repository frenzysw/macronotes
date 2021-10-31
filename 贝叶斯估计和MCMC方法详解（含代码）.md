## 贝叶斯估计和MCMC方法详解（含代码）
>主要参考书籍：**Blake and Mumtaz (2017):** ***Applied Bayesian Econometrics for Central Bankers Updated***

### 1.1 线性回归（自回归）的贝叶斯估计
对于线性回归（自回归）的贝叶斯估计，先验分布的选取主要考虑系数和误差方差：

1. 系数未知，误差方差已知，此时似然函数是一个方差已知，均值未知的正态分布，先验分布应选取正态分布
2. 系数已知，误差方差未知，此时似然函数是一个均值已知，方差未知的正态分布，先验分布应选取倒伽马分布（方差倒数服从伽马分布）
3. 系数和误差方差均未知，此时似然函数是一个均值和方差均未知的正态分布，**先验分布应选取以上两者的乘积，但此时系数的先验分布需以误差方差为条件（方差仍为无条件分布）**

* 在情况3中，最后得到的是一个系数和误差方差的联合后验分布，**需要将其积分为分别对系数和误差方差的边缘后验分布以进行进一步的推断（进行统计推断需要对边缘后验分布的形式进行，如`$f(B|Y_t)$`）**。在线性回归的“**自然共轭先验**”（即与似然函数同形式的先验分布）情况下，这样的积分是有解析解的。系数的边缘后验分布是一个多元t分布，而误差方差（的倒数）的后验分布是一个伽马分布。但是，当先验分布发生变化时（如对系数的先验分布不以误差方差为条件），这样的解析解便无法导出。积分求解析解的困难导致了吉布斯采样的出现和流行。

* 情况1和2只有一个未知变量，不需要进行吉布斯采样，而情况3有两个变量，可能面临难以积分的问题，需要使用吉布斯采样。

* 吉布斯采样最后获得的是每一个变量的**边缘概率密度**（由于用于求解的条件概率密度是以观测值`$Y$`为条件的，最终的样本往往也以`$Y$`为条件）`$f(x_i)$`，如果使用常规方法则需要使用积分的方式求解。

![吉布斯采样过程](https://note.youdao.com/yws/res/2962/CB6850B5C6804F2497BF6685C6C21C83)

* 情况3的吉布斯采样的操作方法中需要两个条件后验分布，即`$f(B\mid Y_t,\cfrac{1}{\sigma^2})$`与`$f(\cfrac{1}{\sigma^2}|Y_t,B)$`，而这两个分布可从1、2的情况中获取（通过假设部分参数已知）。

* 情况1下对系数的迭代方法：
![系数后验估计迭代](https://note.youdao.com/yws/res/530/7794AFFA3ED249FB9555E4C412E85D1A)

* 情况2下对误差方差的迭代方法：
![误差方差后验估计迭代](https://note.youdao.com/yws/res/541/FC0FE0EC90604817AABC877EA245D6E4)

* 进行**模型比较**的常用方法是通过获得的样本进行相应的点估计和区间估计，但也可以通过边缘似然函数`$F(Y)$`进行：通过比较数据集`$Y$`在不同参数（即`$B$`和`$σ^2$`）的模型中的似然函数来进行统计推断（函数值越大说明模型越好），缺点是计算上比较复杂（很多时候没法积分，还需要另外使用吉布斯采样）且当先验分布不理想时这样的判定方法可能出错。
```math
F(Y)=\int {F(Y\mid B,\sigma^2)p(B,\sigma^2)d\Xi}
```

* 具体代码如下（其中中文注释是我写的，下同）：
```Matlab
clear
addpath('functions'); %this line adds functions to take lags etc
%an AR 2 model for US inflation
%load inflation data
Y=xlsread('inflation.xls');
T=rows(Y);
X=[ones(T,1) lag0(Y,1) lag0(Y,2)];

%remove missing obs
Y=Y(3:end);
X=X(3:end,:);
T=rows(X);

%step 1 set priors and starting values
%priors for B
B0=[0;0;0];
Sigma0=eye(3);
%priors for sigma2
T0=1;
D0=0.1;

%starting values
B=B0;
sigma2=1;
reps=25000;   %total numbers of Gibbs iterations
burn=15000;   %percent of burn-in iterations
out1=[];
out2=[];
for i=1:reps                         %注意每一次参数的更新都分解为了包含初始值的参数的部分（如sigma0）和包含上一次迭代的参数的部分（如sigma2）
    %step 2 Sample B conditional on sigma N(M*,V*)
   M=inv(inv(Sigma0)+(1/sigma2)*(X'*X))*(inv(Sigma0)*B0+(1/sigma2)*X'*Y);
    V=inv(inv(Sigma0)+(1/sigma2)*(X'*X));
    chck=-1;
    while chck<0                     %check for stability
        B=M+(randn(1,3)*chol(V))';   %这里的Cholesky分解对应公式中的V1/2
        b=[B(2) B(3);1   0];         %留意一下这种检查单位根的方法（伴随式）
        ee=max(abs(eig(b)));
        if ee<=1
            chck=1;
        end
    end
    
    %step 3 sample sigma2 conditional on B from IG(T1,D1);
    %compute residuals
    resids=Y-X*B;                    %注意这里是X*B，如果将观测值以行排列的话就该使用该形式（书上一开始的式子错了）

    %compute posterior df and scale matrix
    T1=T0+T;
    D1=D0+resids'*resids;
    %draw from IG
    z0=randn(T1,1);
    z0z0=z0'*z0;
    sigma2=D1/z0z0;                  %这里有一个十分值得注意的地方，如果使用MATLAB自带的伽马函数抽取器（gamrnd），要注意其参数并非(T1/2,D1/2)，而是(T1/2,2/D1)
                                     %因为当X~G(X,Y)时，1/X~IG(X,1/Y)，因此语句应为sigma2=1/gamrnd(T1/2,2/D1);
    if i>burn
        out1=[out1;B'];
        out2=[out2;sigma2];
    end
end

%plot marginal posterior distributions
subplot(2,2,1);
hist(out1(:,1),50);
axis tight
title('Constant');
subplot(2,2,2);
hist(out1(:,2),50);
axis tight
title('AR(1) coefficient');
subplot(2,2,3);
hist(out1(:,3),50);
axis tight
title('AR(2) coefficient');
subplot(2,2,4);
hist(out2(:,1),50);
axis tight
title('\sigma^{2}');

%compute mean of the marginal posterior distribution of B
MB=mean(out1);
%compute standard error
VB=std(out1);
%compute 95% error band
EB=prctile(out1,[5 95]);
```
### 1.2 扰动项序列相关的自回归模型的贝叶斯估计
* 接下来考察如何使用吉布斯采样对存在扰动项序列相关的自回归模型进行贝叶斯估计，模型形式如下：
```math
Y_t=\alpha+B_1Y_{t-1}+B_2Y_{t-2}+v_t

v_t=\rho v_{t-1}+\epsilon_t, \epsilon_t\sim N(0,\sigma^2)
```
* 首先考虑ρ已知的情况，此时可以通过变量变换将原模型转化为一个无序列相关的自回归模型，具体如下：

![式3.15](https://i.loli.net/2020/09/09/CDPcs6vmd5epLox.png)


* 再考虑`$\alpha$`、`$B_1$`和`$B_2$`均已知的情况，通过变换`$v_t=Y_t-(\alpha+B_1 Y_{t-1} +B_2 Y_{t-2})$`后，只需要估计方差方程即可。

* 基于上述分析，估计的基本思路如下：
1. 抽取条件分布`$f(α,B_1,B_2|\sigma^2,ρ)$`（基于3.15式）
2. 抽取条件分布`$f(\rho|\sigma^2,\alpha,B_1,B_2)$`（基于变换`$v_t=Y_t-(\alpha+B_1 Y_{t-1} +B_2 Y_{t-2})$`）
3. 抽取条件分布`$f(\sigma^2|\rho,\alpha,B_1,B_2)$`

* 先验分布选定：
1. 条件分布`$f(α,B_1,B_2|\sigma^2,ρ)$`的先验分布是正态分布
2. 条件分布`$f(\rho|\sigma^2,\alpha,B_1,B_2)$`的先验分布是正态分布
3. 条件分布`$f(\sigma^2|\rho,\alpha,B_1,B_2)$`的先验分布是倒伽马分布

* 对三个条件分布的迭代：
1. 对条件分布`$f(α,B_1,B_2|\sigma^2,ρ)$`的迭代参考上一部分自回归中的情况1，需要注意的是在本例中进行了变量替换，即使用`$Y_{t}^{*}$`和`$Y_{t-1}^{*}$`、`$Y_{t-2}^{*}$`（作为自变量）代替上一部分情况1中的`$Y_{t}$`和`$X_{t}$`。
2. 对条件分布`$f(\rho|\sigma^2,\alpha,B_1,B_2)$`的迭代，经过变换`$v_t=Y_t-(\alpha+B_1 Y_{t-1} +B_2 Y_{t-2})$`后，原模型可表述为`$v_t=\rho v_{t-1}+\epsilon_t$`，此时使用`$v_t$`及其滞后项代替上一部分情况1中的`$Y_{t}$`和`$X_{t}$`（`$v_t$`怎么获得？通过`$Y_t-\alpha-B_1Y_{t-1}-B_2Y_{t-2}=v_t$`获得）。
3. 对条件分布`$f(\sigma^2|\rho,\alpha,B_1,B_2)$`的迭代参考上一部分的情况2，同样要注意`$Y_{t}$`和`$X_{t}$`的替换。

* 代码如下（注意`$Y_{t}$`和`$X_{t}$`的选取）：
```Matlab
clear
addpath('functions');
%an AR 2 model for US inflation with autocorrelated AR(1) disturbances
%load inflation data
Y=xlsread('\data\inflation.xls');
T=rows(Y);
X=[ones(T,1) lag0(Y,1) lag0(Y,2)];

%remove missing obs
Y=Y(3:end);
X=X(3:end,:);
T=rows(X);
%step 1 set priors and starting values
%priors for B
B0=[0;0;0];
Sigma0=eye(3);
%priors for sigma2
T0=1;
D0=0.1;
%priors for rho
rho0=0;
Sigma0r=1;

%starting values
B=B0;
rho=rho0;
sigma2=1;

reps=25000;
burn=24000;
out1=[]; %will save the inflation forecast
out2=[];
out3=[];
out4=[];
for i=1:reps
    
    %step 2 Sample B conditional on sigma N(M*,V*)
    %remove serial correlation
    ystar=Y-lag0(Y,1)*rho;  %ystar的生成
    xstar=X-lag0(X,1)*rho;
    ystar=ystar(2:end,:);
    xstar=xstar(2:end,:);
    M=inv(inv(Sigma0)+(1/sigma2)*(xstar'*xstar))*(inv(Sigma0)*B0+(1/sigma2)*xstar'*ystar);
    V=inv(inv(Sigma0)+(1/sigma2)*(xstar'*xstar));
    chck=-1;
    while chck<0
        B=M+(randn(1,3)*chol(V))';
        b=[B(2) B(3);1   0];
        ee=max(abs(eig(b)));
        if ee<=1
            chck=1;
        end
    end
    
    %step 3 compute rho
    y=Y-X*B;            %这里的 y就是模型中的 v（矩阵形式）
    x=lag0(y,1);        %自变量是 v的滞后项
    y=y(2:end);
    x=x(2:end);
    MM=inv(inv(Sigma0r)+(1/sigma2)*(x'*x))*(inv(Sigma0r)*rho0+(1/sigma2)*x'*y);
    VV=inv(inv(Sigma0r)+(1/sigma2)*(x'*x));
    %draw rho but again ensure stationarity
    chck=-1;
    while chck<0
        rho=MM+(randn(1,1)*chol(VV))';
        
        ee=abs(rho);    %这里可以直接求出系数 rho
        if ee<=1
            chck=1;
        end
    end
    
    
    
    %step 3 sample sigma2 conditional on B from IG(T1,D1);
    %compute residuals
    resids=ystar-xstar*B;
    %compute posterior df and scale matrix
    T1=T0+T;
    D1=D0+resids'*resids;
    %draw from IG
    z0=randn(T1,1);
    z0z0=z0'*z0;
    sigma2=D1/z0z0;
    
    if i>burn
        %compute forecast for 12 quarters
        yhat=zeros(14,1);
        vhat=zeros(14,1);
        yhat(1:2)=Y(end-1:end); %starting values
        cfactor=sqrt(sigma2);  %standard deviation of the shocks
        for m=3:14
            vhat(m)=vhat(m-1)*rho+randn(1,1)*cfactor;
            yhat(m)=[1 yhat(m-1) yhat(m-2)]*B+vhat(m);
            
        end
        %save
        out1=[out1 [Y;yhat(3:end)]];
        out2=[out2;B'];
        out3=[out3;rho];
        out4=[out4;sigma2];
    end
end


figure(1)
TT=1947.25:0.25:2012.5;
out2x=prctile(out1',[10 20 30 40 50 60 70 80 90])';
plot(TT,out2x);

xlim([2000 2013])

save output.mat out2 out3 out4
```


