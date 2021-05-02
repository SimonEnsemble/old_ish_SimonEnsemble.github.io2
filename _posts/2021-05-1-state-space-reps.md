---
layout: post
title: state space representation of closed-loop systems with time delays
tags: [control theory]
snippet: how to derive the state space representation of closed-loop systems with time delays.
author: Cory Simon
---

we will learn to convert a closed-loop transfer function with time delays into a state space representation in the time domain. 
then, we will use `DifferentialEquations.jl` in Julia to numerically solve the resulting delay differential equations for the closed-loop dynamics.

# the feedback control system

consider the feedback control system depicted in the block diagram below.

<figure>
    <img src="/images/feedbackloop.png" alt="image" style="width: 100%;">
</figure>

the variables are:
* $Y(s)$: controlled variable
* $Y_{sp}(s)$: set point for controlled variable
* $Y_m(s)$: measurement of controlled variable
* $U(s)$: manipulated variable (output by the controller)
* $E(s)$: error signal

the transfer functions describe the dynamics of:
* $g_u(s)$: the process
* $g_c(s)$: the feedback control law
* $g_m(s)$: the sensor


# the closed-loop transfer function for the servo response

an implicit equation for the servo response (response of controlled variable $Y$ to set point changes $Y_{sp}$) is, from inspecting the block diagram:

$$
\begin{align}
Y(s)&=g_u(s)g_c(s)E(s)  \\
&= g_u(s)g_c(s)[Y_{sp}(s)-Y_m(s)]  \\
&= g_u(s)g_c(s)[Y_{sp}(s)-g_m(s)Y(s)].
\end{align}
$$

solving for $Y$, the closed-loop response to a set point change is:

$$Y(s)=\dfrac{g_u(s)g_c(s)}{1+g_u(s)g_c(s)g_m(s)}Y_{sp}(s)$$

# the state space representation

in the closed-loop servo response, the output is $y(t)$ and the input is $y_{sp}(t)$. we wish to convert the closed-loop transfer function governing the servo response to a state space representation in the time domain:

$$
\begin{align}
\frac{d\mathbf{x}(t)}{dt}&=\mathbf{A}\mathbf{x}(t)+\mathbf{B}\mathbf{x}(t-\theta)+\mathbf{C}y_{sp}(t) \\
y(t)&=\mathbf{D}\mathbf{x}(t) + \mathbf{E} \mathbf{x}(t-\phi)
\end{align}
$$

here, $\mathbf{x}$ is the internal state vector and $\mathbf{A}$, $\mathbf{B}$, $\mathbf{C}$, $\mathbf{D}$, and $\mathbf{E}$ are constant matrices. the time delays here are $\theta$ and $\phi$. converting a transfer function into this form is usually required for using numerical differential equation solvers.


# learn by example

suppose:

* the process is first-order plus time delay,
$$
g_u(s)=\dfrac{3e^{-0.2s}}{5s+1}
$$
* we use a PI controller,
$$
g_c(s)=\dfrac{2s+1}{2s} 
$$
* the sensor has a time delay,
$$
g_m(s)=e^{-0.1s}
$$

### the closed-loop transfer function

the closed-loop transfer function governing the servo response is then:

$$
\begin{align}
Y(s)= \frac{\frac{3e^{-0.2s}}{5s+1}\frac{2s+1}{2s}}{1+\frac{3e^{-0.2s}}{5s+1}\frac{2s+1}{2s}e^{-0.1s}} Y_{sp}^*(s)
\end{align}
$$

simplifying this into a sort-of (ignoring the time delays) rational function:

$$
\begin{align}
Y(s)= \frac{e^{-0.2s}(6s+3)}{10s^2+2s+e^{-0.3s}(6s+3)} Y_{sp}(s)
\end{align}
$$

### naively converting into the time domain

rearranging the transfer function with the output $Y$ and input $Y_{sp}$ on opposite sides:

$$Y(s)[10s^2+2s+e^{-0.3s}(6s+3)]=Y_{sp}(s)[e^{-0.2s}(6s+3)]$$

we see the equivalent delay differential equation is:

$$10y^{\prime \prime }(t)+2y^{\prime }(t)+6y^{\prime }(t-0.3)+3y(t-0.3)=6y_{sp}^\prime(t-0.2)+3y_{sp}(t-0.2)$$

subject to initial conditions $y(0)=y^\prime(0)=0$.

this is not a state-space representation because we see higher than first-order derivatives in both (i) the output $y(t)$ and (ii) the input $y_{sp}(t)$. 
by appropriately defining a state vector $\mathbf{x}(t):=[y(t), y^\prime(t)]$, we can solve problem (i), but problem (ii) requires a different approach.

### trick: introducing an intermediate variable

to handle the derivative of the input $y_{sp}^\prime(t)$, we split the transfer function into two processes in series, introducing an intermediate variable $V(s)$:

$$
\begin{align}
V(s)&:=\frac{1}{10s^2+2s+e^{-0.3s}(6s+3)}Y_{sp}(s)\\
Y(s)&=e^{-0.2s}(6s+3)V(s)
\end{align}
$$

this allows to write a state space representation for the intermediate variable $v(t)$, then write $y(t)$ in terms of $v(t)$ and its derivative. because the numerator is one, we won't see any derivatives of the input $y_{sp}(t)$. :+1:

rearranging to clearly see the derivatives involved:

$$
\begin{align}
V(s)[10s^2+2s+e^{-0.3s}(6s+3)]&=Y_{sp}(s)\\
Y(s)&=e^{-0.2s}(6s+3)V(s)
\end{align}
$$

### converting to state space representation

the equivalent delay differential equations in the time domain is:

$$
\begin{align}
&10v^{\prime\prime}(t)+2v^\prime(t)+6v^\prime(t-0.3)+3v(t-0.3)=y_{sp}(t) \\
&y(t)=6v^\prime(t-0.2)+3v(t-0.2)
\end{align}
$$

voila, we've eliminated derivatives of the input $y_{sp}(t)$ by introducing this intermediate variable $v(t)$! 

the first equation, given the input $y_{sp}(t)$, can be solved for $v(t)$ independently of the second equation. the second equation then expresses the output $y(t)$ in terms of the intermediate variable $v(t)$.

however, we have higher than first-order derivaties of $v(t)$. to fix that, define the state vector $\mathbf{x}(t):=[x_1(t), x_2(t)]$ such that:

$$
\begin{align}
x_1(t)&:=v(t)\\
x_2(t)&:=v^\prime(t)
\end{align}
$$

then:

$$
\begin{align}
x_1^\prime(t)=&x_2(t)\\
x_2^\prime(t)=&-\frac{1}{5}x_2(t)-\frac{3}{5}x_2(t-0.3)-\frac{3}{10}x_1(t-0.3) +\frac{1}{10}y_{sp}(t)
\end{align}
$$

the first equation follows from the definition of $x_1$ and $x_2$. 
the second equation follows from the differential equation.

and,

$$y(t)=6x_2(t-0.2)+3x(t-0.2)$$

finally, writing this in matrix form, we arrive at the state space representation:

$$
\begin{align}
&\frac{d}{dt}\begin{bmatrix}x_1(t)\\x_2(t)\end{bmatrix}=\begin{bmatrix}0 & 1 \\ 0 &-\frac{1}{5}\end{bmatrix}\begin{bmatrix}x_1(t)\\x_2(t)\end{bmatrix}
+\begin{bmatrix} 0 & 0 \\ -\frac{3}{10} &-\frac{3}{5}\end{bmatrix}\begin{bmatrix}x_1(t-0.3)\\x_2(t-0.3)\end{bmatrix} + \begin{bmatrix}0 \\ \frac{1}{10}\end{bmatrix}y_{sp}(t)\\
&y(t)=\begin{bmatrix}3 & 6\end{bmatrix}\begin{bmatrix}x_1(t-2)\\x_2(t-2)\end{bmatrix}
\end{align}
$$

subject to initial condition:

$$
\begin{bmatrix}x_1(0)\\x_2(0)\end{bmatrix} = \begin{bmatrix}0\\0\end{bmatrix}
$$

the state space representation is important because numerical solvers usually require the ODE system to be written in this form.

### numerical solution in `DifferentialEquations.jl`

we will use `DifferentialEquations.jl` to numerically solve the system of delay differential equations above subject to a unit step input (a unit set point change). see the [documentation on solving delay differential equations](https://diffeq.sciml.ai/stable/tutorials/dde_example/) to understand more fully.

```julia
using DifferentialEquations, PyPlot

# define matrices
A = [0    1;
     0  -1/5]
B = [  0      0; 
     -3/10  -3/5]
C = [0, 1/10]
E = [3 6]

# history function.
#   stores history of solution to handle time delay
#   initiate as v, v′ = 0 for t < 0
h(params, t) = zeros(2)

# the input
y_sp(t) = t > 0.0 ? 1.0 : 0.0

# right-hand side of ODE in x
#   h(t-0.3) serves as x(t-0.3), looking up history of solution.
f(x, h, params, t) = A * x + B * h(params, t - 0.3) + C * y_sp(t)

# initial condition
x₀ = [0.0, 0.0]

# numerical solution
final_time = 25.0
prob = DDEProblem(f, x₀, h, (0.0, final_time))
sol = solve(prob, d_discontinuities=[0.0, 0.3])

# extracting solution at an array of times
t = range(0.0, final_time, length=300)
x = sol.(t) # internal state vectors
t = t .+ 0.2 # shift time
y = [(E * x_i)[1] for x_i in x]

# plot
figure()
plot(t, y_sp.(t), color="C1", label=L"$y_{sp}(t)$", linestyle="--")
plot(t, y, color="C0", label=L"$y(t)$")
plot([-1.0, 0.2], [0, 0], color="C0")
plot([-1.0, 0.2], [0, 0], color="C1", linestyle="--")
xlabel(L"$t$")
ylabel(L"$y(t)$")
legend()
```

voila, we see the closed-loop response, under PI control, to the unit set point change.
<figure>
    <img src="/images/time_delay_response.png" alt="image" style="width: 100%;">
</figure>

