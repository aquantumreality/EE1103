# ODE

Landau–Lifshitz–Gilbert equation (LLG equation) describes the magnetic precession in solid.
Two forms of LLG equation were given, one in terms of m and another in terms of θ, ϕ.
We used the latter to solve LLG equation:

dθ/dt = ( γα/(α2+1) ) Hzsinθ, −dφ/dt = θ'αsinθ = (dθ/dt) (1/αsinθ)

Where γ = -1.76e11, α can be varied from 0 to 1,
and Hz = (1/μ0)100e-03 A/m.

This ODE can be solved by Midpoint method, third- and fourth-order Runge-Kutta method (RK3, RK4) and Runge–Kutta–Fehlberg method (RK45) etc.

The aims of this assignment are as follows:
 1. Solve the ODE with RK45 and any other method of your choice and plot error vs no. of iteration curve for a specific step-size.
 2. Vary the step-size and record root-mean-square error for each step-size. Plot step-size vs root-mean-square error.
 3. Vary α and record the switching time Tsw for each α. Plot α vs Tsw.
 4. Plot reduced magnetization vector (m→) in three dimensions.

Take intial values of polar and azimuthal angles as θin = 179∘ and ϕin = 1∘.
