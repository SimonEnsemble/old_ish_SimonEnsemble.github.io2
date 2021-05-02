### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# ╔═╡ b926bbf4-ab68-11eb-136c-952f90c757df
using DifferentialEquations, PyPlot

# ╔═╡ ab5e06ed-f1f3-41a3-99e2-2affda673af8
PyPlot.matplotlib.style.use("Solarize_Light2")

# ╔═╡ 1b8c79c1-d663-42dc-beda-29a881dcbd2a
begin
	# define matrices
	A = [0    1;
		 0  -1/5]
	B = [  0      0; 
		 -3/10  -3/5]
	C = [0, 1/10]
	E = [3 6]
end

# ╔═╡ 3f443782-32d8-44f0-baf3-41645d949573
y_sp(t) = t > 0.0 ? 1.0 : 0.0

# ╔═╡ b9aa78a5-d7bb-4bd8-a38d-46ef023670ca
# right-hand side of ODE in x
f(x, h, params, t) = A * x + B * h(params, t - 0.3) + C * y_sp(t)

# ╔═╡ c8c53be9-39ed-4354-90ba-27ae3f9648a6
# initial condition
x₀ = [0.0, 0.0]

# ╔═╡ 874519d1-b28d-4749-89e8-fa32fa901a15
# history function.
#   stores history of solution to handle time delay
#   initiate as v, v′ = 0 for t < 0
h(params, t) = zeros(2)

# ╔═╡ 805a0a12-3780-49b7-8296-ca76da5ea247
begin
	final_time = 25.0
	prob = DDEProblem(f, x₀, h, (0.0, final_time))
	sol = solve(prob, d_discontinuities=[0.0, 0.3])
end

# ╔═╡ e2e938e4-6cf4-4cbb-95b4-f56beeea3048
begin
	# get solution at an array of times to plot it
	t = range(0.0, final_time, length=300)
	x = sol.(t) # internal state vectors
	# shift time
	t = t .+ 0.2
	# determine y at these times
	y = [(E * x_i)[1] for x_i in x] # y
	
	figure()
	plot(t, y_sp.(t), color="C1", label=L"$y_{sp}(t)$", linestyle="--")
	plot(t, y, color="C0", label=L"$y(t)$")
	plot([-1.0, 0.2], [0, 0], color="C0")
	plot([-1.0, 0.2], [0, 0], color="C1", linestyle="--")
	xlabel(L"$t$")
	ylabel(L"$y(t)$")
	legend()
	tight_layout()
	savefig("../images/time_delay_response.png", format="png", dpi=300)
	gcf()
end

# ╔═╡ 117c23bd-a09b-43f1-a0b1-1a031246c1ca
y

# ╔═╡ Cell order:
# ╠═b926bbf4-ab68-11eb-136c-952f90c757df
# ╠═ab5e06ed-f1f3-41a3-99e2-2affda673af8
# ╠═1b8c79c1-d663-42dc-beda-29a881dcbd2a
# ╠═3f443782-32d8-44f0-baf3-41645d949573
# ╠═b9aa78a5-d7bb-4bd8-a38d-46ef023670ca
# ╠═c8c53be9-39ed-4354-90ba-27ae3f9648a6
# ╠═874519d1-b28d-4749-89e8-fa32fa901a15
# ╠═805a0a12-3780-49b7-8296-ca76da5ea247
# ╠═e2e938e4-6cf4-4cbb-95b4-f56beeea3048
# ╠═117c23bd-a09b-43f1-a0b1-1a031246c1ca
