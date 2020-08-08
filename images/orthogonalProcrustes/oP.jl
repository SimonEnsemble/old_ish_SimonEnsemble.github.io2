### A Pluto.jl notebook ###
# v0.11.3

using Markdown
using InteractiveUtils

# ╔═╡ ed7c391e-d937-11ea-31e0-d5884018ca29
md"# Orthogonal Procrustes"

# ╔═╡ fd1fbc70-d936-11ea-2712-ed19c93c0c9f
begin
	using DelimitedFiles
	using Gadfly, Cairo, Colors
	using LinearAlgebra
	using Statistics
	
	set_default_plot_size(15cm,13cm)
	color_palette = Scale.color_discrete().f(7)[[1, 4, 6]]                                      
	themes = [Theme(default_color=color, 
			        panel_fill="white", 
			        grid_color="lightgray",            
                    background_color=parse(Colorant, "white")
			        ) for color in color_palette];
end

# ╔═╡ 3b0cfa3e-d937-11ea-2d53-255142465b91
B = readdlm("point_cloud.txt")

# ╔═╡ 3d667abc-d937-11ea-338d-e194f0379319
p1 = plot(x=B[1, :], y=B[2, :], Geom.point, themes[1],
         Guide.title("a point set, B"), themes[1],
         Coord.cartesian(fixed=true),
        )

# ╔═╡ 53c23896-d937-11ea-2e06-6b6331e0fe4c
draw(PNG("ptset_B.png", 5inch, 4inch, dpi=300), p1)

# ╔═╡ 7ca103a0-d937-11ea-0d5a-af0e3a9141a6
rotation_matrix2d(θ::Float64) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

# ╔═╡ 83746456-d937-11ea-04a1-f15bf81a17d4
begin
	θ = π * 0.8
	R_known = rotation_matrix2d(θ)
end

# ╔═╡ 8b9acb16-d937-11ea-0e02-7f6bcec57932
begin
	ϵ = 0.01 # noise
	A = R_known * B .+ ϵ * randn(size(B)...)
	
	p2 = plot(x=A[1, :], y=A[2, :], Geom.point, themes[3],
	         Guide.title("point set, A"), themes[3],
	         Coord.cartesian(fixed=true),
	    )
end

# ╔═╡ a85a6086-d937-11ea-3fed-21f714cad837
draw(PNG("ptset_A.png", 5inch, 4inch, dpi=300), p2)

# ╔═╡ b04a9cfa-d937-11ea-38cd-e1a2a9d588e3
p3 = plot(layer(x=A[1, :], y=A[2, :], Geom.point, themes[3]),                            
         layer(x=B[1, :], y=B[2, :], Geom.point, themes[1]),    
         Guide.title("orthogonal Procrustes"), themes[1],
         Coord.cartesian(fixed=true),
         Guide.manual_color_key("point set", ["A", "B"], [color_palette[3], color_palette[1]])
    )  

# ╔═╡ b7f97a90-d937-11ea-1d20-a3752908fe8c
draw(PNG("before_alignment.png", 5inch, 4inch, dpi=300), p3)

# ╔═╡ bcb90346-d937-11ea-1bbb-2f481b298a95
p4 = plot(layer(x=A[1, :], y=A[2, :], Geom.point, themes[3]),                            
         layer(x=B[1, :], y=B[2, :], Geom.point, themes[1]),
         [layer(x=[A[1, i], B[1, i]], y=[A[2, i], B[2, i]], 
            Geom.line, Theme(default_color=colorant"gray")) for i = 1:size(A)[2]]...,
         Guide.title("orthogonal Procrustes"), themes[1],
         Coord.cartesian(fixed=true),
         Guide.manual_color_key("point set", ["A", "B"], [color_palette[3], color_palette[1]])
    )  

# ╔═╡ c17420c0-d937-11ea-3fad-0b62648de42f
draw(PNG("correspondence.png", 5inch, 4inch, dpi=300), p4)

# ╔═╡ c61c4986-d937-11ea-2423-51c2154ee671
F = svd(A * B')

# ╔═╡ c9fae814-d937-11ea-3e2a-8d010deb8f54
R = F.V * F.U'

# ╔═╡ ce49e546-d937-11ea-3c90-2150e15a7ce8
isapprox(R, rotation_matrix2d(-θ), rtol=0.01)

# ╔═╡ d4242594-d937-11ea-1a8e-c1b01e9563f8
A_transformed = R * A

# ╔═╡ d86f9a84-d937-11ea-2807-6de9156292c8
p5 = plot(layer(x=A_transformed[1, :], y=A_transformed[2, :], Geom.point, themes[3]),                            
         layer(x=B[1, :], y=B[2, :], Geom.point, themes[1]),    
         Guide.title("orthogonal Procrustes"), themes[1],
         Coord.cartesian(fixed=true),
         Guide.manual_color_key("point set", ["RA", "B"], [color_palette[3], color_palette[1]])
    ) 

# ╔═╡ e461e8cc-d937-11ea-07f3-a1a8abd49594
draw(PNG("after_alignment.png", 5inch, 4inch, dpi=300), p5)

# ╔═╡ Cell order:
# ╟─ed7c391e-d937-11ea-31e0-d5884018ca29
# ╠═fd1fbc70-d936-11ea-2712-ed19c93c0c9f
# ╠═3b0cfa3e-d937-11ea-2d53-255142465b91
# ╠═3d667abc-d937-11ea-338d-e194f0379319
# ╠═53c23896-d937-11ea-2e06-6b6331e0fe4c
# ╠═7ca103a0-d937-11ea-0d5a-af0e3a9141a6
# ╠═83746456-d937-11ea-04a1-f15bf81a17d4
# ╠═8b9acb16-d937-11ea-0e02-7f6bcec57932
# ╠═a85a6086-d937-11ea-3fed-21f714cad837
# ╠═b04a9cfa-d937-11ea-38cd-e1a2a9d588e3
# ╠═b7f97a90-d937-11ea-1d20-a3752908fe8c
# ╠═bcb90346-d937-11ea-1bbb-2f481b298a95
# ╠═c17420c0-d937-11ea-3fad-0b62648de42f
# ╠═c61c4986-d937-11ea-2423-51c2154ee671
# ╠═c9fae814-d937-11ea-3e2a-8d010deb8f54
# ╠═ce49e546-d937-11ea-3c90-2150e15a7ce8
# ╠═d4242594-d937-11ea-1a8e-c1b01e9563f8
# ╠═d86f9a84-d937-11ea-2807-6de9156292c8
# ╠═e461e8cc-d937-11ea-07f3-a1a8abd49594
