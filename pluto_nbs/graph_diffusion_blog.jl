### A Pluto.jl notebook ###
# v0.19.13

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 175f88fe-f0ab-11eb-36bf-4deaf4f591fd
using Statistics, LightGraphs, Printf, GraphPlot, Random,  PlutoUI, LinearAlgebra, Colors, SimpleWeightedGraphs, StatsBase

# ╔═╡ e3115594-5a34-43f9-9f0e-fc32da23d855
md"# diffusion on graphs
_by Cory Simon_

let's learn about and simulate [using Julia] diffusion [e.g., of heat] on graphs. 🐸 

## simple, undirected, weighted graphs

suppose we have a simple, undirected, weighted graph. such a graph is composed of a set of nodes and a set of edges associated with weights. each node represents an entity, such as a human individual. each edge connects two nodes and represents a relationship between them, such as friendship-- the weight represents the strength of the relationship.
"

# ╔═╡ a98bbd67-2a56-4aea-9430-1e7d5fa63907
md"
!!! example \"example: a simple, undirected, weighted graph\"
	circles represent nodes. gray lines connecting pairs of nodes represent edges---labeled with a weight.
"

# ╔═╡ b46a8178-bd2d-41a3-a3e8-c02c96b4b42f
begin
	Random.seed!(131)
	
	nv = 7  # number of nodes

	# generate the graph
	g = SimpleWeightedGraph(nv)  # or use `SimpleWeightedDiGraph` for directed graphs
	add_edge!(g, 1, 4, round(rand(), digits=2))
	add_edge!(g, 3, 4, round(rand(), digits=2))
	add_edge!(g, 3, 5, round(rand(), digits=2))
	add_edge!(g, 4, 2, round(rand(), digits=2))
	add_edge!(g, 6, 5, round(rand(), digits=2))
	add_edge!(g, 5, 7, round(rand(), digits=2))
	add_edge!(g, 3, 7, round(rand(), digits=2))
	
	w = [ed.weight for ed in edges(g)]
	# draw the graph
	nodelabels = ["v<sub>$i</sub>" for i = 1:nv]
	Random.seed!(131)
	gplot(g, nodefillc=RGB(1.0, 1.0, 1.0), edgelabel=w,
		  nodestrokec=RGB(0.0, 0.0, 0.0), nodestrokelw=0.5, nodelabel=nodelabels)
end

# ╔═╡ 85ba66ec-ae48-443c-8a30-1673377bdd81
md"
b.t.w.
* _undirected_ graph $\implies$ the relationships represented by edges are symmetric
* _simple_ graph $\implies$ (i) at most one edge connects a pair of nodes. (ii) no loops, where an edge connects a node to itself.

🕔 introducing some notation...
* let $\mathcal{V}=\{v_1, v_2, ..., v_n\}$ be the set of nodes, with $n=|\mathcal{V}|$ the number of nodes
* let $\mathcal{E}$ be the set of edges. $\{v_i, v_j\} \in \mathcal{E}$ implies the nodes $v_i$ and $v_j$ are joined by an edge (i.e., they are adjacent).
* let $w_{ij}\in\mathbb{R}$ be the edge weight associated with two adjacent nodes $v_i$ and $v_j$. we take $w_{ij}>0$.
  * so $w_{ij}$ defined only if $\{v_i, v_j\}\in\mathcal{E}$
  * undirected graph $\implies$ symmetry $w_{ij}=w_{ji}$
* let $\mathcal{N}(v_i)=\{v_j : \{v_i, v_j\} \in\mathcal{E}\}$ be the _neighbors_ of node $v_i$---the set of nodes adjacent to node $v_i$.
"

# ╔═╡ 9ba2d787-572a-4d3d-a21f-2bd6f5757046
md"
## matrices characterizing the graph

#### adjacency matrix
the symmetric, $n\times n$ _adjacency matrix_ $A$ of a graph assembles its edge weights into a matrix. entry $(i, j)$ of $A$ is:

```math
A_{ij} :=
        \begin{cases}
        w_{ij}, &\{v_i, v_j\} \in \mathcal{E} \\ 
        0,   &   \{v_i, v_j\} \notin \mathcal{E}. 
        \end{cases}
```

so row/column $i$ in $A$ describes the connectivity of node $v_i$.

the adjacency matrix of a graph completely characterizes it; if we know $A$, we can draw the graph.

!!! example \"example: adjacency matrix\"
	the adjacency matrix of our example graph above is:
"

# ╔═╡ bbdda30e-d5f8-4243-b187-ce356a51063a
Matrix(adjacency_matrix(g))

# ╔═╡ 7e27082c-5e9e-4067-8b76-2a6dcf33410c
md"

#### degree matrix
the diagonal, $n\times n$ _degree matrix_ $D$ of a graph lists the degree of node $i$ on the diagonal entry $D_{ii}$. for a weighted graph, we define the degree of node $v_i$ as:
```math
\begin{equation}
	\deg (v_i)=\sum_{j=1}^nA_{ij} = \sum_{j \, :\, v_j\in\mathcal{N}(v_i)} w_{ij} = D_{ii}
\end{equation}
```

if $w_{ij}=1$ on all edges, the degree of a node is simply the number of neighbors it has, $\deg(v_i)=|\mathcal{N}(v_i)|$.

!!! example \"example: degree matrix\"
	the degree matrix of our example graph above is:
"

# ╔═╡ 1068fb73-5814-4b96-9bc9-a4abb2dd38a9
Matrix(degree_matrix(g))

# ╔═╡ a6df9309-eb68-47e1-b284-ccdd4735f1cb
md"
#### Laplacian matrix

the symmetric, $n\times n$ _Laplacian matrix_ $L$ of a graph is defined as $L:=D-A$. entry $(i, j)$ of $L$ is:

```math
\begin{equation}
L_{ij} :=
        \begin{cases}
        -A_{ij}, & i \neq j \\ 
         \deg(v_i), & i = j
        \end{cases}
\end{equation}
```
or, using the Kronecker delta function $\delta_{ij}$, we can write a one-liner:
```math
\begin{equation}
	L_{ij}:=\delta_{ij}\deg(v_i)-A_{ij}.
\end{equation}
```

later, we'll show $L$ is symmetric, semi-positive definite. this guarantees it will have real, non-negative eigenvalues and orthogonal eigenvectors.

it's not a coincidence that the Laplacian matrix $L$ shares the name with the Laplacian operator $\Delta$ in a continuous space; the Laplacian matrix turns out to be central to diffusion on graphs.

!!! example \"example: Laplacian matrix\"
	the Laplacian matrix of our example graph above is:
"

# ╔═╡ 1df0b874-eda6-4e85-875f-a286793c45ce
Matrix(laplacian_matrix(g))

# ╔═╡ 8df2c136-1aa5-4b19-bdb5-6562e046d284
md"
## graph temperature $\boldsymbol \phi$

now imagine each node $v_i$ is an object associated with a temperature $\phi_i$.

stacking the temperature of each node of the graph into a vector gives the _state vector_ of the graph $\boldsymbol \phi \in\mathbb{R}^n$.
as we allow adjacent nodes to exchange heat, the state vector is dynamic (changes with time), i.e. $\boldsymbol \phi = \boldsymbol \phi(t)$, 

!!! example \"example: visualizing the state vector of a graph\"
	suppose the temperature of node $i$ of our example graph is $\phi_i=i/n$. here we color each node according to its temperature (dark green = hot, white = cold).
"

# ╔═╡ 614b581f-b29c-4609-b45d-500f22bcbcbd
begin
	cmap′ = colormap("Greens", nv)
	Random.seed!(131)
    gplot(g, nodefillc=cmap′, edgelabel=w,
	  nodestrokec=RGB(0.0, 0.0, 0.0), nodestrokelw=0.5, nodelabel=nodelabels)
end

# ╔═╡ 8db987ad-a429-4d55-9483-e01cdd834647
md"## Newton's law of cooling

each node exchanges heat with its neighbors.

_Newton's law of cooling_ specifies that heat continuously (in time) transfers from node $v_i$ to an adjacent node $v_j$ at a rate:

```math
\begin{equation}
\alpha w_{ij} (\phi_i - \phi_j).
\end{equation}
```
so, the rate of heat exchange between node $v_i$ and $v_j$ increases as
* the thermal diffusivity, $\alpha \in \mathbb{R}$, increases
* the strength of the connection between them, $w_{ij}$, increases
* the difference in temperature between them, $|\phi_i-\phi_j|$, increases. 

nodes $v_i$ and $v_j$ exchange heat only if they are adjacent, i.e. only if $\{v_i, v_j\} \in \mathcal{E}$. thus the interesting thing here is that the connectivity of the nodes in the graph will play an important role in the dynamics of $\boldsymbol \phi$.

✔ if node $i$ is hotter than node $j$, then $\phi_i > \phi_j$ and indeed heat transfers from node $i$ to node $j$. if instead node $j$ is hotter, the quantity above is negative and, actually, heat transfers from node $j$ to node $i$.	

## dynamics of $\phi_i$

the differential equation (DE) describing the dynamics of the temperature of node $i$, $\phi_i$, accounts for the heat transfer to/from all of its neighboring nodes:
```math
\begin{equation}
	\frac{d\phi_i}{dt}=-\sum_{j\, :\, v_j\in\mathcal{N}(v_i)}\alpha w_{ij} (\phi_i - \phi_j) \text{ for } i \in \{1, ..., n\}.
\end{equation}
```
i.e., the rate at which $\phi_i$ decreases is given by the _net_ rate of heat transfer from node $i$ to its neighboring nodes.

this is a _system_ of _coupled_ differential equations. e.g., the two differential equations for $\phi_i$ and $\phi_j$ are clearly coupled if $v_j \in \mathcal{N}(v_i)$ because (i) the rate of change of $\phi_i$ depends on the temperaure of node $v_j$ and (ii) the rate of change of $\phi_j$ depends on the temperature of node $v_i$. generally, if there exists a path between node $v_i$ and $v_j$ (i.e., if nodes $v_j$ and $v_j$ are connected) in the graph, then $\phi_i$ and $\phi_j$ are coupled since heat can transfer along this path.

🚀 voila, we've posed diffusion of heat on a graph, resulting in an autonomous dynamical system. the DEs above give the evolution of the graph state (temperature) vector $\boldsymbol \phi= \boldsymbol \phi(t)$.

## understanding the dynamics of $\phi_i$

to gain further understanding, let's write the dynamic model for $\phi_i(t)$ a bit differently.
```math
\begin{align}
	\frac{d\phi_i}{d(\alpha t)}&=
-\phi_i\sum_{j\,:\,v_j\in\mathcal{N}(v_i)}w_{ij} + \sum_{j\,:\,v_j\in\mathcal{N}(v_i)}w_{ij}\phi_j \\
&= \left(\sum_{j\,:\,v_j\in\mathcal{N}(v_i)}w_{ij}\right) \left(\frac{\sum_{j\,:\,v_j\in\mathcal{N}(v_i)}w_{ij}\phi_j}{\sum_{j\,:\,v_j\in\mathcal{N}(v_i)}w_{ij}} - \phi_i \right)\\
&=\deg(v_i)\left( \langle \phi_j \rangle_{\mathcal{N}(v_i)} - \phi_i \right),
\end{align}
```
with $\langle \cdot \rangle_{\mathcal{N}(v_i)}$ the weighted average over the neighboring nodes of $v_i$. this gives us some intuition: the temperature of node $v_i$, $\phi_i$, increases at a rate proportional to the difference between the weighted average of the temperature of its neighbors and its temperature. so, node $i$ heats up if the (weighted) average of the temperature of its neighboring nodes is greater than its temperature. on the other hand, node $v_i$ cools down if the (weighted) average of the temperature of its neighboring nodes is lower than its temperature. 😎

## the dynamics of the state vector $\boldsymbol \phi$

let's write a DE governing the dynamics of the state vector $\boldsymbol \phi=\boldsymbol \phi(t)$ of the graph by converting the system of DEs above into vector form $\frac{d\boldsymbol \phi}{dt}=\cdots$. this facilitates finding the solution.

to achieve this, we write the DE for $\phi_i$ differently.
```math
\begin{align}
	\frac{d\phi_i}{d(\alpha t)}&=-\sum_{j=1}^n A_{ij} (\phi_i - \phi_j)\\
&= -\phi_i\sum_{j=1}^n A_{ij} + \sum_{j=1}^n A_{ij}\phi_j \\
&= -\deg(v_i)\phi_i + \sum_{j=1}^n A_{ij}\phi_j \\
&= -\sum_{j=1}^n \left(\delta_{ij} \deg(v_i)-A_{ij} \right)\phi_j\\

\end{align}
```

🙌 element $(i,j)$ of the Laplacian matrix $L$ appears!

by the definition of matrix multiplication...

!!! note \"the diffusion equation on a graph \" 
	 the following vector differential equation (DE) governs the dynamics of the state vector $\boldsymbol \phi=\boldsymbol \phi(t)$ of the graph giving the temperature of the nodes:
	```math
	\begin{equation}
		\frac{d\boldsymbol \phi}{dt} = - \alpha L \boldsymbol \phi
	\end{equation}
	```
	with $\alpha$ the thermal diffusivity and $L$ the Laplacian matrix of the graph. this DE is subject to an initial condition $\boldsymbol \phi(t=0)=\boldsymbol \phi_0$.

clearly, the dynamics of $\boldsymbol \phi(t)$ depend on the connectivity of the graph, through the Laplacian matrix $L$.
"

# ╔═╡ 748f3aff-e0eb-426e-acd3-826d8035da1e
md"## analogy with Laplacian operator

there is an analogy here, $L::\Delta$ with $\Delta$ the Laplace operator. 
$L$ is an operator on vectors for diffusion on graphs; $\Delta$ is an operator on functions for diffusion on continuous manifolds. the Laplace operator $\Delta$ operates on a function giving the temperature at each point on the manifold. the Laplacian matrix $L$ operates on (i.e., left-multiplies) a vector giving the temperature at each node of the graph. 

let's see how the Laplacian matrix $L$ operates on (linearly transforms) a vector $\boldsymbol \phi \in \mathbb{R}^n$. we're considering the function $f(\boldsymbol \phi) = L \boldsymbol \phi$ that linearly transforms $\boldsymbol \phi$ to a new vector $L\boldsymbol \phi \in \mathbb{R}^n$. element $i$ of $L\boldsymbol \phi$ appears as the right-hand-side of the DE for $\phi_i$ above:

```math
\begin{equation}
	(L \boldsymbol \phi)_i=\deg(v_i)\left( \phi_i- \langle \phi_j \rangle_{\mathcal{N}(v_i)} \right)
\end{equation}
```

so the Laplacian matrix operates on the state vector $\boldsymbol \phi$ of the graph to give a new state vector $L\boldsymbol \phi$ whose element $i$ is given by the degree of node $v_i$ times the difference between the temperature of node $v_i$ and the (weighted) average temperature of its neighbors.

the Laplacian operator on a function defined over a line gives conceptually the same thing. let $x\in\mathbb{R}$ be a point on the line. then given a function $\theta(x)$, the Laplacian operator looks like the following if we do a centered finite difference approximation with $\delta x$ the spatial step:

```math
\begin{equation}
	\Delta \theta = (\Delta \theta)(x) = \frac{d^2\theta}{dx^2}\approx \frac{2}{(\delta x)^2} \left( \frac{\theta(x+\delta x) + \theta(x-\delta x)}{2} - \theta(x)  \right)
\end{equation}
```
again we see a difference between (i) the value of $\theta$ at $x$ and (ii) the average of the temperature at the two points neighboring $x$. and the $2$ in front is the degree of node \"$x$\". I guess we could say that when we discretize the diffusion equation on a line, we convert it to a diffusion problem on a graph---a regular graph where the degree of every node is two (ignoring boundary conditions here). (❓ not sure why the sign is flipped...)
"

# ╔═╡ 87068939-2ec9-4b59-aff0-760a6f8a9fb1
begin
	g3 = SimpleGraph(5)
	for i = 1:4
		add_edge!(g3, i, i+1)
	end
	gplot(g3, nodefillc=RGB(1.0, 1.0, 1.0),
	  nodestrokec=RGB(0.0, 0.0, 0.0), nodestrokelw=0.5, nodelabel=["x - 2δx", "x - δx","x", "x+δx", "x+2δx"])
end

# ╔═╡ 4391b845-b4d0-418c-a293-52b493eb2284
md"
n.b. one more way to see $L$ as an operator is that element $i$ of $L\boldsymbol \phi$ is a sum of the differences in temperature between node $v_i$ and each of its neighbors:
```math
\begin{equation}
	(L \boldsymbol \phi)_i= \sum_{j\, :\, v_j \in \mathcal{N}(v_i)} w_{ij} \left( \phi_i- \phi_j \right)
\end{equation}
```
"

# ╔═╡ 287ecc14-d677-4c73-a2aa-3ed422197eb9
md"## the Laplacian matrix has an eigenvalue of $0$

suppose the temperature of all nodes is the same, $\phi_i = 1/\sqrt{n}$ for all $i \in \{1,...,n\}$. then,
```math
\begin{equation}
	(L \boldsymbol \phi)_i=\deg(v_i)\left( \phi_i- \langle \phi_j \rangle_{\mathcal{N}(v_i)} \right) = 0 \quad \forall i\in \{1, ..., n\}.
\end{equation}
```
that is, $L \boldsymbol \phi = \mathbf{0}$ and the constant vector $\boldsymbol \phi = \mathbf{1} / \sqrt{n}$ is an eigenvector of $L$ with eigenvalue $0$. in other words, $L$ has a non-trivial null space.

if the graph is not connected, the null space of $L$ has additional dimensions. suppose the graph has $k$ connected components. then we can assign each node in component $k$ to have a temperature $\phi_i=c_k$ with $c_k$ a constant. the constants $\{c_1, c_2, ..., c_k\}$ can be distinct, with $L \boldsymbol \phi = \mathbf{0}$ still. one can show that the Laplacian matrix of a graph with $k$ connected components has a zero eigenvalue with multiplicity $k$. very cool IMO.
"

# ╔═╡ a5b602ed-9ff6-4986-bab8-76e78c886d57

md"
## conservation of energy

for kicks, we show that energy $\sum_i \phi_i(t)$ is a conserved quantity in the dynamical system. i.e., the energy does not change with time:
```math
\begin{equation}
	\frac{d}{dt} \sum_{i=1}^n  \phi_i = \sum_{i=1}^n \frac{d\phi_i}{dt} = 0.
\end{equation}
```
thus, the energy must stay at its initial value for all time, $\sum_{i=1}^n  \phi_i(t)=\sum_{i=1}^n  \phi_{i, 0}\, \forall t\geq 0$.

well,
```math
\begin{align}
    \alpha^{-1}\sum_{i=1}^n \frac{d\phi_i}{dt} &= -\sum_{i=1}^n \sum_{j\,:\,v_j\in\mathcal{N}(v_i)} w_{ij} (\phi_i - \phi_j) \\
&= -\sum_{i=1}^n \sum_{j=1}^n  A_{ij} (\phi_i - \phi_j) \\
&= \sum_{i=1}^n \sum_{j=1}^n  A_{ij} \phi_j - \sum_{i=1}^n \sum_{j=1}^n  A_{ij} \phi_i.
\end{align}
```

let's mess with the first term on the right. switch order of sums. use symmetry. swap arbitrary index variables $i$ and $j$.
```math
\begin{equation}
\sum_{i=1}^n \sum_{j=1}^n  A_{ij}\phi_j =\sum_{j=1}^n \sum_{i=1}^n A_{ij}\phi_j= \sum_{j=1}^n \sum_{i=1}^n A_{ji} \phi_j = \sum_{i=1}^n \sum_{j=1}^n A_{ij} \phi_i
\end{equation}
```

aha! the two terms cancel. so we've showed the energy in the graph is conserved. this happens because all of the heat stays in the system. there is no source or sink.

pretty confident we can also show, more strictly, energy is conserved within connected components of the graph...
"

# ╔═╡ 490fb33a-9d9f-401e-8906-801ab82dc18e
md"
## the Laplacian matrix $L$ is symmetric, semi-positive definite

the Laplacian matrix $L$ is symmetric since $L_{ij}=\delta_{ij} \deg(v_i) -A_{ij} = L_{ji}= \delta_{ji} \deg(v_j) -A_{ji}$. we'll also show it's semi-positive definite, so that $\boldsymbol \phi ^\intercal L \boldsymbol \phi \geq 0$ for all possible $\boldsymbol \phi \in \mathbb{R}^n$.

we already figured out $(L\boldsymbol \phi)_i$ so we can just use the definition of the dot product to write:
```math
\begin{align}
	\boldsymbol \phi ^\intercal L \boldsymbol \phi & = \sum_{i=1}^n \phi_i(L\boldsymbol \phi)_i \\
	&= \sum_{i=1}^n \phi_i \sum_{j=1}^n A_{ij} (\phi_i - \phi_j) \\
	&= \sum_{i=1}^n \sum_{j=1}^n A_{ij} (\phi_i^2 - \phi_i \phi_j) \\
	&= \frac{1}{2} \sum_{i=1}^n  \sum_{j=1}^n  A_{ij} (\phi_i^2 - 2\phi_i \phi_j + \phi_i^2) \\
	&= \frac{1}{2} \sum_{i=1}^n  \sum_{j=1}^n A_{ij} (\phi_i^2 - 2\phi_i \phi_j + \phi_j^2) \\
	&= \frac{1}{2} \sum_{i=1}^n  \sum_{j=1}^n A_{ij} (\phi_i- \phi_j)^2 \\
	&= \sum_{\{i, j\} \,:\, \{v_i, v_j\} \in \mathcal{E}} w_{ij} (\phi_i- \phi_j)^2 \geq 0.
\end{align}
```
so the number $\boldsymbol \phi ^\intercal L \boldsymbol \phi$ is the weighted sum of square differences between the temperatures of adjacent nodes on the graph. and, it's greater than or equal to zero. it's equal to zero only when the temperatures of all nodes in the graph have the same temperature as the other nodes in the connected component it belongs to. note, this quantity plays an important role for embedding the nodes of a graph into a continuous space.

to review some basic linear algebra, we next convince ourselves of some properties of $L$ that arise because $L$ is real and symmetric semi-positive-definite.

#### $\implies L$ has real eigenvalues and eigenvectors

if we let $^*$ mean \"take the complex conjugate\", then $L^*=L$ since $L$ is real. and $L=L^\intercal$ since it's symmetric.
suppose $L\mathbf{v}=\lambda \mathbf{v}$. we wish to show $\lambda$ is real. we can do that if we show $\lambda^*=\lambda$. keep in mind $\mathbf{v}$ could, for all we know at this point, have complex entries. 

first, note 
```math
\begin{equation}
	\mathbf{v} ^{\intercal *} (L \mathbf{v}) =\mathbf{v}^{\intercal *} (\lambda \mathbf{v}) = \lambda ||\mathbf{v}||^2.
\end{equation}
```

second, note:
```math
\begin{equation}
(\lambda \mathbf{v})^*=(L\mathbf{v})^* \implies \lambda^*\mathbf{v}^*=L^*\mathbf{v}^*=L\mathbf{v}^* \implies \mathbf{v}^{\intercal *}L^\intercal =  \mathbf{v}^{\intercal *}L=\lambda^*\mathbf{v}^{\intercal *}.
\end{equation}
```
allowing us to write $\mathbf{v} ^{\intercal *} L \mathbf{v}=\lambda^*||\mathbf{v}||^2$. 

third, equating the two ways we're able to write $\mathbf{v} ^{\intercal *} L \mathbf{v}$, we find $\lambda = \lambda^*$.

we can show the eigenvectors are real with a similar argument.

#### $\implies $ the eigenvectors of $L$ are orthogonal

suppose $L\mathbf{v}=\lambda \mathbf{v}$ and $L\mathbf{u}=\mu\mathbf{u}$. we wish to show $\mathbf{v}^\intercal \mathbf{u}=0$ provided $\mu \neq \lambda$.

first:
```math
\begin{equation}
	\mathbf{v}^\intercal [L\mathbf{u}=\mu\mathbf{u}] \implies \mathbf{v}^\intercal L\mathbf{u} = \mu \mathbf{v}^\intercal \mathbf{u}.
\end{equation}
```
second,
```math
\begin{equation}
	\mathbf{u}^\intercal [L\mathbf{v}=\lambda\mathbf{v}] \implies \mathbf{u}^\intercal L\mathbf{v} = \lambda \mathbf{u}^\intercal \mathbf{v}.
\end{equation}
```
and transposing:
```math
\begin{equation}
	\mathbf{v}^\intercal L^\intercal\mathbf{u} =\mathbf{v}^\intercal L\mathbf{u}= \lambda \mathbf{v}^\intercal \mathbf{u}.
\end{equation}
```
third, put the two together:
```math
\begin{equation}
	\mathbf{v}^\intercal L\mathbf{u}= \lambda \mathbf{v}^\intercal \mathbf{u}=\mu \mathbf{v}^\intercal \mathbf{u}.
\end{equation}
```
which implies we must have $\mathbf{v}^\intercal \mathbf{u}=0$ since $\lambda\neq\mu$.

what happens if an eigenvalue is repeated? ❓ my argument here should only partially convince you.

!!! note
	an important and useful (we'll see below) implication is that we can construct an orthonormal basis of $\mathbb{R}^n$ using the eigenvectors of $L$.

#### $\implies L$ has non-negative eigenvalues

suppose $L\mathbf{v}=\lambda \mathbf{v}$. then:
```math
\begin{equation}
0\leq \mathbf{v}^\intercal L\mathbf{v} = \lambda \mathbf{v}^\intercal \mathbf{v} = \lambda ||\mathbf{v}||^2
\end{equation}
```
since $L$ is semi-positive definite. well $||\mathbf{v}||^2>0$ for an eigenvector so we have $\lambda \geq 0$.

!!! example \"Example: eigenvectors and eigenvalues of the Laplacian matrix of our example graph\"
	we print the computed eigenvalues of the Laplacian matrix $L$ of our example graph below. we have one eigenvalue that is zero (not repeated b/c graph is connected), and the rest are positive. 
"

# ╔═╡ 26b7a075-e332-49e6-b1c1-ed45b3af26e3
λs, Q = eigen(Matrix(laplacian_matrix(g)));

# ╔═╡ f0128096-c50c-4e12-a968-df4d4df7684f
with_terminal() do
	for i in vertices(g)
		println("eigenvalue λ_$i = ", round(λs[i], digits=3))
	end
end

# ╔═╡ 4f8be09c-9b38-479d-a617-8a220f290602
md"
## analytical solution to the diffusion equation on graphs

the key to finding the solution to the diffusion equation on graphs is to express the solution $\boldsymbol \phi = \boldsymbol \phi(t)$ at any given time $t$ as a linear combination of the eigenvectors of the Laplacian matrix $L$:
```math
\begin{equation}
	\boldsymbol \phi(t)= \sum_{i=1}^n \beta_i(t) \mathbf{v}_i
\end{equation}
```
with eigenvector $\mathbf{v}_i$ satisfying $L\mathbf{v}_i=\lambda_i\mathbf{v}_i$ and $||\mathbf{v}_i||=1$ and $\beta_i=\beta_i(t)$ time-varying coefficients to be determined. 
since the Laplacian matrix $L$ is symmetric and real, we can construct its set of eigenvectors $\{\mathbf{v}_1, \mathbf{v}_2, ..., \mathbf{v}_n\}$ to form an orthonormal basis for $\mathbb{R}^n$, justifying this expansion. 

substituting this expansion of $\boldsymbol \phi(t)$ into the diffusion equation, we have:
```math
\begin{align}
	\frac{d}{dt} \sum_{i=1}^n \beta_i(t) \mathbf{v}_i &= - \alpha L \left(\sum_{i=1}^n \beta_i(t) \mathbf{v}_i\right) \\
\sum_{i=1}^n \frac{d\beta_i}{dt} \mathbf{v}_i & = -\alpha \sum_{i=1}^n \beta_i(t) L \mathbf{v}_i \\
&= -\alpha \sum_{i=1}^n \beta_i(t) \lambda_i \mathbf{v}_i
\end{align}
```

bringing everything on one side,
```math
\begin{equation}
\sum_{i=1}^n \left(\frac{d\beta_i}{dt} + \alpha \lambda_i \beta_i(t) \right) \mathbf{v}_i 
=0
\end{equation}
```
since the $\mathbf{v}_i$'s are linearly independent, we the following independent ODE must hold for each $\beta_i(t)$
```math
\begin{equation}
	\frac{d\beta_i}{dt} + \alpha \lambda_i \beta_i=0 \text{ for } i \in \{1,...,n\}
\end{equation}
```
whose solution is:
```math
\begin{equation}
	\beta_i(t) = \beta_{i, 0} e^{-\alpha \lambda_i t} \text{ for } i \in \{1,...,n\}
\end{equation}
```
the coefficients $\beta_{i, 0}$ follow from invoking the initial condition $\boldsymbol \phi(t=0)= \boldsymbol  \phi_0$:
```math
\begin{equation}
\boldsymbol \phi(t=0) = \boldsymbol \phi_0 = \sum_{i=1}^n \beta_{i, 0}\mathbf{v}_i. 
\end{equation}
```
multiply the above equation on the left by $\mathbf{v}_i^\intercal$ and consider that the eigenvectors are orthonormal to see:
```math
\begin{equation}
	\beta_{i,0}= \mathbf{v}_i^\intercal \boldsymbol \phi_0.
\end{equation}
```

!!! note \"the solution to the diffusion equation on graphs\"
	the solution to the diffusion equation on a graph is:
	```math
	\begin{equation}
		\boldsymbol \phi(t) = \sum_{i=1}^n \mathbf{v}_i^\intercal \boldsymbol \phi_0 e^{-\alpha \lambda_i t} \mathbf{v}_i,
	\end{equation}
	```
	where $\mathbf{v}_i$ and $\lambda_i$ are eigenvector and eigenvalue $i$ of the Laplacian matrix $L$ of the graph.

since $L$ is symmetric positive semi-definite, its eigenvalues $\lambda_i\geq 0$. thus, the solution is bounded: each exponential response mode associated with a non-zero eigenvalue decays with time, leaving a constant vector. the steady-state (as $t\rightarrow \infty$) solution involves the eigenvectors of the Laplacian matrix associated with the zero eigenvalue, which lie in its null space. if the graph is connected, then the steady state solution is where the temperature of each node is the same, $\phi_i=\sum_{j=1}^n \phi_{j, 0} / n$.
"

# ╔═╡ d112d370-f127-48c8-9349-f44516e31102
md"## visualization of $\boldsymbol \phi(t)$ on our graph

!!! example \"Example: visualization of diffusion on our example graph\"
	suppose the initial state vector $\boldsymbol \phi(t=0)=\boldsymbol \phi_0 = [0, 0, 0, 1, 0, 0, 0]$. the visualizations below show the analytical solution to the diffusion equation on the example graph for $\alpha=1$.
"

# ╔═╡ bec1d199-a318-4deb-be1b-50ac41de4118
ts = [0.0, 0.15, 0.5, 2.0, 100.0];

# ╔═╡ 6a567ee2-0a3b-4f94-bc66-165256441040
with_terminal() do
	print("t = ", ts[1])
end

# ╔═╡ 096cfb52-daf9-4142-9a32-6b85a4320c0b
begin
	ϕ0 = zeros(length(vertices(g)))
	ϕ0[4] = 1.0
	
	cmap = colormap("Greens", 101)
	function ϕ_to_colors(ϕvec)
		ids = ceil.(Int, ϕvec * 100) .+ 1
		ids[ids .< 1] .= 1
		ids[ids .> 101] .= 101
		return cmap[ids]
	end
	
	function ϕ(t)
		c₀ = Q' * ϕ0
		c = zeros(nv)
		for i = 1:nv
			c .+= c₀[i] * Q[:, i] * exp(-λs[i] * t)
		end
		return c
	end
	
	function viz_solution(t)
		Random.seed!(131)
		nodefillc = ϕ_to_colors(ϕ(t))
		gplot(g, nodefillc=nodefillc, edgelabel=w,
			  nodestrokec=RGB(0.0, 0.0, 0.0), nodestrokelw=0.5, 
			  nodelabel=nodelabels)
	end
	
	viz_solution(0.0)
end

# ╔═╡ eeb3a282-f53e-4484-8ddd-82daa43a178d
with_terminal() do
	print("t = ", ts[2])
end

# ╔═╡ 75cfcbbb-7896-4dc8-833f-5d3ce14ed02e
viz_solution(ts[2])

# ╔═╡ 11464714-e089-4e53-b6c0-e61c2526e7ef
with_terminal() do
	print("t = ", ts[3])
end

# ╔═╡ 02a35bd3-3bf2-42bc-8c91-83794391ddc8
viz_solution(ts[3])

# ╔═╡ 23498641-43ad-4948-aede-b7151e6c0087
with_terminal() do
	print("t = ", ts[4])
end

# ╔═╡ aec75d65-2278-402d-ab29-17089ab9f5ed
viz_solution(ts[4])

# ╔═╡ 768908ad-6756-476d-9608-5b2b53e38591
with_terminal() do
	print("t = ", ts[5])
end

# ╔═╡ ca62288b-9b41-47a6-9790-4bce4c016068
viz_solution(ts[5])

# ╔═╡ 61690352-55bc-4b86-ab5b-3eaf3efce93f
md"## diffusion on graph as a gradient flow

we can describe diffusion on a graph as a gradient flow

```math
\begin{equation}
\frac{d\phi}{dt} = - \nabla_\phi f(\phi).
\end{equation}
```
where $f(\phi)$ is a potential function that monotonically decreases as time goes on. this is seen by the differential equation: $\phi$ evolves in the direction of the negative gradient of the function $f(\phi)$. hence, $f(\phi)$ will decrease over time.

above, we showed the diffusion equation on a graph is:
```math
\begin{equation}
\frac{d\phi}{dt} = - \alpha L \phi.
\end{equation}
```

with our matrix calculus skills, and considering $L$ is symmetric, we notice that
```math
\begin{equation}
\nabla_\phi \alpha \phi^\intercal L \phi = 2 \alpha L \phi.
\end{equation}
```
therefore, the potential function here is:
```math
\begin{equation}
f(\phi):= \frac{1}{2}\alpha \phi^\intercal L \phi,
\end{equation}
```
called the _Dirichlet energy_ of the graph. in words, (see above) the Dirichlet energy is the weighted sum of square differences between the temperature of adjacent nodes in the graph. intuitively, the Dirichlet energy decreases over time because heat travels from hot to cold nodes, smoothing out temperature differences between adjacent nodes. we may think of the Dirichlet energy as a metric of smoothness of the temperature function on top of the graph (again, thinking of $\phi$ as a function that assigns a temperature to each node of the graph).
"

# ╔═╡ ffeae924-5ff3-4da9-8d0a-d5933dcf6c92
md"## Julia code for diffusion on graphs
the code for the above figures is hidden in cells in this Pluto.jl notebook, which you can download via the button on the top right. 

here, I'll start a new example showing Julia code for simulating diffusion on graphs. this relies on `LightGraphs.jl`, `SimpleWeightedGraphs.jl`, `StatsBase.jl`, `Colors.jl`, and `GraphPlot.jl`.

#### generate a graph
"

# ╔═╡ cf89067a-2b1b-4995-a5e8-c6d06cbd4004
begin
	nb_nodes = 8 # number of nodes
	nb_edges = 12  # number of edges
	
	# construct a simple, weighted, undirected graph
	graph = SimpleWeightedGraph(nb_nodes)
	while length(edges(graph)) < nb_edges
		# choose random pair of vertices to connect
		v_i, v_j = sample(1:nb_nodes, 2, replace=false)
		# generate random edge weight
		w_ij = round(rand(), digits=2)
		add_edge!(graph, v_i, v_j, w_ij)
	end
	
	# plot graph
	Random.seed!(97330) # so the layout is the same
	gplot(graph, nodefillc=RGB(1.0, 1.0, 1.0), 
		  edgelabel=[ed.weight for ed in edges(graph)],
		  nodelabel=["v<sub>$i</sub>" for i = 1:nb_nodes],
		  nodestrokec=RGB(0.0, 0.0, 0.0), nodestrokelw=0.5)
end

# ╔═╡ fb7859a8-c671-4517-a5fa-4eccdaeb22eb
md"
#### compute eigenvectors and eigenvalues of its Laplacian matrix
"

# ╔═╡ 7c79031a-d566-4844-a6c6-177e7709b571
L = laplacian_matrix(graph) # Laplacian matrix of the graph

# ╔═╡ c5344c65-43d7-4b22-84c6-60ff7cb09e1d
# compute eigenvalues and eigenvectors
eigs = eigen(Matrix(L));

# ╔═╡ ef4745a6-c8bf-4a8a-a8fb-9083134be316
eigs.values # eigenvalues

# ╔═╡ 2dc2202d-9e23-4be4-8a67-c8abe850d51e
eigs.vectors # eigenvectors (in the columns)

# ╔═╡ f8b8f8bd-cfa4-4915-acb1-bc918bb511e1
# orthonormal eigenvectors, told ya
round.(eigs.vectors' * eigs.vectors, digits=10)

# ╔═╡ 0677112e-8c77-4373-a0db-707bb1915c5e
md"
#### construct the solution to the diffusion equation on the graph
"

# ╔═╡ ba43dad8-9db6-4ef2-9f0e-8899006a658f
begin
	# thermal diffusivity
	α = 1.0
	# initial condition
	ϕ₀ = zeros(nb_nodes)
	ϕ₀[rand(1:nb_nodes)] = 1.0 # give all heat to a random node
	
	# returns state vector of graph at time t
	function ϕ_vector(t::Float64)
		# will be state vector of graph
		ϕ = zeros(nb_nodes)
		for i = 1:nb_nodes
			# eigenvector i
			v_i = eigs.vectors[:, i]
			# eigenvalue i
			λ_i = eigs.values[i]
			# add response mode to solution
			ϕ .+= dot(v_i, ϕ₀) * exp(-λ_i * α * t) * v_i
		end
		return ϕ
	end
end

# ╔═╡ 280ec4d5-e352-49e5-8bab-9ea93ec1c628
md"
#### visualize the solution over time

click the clock to run the time!
"

# ╔═╡ cf28c22f-96ad-4e50-a2ff-d6ba1df6e448
# first, color scheme for the nodes
orange_cmap = colormap("Oranges", 100)

# ╔═╡ 02e7d340-6b44-4066-853a-27bc03604d26
# map a number to a color
ϕᵢ_to_color(ϕᵢ::Float64) = orange_cmap[ceil(Int, ϕᵢ * 99) + 1]

# ╔═╡ 74a904f5-0194-4688-b074-74f1d738c62d
@bind j Clock(0.1, true) # set up clock for time

# ╔═╡ e762ac49-5667-4759-b8b2-a2a25c3a80f9
begin
	# compute current time
	Δt = 0.025 # time increment
	t = j * Δt # time
	
	# compute solution
	ϕ_vector_t = ϕ_vector(t)
end

# ╔═╡ 9c8586a5-67da-4d7e-8192-d89cc9370f40
with_terminal() do
	@printf("t = %.2f", t)
	@printf("\tenergy = %.2f", sum(ϕ_vector_t))
end

# ╔═╡ 1f9e1123-de47-4e07-8e47-0595a10358fd
begin
	# visualize solution
	Random.seed!(97330) # so the layout is the same
	gplot(graph, nodefillc=ϕᵢ_to_color.(ϕ_vector_t), 
		  edgelabel=[ed.weight for ed in edges(graph)],
		  nodelabel=["v<sub>$i</sub>" for i = 1:nb_nodes],
		  nodestrokec=RGB(0.0, 0.0, 0.0), nodestrokelw=0.5)
end

# ╔═╡ 0366e3c4-5ec5-4c96-afdf-17407136a691
md"# references

I found the following references pivotal for writing this blog post:
* Radu Horaud. A Short Tutorial on Graph Laplacians, Laplacian Embedding, and Spectral Clustering. [link](https://csustan.csustan.edu/~tom/Clustering/GraphLaplacian-tutorial.pdf)
* Laplacian matrix page on Wikipedia. [link](https://en.wikipedia.org/wiki/Laplacian_matrix)"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
GraphPlot = "a2cc645c-3eea-5389-862e-a155d0052231"
LightGraphs = "093fc24a-ae57-5d10-9952-331d41423f4d"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
SimpleWeightedGraphs = "47aef6b3-ad0c-573a-a1e2-d07658019622"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
Colors = "~0.12.8"
GraphPlot = "~0.4.4"
LightGraphs = "~1.3.5"
PlutoUI = "~0.7.9"
SimpleWeightedGraphs = "~1.1.1"
StatsBase = "~0.33.9"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "f87e559f87a45bece9c9ed97458d3afe98b1ebb9"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.1.0"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "32a2b8af383f11cbb65803883837a149d10dfe8a"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.10.12"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dc7dedc2c2aa9faf59a55c622760a25cbefbe941"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.31.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[Compose]]
deps = ["Base64", "Colors", "DataStructures", "Dates", "IterTools", "JSON", "LinearAlgebra", "Measures", "Printf", "Random", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "c6461fc7c35a4bb8d00905df7adafcff1fe3a6bc"
uuid = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
version = "0.9.2"

[[DataAPI]]
git-tree-sha1 = "ee400abb2298bd13bfc3df1c412ed228061a2385"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.7.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4437b64df1e0adccc3e5d1adbc3ac741095e4677"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.9"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[GraphPlot]]
deps = ["ArnoldiMethod", "ColorTypes", "Colors", "Compose", "DelimitedFiles", "LightGraphs", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "dd8f15128a91b0079dfe3f4a4a1e190e54ac7164"
uuid = "a2cc645c-3eea-5389-862e-a155d0052231"
version = "0.4.4"

[[Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "81690084b6198a2e1da36fcfda16eeca9f9f24e4"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.1"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LightGraphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "432428df5f360964040ed60418dd5601ecd240b6"
uuid = "093fc24a-ae57-5d10-9952-331d41423f4d"
version = "1.3.5"

[[LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "6a8a2a625ab0dea913aba95c11370589e0239ff0"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.6"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "2ca267b08821e86c5ef4376cffed98a46c2cb205"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.1"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "c8abc88faa3f7a3950832ac5d6e690881590d6dc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "1.1.0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "5f6c21241f0f655da3952fd60aa18477cf96c220"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.1.0"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[SimpleWeightedGraphs]]
deps = ["LightGraphs", "LinearAlgebra", "Markdown", "SparseArrays", "Test"]
git-tree-sha1 = "f3f7396c2d5e9d4752357894889a87340262f904"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.1.1"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "5b2f81eeb66bcfe379947c500aae773c85c31033"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.8"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "fed1ec1e65749c4d96fc20dd13bea72b55457e62"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.9"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─175f88fe-f0ab-11eb-36bf-4deaf4f591fd
# ╟─e3115594-5a34-43f9-9f0e-fc32da23d855
# ╟─a98bbd67-2a56-4aea-9430-1e7d5fa63907
# ╟─b46a8178-bd2d-41a3-a3e8-c02c96b4b42f
# ╟─85ba66ec-ae48-443c-8a30-1673377bdd81
# ╟─9ba2d787-572a-4d3d-a21f-2bd6f5757046
# ╟─bbdda30e-d5f8-4243-b187-ce356a51063a
# ╟─7e27082c-5e9e-4067-8b76-2a6dcf33410c
# ╟─1068fb73-5814-4b96-9bc9-a4abb2dd38a9
# ╟─a6df9309-eb68-47e1-b284-ccdd4735f1cb
# ╟─1df0b874-eda6-4e85-875f-a286793c45ce
# ╟─8df2c136-1aa5-4b19-bdb5-6562e046d284
# ╟─614b581f-b29c-4609-b45d-500f22bcbcbd
# ╟─8db987ad-a429-4d55-9483-e01cdd834647
# ╟─748f3aff-e0eb-426e-acd3-826d8035da1e
# ╟─87068939-2ec9-4b59-aff0-760a6f8a9fb1
# ╟─4391b845-b4d0-418c-a293-52b493eb2284
# ╟─287ecc14-d677-4c73-a2aa-3ed422197eb9
# ╟─a5b602ed-9ff6-4986-bab8-76e78c886d57
# ╟─490fb33a-9d9f-401e-8906-801ab82dc18e
# ╟─26b7a075-e332-49e6-b1c1-ed45b3af26e3
# ╟─f0128096-c50c-4e12-a968-df4d4df7684f
# ╟─4f8be09c-9b38-479d-a617-8a220f290602
# ╟─d112d370-f127-48c8-9349-f44516e31102
# ╟─bec1d199-a318-4deb-be1b-50ac41de4118
# ╟─6a567ee2-0a3b-4f94-bc66-165256441040
# ╟─096cfb52-daf9-4142-9a32-6b85a4320c0b
# ╟─eeb3a282-f53e-4484-8ddd-82daa43a178d
# ╟─75cfcbbb-7896-4dc8-833f-5d3ce14ed02e
# ╟─11464714-e089-4e53-b6c0-e61c2526e7ef
# ╟─02a35bd3-3bf2-42bc-8c91-83794391ddc8
# ╟─23498641-43ad-4948-aede-b7151e6c0087
# ╟─aec75d65-2278-402d-ab29-17089ab9f5ed
# ╟─768908ad-6756-476d-9608-5b2b53e38591
# ╟─ca62288b-9b41-47a6-9790-4bce4c016068
# ╟─61690352-55bc-4b86-ab5b-3eaf3efce93f
# ╟─ffeae924-5ff3-4da9-8d0a-d5933dcf6c92
# ╠═cf89067a-2b1b-4995-a5e8-c6d06cbd4004
# ╟─fb7859a8-c671-4517-a5fa-4eccdaeb22eb
# ╠═7c79031a-d566-4844-a6c6-177e7709b571
# ╠═c5344c65-43d7-4b22-84c6-60ff7cb09e1d
# ╠═ef4745a6-c8bf-4a8a-a8fb-9083134be316
# ╠═2dc2202d-9e23-4be4-8a67-c8abe850d51e
# ╠═f8b8f8bd-cfa4-4915-acb1-bc918bb511e1
# ╟─0677112e-8c77-4373-a0db-707bb1915c5e
# ╠═ba43dad8-9db6-4ef2-9f0e-8899006a658f
# ╟─280ec4d5-e352-49e5-8bab-9ea93ec1c628
# ╠═cf28c22f-96ad-4e50-a2ff-d6ba1df6e448
# ╠═02e7d340-6b44-4066-853a-27bc03604d26
# ╠═74a904f5-0194-4688-b074-74f1d738c62d
# ╠═e762ac49-5667-4759-b8b2-a2a25c3a80f9
# ╠═9c8586a5-67da-4d7e-8192-d89cc9370f40
# ╠═1f9e1123-de47-4e07-8e47-0595a10358fd
# ╟─0366e3c4-5ec5-4c96-afdf-17407136a691
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
