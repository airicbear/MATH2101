### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ 2eed39f0-10ca-11ec-174b-09d44d195252
let
	using LinearAlgebra
	using RowEchelon
	using Primes
	
	md"# Linear Algebra In-Class Notes"
end

# ╔═╡ 53bf3df6-2a96-4933-a289-3e492d76b68a
md"""
### Dot product properties

For vectors $\vec{x} = \begin{bmatrix} x_1 \\ ⋮ \\ x_n \end{bmatrix}$ and $\vec{y} = \begin{bmatrix} y_1 \\ ⋮ \\ y_n \end{bmatrix}$ and similarly for vectors $\vec{u}$ and $\vec{v}$ and scalar $a$, the following properties are true:

$\vec{x} ⋅ \vec{y} = x_1 y_1 + ⋯ + x_n y_n = \vec{y} ⋅ \vec{x}$

$\begin{cases}
\vec{x} ⋅ (\vec{u} + \vec{v}) = \vec{x} ⋅ \vec{u} + \vec{x} ⋅ \vec{v} \\
\vec{x} ⋅ (a \vec{u}) = a (x ⋅ \vec{u})
\end{cases}$

$\vec{x} ⋅ \vec{x} ≥ 0 \quad \& \quad \vec{x} ⋅ \vec{x} = 0 \iff \vec{x} = 0$
"""

# ╔═╡ be4d3fdb-e53e-4549-a89a-43c9f3b35844
md"""
### Definition

The **normal** of the vector $\vec{x}$ is

$\|\vec{x}\| = \sqrt{\vec{x} ⋅ \vec{x}}$

##### Observation

$\|\vec{x}\| = 0 \iff \vec{x} = \vec{0}$
"""

# ╔═╡ beed4376-cf5b-487d-9fa8-4da86888eba9
md"""
### Definition

The **distance** between two vectors is

$d(\vec{x},\vec{y}) = \|\vec{x} - \vec{y}\| = \sqrt{(\vec{x} - \vec{y}) ⋅ (\vec{x} - \vec{y})}$

##### Observations

We want the distance to satisfy the following:

$d(\vec{x},\vec{x}) = 0 \text{ and } d(\vec{x},\vec{y}) = 0 \iff \vec{x} = \vec{y}$

$d(\vec{x},\vec{y}) = d(\vec{y},\vec{x})$

$d(\vec{x},\vec{y}) ≥ 0$

$d(\vec{x},\vec{y}) ≤ d(\vec{x},\vec{z}) + d(\vec{z},\vec{y})$
"""

# ╔═╡ fb035399-b22e-4afa-aeab-957825008466
md"""
### Notes

$\|\vec{y} - \vec{x}\| = d(\vec{x},\vec{y})$

$|\vec{x} ⋅ \vec{y}| ≤ \|\vec{x}\| \|\vec{y}\|$

$\begin{align*}
\|x - y\|^2 &= (x - y) ⋅ (x - y) \\
&= (x - y) ⋅ x - (x - y) ⋅ y \\
&= x ⋅ x - y ⋅ x - x ⋅ y + y ⋅ y \\
&≤ \|x\|^2 + 2|x⋅y| + \|y\|^2 \\
&≤ \|x\|^2 + 2 \|x\| \|y\| + \|y\|^2 \\
&≤ (\|x\| + \|y\|)^2
\end{align*}$

$\|x - y\| ≤ \|x\| + \|y\|$
"""

# ╔═╡ aca0cead-df66-4fee-a78f-54fcef15fc36
md"""
### Cauchy-Schwarz inequality

$\|\vec{x} - \vec{y}\| ≤ \|\vec{x}\| \|\vec{y}\|$

The Triangle inequality holds.

##### Proof

If $\vec{y} = \vec{0}$, then $\vec{x} ⋅ \vec{y} = 0$.

$|\vec{x} ⋅ \vec{y}| = 0 ≤ \|\vec{x}\| ⋅ \|\vec{y}\| = 0$

Assume $\vec{y} \neq \vec{0}$

Consider the function

$P_{\vec{x},\vec{y}} : ℝ → ℝ$

$\begin{align*}
P_{\vec{x},\vec{y}}(t) &= \|\vec{x} + \vec{y}\|^2 \\
&= (\vec{x} + t\vec{y})(\vec{x} + t\vec{y}) \\
&= \vec{x} ⋅ \vec{x} + \vec{x} ⋅ \vec{y} + t \vec{y} ⋅ \vec{x} + t^2 ⋅ \vec{y} ⋅ \vec{y} \\
&= \|\vec{x}\|^2 + 2t ⋅ \vec{x} ⋅ \vec{y} + t^2 ⋅ \|\vec{y}\|^2
\end{align*}$

Recall Bhaskara's formula:

If $a ∈ ℝ$, $a \neq 0$, then $p(t) = at^2 + bt + c$ has its zeros given by

$t_± = \frac{-b ± \sqrt{b^2 - 4ac}}{2a}$

By the function $P_{\vec{x},\vec{y}}$ and Bhaskara's formula, we get that

$a = \|\vec{y}\|^2$

$b = 2\vec{x}\vec{y}$

$c = \|\vec{x}\|^2$

Because implies $P_{\vec{x},\vec{y}} (t) ≥ 0$, we get from

$b^2 - 4ac ≤ 0 \implies b^2 ≤ 4ac \implies 4(x⋅y)^2 ≤ 4\|\vec{x}\|^2 ≤4\|\vec{x}\|^2\|\vec{y}\|^2$

$\implies |x ⋅ y| ≤ \|\vec{x}\| \|\vec{y}\|$
"""

# ╔═╡ fb9fbd9b-0cfa-4a7c-80e2-5f0ad23758cb
md"""
### Notes

Let us explore $P_{\vec{x},\vec{y}} (t)$ whenever $\vec{y} ≠ \vec{0}$

${P_{\vec{x},\vec{y}}}'(t) = 2\vec{x} ⋅ \vec{y} + 2t\|\vec{y}\|^2$

${P_{\vec{x},\vec{y}}}''(t) = 2\|\vec{y}\|^2 > 0$

The global minimum is attained provided

${P_{\vec{x},\vec{y}}}'(t_{\min}) = 0$

$2\vec{x}⋅\vec{y} + 2t_{\min} \|\vec{y}\|^2 = 0$

$t_{\min} = \frac{-\vec{x} ⋅ \vec{y}}{\|\vec{y}\|^2}$

Our computation says that the vector $t_{\min} ⋅ \vec{y}$ is such that

$d(\vec{x},-t_{\min}\vec{y}) = \|\vec{x} + t_{\min} \vec{y}\| = \sqrt{P_{\vec{x},\vec{y}}(t_{\min})} ≤ d(\vec{x}, t\vec{y}), \quad ∀t ∈ ℝ$
"""

# ╔═╡ b644db36-0781-46bf-b09e-e35362ce6510
md"""
### Definition

Given $\vec{x}$, $\vec{y}$; $\; \vec{y} ≠ \vec{0}$.
Then

$\text{comp}_{\vec{y}}(\vec{x}) = \frac{\vec{x} ⋅ \vec{y}}{\|\vec{y}\|^2}$

$\text{proj}_{\vec{y}}(\vec{x}) = \text{comp}_{\vec{y}}(\vec{x}) ⋅ \vec{y} = \frac{(\vec{x} ⋅ \vec{y})}{\|y\|^2} \vec{y}$

If $\vec{x} ⋅ \vec{y} = 0$ then the projection is $\vec{0}$ and we say that $\vec{x}$ and $\vec{y}$ are orthogonal.
The vector $\vec{0}$ is orthogonal to every vector.
"""

# ╔═╡ fb4afe70-d51e-41d3-aeaf-d6a76a17aca1
md"""
### Example

$A = \begin{bmatrix}
-1 & -1 & 2 & -1 \\
-1 & 3 & -1 & -1 \\
-2 & 3 & 0 & -2
\end{bmatrix}$

$A^t = \begin{bmatrix}
-1 & -1 & -2 \\
-1 & 3 & 3 \\
2 & -1 & 0 \\
-1 & -1 & -2
\end{bmatrix}$

$AA^t = \begin{bmatrix}
7 & -3 & 1 \\
-3 & 12 & 13 \\
1 & 13 & 17
\end{bmatrix}$

$(AA^t)^{-1} = \begin{bmatrix}
\frac{35}{2} & 32 & -\frac{51}{2} \\
32 & 59 & -47 \\
-\frac{51}{2} & -47 & -\frac{75}{2}
\end{bmatrix} = Q$

$B = A^t Q$

$AB = AA^t Q = AA^t (AA^t)^{-1} = I$
"""

# ╔═╡ 168c44b5-0616-4b84-b2d7-18f979caa48b
let
	A = [-1 -1  2 -1
		 -1  3 -1 -1
		 -2  3  0 -2]
	
	Q = inv(A * A')
	
	B = A' * Q
	AB = A * A' * Q
end

# ╔═╡ 68921ad6-100a-408b-9d33-4e590cf748ca
md"""
### Note

$\begin{align*}
\left\langle \begin{bmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix}, \begin{bmatrix} y_1 \\ y_2 \end{bmatrix} \right\rangle
&= \left\langle \begin{bmatrix} a_{11} x_1 + a_{12} x_2 + a_{13} x_3 \\ a_{21} x_1 + a_{22} x_2 + a_{23} x_3 \end{bmatrix}, \begin{bmatrix} y_1 \\ y_2 \end{bmatrix} \right\rangle \\
&= (a_{11} x_1 + a_{12} x_2 + a_{13} x_3) y_1 + (a_{21} x_1 + a_{22} x_2 + a_{23} x_3) y_2 \\
&= x_1 (a_{11} y_1 + a_{21} y_2) + x_2 (a_{12} y_1 + a_{22} y_2) + x_3 (a_{13} y_1 + a_{23} y_2) \\
&= \left\langle \begin{bmatrix} \end{bmatrix} \right\rangle
\end{align*}$
"""

# ╔═╡ cf50bb73-d2ec-4f96-9707-4dbc6eb0f4cc
md"""

$Ax = b$

``A`` has rank $m$

For all $b$ there is $x$ such that $Ax = b$

1) Suppose $y_0$ such that $\langle Ax, y_0 \rangle = 0$ for all $x$ then $y_0 = 0$. (Why? Let $x_0$ solve $Ax_0 = y_0$. With these

$0 = \langle Ax_0, y_0 \rangle = \langle y_0, y_0 \rangle = \| y_0 \|^2,\; y_0 y_0 = 0$

To say that $\langle Ax, y_0\rangle = 0$ for all $x$ is the same as
$\langle x, A^k y_0 \rangle$ for all $x$ implies $y_0 = 0$
"""

# ╔═╡ 5826d5ba-0594-4799-b548-a2d1df4581cb
md"""
If you can always solve $Ax = b$ (for any $b$) then there is $B$ such that $AB = I$
"""

# ╔═╡ dbea1a16-bf27-4799-88e0-0c57852e1011
md"""
$AB = I$

$I = CA$

$B = (CA)B = C(AB) = CI = C$
"""

# ╔═╡ aa3d3c00-7437-477d-90b9-eee3e283d7a7
md"""
### Example

Example of no solution.
"""

# ╔═╡ 15f7bf64-3c3d-4697-b6e5-7b8c1988d967
let
	A = [-1 -1  2 -1
		 -1  3 -1 -1
		 -2  3  0 -2]
	
	b = [-1
		  0
		  0
		  1]
	
	rref([A' b])
end

# ╔═╡ 6394ab68-aa67-4500-8bf1-24cc44a7ebba
md"""
### Definition

``S`` is a set of elements of $ℝ^n$.

``\text{span}(S) = \text{set of all linear combinations of elements of } S.``
"""

# ╔═╡ fca1ed81-e9fb-4ad6-8c28-ce05cf65e55b
md"""
### Example

$S = \{\vec{u}_1, …, \vec{u}_m\} \qquad \vec{u}_1, …, \vec{u}_m \text{ in } ℝ^m$

$\text{span} = \begin{cases} \vec{u}, & \text{where }\vec{u} = a_1 \vec{u}_1 + ⋯ + a_n \vec{u}_n \text{ for } a_1 … a_n ∈ ℝ\end{cases}$
"""

# ╔═╡ 64e33045-3e3e-41c1-af6e-2f8a7cfaf293
md"""
$V = \text{span}(S)$

1. ``0 ∈ V``

2. ``\vec{v}, \vec{w}`` in $V$ then $\vec{v} + \vec{w}$ in $V$
"""

# ╔═╡ ffff9c05-9548-4668-9ec0-216c8a083628
md"""
$C = \begin{bmatrix} a & b \\ c & d \end{bmatrix}$

$Δ = ad - bc ≠ 0$

$C^{-1} = \frac{1}{Δ} \begin{bmatrix} d & -b \\ -c & a \end{bmatrix}$
"""

# ╔═╡ 27cdc9e1-d3ad-4b89-8d5c-800cda9c218c
md"""
$\text{rg}(T) = \text{span of the columns of } A \text{ is a vector subspace of } ℝ^n$
Also: the vector $\vec{x}$ is $Ax = 0$ is a vector subspace of $ℝ^n$
"""

# ╔═╡ a0fa2bf7-5918-4590-b34b-70fa34f0c629
md"""
$v_1, …, v_s \qquad \text{vectors in } ℝ^n$

If $u_1, …, u_r$ in $\text{span}(v_1, …, v_s)$ and $r > s$, then $u_1, …, u_r$ is linearly dependent.

``V = \text{span}(v_1, …, v_s)`` and $v_1, …, v_s$ is independent.
$V = \text{span}(u_1, …, u_r)$, also $u_1, …, u_r$ is independent.
Then $s = r$.

If not, then, any $s < r$.
Then $u_1, …, u_r$ is linearly dependent.
But that's not true.

$r < s$

Likewise $r < s$ cannot be true, so $s = r$

If $\{u_1, …, u_s\}$ span $V$ and are linearly independent, then $\{u_1, …, u_r\}$ is called a basis, and $r$ is called the *dimension* of $V$ ($\dim{V} = r$)
"""

# ╔═╡ 24e93099-2d37-4869-8ba7-6b88a3a3d58b
md"""
Let $u_1, …, u_n$ be vectors in $ℝ^n$.
Let $v_1, …, v_r$ be linear combinations of the $u_i$.
If $r > s$, then $v_1, …, v_r$ are linearly independent.
(There are $x_1, …, x_r$, not all 0, such that $x_1 v_1 + x_2 v_2 + ⋯ x_r v_r = 0$)

Why can I find $x_j$ not all zero such that $x_1 v_1 + ⋯ x_s v_s = 0$

$\begin{align*}
v_1 &= a_{11} u_1 + a_{21} u_2 + ⋯ + a_{s1} u_s \\
v_2 &= a_{12} u_1 + a_{22} u_2 + ⋯ + a_{s2} u_s \\
&\;\;⋮ \\
v_r &= a_{1r} u_1 + a_{2r} u_2 + ⋯ + a_{sr} u_s \\
\end{align*}$

$\begin{align*}
x_1 v_1 + ⋯ x_r v_r &= (x_1 a_{11} + x_2 a_{12} + ⋯ + x_r a_{1r}) u_1 \\
&+ (x_1 a_{11} + x_2 a_{22} + ⋯ + x_r a_{2r}) u_2 \\
&+ ⋯ \\
&+ (x_1 a_{s1} + x_2 a_{s2} + ⋯ + x_r a_{sr}) u_s
\end{align*}$
"""

# ╔═╡ 21d5f462-2958-4b10-8e67-a5b68596e490
md"""
### Example (similar to Homework #2)

$V = \{a_0 + a_1 t + a_2 t^2 : a_1, a_2, a_3 ∈ ℝ\}$

$\begin{align*}
t^0 && t && t^2 \\
⋮ && ⋮ && ⋮ \\
v_0 && v_1 && v_2
\end{align*}$

$T : V → ℝ^3$

$T(v) = \begin{bmatrix} v(-3) \\ v(1) \\ v(2) \end{bmatrix}$

$\begin{align*}
T(t^0) &= \begin{bmatrix} 1 \\ 1 \\ 1 \end{bmatrix} & T(t^1) &= \begin{bmatrix} -3 \\ 1 \\ 2 \end{bmatrix} & T(t^2) &= \begin{bmatrix} 9 \\ 1 \\ 4 \end{bmatrix} \\
&= 𝐞_1 + 𝐞_2 + 𝐞_3 &&= -3𝐞_1 + 𝐞_2 + 2𝐞_3 &&= 9𝐞_1 + 𝐞_2 + 4𝐞_3
\end{align*}$

$𝐞_1 = \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix} \qquad 𝐞_2 = \begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix} \qquad 𝐞_3 = \begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix}$

Matrix of $T$ is the basis of $v_1, v_2, v_3$ and $𝐞_1, 𝐞_2, 𝐞_3$ of $ℝ^3$ is:

$A = \begin{bmatrix} 1 & -3 & 9 \\ 1 & 1 & 1 \\ 1 & 2 & 4 \end{bmatrix}$
"""

# ╔═╡ f1b82a35-9de1-49bd-808f-4c27a8830af3
md"""
$\begin{bmatrix}
1 & -3 & 9 & 1 & 0 & 0 \\
1 & 1 & 1 & 0 & 1 & 0 \\
1 & 2 & 4 & 0 & 0 & 4
\end{bmatrix} ≃ \begin{bmatrix}
1 & 0 & 0 & \frac{1}{10} & \frac{3}{2} & -\frac{3}{5} \\
0 & 1 & 0 & -\frac{3}{20} & -\frac{1}{4} & \frac{2}{5} \\
0 & 0 & 1 & \frac{1}{20} & -\frac{1}{4} & \frac{1}{5}
\end{bmatrix}$

$B = A^{-1} = \begin{bmatrix}
\frac{1}{10} & \frac{3}{2} & -\frac{3}{5} \\
-\frac{3}{20} & -\frac{1}{4} & \frac{2}{5} \\
\frac{1}{20} & -\frac{1}{4} & \frac{1}{5}
\end{bmatrix}$

$S : ℝ^3 → V$

$S(𝐞_1) = \frac{1}{10} v_0 - \frac{3}{20} v_1 + \frac{1}{20} v_2$

$S(𝐞_2) = \frac{3}{2} v_0 - \frac{1}{4} v_1 - \frac{1}{4} v_2$

$S(𝐞_3) = -\frac{3}{5} v_0 + \frac{2}{5} v_1 + \frac{1}{5} v_2$
"""

# ╔═╡ ed51563d-1d21-4d66-a0dc-8dce5d011055
md"## Determinants"

# ╔═╡ a2b18ba2-f602-4c55-8b6a-9b9e58e0dcd8
let
	A = [1  3
		 -2 5]
	B = [5 -3
		 2 1]
	
	A * B / det(A)
end

# ╔═╡ a51849d1-a595-410e-b121-fbdfbd2e5ed0
let
	A = [2 -1 1
		 1  2 0
		 -3 2 1]
	
	E1 = [0 1 0
		  1 0 0
		  0 0 1] # => (* -1)
	
	E2 = [ 1 0 0
		  -2 1 0
		   0 0 1]
	
	E3 = [1 0 0
		  0 1 0
		  3 0 1]
	
	E4 = [1 0    0
	      0 1    0
		  0 8//5 1]
	
	E5 = [1 0 0
		  0 1 0
		  0 0 5//13] # => (* 5/13)
	
	E6 = [1 0 0
		  0 1 -1
		  0 0 1]
	
	E7 = [1 0     0
		  0 -1//5 0
		  0 0     1] # => (* -1/5)
	
	E8 = [1 -2 0
		  0  1 0
		  0  0 1]
	
	Int.(E8 * E7 * E6 * E5 * E4 * E3 * E2 * E1 * A)
	
	detA = 1 / (-1 * (5/13) * (-1/5))
end

# ╔═╡ 60df75bd-993c-447f-8c56-40f10ffbf771
md"""
### Eigenvalues and eigenvectors

Let $A$ be an $n × n$ matrix.
Let $\lambda$ be a number.
If you can find a nonzero vector $v$ such that $A𝐯 = \lambda 𝐯$ then we say that $\lambda$ is an eigenvalue of $A$ and $𝐯$ an eigenvector.

$T : V → V$

With some basis of $V$, get the matrix of $T$, $A$ and $A$ happens to have a basis of eigenvectors

$\sum_j q_{ij} w_j$

$T(u_1) = \lambda_1 w_1$

$(A - \lambda I)𝐯 = 0$


"""

# ╔═╡ c160a8f5-a0d5-4366-8673-6d1b214cdaa7
md"""
### Example

$A = \begin{bmatrix} -13 & 5 \\ -25 & 17 \end{bmatrix}$

$p_A(\lambda) = \lambda^2 - 4\lambda - 96 = (\lambda - 12)(\lambda + 8)$

$(A - 12I) \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = 0$

$\begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} \frac{1}{5} \\ 1 \end{bmatrix}$

$(A + 8I) \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = 0$

$\begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} 1 \\ 1 \end{bmatrix}$

Eigenvectors corresponding to different eigenvalues are linearly independent.

A same matrix $\lambda_1, \lambda_2, …, \lambda_k$ different eigenvalues $v_1, …, v_n$ are eigenvectors

If $x_1 v_1 + ⋯ x_n v_n = 0$ then $x_1 = ⋯ = x_n = 0$
"""

# ╔═╡ dfe7d42d-b9e8-4ac9-8b72-d1c132b28123
md"""
### Example

$A = \begin{bmatrix} 0 & 0 & 1 & 1 \\ 1 & -3 & -9 & 7 \\ 1 & 0 & 0 & 1 \\ -2 & -4 & -6 & 9 \end{bmatrix}$

$A - \lambda I = \begin{bmatrix} -\lambda & 0 & 1 & 1 \\ 1 & -3 - \lambda & -9 & 7 \\ 1 & 0 & -\lambda & 1 \\ -2 & -4 & -6 & 9 - \lambda \end{bmatrix}$

$\begin{align*}
&\det(A - \lambda I) \\&= -\lambda \begin{vmatrix} -3 - \lambda & -9 & 7 \\ 0 & -\lambda & 1 \\ -4 & -6 & 9 - \lambda \end{vmatrix} + \begin{vmatrix} 1 & -3 - \lambda & 7 \\ 1 & 0 & 1 \\ -2 & -4 & 9 - \lambda \end{vmatrix} + \begin{vmatrix} 1 & -3 - \lambda & -9 \\ 1 & 0 & -\lambda \\ -2 & -4 & -6 \end{vmatrix} \\
&= -\lambda(-\lambda^3 + 6\lambda^2 - 7\lambda + 18) + (9 + 81 - \lambda^2) - (-2\lambda^2 - 16\lambda + 18) \\
&= \lambda^4 - 6\lambda^3 + 8\lambda^2 + 6\lambda - 9
\end{align*}$

$\implies \lambda_1 = -1, \lambda_2 = 1, \lambda_3 = 3$
"""

# ╔═╡ ab4bd76b-421d-4231-9c40-de8d77432c85
let
	A = [0 0 1 1; 1 -3 -9 7; 1 0 0 1; -2 -4 -6 9]
	[rref(A + I), rref(A - I), rref(A - 3I)]
end

# ╔═╡ b219682a-6940-4306-9226-7e5c848c4f3d
md"""
UPenn lecture on generalized eigenvectors:

$A = \begin{bmatrix} 1 & 1 \\ 0 & 1 \end{bmatrix}$

$\begin{align*}
\det(A - \lambda I) &= (1 - \lambda)^2 \implies \lambda = 1  \\
(A - \lambda I)𝐯 = 𝟎 &\implies \begin{bmatrix} 0 & 1 & 0 \\ 0 & 0 & 0 \end{bmatrix} \\
&\implies 𝐯 = \begin{bmatrix} 1 \\ 0 \end{bmatrix}
\end{align*}$
"""

# ╔═╡ a1addfcc-23eb-41ed-860c-622269b78488
md"""
UPenn lecture on generalized eigenvectors:

$A = \begin{bmatrix} 1 & 2 & 0 \\ 1 & 1 & 2 \\ 0 & -1 & 1 \end{bmatrix}$

$A - \lambda I = \begin{bmatrix} 1 - \lambda & 2 & 0 \\ 1 & 1 - \lambda & 2 \\ 0 & -1 & 1 - \lambda \end{bmatrix}$

$\begin{align*}
\det(A - \lambda I) &= (1 - \lambda) \begin{vmatrix} 1 - \lambda & 2 \\ -1 & 1 - \lambda \end{vmatrix} - 2 \begin{vmatrix} 1 & 2 \\ 0 & 1 - \lambda \end{vmatrix} \\
&= (1 - \lambda) \left[(1 - \lambda)(1 - \lambda) + 2\right] - 2 (1 - \lambda) \\
&= (1 - \lambda) \left(1 - 2 \lambda + \lambda^2 + 2\right) - 2 + 2 \lambda \\
&= 1 - 2\lambda + \lambda^2 + 2 - \lambda + 2\lambda^2 - \lambda^3 - 2\lambda - 2 + 2 \lambda \\
&= 1 - 3\lambda + 3\lambda^2 - \lambda^3 \\
&= -(\lambda^3 - 3\lambda^2 + 3\lambda - 1) \\
&= -(\lambda^2 - 2\lambda + 1)(\lambda - 1) \\
&= -(\lambda - 1)^3 \\
&\implies \lambda = 1
\end{align*}$

$\begin{array}{c | c c c c}
1 & 1 & -3 & 3 & -1 \\
&  & 1 & -2 & 1 \\
\hline
& 1 & -2 & 1 & 0
\end{array}$
"""

# ╔═╡ 69702313-8e8c-420f-985c-4e960eff5c0e
md"""
$\begin{align*}
(A - \lambda I)𝐯_1 = 𝟎 &\implies \begin{bmatrix} 0 & 2 & 0 & 0 \\ 1 & 0 & 2 & 0 \\ 0 & -1 & 0 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & 0 & 2 & 0 \\ 0 & 2 & 0 & 0 \\ 0 & -1 & 0 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & 0 & 2 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{bmatrix} \\
&\implies 𝐯_1 = \begin{bmatrix} -2 \\ 0 \\ 1 \end{bmatrix}
\end{align*}$

$\begin{align*}
(A - \lambda I)𝐯_2 = 𝐯_1 &\implies \begin{bmatrix} 0 & 2 & 0 & -2 \\ 1 & 0 & 2 & 0 \\ 0 & -1 & 0 & 1 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & 0 & 2 & 0 \\ 0 & 2 & 0 & -2 \\ 0 & -1 & 0 & 1 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & 0 & 2 & 0 \\ 0 & 1 & 0 & -1 \\ 0 & 0 & 0 & 0 \end{bmatrix} \\
&\implies 𝐯_2 = \begin{bmatrix} -2v_3 \\ 0 \\ v_3 \end{bmatrix} + \begin{bmatrix} 0 \\ -1 \\ 0 \end{bmatrix} = \begin{bmatrix} 0 \\ -1 \\ 0 \end{bmatrix} \text{ for } v_3 = 0
\end{align*}$

$\begin{align*}
(A - \lambda I)𝐯_3 = 𝐯_2 &\implies \begin{bmatrix} 0 & 2 & 0 & 0 \\ 1 & 0 & 2 & -1 \\ 0 & -1 & 0 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & 0 & 2 & -1 \\ 0 & 2 & 0 & 0 \\ 0 & -1 & 0 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & 0 & 2 & -1 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{bmatrix} \\
&\implies 𝐯_2 = \begin{bmatrix} -1-2v_3 \\ 0 \\ v_3 \end{bmatrix} = \begin{bmatrix} -1 \\ 0 \\ 0 \end{bmatrix} \text{ for } v_3 = 0
\end{align*}$
"""

# ╔═╡ 6a21e355-79d1-4cf9-a301-39f48603e3be
md"""
### Definition

For a linear operator $T : V → V$, a nonzero vector $𝐯$ satisfying $(A - \lambda I)^k 𝐯 = 𝟎$ for some positive integer $k$ and some scalar $\lambda$ is called a **generalized eigenvector** of $T$.
"""

# ╔═╡ ca03d7bf-e557-44fd-a6e0-440ffdf0ad38
md"""
### Example

Show that $𝐯 = \begin{bmatrix} 4 \\ 1 \end{bmatrix}$ is a generalized 2-eigenvector for $A = \begin{bmatrix} 1 & -1 \\ 1 & 3 \end{bmatrix}$ that is not a (regular) 2-eigenvector.

- We compute $(A - 2I)𝐯 = \begin{bmatrix} 1 & 1 \\ -1 & -1 \end{bmatrix} \begin{bmatrix} 4 \\ 1 \end{bmatrix} = \begin{bmatrix} 5 \\ -5 \end{bmatrix}$, and since this is not zero, $𝐯$ is not a 2-eigenvector.

- However, $(A - 2I)^2 𝐯 = \begin{bmatrix} 1 & 1 \\ -1 & -1 \end{bmatrix} \begin{bmatrix} 5 \\ -5 \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \end{bmatrix}$, and so $𝐯$ is a generalized 2-eigenvector, with $k = 2$.
"""

# ╔═╡ a4e857e7-c39e-4a9c-b675-62774642033d
md"""
### Proposition (Generalized Eigenspaces)

For a linear operator $T : V → V$, the set of vectors $𝐯$ satisfying $(T - \lambda I)^k 𝐯 = 𝟎$ for some positive integer $k$ is a subspace of $V$.
This subspace is called the **generalized $\lambda$-eigenspace** of $T$.
"""

# ╔═╡ cbcc37dd-c176-4d25-802e-313031f8403e
md"""
### Proposition (Eigenvalues for Generalized Eigenvectors)

If $T : V → V$ is a linear operator and $𝐯$ is a nonzero vector satisfying $(T - \lambda I)^k 𝐯 = 𝟎$ for some positive integer $k$ and some scalar $\lambda$, then $\lambda$ is an eigenvalue of $T$.
Furthermore, the eigenvalue associated to a generalized eigenvector is unique.
"""

# ╔═╡ 09ba27d4-576c-448c-9b2f-d9162323ea72
md"""
### Theorem (Computing Generalized Eigenspaces)

If $T : V → V$ is a linear operator and $V$ is finite-dimensional, then the generalized $\lambda$-eigenspace of $T$ is equal to $\ker(T - \lambda I)^{\dim(V)}$.
In other words, if $(T - \lambda I)^k 𝐯 = 𝟎$ for some positive integer $k$, then in fact $(T - \lambda I)^{\dim(V)} 𝐯 = 𝟎$.
"""

# ╔═╡ 0e05d296-64c6-4fd1-8395-893d16954061
md"""
### Example

Find the generalized eigenspaces of $A = \begin{bmatrix} 2 & 0 & 0 \\ -1 & 2 & 1 \\ 1 & -1 & 0 \end{bmatrix}$.

- The characteristic polynomial is $\det(\lambda I - A) = (\lambda - 1)^2 (\lambda - 2)$ so the eigenvalues are $\lambda = 1, 1, 2$.

- For the generalized 1-eigenspace, we must compute the nullspace of $(A - I)^3 = \begin{bmatrix} 1 & 0 & 0 \\ -1 & 0 & 0 \\ 1 & 0 & 0 \end{bmatrix}$.

- Upon row-reducing, we see that the generalized 1-eigenspace has dimension 2 and is spanned by (0, 1, 0) and (0, 0, 1).

- Note here that neither of the generalized 1-eigenvectors is a 1-eigenvector, and (in fact) the 1-eigenspace of $A$ is only 1-dimensional.
  This means $A$ is not diagonalizable.
"""

# ╔═╡ c8915158-e05a-4125-a911-214a145bf0dd
md"""
### Example

Find the generalized eigenspaces of $A = \begin{bmatrix} 2 & 0 & 0 \\ -1 & 2 & 1 \\ 1 & -1 & 0 \end{bmatrix}$.

- For the generalized 2-eigenspace, we must compute the nullspace of $(A - 2I)^3 = \begin{bmatrix} 0 & 0 & 0 \\ -1 & 2 & 3 \\ 1 & -3 & -4 \end{bmatrix}$.

- Upon row-reducing, we see that the generalized 2-eigenspace has dimension 1 and is spanned by (1, -1, 1).

In the example, observe that $V$ does not have a basis of eigenvectors of $A$ since 1-eigenspace is only 1-dimensional.

Nonetheless, $V$ does possess a basis of *generalized* eigenvectors.
"""

# ╔═╡ 8ee8ea63-578a-4992-be46-dee8c64b606c
md"""
### Theorem (Independence of Generalized Eigenvectors)

If $𝐯_1, 𝐯_2, …, 𝐯_n$ are generalized eigenvectors of $T$ associated to distinct eigenvalues $\lambda_1, \lambda_2, …, \lambda_n$, then $𝐯_1, 𝐯_2, …, 𝐯_n$ are linearly independent.
"""

# ╔═╡ 33623f80-2cce-46ac-83c5-ab7b2b7e02b5
md"""
### Theorem (Upper-Triangular Associated Matrix)

Suppose $T : V → V$ is a linear operator on a finite-dimensional vector space such that the scalar field of $V$ contains all eigenvalues of $T$.
If $\lambda$ is an eigenvalue of $T$ having multiplicity $d$, then there exists a basis $\beta$ of $V$ such that $[T]_\beta^\beta$ is upper-triangular and the last $d$ entries on the diagonal are equal to $\lambda$.
"""

# ╔═╡ 2ef13edf-2c14-422c-a044-e064b1cb98eb
md"""
### Theorem (Dimension of Generalized Eigenspace)

If $V$ is finite-dimensional, $T : V → V$ is linear, and $\lambda$ is a scalar, then the dimension of the generalized $\lambda$-eigenspace is equal to the multiplicity $d$ of $\lambda$ as a root of the characteristic polynomial of $T$, and in fact the generalized $\lambda$-eigenspace is the kernel of $(T - \lambda I)^d$.
"""

# ╔═╡ c9e7ba69-15e2-40fe-9f9c-501ac8e378bb
md"""
### Example

Suppose the characteristic polynomial of $T$ is $p(t) = t^3(t - 2)^2$.
Then:

- The generalized 0-eigenspace has dimension 3 and is $\ker(T^3)$.

- The generalized 2-eigenspace has dimension 2 and is $\ker(T - 2I)^2$.
"""

# ╔═╡ 72ff60d6-5825-4829-92b2-e9b7822a3204
let
	A = [0 0 1 0; 0 2 -3 1; 0 1 -2 1; 0 0 -1 1]
	rref(I - A)
end

# ╔═╡ 8b93c8e1-34c9-45aa-a8c6-f2a86417b632
md"""
### Example

Find the dimensions of the generalized eigenspaces of $A = \begin{bmatrix} 0 & 0 & 1 & 0 \\ 0 & 2 & -3 & 1 \\ 0 & 1 & -2 & 1 \\ 0 & 0 & -1 & 1 \end{bmatrix}$, and then verify the result by finding a basis for each generalized eigenspace.
Also, decide whether or not $A$ is diagonalizable.

- Some computation produces $\det(\lambda I - A) = \lambda^3 (\lambda - 1)$.

- By the Theorem, the dimension of the generalized 0-eigenspace is 3 and the dimension of the generalized 1-eigenspace is 1.

- For the generalized 0-eigenspace, the nullspace of $A^3 = \begin{bmatrix} 0 & 0 & 0 & 0 \\ 0 & 1 & -1 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & -1 & 1 & 0 \end{bmatrix}$ has basis $\begin{bmatrix} 1 \\ 0 \\ 0 \\ 0 \end{bmatrix}, \begin{bmatrix} 0 \\ 1 \\ 1 \\ 0 \end{bmatrix}, \begin{bmatrix} 0 \\ 0 \\ 0 \\ 1 \end{bmatrix}$.

- For the generalized 1-eigenspace, the nullspace of $I - A = \begin{bmatrix} 1 & 0 & -1 & 1 \\ 0 & -1 & 3 & -1 \\ 0 & -1 & 3 & -1 \\ 0 & 0 & 1 & 0 \end{bmatrix}$ has basis vector $\begin{bmatrix} 0 \\ 1 \\ 0 \\ -1 \end{bmatrix}$.

- The matrix of $A$ is not diagonalizable because there is not a basis of (regular) eigenvectors, as the 0-eigenspace only has dimension 1 (it is spanned by (1, 0, 0, 0), as can be seen by row-reducing $A$).
"""

# ╔═╡ b4dada21-2a39-4800-bab3-00830e7fe1f4
md"""
### Theorem (Spectral Decomposition)

If $V$ is finite-dimensional, $T : V → V$ is linear, and all eigenvalues of $T$ lie in the scalar field of $V$, then $V$ has a basis of generalized eigenvectors of $T$.
"""

# ╔═╡ 248fadd2-7321-48fa-b24c-78d4c1f65862
md"""
### PS24 #2

$A = \begin{bmatrix} 4 & 1 & 1 \\ -1 & 6 & 2 \\ 0 & 0 & 5 \end{bmatrix}$

$A - \lambda I = \begin{bmatrix} 4 - \lambda & 1 & 1 \\ -1 & 6 - \lambda & 2 \\ 0 & 0 & 5 - \lambda \end{bmatrix}$

$\begin{align*}
\det(A - \lambda I) &= (5 - \lambda) \left[(4 - \lambda)(6 - \lambda) - (1)(-1)\right] \\
&= (5 - \lambda)(\lambda - 10 \lambda + 25) \\
&= (5 - \lambda)(\lambda - 5)^2 \\
&= -(\lambda - 5)^3
\end{align*}$
"""

# ╔═╡ 4e494279-f03c-4e55-a0e8-732c4134f245
let
	A = [4 1 1; -1 6 2; 0 0 5]
	nullspace((A - 5I)^3)
end

# ╔═╡ 9281bcee-8528-41bb-bce2-04de75f6d6a4
md"""
### Example

$P = \{a_0 t^0 + a_1 t^1 + a_2 t^2 : a_0, a_1, a_2 ∈ ℝ\}$

$p, q ∈ P$

$⟨p, q⟩ = \frac{1}{2} ∫_{-1}^1 p(t) q(t) \;dt$
"""

# ╔═╡ 2e0ec81b-d661-4510-81d1-dd9dfbf77435
md"""
### Example

$⟨p, q⟩' = a_0$

$p = a_0 t^0 + a_1 b^1 + a_1 t^2$

$q = b_0 t^0 + b_1 t^1 + b_2 t^2$

$\left.p\left(\frac{d}{dt}\right)\right|_{t = 0} = \left(a_0 \frac{d^0}{dt^0} + a_1 \frac{d^1}{dt^1} + a_2 \frac{d^2}{dt^2}\right)$

$\begin{align*}
p\left(\frac{d}{dt}\right) &= \left(a_0 \frac{d^0}{dt^0} + a_1 \frac{d^1}{dt^1} + a_2 \frac{d^2}{dt^2}\right)_{t = 0} \\
&= a_0 ⋅ 0 + a_1 ⋅ 1 + a_2 ⋅ 0 ⋅ 1 = a_1
\end{align*}$
"""

# ╔═╡ 7397679e-dc18-4e5e-abce-4f6e1d3f5b55
md"""
$V = \{ f : [0,x] → ℝ : f \text{ is differentiable to any order} \}$

$T : V → V$

$T(f) = -f''$

$⟨f, g⟩ = ∫_0^π f(x) g(x) \;dx$

$T(f) = \lambda f$

$-f'' = \lambda f$

$f'' + \lambda f = 0$
"""

# ╔═╡ 2361bffd-658f-418c-bf90-2a562a525ce8
md"""
### Symmetric Transformations

$⟨Tv, w⟩ = ⟨v, Tw⟩ \text{ for all } v, w$
"""

# ╔═╡ 952dca44-a6df-4377-adf7-464fde9a1fcc
md"""
### Example

$f, g ∈ V$

$\begin{align*}
∫_0^π -f''(x) g(x) \;dx &= ∫_0^π - (f')' g(x) \;dx \\
&= -\left[f' g\;\Big|_0^π - ∫_0^π f'(x) g'(x)\right] \\
&= ∫_0^π f'(x) g'(x) \;dx
\end{align*}$

$\begin{align*}
∫_0^π -g''(x) f(x) \;dx &= ∫_0^π - (g')' f(x) \;dx \\
&= -\left[g' f\;\Big|_0^π - ∫_0^π g'(x) f'(x)\right] \\
&= ∫_0^π f'(x) g'(x) \;dx
\end{align*}$

So:

$⟨Tf, g⟩ = ⟨f, Tg⟩$

$⟨Tf, f⟩ = ⟨f', f'⟩ = \|f\|^2 ≥ 0$
"""

# ╔═╡ 2abf9ece-4419-4f00-88fa-3cc1b2b80deb
md"""
$V = \{f : [0,\pi] → ℝ : f \text{ is } ∞ \text{ diff, } f(0) = f(\pi)\}$

$⟨f, g⟩ = ∫_0^\pi f(t) g(t) \;dt$

$T : V → V$

$T(f) = -f''$

$⟨T(f), g⟩ = ∫_0^\pi f'(t) g'(t) \;dt$

$⟨f, T(g)⟩ = ∫_0^\pi f'(t) g'(t) \;dt$

$∫_a^b u' v \;dt = \left[uv\right]_a^b - ∫_a^b u v' \;dt$

$f_k(t) = \sin(kt)\; \quad k = 1,2, …$

$T(f_k) = k^2 f_k$

$⟨f_k, f_ℓ⟩ = 0\;\; \text{ if } k ≠ ℓ$

$\cos{(\alpha + \beta)} = \cos{\alpha} \cos{\beta} - \sin{\alpha} \sin{\beta}$

$\cos{(\alpha - \beta)} = \cos{\alpha} \cos{\beta} + \sin{\alpha} \sin{\beta}$

$\cos{(\alpha - \beta)} - \cos{(\alpha + \beta)} = 2(\sin{\alpha} \sin{\beta})$

$\begin{align*}
⟨f_k, f_ℓ⟩ &= ∫_0^\pi \sin(kt) \sin(ℓt) \;dt \\
&= \frac{1}{2} ∫_0^\pi \left[\cos{((k - ℓ) t)} - \cos{((k - ℓ) t)}\right] \;dt \\
&= \frac{1}{2} \left[\frac{1}{k - ℓ} \sin{((k - ℓ) t)} - \frac{1}{k - ℓ} \sin{((k - ℓ) t)}\right]_0^\pi \\
&= 0
\end{align*}$
"""

# ╔═╡ 6787b7a5-3200-4005-9e98-d2a125a858ed
md"""
### Example

``V`` inner product space, finite dimension

``T : V → V``

``⟨Tu, v⟩ = ⟨u, Tv⟩`` for all ``u,v`` in ``V``

``v_1,…,v_n`` ``w`` ``A`` vectors of ``T``

get ``p_A(\lambda)`` let ``\lambda_1`` be a root,

Find an eigenvector for ``\lambda_1``

``Aw_1 = \lambda_1 w_1`` where ``w_1`` in ``ℝ^n``

``w_1 = \begin{bmatrix} \alpha_{11} \\ \alpha_{12} \\ ⋮ \\ \alpha_{1n} \end{bmatrix}``

``u_1 = \alpha_{11} v_1 + \alpha_{12} v_2 + ⋯ + \alpha_{1n} v_n``

``T(u_1) = \lambda_1 u_1``

``V_1 = \{v ∈ V : v ⟂ u_1\}``

If ``v ∈ V_1``, then ``T(v) ∈ V_1``

``v ∈ V_1`` ``⟨v, u_1⟩ = 0``

``⟨Tv, u_1⟩ = ⟨v, Tu_1⟩ = ⟨v, \lambda_1 u_1⟩ = \lambda ⟨v, u_1⟩ = \lambda_1 0 = 0``
"""

# ╔═╡ 955a2db6-d31a-4bb3-a9e0-476c7359189a
md"""
### Example

$\begin{align*}
x &= x_1 v_1 + ⋯ + x_n v_n \\
y &= y_1 v_1 + ⋯ + y_n v_n
\end{align*}$

$\begin{align*}
⟨x, y⟩ &= ⟨x_1 v_1 + x_2 v_2 + x_3 v_3, y_1 v_1 + y_2 v_2 + y_3 v_3⟩ \\
&= x_1 ⟨v_1, y⟩ + x_2 ⟨v_2, y⟩ + x_3 ⟨v_3, y⟩ \\
&= x_1 ⟨v_1, v_1⟩ y_1 + ⋯ + x_3 ⟨v_3, v_3⟩ y_3 \\
&= \begin{bmatrix} x_1 & x_2 & x_3 \end{bmatrix}
\begin{bmatrix}
⟨v_1, v_1⟩ & ⟨v_1, v_2⟩ & ⟨v_1, v_3⟩ \\
⟨v_2, v_1⟩ & ⟨v_2, v_2⟩ & ⟨v_2, v_3⟩ \\
⟨v_3, v_1⟩ & ⟨v_3, v_2⟩ & ⟨v_3, v_3⟩
\end{bmatrix}
\end{align*}$

$\begin{align*}
⟨v_1, v_1⟩ &= \lambda_1 \\
⟨v_1, v_2⟩ = ⟨v_2, v_1⟩ &= 0 \\
⟨v_2, v_2⟩ &= \lambda_2 \\
⟨v_3, v_3⟩ &= \lambda_3
\end{align*}$
"""

# ╔═╡ 362a937a-74b5-4e99-9c69-af3657a0b2a8
md"""
### Example

``V`` polynomials of degree ≤ 2

``⟨p, q⟩ ∫_0^1 p(t) q(t) \;dt``

``t^0, t^1, t^2``
"""

# ╔═╡ 8612b1d3-c521-471e-b90b-b610467ea7ed
md"""
### PS29 #22

$(A - \lambda_2 I) v = (\lambda_1 - \lambda_2) v$

$(A - \lambda_2 I) (A - \lambda_2 I) v = (\lambda_1 - \lambda_2)^2 v$

In general,

$(A - \lambda_2 I)^{m_2} v = (\lambda_1 - \lambda_2)^{m_2} v$

Can it be that $(A - \lambda_2)^{m_2} v = 0$?

No because $\lambda_1 ≠ \lambda_2$.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Primes = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
RowEchelon = "af85af4c-bcd5-5d23-b03a-a909639aa875"

[compat]
Primes = "~0.5.0"
RowEchelon = "~0.2.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Primes]]
git-tree-sha1 = "afccf037da52fa596223e5a0e331ff752e0e845c"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.0"

[[RowEchelon]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f479526c4f6efcbf01e7a8f4223d62cfe801c974"
uuid = "af85af4c-bcd5-5d23-b03a-a909639aa875"
version = "0.2.1"
"""

# ╔═╡ Cell order:
# ╟─2eed39f0-10ca-11ec-174b-09d44d195252
# ╟─53bf3df6-2a96-4933-a289-3e492d76b68a
# ╟─be4d3fdb-e53e-4549-a89a-43c9f3b35844
# ╟─beed4376-cf5b-487d-9fa8-4da86888eba9
# ╟─fb035399-b22e-4afa-aeab-957825008466
# ╟─aca0cead-df66-4fee-a78f-54fcef15fc36
# ╟─fb9fbd9b-0cfa-4a7c-80e2-5f0ad23758cb
# ╟─b644db36-0781-46bf-b09e-e35362ce6510
# ╟─fb4afe70-d51e-41d3-aeaf-d6a76a17aca1
# ╠═168c44b5-0616-4b84-b2d7-18f979caa48b
# ╟─68921ad6-100a-408b-9d33-4e590cf748ca
# ╟─cf50bb73-d2ec-4f96-9707-4dbc6eb0f4cc
# ╟─5826d5ba-0594-4799-b548-a2d1df4581cb
# ╟─dbea1a16-bf27-4799-88e0-0c57852e1011
# ╟─aa3d3c00-7437-477d-90b9-eee3e283d7a7
# ╠═15f7bf64-3c3d-4697-b6e5-7b8c1988d967
# ╟─6394ab68-aa67-4500-8bf1-24cc44a7ebba
# ╟─fca1ed81-e9fb-4ad6-8c28-ce05cf65e55b
# ╟─64e33045-3e3e-41c1-af6e-2f8a7cfaf293
# ╟─ffff9c05-9548-4668-9ec0-216c8a083628
# ╟─27cdc9e1-d3ad-4b89-8d5c-800cda9c218c
# ╟─a0fa2bf7-5918-4590-b34b-70fa34f0c629
# ╟─24e93099-2d37-4869-8ba7-6b88a3a3d58b
# ╟─21d5f462-2958-4b10-8e67-a5b68596e490
# ╟─f1b82a35-9de1-49bd-808f-4c27a8830af3
# ╟─ed51563d-1d21-4d66-a0dc-8dce5d011055
# ╠═a2b18ba2-f602-4c55-8b6a-9b9e58e0dcd8
# ╠═a51849d1-a595-410e-b121-fbdfbd2e5ed0
# ╟─60df75bd-993c-447f-8c56-40f10ffbf771
# ╟─c160a8f5-a0d5-4366-8673-6d1b214cdaa7
# ╟─dfe7d42d-b9e8-4ac9-8b72-d1c132b28123
# ╠═ab4bd76b-421d-4231-9c40-de8d77432c85
# ╟─b219682a-6940-4306-9226-7e5c848c4f3d
# ╟─a1addfcc-23eb-41ed-860c-622269b78488
# ╟─69702313-8e8c-420f-985c-4e960eff5c0e
# ╟─6a21e355-79d1-4cf9-a301-39f48603e3be
# ╟─ca03d7bf-e557-44fd-a6e0-440ffdf0ad38
# ╟─a4e857e7-c39e-4a9c-b675-62774642033d
# ╟─cbcc37dd-c176-4d25-802e-313031f8403e
# ╟─09ba27d4-576c-448c-9b2f-d9162323ea72
# ╟─0e05d296-64c6-4fd1-8395-893d16954061
# ╟─c8915158-e05a-4125-a911-214a145bf0dd
# ╟─8ee8ea63-578a-4992-be46-dee8c64b606c
# ╟─33623f80-2cce-46ac-83c5-ab7b2b7e02b5
# ╟─2ef13edf-2c14-422c-a044-e064b1cb98eb
# ╟─c9e7ba69-15e2-40fe-9f9c-501ac8e378bb
# ╠═72ff60d6-5825-4829-92b2-e9b7822a3204
# ╟─8b93c8e1-34c9-45aa-a8c6-f2a86417b632
# ╟─b4dada21-2a39-4800-bab3-00830e7fe1f4
# ╟─248fadd2-7321-48fa-b24c-78d4c1f65862
# ╠═4e494279-f03c-4e55-a0e8-732c4134f245
# ╟─9281bcee-8528-41bb-bce2-04de75f6d6a4
# ╟─2e0ec81b-d661-4510-81d1-dd9dfbf77435
# ╟─7397679e-dc18-4e5e-abce-4f6e1d3f5b55
# ╟─2361bffd-658f-418c-bf90-2a562a525ce8
# ╟─952dca44-a6df-4377-adf7-464fde9a1fcc
# ╟─2abf9ece-4419-4f00-88fa-3cc1b2b80deb
# ╟─6787b7a5-3200-4005-9e98-d2a125a858ed
# ╟─955a2db6-d31a-4bb3-a9e0-476c7359189a
# ╟─362a937a-74b5-4e99-9c69-af3657a0b2a8
# ╟─8612b1d3-c521-471e-b90b-b610467ea7ed
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
