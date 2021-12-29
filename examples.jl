### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ c8e9d359-45b8-47c6-80aa-2eed394d1bb3
begin
	using LinearAlgebra
	using RowEchelon
	
	md"# Linear Algebra Examples"
end

# ╔═╡ c458e6c0-1d48-4fdf-ab1e-16a030dc7cc2
md"""
### Example 8.7: Finding the eigenvalues

Find the eigenvalues of the matrix

$A = \begin{bmatrix} -5 & 2 \\ -7 & 4 \end{bmatrix}.$

**Solution.**
By Theorem 8.6, a scalar $\lambda$ is an eigenvalue of $A$ if and only if $\det(A - \lambda I) = 0$.
We calculate the determinant:

$\begin{align*}
\det(A - \lambda I) &= \begin{vmatrix} -5 - \lambda & 2 \\ -7 & 4 - \lambda \end{vmatrix} \\
&= (-5 - \lambda)(4 - \lambda) + 14 \\
&= \lambda^2 + \lambda - 6
\end{align*}$

Therefore, $\lambda$ is an eigenvalue if and only if $\lambda^2 + \lambda - 6 = 0$.
We can find the roots of this equation using the quadratic formula, or equivalently, by factoring the left-hand side:

$\lambda^2 + \lambda - 6 = 0 \iff (\lambda + 3)(\lambda - 2) = 0.$

Therefore, the eigenvalues are $\lambda = -3$ and $\lambda = 2$.
"""

# ╔═╡ 95eaff6e-30c3-4303-a0c6-15876bf70584
md"""
### Example 8.8: Finding the eigenvalues

Find the eigenvalues of the matrix

$A = \begin{bmatrix} 5 & -4 & 4 \\ 2 & -1 & 2 \\ 0 & 0 & 2 \end{bmatrix}.$

**Solution.**
Once again, we calculate $\det(A - \lambda I)$:

$\begin{align*}
\det(A - \lambda I) &= \begin{vmatrix} 5 - \lambda & -4 & 4 \\ 2 & -1 - \lambda & 2 \\ 0 & 0 & 2 - \lambda \end{vmatrix} \\
&= (5 - \lambda)(-1 - \lambda)(2 - \lambda) - 2(-4)(2 - \lambda) \\
&= -\lambda^3 + 6\lambda^2 - 11\lambda + 6 \\
&= (3 - \lambda)(1 - \lambda)(2 - \lambda)
\end{align*}$

The eigenvalues are the roots of this polynomial, i.e., the solutions of the equation $(\lambda - 3)(\lambda - 1)(2 - \lambda) = 0$.
Therefore, the eigenvalues of $A$ are $\lambda = 1$, $\lambda = 2$, and $\lambda = 3$.
"""

# ╔═╡ 01d0fc2a-fd0a-4eb8-8ff6-d47011eeb8bf
md"""
### Example 8.9: No real eigenvalue

Find the eigenvalues of the matrix

$A = \begin{bmatrix} 0 & -1 \\ 1 & 0 \end{bmatrix}$

**Solution.**
We have

$\det(A - \lambda I) = \begin{vmatrix} -\lambda & -1 \\ 1 & -\lambda \end{vmatrix} = \lambda^2 + 1.$

Since $\lambda^2 + 1 = 0$ does not have any solutions in the real numbers, the matrix $A$ has no real eigenvalues.
(However, if we were working over the field of complex numbers rather than real numbers, this matrix would have eigenvalues $\lambda = ±i$).
"""

# ╔═╡ f808acd5-e84b-4aee-a7b2-5a4a6f243e7a
md"""
### Example 8.15: Eigenvalues of triangular matrix

Find the eigenvalues of

$A = \begin{bmatrix} 1 & 2 & 4 \\ 0 & 4 & 7 \\ 0 & 0 & 6 \end{bmatrix}.$

**Solution:**
We calculate $\det(A - \lambda I) = 0$ as follows:

$\det(A - \lambda I) = \det{\begin{bmatrix} 1 - \lambda & 2 & 4 \\ 0 & 4 - \lambda & 7 \\ 0 & 0 & 6 - \lambda \end{bmatrix}} = (1 - \lambda)(4 - \lambda)(6 - \lambda)$

Solving the equation $(1 - \lambda)(4 - \lambda)(6 - \lambda) = 0$ results in the eigenvalues $\lambda_1 = 1$, $\lambda_2 = 4$, and $\lambda_3 = 6$.
Thus the eigenvalues are the entries on the main diagonal of $A$.
"""

# ╔═╡ b561b2aa-3b1a-11ec-0989-7b5de473f175
md"""
### Example 8.22: Diagonalizing a matrix

Diagonalize the matrix

$A = \begin{bmatrix} 3 & 0 & 2 \\ 6 & 4 & 3 \\ -4 & 0 & -3 \end{bmatrix}.$

In other words, find an invertible matrix $P$ and a diagonal matrix $D$ such that $P^{-1}AP = D$.

**Solution.**
By Theorem 8.2.1, we use the eigenvectors of $A$ as the columns of $P$ and the corresponding eigenvalues as the diagonal entries of $D$.
We already found the eigenvectors and -values of $A$ in Example 8.13.
They were

$𝐯_1 = \begin{bmatrix} -1 \\ 1 \\ 1 \end{bmatrix}, \quad 𝐯_2 = \begin{bmatrix} -1 \\ 0 \\ 2 \end{bmatrix}, \; \text{ and } \; 𝐯_3 = \begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix}$

with corresponding eigenvalues $\lambda_1 = 1$, $\lambda_2 = -1$, and $\lambda_3 = 4$.
Therefore we can use

$P = \begin{bmatrix} -1 & -1 & 0 \\ 1 & 0 & 1 \\ 1 & 2 & 0 \end{bmatrix} \; \text{ and } \; D = \begin{bmatrix} 1 & 0 & 0 \\ 0 & -1 & 0 \\ 0 & 0 & 4 \end{bmatrix}.$

To double-check that $P^{-1}AP$ is indeed equal to $D$, we first compute the inverse of $P$:

$P^{-1} = \begin{bmatrix} -2 & 0 & -1 \\ 1 & 0 & 1 \\ 2 & 1 & 1 \end{bmatrix}.$

Then

$P^{-1}AP = \begin{bmatrix} -2 & 0 & -1 \\ 1 & 0 & 1 \\ 2 & 1 & 1 \end{bmatrix} \begin{bmatrix} 3 & 0 & 2 \\ 6 & 4 & 3 \\ -4 & 0 & -3 \end{bmatrix} \begin{bmatrix} -1 & -1 & 0 \\ 1 & 0 & 1 \\ 1 & 2 & 0 \end{bmatrix} = \begin{bmatrix} 1 & 0 & 0 \\ 0 & -1 & 0 \\ 0 & 0 & 4 \end{bmatrix} = D.$

Alternatively, we could have checked that $AP = PD$, which would not have required computing $P^{-1}$.
"""

# ╔═╡ c2467751-6b99-4d73-8bb7-0578eb225f07
md"""
### Example 8.24: A matrix that cannot be diagonalized

Show that the matrix $A = \begin{bmatrix} 1 & 1 \\ 0 & 1 \end{bmatrix}$ cannot be diagonalized.

**Solution.**
Through the usual procedure, we find that the characteristic polynomial is $(1 - \lambda)^2$, and therefore the only eigenvalue is $\lambda = 1$.
To find the eigenvectors, we solve the equation $(A - I)𝐯 = 𝟎$:

$\left[\begin{array}{cc|c}
0 & 1 & 0 \\ 0 & 0 & 0
\end{array}\right]$

The general solution is

$𝐯 = t\begin{bmatrix} 1 \\ 0 \end{bmatrix}.$

Because the solution space is 1-dimensional, there is only one basic eigenvector:

$𝐯_1 = \begin{bmatrix} 1 \\ 0 \end{bmatrix}.$

Since the matrix $A$ has only one basic eigenvector, we cannot find two linearly independent eigenvectors.
Therefore, by Theorem 8.21, $A$ cannot be diagonalized.
"""

# ╔═╡ bcf89321-2693-4c21-8f57-d86902688918
md"""
### Example 8.25: Raising a matrix to a high power.

Let $A = \begin{bmatrix} 2 & 1 & 0 \\ 0 & 1 & 0 \\ -1 & -1 & 1 \end{bmatrix}.$
Find $A^{50}$.

**Solution:**
First, we will diagonalize $A$.
Following the usual steps, we find that the eigenvalues are $\lambda = 1$ and $\lambda = 2$.
The basic eigenvectors corresponding to $\lambda = 1$ are

$𝐯_1 = \begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix} \; \text{ and } \; 𝐯_2 = \begin{bmatrix} -1 \\ 1 \\ 0 \end{bmatrix}$

and the basic eigenvector corresponding to $\lambda = 2$ is

$𝐯_3 = \begin{bmatrix} -1 \\ 0 \\ 1 \end{bmatrix}$

Now we construct $P$ by using the basic eigenvectors of $A$ as the columns of $P$.
Thus

$P = \begin{bmatrix} 0 & -1 & -1 \\ 0 & 1 & 0 \\ 1 & 0 & 1 \end{bmatrix}$

The inverse of $P$ is

$P^{-1} = \begin{bmatrix} 1 & 1 & 1 \\ 0 & 1 & 0 \\ -1 & -1 & 0 \end{bmatrix}.$

Then

$P^{-1} AP = \begin{bmatrix} 1 & 1 & 1 \\ 0 & 1 & 0 \\ -1 & -1 & 0 \end{bmatrix} \begin{bmatrix} 2 & 1 & 0 \\ 0 & 1 & 0 \\ -1 & -1 & 1 \end{bmatrix} \begin{bmatrix} 0 & -1 & -1 \\ 0 & 1 & 0 \\ 1 & 0 & 1 \end{bmatrix} = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 2 \end{bmatrix} = D.$

Now it follows by rearranging the equation that $A = PDP^{-1}$, and therefore, as noted above,

$\begin{align*}
A^{50} = PD^{50} P^{-1} &= \begin{bmatrix} 0 & -1 & -1 \\ 0 & 1 & 0 \\ 1 & 0 & 1 \end{bmatrix} \begin{bmatrix} 1^{50} & 0 & 0 \\ 0 & 1^{50} & 0 \\ 0 & 0 & 2^{50} \end{bmatrix} \begin{bmatrix} 1 & 1 & 1 \\ 0 & 1 & 0 \\ -1 & -1 & 0 \end{bmatrix} \\
&= \begin{bmatrix} 2^{50} & -1 + 2^{50} & 0 \\ 0 & 1 & 0 \\ 1 - 2^{50} & 1 - 2^{50} & 1 \end{bmatrix}.
\end{align*}$

Thus, through diagonalization, we have efficiently computed a high power of $A$.
The following example shows that we can also use the same technique for finding a square root of a matrix.
"""

# ╔═╡ e66def3f-4560-42ca-a359-3da470459b5c
md"""
### Example 8.41: Algebraic and geometric multiplicity

Let

$A = \begin{bmatrix} 3 & 1 & 0 & 0 & 0 \\ 0 & 3 & 0 & 0 & 0 \\ 0 & 0 & 4 & 0 & 0 \\ 0 & 0 & 0 & 4 & 0 \\ 0 & 0 & 0 & 0 & 5 \end{bmatrix}.$

Find the algebraic and geometric multiplicity of each eigenvalue of $A$.
"""

# ╔═╡ 04dd5528-f671-47d0-a2d9-94a6ff14d581
let
	A = [3 1 0 0 0
	     0 3 0 0 0
	     0 0 4 0 0
	     0 0 0 4 0
		 0 0 0 0 5]

	
	eigvals(A), nullspace(A - 3I), nullspace(A - 4I), nullspace(A - 5I)
end

# ╔═╡ 8d49476d-2523-46c6-aea8-85d8c88c0f3b
md"""
### Example 11.2: ℝⁿ with the usual dot product

``V = ℝ^n``, with the usual dot product

$⟨𝐮, 𝐯⟩ = 𝐮 ⋅ 𝐯 = u_1 v_1 + … + u_n v_n$

is an inner product space.
"""

# ╔═╡ b097dba0-45fc-452d-b144-478f660edb2f
md"""
### Example 11.3: ℝ² with a non-standard inner product.

Let $V = ℝ^2$, and consider the matrix $A = \begin{bmatrix} 1 & 1 \\ 1 & 2 \end{bmatrix}$.
Define an operation $⟨—,—⟩$ by

$⟨𝐮, 𝐯⟩ = 𝐮^T A𝐯 = \begin{bmatrix} u_1 & u_2 \end{bmatrix} \begin{bmatrix} 1 & 1 \\ 1 & 2 \end{bmatrix} \begin{bmatrix} v_1 \\ v_2 \end{bmatrix} = u_1 v_1 + u_1 v_2 + u_2 v_1 + 2u_2 v_2.$

Then $ℝ^2$, with this operation, is an inner product space.
"""

# ╔═╡ c8ee8fb8-ca4e-4ee4-bb19-42db3a0234ef
md"""
### Example 11.4: Continuous functions on an interval

Let $a < b$ be real numbers, and let

$[a,b] = \{x ∣ a ≤ x ≤ b\}$

be the closed interval.
Let $V = C[a,b]$ be the vector space of all continuous functions $f : [a,b] → ℝ$.
Given two functions $f,g ∈ C[a,b]$, define

$⟨f, g⟩ = ∫_a^b f(x) g(x) \; dx.$

This is an inner product.
"""

# ╔═╡ 9e726f96-0b60-4982-83f3-8dd30d697095
md"""
### Example 11.5: Polynomials

Let $𝐏$ be the vector space of polynomials (of any degree) with real coefficients.
Let $a < b$ be real numbers, and consider the inner product defined by

$⟨p, q⟩ = ∫_a^b p(x) q(x) \;dx.$

This is an inner product space.
"""

# ╔═╡ e0e6115c-c057-49dd-9160-fcc092832566
md"""
### Example 11.6: Real Hilbert space

Let $\textbf{Hilb}_ℝ$ be the vector space of all infinite sequences of real numbers $a = (a_0, a_1, a_2, …)$ satisfying

${a_0}^2 + {a_1}^2 + {a_2}^2 + … < ∞.$

These are called the **square summable** sequences.
(One needs to do a little bit of work to show that it is indeed a vector space; in particular, to show that the sum of two square summable sequences is square summable.)
On this space, we can define an inner product as follows:

$⟨a, b⟩ = a_0 b_0 + a_1 b_1 + a_2 b_2 + ….$

The details will be worked out in Exercise 11.1.6.
This inner product space is called **(real) Hilbert space**.
"""

# ╔═╡ c5609164-72e2-475d-8ab9-be9ba6605e75
md"""
### Example 11.10: Norm in C[-1,1]

Calculate the norm of $f(x) = x^2$ in $C[-1, 1]$.

**Solution.**
In the vector space $C[-1,1]$, the inner product is defined as

$⟨f, g⟩ = ∫_{-1}^1 f(x) g(x) \;dx.$

We therefore have

$⟨f, f⟩ = ∫_{-1}^1 f(x)^2 \;dx = ∫_{-1}^1 x^4 \;dx = \left[\frac{1}{5}x^5\right]_{-1}^1 = \frac{2}{5}.$

Therefore, $\|f\| = \sqrt{⟨f, f⟩} = \sqrt{\frac{2}{5}}.$

A vector $𝐮$ in an inner product space is called **normalized** or a **unit vector** if $\|𝐮\| = 1$.
We note that if $𝐯$ is any non-zero vector in an inner product space, then

$𝐮 = \frac{1}{\|𝐯\|} 𝐯$

is normalized.
"""

# ╔═╡ e4394bdd-5946-4c98-83ff-62bb9907649b
md"""
### Example 11.14: Find the angle between two vectors

Find the angle between 1 and $x^2$ in $C[-1,1]$.

**Solution.**
We have

$\begin{align*}
⟨1, 1⟩ &= ∫_{-1}^1 1 ⋅ 1 \;dx = 2, \\
\left\langle x^2, x^2 \right\rangle &= ∫_{-1}^1 x^2 ⋅ x^2 \;dx = \frac{2}{5}, \\
\left\langle 1, x^2 \right\rangle &= ∫_{-1}^1 1 ⋅ x^2 \;dx = \frac{2}{3}.
\end{align*}$

Therefore

$\cos{\theta} = \frac{\left\langle 1, x^2 \right\rangle}{\|1\|\|x^2\|} = \frac{\frac{2}{3}}{\sqrt{2} \sqrt{\frac{2}{5}}} = \frac{\sqrt{5}}{3}.$

The angle $\theta$ is $\cos^{-1}\left(\frac{\sqrt{5}}{3}\right)$, which is approximately 0.7297 radians or 41.81 degrees.
"""

# ╔═╡ 520768e2-89c7-40bb-9550-9011c6e50eb4
md"""
### Example 11.18: Orthogonal complement

Consider the inner product space $𝐏_3$ of polynomials of degree at most 3, with the inner product defined by

$⟨f, g⟩ = ∫_{-1}^1 f(x) g(x) \;dx$

Find the orthogonal complement of $\{x^2\}$.

**Solution.**
We have to compute the set of all polynomials of the form $p(x) = ax^3 + bx^2 + cx + d$ that are orthogonal to $x^2$.
So let us compute the inner product:

$\begin{align*}
⟨p(x), x^2⟩ &= ∫_{-1}^1 (ax^3 + bx^2 + cx + d) x^2 \;dx \\
&= ∫_{-1}^1 ax^5 + bx^4 + cx^3 + dx^2 \;dx \\
&= 0a + \frac{2}{5}b + 0c + \frac{2}{3} d.
\end{align*}$

Setting this equal to 0, we see that $⟨p(x), x^2⟩ = 0$ if and only if $\frac{2}{5} b + \frac{2}{3} d = 0$, or equivalently, $3b + 5d = 0$.
The basic solutions are

$\begin{bmatrix} a \\ b \\ c \\ d \end{bmatrix} = \begin{bmatrix} 1 \\ 0 \\ 0 \\ 0 \end{bmatrix}, \quad \begin{bmatrix} a \\ b \\ c \\ d \end{bmatrix} = \begin{bmatrix} 0 \\ 5 \\ 0 \\ -3 \end{bmatrix}, \quad \begin{bmatrix} a \\ b \\ c \\ d \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \\ 1 \\ 0 \end{bmatrix}$

giving the following basis for the space of polynomials orthogonal to $x^2$:

$\{x^3, 5x^2 - 3, x\}.$
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
RowEchelon = "af85af4c-bcd5-5d23-b03a-a909639aa875"

[compat]
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

[[RowEchelon]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f479526c4f6efcbf01e7a8f4223d62cfe801c974"
uuid = "af85af4c-bcd5-5d23-b03a-a909639aa875"
version = "0.2.1"
"""

# ╔═╡ Cell order:
# ╟─c8e9d359-45b8-47c6-80aa-2eed394d1bb3
# ╟─c458e6c0-1d48-4fdf-ab1e-16a030dc7cc2
# ╟─95eaff6e-30c3-4303-a0c6-15876bf70584
# ╟─01d0fc2a-fd0a-4eb8-8ff6-d47011eeb8bf
# ╟─f808acd5-e84b-4aee-a7b2-5a4a6f243e7a
# ╟─b561b2aa-3b1a-11ec-0989-7b5de473f175
# ╟─c2467751-6b99-4d73-8bb7-0578eb225f07
# ╟─bcf89321-2693-4c21-8f57-d86902688918
# ╟─e66def3f-4560-42ca-a359-3da470459b5c
# ╠═04dd5528-f671-47d0-a2d9-94a6ff14d581
# ╟─8d49476d-2523-46c6-aea8-85d8c88c0f3b
# ╟─b097dba0-45fc-452d-b144-478f660edb2f
# ╟─c8ee8fb8-ca4e-4ee4-bb19-42db3a0234ef
# ╟─9e726f96-0b60-4982-83f3-8dd30d697095
# ╟─e0e6115c-c057-49dd-9160-fcc092832566
# ╟─c5609164-72e2-475d-8ab9-be9ba6605e75
# ╟─e4394bdd-5946-4c98-83ff-62bb9907649b
# ╟─520768e2-89c7-40bb-9550-9011c6e50eb4
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
