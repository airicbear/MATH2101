### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ fba59878-c8c1-4f42-8214-6d99568d8de1
begin
	using LinearAlgebra
	using RowEchelon
	using Polynomials
	using Symbolics
	
	md"# Linear Algebra Homework"
end

# ╔═╡ 89066d76-248e-11ec-3893-757003b33e14
md"## Homework #1"

# ╔═╡ 6fce260e-d916-4ac1-902f-8366d179c9df
md"""
### HW1 #1

Let

$A = \begin{bmatrix} 2 & -1 & 3 & 0 \\ -3 & 1 & 0 & 2 \end{bmatrix}$

Find $B$ such that $AB = I$ where $I = \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix}$.
"""

# ╔═╡ 58d62bd7-a0e7-4746-a2a2-fb09d053899b
let
	A = [ 2 -1 3 0
		 -3  1 0 2]
	
	C = rref([A*A' I])[:,3:end]

	B = rationalize.(A' * C)
end

# ╔═╡ 0fbade02-25e8-417d-b4c7-4ce30449ac62
md"""
$A^t = \begin{bmatrix} 2 & -3 \\ -1 & 1 \\ 3 & 0 \\ 0 & 2 \end{bmatrix}$

$AA^t = \begin{bmatrix} 4 + 1 + 9 & -6 - 1 \\ -6 - 1 & 9 + 1 + 4 \end{bmatrix}
= \begin{bmatrix} 14 & -7 \\ -7 & 14 \end{bmatrix}$

$\begin{bmatrix} AA^t & I \end{bmatrix} = \begin{bmatrix} 14 & -7 & 1 & 0 \\ -7 & 14 & 0 & 1 \end{bmatrix}$

$\begin{align*}
R_1 \gets \frac{1}{14} R_1, \; R_2 \gets -\frac{1}{7} R_2 &\implies \begin{bmatrix} 1 & -\frac{1}{2} & \frac{1}{14} & 0 \\ 1 & -2 & 0 & -\frac{1}{7} \end{bmatrix} \\
R_2 \gets R_2 - R_1 &\implies \begin{bmatrix} 1 & -\frac{1}{2} & \frac{1}{14} & 0 \\ 0 & -\frac{3}{2} & -\frac{1}{14} & -\frac{1}{7} \end{bmatrix} \\
R_2 \gets -\frac{2}{3} R_2 &\implies \begin{bmatrix} 1 & -\frac{1}{2} & \frac{1}{14} & 0 \\ 0 & 1 & \frac{1}{21} & \frac{2}{21} \end{bmatrix} \\
R_1 \gets R_1 + \frac{1}{2} R_2 &\implies \begin{bmatrix} 1 & 0 & \frac{2}{21} & \frac{1}{21} \\ 0 & 1 & \frac{1}{21} & \frac{2}{21} \end{bmatrix}
\end{align*}$

$C = \begin{bmatrix} \frac{2}{21} & \frac{1}{21} \\ \frac{1}{21} & \frac{2}{21} \end{bmatrix}$
"""

# ╔═╡ c6fb02f8-d4e1-42e3-b7c0-3a59402f1dbb
md"""
$\begin{align*}
B = A^t C &= \begin{bmatrix} 2 & -3 \\ -1 & 1 \\ 3 & 0 \\ 0 & 2 \end{bmatrix} \begin{bmatrix} \frac{2}{21} & \frac{1}{21} \\ \frac{1}{21} & \frac{2}{21} \end{bmatrix} \\
&= \begin{bmatrix} \frac{1}{21} & -\frac{4}{21} \\ -\frac{1}{21} & \frac{1}{21} \\ \frac{2}{7} & \frac{1}{7} \\ \frac{2}{21} & \frac{4}{21} \end{bmatrix}
\end{align*}$
"""

# ╔═╡ 13ca6959-333b-42ac-b879-566be5ef2615
md"""
### HW1 #2

Let

$A = \begin{bmatrix} 1 & -2 & 0 & 1 \\ 3 & -2 & 4 & -1 \\ -1 & 0 & -2 & 1 \end{bmatrix}$

The columns of $A$ span a vector subspace $V$ of $ℝ^3$.

(a) Find the dimension of $V$

Let

$W = \left\{ \vec{y} = \begin{bmatrix} y_1 \\ y_2 \\ y_3 \end{bmatrix} : A^t y = 0 \right\}.$

This is another subspace of $ℝ^3$

(b) Find the dimension of $W$.

(c) Let $\vec{z}$ be an element of $V$, let $\vec{y}$ be an element of $W$.
Compute

$\vec{z} ⋅ \vec{y}$

Hint: For (c) use that if $\vec{v}_1, …, \vec{v}_4$ are the columns of $A$, then $\vec{z} = x_1 \vec{v}_1 + ⋯ + x_4 \vec{v}_4$ for some numbers $x_1, …, x_4$.
"""

# ╔═╡ f8eb8ef6-b8a6-4b4f-b8eb-3cc3125796e5
let
	A = [ 1 -2  0  1
		  3 -2  4 -1
		 -1  0 -2  1]
	
	rref(A)
end

# ╔═╡ d5523a14-f49e-4a66-a009-d84578179de7
md"(a) $\text{dim}(V) = 2$"

# ╔═╡ c86dce76-be3a-44d5-ad62-5ae915349403
let
	A = [ 1 -2  0  1
		  3 -2  4 -1
		 -1  0 -2  1]
	
	rref(A')
end

# ╔═╡ c87c49d4-1139-46cb-bf43-17d6fd09c993
md"(b) $\dim(W) = 2$"

# ╔═╡ bb4f86cd-b207-4478-a7b9-e4ff7146fb8c
md"""
### HW1 #3

Let $\vec{w}_1, \vec{w}_2, \vec{w}_3$ be linearly independent vectors in $ℝ^n$, let

$\begin{align*}
\vec{v}_1 &= 2\vec{w}_2 - 2\vec{w}_3 \\
\vec{v}_2 &= -\vec{w}_1 + 2\vec{w}_2 - 3\vec{w}_3 \\
\vec{v}_3 &= -3\vec{w}_1 - \vec{w}_2
\end{align*}$

Are $\vec{v}_1, \vec{v}_2, \vec{v}_3$ linearly independent?
Explain.
"""

# ╔═╡ dc73982a-8dab-4173-8c70-e180404e7aec
let
	A = [ 0  2  2
		 -1  2 -3
		 -3 -1  0]
	
	rref(A)
end

# ╔═╡ 370ad756-c609-460e-a9f3-2d5305e594c6
md"Yes. The rank is equal to the number of variables (r = 3)."

# ╔═╡ fb839ada-ae51-41dd-bf95-45a1cc60d6f0
md"""
### HW1 #4

Let

$\vec{u}_1 = \begin{bmatrix} 1 & 3 & -2 & 2 \end{bmatrix},
\vec{u}_2 = \begin{bmatrix} -1 & -3 & 5 & 0 \end{bmatrix},
\vec{u}_3 = \begin{bmatrix} 2 & 6 & 2 & 8 \end{bmatrix}$

These row vectors (i.e., $M^{1 × 4}$ matrices) span a vector space $U$.
Find a basis of $U$.
Given your answer, what is the dimension of $U$?
"""

# ╔═╡ 774cfe46-d41b-406a-ab40-8da0707ffe09
let
	u1 = [1 3 -2 2]
	u2 = [-1 -3 5 0]
	u3 = [2 6 2 8]
	
	M = [u1; u2; u3]
	
	rationalize.(rref(M))
end

# ╔═╡ f3fddbe8-67c4-46e2-802c-5bf292f69f11
md"""
A basis of $U$ is $\left\{\begin{bmatrix} 1 & 3 & 0 & \frac{10}{3} \end{bmatrix}, \begin{bmatrix} 0 & 0 & 1 & \frac{2}{3} \end{bmatrix}\right\}$.
The dimension of $U$ is 2.
"""

# ╔═╡ 34582ac8-c4e4-483e-ba3c-b7b83c44d199
md"""
### Homework #2

Let $V$ be the set of polynomials of degree at most two in the variable $t$:

$V = \{a_0 + a_1 t + a_2 t^2 : a_0, a_1, a_2 ∈ ℝ\}$

This is the set of linear combinations of $t^0$, $t^1$, and $t^2$ (the term $a_0$ is really $a_0 t^0$).
Define $T : V → ℝ$ by:

$T(v) = \begin{bmatrix} v(-1) \\ v(0) \\ v(1) \end{bmatrix},$

This is a linear operator. (You should convince yourself this is true.)

(a) Write the matrix of $T$ with respect to the bases

$\vec{v}_1 = t^0, \quad \vec{v}_1 = t, \quad \vec{v}_1 = t^2$

of $V$ and

$𝐞_1 = \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix}, \quad 𝐞_2 = \begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix}, \quad 𝐞_3 = \begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix}$

of $ℝ^3$.
Call this matrix $A$

(b) The matrix $A$ is invertible, its inverse is

$B = \begin{bmatrix} 0 & 1 & 0 \\ -\frac{1}{2} & 0 & \frac{1}{2} \\ \frac{1}{2} & -1 & \frac{1}{2} \end{bmatrix}.$

This matrix corresponds to an operator $S : ℝ^3 → V$,

$S(\vec{r}) = \alpha_0(r_1, r_2, r_3) t^0 + \alpha_1(r_1, r_2, r_3) t + \alpha_2(r_1, r_2, r_3) t^2.$

Write down a formula for the coefficients $\alpha_j$.

(c) Find $T(S(\vec{r}))$ and $S(T(v))$ when $v = a_0 + a_1 t + a_2 t^2$.
"""

# ╔═╡ aba0a023-be9b-4f41-8728-6d162ef8a558
md"""
**(a)**

$v(t) = a_0 + a_1 t + a_2 t^2$

$v(-1) = a_0 - a_1 + a_2$

$v(0) = a_0$

$v(1) = a_0 + a_1 + a_2$

$T(v) = \begin{bmatrix} v(-1) \\ v(0) \\ v(1) \end{bmatrix} = \begin{bmatrix} a_0 - a_1 + a_2 \\ a_0 \\ a_0 + a_1 + a_2 \end{bmatrix}$

$T(𝐞_1) = \begin{bmatrix} 1 - 0 + 0 \\ 1 \\ 1 + 0 + 0 \end{bmatrix} = \begin{bmatrix} 1 \\ 1 \\ 1 \end{bmatrix}$

$T(𝐞_2) = \begin{bmatrix} 0 - 1 + 0 \\ 0 \\ 0 + 1 + 0 \end{bmatrix} = \begin{bmatrix} -1 \\ 0 \\ 1 \end{bmatrix}$

$T(𝐞_3) = \begin{bmatrix} 0 - 0 + 1 \\ 0 \\ 0 + 0 + 1 \end{bmatrix} = \begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix}$

$A = \begin{bmatrix}
1 & -1 & 1 \\
1 & 0 & 0 \\
1 & 1 & 1
\end{bmatrix}$
"""

# ╔═╡ cdd75ddd-39e2-41da-9036-e0ff29d72aca
md"""
**(b)**

$B\vec{r} = \begin{bmatrix} 0 & 1 & 0 \\ -\frac{1}{2} & 0 & \frac{1}{2} \\ \frac{1}{2} & -1 & \frac{1}{2} \end{bmatrix} \begin{bmatrix} r_1 \\ r_2 \\ r_3 \end{bmatrix} = \begin{bmatrix} r_2 \\ -\frac{1}{2} r_1 + \frac{1}{2} r_3 \\ \frac{1}{2} r_1 - r_2 + \frac{1}{2} r_3 \end{bmatrix}$

So

$\begin{align*}
\alpha_0(r_1, r_2, r_3) &= r_2 \\
\alpha_1(r_1, r_2, r_3) &= -\frac{1}{2} r_1 + \frac{1}{2} r_3 \\
\alpha_2(r_1, r_2, r_3) &= \frac{1}{2} r_1 - r_2 + \frac{1}{2} r_3
\end{align*}$

In general,

$\alpha_j(r_1, r_2, r_3) = B_{j,1} + B_{j,2} + B_{j,3}$
"""

# ╔═╡ 6cff14e0-548f-47eb-a34f-7f373b9efe63
md"""
**(c)**

$\begin{align*}
T(S(\vec{r})) &= \begin{bmatrix} 1 & -1 & 1 \\ 1 & 0 & 0 \\ 1 & 1 & 1 \end{bmatrix} \left(\begin{bmatrix} 0 & 1 & 0 \\ -\frac{1}{2} & 0 & \frac{1}{2} \\ \frac{1}{2} & -1 & \frac{1}{2} \end{bmatrix} \begin{bmatrix} r_1 \\ r_2 \\ r_3 \end{bmatrix}\right) \\
&= \begin{bmatrix} 1 & -1 & 1 \\ 1 & 0 & 0 \\ 1 & 1 & 1 \end{bmatrix} \begin{bmatrix} r_2 \\ -\frac{1}{2} r_1 + \frac{1}{2} r_3 \\ \frac{1}{2} r_1 - r_2 + \frac{1}{2} r_3 \end{bmatrix} \\
&= \begin{bmatrix}
r_2 + \frac{1}{2} r_1 - \frac{1}{2} r_3 + \frac{1}{2} r_1 - r_2 + \frac{1}{2} r_3 \\
r_2 \\
r_2 - \frac{1}{2} r_1 + \frac{1}{2} r_3 + \frac{1}{2} r_1 - r_2 + \frac{1}{2} r_3
\end{bmatrix} \\
&= \begin{bmatrix} r_1 \\ r_2 \\ r_3 \end{bmatrix}
\end{align*}$

$\begin{align*}
S(T(v)) &= \begin{bmatrix} 0 & 1 & 0 \\ -\frac{1}{2} & 0 & \frac{1}{2} \\ \frac{1}{2} & -1 & \frac{1}{2} \end{bmatrix} \begin{bmatrix} a_0 - a_1 + a_2 \\ a_0 \\ a_0 + a_1 + a_2 \end{bmatrix} \\
&= \begin{bmatrix}
a_0 \\
-\frac{1}{2} a_0 + \frac{1}{2} a_1 - \frac{1}{2} a_2 + \frac{1}{2} a_0 + \frac{1}{2} a_1 + \frac{1}{2} a_2 \\
\frac{1}{2} a_0 - \frac{1}{2} a_1 + \frac{1}{2} a_2 - a_0 + \frac{1}{2} a_0 + \frac{1}{2} a_1 + \frac{1}{2} a_2
\end{bmatrix} \\
&= \begin{bmatrix} a_0 \\ a_1 \\ a_2 \end{bmatrix}
\end{align*}$
"""

# ╔═╡ c42f0ad9-08bb-42ce-8b31-63ce79651a79
md"## Homework #3"

# ╔═╡ a7b4e7de-1765-4e9e-bf26-bd5678334717
md"""
### HW3 #1

The matrix

$A' = \begin{bmatrix}
1 & -19 & 0 & -3 \\
2 & -4 & 2 & 0 \\
2 & 6 & 5 & 1 \\
5 & 6 & 3 & 4
\end{bmatrix}$

was obtained from the matrix

$A = \begin{bmatrix}
1 & -1 & 0 & -3 \\
2 & 2 & 2 & 0 \\
2 & 7 & 5 & 1 \\
5 & 1 & 3 & 4
\end{bmatrix}$

by replacing the second *column* of $A$ by $A_2 - 3A_1 + 5A_4$ ($A_i$ is the $i$-th column of $A$).
Given that $\det(A) = 20$, find $\det(A')$ using the properties listed on the sidebar.
"""

# ╔═╡ 24dae75e-90f1-4d64-968b-addac182a164
md"""
$A' = \begin{bmatrix} A_1 & A_2 - 3A_1 + 5A_4 & A_3 & A_4 \end{bmatrix}$

$\det(A')$

$\begin{align*}
&= \det(\begin{bmatrix} A_1 & A_2 - 3A_1 + 5A_4 & A_3 & A_4 \end{bmatrix} ) \\
&= -\det(\begin{bmatrix} A_2 - 3A_1 + 5A_4 & A_1 & A_3 & A_4 \end{bmatrix}) \\
&= -\det(\begin{bmatrix} A_2 & A_1 & A_3 & A_4 \end{bmatrix}) - \det(\begin{bmatrix} -3A_1 & A_1 & A_3 & A_4 \end{bmatrix}) - \det(\begin{bmatrix} 5A_4 & A_1 & A_3 & A_4 \end{bmatrix}) \\
&= \det(\begin{bmatrix} A_1 & A_2 & A_3 & A_4 \end{bmatrix}) + 3 \det(\begin{bmatrix} A_1 & A_1 & A_3 & A_4 \end{bmatrix}) - 5 \det(\begin{bmatrix} A_4 & A_1 & A_3 & A_4 \end{bmatrix}) \\
&= \det(A) + 3(0) - 5(0) \\
&= \det(A) \\
&= 20
\end{align*}$
"""

# ╔═╡ e50a0ec8-4336-43fe-9ebe-6fb82fefbf89
md"""
### HW3 #2

Let

$A = \begin{bmatrix} 1 & 3 & x_1 \\ 2 & 4 & x_2 \\ -1 & 0 & x_3 \end{bmatrix}$

(a) Find the solutions of $\det(A) = 0$ (the solutions are column vectors).

(b) Verify that every solution in (a) is a linear combination of the first two columns of $A$.
"""

# ╔═╡ 8c06c619-f214-42da-8fe2-19b2ec9de635
md"""
**(a)**

$\begin{align*}
\det(A) &= x_1 \begin{vmatrix} 2 & 4 \\ -1 & 0 \end{vmatrix} - x_2 \begin{vmatrix} 1 & 3 \\ -1 & 0 \end{vmatrix} + x_3 \begin{vmatrix} 1 & 3 \\ 2 & 4 \end{vmatrix} \\
&= x_1 (0 + 4) - x_2 (0 + 3) + x_3 (4 - 6) \\
&= 4x_1 - 3x_2 - 2x_3 = 0
\end{align*}$

$\begin{align*}
x_1 &= \frac{3}{4}x_2 + \frac{1}{2}x_3 \\
x_2 &= x_2 \\
x_3 &= x_3
\end{align*}$

$\begin{align*}
x_2 &= s\\
x_3 &= t
\end{align*}$

$\begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} = \begin{bmatrix}  \frac{3}{4} s + \frac{1}{2} t \\ s \\ t \end{bmatrix} = s \begin{bmatrix} \frac{3}{4} \\ 1 \\ 0 \end{bmatrix} + t \begin{bmatrix} \frac{1}{2} \\ 0 \\ 1 \end{bmatrix}$
"""

# ╔═╡ 348ea2aa-aef7-4464-92f0-9c7317552dda
md"""
**(b)**

We want to show that the span of the first two columns of $A$ is equal to the basis for $\text{null}(A)$:

$\left\{\begin{bmatrix} 1 \\ 2 \\ -1 \end{bmatrix}, \begin{bmatrix} 3 \\ 4 \\ 0 \end{bmatrix}\right\} = \left\{\begin{bmatrix} \frac{3}{4} \\ 1 \\ 0 \end{bmatrix}, \begin{bmatrix} \frac{1}{2} \\ 0 \\ 1 \end{bmatrix}\right\}$

Let $\vec{v}_1, \vec{v}_2$ represent the first two columns and $\vec{u}_1, \vec{u}_2$ represent the columns in the basis of the nullspace of $A$.

$\left[\begin{array}{cc|cc}
1 & 3 & \frac{3}{4} & \frac{1}{2} \\
2 & 4 & 1 & 0 \\
-1 & 0 & 0 & 1
\end{array}\right]
\leadsto
\left[\begin{array}{cc|cc}
1 & 3 & \frac{3}{4} & \frac{1}{2} \\
1 & 2 & \frac{1}{2} & 0 \\
0 & 0 & \frac{3}{4} & \frac{3}{2}
\end{array}\right]
\leadsto
\left[\begin{array}{cc|cc}
1 & 3 & \frac{3}{4} & \frac{1}{2} \\
0 & 1 & \frac{1}{4} & \frac{1}{2} \\
0 & 0 & \frac{3}{4} & \frac{3}{2}
\end{array}\right]$
$\leadsto
\left[\begin{array}{cc|cc}
1 & 3 & \frac{3}{4} & \frac{1}{2} \\
0 & 1 & \frac{1}{4} & \frac{1}{2} \\
0 & 0 & 0 & 0
\end{array}\right]
\leadsto
\left[\begin{array}{cc|cc}
1 & 0 & 0 & -1 \\
0 & 1 & \frac{1}{4} & \frac{1}{2} \\
0 & 0 & 0 & 0
\end{array}\right]$

$\begin{align*}
&\implies \frac{1}{4} \vec{v}_2 = \vec{u}_1 \\
&\implies -\vec{v}_1 + \frac{1}{2} \vec{v}_2 = \vec{u}_2
\end{align*}$
"""

# ╔═╡ 06f5fec0-595b-4d29-8530-0175155683c4
md"## Homework #4"

# ╔═╡ 9bdb9179-b5e4-4dbf-9277-97c5b670e81d
md"""
### HW4 #1

Let

$𝐮 = \begin{bmatrix} \frac{\sqrt{3}}{2 \sqrt{2}} \\ \frac{1}{\sqrt{2}} \\ \frac{1}{2 \sqrt{2}} \end{bmatrix}$

Let $T(𝐯) = \text{proj}_𝐮(𝐯)$.
Write down the matrix of $T$ with respect to the canonical basis of $ℝ^3$, then find all its eigenvalues and a basis for each of its eigenspaces.
"""

# ╔═╡ 0d502224-dba6-4d35-a0a5-11ce171bbc02
md"""
$\text{proj}_𝐮\left(𝐯\right) = \frac{𝐮 ⋅ 𝐯}{\|𝐮\|^2} 𝐮$

$𝐮 ⋅ 𝐯 = \frac{\sqrt{3}}{2 \sqrt{2}} v_1 + \frac{1}{\sqrt{2}} v_2 + \frac{1}{2 \sqrt{2}} v_3$

$\|𝐮\|^2 = \left(\frac{\sqrt{3}}{2\sqrt{2}}\right)^2 + \left(\frac{1}{\sqrt{2}}\right)^2 + \left(\frac{1}{2\sqrt{2}}\right)^2 = \frac{3}{8} + \frac{4}{8} + \frac{1}{8} = 1$

$\begin{align*}
T(𝐯) = \text{proj}_𝐮(𝐯) &= \begin{bmatrix} \left(\frac{\sqrt{3}}{2 \sqrt{2}} v_1 + \frac{1}{\sqrt{2}} v_2 + \frac{1}{2 \sqrt{2}} v_3\right) \frac{\sqrt{3}}{2 \sqrt{2}} \\ \left(\frac{\sqrt{3}}{2 \sqrt{2}} v_1 + \frac{1}{\sqrt{2}} v_2 + \frac{1}{2 \sqrt{2}} v_3\right) \frac{1}{\sqrt{2}} \\ \left(\frac{\sqrt{3}}{2 \sqrt{2}} v_1 + \frac{1}{\sqrt{2}} v_2 + \frac{1}{2 \sqrt{2}} v_3\right) \frac{1}{2 \sqrt{2}} \end{bmatrix} \\
&= \begin{bmatrix} \frac{3}{8} v_1 + \frac{\sqrt{3}}{4} v_2 + \frac{\sqrt{3}}{8} v_3 \\
\frac{\sqrt{3}}{4} v_1 + \frac{1}{2} v_2 + \frac{1}{4} v_3 \\ \frac{\sqrt{3}}{8} v_1 + \frac{1}{4} v_2 + \frac{1}{8} v_3 \end{bmatrix}
\end{align*}$
"""

# ╔═╡ 30dd5e24-eec8-4b5f-bd13-50fe1c3162c9
let
	u = [sqrt(3) / 2sqrt(2)
		 1 / sqrt(2)
		 1 / 2sqrt(2)]
	
	[u*u[1] u*u[2] u*u[3]]
end

# ╔═╡ fb879466-d6b8-4200-86fc-ce95b58d7794
md"""
$T(𝐞_1) = \begin{bmatrix} \frac{3}{8} \\ \frac{\sqrt{3}}{4} \\ \frac{\sqrt{3}}{8} \end{bmatrix}
\quad T(𝐞_2)  = \begin{bmatrix} \frac{\sqrt{3}}{4} \\ \frac{1}{2} \\ \frac{1}{4} \end{bmatrix}
\quad T(𝐞_3) = \begin{bmatrix} \frac{\sqrt{3}}{8} \\ \frac{1}{4} \\ \frac{1}{8} \end{bmatrix}$

$A = \begin{bmatrix} \frac{3}{8} & \frac{\sqrt{3}}{4} & \frac{\sqrt{3}}{8} \\
\frac{\sqrt{3}}{4} & \frac{1}{2} & \frac{1}{4} \\ \frac{\sqrt{3}}{8} & \frac{1}{4} & \frac{1}{8} \end{bmatrix}$

$\lambda I - A = \begin{bmatrix} \lambda - \frac{3}{8} & \frac{\sqrt{3}}{4} & \frac{\sqrt{3}}{8} \\
\frac{\sqrt{3}}{4} & \lambda - \frac{1}{2} & \frac{1}{4} \\ \frac{\sqrt{3}}{8} & \frac{1}{4} & \lambda - \frac{1}{8} \end{bmatrix}$

$A - \lambda I = \begin{bmatrix} \frac{3}{8} - \lambda & \frac{\sqrt{3}}{4} & \frac{\sqrt{3}}{8} \\
\frac{\sqrt{3}}{4} & \frac{1}{2} - \lambda & \frac{1}{4} \\ \frac{\sqrt{3}}{8} & \frac{1}{4} & \frac{1}{8} - \lambda \end{bmatrix}$
"""

# ╔═╡ 21779de1-f5a0-4908-800b-47b285255cc6
md"""
$\begin{align*}
&\det(\lambda I - T) \\
&= \left(\lambda - \frac{3}{8}\right) \begin{vmatrix} \lambda - \frac{1}{2} & \frac{1}{4} \\ \frac{1}{4} & \lambda - \frac{1}{8} \end{vmatrix} - \frac{\sqrt{3}}{4} \begin{vmatrix} \frac{\sqrt{3}}{4} & \frac{\sqrt{3}}{8} \\ \frac{1}{4} & \lambda - \frac{1}{8} \end{vmatrix} + \frac{\sqrt{3}}{8} \begin{vmatrix} \frac{\sqrt{3}}{4} & \frac{\sqrt{3}}{8} \\ \lambda - \frac{1}{2} & \frac{1}{4} \end{vmatrix} \\
&= \left(\lambda - \frac{3}{8}\right) \left[\left(\lambda - \frac{1}{2}\right)\left(\lambda - \frac{1}{8}\right) - \left(\frac{1}{4}\right)\left(\frac{1}{4}\right)\right] \\
&\quad - \left(\frac{\sqrt{3}}{4}\right) \left[\left(\frac{\sqrt{3}}{4}\right)\left(\lambda - \frac{1}{8}\right) - \left(\frac{\sqrt{3}}{8}\right)\left(\frac{1}{4}\right)\right] \\
&\quad + \left(\frac{\sqrt{3}}{8}\right) \left[\left(\frac{\sqrt{3}}{4}\right)\left(\frac{1}{4}\right) - \left(\frac{\sqrt{3}}{8}\right)\left(\lambda - \frac{1}{2}\right)\right] \\
&= \left(\lambda - \frac{3}{8}\right) \left(\lambda^2 - \frac{5}{8} \lambda\right) - \frac{\sqrt{3}}{4} \left(\frac{\sqrt{3}}{4} \lambda - \frac{\sqrt{3}}{16}\right) + \frac{\sqrt{3}}{8} \left(\frac{\sqrt{3}}{8} - \frac{\sqrt{3}}{8} \lambda\right) \\
&= \lambda^3 - \lambda^2 + \frac{15}{64} \lambda - \frac{3}{16} \lambda + \frac{3}{64} + \frac{3}{64} - \frac{3}{64} \lambda \\
&= \lambda^3 - \lambda^2 + \frac{6}{64}
\end{align*}$
"""

# ╔═╡ 48542abe-5050-41a3-8a31-3f7ad5fef612
md"""
$\begin{align*}
&\det(T - \lambda I) \\
&= \left(\frac{3}{8} - \lambda\right) \begin{vmatrix} \frac{1}{2} - \lambda & \frac{1}{4} \\ \frac{1}{4} & \frac{1}{8} - \lambda \end{vmatrix} - \frac{\sqrt{3}}{4} \begin{vmatrix} \frac{\sqrt{3}}{4} & \frac{\sqrt{3}}{8} \\ \frac{1}{4} & \frac{1}{8} - \lambda \end{vmatrix} + \frac{\sqrt{3}}{8} \begin{vmatrix} \frac{\sqrt{3}}{4} & \frac{\sqrt{3}}{8} \\ \frac{1}{2} - \lambda & \frac{1}{4} \end{vmatrix} \\
&= \left(\frac{3}{8} - \lambda\right) \left[\left(\frac{1}{2} - \lambda\right)\left(\frac{1}{8} - \lambda\right) - \left(\frac{1}{4}\right)\left(\frac{1}{4}\right)\right] \\
&\quad - \left(\frac{\sqrt{3}}{4}\right) \left[\left(\frac{\sqrt{3}}{4}\right)\left(\frac{1}{8} - \lambda\right) - \left(\frac{\sqrt{3}}{8}\right)\left(\frac{1}{4}\right)\right] \\
&\quad + \left(\frac{\sqrt{3}}{8}\right) \left[\left(\frac{\sqrt{3}}{4}\right)\left(\frac{1}{4}\right) - \left(\frac{\sqrt{3}}{8}\right)\left(\frac{1}{2} - \lambda\right)\right] \\
&= \left(\frac{3}{8} - \lambda\right) \left(\lambda^2 - \frac{5}{8} \lambda\right) - \frac{\sqrt{3}}{4} \left(-\frac{\sqrt{3}}{4} \lambda\right) + \frac{\sqrt{3}}{8} \left(\frac{\sqrt{3}}{8} \lambda\right) \\
&= \frac{3}{8} \lambda^2 - \frac{15}{64} \lambda - \lambda^3 + \frac{5}{8} \lambda^2 + \frac{3}{16} \lambda + \frac{3}{64} \lambda \\
&= -\lambda^3 + \lambda^2 \\
&= -\lambda^2 (\lambda - 1)
\end{align*}$

$\implies \lambda_1 = 0, \lambda_2 = 1$
"""

# ╔═╡ 00237bfe-45af-4d19-a61c-2da6ca08ec03
md"""
$2\sqrt{3} + 3x = 0 \implies x = -\frac{2\sqrt{3}}{3}$
$\sqrt{3} + 3x = 0 \implies x = -\frac{\sqrt{3}}{3}$

$\begin{align*}
(T - \lambda_1 I)𝐯 = T𝐯 = 𝟎 &\implies \begin{bmatrix} \frac{3}{8} & \frac{\sqrt{3}}{4} & \frac{\sqrt{3}}{8} & 0 \\ \frac{\sqrt{3}}{4} & \frac{1}{2} & \frac{1}{4} & 0 \\ \frac{\sqrt{3}}{8} & \frac{1}{4} & \frac{1}{8} & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 3 & 2\sqrt{3} & \sqrt{3} & 0 \\ 2\sqrt{3} & 4 & 2 & 0 \\ \sqrt{3} & 2 & 1 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 3 & 2\sqrt{3} & \sqrt{3} & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & \frac{2\sqrt{3}}{3} & \frac{\sqrt{3}}{3} & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{bmatrix} \\
\end{align*}$

$\text{Basis of } E_0 = \left\{\begin{bmatrix} -\frac{2\sqrt{3}}{3} \\ 1 \\ 0 \end{bmatrix}, \begin{bmatrix} -\frac{\sqrt{3}}{3} \\ 0 \\ 1 \end{bmatrix}\right\}$
"""

# ╔═╡ 97ead013-1341-4db6-9bb3-9dad84921512
md"""
$2\sqrt{3} - 5x = 0 \implies x = \frac{2\sqrt{3}}{5}$
$\sqrt{3} - 5x = 0 \implies x = \frac{\sqrt{3}}{5}$

$\begin{align*}
(T - \lambda_2 I)𝐯 = (T - I)𝐯 = 𝟎 &\implies \begin{bmatrix} -\frac{5}{8} & \frac{\sqrt{3}}{4} & \frac{\sqrt{3}}{8} & 0 \\ \frac{\sqrt{3}}{4} & -\frac{1}{2} & \frac{1}{4} & 0 \\ \frac{\sqrt{3}}{8} & \frac{1}{4} & -\frac{7}{8} & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} -5 & 2\sqrt{3} & \sqrt{3} & 0 \\ 2\sqrt{3} & -4 & 2 & 0 \\ \sqrt{3} & 2 & -7 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} -5 & 2\sqrt{3} & \sqrt{3} & 0 \\ 0 & -\frac{8}{5} & \frac{16}{5} & 0 \\ 0 & \frac{16}{5} & -\frac{32}{5} & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & -\frac{2\sqrt{3}}{5} & -\frac{\sqrt{3}}{5} & 0 \\ 0 & 1 & -2 & 0 \\ 0 & 1 & -2 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & 0 & -\sqrt{3} & 0 \\ 0 & 1 & -2 & 0 \\ 0 & 0 & 0 & 0 \end{bmatrix} \\
\end{align*}$

$\text{Basis of } E_1 = \left\{\begin{bmatrix} \sqrt{3} \\ 2 \\ 1 \end{bmatrix}\right\}$
"""

# ╔═╡ 37becf27-e980-4d6e-9d70-e96fe93adebd
let
	T = 1/8 * [3        2sqrt(3) sqrt(3)
			   2sqrt(3) 4        2
		 	   sqrt(3)  2        1]
	rref(T - I)
end

# ╔═╡ a6138b07-9ed8-4210-85f9-18d71026c8e9
md"""
### HW4 #2

Let $A$ be an $n × n$ matrix, and let

$q(\lambda) = q_3 \lambda^3 + q_2 \lambda^2 + q_1 \lambda + q_0 \lambda^0$

be some polynomial of degree 3 (the degree is not important).
Write down $q(A)$ by replacing every instance of $\lambda$ in $q$ by $A$ (use that $A^0$ is the identity matrix: $A^0 = I$).
Let $\lambda_0$ be an eigenvalue of $A$.
Why is $q(\lambda_0)$ an eigenvalue of $q(A)$?
"""

# ╔═╡ c8138ab3-dada-4b16-b495-834af5546e6a
md"""
$q(A) = q_3 A^3 + q_2 A^2 + q_1 A + q_0$

Let $\vec{v}$ be an eigenvector of $A$.
To show that $q(\lambda_0)$ is an eigenvalue of $q(A)$ then we need to show

$q(A) \vec{v} = q(\lambda_0) \vec{v}.$

Using

$A\vec{v} = \lambda_0 \vec{v} \implies A^2 \vec{v} = {\lambda_0}^2 \vec{v} \implies A^3 \vec{v} = {\lambda_0}^3 \vec{v}$

$\begin{align*}
q(A) \vec{v} &= (q_3 A^3 + q_2 A^2 + q_1 A + q_0) \vec{v} \\
&= q_3 A^3 \vec{v} + q_2 A^2 \vec{v} + q_1 A \vec{v} + q_0 \vec{v} \\
&= q_3 {\lambda_0}^3 \vec{v} + q_2 {\lambda_0}^2 \vec{v} + q_1 \lambda_0 \vec{v} + q_0 \vec{v} \\
&= (q_3 {\lambda_0}^3 + q_2 {\lambda_0}^2 + q_1 \lambda_0 + q_0) \vec{v} \\
&= q(\lambda_0) \vec{v}
\end{align*}$

Hence $q(\lambda_0)$ is an eigenvalue of $q(A)$ since $q(A) \vec{v} = q(\lambda_0) \vec{v}$.
"""

# ╔═╡ 74ba90f1-b990-40f8-bf44-adf7f127d449
md"## Homework #5"

# ╔═╡ 9d7c0449-5b2b-4e77-a2d9-8c5a8734307b
md"""
### HW5 #1

The characteristic polynomial of

$A = \begin{bmatrix} 3 & 0 & -1 & 0 & 0 \\ 1 & 3 & 0 & 0 & 1 \\ 0 & -1 & 1 & 0 & -1 \\ 3 & 2 & 0 & 2 & 3 \\ -1 & 0 & 1 & 0 & 2 \end{bmatrix}$

happens to be $p(\lambda) = (\lambda - 2)^4 (\lambda - 3)$
"""

# ╔═╡ 3fca3dec-6e40-4c4d-8782-6738d82ee5d2
let
	@variables x
	A = [3 0 -1 0 0; 1 3 0 0 1; 0 -1 1 0 -1; 3 2 0 2 3; -1 0 1 0 2]
	nullspace(A - 2I), nullspace((A - 2I)^2), nullspace((A - 2I)^4)
end

# ╔═╡ b25706e5-bfee-4a21-92d7-09b2c77484a4
md"""
### HW5 #2

Exactly one of the four matrices

$\begin{bmatrix} 2 & 0 & 0 & 0 & 0 \\ 0 & 2 & 0 & 0 & 0 \\ 0 & 0 & 2 & 0 & 0 \\ 0 & 0 & 0 & 2 & 0 \\ 0 & 0 & 0 & 0 & 3 \end{bmatrix}, \begin{bmatrix} 2 & 1 & 0 & 0 & 0 \\ 0 & 2 & 0 & 0 & 0 \\ 0 & 0 & 2 & 0 & 0 \\ 0 & 0 & 0 & 2 & 0 \\ 0 & 0 & 0 & 0 & 3 \end{bmatrix}, \begin{bmatrix} 2 & 1 & 0 & 0 & 0 \\ 0 & 2 & 1 & 0 & 0 \\ 0 & 0 & 2 & 0 & 0 \\ 0 & 0 & 0 & 2 & 0 \\ 0 & 0 & 0 & 0 & 3 \end{bmatrix}, \begin{bmatrix} 2 & 1 & 0 & 0 & 0 \\ 0 & 2 & 0 & 0 & 0 \\ 0 & 0 & 2 & 1 & 0 \\ 0 & 0 & 0 & 2 & 0 \\ 0 & 0 & 0 & 0 & 3 \end{bmatrix}$

has the same general behavior as the matrix of Problem 1:
All have the same characteristic polynomial as ``A``, but only one of them, call it ``A'``, satisfies

$\dim(\text{Null}(A' - 2I)) = \dim(\text{Null}(A - 2I)) \text{ and Null}((A' - 2I)^2) = \text{Null}((A - 2I)^4).$

Find it.
"""

# ╔═╡ ddfec193-2c2e-42d8-ba0f-af47a02b3ddc
let
	A = [3 0 -1 0 0; 1 3 0 0 1; 0 -1 1 0 -1; 3 2 0 2 3; -1 0 1 0 2]
	
	nullspace.([(A - 2I)^2, (A - 2I)^4])
end

# ╔═╡ cb942a73-1688-46b9-be93-836ff286f912
let
	A1 = [2 0 0 0 0; 0 2 0 0 0; 0 0 2 0 0; 0 0 0 2 0; 0 0 0 0 3]
	A2 = [2 1 0 0 0; 0 2 0 0 0; 0 0 2 0 0; 0 0 0 2 0; 0 0 0 0 3]
	A3 = [2 1 0 0 0; 0 2 1 0 0; 0 0 2 0 0; 0 0 0 2 0; 0 0 0 0 3]
	A4 = [2 1 0 0 0; 0 2 0 0 0; 0 0 2 1 0; 0 0 0 2 0; 0 0 0 0 3]

	rref.([A1 - 2I, A2 - 2I, A3 - 2I, A4 - 2I])
end

# ╔═╡ a797996c-852d-472e-a16b-927480be9cc3
let
	A1 = [2 0 0 0 0; 0 2 0 0 0; 0 0 2 0 0; 0 0 0 2 0; 0 0 0 0 3]
	nullspace.([(A1 - 2I), (A1 - 2I)^2])
	(A1 - 2I)
end

# ╔═╡ 1cf7d2ce-ae44-48b0-ac2f-72801923ad21
let
	A2 = [2 1 0 0 0; 0 2 0 0 0; 0 0 2 0 0; 0 0 0 2 0; 0 0 0 0 3]
	nullspace.([(A2 - 2I), (A2 - 2I)^2])
end

# ╔═╡ f8cb5e64-d615-47b4-8f97-de988bd3c144
let
	A3 = [2 1 0 0 0; 0 2 1 0 0; 0 0 2 0 0; 0 0 0 2 0; 0 0 0 0 3]
	# nullspace.([(A3 - 2I), (A3 - 2I)^2])
	(A3 - 2I)^2
end

# ╔═╡ a916da5f-d5fe-4001-aed0-64bedaf02a68
let
	A4 = [2 1 0 0 0; 0 2 0 0 0; 0 0 2 1 0; 0 0 0 2 0; 0 0 0 0 3]
	nullspace.([(A4 - 2I), (A4 - 2I)^2])
end

# ╔═╡ ecf7f160-daa9-47e7-9b25-17ac36415438
md"""
### HW5 #3

Let ``v_1, v_2, v_3`` be a basis of some three-dimensional inner product vector space.
Suppose the inner product ``⟨,⟩`` satisfies

$\begin{bmatrix} ⟨v_1,v_1⟩ & ⟨v_1,v_2⟩ & ⟨v_1,v_3⟩ \\ ⟨v_2,v_1⟩ & ⟨v_2,v_2⟩ & ⟨v_2,v_3⟩ \\ ⟨v_3,v_1⟩ & ⟨v_3,v_2⟩ & ⟨v_3,v_3⟩ \end{bmatrix} = \begin{bmatrix} 4 & 2 & 1 \\ 2 & 5 & 2 \\ 1 & 2 & 4 \end{bmatrix}.$

Apply Gram-Schmidt to the basis ``v_1,v_2,v_3,`` call the resulting vectors ``u_1,u_2,u_3``.
The ``u_i`` are certain linear combinations of the ``v_i``:

$u_1 = p_{11} v_1 + p_{21} v_2 + p_{31} v_3, \; u_2 = p_{12} v_1 + p_{22} v_2 + p_{32} v_3, \; u_3 = p_{13} v_1 + p_{23} v_2 + p_{33} v_3.$

Form the matrix ``P = [p_{ij}]`` (``i`` indexes the rows), find ``\det{P}``.
"""

# ╔═╡ 51fdc40e-2e75-40e2-8bf2-a46055fe2c4d
md"""
$\begin{align*}
u_1 &= v_1 \\
u_2 &= v_2 - \frac{⟨u_1, v_2⟩}{⟨u_1,u_1⟩} u_1 \\
&= v_2 - \frac{1}{2} v_1 \\
&= -\frac{1}{2} v_1 + v_2 \\
u_3 &= v_3 - \frac{⟨u_1, v_3⟩}{⟨u_1,u_1⟩} u_1 - \frac{⟨u_2, v_3⟩}{⟨u_2,u_2⟩} u_2 
\\
&= v_3 - \frac{1}{4} v_1 - \frac{⟨v_2 - \frac{1}{2} v_1, v_3⟩}{⟨v_2 - \frac{1}{2} v_1,v_2 - \frac{1}{2} v_1⟩} (v_2 - \frac{1}{2} v_1) \\
&= v_3 - \frac{1}{4} v_1 - \frac{⟨v_2, v_3⟩ - \frac{1}{2} ⟨v_1, v_3⟩}{⟨v_2, v_2 - \frac{1}{2} v_1⟩ - \frac{1}{2} ⟨v_1,v_2 - \frac{1}{2} v_1⟩} (v_2 - \frac{1}{2} v_1) \\
&= v_3 - \frac{1}{4} v_1 - \frac{⟨v_2, v_3⟩ - \frac{1}{2} ⟨v_1, v_3⟩}{(⟨v_2, v_2⟩ - \frac{1}{2} ⟨v_2, v_1⟩) - \frac{1}{2} (⟨v_1,v_2⟩ - \frac{1}{2} ⟨v_1,v_1⟩)} (v_2 - \frac{1}{2} v_1) \\
&= v_3 - \frac{1}{4} v_1 - \frac{2 - \frac{1}{2}}{(5 - 1) - \frac{1}{2} (2 - 2)} (v_2 - \frac{1}{2} v_1) \\
&= v_3 - \frac{1}{4} v_1 - \frac{3}{8} (v_2 - \frac{1}{2} v_1) \\
&= v_3 - \frac{1}{4} v_1 - \frac{3}{8} v_2 - \frac{3}{16} v_1 \\
&= -\frac{7}{16} v_1 - \frac{3}{8} v_2 + v_3
\end{align*}$
"""

# ╔═╡ 1f609164-a435-4ca7-b41c-9064397ffaf0
md"""
$P = \begin{bmatrix} 1 & -\frac{1}{2} & -\frac{7}{16} \\ 0 & 1 & -\frac{3}{8} \\ 0 & 0 & 1 \end{bmatrix}$

$\det(P) = 1$
"""

# ╔═╡ d9f4fade-7503-468e-9b8f-60e51b2bad02
md"""
### HW5 #4

Let ``V``, the basis, and the inner product be as in Problem 3.
The subspace of ``V`` orthogonal to ``v_3`` consists of all vectors ``v = a_1 v_1 + a_2 v_2 + a_3 v_3`` such that

$⟨a_1 v_1 + a_2 v_2 + a_3 v_3, v_3⟩ = 0.$

This is a condition (a linear equation of rank 1) on the coefficients:
one equation in three unknowns.
Find two independent solutions and use them to write a basis of ``\{v_3\}^\perp``
"""

# ╔═╡ 19913efa-5aaa-402f-b85d-058296e592db
md"""
$\begin{align*}
⟨a_1 v_1 + a_2 v_2 + a_3 v_3, v_3⟩ &= 0 \\
a_1 ⟨v_1, v_3⟩ + a_2 ⟨v_2, v_3⟩ + a_3 ⟨v_3, v_3⟩ &= 0 \\
a_1 + 2a_2 + 4a_3 &= 0
\end{align*}$

$\begin{bmatrix} 1 & 2 & 4 \end{bmatrix}$

$x = \begin{bmatrix} -2a_2 - 4a_3 \\ a_2 \\ a_3 \end{bmatrix} = \begin{bmatrix} -2a_2 \\ a_2 \\ 0 \end{bmatrix} + \begin{bmatrix} -4a_3 \\ 0 \\ a_3 \end{bmatrix} = a_2 \begin{bmatrix} -2 \\ 1 \\ 0 \end{bmatrix} + a_3 \begin{bmatrix} -4 \\ 0 \\ 1 \end{bmatrix}$

$\left\{\begin{bmatrix} -2 \\ 1 \\ 0 \end{bmatrix}, \begin{bmatrix} -4 \\ 0 \\ 1 \end{bmatrix}\right\}$
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Polynomials = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
RowEchelon = "af85af4c-bcd5-5d23-b03a-a909639aa875"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
Polynomials = "~2.0.17"
RowEchelon = "~0.2.1"
Symbolics = "~4.1.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgCheck]]
git-tree-sha1 = "dedbbb2ddb876f899585c4ec4433265e3017215a"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.1.0"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "265b06e2b1f6a216e0e8f183d28e4d354eab3220"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.2.1"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[AutoHashEquals]]
git-tree-sha1 = "45bb6705d93be619b81451bb2006b7ee5d4e4453"
uuid = "15f4f7f2-30c1-5605-9d31-71845cf9641f"
version = "0.2.0"

[[BangBang]]
deps = ["Compat", "ConstructionBase", "Future", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables", "ZygoteRules"]
git-tree-sha1 = "0ad226aa72d8671f20d0316e03028f0ba1624307"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.32"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[Bijections]]
git-tree-sha1 = "705e7822597b432ebe152baa844b49f8026df090"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.3"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4c26b4e9e91ca528ea212927326ece5918a04b47"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.2"

[[ChangesOfVariables]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "9a1d594397670492219635b35a3d830b04730d62"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.1"

[[Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dce3e3fea680869eaa0b774b2e8343e9ff442313"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.40.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[CompositeTypes]]
git-tree-sha1 = "d5b014b216dc891e81fea299638e4c10c657b582"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.2"

[[CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DefineSingletons]]
git-tree-sha1 = "77b4ca280084423b728662fe040e5ff8819347c5"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.1"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[DiffRules]]
deps = ["LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "d8f468c5cd4d94e86816603f7d18ece910b4aaf1"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.5.0"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "7f3bec11f4bcd01bc1f507ebce5eadf1b0a78f47"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.34"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays", "Statistics"]
git-tree-sha1 = "5f5f0b750ac576bcf2ab1d7782959894b304923e"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.9"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "1b4665a7e303eaa7e03542cfaef0730cb056cb00"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.3.21"

[[EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "3fe985505b4b667e1ae303c9ca64d181f09d5c05"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.1.3"

[[ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "8756f9935b7ccc9064c6eef0bff0ad643df733a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.7"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[InitialValues]]
git-tree-sha1 = "7f6a4508b4a6f46db5ccd9799a3fc71ef5cad6e6"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.2.11"

[[InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "19cb49649f8c41de7fea32d089d37de917b553da"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.0.1"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "3cc368af3f110a767ac786560045dceddfc16758"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.3"

[[Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "323a38ed1952d30586d0fe03412cde9399d3618b"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.5.0"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "LinearAlgebra", "MacroTools", "StaticArrays"]
git-tree-sha1 = "3609bbf5feba7b22fb35fe7cb207c8c8d2e2fc5b"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.6.7"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "be9eef9f9d78cecb6f262f3c10da151a6c5ab827"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.5"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Metatheory]]
deps = ["AutoHashEquals", "DataStructures", "Dates", "DocStringExtensions", "Parameters", "Reexport", "TermInterface", "ThreadsX", "TimerOutputs"]
git-tree-sha1 = "0d3b2feb3168e4deb78361d3b5bb5c2e51ea5271"
uuid = "e9d8d322-4543-424a-9be4-0cc815abe26c"
version = "1.3.2"

[[MicroCollections]]
deps = ["BangBang", "Setfield"]
git-tree-sha1 = "4f65bdbbe93475f6ff9ea6969b21532f88d359be"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "29714d0a7a8083bba8427a4fbfb00a540c681ce7"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.3"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "45c9940cec79dedcdccc73cc6dd09ea8b8ab142c"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.3.18"

[[MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "8d9496b2339095901106961f44718920732616bb"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "0.2.22"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "ee26b350276c51697c9c2d88a072b339f9f03d73"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.5"

[[Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "d911b6a12ba974dabe2291c6d450094a7226b372"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.1"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[Polynomials]]
deps = ["Intervals", "LinearAlgebra", "MutableArithmetics", "RecipesBase"]
git-tree-sha1 = "7499556d31417baeabaa55d266a449ffe4ec5a3e"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "2.0.17"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[RecursiveArrayTools]]
deps = ["ArrayInterface", "ChainRulesCore", "DocStringExtensions", "FillArrays", "LinearAlgebra", "RecipesBase", "Requires", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "c944fa4adbb47be43376359811c0a14757bdc8a8"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.20.0"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Referenceables]]
deps = ["Adapt"]
git-tree-sha1 = "e681d3bfa49cd46c3c161505caddf20f0e62aaa9"
uuid = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"
version = "0.1.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[RowEchelon]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f479526c4f6efcbf01e7a8f4223d62cfe801c974"
uuid = "af85af4c-bcd5-5d23-b03a-a909639aa875"
version = "0.2.1"

[[RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "cdc1e4278e91a6ad530770ebb327f9ed83cf10c4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.3"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "RecipesBase", "RecursiveArrayTools", "StaticArrays", "Statistics", "Tables", "TreeViews"]
git-tree-sha1 = "b3d23aa4e5f621b574b3b0d41c62c8624d27192a"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.19.5"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "0afd9e6c623e379f593da01f20590bacc26d1d14"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.1"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

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

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f0bccf98e16759818ffc5d97ac3ebf87eb950150"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.8.1"

[[SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "39c9f91521de844bad65049efd4f9223e7ed43f9"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.14"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "e7bc80dc93f50857a5d1e3c8121495852f407e6a"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.4.0"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "0f2aa8e32d511f758a2ce49208181f7733a0936a"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.1.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "2bb0cb32026a66037360606510fca5984ccc6b75"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.13"

[[StatsFuns]]
deps = ["ChainRulesCore", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "bedb3e17cc1d94ce0e6e66d3afa47157978ba404"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.14"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "Metatheory", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TermInterface", "TimerOutputs"]
git-tree-sha1 = "5255e65d129c8edbde92fd2ede515e61098f93df"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "0.18.1"

[[Symbolics]]
deps = ["ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "IfElse", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Metatheory", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TermInterface", "TreeViews"]
git-tree-sha1 = "1f738ebade1567d461c4db9ef15470a1fdf0b9b5"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "4.1.1"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "fed34d0e71b91734bf0a7e10eb1bb05296ddbcd0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[TermInterface]]
git-tree-sha1 = "7aa601f12708243987b88d1b453541a75e3d8c7a"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "0.2.3"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[ThreadsX]]
deps = ["ArgCheck", "BangBang", "ConstructionBase", "InitialValues", "MicroCollections", "Referenceables", "Setfield", "SplittablesBase", "Transducers"]
git-tree-sha1 = "abcff3ac31c7894550566be533b512f8b059104f"
uuid = "ac1d9e8a-700a-412c-b207-f0111f4b6c0d"
version = "0.1.8"

[[TimeZones]]
deps = ["Dates", "Downloads", "InlineStrings", "LazyArtifacts", "Mocking", "Pkg", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "8de32288505b7db196f36d27d7236464ef50dba1"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.6.2"

[[TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "7cb456f358e8f9d102a8b25e8dfedf58fa5689bc"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.13"

[[Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "bccb153150744d476a6a8d4facf5299325d5a442"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.67"

[[TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─fba59878-c8c1-4f42-8214-6d99568d8de1
# ╟─89066d76-248e-11ec-3893-757003b33e14
# ╟─6fce260e-d916-4ac1-902f-8366d179c9df
# ╠═58d62bd7-a0e7-4746-a2a2-fb09d053899b
# ╟─0fbade02-25e8-417d-b4c7-4ce30449ac62
# ╟─c6fb02f8-d4e1-42e3-b7c0-3a59402f1dbb
# ╟─13ca6959-333b-42ac-b879-566be5ef2615
# ╠═f8eb8ef6-b8a6-4b4f-b8eb-3cc3125796e5
# ╟─d5523a14-f49e-4a66-a009-d84578179de7
# ╠═c86dce76-be3a-44d5-ad62-5ae915349403
# ╟─c87c49d4-1139-46cb-bf43-17d6fd09c993
# ╟─bb4f86cd-b207-4478-a7b9-e4ff7146fb8c
# ╠═dc73982a-8dab-4173-8c70-e180404e7aec
# ╟─370ad756-c609-460e-a9f3-2d5305e594c6
# ╟─fb839ada-ae51-41dd-bf95-45a1cc60d6f0
# ╠═774cfe46-d41b-406a-ab40-8da0707ffe09
# ╟─f3fddbe8-67c4-46e2-802c-5bf292f69f11
# ╟─34582ac8-c4e4-483e-ba3c-b7b83c44d199
# ╟─aba0a023-be9b-4f41-8728-6d162ef8a558
# ╟─cdd75ddd-39e2-41da-9036-e0ff29d72aca
# ╟─6cff14e0-548f-47eb-a34f-7f373b9efe63
# ╟─c42f0ad9-08bb-42ce-8b31-63ce79651a79
# ╟─a7b4e7de-1765-4e9e-bf26-bd5678334717
# ╟─24dae75e-90f1-4d64-968b-addac182a164
# ╟─e50a0ec8-4336-43fe-9ebe-6fb82fefbf89
# ╟─8c06c619-f214-42da-8fe2-19b2ec9de635
# ╟─348ea2aa-aef7-4464-92f0-9c7317552dda
# ╟─06f5fec0-595b-4d29-8530-0175155683c4
# ╟─9bdb9179-b5e4-4dbf-9277-97c5b670e81d
# ╟─0d502224-dba6-4d35-a0a5-11ce171bbc02
# ╠═30dd5e24-eec8-4b5f-bd13-50fe1c3162c9
# ╟─fb879466-d6b8-4200-86fc-ce95b58d7794
# ╟─21779de1-f5a0-4908-800b-47b285255cc6
# ╟─48542abe-5050-41a3-8a31-3f7ad5fef612
# ╟─00237bfe-45af-4d19-a61c-2da6ca08ec03
# ╟─97ead013-1341-4db6-9bb3-9dad84921512
# ╠═37becf27-e980-4d6e-9d70-e96fe93adebd
# ╟─a6138b07-9ed8-4210-85f9-18d71026c8e9
# ╟─c8138ab3-dada-4b16-b495-834af5546e6a
# ╟─74ba90f1-b990-40f8-bf44-adf7f127d449
# ╟─9d7c0449-5b2b-4e77-a2d9-8c5a8734307b
# ╠═3fca3dec-6e40-4c4d-8782-6738d82ee5d2
# ╟─b25706e5-bfee-4a21-92d7-09b2c77484a4
# ╠═ddfec193-2c2e-42d8-ba0f-af47a02b3ddc
# ╠═cb942a73-1688-46b9-be93-836ff286f912
# ╠═a797996c-852d-472e-a16b-927480be9cc3
# ╠═1cf7d2ce-ae44-48b0-ac2f-72801923ad21
# ╠═f8cb5e64-d615-47b4-8f97-de988bd3c144
# ╠═a916da5f-d5fe-4001-aed0-64bedaf02a68
# ╟─ecf7f160-daa9-47e7-9b25-17ac36415438
# ╟─51fdc40e-2e75-40e2-8bf2-a46055fe2c4d
# ╟─1f609164-a435-4ca7-b41c-9064397ffaf0
# ╟─d9f4fade-7503-468e-9b8f-60e51b2bad02
# ╟─19913efa-5aaa-402f-b85d-058296e592db
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
