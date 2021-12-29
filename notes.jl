### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ 7b8bf49e-fb76-4034-a400-22f90e1983a3
let
	using PlutoUI
	using Plots
	
	function Image(filepath::String)
		@assert !isempty(filepath)
		PlutoUI.LocalResource(joinpath(split(@__FILE__, '#')[1] * ".assets", filepath))
	end

	
	md"# Linear Algebra"
end

# ╔═╡ 542aec7c-3740-4bba-a572-87971b5614fe
PlutoUI.TableOfContents()

# ╔═╡ a5c5f744-82df-4db9-a958-1a9f50005072
md"# Systems of linear equations"

# ╔═╡ 73440efb-2a36-414f-9f53-23d07b792574
begin
	p1 = plot([x -> x, x -> -x + 2], title = "One solution")
	scatter!(p1, (1, 1), color = :red)
	p2 = plot([x -> x, x -> x + 2], title = "No solutions")
	p3 = plot([x -> x, x -> x], title = "Infinitely many solutions")
	plot(p1, p2, p3, legend = false, ticks = nothing)
end

# ╔═╡ 830a5a05-7151-4a70-ae28-6353d50290bb
md"## Geometric view of systems of equations"

# ╔═╡ f13a313a-0359-4ae7-b0aa-c0616df64234
md"## Algebraic view of systems of equations"

# ╔═╡ 459a45a1-1013-4b5a-afbf-a7c9e9556887
md"## Elementary operations"

# ╔═╡ 6e46f7d2-ae75-48c8-b4ce-131b52c0f795
md"## Gaussian elimination"

# ╔═╡ 46d7cf27-29b2-48a4-91e5-808bbae92c68
md"## Gauss-Jordan elimination"

# ╔═╡ d15db133-0b81-401c-9935-5d81bfc6ba9f
md"## Homogeneous systems"

# ╔═╡ 724b61b0-117a-4132-9431-659dd8e8d5ad
md"## Uniqueness of the reduced echelon form"

# ╔═╡ d5bdb6e7-185d-4309-a00d-ce3fa544c4c5
md"## Fields"

# ╔═╡ aa81baa3-fde8-441a-8bdc-7fcdd1d7ff96
md"## Application: Balancing chemical reactions"

# ╔═╡ 90873450-e0c4-4e5a-ad1c-a260e6caf825
md"## Application: Dimensionless variables"

# ╔═╡ 7c17981e-ee11-4e39-bf7d-7cdf94dbb5cf
md"## Application: Resistor networks"

# ╔═╡ a89b4ad7-31fc-418c-bb95-36187aadff22
md"# Vectors in ℝⁿ"

# ╔═╡ 091bb6cb-e600-4df4-9c5a-58c6de4e48cd
md"## Points and vectors"

# ╔═╡ be63bec1-3317-457d-b17e-4f3b232dbd08
md"## Addition"

# ╔═╡ 5e0eac80-5837-4a6f-af5f-f587270751fc
md"## Scalar multiplication"

# ╔═╡ 53f07204-8912-40cb-be72-6678635c731d
md"## Linear combinations"

# ╔═╡ 5b3752a9-9391-4c8a-8daa-32a94c75f9ae
md"## Length of a vector"

# ╔═╡ f8a521e9-7d13-419f-8fd9-b8f98e7614f2
md"## The dot product"

# ╔═╡ 2017c2ed-ebe5-4eb7-8825-44751abe6076
md"""
### Outcomes

A. Compute the dot product of vectors geometrically and algebraically.

B. Use properties of the dot product, including the Cauchy-Schwarz inequality and the triangle inequality, to prove further equalities and inequalities.

C. Determine whether two vectors are orthogonal.

D. Compute the scalar and vector projection of one vector onto another.

E. Decompose a vector into orthogonal components.
"""

# ╔═╡ 41e6ca5e-a95e-42e6-8ab4-b7721028eb2f
md"""
### Proposition 2.27: Cauchy-Schwarz inequality

The dot product satisfies the inequality

$|𝐮 ⋅ 𝐯| ≤ \|𝐮\|\|𝐯\|$

Furthermore equality is obtained if and only if one of $𝐮$ or $𝐯$ is a scalar multiple of the other.
"""

# ╔═╡ 3d20ee51-69f8-4ba1-b005-383f8410b5d5
md"""
### Proposition 2.28: Triangle inequality

For $𝐮, 𝐯 ∈ ℝ^n$, we have

$\|𝐮 + 𝐯\| ≤ \|𝐮\| + \|𝐯\|.$
"""

# ╔═╡ 91baaf1c-479a-46a4-9f75-cca70d08c80a
md"""
### Definition 2.35: Vector projection

Let $𝐮$ be a non-zero vector and $𝐯$ any vector.
Then the **component of $𝐯$ in the direction of $𝐮$** is defined to be the scalar

$\text{comp}_𝐮(𝐯) = \frac{𝐮 ⋅ 𝐯}{\|𝐮\|}.$

The **projection of $𝐯$ onto $𝐮$** is defined to be the vector

$\text{proj}_𝐮(𝐯) = \frac{𝐮 ⋅ 𝐯}{𝐮 ⋅ 𝐮} 𝐮 = \frac{𝐮 ⋅ 𝐯}{\|𝐮\|^2} 𝐮.$

These two operations are called the **scalar projection** and **vector projection**, respectively.
"""

# ╔═╡ b6d18658-7b6b-4f6a-b9b0-ece83963a470
md"## The cross product"

# ╔═╡ 149b568a-70a9-40bd-9628-e3d1485f3216
md"# Lines and planes in ℝⁿ"

# ╔═╡ afe0d69e-3d2a-4fab-9daa-6502918cfc85
md"## Lines"

# ╔═╡ a2dfbd75-53b1-4d1e-9abd-52ca3f165e5d
md"## Planes"

# ╔═╡ b1158480-08c9-48cc-8f01-30a84437b9e2
md"# Matrices"

# ╔═╡ ffb118ad-8baf-4556-9062-0480de6459f3
md"## Definition and equality"

# ╔═╡ 94d46bac-23fc-413e-aaf8-ee2940233eb6
md"## Addition"

# ╔═╡ 9c28afd7-7203-4d9f-91cb-ca7ae1cfa823
md"## Scalar multiplication"

# ╔═╡ 203cc5d2-45be-4704-9bfe-9410e2efd3bf
md"## Matrix multiplication"

# ╔═╡ da3046e7-7436-4e0e-9fe3-bbf28b15c267
md"## Matrix inverses"

# ╔═╡ 3b3bc7a5-6455-41e5-a872-8b7a252e37c6
md"## Elementary matrices"

# ╔═╡ bd75f8f7-c31e-450f-b544-be72ddc60b9a
md"## The transpose"

# ╔═╡ 4d65dbd5-e6c8-482d-a119-53c8f59676c8
md"## Matrix arithmetic modulo 𝑝"

# ╔═╡ db7ec3f6-6e83-4155-83fc-894d1ba9b521
md"## Application: Cryptography"

# ╔═╡ 8bded7bf-2fa5-42d3-ba49-650e78df4d8a
md"# Spans, linear independence, and bases in ℝⁿ"

# ╔═╡ 63fe60cb-df4c-4ba8-9f9c-fb52604a6bd7
md"## Spans"

# ╔═╡ 45ae2f7b-d464-45e3-b416-531dbb78627d
md"## Linear independence"

# ╔═╡ 09abcb81-943c-42cb-b3e0-492353415fe8
md"## Subspaces of ℝⁿ"

# ╔═╡ 25d798ad-879c-4bb4-a538-c26d28f3f486
md"## Basis and dimension"

# ╔═╡ 3686f34f-bb50-4e86-be62-62d4b169ca5b
md"## Column space, row space, and null space of a matrix"

# ╔═╡ 8b863a5c-f515-4669-94c4-d64dbe928bdf
md"# Linear transformations in ℝⁿ"

# ╔═╡ 1d3a1b05-345e-4189-ad4e-fb47c4871804
md"## Linear transformations"

# ╔═╡ 8afc742e-09ab-4ea3-985c-84a7508c27b9
md"""
### Outcomes

A. Determine whether a vector function $T : ℝ^n → ℝ^m$ is a linear transformation.
"""

# ╔═╡ d32a3b9d-d156-4bcb-9072-d7e34b6c600c
md"""
### Definition 6.1: Linear transformation

A vector function $T : ℝ^n → ℝ^m$ is called a **linear transformation**, or simply **linear**, if it satisfies the following two conditions:

1. ``T`` preserves addition, i.e., for all $𝐯, 𝐰 ∈ ℝ^n$, we have $T(𝐯 + 𝐰) = T(𝐯) + T(𝐰)$;

2. ``T`` preserves scalar multiplication, i.e., for all $𝐯 ∈ ℝ^n$ and scalars $k$, we have $T(k𝐯) = kT(𝐯)$.
"""

# ╔═╡ 48575c01-ef6d-4db2-9144-6ec4850044e4
md"""
### Proposition 6.3: Alternative characterization of linear transformations

A vector function $T : ℝ^n → ℝ^m$ is linear if and only if it satisfies the following condition, for all $𝐯, 𝐰 ∈ ℝ^n$ and scalars $a, b$:

$T(a𝐯 + b𝐰) = aT(𝐯) + bT(𝐰).$
"""

# ╔═╡ c42b1c38-6cca-4138-8f54-bddde6b3ad5f
md"## The matrix of a linear transformation"

# ╔═╡ fb7cc1b2-6573-4cf9-85b2-ebf6ab59a274
md"""
### Outcomes

A. Find the matrix corresponding to a linear transformation $T : ℝ^n → ℝ^m$.
"""

# ╔═╡ 508de533-96d5-4817-bbbc-3fc2c41fa848
md"""
### Proposition 6.4: Matrix transformations are linear transformations

Let $A$ be an $m × n$-matrix, and consider the vector function $T : ℝ^n → ℝ^m$ defined by $T(𝐯) = A𝐯$.
Then $T$ is a linear transformation.
"""

# ╔═╡ 033cca18-dcb6-4720-bf13-5379cd334df2
md"""
### Theorem 6.5: Linear transformations are matrix transformations

Let $T : ℝ^n → ℝ^m$ be any linear transformation.
Then there exists an $m × n$-matrix $A$ such that for all $𝐯 ∈ ℝ^n$,

$T(𝐯) = A𝐯.$

In other words, $T$ is a matrix transformation.
"""

# ╔═╡ bdb1858d-f3ab-44f9-96da-759c31998d6f
md"## Geometric interpretation of linear transformations"

# ╔═╡ 48ac0870-e5e3-4581-94ee-182c87b4e78e
md"## Properties of linear transformations"

# ╔═╡ 3cde46c8-3576-4bd4-b5e4-aa8040ab7037
md"""
### Outcomes

A. Use properties of linear transformations to solve problems.

B. Find the composite of transformations and the inverse of a transformation.
"""

# ╔═╡ dbc16007-9c67-4b5c-9c88-c3680294138a
md"""
### Proposition 6.14: Properties of linear transformations

Let $T : ℝ^n → ℝ^m$ be a linear transformation.
Then

- ``T`` preserves the zero vector $T(𝟎) = 𝟎$.

- ``T`` preserves negation: $T(-𝐯) = -T(𝐯)$.

- ``T`` preserves linear combinations:

  $T(a_1 𝐯_1 + a_2 𝐯_2 + … + a_k 𝐯_k) = a_1 T(𝐯_1) + a_2 T(𝐯_2) + … + a_k T(𝐯_k).$
"""

# ╔═╡ b980e504-1e37-4d17-a246-685fbb9a346e
md"""
### Definition 6.16: Composition of linear transformations

Let $S : ℝ^k → ℝ^n$ and $T : ℝ^n → ℝ^m$ be linear transformations.
Then the **composition** of $S$ and $T$ (also called the **composite transformation** of $S$ and $T$) is the linear transformation

$T ∘ S : ℝ^k → ℝ^m$

that is defined by

$(T ∘ S)(𝐯) = T(S(𝐯)).$

for all $𝐯 ∈ ℝ^k$.
"""

# ╔═╡ 907c43a9-d4e6-4043-bca9-891f4a81ac1d
md"""
### Theorem 6.17: Matrix of a composite transformation

Let $S : ℝ^k → ℝ^n$ and $T : ℝ^n → ℝ^m$ be linear transformations.
Let $A$ be the matrix corresponding to $S$, and let $B$ be the matrix corresponding to $T$.
Then the matrix corresponding to the composite linear transformation $T ∘ S$ is $BA$.
"""

# ╔═╡ 8e304446-19d7-4d39-bdb0-8057b3b4b057
md"""
### Definition 6.20: Inverse of a transformation

Let $T, S : ℝ^n → ℝ^n$ be linear transformations.
Suppose that for each $𝐯 ∈ ℝ^n$,

$(S ∘ T)(𝐯) = 𝐯$

and

$(T ∘ S)(𝐯) = 𝐯.$

Then $S$ is called the **inverse** of $T$, and we write $S = T^{-1}$.
"""

# ╔═╡ 60cf4ac0-7825-4385-99fb-0d1c1811a1e0
md"""
### Theorem 6.22: Matrix of the inverse transformation

Let $T : ℝ^n → ℝ^n$ be a linear transformation and let $A$ be the corresponding $n × n$-matrix.
Then $T$ has an inverse if and only if the matrix $A$ is invertible.
In this case, the matrix of $T^{-1}$ is $A^{-1}$.
"""

# ╔═╡ 960b204f-5721-4c76-a91a-e61b64c22224
md"## Application: Perspective rendering"

# ╔═╡ b28fe842-847f-45d8-aa8a-5a5a66467578
md"# Determinants"

# ╔═╡ 243bc5e3-3ce9-404a-920b-25ac294f92b9
md"## Determinants of 2 × 2- and 3 × 3-matrices"

# ╔═╡ 075f38de-966e-405a-ae01-b4a3b179a502
md"""
### Outcomes

A. Calculate the determinant of $2 × 2$-matrices and $3 × 3$-matrices.
"""

# ╔═╡ e32b08fb-def4-43a6-b031-f28ae5b56028
md"""
### Definition 7.1: Determinant of a $2 × 2$-matrix

Let $A = \begin{bmatrix} a & b \\ c & d \end{bmatrix}$. Then

$\det(A) = ad - bc.$
"""

# ╔═╡ 432b500e-368b-4a45-a888-5390efec5508
md"""
### Definition 7.3: Determinant of a $3 × 3$-matrix

Let $A = \begin{bmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{bmatrix}$.
Then

$\det(A) = a_{11} a_{22} a_{33} + a_{12} a_{23} a_{31} + a_{13} a_{21} a_{32} - a_{31} a_{22} a_{13} - a_{32} a_{23} a_{11} - a_{33} a_{21} a_{12}.$
"""

# ╔═╡ 0889fdf6-32f5-4566-a9d5-bcabed0bf24f
md"## Minors and cofactors"

# ╔═╡ 4260cefc-050f-4da1-8257-c4eee5be38fa
md"""
### Outcomes

A. Compute minors and cofactors of matrices.

B. Use cofactor expansion to compute the determinant of an $n × n$-matrix.
"""

# ╔═╡ 59c4f604-abf7-4589-802b-eb12f3bc0553
md"""
### Definition 7.5: The $ij^{\text{th}}$ minor of a matrix

Let $A$ be an $n × n$-matrix.
The $ij^{\text{th}}$ **minor** of $A$, denoted by $M_{ij}$, is the determinant of the $(n - 1) × (n - 1)$-matrix that is obtained by deleting the $i^{\text{th}}$ row and the $j^{\text{th}}$ column of $A$.
"""

# ╔═╡ 959dd506-e742-4f4c-a9b1-450032b9c3c7
md"""
### Definition 7.7: The $ij^{\text{th}}$ cofactor of a matrix

Suppose $A$ is an $n × n$-matrix.
The $ij^{\text{th}}$  **cofactor**, denoted by $C_{ij}$, is defined to be

$C_{ij} = (-1)^{i + j} M_{ij}$
"""

# ╔═╡ 024b87e3-f48f-46e5-bdca-84f383146d47
md"""
### Definition 7.9: The determinant of an $n × n$-matrix

Let $A$ be an $n × n$-matrix.
Then $\det(A)$ is calculated by picking a row (or column) and taking the product of each entry in that row (column) with its cofactor and adding these products together.
This process is known as **expanding along the $i^{\text{th}}$ row (or column)**.

In formulas, the process of expanding along the $i^{\text{th}}$ row is given as follow:

$\det(A) = a_{i1} C_{i1} + a_{i2} C_{i2} + … + a_{in} C_{in}.$

Similarly, the process of expanding along the $j^{\text{th}}$ column is:

$\det(A) = a_{1j} C_{1j} + a_{2j} C_{2j} + … + a_{nj} C_{nj}.$
"""

# ╔═╡ 56778fde-510b-4877-8c0d-7f4a2e39cea3
md"""
### Theorem 7.11: The determinant is well-defined

Expanding an $n × n$-matrix along any row or column always gives the same answer, which is the determinant.
"""

# ╔═╡ 27a52480-d576-4bdf-a063-7e005fae16c9
md"## The determinant of a triangular matrix"

# ╔═╡ f7de8cf8-ef6a-461c-87a6-046618ea9ec3
md"""
### Outcomes

A. Calculate the determinant of an upper or lower triangular matrix.
"""

# ╔═╡ 61cd090b-868f-4ca8-b3a3-198bfaf600a7
md"""
### Definition 7.13: Triangular matrices

A square matrix $A$ is **upper triangular** if $a_{ij} = 0$ whenever $i > j$.
In other words, a matrix is upper triangular if the entries below the main diagonal are 0.
Thus, an upper triangular matrix looks as follows, where $∗$ refers to any non-zero number:

$\begin{bmatrix}
∗ & ∗ & ⋯ & ∗ & ∗ \\
0 & ∗ & ⋯ & ∗ & ∗ \\
⋮ & ⋮ & ⋱ & ⋮ & ⋮ \\
0 & 0 & ⋯ & ∗ & ∗ \\
0 & 0 & ⋯ & 0 & ∗
\end{bmatrix}$

Similarly, a square matrix is **lower triangular** if entries above the main diagonal are 0.
"""

# ╔═╡ 4736289b-2a7f-415d-b2ac-ae7c15499197
md"""
### Theorem 7.14: Determinant of a triangular matrix

Let $A$ be an upper or lower triangular matrix.
Then $\det(A)$ is equal to the product of the entries on the main diagonal.
Written as a formula, we have

$\det(A) = a_{11} a_{22} ⋯ a_{nn}.$
"""

# ╔═╡ 51ecfe8d-c4a3-4ba1-a6f1-7f792f2e0933
md"## Determinants and row operations"

# ╔═╡ ca305f42-fb86-41a0-ba71-189de923346a
md"""
### Outcomes

A. Determine the effect of a row operation on the determinant of a matrix.

B. Use row operations to calculate a determinant.
"""

# ╔═╡ 4354eea7-5e35-4979-8c44-53a34c1ddffb
md"""
### Theorem 7.16: Effect of row operations on the determinant

Let $A$ be an $n × n$-matrix.

1. If $B$ is obtained from $A$ by switching two rows, then

$\det(B) = -\det(A).$

2. If $B$ is obtained from $A$ by multiplying one row by a non-zero scalar $k$, then

$\det(B) = k \det(A).$

3. If $B$ is obtained from $A$ by adding a multiple of one row to another row, then

$\det(B) = \det(A).$
"""

# ╔═╡ 928d92a1-baa6-4b36-824d-7c56f85c857c
md"## Properties of determinants"

# ╔═╡ b201693c-d02e-4e2c-86b4-3b2db1c5fa8c
md"""
### Outcomes

A. Use the determinant of a square matrix to decide whether the matrix is invertible.

B. From the determinants of two matrices, calculate the determinant of their product.

C. From the determinant of a matrix, calculate the determinant of its inverse.

D. From the determinant of a matrix, calculate the determinant of its transpose.

E. Calculate the determinant of $kA$, if the determinant of $A$ is known.

F. Without calculation, find the determinant of a matrix containing a row or column of Eros, or a matrix containing a row (or column) that is a scalar multiple of another row (or column).

G. Use algebraic properties to reason about determinants.
"""

# ╔═╡ aa765c1b-1f13-45ef-bbef-851781a64d29
md"""
### Theorem 7.20: Determinants and invertible matrices

Let $A$ be an $n × n$-matrix.
Then $A$ is invertible if and only if $\det(A) ≠ 0$.
"""

# ╔═╡ f819b979-07b8-47b2-866c-6d77889f3b6d
md"""
### Corollary 7.22: Determinants and homogenous systems

Let $A$ be an $n × n$-matrix.
Then the homogenous system $A𝐯 = 𝟎$ has non-trivial solutions if and only if $\det(A) = 0$.
"""

# ╔═╡ 93310b8c-09b7-4f3c-936f-d95586c5f736
md"""
### Theorem 7.23: Determinant of a product

Let $A$ and $B$ be $n × n$-matrices.
Then

$\det(AB) = \det(A) \det(B)$
"""

# ╔═╡ 6a7d7172-34d3-420a-8d36-fb51739e98f4
md"""
### Proposition 7.25: Properties of determinants

Let $A, B$ be $n × n$-matrices.
Then:

1. ``\det(AB) = \det(A) \det(B)``.

2. ``\det(I) = 1``.

3. ``A`` is invertible if and only if $\det(A) ≠ 0$.
   Moreover, if this is the case, then

   $\det(A^{-1}) = \frac{1}{\det(A)}.$

4. ``\det(kA) = k^n \det(A)``.

5. ``\det(A^T) = \det(A)``.
"""

# ╔═╡ a44deee1-5ecc-47b3-81b3-c30f0b53b2f5
md"""
### Theorem 7.26: Special matrices with zero determinant

Let $A$ be an $n × n$-matrix.

1. If $A$ has a row consisting only of zeros, or a column consisting only of zeros, then $\det(A) = 0$.

2. If $A$ has a row that is a scalar multiple of another row, or a column that is a scalar multiple of another column, then $\det(A) = 0$.
"""

# ╔═╡ 1f7f2aa2-5913-4b4d-8c2d-94d76487c690
md"## Application: A formula for the inverse of a matrix"

# ╔═╡ a573916f-280b-4055-9a30-d9aa663d7154
md"""
### Outcomes

A. Find the cofactor matrix and the adjugate of a matrix.

B. Find the inverse of a matrix using the adjugate formula.
"""

# ╔═╡ 6869f242-3f47-4c88-a669-b525bd101e8c
md"""
### Definition 7.27: The cofactor matrix

Let $A$ be an $n × n$-matrix.
Then the **cofactor matrix of $A$**, detnoed $\text{cof}(A)$, is defined by

$\text{cof}(A) = \begin{bmatrix} C_{ij} \end{bmatrix} = \begin{bmatrix} C_{11} & C_{12} & ⋯ & C_{1n} \\ C_{21} & C_{22} & ⋯ & C_{2n} \\ ⋮ & ⋮ & ⋱ & ⋮ \\ C_{n1} & C_{n2} & ⋯ & C_{nn} \end{bmatrix}$

where $C_{ij}$ is the $ij^{\text{th}}$ cofactor of $A$.
"""

# ╔═╡ 4ce46738-1c6a-479d-9fec-581c32cb5380
md"""
### Theorem 7.29: Formula for the inverse

Let $A$ be an $n × n$-matrix.
Then

$A \; \text{adj}(A) = \text{adj}(A) \;A = \det(A) I.$

Moreover, $A$ is invertible if and only if $\det(A) ≠ 0$.
In this case, we have:

$A^{-1} = \frac{1}{\det(A)} \text{adj}(A).$

We call this the **adjugate formula** for the matrix inverse.
"""

# ╔═╡ a273bbe6-f563-4bdd-b3dc-8bb5311ef4cc
md"## Application: Cramer's rule"

# ╔═╡ 30696214-55a8-4577-a38d-0be73eb85d49
md"""
### Outcomes

A. Use Cramer's rule to solve a system of equations with invertible coefficient matrix.
"""

# ╔═╡ ef5fb6e2-0761-4d4e-8463-8f8059b0eb2f
md"""
### Theorem 7.34: Cramer's rule

Suppose $A$ is an invertible $n × n$-matrix and we wish to solve the system $A𝐱 = 𝐛$, where $𝐱 = \begin{bmatrix} x_1, …, x_n \end{bmatrix}$.
Then $x_i$ can be computed by the rule

$x_i = \frac{\det(A_i)}{\det(A)},$

where $A_i$ is the matrix obtained by replacing the $i^{\text{th}}$ column of $A$ with $𝐛$.
"""

# ╔═╡ d03356a5-706b-447e-8fc8-d8dbde247277
md"# Eigenvalues, eigenvectors, and diagonalization"

# ╔═╡ 6c19d1a2-4012-4925-8145-93cf8d88cd01
md"## Eigenvectors and eigenvalues"

# ╔═╡ 38574c31-f3d9-4023-8bfd-58aa7ba2b1d1
md"""
### Outcomes

A. Determine whether a vector is an eigenvector of a matrix.

B. Given an eigenvector, find the corresponding eigenvalue.

C. Given an eigenvalue, find the corresponding eigenvectors.

D. Find a basis for the eigenspace of a given eigenvalue.
"""

# ╔═╡ b1cf5cf6-cd9e-43e0-8a50-4bbbb141a50c
md"""
### Definition 8.1: Eigenvalues and eigenvectors

Let $A$ be an $n × n$-matrix.
Suppose that $𝐯 ∈ ℝ^n$ is a non-zero vector such that $A𝐯$ is a scalar multiple of $𝐯$.
In other words, suppose that there exists a scalar $\lambda$ such that

$A𝐯 = \lambda 𝐯$

Then $𝐯$ is called an **eigenvector** of $A$, and $\lambda$ is called the corresponding **eigenvalue**.
"""

# ╔═╡ 11a03a13-e039-47ef-8434-504751b95fcb
md"""
### Definition 8.4: Eigenspace

Let $A$ be an $n × n$-matrix, and let $\lambda$ be an eigenvalue of $A$.
The **eigenspace** of $\lambda$ is the set

$E_\lambda = \{𝐯 ∣ A𝐯 = \lambda 𝐯\}$

It is a subspace of $ℝ^n$.
"""

# ╔═╡ 7f6cda4c-6da5-497b-b82c-6807a8aee2b5
md"## Finding eigenvalues"

# ╔═╡ 2ae68b7f-c082-4a53-a3f3-038ef191ba17
md"""
### Outcomes

A. Find the characteristic polynomial, eigenvalues, and eigenvectors of a matrix.

B. Find the eigenvalues of a triangular matrix.
"""

# ╔═╡ 3b9596ec-45b7-4960-a38f-7406d12ee53a
md"""
### Theorem 8.6: Eigenvalues

Let $A$ be a square matrix, and let $\lambda$ be a scalar.
Then $\lambda$ is an eigenvalue of $A$ if and only if

$\det(A - \lambda I) = 0.$
"""

# ╔═╡ 0b599174-38b3-4140-b05e-7e0f8461b284
md"""
### Definition 8.10: Characteristic polynomial

Let $A$ be a square matrix.
The expression

$p(\lambda) = \det(A - \lambda I)$

is called the **characteristic polynomial** of $A$.
"""

# ╔═╡ fd5407f1-ea60-4b63-abf8-ee7bbdbcfc70
md"""
### Procedure 8.12: Finding eigenvalues and eigenvectors

Let $A$ be an $n × n$-matrix.
To find the eigenvalues and eigenvectors of $A$:

1. Calculate the characteristic polynomial $\det(A - \lambda I)$.

2. The eigenvalues are the roots of the characteristic polynomial.

3. For each eigenvalue $\lambda$, find a basis for the eigenvectors by solving the homogeneous system

$(A - \lambda I)𝐯 = 𝟎.$

To double-check your work, make sure that $A𝐯 = \lambda 𝐯$ for each eigenvalue $\lambda$ and associated eigenvector $𝐯$.
"""

# ╔═╡ 1d0c2282-23c6-4e50-a4b3-e8e9e87f580b
md"""
### Proposition 8.16: Eigenvalues of a triangular matrix

Let $A$ be an upper or lower triangular matrix.
Then the eigenvalues of $A$ are the entries on the main diagonal.
"""

# ╔═╡ b6cdc189-6910-4152-9607-540bdfd1dffc
md"## Geometric interpretation of eigenvectors"

# ╔═╡ 153be9b2-cab0-4bbf-8769-ee539c12eaf6
md"""
### Outcomes

A. Visualize the effect of a linear transformation by considering its eigenvectors and eigenvalues.
"""

# ╔═╡ d65154ad-8298-4e18-aa17-50958e4b9a05
md"## Diagonalization"

# ╔═╡ 82aa67fa-95da-406c-9e3c-a3ed4b4bac16
md"""
### Outcomes

A. Compute sums, products, and powers of diagonal matrices.

B. Determine whether a square matrix is diagonalizable, and diagonalize it if possible.
"""

# ╔═╡ 06b27e13-053e-4b12-bb1d-c19e696c5571
md"""
### Diagonal matrix

A square matrix $D$ is called a **diagonal matrix** if all entries except those on the main diagonal are zero.
Such matrices look like the following:

$\begin{bmatrix} d_{11} & 0 & ⋯ & 0 \\ 0 & d_{22} & ⋯ & 0 \\ ⋮ & ⋮ & ⋱ & ⋮ \\ 0 & 0 & ⋯ & d_{nn} \end{bmatrix}$
"""

# ╔═╡ 1ba30c58-59b6-4b72-8473-89836786fb3e
md"""
### Definition 8.20: Diagonalizable matrix

Let $A$ be an $n × n$-matrix.
Then $A$ is said to be **diagonalizable** if there exists an invertible matrix $P$ and a diagonal matrix $D$ such that

$P^{-1}AP = D.$
"""

# ╔═╡ b446e558-e495-4b41-ad39-311d3284d140
md"""
### Theorem 8.21: Diagonalization and eigenvectors

An $n × n$-matrix $A$ is diagonalizable if and only if $A$ has $n$ linearly independent eigenvectors.

Moreover, in this case, let $P$ be the invertible matrix whose columns are $n$ linearly independent eigenvectors of $A$, and let $D$ be the diagonal matrix whose diagonal entries are the corresponding eigenvalues.
Then $P^{-1} AP = D$.
"""

# ╔═╡ f96db188-18fc-4310-a9a8-8cefdce80b91
md"## Application: Matrix powers"

# ╔═╡ 2df1c6f0-0e5e-4db6-ac99-568b5a55e575
md"""
### Outcomes

A. Use diagonalization to raise a matrix to a high power.

B. Use diagonalization to compute a square root of a matrix.
"""

# ╔═╡ 7563116c-78e8-4d38-8f49-fed240c130d3
md"## Application: Solving recurrences"

# ╔═╡ 9528b092-0411-48af-88d0-15b12fd6edca
md"""
### Outcomes

A. Solve a linear recurrence relation using diagonalization
"""

# ╔═╡ 1bf608a5-90f4-46c6-b248-450a959aa671
md"## Application: Systems of linear differential equations"

# ╔═╡ 6a9232c8-0198-4e94-a482-b6c8fe163872
md"""
### Outcomes

A. Solve a system of second order linear differential equations.
"""

# ╔═╡ 73b31e91-4ef2-4fa8-8197-5cd78f86d569
md"""
### Theorem 8.32: Solutions of $y'' = qy$

Let $q$ be a positive real number.
The differential equation

$\begin{align*}
y'' &= -qy &\text{ has basic solutions }&& y &= \sin(\sqrt{q}x) &\text{ and }&& y &= \cos(\sqrt{q} x), \\
y'' &= 0 &\text{ has basic solutions }&& y &= 1 &\text{ and }&& y &= x, \\
y'' &= qy &\text{ has basic solutions }&& y &= e^{\sqrt{q} x} &\text{ and }&& y &= e^{-\sqrt{q} x}, \\
\end{align*}$
"""

# ╔═╡ 9263dc6b-a174-45a5-8064-c1b7e7c1583d
md"## Application: The matrix exponential"

# ╔═╡ bcee0296-f416-4758-aa67-357f9ac7f29e
md"""
### Outcomes

A. Compute $e^A$, $\sin{A}$, and $\cos{A}$, for a square matrix $A$.

B. Apply any analytic function to a square matrix.
"""

# ╔═╡ e37d13f7-94b0-4ea8-ad9e-422f32934df3
md"""
### Theorem 8.36: Matrix functions by diagonalization

Suppose $A$ is a diagonalizable square matrix, with $A = PDP^{-1}$.
Then

$\begin{align*}
e^A &= Pe^D P^{-1}, \\
\sin{A} &= P(\sin{D}) P^{-1}, \\
\cos{A} &= P(\cos{D}) P^{-1}.
\end{align*}$
"""

# ╔═╡ e3b87191-a4e3-4d76-be69-19254fef192b
md"## Properties of eigenvectors and eigenvalues"

# ╔═╡ d9022444-9aef-42da-8b0d-82ba37d0230d
md"""
### Outcomes

A. Know that eigenvectors corresponding to distinct eigenvalues are linearly independent.

B. Compute the algebraic and geometric multiplicity of an eigenvalue.

C. Determine whether a matrix is diagonalizable from the geometric multiplicities of its eigenvalues.
"""

# ╔═╡ f6ac9d37-8369-43e8-b4e4-377e62b90ab3
md"""
### Proposition 8.38: Eigenvectors for different eigenvalues are linearly independent

Let $A$ be a square matrix, and suppose that $A$ distinct eigenvalues $\lambda_1, …, \lambda_k$ with corresponding eigenvectors $𝐯_1, …, 𝐯_k$.
Then $𝐯_1, …, 𝐯_k$ are linearly independent.
"""

# ╔═╡ 91fc3fa2-9ce6-4f13-9286-591b4ae0d98b
md"""
### Corollary 8.39: Distinct eigenvalues

Let $A$ be an $n × n$-matrix and suppose it has $n$ distinct eigenvalues.
Then $A$ is diagonalizable.
"""

# ╔═╡ 6b917b12-8020-4e2e-9a75-46d2c90a8172
md"""
### Definition 8.40: Algebraic and geometric multiplicity

Let $\hat{\lambda}$ be an eigenvalue of a square matrix $A$.
Then the **algebraic multiplicity** of $\hat{\lambda}$ is power such that $(\hat{\lambda} - \lambda)^k$ is a factor of the characteristic polynomial.
The **geometric multiplicity** of $\hat{\lambda}$ is the dimension of its eigenspace $E_{\lambda}$.
"""

# ╔═╡ 24108eeb-4dfd-4cd9-8857-3053523e4faf
md"""
### Proposition 8.42: Algebraic and geometric multiplicity

Let $\hat{\lambda}$ be an eigenvalue of matrix $A$, with algebraic multiplicity $k$ and geometric multiplicity $m$.
Then

$1 ≤ m ≤ k.$
"""

# ╔═╡ 800d136c-23bb-416d-9fe9-f832470a15be
md"""
### Proposition 8.43: Geometric multiplicity and diagonalization

An $n × n$-matrix $A$ is diagonalizable if and only if the sum of the geometric multiplicities of all the eigenvalues of $A$ is $n$.
"""

# ╔═╡ d6c26b03-76a5-4282-8ba2-e84b03ccc77b
md"## The Cayley-Hamilton Theorem"

# ╔═╡ f3540f07-1cca-4799-b44c-d42a7eb6224f
md"""
### Outcomes

A. For a square matrix $A$, find a polynomial $p(x)$ such that $p(A) = 0$.
"""

# ╔═╡ d8e684b6-27cb-49a7-9abf-71d3c3d99887
md"""
### Theorem 8.44: Cayley-Hamilton theorem

Let $A$ be a square matrix and let $p(\lambda) = \det(A - \lambda I)$ be its characteristic polynomial.
Then $p(A) = 0$.
"""

# ╔═╡ f9bac3a0-c4df-4e31-b0b3-7820fd90cce8
md"""
### Lemma 8.46: Polynomials with matrix coefficients

Let $A_0, …, A_m$ be $n × n$-matrices and assume that for all scalars $\lambda$,

$A_0 + A_1 \lambda + … + A_m \lambda^m = 0.$

Then each $A_i = 0$.
"""

# ╔═╡ 7ab090a5-0a68-4f21-9414-744246d8ed07
md"""
### Corollary 8.47:

Let $A_i$ and $B_i$ be $n × n$-matrices and suppose that

$A_0 + A_1 \lambda + … + A_m \lambda^m = B_0 + B_1 \lambda + … + B_m \lambda^m$

for all $\lambda$.
Then for any $n × n$-matrix $C$,

$A_0 + A_1 C + … + A_m C^m = B_0 + B_1 C + … + B_m C^m.$
"""

# ╔═╡ ae1d7bfd-92d6-4972-b784-d23886748098
md"## Complex eigenvalues and eigenvectors"

# ╔═╡ 42705104-3f52-483d-b333-f75dbd5fb109
md"""
### Outcomes

A. Find the complex eigenvalues and eigenvectors of a matrix.

B. Diagonalize a matrix over the complex numbers.
"""

# ╔═╡ 90913ec4-40e9-410d-8a1e-f860cc17ebb1
md"""
### Proposition 8.51: Complex conjugate eigenvalues

Let $A$ be a square matrix whose entries are real numbers.
If $\lambda$ is an eigenvalue of $A$, then so is $\bar{\lambda}$.
"""

# ╔═╡ fe143e38-ce38-4086-af09-0b0100d7a447
md"""
### Proposition 8.52: Diagonalizabilty criterion

A square matrix $A$ is diagonalizable over the complex numbers if and only if the geometric multiplicity of each eigenvalue is equal to its algebraic multiplicity.
"""

# ╔═╡ 84d76708-e6a3-4cf0-ab47-d11c563b1ead
md"# Vector spaces"

# ╔═╡ 2edd043b-84ef-4b53-b9d5-37756de1f7a9
md"## Definition of vector spaces"

# ╔═╡ f9d868bb-eecd-468d-914e-d4879fa39e3d
md"""
### Outcomes

A. Develop the concept of a vector space through axioms.

B. Use the vector space axioms to determine if a set and its operations constitute a vector space.

C. Encounter several examples of vector spaces.
"""

# ╔═╡ 551ab18d-46a1-4fcc-ad53-438a3fe48b08
md"""
### Definition 9.1: Vector space

Let $K$ be a field.
A **vector space** over $K$ is a set $V$ equipped with two operations of **addition** and **scalar multiplication**, such that the following properties hold:

- (A1) Commutative law of addition: $𝐮 + 𝐯 = 𝐯 + 𝐮$.

- (A2) Associative law of addition: $(𝐮 + 𝐯) + 𝐰 = 𝐮 + (𝐯 + 𝐰)$.

- (A3) The existence of an additive unit: there exists an element $𝟎 ∈ V$ such that for all $𝐮, $𝐮 + 𝟎 = 𝐮$.

- (A4) The law of additive inverses: $𝐮 + (-𝐮) = 𝟎$.

- (SM1) The distributive law over vector addition: $k(𝐮 + 𝐯) = k𝐮 + k𝐯$.

- (SM2) The distributive law over scalar addition: $(k + ℓ)𝐮 = k𝐮 + ℓ𝐮$.

- (SM3) The associative law for scalar multiplication: $k(ℓ𝐮) = (kℓ)𝐮$.

- (SM4) The rule for multiplication by one: $1𝐮 = 𝐮$.
"""

# ╔═╡ 88d3aeb5-9d94-4a27-97e4-e75fea240b66
md"""
### Proposition 9.9: Elementary consequences of the vector space axioms

In any vector space, the following are true:

(a) The additive unit is unique. In other words, whenever $𝐮 + 𝐯 = 𝐮$, then $𝐯 = 𝟎$.

(b) Additive inverses are unique. In other words, whenever $𝐮 + 𝐯 = 𝟎$, then $𝐯 = -𝐮$.

(c) $0𝐮 = 𝟎$ for all vectors $𝐮$.

(d) The following **cancellation law** holds: if $𝐮 + 𝐰 = 𝐯 + 𝐰$, then $𝐮 = 𝐯$.
"""

# ╔═╡ 9ee84712-c42b-42c1-8741-1503d8fe67a0
md"## Linear combinations, span, and linear independence"

# ╔═╡ ff1e5919-3bda-4ad2-bde2-309304e35a04
md"""
### Outcomes

A. Determine if a vector is within a given span.

B. Determine if a set is spanning.

C. Determine if a set is linearly independent.
"""

# ╔═╡ 8424d4eb-a0d3-40b3-892e-97788c073514
md"""
### Definition 9.10: Linear combination

Let $V$ be a vector space over a field $K$.
Let $𝐮_1, …, 𝐮_n ∈ V$.
A vector $𝐯 ∈ V$ is called a **linear combination** of $𝐮_1, …, 𝐮_n$ if there exists scalars $a_1, …, a_n ∈ K$ such that

$𝐯 = a_1 𝐮_1 + … + a_n 𝐮_n.$
"""

# ╔═╡ d19e2ad3-d0ce-49c8-8f7e-28f0520dcb17
md"""
### Definition 9.13: Span of a set of vectors

Let $V$ be a vector space over some field $K$, and let $S$ be a set of vectors (i.e., a subset of $V$).
The **span** of $S$ is the set of all linear combinations of elements of $S$.
In symbols, we have

$\text{span } S = \{a_1 𝐮_1 + … + a_k 𝐮_k ∣ 𝐮_1, …, 𝐮_k ∈ S \text{ and } a_1, …, a_k ∈ K\}.$
"""

# ╔═╡ 1f5f99b8-e8e4-4d07-ad41-ff30692a9033
md"""
### Definition 9.17: Linear independence

Let $V$ be a vector space over some field $K$.
A finite set of vectors $\{𝐮_1, …, 𝐮_k\}$ is called **linearly independent** if the equation

$a_1 𝐮_1 + … + a_k 𝐮_k = 𝟎$

has only the trivial solution $a_1, …, a_k = 0$.
An infinite set $S$ of vectors is called linearly independent if every finite subset of $S$ is linearly independent.
A set of vectors is called **linearly dependent** if it is not linearly independent.
"""

# ╔═╡ 4ce6a89d-2139-4c0a-9b8b-7efa0e4bf7cb
md"""
### Proposition 9.22: Linear dependence and redundant vectors

Let $V$ be a vector space, and let $𝐮_1, 𝐮_2, …$ be a (finite or infinite) sequence of vectors in $V$.
If $𝐮_1, 𝐮_2, …$ are linearly dependent, then at least one of the vectors can be written as a linear combination of earlier vectors in the sequence:

$𝐮_j = a_1 𝐮_1 + a_2 𝐮_2 + … + a_{j - 1} 𝐮_{j - 1},$

for some $j$.
"""

# ╔═╡ d78405ce-f2ac-46a6-9086-b8a1d18d532f
md"""
### Proposition 9.24: Adding to a linearly independent set

Suppose $\{𝐮_1, …, 𝐮_k\}$ is linearly independent and $𝐯 ∉ \text{span} \{𝐮_1, …, 𝐮_k\}$.
Then the set

$\{𝐮_1, …, 𝐮_k, 𝐯\}$

is also linearly independent.
"""

# ╔═╡ 8f16fcc6-df0c-4c98-a7a1-2f44755258d2
md"## Subspaces"

# ╔═╡ e4b6cbcc-bc7c-4005-9bc1-1d84fe6cf28f
md"""
### Outcomes

A. Determine whether a set of vectors is a subspace of a given vector space.

B. Determine whether two sets of vectors span the same subspace.
"""

# ╔═╡ 075cd9b2-e6de-48aa-ad23-a7d18025ccec
md"""
### Definition 9.25: Subspace

Let $V$ be a vector space over a field $K$.
A subset $W ⊆ V$ is said to be a **subspace** of $V$ if the following conditions hold:

1. ``𝟎 ∈ W``, where $𝟎$ is the additive unit of $V$.

2. ``W`` is **closed under addition**: Whenever $𝐮, 𝐯 ∈ W$, then $𝐮 + 𝐯 ∈ W$.

3. ``W`` is **closed under scalar multiplication**: Whenever $k ∈ K$ and $𝐮 ∈ W$, then $k𝐮 ∈ W$.
"""

# ╔═╡ d77fea95-8d56-44e5-a419-e118df68d8d4
md"""
### Proposition 9.33: Subspaces are vector spaces

Let $W$ be a subspace of a vector space $V$.
Then $W$ satisfies the vector space axioms (A1)--(A4) and (SM1)--(SM4), with respect to the same operations (addition and scalar multiplication) as those defined on $V$.
"""

# ╔═╡ 702e7172-b8cf-40bd-bba6-3aba7e5cd327
md"""
### Proposition 9.34: Span is smallest subspace containing given vectors

Let $V$ be a vector space over some field $K$, and consider a set of vectors $S ⊆ V$.
Then $\text{span } S$ is the smallest subspace of $V$ containing $S$.
More explicitly, we have:

(a) The set $\text{span } S$ is a subspace of $V$, and. $S ⊆ \text{span } S$.

(b) If $W$ is any other subspace of $V$ such that $S ⊆ W$, then $\text{span } S ⊆ W$.
"""

# ╔═╡ f1b83797-fca2-490b-b4d2-63368920696d
md"## Basis and dimension"

# ╔═╡ 2122657c-1eca-44bc-b771-272c14b883db
md"""
### Outcomes

A. Find a basis of a given vector space.

B. Determine the dimension of a vector space.

C. Extend a linearly independent set of vectors to a basis.

D. Shrink a spanning set of vectors to a basis.
"""

# ╔═╡ de7b4001-c21d-4aff-aad4-a19f821c1722
md"""
### Definition 9.36: Basis

Let $V$ be a vector space.
A set $B$ of vectors is called **basis** of $V$ if

1. ``B`` is a spanning set for $V$, and

2. ``B`` is linearly independent.
"""

# ╔═╡ fa9315e0-3f85-4f6e-ab08-1956ce173a4d
md"""
### Theorem 9.40: Existence of bases

Every vector space has a basis.
"""

# ╔═╡ 0cafaf89-4a78-4d10-b97b-5b9ba67afda0
md"""
### Lemma 9.41: Exchange Lemma

Let $V$ be a vector space over a field $K$.
Suppose $𝐮_1, …, 𝐮_r$ are linearly independent elements of $\text{span} \{𝐯_1, …, 𝐯_s\}$.
Then $r ≤ s$.
"""

# ╔═╡ c9e48618-a788-4274-894f-a4f5db700991
md"""
### Theorem 9.42: Bases are of the same size

Let $V$ be a vector space over some field $K$, and let $B_1$ and $B_2$ be bases of $V$.
Then either $B_1$ and $B_2$ are both finite and have the same number of elements, or else $B_1$ and $B_2$ are both infinite.
"""

# ╔═╡ 16dd29dd-7b63-48b0-956f-3f1ae2fbb405
md"""
### Definition 9.43: Dimension

Let $V$ be a vector space over a field $K$.
If $V$ has a basis consisting of $n$ vectors, we say that $V$ has **dimension** $n$, and we write $\dim(V) = n$.
In this case we also say that $V$ is **finite-dimensional**.
If $V$ has an infinite basis, we say that $V$ is **infinite-dimensional**, and we write $\dim(V) = ∞$.
"""

# ╔═╡ a16acb0e-da2a-404c-8468-74a0f64cac09
md"""
### Proposition 9.46: Extending a linearly independent set to a basis

Let $V$ be a vector space, and let $S ⊆ V$ be a linearly independent set of vectors.
Then $S$ can be extended to a basis of $V$, i.e., there exists a basis of $B$ of $V$ such that $S ⊆ B$.
"""

# ╔═╡ 87cd5b60-f6d3-46fb-93dd-6fc4a4c8b8cd
md"""
### Proposition 9.47: Shrinking a spanning set to a basis

Let $V$ be a vector space, and let $S ⊆ V$ be a spanning set of $V$.
Then $S$ can be shrunk tot a basis, i.e., there exists a basis $B$ of $V$ such that $B ⊆ S$.
"""

# ╔═╡ 5caf1942-3af4-4f8a-a135-2ecce1c3968a
md"## Application: Error correcting codes"

# ╔═╡ 0fbe99d6-5596-460c-b080-727fb5b27bec
md"""
### Outcomes

A. Determine the block length, message length, Hamming distance, and rate of a code.

B. Determine whether a code is $m$-error detecting or $m$-error correcting.

C. Find generator and check matrices for a linear code.

D. Use the syndrome method to correct errors in code blocks.

E. Construct and use Hamming codes.
"""

# ╔═╡ 4f3a9db6-9114-42fc-9fca-0ab061855756
md"""
### Definition 9.51: Code

Let $n$ and $k$ be positive integers with $k ≤ n$.
A **code** with **message length** $k$ and **block length** $n$ is a set $C$ of $2^k$ different vectors in $ℤ_2^n$.
We call the elements of $C$ the **code blocks**.
Since there are $2^k$ different code blocks, they can be used to encode **message blocks** of length $k$.
We say that the **rate** of the code is $\frac{k}{n}$, because for every $n$ bits of encoded data transmitted, $k$ bits of decoded data are obtained.
"""

# ╔═╡ 18b7af82-28b7-4ed6-89ba-66d86429e4a3
md"""
### Definition 9.53: Hamming weight and Hamming distance

- Let $𝐯$ be a vector in $ℤ_2^n$.
  The **Hamming weight** of $𝐯$, denoted $W(𝐯)$, is the number of components of $𝐯$ that are equal to 1.

- Let $𝐯, 𝐰$ be two vectors in $ℤ_2^n$. The **Hamming distance** between $𝐯$ and $𝐰$, denoted $D(𝐯, 𝐰)$, is the number of components where $𝐯$ and $𝐰$ differ.
  We can also express this as the Hamming weight of $𝐯 - 𝐰$, i.e., $D(𝐯, 𝐰) = W(𝐯 - 𝐰)$.

- Finally, we say that the **Hamming distance of a code** is equal to the smallest Hamming distance between any two code blocks.

A code with block length $n$, message length $k$, and Hamming distance $d$ is also called an **$(n, k, d)$-code**.
"""

# ╔═╡ 48959034-7519-4824-a3b8-e6688a393996
md"""
### Definition 9.55: Error detection and error correction

- We say that a code $C$ is **$m$-error detecting** if for all valid code blocks $𝐯$, whenever $𝐰$ is obtained from $𝐯$ by introducing up to $m$ bit errors, then $𝐰$ is not a valid code block.

- We say that a code $C$ is **$m$-error correcting** if for all valid code blocks $𝐯$, whenever $𝐰$ is obtained from $𝐯$ by introducing up to $m$ bit errors, then $𝐯$ is the only valid code block within Hamming distance $m$ of $𝐰$.
"""

# ╔═╡ 1552bc09-d075-49ee-9f01-229b699f6a66
md"""
### Proposition 9.56: Error detection and error correction

Consider a code with Hamming distance $d$.
Then the code is:

- ``m``-error detecting if $m ≤ d - 1$;

- ``m``-error correcting if $2m ≤ d - 1$.
"""

# ╔═╡ 2fa1aac1-92f4-4997-9fce-27122ba7afa8
md"""
### Definition 9.59: Linear code

A **linear code** is a subspace of $ℤ_2^n$.
If the subspace is $k$-dimensional, then the code has message length $k$ and block length $n$.
"""

# ╔═╡ fab4aed3-5569-4a00-a539-22ffd6ec944f
md"""
### Definition 9.61: Generator matrix and check matrix

Consider a linear code $C$ with block length $n$ and message length $k$, i.e., a $k$-dimensional subspace of $ℤ_2^n$.

- A **generator matrix** for the code is an $n × k$-matrix $G$ such that $C$ is the column space of $G$.

- A **check matrix** for the code is an $(n - k) × n$-matrix $H$ such that $C$ is the null space of $H$.
"""

# ╔═╡ 7778a9ff-7a41-4fe1-8ea0-c7b985e07da4
md"""
### Theorem 9.66: Check matrices of 1-error correcting codes

Consider a linear code with check matrix $H$.
Then the code is 1-error correcting if and only if all columns of the check matrix are non-zero and distinct.
"""

# ╔═╡ ac068c9d-b2f9-47f3-83f3-8e1cb80836fd
md"""
### Definition 9.67: Hamming code

Let $r ≥ 2$.
The $r^{\text{th}}$ **Hamming code** is a linear code whose check matrix $H$ is an $r × (2^r - 1)$-matrix that has all possible non-zero column vectors as its columns.
"""

# ╔═╡ cda437ab-80fd-4ae2-b7ce-344dadec1a0f
md"# Linear transformation of vector spaces"

# ╔═╡ 4f976b2c-f0d9-4b20-b62f-ce11cddb9648
md"## Definition and examples"

# ╔═╡ bdc9319f-244f-4ea6-8ce8-20de224977e7
md"## The algebra of linear transformations"

# ╔═╡ e55fe58e-a849-441a-975a-12e856c3e9f0
md"## Linear transformations defined on a basis"

# ╔═╡ fe1d9c14-f23e-42e1-b15d-7387271ef1a1
md"## The matrix of a linear transformation"

# ╔═╡ d11828b5-7470-4c00-b29e-fb4af280db91
md"# Inner product spaces"

# ╔═╡ 5d826c7a-1ae2-4b79-864e-271a485c0536
md"## Real inner product spaces"

# ╔═╡ 523d9311-5d10-452d-8148-b20c10a4b40e
md"""
### Outcomes

A. Check whether an operation is an inner product.

B. Give an example of a vector space on which more than one inner product can be defined.

C. Calculate the inner product of vectors in various examples of inner product spaces.

D. Calculate the norm of a vector and the angle between two vectors.

E. Use the Cauchy-Schwarz inequality and the triangle inequality to reason about the size of inner products and norms of vectors.
"""

# ╔═╡ 5d8c8815-6484-425c-b623-f67bce72953c
md"""
### Definition 11.1: Real Inner product space

A **(real) inner product space** is a real vector space $V$ equipped with an operation that assigns to any pair of vectors $𝐮, 𝐯 ∈ V$ a real number $⟨𝐮, 𝐯⟩$, called the **inner product** of $𝐮$ and $𝐯$.
This operation must satisfy the following properties:

1. Symmetry: $⟨𝐮, 𝐯⟩ = ⟨𝐯, 𝐮⟩$.

2. Linearity: $⟨𝐮, k𝐯 + ℓ𝐰⟩ = k⟨𝐮, 𝐯⟩ + ℓ⟨𝐮, 𝐰⟩$.

3. The positive definite property: $⟨𝐮, 𝐮⟩ ≥ 0$, and moreover, $⟨𝐮, 𝐮⟩ = 0$ if and only if $𝐮 = 𝟎$.
"""

# ╔═╡ fb1b1273-1eda-443a-8aca-8f62f74a372e
md"""
### Proposition 11.7: Left linearity

Let $V$ be an inner product space.
Then $⟨k𝐮 + ℓ𝐯, 𝐰⟩ = k⟨𝐮, 𝐰⟩ + ℓ⟨𝐯, 𝐰⟩$.
"""

# ╔═╡ bb565b5d-c912-4d89-aa56-24b2c29343c4
md"""
### Theorem 11.8: Cauchy-Schwarz inequality

Let $𝐮, 𝐯$ be vector in an inner product space.
Then

$⟨𝐮, 𝐯⟩^2 ≤ ⟨𝐮, 𝐮⟩ ⋅ ⟨𝐯, 𝐯⟩.$
"""

# ╔═╡ b0e78f06-d995-404b-9e0c-00b193b607bf
md"""
### Definition 11.9: Norm

Let $𝐮$ be a vector in an inner product space.
Then the **norm** of $𝐮$ is defined as

$\|𝐮\| = \sqrt{⟨𝐮, 𝐮⟩}.$
"""

# ╔═╡ 222808aa-3620-4bbd-9ad3-878677c2ca32
md"""
### Proposition 11.11: Inner product and norm

Let $𝐮, 𝐯$ be vectors in an inner product space.
Then

$|⟨𝐮, 𝐯⟩| ≤ \|𝐮\|\|𝐯\|.$
"""

# ╔═╡ 89ac4300-da7e-4ad3-8488-62ca2d85bafd
md"""
### Proposition 11.12: Triangle inequality

Let $𝐮, 𝐯$ be vectors in an inner product space.
Then

$\|𝐮 + 𝐯\| ≤ \|𝐮\| + \|𝐯\|.$

This is called the **triangle inequality**.
"""

# ╔═╡ e87ddc4a-599c-43b2-8d03-9ea8250326f5
md"""
### Definition 11.13: Angle between two vectors

Let $𝐮, 𝐯$ be vectors in an inner product space.
The **angle** between $𝐮$ and $𝐯$ is, by definition, the unique $\theta$ such that $0 ≤ \theta ≤ \pi$ and

$\cos{\theta} = \frac{⟨𝐮, 𝐯⟩}{\|𝐮\|\|𝐯\|}.$
"""

# ╔═╡ 81fc1a03-7eeb-46a7-98d2-667e20414750
md"## Orthogonality"

# ╔═╡ 0b08309f-7437-4b86-bd8b-caa380f7e56e
md"""
### Outcomes

A. Determine whether two vectors in an inner product space are orthogonal.

B. Find the orthogonal complement of a set of vectors.

C. Determine whether a set of vectors is orthogonal and/or orthonormal.

D. Check whether a basis is orthogonal and/or orthonormal.

E. Calculate the Fourier coefficients of a vector with respect to an orthogonal basis.
"""

# ╔═╡ 0ee828dd-3a3d-4572-89d3-9ae8fc2e6522
md"""
### Definition 11.15: Orthogonality

Let $𝐮$ and $𝐯$ be vectors in an inner product space.
We say that $𝐮$ and $𝐯$ are **orthogonal** if

$⟨𝐮, 𝐯⟩ = 0.$

We also write $𝐮 \perp 𝐯$ to indicate that $𝐮$ and $𝐯$ are orthogonal.
"""

# ╔═╡ 90642e5a-e33f-40d8-9b42-e1ad1bc16bd6
md"""
### Definition 11.16: Orthogonal complement

Let $S$ be a subset of an inner product space $V$.
The **orthogonal complement** of $S$ is the set

$S^{\perp} = \{𝐯 ∈ V ∣ ⟨𝐯, 𝐰⟩ = 0 \text{ for all } 𝐰 ∈ S\}$
"""

# ╔═╡ cdc19a24-fd40-4734-804c-644b66f2480a
md"""
### Proposition 11.17: Orthogonal complement

If $S$ is any subset of an inner product space $V$, then $S^{\perp}$ is a subspace.
"""

# ╔═╡ d9a093fc-b519-4946-950a-a38ab65fd0c7
md"""
### Definition 11.19: Orthogonal and orthonormal sets of vectors

A set of vectors $\{𝐮_1, …, 𝐮_k\}$ in an inner product space is called an **orthogonal set** if the vectors are non-zero and pairwise orthogonal, i.e., for all $i$, $𝐮_i ≠ 𝟎$ and for all $i ≠ j$, $𝐮_i \perp 𝐮_j$.

Moreover, the set of vectors $\{𝐮_1, …, 𝐮_k\}$ is called **orthonormal** if it is orthogonal and each vector is normalized, i.e., $\|𝐮_i\| = 1$.
"""

# ╔═╡ 4751700f-5abd-4366-89d8-34ece5ac6028
md"""
### Proposition 11.20: Orthogonal set is linearly independent

If $\{𝐮_1, …, 𝐮_k\}$ is an orthogonal set of vectors, then $𝐮_1, …, 𝐮_k$ are linearly independent.
"""

# ╔═╡ 729bdbe9-79f2-426b-9982-a1cf5e9232a4
md"""
### Definition 11.22: Orthogonal and orthonormal bases

Let $V$ be an inner product space, and let $W$ be a subspace of $V$.
We say that $\{𝐮_1, …, 𝐮_k\}$ is an **orthogonal basis** for $W$ if it is an orthogonal set and spans $W$.

If, moreover, each $𝐮_i$ is normalized, we say that $\{𝐮_1, …, 𝐮_k\}$ is an **orthonormal basis** for $W$.
"""

# ╔═╡ bad3a9e1-4b57-48ef-ad54-81a4008691dd
md"""
### Proposition 11.24: Fourier coefficients

Let $B = \{𝐮_1, …, 𝐮_n\}$ be an orthogonal basis of some space $W$, and suppose $𝐯 ∈ W$.
Then

$𝐯 = a_1 𝐮_1 + … + a_n 𝐮_n,$

where

$a_i = \frac{⟨𝐮_1, 𝐯⟩}{⟨𝐮_i, 𝐮_i⟩}.$

In this situation, the coordinates $a_1, …, a_n$ are also called the **Fourier coefficients** of $𝐯$ (with respect to the orthogonal basis $B$).

In case $B$ is an orthonormal basis, the formula is even simpler.
In that case:

$a_i = ⟨𝐮_i, 𝐯⟩.$
"""

# ╔═╡ 9e049b8e-6fa9-4ec3-9a79-8525b3e9f5b7
md"""
### Proposition 11.27: Norm of orthogonal linear combination

Let $V$ be an inner product space, and suppose that $𝐮_1, …, 𝐮_k$ are orthogonal.
Then

$\|a_1 𝐮_1 + … + a_k 𝐮_k\|^2 = {a_1}^2 \|𝐮_1\|^2 + … + {a_k}^2 \|𝐮_k\|^2.$
"""

# ╔═╡ da2c86a7-db2b-4859-92b4-ac3af222725a
md"""
### Proposition 11.28: Inner product and dot product

Let $V$ be an inner product space, and suppose that $B = \{𝐮_1, …, 𝐮_n\}$ is an orthonormal basis.
Then for any pair of vectors $𝐯, 𝐰 ∈ V$, their inner product is equal to the dot product of their coordinate vectors with respect to the basis $B$, i.e.

$⟨𝐯, 𝐰⟩ = [𝐯]_B ⋅ [𝐰]_B.$
"""

# ╔═╡ 06519655-75c4-4330-903a-dff5f0b8d935
md"## The Gram-Schmidt orthogonalization procedure"

# ╔═╡ d3105f5f-d243-4028-8ee6-b62a6c6a72fd
md"""
### Outcomes

A. Use the Gram-Schmidt procedure to find an orthogonal basis of a subspace of an inner product space.

B. Find an orthonormal basis of a subspace.
"""

# ╔═╡ c784a78d-fae9-4930-abb7-c8818cf5f49f
md"""
### Proposition 11.29: Gram-Schmidt orthogonalization procedure for 2 vectors

Let $\{𝐯_1, 𝐯_2\}$ be a basis for some subspace $W$ of an inner product space $V$.
Define vectors $𝐮_1, 𝐮_2$ as follows:

$\begin{align*}
𝐮_1 &= 𝐯_1, \\
𝐮_2 &= 𝐯_2 - \frac{⟨𝐮_1, 𝐯_2⟩}{⟨𝐮_1, 𝐮_1⟩} 𝐮_1.
\end{align*}$

Then $\{𝐮_1, 𝐮_2\}$ is an orthogonal basis of $W$.
"""

# ╔═╡ 1feb211b-0b85-4100-b9d2-9421585d990f
md"""
### Proposition 11.31: Gram-Schmidt orthogonalization procedure for $k$ vectors

Let $\{𝐯_1, …, 𝐯_k\}$ be a basis for some subspace $W$ of an inner product space $V$.
Define vectors $𝐮_1, …, 𝐮_k$ as follows:

$\begin{align*}
𝐮_1 &= 𝐯_1, \\
𝐮_2 &= 𝐯_2 - \frac{⟨𝐮_1, 𝐯_2⟩}{⟨𝐮_1, 𝐮_1⟩} 𝐮_1, \\
𝐮_3 &= 𝐯_3 - \frac{⟨𝐮_1, 𝐯_3⟩}{⟨𝐮_1, 𝐮_1⟩} 𝐮_1 - \frac{⟨𝐮_2, 𝐯_3⟩}{⟨𝐮_2, 𝐮_2⟩} 𝐮_2, \\
&⋮ \\
𝐮_k &= 𝐯_k - \frac{⟨𝐮_1, 𝐯_k⟩}{⟨𝐮_1, 𝐮_1⟩} 𝐮_1 - \frac{⟨𝐮_2, 𝐯_k⟩}{⟨𝐮_2, 𝐮_2⟩} 𝐮_2 - … - \frac{⟨𝐮_{k - 1}, 𝐯_k⟩}{⟨𝐮_{k - 1}, 𝐮_{k - 1}⟩} 𝐮_{k - 1}, \\
\end{align*}$

Then $\{𝐮_1, …, 𝐮_k\}$ is an orthogonal basis of $W$.
"""

# ╔═╡ da118ab7-1e2e-40b4-b1b2-9a1696cd125a
md"## Orthogonal projections and Fourier series"

# ╔═╡ 9e10f35a-db41-4f11-b350-28882237973a
md"## Application: Least squares approximations and curve fitting"

# ╔═╡ c1899059-6569-4b36-b6ef-218c0f4bd8c0
md"## Orthogonal functions and orthogonal matrices"

# ╔═╡ 7c082dfb-d073-4370-8242-f1cca35f2557
md"## Diagonalization of symmetric matrices"

# ╔═╡ 89463c44-2ece-4556-a6d1-d381155a84e4
md"## Positive semidefinite and positive definite matrices"

# ╔═╡ 5642769f-3467-411b-906d-dc7b49cbd1a1
md"## Application: Simplification of quadratic forms"

# ╔═╡ 2c625a46-6fd3-43da-a577-fbc658b1a7a7
md"## Complex inner product spaces"

# ╔═╡ 20bb4627-80c7-4dd6-8047-bdd6ecccaf33
md"## Application: Principal component analysis"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Plots = "~1.23.2"
PlutoUI = "~0.7.17"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "3533f5a691e60601fe60c90d8bc47a27aa2907ec"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.0"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "a851fec56cb73cfdf43762999ec72eff5b86882a"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.15.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dce3e3fea680869eaa0b774b2e8343e9ff442313"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.40.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

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

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "0c603255764a1fa0b61752d2bec14cfbd18f7fe8"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+1"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "d189c6d2004f63fd3c91748c458b09f26de0efaa"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.61.0"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fd75fa3a2080109a2c0ec9864a6e14c60cca3866"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.62.0+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "14eece7a3308b4d8be910e265c724a6ba51a9798"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.16"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "5efcf53d798efede8fee5b2c8b09284be359bf24"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.2"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "f0c6489b12d28fb4c2103073ec7452f3423bd308"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.1"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

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

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "6193c3815f13ba1b78a51ce391db8be016ae9214"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.4"

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

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "d911b6a12ba974dabe2291c6d450094a7226b372"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.1"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "b084324b4af5a438cd63619fd006614b3b20b87b"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.15"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "ca7d534a27b1c279f05cd094196cb70c35e3d892"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.23.2"

[[PlutoUI]]
deps = ["Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "615f3a1eff94add4bca9476ded096de60b46443b"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.17"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

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

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "7ad0dfa8d03b7bcf8c597f59f5292801730c55b8"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.4.1"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

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
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "eb35dcc66558b2dda84079b9a1be17557d32091a"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.12"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "2ce41e0d042c60ecd131e9fb7154a3bfadbf50d3"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.3"

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

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll"]
git-tree-sha1 = "2839f1c1296940218e35df0bbb220f2a79686670"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.18.0+4"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─7b8bf49e-fb76-4034-a400-22f90e1983a3
# ╟─542aec7c-3740-4bba-a572-87971b5614fe
# ╟─a5c5f744-82df-4db9-a958-1a9f50005072
# ╠═73440efb-2a36-414f-9f53-23d07b792574
# ╟─830a5a05-7151-4a70-ae28-6353d50290bb
# ╟─f13a313a-0359-4ae7-b0aa-c0616df64234
# ╟─459a45a1-1013-4b5a-afbf-a7c9e9556887
# ╟─6e46f7d2-ae75-48c8-b4ce-131b52c0f795
# ╟─46d7cf27-29b2-48a4-91e5-808bbae92c68
# ╟─d15db133-0b81-401c-9935-5d81bfc6ba9f
# ╟─724b61b0-117a-4132-9431-659dd8e8d5ad
# ╟─d5bdb6e7-185d-4309-a00d-ce3fa544c4c5
# ╟─aa81baa3-fde8-441a-8bdc-7fcdd1d7ff96
# ╟─90873450-e0c4-4e5a-ad1c-a260e6caf825
# ╟─7c17981e-ee11-4e39-bf7d-7cdf94dbb5cf
# ╟─a89b4ad7-31fc-418c-bb95-36187aadff22
# ╟─091bb6cb-e600-4df4-9c5a-58c6de4e48cd
# ╟─be63bec1-3317-457d-b17e-4f3b232dbd08
# ╟─5e0eac80-5837-4a6f-af5f-f587270751fc
# ╟─53f07204-8912-40cb-be72-6678635c731d
# ╟─5b3752a9-9391-4c8a-8daa-32a94c75f9ae
# ╟─f8a521e9-7d13-419f-8fd9-b8f98e7614f2
# ╟─2017c2ed-ebe5-4eb7-8825-44751abe6076
# ╟─41e6ca5e-a95e-42e6-8ab4-b7721028eb2f
# ╟─3d20ee51-69f8-4ba1-b005-383f8410b5d5
# ╟─91baaf1c-479a-46a4-9f75-cca70d08c80a
# ╟─b6d18658-7b6b-4f6a-b9b0-ece83963a470
# ╟─149b568a-70a9-40bd-9628-e3d1485f3216
# ╟─afe0d69e-3d2a-4fab-9daa-6502918cfc85
# ╟─a2dfbd75-53b1-4d1e-9abd-52ca3f165e5d
# ╟─b1158480-08c9-48cc-8f01-30a84437b9e2
# ╟─ffb118ad-8baf-4556-9062-0480de6459f3
# ╟─94d46bac-23fc-413e-aaf8-ee2940233eb6
# ╟─9c28afd7-7203-4d9f-91cb-ca7ae1cfa823
# ╟─203cc5d2-45be-4704-9bfe-9410e2efd3bf
# ╟─da3046e7-7436-4e0e-9fe3-bbf28b15c267
# ╟─3b3bc7a5-6455-41e5-a872-8b7a252e37c6
# ╟─bd75f8f7-c31e-450f-b544-be72ddc60b9a
# ╟─4d65dbd5-e6c8-482d-a119-53c8f59676c8
# ╟─db7ec3f6-6e83-4155-83fc-894d1ba9b521
# ╟─8bded7bf-2fa5-42d3-ba49-650e78df4d8a
# ╟─63fe60cb-df4c-4ba8-9f9c-fb52604a6bd7
# ╟─45ae2f7b-d464-45e3-b416-531dbb78627d
# ╟─09abcb81-943c-42cb-b3e0-492353415fe8
# ╟─25d798ad-879c-4bb4-a538-c26d28f3f486
# ╟─3686f34f-bb50-4e86-be62-62d4b169ca5b
# ╟─8b863a5c-f515-4669-94c4-d64dbe928bdf
# ╟─1d3a1b05-345e-4189-ad4e-fb47c4871804
# ╟─8afc742e-09ab-4ea3-985c-84a7508c27b9
# ╟─d32a3b9d-d156-4bcb-9072-d7e34b6c600c
# ╟─48575c01-ef6d-4db2-9144-6ec4850044e4
# ╟─c42b1c38-6cca-4138-8f54-bddde6b3ad5f
# ╟─fb7cc1b2-6573-4cf9-85b2-ebf6ab59a274
# ╟─508de533-96d5-4817-bbbc-3fc2c41fa848
# ╟─033cca18-dcb6-4720-bf13-5379cd334df2
# ╟─bdb1858d-f3ab-44f9-96da-759c31998d6f
# ╟─48ac0870-e5e3-4581-94ee-182c87b4e78e
# ╟─3cde46c8-3576-4bd4-b5e4-aa8040ab7037
# ╟─dbc16007-9c67-4b5c-9c88-c3680294138a
# ╟─b980e504-1e37-4d17-a246-685fbb9a346e
# ╟─907c43a9-d4e6-4043-bca9-891f4a81ac1d
# ╟─8e304446-19d7-4d39-bdb0-8057b3b4b057
# ╟─60cf4ac0-7825-4385-99fb-0d1c1811a1e0
# ╟─960b204f-5721-4c76-a91a-e61b64c22224
# ╟─b28fe842-847f-45d8-aa8a-5a5a66467578
# ╟─243bc5e3-3ce9-404a-920b-25ac294f92b9
# ╟─075f38de-966e-405a-ae01-b4a3b179a502
# ╟─e32b08fb-def4-43a6-b031-f28ae5b56028
# ╟─432b500e-368b-4a45-a888-5390efec5508
# ╟─0889fdf6-32f5-4566-a9d5-bcabed0bf24f
# ╟─4260cefc-050f-4da1-8257-c4eee5be38fa
# ╟─59c4f604-abf7-4589-802b-eb12f3bc0553
# ╟─959dd506-e742-4f4c-a9b1-450032b9c3c7
# ╟─024b87e3-f48f-46e5-bdca-84f383146d47
# ╟─56778fde-510b-4877-8c0d-7f4a2e39cea3
# ╟─27a52480-d576-4bdf-a063-7e005fae16c9
# ╟─f7de8cf8-ef6a-461c-87a6-046618ea9ec3
# ╟─61cd090b-868f-4ca8-b3a3-198bfaf600a7
# ╟─4736289b-2a7f-415d-b2ac-ae7c15499197
# ╟─51ecfe8d-c4a3-4ba1-a6f1-7f792f2e0933
# ╟─ca305f42-fb86-41a0-ba71-189de923346a
# ╟─4354eea7-5e35-4979-8c44-53a34c1ddffb
# ╟─928d92a1-baa6-4b36-824d-7c56f85c857c
# ╟─b201693c-d02e-4e2c-86b4-3b2db1c5fa8c
# ╟─aa765c1b-1f13-45ef-bbef-851781a64d29
# ╟─f819b979-07b8-47b2-866c-6d77889f3b6d
# ╟─93310b8c-09b7-4f3c-936f-d95586c5f736
# ╟─6a7d7172-34d3-420a-8d36-fb51739e98f4
# ╟─a44deee1-5ecc-47b3-81b3-c30f0b53b2f5
# ╟─1f7f2aa2-5913-4b4d-8c2d-94d76487c690
# ╟─a573916f-280b-4055-9a30-d9aa663d7154
# ╟─6869f242-3f47-4c88-a669-b525bd101e8c
# ╟─4ce46738-1c6a-479d-9fec-581c32cb5380
# ╟─a273bbe6-f563-4bdd-b3dc-8bb5311ef4cc
# ╟─30696214-55a8-4577-a38d-0be73eb85d49
# ╟─ef5fb6e2-0761-4d4e-8463-8f8059b0eb2f
# ╟─d03356a5-706b-447e-8fc8-d8dbde247277
# ╟─6c19d1a2-4012-4925-8145-93cf8d88cd01
# ╟─38574c31-f3d9-4023-8bfd-58aa7ba2b1d1
# ╟─b1cf5cf6-cd9e-43e0-8a50-4bbbb141a50c
# ╟─11a03a13-e039-47ef-8434-504751b95fcb
# ╟─7f6cda4c-6da5-497b-b82c-6807a8aee2b5
# ╟─2ae68b7f-c082-4a53-a3f3-038ef191ba17
# ╟─3b9596ec-45b7-4960-a38f-7406d12ee53a
# ╟─0b599174-38b3-4140-b05e-7e0f8461b284
# ╟─fd5407f1-ea60-4b63-abf8-ee7bbdbcfc70
# ╟─1d0c2282-23c6-4e50-a4b3-e8e9e87f580b
# ╟─b6cdc189-6910-4152-9607-540bdfd1dffc
# ╟─153be9b2-cab0-4bbf-8769-ee539c12eaf6
# ╟─d65154ad-8298-4e18-aa17-50958e4b9a05
# ╟─82aa67fa-95da-406c-9e3c-a3ed4b4bac16
# ╟─06b27e13-053e-4b12-bb1d-c19e696c5571
# ╟─1ba30c58-59b6-4b72-8473-89836786fb3e
# ╟─b446e558-e495-4b41-ad39-311d3284d140
# ╟─f96db188-18fc-4310-a9a8-8cefdce80b91
# ╟─2df1c6f0-0e5e-4db6-ac99-568b5a55e575
# ╟─7563116c-78e8-4d38-8f49-fed240c130d3
# ╟─9528b092-0411-48af-88d0-15b12fd6edca
# ╟─1bf608a5-90f4-46c6-b248-450a959aa671
# ╟─6a9232c8-0198-4e94-a482-b6c8fe163872
# ╟─73b31e91-4ef2-4fa8-8197-5cd78f86d569
# ╟─9263dc6b-a174-45a5-8064-c1b7e7c1583d
# ╟─bcee0296-f416-4758-aa67-357f9ac7f29e
# ╟─e37d13f7-94b0-4ea8-ad9e-422f32934df3
# ╟─e3b87191-a4e3-4d76-be69-19254fef192b
# ╟─d9022444-9aef-42da-8b0d-82ba37d0230d
# ╟─f6ac9d37-8369-43e8-b4e4-377e62b90ab3
# ╟─91fc3fa2-9ce6-4f13-9286-591b4ae0d98b
# ╟─6b917b12-8020-4e2e-9a75-46d2c90a8172
# ╟─24108eeb-4dfd-4cd9-8857-3053523e4faf
# ╟─800d136c-23bb-416d-9fe9-f832470a15be
# ╟─d6c26b03-76a5-4282-8ba2-e84b03ccc77b
# ╟─f3540f07-1cca-4799-b44c-d42a7eb6224f
# ╟─d8e684b6-27cb-49a7-9abf-71d3c3d99887
# ╟─f9bac3a0-c4df-4e31-b0b3-7820fd90cce8
# ╟─7ab090a5-0a68-4f21-9414-744246d8ed07
# ╟─ae1d7bfd-92d6-4972-b784-d23886748098
# ╟─42705104-3f52-483d-b333-f75dbd5fb109
# ╟─90913ec4-40e9-410d-8a1e-f860cc17ebb1
# ╟─fe143e38-ce38-4086-af09-0b0100d7a447
# ╟─84d76708-e6a3-4cf0-ab47-d11c563b1ead
# ╟─2edd043b-84ef-4b53-b9d5-37756de1f7a9
# ╟─f9d868bb-eecd-468d-914e-d4879fa39e3d
# ╟─551ab18d-46a1-4fcc-ad53-438a3fe48b08
# ╟─88d3aeb5-9d94-4a27-97e4-e75fea240b66
# ╟─9ee84712-c42b-42c1-8741-1503d8fe67a0
# ╟─ff1e5919-3bda-4ad2-bde2-309304e35a04
# ╟─8424d4eb-a0d3-40b3-892e-97788c073514
# ╟─d19e2ad3-d0ce-49c8-8f7e-28f0520dcb17
# ╟─1f5f99b8-e8e4-4d07-ad41-ff30692a9033
# ╟─4ce6a89d-2139-4c0a-9b8b-7efa0e4bf7cb
# ╟─d78405ce-f2ac-46a6-9086-b8a1d18d532f
# ╟─8f16fcc6-df0c-4c98-a7a1-2f44755258d2
# ╟─e4b6cbcc-bc7c-4005-9bc1-1d84fe6cf28f
# ╟─075cd9b2-e6de-48aa-ad23-a7d18025ccec
# ╟─d77fea95-8d56-44e5-a419-e118df68d8d4
# ╟─702e7172-b8cf-40bd-bba6-3aba7e5cd327
# ╟─f1b83797-fca2-490b-b4d2-63368920696d
# ╟─2122657c-1eca-44bc-b771-272c14b883db
# ╟─de7b4001-c21d-4aff-aad4-a19f821c1722
# ╟─fa9315e0-3f85-4f6e-ab08-1956ce173a4d
# ╟─0cafaf89-4a78-4d10-b97b-5b9ba67afda0
# ╟─c9e48618-a788-4274-894f-a4f5db700991
# ╟─16dd29dd-7b63-48b0-956f-3f1ae2fbb405
# ╟─a16acb0e-da2a-404c-8468-74a0f64cac09
# ╟─87cd5b60-f6d3-46fb-93dd-6fc4a4c8b8cd
# ╟─5caf1942-3af4-4f8a-a135-2ecce1c3968a
# ╟─0fbe99d6-5596-460c-b080-727fb5b27bec
# ╟─4f3a9db6-9114-42fc-9fca-0ab061855756
# ╟─18b7af82-28b7-4ed6-89ba-66d86429e4a3
# ╟─48959034-7519-4824-a3b8-e6688a393996
# ╟─1552bc09-d075-49ee-9f01-229b699f6a66
# ╟─2fa1aac1-92f4-4997-9fce-27122ba7afa8
# ╟─fab4aed3-5569-4a00-a539-22ffd6ec944f
# ╟─7778a9ff-7a41-4fe1-8ea0-c7b985e07da4
# ╟─ac068c9d-b2f9-47f3-83f3-8e1cb80836fd
# ╟─cda437ab-80fd-4ae2-b7ce-344dadec1a0f
# ╟─4f976b2c-f0d9-4b20-b62f-ce11cddb9648
# ╟─bdc9319f-244f-4ea6-8ce8-20de224977e7
# ╟─e55fe58e-a849-441a-975a-12e856c3e9f0
# ╟─fe1d9c14-f23e-42e1-b15d-7387271ef1a1
# ╟─d11828b5-7470-4c00-b29e-fb4af280db91
# ╟─5d826c7a-1ae2-4b79-864e-271a485c0536
# ╟─523d9311-5d10-452d-8148-b20c10a4b40e
# ╟─5d8c8815-6484-425c-b623-f67bce72953c
# ╟─fb1b1273-1eda-443a-8aca-8f62f74a372e
# ╟─bb565b5d-c912-4d89-aa56-24b2c29343c4
# ╟─b0e78f06-d995-404b-9e0c-00b193b607bf
# ╟─222808aa-3620-4bbd-9ad3-878677c2ca32
# ╟─89ac4300-da7e-4ad3-8488-62ca2d85bafd
# ╟─e87ddc4a-599c-43b2-8d03-9ea8250326f5
# ╟─81fc1a03-7eeb-46a7-98d2-667e20414750
# ╟─0b08309f-7437-4b86-bd8b-caa380f7e56e
# ╟─0ee828dd-3a3d-4572-89d3-9ae8fc2e6522
# ╟─90642e5a-e33f-40d8-9b42-e1ad1bc16bd6
# ╟─cdc19a24-fd40-4734-804c-644b66f2480a
# ╟─d9a093fc-b519-4946-950a-a38ab65fd0c7
# ╟─4751700f-5abd-4366-89d8-34ece5ac6028
# ╟─729bdbe9-79f2-426b-9982-a1cf5e9232a4
# ╟─bad3a9e1-4b57-48ef-ad54-81a4008691dd
# ╟─9e049b8e-6fa9-4ec3-9a79-8525b3e9f5b7
# ╟─da2c86a7-db2b-4859-92b4-ac3af222725a
# ╟─06519655-75c4-4330-903a-dff5f0b8d935
# ╟─d3105f5f-d243-4028-8ee6-b62a6c6a72fd
# ╟─c784a78d-fae9-4930-abb7-c8818cf5f49f
# ╟─1feb211b-0b85-4100-b9d2-9421585d990f
# ╟─da118ab7-1e2e-40b4-b1b2-9a1696cd125a
# ╟─9e10f35a-db41-4f11-b350-28882237973a
# ╟─c1899059-6569-4b36-b6ef-218c0f4bd8c0
# ╟─7c082dfb-d073-4370-8242-f1cca35f2557
# ╟─89463c44-2ece-4556-a6d1-d381155a84e4
# ╟─5642769f-3467-411b-906d-dc7b49cbd1a1
# ╟─2c625a46-6fd3-43da-a577-fbc658b1a7a7
# ╟─20bb4627-80c7-4dd6-8047-bdd6ecccaf33
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
