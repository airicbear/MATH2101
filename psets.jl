### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ 35bb1e85-5c56-48e5-a8ea-92ffb0ca6148
begin
	using RowEchelon
	using LinearAlgebra
	using Symbolics
	
	R1(M::Matrix) = M[1:1,:]
	R2(M::Matrix) = M[2:2,:]
	R3(M::Matrix) = M[3:3,:]
	R(r::UnitRange, M::Matrix) = M[r,:]
	
	minor(M::Matrix, i::Int, j::Int) = M[1:end .!= i, 1:end .!= j]
	minor(M::Matrix) = reduce(hcat, [[det(minor(M, i, j)) for i ∈ 1:size(M, 1)] for j ∈ 1:size(M, 2)])
	minors(M::Matrix) = [minor(M, i, j) for i ∈ 1:size(M, 1) for j ∈ 1:size(M, 2)]
	
	cof(M::Matrix, i::Int, j::Int) = (-1)^(i + j) * det(minor(M, i, j))
	cof(M::Matrix) = reduce(hcat, [[cof(M, i, j) for i ∈ 1:size(M, 1)] for j ∈ 1:size(M, 2)])

	md"# Linear Algebra Problem Sets"
end

# ╔═╡ 0727e9d6-c906-410c-b7c3-b7f2922b3309
md"""
## Problem Set 1

1. Review pages 1--23.

2. Find an echelon form of the given matrix (Gaussian elimination):

$\begin{align*}
(\text{A}) &&\qquad \begin{bmatrix}
2 & 3 & -3 & 1 \\
1 & -3 & 2 & 1 \\
-5 & 1 & 1 & -3
\end{bmatrix} \\

(\text{B}) &&\qquad \begin{bmatrix}
2 & 3 & -3 & 1 \\
1 & -3 & 2 & 1 \\
-3 & -9 & 8 & -4
\end{bmatrix} \\

(\text{C}) &&\qquad \begin{bmatrix}
2 & 1 & -3 & 2 \\
1 & -3 & 2 & 1
\end{bmatrix}
\end{align*}$

3. Each of these matrices is the augmented matrix of a linear system of equations in the unknowns $x$, $y$, $z$.
Write down the three systems, and, when possible, use the information gained from the echelon forms to write down solutions.
(Note: One system has exactly one solution, one has many, one has none.)
"""

# ╔═╡ 2262a111-83be-4e51-8518-ef50b26fe450
md"""
### PS1 #2 (A)

$\begin{bmatrix}
2 & 3 & -3 & 1 \\
1 & -3 & 2 & 1 \\
-5 & 1 & 1 & -3
\end{bmatrix}$

$\begin{align*}
R_1 = \frac{1}{2} R_1 &\implies 
\begin{bmatrix}
1 & \frac{3}{2} & -\frac{3}{2} & \frac{1}{2} \\
1 & -3 & 2 & 1 \\
-5 & 1 & 1 & -3
\end{bmatrix} \\

R_2 = R_2 - R_1 &\implies
\begin{bmatrix}
1 & \frac{3}{2} & -\frac{3}{2} & \frac{1}{2} \\
0 & -\frac{9}{2} & \frac{7}{2} & \frac{1}{2} \\
-5 & 1 & 1 & -3
\end{bmatrix} \\

R_3 = R_3 + 5R_1 &\implies
\begin{bmatrix}
1 & \frac{3}{2} & -\frac{3}{2} & \frac{1}{2} \\
0 & -\frac{9}{2} & \frac{7}{2} & \frac{1}{2} \\
0 & \frac{17}{2} & -\frac{13}{2} & -\frac{1}{2}
\end{bmatrix} \\

R_2 = -\frac{2}{9}R_2 &\implies
\begin{bmatrix}
1 & \frac{3}{2} & -\frac{3}{2} & \frac{1}{2} \\
0 & 1 & -\frac{7}{9} & -\frac{1}{9} \\
0 & \frac{17}{2} & -\frac{13}{2} & -\frac{1}{2}
\end{bmatrix} \\

R_3 = R_3 - \frac{17}{2} R_2 &\implies
\begin{bmatrix}
1 & \frac{3}{2} & -\frac{3}{2} & \frac{1}{2} \\
0 & 1 & -\frac{7}{9} & -\frac{1}{9} \\
0 & 0 & \frac{1}{9} & \frac{4}{9}
\end{bmatrix} \\

R_3 = 9R_3 &\implies
\begin{bmatrix}
1 & \frac{3}{2} & -\frac{3}{2} & \frac{1}{2} \\
0 & 1 & -\frac{7}{9} & -\frac{1}{9} \\
0 & 0 & 1 & 4
\end{bmatrix}
\end{align*}$
"""

# ╔═╡ 854c658c-c1ba-4dac-8982-17513b89f1df
let
	A = [2 3 -3
		 1 -3 2
		-5 1 1]
	
	b = [1
		 1
		-3]
	
	M = [A b]

	M = [(1//2)R1(M)
    	 R(2:3,M)]

	M = [R1(M)
		 R2(M) - R1(M)
		 R3(M)]
	
	M = [R(1:2,M)
		 R3(M) + 5R1(M)]
	
	M = [R1(M)
		 (-2//9)R2(M)
		 R3(M)]
	
	M = [R(1:2,M)
		 R3(M) - (17//2)R2(M)]
	
	M = [R(1:2,M)
		 9R3(M)]
end

# ╔═╡ 2468ff5b-b17b-41b5-983d-46d4229c2dc8
md"""
### PS1 #2 (B)

$\begin{bmatrix}
2 & 3 & -3 & 1 \\
1 & -3 & 2 & 1 \\
-3 & -9 & 8 & -4
\end{bmatrix}$

$\begin{align*}
R_1 = \frac{1}{2} R_1 &\implies
\begin{bmatrix}
1 & \frac{3}{2} & -\frac{3}{2} & \frac{1}{2} \\
1 & -3 & 2 & 1 \\
-3 & -9 & 8 & -4
\end{bmatrix} \\

R_2 = R_2 - R_1 &\implies
\begin{bmatrix}
1 & \frac{3}{2} & -\frac{3}{2} & \frac{1}{2} \\
0 & -\frac{9}{2} & \frac{7}{2} & \frac{1}{2} \\
-3 & -9 & 8 & -4
\end{bmatrix} \\

R_3 = R_3 + 3R_1 &\implies
\begin{bmatrix}
1 & \frac{3}{2} & -\frac{3}{2} & \frac{1}{2} \\
0 & -\frac{9}{2} & \frac{7}{2} & \frac{1}{2} \\
0 & -\frac{9}{2} & \frac{7}{2} & -\frac{5}{2}
\end{bmatrix} \\

R_3 = R_3 - R_2 &\implies
\begin{bmatrix}
1 & \frac{3}{2} & -\frac{3}{2} & \frac{1}{2} \\
0 & -\frac{9}{2} & \frac{7}{2} & \frac{1}{2} \\
0 & 0 & 0 & -3
\end{bmatrix} \\
\end{align*}$
"""

# ╔═╡ db250a73-76e6-4288-945f-250e0a200df3
let
	A = [2 3 -3
		 1 -3 2
		-3 -9 8]
	
	b = [1
		 1
		-4]
	
	M = [A b]
	
	M = [(1//2)R1(M)
		R(2:3,M)]
	
	M = [R1(M)
		 R2(M) - R1(M)
		 R3(M)]
	
	M = [R(1:2,M)
		 R3(M) + 3R1(M)]
	
	M = [R(1:2,M)
		 R3(M) - R2(M)]
end

# ╔═╡ 965d09db-1fed-4706-9d7e-2b862a3ce44a
md"""
### PS1 #2 (C)

$\begin{bmatrix}
2 & 1 & -3 & 2 \\
1 & -3 & 2 & 1
\end{bmatrix}$

$\begin{align*}
R_1 = \frac{1}{2} R_1 &\implies
\begin{bmatrix}
1 & \frac{1}{2} & -\frac{3}{2} & 1 \\
1 & -3 & 2 & 1
\end{bmatrix} \\

R_2 = R_2 - R_1 &\implies
\begin{bmatrix}
1 & \frac{1}{2} & -\frac{3}{2} & 1 \\
0 & -\frac{7}{2} & \frac{7}{2} & 0
\end{bmatrix} \\

R_2 = -\frac{2}{7}R_2 &\implies
\begin{bmatrix}
1 & \frac{1}{2} & -\frac{3}{2} & 1 \\
0 & 1 & -1 & 0
\end{bmatrix}
\end{align*}$
"""

# ╔═╡ 7db67aff-3077-43d4-8177-6f5688b86eae
let
	A = [2 1 -3
		 1 -3 2]
	
	b = [2
		 1]
	
	M = [A b]
	
	M = [(1//2)R1(M)
		R2(M)]
	
	M = [R1(M)
		 R2(M) - R1(M)]
	
	M = [R1(M)
		 (-2//7)R2(M)]
end

# ╔═╡ 0df73e47-63fb-46ba-ab1f-68e6f504cab3
md"### PS1 #3"

# ╔═╡ 6b00b2de-35da-452b-953d-810f89e5fc3f
md"""
$(\text{A}) \qquad \begin{cases}
x + \frac{3}{2} y - \frac{3}{2} z = \frac{1}{2} \\
\qquad\;\; y - \frac{7}{9} z = -\frac{1}{9} \\
\qquad\qquad\quad z = 4
\end{cases}$

$\begin{align*}
y - \frac{7}{9} (4) &= -\frac{1}{9} \\
y - \frac{28}{9} &= -\frac{1}{9} \\
y &= 3
\end{align*}$

$\begin{align*}
x + \frac{3}{2} (3) - \frac{3}{2} (4) &= \frac{1}{2} \\
x + \frac{9}{2} - \frac{12}{2} &= \frac{1}{2} \\
x - \frac{3}{2} &= \frac{1}{2} \\
x &= 2
\end{align*}$

$\begin{bmatrix} x \\ y \\ z \end{bmatrix} =
\begin{bmatrix} 2 \\ 3 \\ 4 \end{bmatrix}$
"""

# ╔═╡ 6ac56c3d-8d4c-4a1c-a182-2c97e020ddc2
md"""

$(\text{B}) \qquad \begin{cases}
x + \frac{3}{2} y - \frac{3}{2} z = \frac{1}{2} \\
\;\; -\frac{9}{2}y + \frac{7}{9} z = \frac{1}{2} \\
\qquad\quad\quad\quad 0 = -3
\end{cases}$

$\text{No solution. } (0 \neq -3)$
"""

# ╔═╡ 4a056063-57d1-4596-b6c3-8892c6a78054
md"""
$(\text{C}) \qquad \begin{cases}
x + \frac{1}{2} y - \frac{3}{2} z = 1 \\
\qquad\;\; y - \;\;\; z = 0
\end{cases}$

$y = z$

$\begin{align*}
x + \frac{1}{2} y - \frac{3}{2} y &= 1 \\
-y &= 1 - x \\
y &= x - 1
\end{align*}$

$\begin{cases}
x = y + 1 \\
y = z
\end{cases}$

$\text{Infinitely many solutions.}$
"""

# ╔═╡ 0a904cfb-c4a9-48f8-8829-64a45fe362d9
md"""
## Problem Set 2

1. "Solve" (not really)

$\begin{align*}
x_1 - \frac{1}{2} x_2 - \frac{1}{4} x_4 &= -3/2 \\
-\frac{1}{4} x_1 + x_2 - \frac{1}{4} x_3 - \frac{1}{4} x_5 &= -1/2 \\
-\frac{1}{4} x_2 + x_3 - \frac{1}{4} x_6 &= \frac{5}{4} \\
-\frac{1}{4} x_1 + x_4 - \frac{1}{4} x_5 &= 0 \\
-\frac{1}{4} x_2 - \frac{1}{4} x_4 + x_5 - \frac{1}{4} x_6 &= \frac{1}{4} \\
-\frac{1}{4} x_3 - \frac{1}{4} x_5 + x_6 &= \frac{1}{4}
\end{align*}$

(a) Write down the augmented matrix of the system, call it $A$.

(b) After some work towards finding an echelon form for the augmented matrix, we reach

$\begin{bmatrix}
4 & -1 & 0 & -1 & 0 & 0 & -6 \\

0 & 15 & -4 & -1 & -4 & 0 & -14 \\

0 & 0 & 56 & -1 & -4 & -15 & 61 \\

0 & 0 & 0 & 209 & -60 & -1 & -93 \\

0 & 0 & 0 & 0 & 712 & -225 & -25 \\

0 & 0 & 0 & 0 & -225 & 780 & 435
\end{bmatrix}$

Use this to find $x_6$ as an exact fraction (as integer/integer).

*Note: You only need to work with the last three columns of the last two rows, which is the augmented matrix of a system of two linear equations in two unknowns.*
*Using row operations with this matrix with two rows and three columns, you will end up with some fractions, one of which should be $780 - 225^2 / 712$.*
*Add it as $(780 \times 712 - 225^2) / 712$ and use a calculator to compute the numerator, but not the quotient (e.g., if it were $1/52$, write this, not an approximate value like 0.0192308...).*
*Work your way through the problem like this, it's easier and cleaner.*
*The problem asks for an exact fraction.*
*Check your answer against the approximate value 0.602484.*

(c) In the process of solving (b) you will get an echelon form for $A$.
Use it to find the rank of the system of equations.

2. Let $B$ be the matrix

$\begin{bmatrix}
2 & 1 & 3 & -5 & 1 \\
-2 & 0 & -3 & 1 & 0 \\
1 & 0 & 3 & 4 & 2
\end{bmatrix}.$

(a) Find the reduced echelon form of $B$

(b) Solve each of the following two systems:

$\begin{equation}
\begin{aligned}[c]
2x_1 + x_2 + 3x_3 &= -5 \\
-2x_1 - 3x_3 &= 1 \\
x_1 + 3x_3 &= 4
\end{aligned}

\qquad\qquad

\begin{aligned}[c]
2x_1 + x_2 + 3x_3 &= 1 \\
-2x_1 - 3x_3 &= 0 \\
x_1 + 3x_3 &= 2
\end{aligned}
\end{equation}$

*Note: These are two separate systems of linear equations.*
*The last two columns of $B$ are the constant terms of the systems.*
*See pg. 32.*
"""

# ╔═╡ 9c0518db-b2cd-4499-82bc-3735826373e7
md"### PS2 #1"

# ╔═╡ f2a2e166-6257-4efc-be59-636099c68740
md"""
**(a)**

$\begin{bmatrix}
1 & -\frac{1}{4} & 0 & -\frac{1}{4} & 0 & 0 & -\frac{3}{2} \\
-\frac{1}{4} & 1 & -\frac{1}{4} & 0 & -\frac{1}{4} & 0 & -\frac{1}{2} \\
0 & -\frac{1}{4} & 1 & 0 & 0 & -\frac{1}{4} & \frac{5}{4} \\
-\frac{1}{4} & 0 & 0 & 1 & -\frac{1}{4} & 0 & 0 \\
0 & -\frac{1}{4} & 0 & -\frac{1}{4} & 1 & -\frac{1}{4} & \frac{1}{4} \\
0 & 0 & -\frac{1}{4} & 0 & -\frac{1}{4} & 1 & \frac{1}{4}
\end{bmatrix}$
"""

# ╔═╡ 9c86875e-6070-4601-b24b-3da1a0e79136
md"""
**(b)**

$\begin{bmatrix}
4 & -1 & 0 & -1 & 0 & 0 & -6 \\

0 & 15 & -4 & -1 & -4 & 0 & -14 \\

0 & 0 & 56 & -1 & -4 & -15 & 61 \\

0 & 0 & 0 & 209 & -60 & -1 & -93 \\

0 & 0 & 0 & 0 & 712 & -225 & -25 \\

0 & 0 & 0 & 0 & -225 & 780 & 435
\end{bmatrix}$

$\begin{align*}
R_5 \Longleftrightarrow R_6 &\Longrightarrow
\begin{bmatrix}
4 & -1 & 0 & -1 & 0 & 0 & -6 \\

0 & 15 & -4 & -1 & -4 & 0 & -14 \\

0 & 0 & 56 & -1 & -4 & -15 & 61 \\

0 & 0 & 0 & 209 & -60 & -1 & -93 \\

0 & 0 & 0 & 0 & -225 & 780 & 435 \\

0 & 0 & 0 & 0 & 712 & -225 & -25
\end{bmatrix} \\

R_5 = -\frac{1}{225} R_5 &\Longrightarrow
\begin{bmatrix}
4 & -1 & 0 & -1 & 0 & 0 & -6 \\

0 & 15 & -4 & -1 & -4 & 0 & -14 \\

0 & 0 & 56 & -1 & -4 & -15 & 61 \\

0 & 0 & 0 & 209 & -60 & -1 & -93 \\

0 & 0 & 0 & 0 & 1 & -\frac{52}{15} & -\frac{29}{15} \\

0 & 0 & 0 & 0 & 712 & -225 & -25
\end{bmatrix} \\

R_6 = R_6 - 712 R_5 &\Longrightarrow
\begin{bmatrix}
4 & -1 & 0 & -1 & 0 & 0 & -6 \\

0 & 15 & -4 & -1 & -4 & 0 & -14 \\

0 & 0 & 56 & -1 & -4 & -15 & 61 \\

0 & 0 & 0 & 209 & -60 & -1 & -93 \\

0 & 0 & 0 & 0 & 1 & -\frac{52}{15} & -\frac{29}{15} \\

0 & 0 & 0 & 0 & 0 & \frac{33649}{15} & \frac{20273}{15}
\end{bmatrix} \\

R_6 = \frac{15}{33649} R_6 &\Longrightarrow
\begin{bmatrix}
4 & -1 & 0 & -1 & 0 & 0 & -6 \\

0 & 15 & -4 & -1 & -4 & 0 & -14 \\

0 & 0 & 56 & -1 & -4 & -15 & 61 \\

0 & 0 & 0 & 209 & -60 & -1 & -93 \\

0 & 0 & 0 & 0 & 1 & -\frac{52}{15} & -\frac{29}{15} \\

0 & 0 & 0 & 0 & 0 & 1 & \frac{97}{161}
\end{bmatrix}

x_6 = \frac{97}{161} \approx 0.602
\end{align*}$
"""

# ╔═╡ 48b71ef8-d6ce-45ad-928a-2a2ca0e8ab5f
float(last(([712 -225 -25] - 712*[1 -52//15 -29//15]) * (15//33649)))

# ╔═╡ 8635836b-5e67-4f51-8501-881c358f08a4
md"""
**(b)**

The rank is 6.
"""

# ╔═╡ 3f1cf813-9ea9-41fa-b5a5-a13343f3b20c
md"### PS2 #2"

# ╔═╡ 9a1450a0-d356-42dd-806a-0bd6fc5e40de
md"""
**(a)**

$B = \begin{bmatrix}
2 & 1 & 3 & -5 & 1 \\
-2 & 0 & -3 & 1 & 0 \\
1 & 0 & 3 & 4 & 2
\end{bmatrix}$

$\begin{align*}
R_1 \Longleftrightarrow R_3 &\Longrightarrow
\begin{bmatrix}
1 & 0 & 3 & 4 & 2 \\
-2 & 0 & -3 & 1 & 0 \\
2 & 1 & 3 & -5 & 1
\end{bmatrix} \\

R_2 \Longleftrightarrow R_3 &\Longrightarrow
\begin{bmatrix}
1 & 0 & 3 & 4 & 2 \\
2 & 1 & 3 & -5 & 1 \\
-2 & 0 & -3 & 1 & 0
\end{bmatrix} \\

R_2 = R_2 - 2R_1 &\Longrightarrow
\begin{bmatrix}
1 & 0 & 3 & 4 & 2 \\
0 & 1 & -3 & -13 & -3 \\
-2 & 0 & -3 & 1 & 0
\end{bmatrix} \\

R_3 = R_3 + 2R_1 &\Longrightarrow
\begin{bmatrix}
1 & 0 & 3 & 4 & 2 \\
0 & 1 & -3 & -13 & -3 \\
0 & 0 & 3 & 9 & 4
\end{bmatrix} \\

R_3 = \frac{1}{3} R_3 &\Longrightarrow
\begin{bmatrix}
1 & 0 & 3 & 4 & 2 \\
0 & 1 & -3 & -13 & -3 \\
0 & 0 & 1 & 3 & \frac{4}{3}
\end{bmatrix} \\

R_1 = R_1 - 3R_3 &\Longrightarrow
\begin{bmatrix}
1 & 0 & 0 & -5 & -2 \\
0 & 1 & -3 & -13 & -3 \\
0 & 0 & 1 & 3 & \frac{4}{3}
\end{bmatrix} \\

R_2 = R_2 + 3R_3 &\Longrightarrow
\begin{bmatrix}
1 & 0 & 0 & -5 & -2 \\
0 & 1 & 0 & -4 & 1 \\
0 & 0 & 1 & 3 & \frac{4}{3}
\end{bmatrix}
\end{align*}$
"""

# ╔═╡ 12dfce5c-39ce-46e3-878e-324ec03f73af
let
	A = [2 1 3 -5
		-2 0 -3 1
		 1 0 3  4]
	
	b = [1
		 0
		 2]
	
	M = [A b]
	
	M = [R3(M)
		 R2(M)
		 R1(M)]
	
	M = [R1(M)
		 R3(M)
		 R2(M)]
	
	M = [R1(M)
		R2(M) - 2R1(M)
		R3(M)]
	
	M = [R(1:2,M)
		 R3(M) + 2R1(M)]
	
	M = [R(1:2,M)
		 (1/3)*R3(M)]
	
	M = [R1(M) - 3R3(M)
		 R(2:3,M)]
	
	M = [R1(M)
		 R2(M) + 3R3(M)
		 R3(M)]
	
	M, rref(M)
end

# ╔═╡ 05811194-b5d5-4cf8-af9d-37df311633cc
md"""
**(b)**

For the first system, the solution is

$\begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix}
= \begin{bmatrix} -5 \\ -4 \\ 3 \end{bmatrix}$


For the second system, the solution is

$\begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix}
= \begin{bmatrix} -2 \\ 1 \\ \frac{4}{3} \end{bmatrix}$
"""

# ╔═╡ 20e520b8-bcdb-4e3e-aa22-21042e6ae15d
# First system
let
	x1, x2, x3 = [-5 -4 3]
	[
		2x1 + x2 + 3x3 == -5
		-2x1 - 3x3 == 1
		x1 + 3x3 == 4
	]
end

# ╔═╡ d5bf5f56-3d4c-4afd-9c9d-b9f6a07ec145
# Second system
let
	x1, x2, x3 = [-2 1 4/3]
	[
		2x1 + x2 + 3x3 == 1
		-2x1 - 3x3 == 0
		x1 + 3x3 == 2
	]
end

# ╔═╡ 4f02f2d3-f72b-433f-ba23-6291ed773ddc
md"## Problem Set 3"

# ╔═╡ 59a6c778-3bc8-4462-abc7-1e5bdc1fd8c4
md"""
### PS3 #1

The matrix

$A = \begin{bmatrix}
1 & 2 & 0 \\
0 & 3 & 1 \\
1 & 2 & 2
\end{bmatrix}$

is the matrix of a certain system of equations (just the matrix of the system, not the augmented matrix).
Let $R_1$, $R_2$, $R_3$ be its rows.
Applying the sequence of elementary row operations

$\begin{align*}
{R_3}' &\gets R_3 - R_1 \\
{R_2}' &\gets R_2 - \frac{1}{2} R_3 \\
{R_1}' &\gets R_1 - \frac{2}{3} R_2
\end{align*}$

in sequence starting with $A$ yields

$\begin{bmatrix}
1 & 0 & 0 \\
0 & 3 & 0 \\
0 & 0 & 2
\end{bmatrix} \;.$

Let

$v_1 = 
\begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix}, \quad
v_2 = 
\begin{bmatrix} 2 \\ 3 \\ 2 \end{bmatrix}, \quad
v_3 =
\begin{bmatrix} 0 \\ 1 \\ 2 \end{bmatrix}, \quad
b = 
\begin{bmatrix} -1 \\ -2 \\ 1 \end{bmatrix}$

There are numbers $x$, $y$, $z$ such that

$xv_1 + yv_2 + zv_3 = b$

because $A$ has rank 3.
Find them.

To set up the problem, see/study Section 2.4 (pg. 71) of our textbook.
To solve it, think what happens to the last column in the augmented matrix of a system when you apply elementary row operations?
Following the row operations above leads to the next to last step in solving this problem.
"""

# ╔═╡ b3d9f2fe-afe0-4e85-9216-354386bbdc51
md"""
$\begin{align*}
R_2 = \frac{1}{3} R_2 &\implies v_2 = \begin{bmatrix} 2 \\ 1 \\ 2 \end{bmatrix}, \quad v_3 = \begin{bmatrix} 0 \\ \frac{1}{3} \\ \frac{2}{3} \end{bmatrix} \\
R_3 = \frac{1}{2} R_3 &\implies v_1 = \begin{bmatrix} 1 \\ 0 \\ \frac{1}{2} \end{bmatrix}, \quad v_2 = \begin{bmatrix} 2 \\ 1 \\ 1 \end{bmatrix}, \quad v_3 = \begin{bmatrix} 0 \\ \frac{1}{3} \\ \frac{1}{3} \end{bmatrix}
\end{align*}$
"""

# ╔═╡ cb72d5be-289d-4e5a-aad1-360c716026b3
md"""
$\begin{align*}
\begin{bmatrix} x \\ y \\ z \end{bmatrix} &= (-1) \begin{bmatrix} 1 \\ 0 \\ \frac{1}{2} \end{bmatrix} + (-2) \begin{bmatrix} 2 \\ 1 \\ 1 \end{bmatrix} + (1) \begin{bmatrix} 0 \\ 1 \\ 2 \end{bmatrix} \\
&= \begin{bmatrix} -1 \\ 0 \\ -\frac{1}{2} \end{bmatrix} + \begin{bmatrix} -4 \\ -2 \\ -2 \end{bmatrix} + \begin{bmatrix} 0 \\ 1 \\ 2 \end{bmatrix} \\
&= \begin{bmatrix} -5 \\ -1 \\ -\frac{1}{2} \end{bmatrix}
\end{align*}$
"""

# ╔═╡ 59df2f01-4912-4da7-a349-35447eac202b
let
	A = [1 2 0 1 0 0
		 0 3 1 0 1 0
		 1 2 2 0 0 1]
	
	A = [R(1:2,A)
		 R3(A) - R1(A)]
	
	A = [R1(A)
		 R2(A) - (1//2)R3(A)
		 R3(A)]
	
	A = [R1(A) - (2//3)R2(A)
		 R(2:3, A)]
	
	A = [R1(A)
		 (1//3)R2(A)
		 (1//2)R3(A)]
	
	b = [-1
		 -2
		  1]
	
	b[1] * A[:,4] + b[2] * A[:,5] + b[3] * A[:,6]
end

# ╔═╡ 1927374f-49d4-4d8c-97f3-e2ca80c945e3
md"## Problem Set 4"

# ╔═╡ bd563257-69ea-41ce-bcc6-3c1cc62a0d8d
md"""
### PS4 #1

(a) Find the rank of the system (*)

$\left\{\begin{align*}
2x_1 + 2x_2 - x_3 &= -1 \\
-x_1 + x_3 + x_4 &= -1 \\
-x_1 + 3x_2 + 3x_3 + x_4 &= 0
\end{align*}\right.$

(b) This system has exactly one solution with $x_4 = 1$.
    Find it, write it as a column vector, call it $X_0$.

(c) Let the column vector $X_h$ be the general solution of the associated homogeneous system

$\left\{\begin{align*}
2x_1 + 2x_2 - x_3 &= 0 \\
-x_1 + x_3 + x_4 &= 0 \\
-x_1 + 3x_2 + 3x_3 + x_4 &= 0
\end{align*}\right.$

(The entries of $X_h$ are the solutions $x_1$, $x_2$, $x_3$, $x_4$ of this system.)

(d) Let $X = X_0 + X_h$. Verify that $X$ solves (*) by directly replacing in the equations.
"""

# ╔═╡ e7dc3822-653a-4462-9e7a-1fde12cf638d
let
	rref([2 2 -1 0 -1
		 -1 0  1 1 -1
		 -1 3  3 1  0]),
	rref([2 2 -1 0 0
		 -1 0  1 1 0
		 -1 3  3 1 0])
end

# ╔═╡ f6bc8530-1822-44a7-b0ca-d39b1efd9225
let
	md"""
	(a) Represent it as an augmented matrix:

	$\begin{bmatrix}
	2 & 2 & -1 & 0 & -1 \\
	-1 & 0 & 1 & 1 & -1 \\
	-1 & 3 & 3 & 1 & 0
	\end{bmatrix}$
	
	Reduced row echelon form:
	
	$\begin{align*}
	R_2 \Longleftrightarrow R_3 &\implies \begin{bmatrix}
	2 & 2 & -1 & 0 & -1 \\
	-1 & 3 & 3 & 1 & 0 \\
	-1 & 0 & 1 & 1 & -1
	\end{bmatrix} \\
	R_1 \gets \frac{1}{2} R_1 &\implies \begin{bmatrix}
	1 & 1 & -\frac{1}{2} & 0 & -\frac{1}{2} \\
	-1 & 3 & 3 & 1 & 0 \\
	-1 & 0 & 1 & 1 & -1
	\end{bmatrix} \\
	R_2 \gets R_2 + R_1 &\implies \begin{bmatrix}
	1 & 1 & -\frac{1}{2} & 0 & -\frac{1}{2} \\
	0 & 4 & \frac{5}{2} & 1 & -\frac{1}{2} \\
	-1 & 0 & 1 & 1 & -1
	\end{bmatrix} \\
	R_3 \gets R_3 + R_1 &\implies \begin{bmatrix}
	1 & 1 & -\frac{1}{2} & 0 & -\frac{1}{2} \\
	0 & 4 & \frac{5}{2} & 1 & -\frac{1}{2} \\
	0 & 1 & \frac{1}{2} & 1 & -\frac{3}{2}
	\end{bmatrix} \\
	R_2 \gets \frac{1}{4} R_2 &\implies \begin{bmatrix}
	1 & 1 & -\frac{1}{2} & 0 & -\frac{1}{2} \\
	0 & 1 & \frac{5}{8} & \frac{1}{4} & -\frac{1}{8} \\
	0 & 1 & \frac{1}{2} & 1 & -\frac{3}{2}
	\end{bmatrix} \\
	R_3 \gets R_3 - R_2 &\implies \begin{bmatrix}
	1 & 1 & -\frac{1}{2} & 0 & -\frac{1}{2} \\
	0 & 1 & \frac{5}{8} & \frac{1}{4} & -\frac{1}{8} \\
	0 & 0 & -\frac{1}{8} & \frac{3}{4} & -\frac{11}{8}
	\end{bmatrix} \\
	&\implies \begin{bmatrix}
	1 & 0 & 0 & -7 & 12 \\
	0 & 1 & 0 & 4 & -7 \\
	0 & 0 & 1 & -6 & 11
	\end{bmatrix}\end{align*}$
	
	The rank of the system is 3.
	"""
end

# ╔═╡ e42f0f08-fbbe-4547-a970-cc4204f447b1
md"""
(b) 

$x_1 = 12 + 7t = 19$
$x_2 = -7 - 4t = -11$
$x_3 = 11 + 6t = 17$
$x_4 = t$

$X_0 = \begin{bmatrix} 19 \\ -11 \\ 17 \\ 1 \end{bmatrix}$
"""

# ╔═╡ 9f527f88-244a-4f60-8c45-fcf37b760656
md"""
(c)

$\begin{bmatrix}
2 & 2 & -1 & 0 & 0 \\
-1 & 0 & 1 & 1 & 0 \\
-1 & 3 & 3 & 1 & 0
\end{bmatrix}$

$\begin{align*}
R_1 \Longleftrightarrow R_2 &\implies \begin{bmatrix}
-1 & 0 & 1 & 1 & 0 \\
2 & 2 & -1 & 0 & 0 \\
-1 & 3 & 3 & 1 & 0
\end{bmatrix} \\
R_1 \gets -R_1 &\implies \begin{bmatrix}
1 & 0 & -1 & -1 & 0 \\
2 & 2 & -1 & 0 & 0 \\
-1 & 3 & 3 & 1 & 0
\end{bmatrix} \\
R_2 \gets R_2 - 2R_1 &\implies \begin{bmatrix}
1 & 0 & -1 & -1 & 0 \\
0 & 2 & 1 & 2 & 0 \\
-1 & 3 & 3 & 1 & 0
\end{bmatrix} \\
R_3 \gets R_3 + R_1 &\implies \begin{bmatrix}
1 & 0 & -1 & -1 & 0 \\
0 & 2 & 1 & 2 & 0 \\
0 & 3 & 2 & 0 & 0
\end{bmatrix} \\
R_2 \gets \frac{1}{2} R_2 &\implies \begin{bmatrix}
1 & 0 & -1 & -1 & 0 \\
0 & 1 & \frac{1}{2} & 1 & 0 \\
0 & 3 & 2 & 0 & 0
\end{bmatrix} \\
R_3 \gets R_3 - 3R_2 &\implies \begin{bmatrix}
1 & 0 & -1 & -1 & 0 \\
0 & 1 & \frac{1}{2} & 1 & 0 \\
0 & 0 & \frac{1}{2} & -3 & 0
\end{bmatrix} \\
R_3 \gets 2R_3 &\implies \begin{bmatrix}
1 & 0 & -1 & -1 & 0 \\
0 & 1 & \frac{1}{2} & 1 & 0 \\
0 & 0 & 1 & -6 & 0
\end{bmatrix} \\
R_2 \gets R_2 - \frac{1}{2}R_3 &\implies \begin{bmatrix}
1 & 0 & -1 & -1 & 0 \\
0 & 1 & 0 & 4 & 0 \\
0 & 0 & 1 & -6 & 0
\end{bmatrix} \\
R_1 \gets R_1 + R_3 &\implies \begin{bmatrix}
1 & 0 & 0 & -7 & 0 \\
0 & 1 & 0 & 4 & 0 \\
0 & 0 & 1 & -6 & 0
\end{bmatrix} \\
\end{align*}$

$X_h = t \begin{bmatrix} 7 \\ -4 \\ 6 \\ 1 \end{bmatrix}$
"""

# ╔═╡ 8b056c30-1e7c-4882-9a0c-f23427509f97
md"""
$X = X_0 + X_h$

$X = \begin{bmatrix} 19 \\ -11 \\ 17 \\ 1 \end{bmatrix} + t \begin{bmatrix} 7 \\ -4 \\ 6 \\ 1 \end{bmatrix} = \begin{bmatrix} 19 + 7t \\ -11 - 4t \\ 17 + 6t \\ 1 + t
\end{bmatrix}$

First equation:

$\begin{align*}
2(19 + 7t) + 2(-11 - 4t) - (17 + 6t) &= -1 \\
38 + 14t - 22 - 8t - 17 - 6t &= -1 \\
0 - 1 &= -1 \\
-1 &= -1
\end{align*}$

Second equation:

$\begin{align*}
-(19 + 7t) + (17 + 6t) + (1 + t) &= -1 \\
-19 - 7t + 17 + 6t + 1 + t &= -1 \\
0 - 1 &= -1 \\
-1 &= -1
\end{align*}$

Third equation:

$\begin{align*}
-(19 + 7t) + 3(-11 - 4t) + 3(17 + 6t) + (1 + t) &= 0 \\
-19 - 7t - 33 - 12t + 51 + 18t + 1 + t &= 0 \\
0 &= 0
\end{align*}$
"""

# ╔═╡ 0ebc87b1-5187-4106-9c95-6cfb75fc97a7
md"""
### PS4 #2

(a) Find the rank of the system

$\left\{\begin{align*}
x_1 + 4x_2 + 2x_3 + 4x_4 &= 0 \\
4x_1 + 4x_2 - x_3 - x_4 &= 0
\end{align*}\right.$

(b) Find all solutions of this system.
"""

# ╔═╡ 07d9902a-af9f-474a-a015-56a7ce5b2476
rref([1 4 2 4 0
	  4 4 -1 -1 0])

# ╔═╡ 4b780acb-6fe8-446e-a0d8-293f520dd98c
md"""
(a) Represent the system as an augmented matrix:

$\begin{bmatrix}
1 & 4 & 2 & 4 & 0 \\
4 & 4 & -1 & -1 & 0
\end{bmatrix}$

Reduce to row echelon form.

$\begin{align*}
R_2 \gets R_2 - 4R_1 &\implies \begin{bmatrix}
1 & 4 & 2 & 4 & 0 \\
0 & -12 & -9 & -17 & 0
\end{bmatrix} \\
R_2 \gets -\frac{1}{12}R_2 &\implies \begin{bmatrix}
1 & 4 & 2 & 4 & 0 \\
0 & 1 & \frac{3}{4} & \frac{17}{12} & 0
\end{bmatrix} \\
\end{align*}$

The system's rank is 2.
"""

# ╔═╡ 3c3a4efd-21cc-4d78-b86b-da0118a9eaff
md"""
(b) Reduced row echelon form:

$\begin{align*}
R_1 \gets R_1 - 4R_2 &\implies \begin{bmatrix}
1 & 0 & -1 & -\frac{5}{3} & 0 \\
0 & 1 & \frac{3}{4} & \frac{17}{12} & 0
\end{bmatrix}
\end{align*}$

Solution:

$\begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} \frac{5}{3} x_4 + x_3 \\ -\frac{17}{12} x_4 - \frac{3}{4} x_3 \end{bmatrix}$
"""

# ╔═╡ 97ca197b-7d1f-4675-acd7-cd622946a24e
md"""
### PS4 #3

Let

$\vec{v}_1 = \begin{bmatrix} 1 \\ -1 \\ 7 \\ -6 \end{bmatrix}, \qquad \vec{v}_2 = \begin{bmatrix} 0 \\ 2 \\ 11 \\ -10 \end{bmatrix}$

(a) Verify that every linear combination of these two vectors is a solution of

$\left\{\begin{align*}
3x_1 + x_2 - 2x_3 - 2x_4 &= 0 \\
5x_1 + 3x_2 + 4x_3 + 5x_4 &= 0
\end{align*}\right.$

(b) Find the rank of this system and its general solution.
"""

# ╔═╡ 53fdf8ed-9a96-41d8-81a8-7c5fa91bb00c
rref([3 1 -2 -2 0
	  5 3  4  5 0])

# ╔═╡ 12d3e4b5-b0df-4b90-baf7-3ec9a319d153
md"""
Augmented matrix:

$\begin{bmatrix}
3 & 1 & -2 & -2 & 0 \\
5 & 3 & 4 & 5 & 0
\end{bmatrix}$

Reduced row echelon form:

$\begin{align*}
R_1 \gets \frac{1}{3} R_1 &\implies \begin{bmatrix}
1 & \frac{1}{3} & -\frac{2}{3} & -\frac{2}{3} & 0 \\
5 & 3 & 4 & 5 & 0
\end{bmatrix} \\
R_2 \gets R_2 - 5R_1 &\implies \begin{bmatrix}
1 & \frac{1}{3} & -\frac{2}{3} & -\frac{2}{3} & 0 \\
0 & \frac{4}{3} & \frac{22}{3} & \frac{25}{3} & 0
\end{bmatrix} \\
R_2 \gets \frac{3}{4} R_2 &\implies \begin{bmatrix}
1 & \frac{1}{3} & -\frac{2}{3} & -\frac{2}{3} & 0 \\
0 & 1 & \frac{11}{2} & \frac{25}{4} & 0
\end{bmatrix} \\
R_1 \gets R_1 - \frac{1}{3} R_2 &\implies \begin{bmatrix}
1 & 0 & -\frac{5}{2} & -\frac{11}{4} & 0 \\
0 & 1 & \frac{11}{2} & \frac{25}{4} & 0
\end{bmatrix} \\
\end{align*}$

The rank of the system is 2.

Solution:

$\begin{bmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{bmatrix} = \begin{bmatrix} \frac{11}{4} x_4 + \frac{5}{2} x_3 \\ -\frac{25}{4} x_4 - \frac{11}{2} x_3 \\ x_3 \\ x_4 \end{bmatrix}$
"""

# ╔═╡ 133cf267-aa8d-46a6-a1f7-4d4a833fad75
md"""
For $x_3 = 7$ and $x_4 = -6$:

$\begin{bmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{bmatrix} = \begin{bmatrix} 1 \\ -1 \\ 7 \\ -6 \end{bmatrix}$

For $x_3 = 11$ and $x_4 = -10$:

$\begin{bmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{bmatrix} = \begin{bmatrix} 0 \\ 2 \\ 11 \\ -10 \end{bmatrix}$
"""

# ╔═╡ 50a6df63-9753-46fc-bb20-d6b489108197
md"""
### PS4 #4

$\vec{a} = \begin{bmatrix} 1 \\ -2 \\ 3 \end{bmatrix}, \qquad \vec{b} = \begin{bmatrix} 0 \\ 2 \\ -2 \end{bmatrix}$

Find all vectors $\vec{x} = \begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix}$ such that

$\vec{x} ⋅ \vec{a} = 0, \qquad \vec{x} ⋅ \vec{b} = 0.$
"""

# ╔═╡ a1409499-b628-4506-a5de-c3947f965433
md"""
$\left\{\begin{align*}
x_1 - 2x_2 + 3x_3 &= 0 \\
2x_2 - 2x_3 &= 0
\end{align*}\right.$

Augmented matrix:

$\begin{bmatrix}
1 & -2 & 3 & 0 \\
0 & 2 & -2 & 0
\end{bmatrix}$

Reduced row echelon form:

$\begin{align*}
R_2 \gets \frac{1}{2} R_2 &\implies \begin{bmatrix}
1 & -2 & 3 & 0 \\
0 & 1 & -1 & 0 
\end{bmatrix} \\
R_1 \gets R_1 + 2R_2 &\implies \begin{bmatrix}
1 & 0 & 1 & 0 \\
0 & 1 & -1 & 0 
\end{bmatrix} \\
\end{align*}$

Solution:

$\begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} = \begin{bmatrix} -x_3 \\ x_3 \\ x_3 \end{bmatrix} = \begin{bmatrix} -1 \\ 1 \\ 1 \end{bmatrix} x_3$
"""

# ╔═╡ 5bd95da9-b19b-4502-8f86-e81ec8af4acf
md"""
### PS4 #5

Let

$\vec{a} = \begin{bmatrix} 3 \sqrt{3} \\ 6 \\ 3 \end{bmatrix}, \qquad \begin{bmatrix} -2\sqrt{3} \\ 4 \\ -2 \end{bmatrix}$

Find: (a) $\|\vec{a}\|^2$; (b) $\|\vec{b}\|^2$; (c) $\|\vec{a} + \vec{b}\|^2$; (d) $\vec{a} ⋅ \vec{b}$.

What is $\|\vec{a}\|^2 + \|\vec{b}\|^2$?
"""

# ╔═╡ 904607da-67c6-4f45-a2bb-2e8d58514830
md"""
$\|\vec{a}\|^2 = (3\sqrt{3})^2 + 6^2 + 3^2 = 27 + 36 + 9 = 72$

$\|\vec{b}\|^2 = (-2\sqrt{3})^2 + 4^2 + (-2)^2 = 12 + 16 + 4 = 32$

$\vec{a} + \vec{b} = \begin{bmatrix} \sqrt{3} \\ 10 \\ 1 \end{bmatrix}$

$\|\vec{a} + \vec{b}\|^2 = 3 + 100 + 1 = 104$

$\vec{a} ⋅ \vec{b} = (3\sqrt{3})(-2\sqrt{3}) + (6)(4) + (3)(-2) = -18 + 24 - 6 = 0$

$\|\vec{a}\|^2 + \|\vec{b}\|^2 = 72 + 32 = 104$
"""

# ╔═╡ 509a9d05-fb66-49a8-945a-9b020627f745
md"## Problem Set #5"

# ╔═╡ eb31d37b-ac26-4d04-8598-097dd16b332e
md"""
### PS5 #1

Let $\theta$ be some number, let $\vec{v} = \begin{bmatrix} \cos{\theta} \\ \sin{\theta} \end{bmatrix}$.

(a) Find $\|\vec{v}\|$.

(b) Find a vector $\vec{w}$ orthogonal to $\|\vec{v}\|$ and of the same norm (magnitude) as $\vec{v}$. (How many such vectors are there?)
"""

# ╔═╡ a25c703f-537b-4cfb-ad60-559d027669a7
md"""
**(a)**

$\|\vec{v}\| = \cos^2{\theta} + \sin^2{\theta} = 1$
"""

# ╔═╡ 080dd720-d75f-4108-aa22-cfe946675b41
md"""
**(b)**

$\vec{v} ⋅ \vec{w} = \vec{0}$
"""

# ╔═╡ e9b8568e-a895-48fc-aa32-fabb03e574f2
md"""
### PS5 #2

Each of the vectors

$\vec{a} = \begin{bmatrix} -4 \\ -5 \\ 2 \\ 1 \end{bmatrix}, \qquad \vec{b} = \begin{bmatrix} -1 \\ -\frac{3}{2} \\ 0 \\ 1 \end{bmatrix}, \qquad \vec{c} = \begin{bmatrix} -\frac{5}{2} \\ -\frac{13}{4} \\ 1 \\ 1 \end{bmatrix}$

can be written as a linear combination of the other two.
Find them.
"""

# ╔═╡ 01f57e22-3edb-4a13-8d8d-24ebc30d4878
md"""
### PS5 #3

Let $\vec{a}$ be the solution of the system

$\begin{cases}\begin{align*}
-x_1 + 2x_2 + 2x_3 + 2x_4 = 0 \\
2x_1 + 3x_2 + 2x_4 = 0
\end{align*}\end{cases}$

with $x_3 = 1$ and $x_4 = 0$, let $\vec{b}$ be the solution with $x_3 = 0$ and $x_4 = 1$.
Find two non-collinear vectors $\vec{u}$ and $\vec{v}$ which are orthogonal to all linear combinations of $\vec{a}$ and $\vec{b}$.
"""

# ╔═╡ 4cb458b2-67a6-433c-8286-4b276f52a183
md"""
### PS5 #4

(a) Find a nonzero vector $\vec{v}$ which is orthogonal to every solution of

$2x + 3y + 4z = 1\,.$

Then find a unit vector $\vec{u}$ which is orthogonal to every solution of the equation.
How many such unit vectors are there?

(b) Are the vectors $\vec{v}$ of the first part also orthogonal to all solutions of

$2x + 3y + 4z = 2?$
"""

# ╔═╡ 6487ad98-683c-429b-a0d2-a370cc0b0c83
md"""
**(a)** Let $x = s$, $y = t$.

$4z = 1 - 2s - 3t$

$\begin{bmatrix} s \\ t \\ \frac{1}{4} (1 - 2s - 3t)\end{bmatrix} = \begin{bmatrix} x(s,t) \\ y(s,t) \\ z(s,t) \end{bmatrix} = \vec{s}$

Look for $\begin{bmatrix} c_1 \\ c_2 \\ c_3 \end{bmatrix}$ s.t. $\begin{bmatrix} c_1 \\ c_2 \\ c_3 \end{bmatrix} ⋅ \vec{s} = 0$

$\begin{align*}
\begin{bmatrix} c_1 \\ c_2 \\ c_3 \end{bmatrix} ⋅ \begin{bmatrix} s \\ t \\ \frac{1}{4}(1 - 2s - 3t) \end{bmatrix} &= c_1 s + c_2 t + c_3 \frac{1}{4} (1 - 2s - 3t) \\
&= \left(c_1 - \frac{c_3}{2}\right)s + \left(c_2 - \frac{3}{4} c_3\right) t + \frac{c_3}{4}
\end{align*}$

$c_1 - \frac{c_3}{2} = 0 \implies c_1 = \frac{1}{2} c_3$

$c_2 - \frac{3}{4} c_3 = 0 \implies c_2 = \frac{3}{4} c_3$

$\begin{bmatrix} \frac{1}{2} c_3 \\ \frac{3}{4} c_3 \\ c_3 \end{bmatrix}$
"""

# ╔═╡ 7051900e-7863-45cd-b1dd-7066793b8360
md"""
### PS5 #5

Let $\vec{v} = \begin{bmatrix} \cos{\theta} \\ \sin{\theta} \end{bmatrix}$ (as in Problem 1).
Let $\vec{w} = \begin{bmatrix} 1 \\ 3 \end{bmatrix}$

(a) Find

$\vec{w} - \text{proj}_{\vec{v}}{\vec{w}}\,.$

(b) Evaluate

$(\vec{w} - \text{proj}_{\vec{v}}{\vec{w}}) ⋅ \vec{v}\,.$
"""

# ╔═╡ 39204db3-380a-443f-b974-8bc652118ddf
md"""
### PS5 #6

Let $\vec{u}$ be a nonzero vector.
Let $\vec{w}$ be arbitrary.
Simplify

$\|\text{proj}_{\vec{u}}{\vec{w}}\|^2 + \|w - \text{proj}_{\vec{u}}{\vec{w}}\|^2$
"""

# ╔═╡ f976f042-3cc6-4f50-88c4-7dcd488d7cd5
md"""
### PS5 #7

The angle between certain two nonzero vectors $\vec{v}$ and $\vec{w}$ is $5 \pi / 6$.

(a) Find the angle between $3\vec{v}$ and $7\vec{w}$.

(b) Find the angle between $\vec{v}$ and $-\vec{w}$.
"""

# ╔═╡ 2bc90945-8e1f-47ca-a519-f8d93ddacbcf
md"""
### PS5 #8

Find the angle between the given vectors:

(1) $⟨3,0⟩$, $\left\langle\frac{3}{2}, \frac{3\sqrt{3}}{2}\right\rangle$

(2) $\left\langle2, 2\sqrt{3}\right\rangle$, $\left\langle-\sqrt{3},1\right\rangle$

(3) $⟨0,-3⟩$, $\left\langle 2, -2\sqrt{3} \right\rangle$

(4) $⟨-1,0⟩$, $\left\langle \frac{5}{2}, -\frac{1}{2} \left(5\sqrt{3}\right)\right\rangle$

(5) $\left\langle -2\sqrt{3}, 2 \right\rangle, \left\langle\frac{3 \left(\sqrt{3}-1\right)}{2 \sqrt{2}},-\frac{\left(\sqrt{3}+1\right)}{2\sqrt{2}}\right\rangle$

(6) $\left\langle\frac{3}{\sqrt{2}}, \frac{3}{\sqrt{2}}\right\rangle$, $\left\langle-\frac{1}{\sqrt{2}},-\frac{1}{\sqrt{2}}\right\rangle$
"""

# ╔═╡ 26dbb5a5-0005-40ea-bb88-963b4b995520
md"## Problem Set #6"

# ╔═╡ 7888fa3b-3074-4ccd-8f64-f7022e32aec7
md"""
### PS6 #1

Let

$P = \begin{bmatrix} 1 & 3 & 5 \\ 1 & -2 & 1 \end{bmatrix}, \qquad Q = \begin{bmatrix} -3 & 1 & 2 \\ -5 & -1 & 4 \\ 2 & 9 & 0 \end{bmatrix}\;.$

Find $PQ$.
Does $QP$ make sense?.
"""

# ╔═╡ b5c21a12-d081-4686-ba67-6cd9256cc09e
let
	P = [1 3 5
		 1 -2 1]
	
	Q = [-3 1 2
		 -5 -1 4
		  2 9 0]
	
	R1 = P[1,:]
	R2 = P[2,:]
	
	C1 = Q[:,1]
	C2 = Q[:,2]
	C3 = Q[:,3]
	
	[
		(R1 ⋅ C1) (R1 ⋅ C2) (R1 ⋅ C3)
		(R2 ⋅ C1) (R2 ⋅ C2) (R2 ⋅ C3)
	]
end

# ╔═╡ 5135f07a-a843-4751-9c54-35712ab2ac04
md"``QP`` does not make sense since you cannot multiply a 3x3 matrix by a 2x3 matrix."

# ╔═╡ 57bd2d47-b14c-4459-9c10-7d8d71705254
md"""
### PS6 #2

Let $A$ be the matrix formed by the first two columns of

$E = \begin{bmatrix} 3 & 1 & 1 & 0 \\ -1 & 2 & 0 & 1 \end{bmatrix}$

Bring $E$ to reduced echelon form (call it $E'$), then let $B$ be the matrix formed by the last two columns of $E'$.
Find $AB$ and $BA$
"""

# ╔═╡ 4c41085c-3443-4716-bba2-7839ea337c40
let
	E = [3 1 1 0
		-1 2 0 1]
	
	A = E[:, 1:2]
	
	E′ = rref(E)
	
	B = E′[:, 3:4]
	
	A*B, B*A
end

# ╔═╡ ea7c75d2-eb5e-4d51-8ee9-dee330304b6b
md"""
### PS6 #3

Let

$A = \begin{bmatrix} 2 & 1 & 2 \\ 2 & -3 & -1 \\ 0 & -3 & 0 \end{bmatrix},
\qquad I = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix}\;.$

Let $C$ be the $3 × 6$ matrix whose first three columns are the columns of $A$ (in the same order) and whose last three columns are those of $I$ (also in order).
Use Gauss-Jordan elimination to transform $C$ to reduced echelon form, call that $C'$.
Let $B$ be the $3 × 3$ matrix whose columns are the last three columns of $C'$.
Verify that $BA = I$ and $AB = I$.
"""

# ╔═╡ fb5b591a-f32c-4bfe-a2cb-8183ba942a65
let
	A = [2 1 2
		 2 -3 -1
	     0 -3 0]
	
	C = [A I(3)]
	
	C′ = rref(C)
	
	B = C′[:,4:end]
	
	A*B ≈ I(3), B*A ≈ I(3)
end

# ╔═╡ 8677382d-34c3-4de1-8894-6028380d5b02
md"""
### PS6 #4

Let $A$, $B$, and $I$ be the matrices in problem 3.
Let $b = \begin{bmatrix} -1 \\ 2 \\ 3 \end{bmatrix}$.
Let $x = Bb$.
Find $Ax$.
"""

# ╔═╡ 698abde3-1f79-4b28-9ffa-224972ba388a
let	
	A = [2 1 2
		 2 -3 -1
	     0 -3 0]
	
	C = [A I(3)]
	
	C′ = rref(C)
	
	B = C′[:,4:end]
	
	b = [-1
		 2
		 3]
	
	x = B * b
	
	A * x
end

# ╔═╡ bd10b677-7b37-41ef-a4a0-3ad1017e5acd
md"""
### PS6 #5

Let

$A = \begin{bmatrix}
a_{11} & a_{12} & a_{13} \\
a_{21} & a_{22} & a_{23} \\
a_{31} & a_{32} & a_{33}
\end{bmatrix}, \qquad E = \begin{bmatrix}
1 & 0 & 0 \\
0 & 0 & 1 \\
0 & 1 & 0
\end{bmatrix}\;.$

Find $EA$.
"""

# ╔═╡ 88d2cf92-46fa-4f9b-aab9-9a2fb8111d1c
md"""
$\begin{align*}
EA &= \begin{bmatrix}
1 & 0 & 0 \\
0 & 0 & 1 \\
0 & 1 & 0
\end{bmatrix} \begin{bmatrix}
a_{11} & a_{12} & a_{13} \\
a_{21} & a_{22} & a_{23} \\
a_{31} & a_{32} & a_{33}
\end{bmatrix} \\
&= \begin{bmatrix}
\left(\begin{bmatrix} a_{11} \\ a_{21} \\ a_{31} \end{bmatrix} ⋅ \begin{bmatrix}1 \\ 0 \\ 0\end{bmatrix}\right) & \left(\begin{bmatrix} a_{12} \\ a_{22} \\ a_{32} \end{bmatrix} ⋅ \begin{bmatrix}1 \\ 0 \\ 0\end{bmatrix}\right) & \left(\begin{bmatrix} a_{13} \\ a_{23} \\ a_{33} \end{bmatrix} ⋅ \begin{bmatrix}1 \\ 0 \\ 0\end{bmatrix}\right) \\
\left(\begin{bmatrix} a_{11} \\ a_{21} \\ a_{31} \end{bmatrix} ⋅ \begin{bmatrix}0 \\ 0 \\ 1\end{bmatrix}\right) & \left(\begin{bmatrix} a_{12} \\ a_{22} \\ a_{32} \end{bmatrix} ⋅ \begin{bmatrix}0 \\ 0 \\ 1\end{bmatrix}\right) & \left(\begin{bmatrix} a_{13} \\ a_{23} \\ a_{33} \end{bmatrix} ⋅ \begin{bmatrix}0 \\ 0 \\ 1\end{bmatrix}\right) \\
\left(\begin{bmatrix} a_{11} \\ a_{21} \\ a_{31} \end{bmatrix} ⋅ \begin{bmatrix}0 \\ 1 \\ 0\end{bmatrix}\right) & \left(\begin{bmatrix} a_{12} \\ a_{22} \\ a_{32} \end{bmatrix} ⋅ \begin{bmatrix}0 \\ 1 \\ 0\end{bmatrix}\right) & \left(\begin{bmatrix} a_{13} \\ a_{23} \\ a_{33} \end{bmatrix} ⋅ \begin{bmatrix}0 \\ 1 \\ 0\end{bmatrix}\right) 
\end{bmatrix} \\
&= \begin{bmatrix}
a_{11} & a_{12} & a_{13} \\
a_{31} & a_{32} & a_{33} \\
a_{21} & a_{22} & a_{23}
\end{bmatrix}
\end{align*}$
"""

# ╔═╡ 329ecc56-67a8-4ffb-8044-7dfca93a8116
md"## Problem Set #7"

# ╔═╡ a29e9713-6550-43ee-95eb-eb35d46e51cc
md"""
### PS7 #1

The task is to solve $(^*)$

$\begin{cases}\begin{align*}
-x_1 + 2x_2 + x_3 &= -1 \\
3x_1 - x_2 - x_3 &= 0 \\
-x_1 + 3x_2 + 2x_3 &= -1
\end{align*}\end{cases}$

Applying Gauss-Jordan elimination to

$\begin{bmatrix}
-1 & 2 & 1 & 1 & 0 & 0 \\
3 & -1 & -1 & 0 & 1 & 0 \\
-1 & 3 & 2 & 0 & 0 & 1
\end{bmatrix}$

leads to

$\begin{bmatrix}
1 & 0 & 0 & -\frac{1}{3} & \frac{1}{3} & \frac{1}{3} \\
0 & 1 & 0 & \frac{5}{3} & \frac{1}{3} & -\frac{2}{3} \\
0 & 0 & 1 & -\frac{8}{3} & -\frac{1}{3} & \frac{5}{3}
\end{bmatrix}$

Use this information to solve $(^*)$.
"""

# ╔═╡ 2415a79f-b72c-41d4-9784-0dec68ee4423
md"""
$Ax = b \implies x = A^{-1} b$

$A^{-1} = 
\begin{bmatrix}
-\frac{1}{3} & \frac{1}{3} & \frac{1}{3} \\
\frac{5}{3} & \frac{1}{3} & -\frac{2}{3} \\
-\frac{8}{3} & -\frac{1}{3} & \frac{5}{3}
\end{bmatrix}$

$b = \begin{bmatrix} -1 \\ 0 \\ -1 \end{bmatrix}$

$\begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} = \begin{bmatrix} 0 \\ -1 \\ 1 \end{bmatrix}$
"""

# ╔═╡ f25e481c-3639-4980-947a-83128a148a81
let
	A = [-1 2 1
		3 -1 -1
		-1 3 2]
	
	B = [-1//3 1//3 1//3
		5//3 1//3 -2//3
		-8//3 -1//3 5//3]
	
	b = [-1; 0; -1]
	
	float.(B)
	
	float.(B) * b
	
	inv(A) * b
end

# ╔═╡ cf5ac8d4-d741-4946-a7a6-1a3bed5304d8
md"""
### PS7 #2

Find the general solution of $(^*)$

$\begin{cases}\begin{align*}
4x_1 + 4x_2 + 4x_3 + x_4 &= 1 \\
3x_1 + 2x_2 + 4x_3 + x_4 &= 1 \\
-x_1 + 2x_2 + 3x_3 + 2x_4 &= 1 \\
\end{align*}\end{cases}$

$A = \begin{bmatrix}
4 & 4 & 4 & 1 & 1 \\
3 & 2 & 4 & 1 & -2 \\
-1 & 2 & 3 & 2 & 2
\end{bmatrix}$

$A' = \begin{bmatrix}
1 & 0 & 0 & -\frac{5}{14} & -\frac{11}{14} \\
0 & 1 & 0 & \frac{5}{28} & \frac{53}{28} \\
0 & 0 & 1 & \frac{3}{7} & -\frac{6}{7}
\end{bmatrix}$

$\begin{align*}
x_1 - \frac{5}{14} x_4 &= -\frac{11}{14} \\
x_2 + \frac{5}{28} x_4 &= \frac{53}{28} \\
x_3 + \frac{3}{7} x_4 &= -\frac{6}{7}
\end{align*}$

$\begin{bmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{bmatrix} = \begin{bmatrix} -\frac{11}{14} \\ \frac{53}{28} \\ -\frac{6}{7} \\ 0 \end{bmatrix} + \begin{bmatrix} \frac{5}{14} \\ -\frac{5}{28} \\ -\frac{3}{7} \\ 1 \end{bmatrix} x_4$
"""

# ╔═╡ ad6a34d8-9bd1-4105-8afd-f19e5328933a
md"""
### PS7 #3

Find the rank of the system
"""

# ╔═╡ 29860b91-212e-4b88-b199-af479a132c01
let
	A = [-6  2 -8 -4 -2
		  3 -1  5  2  2
		  9 -3  9  7 -1
		  3 -1  4  1 -3]
	
	ceil.(Int, rref(A[1:4,1:4])), ceil.(Int, rref(A))
end

# ╔═╡ b0927f96-ec7b-4a96-8a10-d69b55ec2c39
md"r = 4"

# ╔═╡ 590ba8f2-327c-43cf-8ea6-79bd2bccbbfe
md"""
### PS7 #4
"""

# ╔═╡ b5de1268-62a9-46c0-84b4-fb91160d729f
let
	A1 = [ 1 -2
		  -4  3]
	
	A2 = [3 1
		  1 2]
	
	B = [ 9 -4
		 -10 13]
	
	rref(
	[A1[1] A2[1] B[1]
	 A1[2] A2[2] B[2]
	 A1[3] A2[3] B[3]
	 A1[4] A2[4] B[4]])
end

# ╔═╡ 6aab5a5e-2e68-4e1b-b4aa-daaf6d4482bf
md"``x_1 = 3`` and ``x_2 = 2``"

# ╔═╡ 06cbfac7-ad60-4e22-a886-923ffebf04ad
md"""
### PS7 #5
"""

# ╔═╡ 820b0411-9647-455d-badf-804886f89dc2
let
	E1 = [0 -1; 1 0]
	E2 = [1 0; -2 1]
	E3 = [1 0; 0 -1]
	E4 = [1 -3; 0 1]
	C = [2 5 1 0; -1 -3 0 1]
	E1 * C, E2 * E1 * C, E3 * E2 * E1 * C, E4 * (E3 * E2 * E1 * C)
end

# ╔═╡ 26c08822-530b-4d78-9a27-381baa1f4bea
let
	E1 = [0 -1; 1 0]
	E2 = [1 0; -2 1]
	E3 = [1 0; 0 -1]
	E4 = [1 -3; 0 1]
	A = [2 5; -1 -3]
	B = E4 * E3 * E2 * E1
	B * A
end

# ╔═╡ 901848a7-6f10-4730-91f8-aaeadb132b2a
md"""
### PS7 #6

Norm (magnitude) $\|u\| = \sqrt{\vec{u} ⋅ \vec{u}}$

Projection of $\vec{x}$ on $\vec{y}$ is

$\frac{(\vec{x} ⋅ \vec{y})}{\|\vec{y}\|^2} \vec{y}$

Orthogonality is when

$\vec{x} ⋅ \vec{y} = 0$
"""

# ╔═╡ def19b3e-9ec8-4742-b378-415a3df7107d
md"## Problem Set 8"

# ╔═╡ a96cd412-dccf-45bb-89cf-05b24774b9c2
md"""
### PS8 #1

Let

$A = \begin{bmatrix} 1 & -3 & -3 \\ 3 & -1 & 0 \end{bmatrix}$

Find $B$ such that $AB = I$ where $I = \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix}$.
(This is possible because the rank of $A$ is equal to the number of rows.)
"""

# ╔═╡ 88dee7bb-53c3-4041-9712-ca66e6624501
let
	A = [1 -3 -3
		 3 -1  0]
	
	M = [(A * A') I]
	
	C = rref(M)[:,3:4]
	
	B = A' * C
	
	rationalize.(B)
end

# ╔═╡ b8744e45-4c20-4fcc-911a-e5c57c38c2ef
md"""
$A^t = \begin{bmatrix} 1 & 3 \\ -3 & -1 \\ -3 & 0 \end{bmatrix}$

$\begin{align*}
AA^t &= \begin{bmatrix} 1 & -3 & -3 \\ 3 & -1 & 0 \end{bmatrix} \begin{bmatrix} 1 & 3 \\ -3 & -1 \\ -3 & 0 \end{bmatrix} \\
&= \begin{bmatrix}
(1)(1) + (-3)(-3) + (-3)(-3) & (1)(3) + (-3)(-1) + (-3)(0) \\
(3)(1) + (-1)(-3) + (0)(-3) & (3)(3) + (-1)(-1) + (0)(0)
\end{bmatrix} \\
&= \begin{bmatrix}
19 & 6 \\
6 & 10
\end{bmatrix}
\end{align*}$

$\begin{align*}\begin{bmatrix} AA^t & I\end{bmatrix}
&= \begin{bmatrix} 19 & 6 & 1 & 0 \\ 6 & 10 & 0 & 1 \end{bmatrix} \\
R_1 \gets \frac{1}{19} R_1 &\implies \begin{bmatrix} 1 & \frac{6}{19} & \frac{1}{19} & 0 \\ 6 & 10 & 0 & 1 \end{bmatrix} \\
R_2 \gets R_2 - 6 R_1 &\implies \begin{bmatrix} 1 & \frac{6}{19} & \frac{1}{19} & 0 \\ 0 & \frac{154}{19} & -\frac{6}{19} & 1 \end{bmatrix} \\
R_2 \gets \frac{19}{154} R_2 &\implies \begin{bmatrix} 1 & \frac{6}{19} & \frac{1}{19} & 0 \\ 0 & 1 & -\frac{3}{77} & \frac{19}{154} \end{bmatrix} \\
R_1 \gets R_1 - \frac{6}{19} R_2 &\implies \begin{bmatrix} 1 & 0 & \frac{5}{77} & -\frac{3}{77} \\ 0 & 1 & -\frac{3}{77} & \frac{19}{154} \end{bmatrix} \\
\end{align*}$
"""

# ╔═╡ 49a8acd9-0db5-4cb2-9cf5-0437483dbe71
md"""
$C = \begin{bmatrix} \frac{5}{77} & -\frac{3}{77} \\ -\frac{3}{77} & \frac{19}{154} \end{bmatrix}$

$\begin{align*}
B = A^t C &= \begin{bmatrix} 1 & 3 \\ -3 & -1 \\ -3 & 0 \end{bmatrix} \begin{bmatrix} \frac{5}{77} & -\frac{3}{77} \\ -\frac{3}{77} & \frac{19}{154} \end{bmatrix} \\
&= \begin{bmatrix}
(1) \left(\frac{5}{77}\right) + (3) \left(-\frac{3}{77}\right) & (1) \left(-\frac{3}{77}\right) + (3) \left(\frac{19}{154}\right) \\
(-3) \left(\frac{5}{77}\right) + (-1) \left(-\frac{3}{77}\right) & (-3) \left(-\frac{3}{77}\right) + (-1) \left(\frac{19}{154}\right) \\
(-3) \left(\frac{5}{77}\right) + (0) \left(-\frac{3}{77}\right) & (-3) \left(-\frac{3}{77}\right) + (0) \left(\frac{19}{154}\right)
\end{bmatrix} \\
&= \begin{bmatrix}
-\frac{4}{77} & \frac{51}{154} \\
-\frac{12}{77} & -\frac{1}{154} \\
-\frac{15}{77} & \frac{9}{77}
\end{bmatrix}
\end{align*}$
"""

# ╔═╡ 4c1e7119-db6d-413c-ad60-9fdfa0f58d18
3*19//154

# ╔═╡ 12693645-0282-46f3-9a4a-02c122f597fc
md"""
$\begin{align*}
AB &= \begin{bmatrix} 1 & -3 & -3 \\ 3 & -1 & 0 \end{bmatrix} \begin{bmatrix}
-\frac{4}{77} & \frac{51}{154} \\
-\frac{12}{77} & -\frac{1}{154} \\
-\frac{15}{77} & \frac{9}{77}
\end{bmatrix} \\
&= \begin{bmatrix}
(1) \left(-\frac{4}{77}\right) + (-3) \left(-\frac{12}{77}\right) + (-3) \left(-\frac{15}{77}\right) & (1) \left(\frac{51}{154}\right) + (-3) \left(-\frac{1}{154}\right) + (-3) \left(\frac{9}{77}\right) \\
(3) \left(-\frac{4}{77}\right) + (-1) \left(-\frac{12}{77}\right) + (0) \left(-\frac{15}{77}\right) & (3) \left(\frac{51}{154}\right) + (-1) \left(-\frac{1}{154}\right) + (0) \left(\frac{9}{77}\right)
\end{bmatrix} \\
&= \begin{bmatrix}
1 & 0 \\
0 & 1
\end{bmatrix}
\end{align*}$
"""

# ╔═╡ 16e43f06-127b-40d4-a26a-3e5a1e12a847
md"""
### PS8 #2

Let

$A = \begin{bmatrix}
-4 & -4 & 4 \\
0 & 1 & 1
\end{bmatrix}$

Find $B$ such that $AB = I$ where $I$ is as above.
"""

# ╔═╡ 37b94d4a-bd8b-412b-b20b-b4c8453b6d41
let
	A = [-4 -4 4
		  0  1 1]
	
	C = rref([A * A' I])[:,3:4]
	
	B = A' * C
	
	rationalize.(B)
end

# ╔═╡ b18830f5-dc83-4160-b87e-4212fd752175
md"""
$A^t = \begin{bmatrix}
-4 & 0 \\
-4 & 1 \\
4 & 1
\end{bmatrix}$

$\begin{align*}
AA^t &= \begin{bmatrix} -4 & -4 & 4 \\ 0 & 1 & 1 \end{bmatrix}
\begin{bmatrix} -4 & 0 \\ -4 & 1 \\ 4 & 1 \end{bmatrix} \\
&= \begin{bmatrix}
(-4)(-4) + (-4)(-4) + (4)(4) & (-4)(0) + (-4)(1) + (4)(1) \\
(0)(-4) + (1)(-4) + (1)(4) & (0)(0) + (1)(1) + (1)(1)
\end{bmatrix} \\
&= \begin{bmatrix}
16 + 16 + 16 & -4 + 4 \\
-4 + 4 & 1 + 1
\end{bmatrix} \\
&= \begin{bmatrix}
48 & 0 \\
0 & 2
\end{bmatrix}
\end{align*}$

$\begin{align*}
\begin{bmatrix} AA^t & I \end{bmatrix} &= \begin{bmatrix} 48 & 0 & 1 & 0 \\ 0 & 2 & 0 & 1 \end{bmatrix} \\
R_1 \gets \frac{1}{48} R_1, \quad R_2 \gets \frac{1}{2} R_2 &\implies \begin{bmatrix} 1 & 0 & \frac{1}{48} & 0 \\ 0 & 1 & 0 & \frac{1}{2} \end{bmatrix}
\end{align*}$
"""

# ╔═╡ 63db7ceb-b710-4830-acf0-61ede14c2fe5
md"""
$C = \begin{bmatrix} \frac{1}{48} & 0 \\ 0 & \frac{1}{2} \end{bmatrix}$

$\begin{align*}
B = A^t C &= \begin{bmatrix} -4 & 0 \\ -4 & 1 \\ 4 & 1 \end{bmatrix} \begin{bmatrix} \frac{1}{48} & 0 \\ 0 & \frac{1}{2} \end{bmatrix} \\
&= \begin{bmatrix} (-4) \left(\frac{1}{48}\right) & 0 \\ (-4) \left(\frac{1}{48}\right) & \frac{1}{2} \\ \frac{4}{48} & \frac{1}{2} \end{bmatrix} \\
&= \begin{bmatrix} -\frac{1}{12} & 0 \\ -\frac{1}{12} & \frac{1}{2} \\ \frac{1}{12} & \frac{1}{2} \end{bmatrix}
\end{align*}$
"""

# ╔═╡ 3eeb4786-1d9c-4cd3-8d2c-7910418f71bd
md"""
$\begin{align*}
AB &= \begin{bmatrix} -4 & -4 & 4 \\ 0 & 1 & 1 \end{bmatrix} \begin{bmatrix} -\frac{1}{12} & 0 \\ -\frac{1}{12} & \frac{1}{2} \\ \frac{1}{12} & \frac{1}{2} \end{bmatrix} \\
&= \begin{bmatrix} \frac{4}{12} + \frac{4}{12} + \frac{4}{12} & -\frac{4}{2} + \frac{4}{2} \\ -\frac{1}{12} + \frac{1}{12} & \frac{1}{2} + \frac{1}{2} \end{bmatrix} \\
&= \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix}
\end{align*}$
"""

# ╔═╡ 77994a8b-d5e0-4673-899f-d16adefba744
md"""
### PS8 #3

This time $A$ is $2 × 4$

$A = \begin{bmatrix} 1 & -2 & 0 & 1 \\ -1 & 0 & -2 & 1 \end{bmatrix}$

Find $B$ such that $AB = I$.
"""

# ╔═╡ 9b0bf121-c4a7-4dfc-a60e-29b8f616940b
md"""
$A^t = \begin{bmatrix} 1 & -1 \\ -2 & 0 \\ 0 & -2 \\ 1 & 1 \end{bmatrix}$

$\begin{align*}
AA^t &= \begin{bmatrix} 1 & -2 & 0 & 1 \\ -1 & 0 & -2 & 1 \end{bmatrix} \begin{bmatrix} 1 & -1 \\ -2 & 0 \\ 0 & -2 \\ 1 & 1 \end{bmatrix} \\
&= \begin{bmatrix} (1)(1) + (-2)(-2) + (0)(0) + (1)(1) & (1)(-1) + (-2)(0) + (0)(-2) + (1)(1) \\ (-1)(1) + (0)(-2) + (-2)(0) + (1)(1) & (-1)(-1) + (0)(0) + (-2)(-2) + (1)(1) \end{bmatrix} \\
&= \begin{bmatrix} 1 + 4 + 1 & -1 + 1 \\ -1 + 1 & 1 + 4 + 1 \end{bmatrix} \\
&= \begin{bmatrix} 6 & 0 \\ 0 & 6 \end{bmatrix}
\end{align*}$

$\begin{align*}
\begin{bmatrix} AA^t & I \end{bmatrix}
&= \begin{bmatrix} 6 & 0 & 1 & 0 \\ 0 & 6 & 0 & 1 \end{bmatrix} \\
R_1 \gets \frac{1}{6} R_1, \quad R_2 \gets \frac{1}{6} R_2 &\implies \begin{bmatrix} 1 & 0 & \frac{1}{6} & 0 \\ 0 & 1 & 0 & \frac{1}{6} \end{bmatrix}
\end{align*}$
"""

# ╔═╡ 0dc779f0-3d01-4252-a5fd-3ba9d8d0ddd4
md"""
$C = \begin{bmatrix} \frac{1}{6} & 0 \\ 0 & \frac{1}{6} \end{bmatrix}$

$\begin{align*}
B = A^t C &= \begin{bmatrix} 1 & -1 \\ -2 & 0 \\ 0 & -2 \\ 1 & 1 \end{bmatrix} \begin{bmatrix} \frac{1}{6} & 0 \\ 0 & \frac{1}{6} \end{bmatrix} \\
&= \begin{bmatrix} (1) \left(\frac{1}{6}\right) + (-1)(0) & (1)(0) + (-1) \left(\frac{1}{6}\right) \\
(-2) \left(\frac{1}{6}\right) + (0)(0) & (-2)(0) + (0)\left(\frac{1}{6}\right) \\
(0) \left(\frac{1}{6}\right) + (-2)(0) & (0)(0) + (-2)\left(\frac{1}{6}\right) \\
(1) \left(\frac{1}{6}\right) + (1)(0) & (1)(0) + (1)\left(\frac{1}{6}\right)
\end{bmatrix} \\
&= \begin{bmatrix} \frac{1}{6} & -\frac{1}{6} \\ -\frac{1}{3} & 0 \\ 0 & -\frac{1}{3} \\ \frac{1}{6} & \frac{1}{6} \end{bmatrix}
\end{align*}$
"""

# ╔═╡ 8a798a6e-46b7-498e-be2f-253a45b13900
md"""
$\begin{align*}
AB &= \begin{bmatrix} 1 & -2 & 0 & 1 \\ -1 & 0 & -2 & 1 \end{bmatrix} \begin{bmatrix} \frac{1}{6} & -\frac{1}{6} \\ -\frac{1}{3} & 0 \\ 0 & -\frac{1}{3} \\ \frac{1}{6} & \frac{1}{6} \end{bmatrix} \\
&= \begin{bmatrix}
(1) \left(\frac{1}{6}\right) + (-2) \left(-\frac{1}{3}\right) + (0)(0) + (1) \left(\frac{1}{6}\right) & (1) \left(-\frac{1}{6}\right) + (-2)(0) + (0) \left(-\frac{1}{3}\right) + (1) \left(\frac{1}{6}\right) \\
(-1) \left(\frac{1}{6}\right) + (0) \left(-\frac{1}{3}\right) + (-2)(0) + (1) \left(\frac{1}{6}\right) & (-1) \left(-\frac{1}{6}\right) + (0)(0) + (-2) \left(-\frac{1}{3}\right) + (1) \left(\frac{1}{6}\right)
\end{bmatrix} \\
&= \begin{bmatrix} \frac{1}{6} + \frac{2}{3} + \frac{1}{6} & -\frac{1}{6} + \frac{1}{6} \\
-\frac{1}{6} + \frac{1}{6} & \frac{1}{6} + \frac{2}{3} + \frac{1}{6} \end{bmatrix} \\
&= \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix}
\end{align*}$
"""

# ╔═╡ b1a4e8f3-a5a0-478a-bd85-137e8cf4cd59
md"""
### PS8 #4

The rank of the $3 × 4$ matrix

$A = \begin{bmatrix} 1 & -2 & 0 & 1 \\ 3 & -2 & 4 & -1 \\ -1 & 0 & -2 & 1 \end{bmatrix}$

is 2, so something will fail here.
Follow the steps indicated on the side bar for Problem 1 until the method crashes.
(Explain why it does not work.)
"""

# ╔═╡ 32662135-aad8-47e1-bc33-976a3bf5445f
md"""
$A^t = \begin{bmatrix} 1 & 3 & -1 \\ -2 & -2 & 0 \\ 0 & 4 & -2 \\ 1 & -1 & 1 \end{bmatrix}$

$\begin{align*}
AA^t &= \begin{bmatrix} 1 & -2 & 0 & 1 \\ 3 & -2 & 4 & -1 \\ -1 & 0 & -2 & 1 \end{bmatrix} \begin{bmatrix} 1 & 3 & -1 \\ -2 & -2 & 0 \\ 0 & 4 & -2 \\ 1 & -1 & 1 \end{bmatrix} \\
&= \begin{bmatrix} 6 & 6 & 0 \\ 6 & 30 & -12 \\ 0 & -12 & 6 \end{bmatrix}
\end{align*}$

$\begin{bmatrix} AA^t & I \end{bmatrix} = \begin{bmatrix} 6 & 6 & 0 & 1 & 0 & 0 \\ 6 & 30 & -12 & 0 & 1 & 0 \\
0 & -12 & 6 & 0 & 0 & 1 \end{bmatrix}$

$\text{rref}\left(\begin{bmatrix} AA^t & I \end{bmatrix}\right) = \begin{bmatrix} 1 & 0 & \frac{1}{2} & 0 & \frac{1}{6} & \frac{5}{12} \\
0 & 1 & -\frac{1}{2} & 0 & 0 & -\frac{1}{12} \\
0 & 0 & 0 & 1 & -1 & -2
\end{bmatrix}$
"""

# ╔═╡ eadcfe72-cf6c-4d97-855e-be9dac1a0b5b
md"""
The method crashes because $\text{rref}\left(\begin{bmatrix} AA^t & I \end{bmatrix}\right)$ has rank 2 so there is no way to obtain an identity matrix $I$ on the left side of the augmented matrix.
"""

# ╔═╡ 2775318a-9d4f-412a-9c5c-fbb698174070
md"## Problem Set 9"

# ╔═╡ e985723e-144d-480a-bbd3-1eebe1331347
md"""
### PS9 #1
"""

# ╔═╡ a1ce5603-defc-4a95-9943-f38fbf40dfb5
let
	A = [ 0  3
		  6  2
		 -2 -1]
	
	C = A' * A
	
	Δ = C[1,1] * C[2,2] - C[1,2] * C[2,1]
	
	Cⁱ = (1 / Δ) * [C[2,2] -C[1,2]; -C[2,1] C[1,1]]
	
	B = Cⁱ * A'
	
	B * A
end

# ╔═╡ 4455c750-19f2-4628-9201-60c18aecb291
md"""
(a)

$\begin{align*}
B &= C^{-1} A^t \\
BA &= C^{-1} A^t A \\
BA &= C^{-1} C \\
BA &= I
\end{align*}$
"""

# ╔═╡ 63020d78-c97a-4212-b729-87305143d597
md"""
(b)

$\begin{align*}
Ax &= 0 \\
B(Ax) &= 0 \\
(BA) x &= 0
\end{align*}$
"""

# ╔═╡ 67aa6210-dd6b-4437-be6f-91e021eb75ef
md"""
### PS9 #2
"""

# ╔═╡ cf46c827-b872-42ac-b06b-84d93ce28970
let
	u1 = [1; 3; -2]
	u2 = [0; 1; 1]
	u3 = [2; 5; -5]
	
	M = [u1 u2 u3]
	
	rref(M)
end

# ╔═╡ f1ec3d56-f091-4d3f-a992-afb994e0f2a7
md"""
### PS9 #3
"""

# ╔═╡ 358f9c4a-ef4b-478b-b1ae-3f1f42d1da78
let
	v1 = [2; 1; 1]
	v2 = [1; 0; 1]
	v3 = [-2; 1; 1]
	
	M = [v1 v2 v3 zeros(3)]
	
	rref(M)
end

# ╔═╡ 75098ed4-b3b6-41df-b3ce-b4c92bf1f091
md"""
### PS9 #4
"""

# ╔═╡ 253187a1-8942-4ce3-ba7c-c3e1b5275828
let
	u1 = [1; 3; -2]
	u2 = [0; 1; 1]
	u3 = [2; 5; -5]
	
	M = [u1 u2 u3]
	
	rref(M)
	
	2u1 - u2
end

# ╔═╡ f685369a-9768-496a-be6f-f006b70524c5
md"## Problem Set 10"

# ╔═╡ 1ff64173-6746-4029-8014-7afcaddb2e71
md"""
### PS10 #1
"""

# ╔═╡ 29da4f6f-ce23-4092-8aaf-a6af74856b6b
let
	u1 = [1; -1]
	u2 = [-4; 4]
	u3 = [0; 1]
	u4 = [2; -3]
	u5 = [-1; 4]
	A = [u1 u2 u3 u4 u5 zeros(2)]
	rref(A), rref([u1 u3])
end

# ╔═╡ c5e60741-c889-41e9-86ae-f3d301f6bb21
md"The vectors $\vec{u}_2$, $\vec{u}_4$, and $\vec{u}_5$ are redundant. $\vec{u}_1$ and $\vec{u}_3$ are linearly independent."

# ╔═╡ d9970466-8e57-4e98-86ac-114a9450df7b
md"""
### PS10 #2
"""

# ╔═╡ a234609f-b53d-492c-a34c-c264d1e14bab
let
	w1 = [-2; -1; 1]
	w2 = [-3; -6; -3]
	w3 = [0; 2; 2]
	w4 = [2; 4; 2]
	A = [w1 w2 w3 w4 zeros(3)]
	rationalize.(rref(A)), rref([w1 w2])
end

# ╔═╡ 58c6e499-1ff5-42c6-9fb8-02708f56427b
md"The vectors $\vec{w}_3$ and $\vec{w}_4$ are redundant. $\vec{w}_1$ and $\vec{w}_2$ are not linearly independent."

# ╔═╡ 6a8cefbb-60d8-45d4-8d91-3f6896ed6ffa
md"""
### PS10 #3
"""

# ╔═╡ 30993fa2-4fe8-4a65-b642-d53dbd0fcbc5
let
	v1 = [1; -1; -2; 0; 1]
	v2 = [4; -4; -8; 0; 4]
	v3 = [-1; 2; 1; 2; -1]
	v4 = [3; -4; -5; -2; 3]
	v5 = [1; -2; -1; 0; 2]
	A = [v1 v2 v3 v4 v5 zeros(5)]
	Int.(rref(A))
end

# ╔═╡ ac6f7409-9feb-426e-b82a-27333ffaeecd
md"The vectors $\vec{v}_2$ and $\vec{v}_4$ are redundant.
The linearly independent set is $\{\vec{v}_1, \vec{v}_3, \vec{v}_5\}$"

# ╔═╡ d37e2898-05f0-4855-a81d-b06ab113c541
md"""
### PS10 #4
"""

# ╔═╡ b8365990-a6a1-433e-83bf-4950d88c61e6
let
	v1 = [1; 2; 0]
	v2 = [4; 0; 1]
	v3 = [0; 1; 0]
	A = [v1 v2 v3]
	Int.(rref(A))
end

# ╔═╡ 363a9d0f-4e24-4a48-ba34-463e50a5078b
md"No, the rank is 3"

# ╔═╡ b3fb8b64-83eb-44b0-b721-63d6d28637ab
md"## Problem Set 11"

# ╔═╡ 1f2c250c-fd6c-435d-9548-99a11e8d354a
md"""
### PS11 #1

Let

$\vec{e}_1 = \begin{bmatrix} 1 \\ 0 \\ 0 \\ 0 \\ 0 \end{bmatrix}
, \qquad \vec{e}_2 = \begin{bmatrix} 0 \\ 1 \\ 0 \\ 0 \\ 0 \end{bmatrix}
, \qquad \vec{e}_3 = \begin{bmatrix} 0 \\ 0 \\ 1 \\ 0 \\ 0 \end{bmatrix}
, \qquad \vec{e}_4 = \begin{bmatrix} 0 \\ 0 \\ 0 \\ 1 \\ 0 \end{bmatrix}
, \qquad \vec{e}_5 = \begin{bmatrix} 0 \\ 0 \\ 0 \\ 0 \\ 1 \end{bmatrix}$

Explain why $\{e_1, …, e_5\}$ is a basis of $ℝ^5$
"""

# ╔═╡ cdc02ef6-fe47-45aa-a0f3-32a2f8a9f952
md"""
``\{e_1, …, e_5\}`` is the **standard basis** of $ℝ^5$.

Clearly it satisfies the two conditions:

1. ``\text{span} \{e_1, …, e_5\} = V`` where $V$ is a subspace of $ℝ^5$

2. ``e_1, …, e_5`` are linearly independent.
"""

# ╔═╡ c924378b-0d35-4a50-a804-6e03adc1bde6
md"""
### PS11 #2

Let $V$ be a vector subspace of $ℝ^n$.
What is the biggest dimension $V$ can have?
"""

# ╔═╡ bb6f8b89-c199-4393-a7ee-ab5b06d6ad6a
md"""
The biggest dimension $V$ can have is $n$.
"""

# ╔═╡ 6e6fdde8-36eb-44a1-855e-21b676c6ece8
md"""
### PS11 #3

Let $\vec{w}_1, \vec{w}_2, \vec{w}_3$ be linearly independent vectors in $ℝ^n$,
let

$\begin{align*}
\vec{v}_1 &= 3\vec{w}_2 + 2\vec{w}_3 \\
\vec{v}_2 &= 3\vec{w}_1 + 7\vec{w}_2 - \vec{w}_3 \\
\vec{v}_3 &= -3\vec{w}_1 - 3\vec{w}_2 + \vec{w}_3
\end{align*}$

Argue that $v_1$, $v_2$ and $v_3$ are linearly independent.
"""

# ╔═╡ 065df513-1ffa-47eb-88ae-a69a834dcac5
md"""
$x_1 \vec{v}_1 + x_2 \vec{v}_2 + x_3 \vec{v}_3 = 0$

if and only if $x_1, x_2, x_3 = 0$.
To do this we substitute for each $\vec{v}_i$:

$\begin{align*}
x_1 (3\vec{w}_2 + 2\vec{w}_3) + x_2 (3\vec{w}_1 + 7\vec{w}_2 - \vec{w}_3) + x_3(-3\vec{w}_1 - 3\vec{w}_2 + \vec{w}_3) &= 0 \\
(3x_1 \vec{w}_2 + 2x_1 \vec{w}_3) + (3x_2 \vec{w}_1 + 7x_2 \vec{w}_2 - x_2 \vec{w}_3) + (-3x_3 \vec{w}_1 - 3x_3 \vec{w}_2 + x_3 \vec{w}_3) &= 0 \\
(3x_2 - 3x_3) \vec{w}_1 + (3x_1 + 7x_2 - 3x_3) \vec{w}_2 + (2x_1 - x_2 + x_3) \vec{w}_3 &= 0
\end{align*}$

We get the system of equations

$\begin{align*}
3x_2 - 3x_2 &= 0 \\
3x_1 + 7x_2 - 3x_3 &= 0 \\
2x_1 - x_2 + x_3 &= 0
\end{align*}$

The augmented matrix is

$\begin{bmatrix} 0 & 3 & -3 \\ 3 & 7 & -3 \\ 2 & -1 & 1 \end{bmatrix}$
"""

# ╔═╡ d070d235-c7fc-4082-b1f5-0d2cd4c59ea0
let
	A = [ 0  3  -3
		  3  7 -3
		  2 -1  1]
	
	rref(A)
end

# ╔═╡ 4cfd479d-4c64-4351-9324-9db91c5a3178
md"""
### PS11 #4

Let

$A = \begin{bmatrix}
1 & -2 & 3 \\
2 & -4 & 6 \\
-3 & 6 & -9
\end{bmatrix}$

The set of solutions to $Ax = 0$ is a vector subspace of $ℝ^3$.
Find a basis of this set.
"""

# ╔═╡ bf19bf1a-7801-4ffd-96fa-b6d5c281b55d
let
	A = [ 1 -2  3
		  2 -4  6
		 -3  6 -9]
	
	rref([A zeros(3)])
end

# ╔═╡ 52541443-187c-4ef3-a24e-3142d7ad9bbf
md"""
$\begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} = x_2 \begin{bmatrix} 2 \\ 1 \\ 0 \end{bmatrix} + x_3 \begin{bmatrix} -3 \\ 0 \\ 1 \end{bmatrix}$

The basis is

$\left\{\begin{bmatrix} 2 \\ 1 \\ 0 \end{bmatrix}, \begin{bmatrix} -3 \\ 0 \\ 1 \end{bmatrix}\right\}$
"""

# ╔═╡ c2532725-32ea-40e6-9d24-edd0c4b9bda2
md"""
### PS11 #5

View the *rows* of

$A = \begin{bmatrix}
1 & 4 & -1 & 3 & 1 \\
-1 & -4 & 2 & -4 & -2 \\
-2 & -8 & 1 & -5 & -1 \\
0 & 0 & 2 & -2 & 0 \\
1 & 4 & -1 & 3 & 2
\end{bmatrix}$

as elements of $ℝ^5$.
As such they span a vector subspace $V$ of $ℝ^5$.
Find the dimension of that subspace.
"""

# ╔═╡ eef85d24-92ce-492d-87d5-40674a28d2ee
let
	A = [1 4 -1 3 1
		 -1 -4 2 -4 -2
		 -2 -8 1 -5 -1
		 0 0 2 -2 0
		 1 4 -1 3 2]
	
	rref([A' zeros(5)])
end

# ╔═╡ a6ec0c87-245d-4548-8199-5438f3911433
md"""
$\dim(V) = 3$

$\begin{bmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \\ x_5 \end{bmatrix} = \begin{bmatrix} 3x_3 \\ x_3 + x_5 \\ x_3 \\ -\frac{1}{2}x_5 \\ x_5 \end{bmatrix} = \begin{bmatrix} 3x_3 \\ x_3 \\ x_3 \\ 0 \\ 0\end{bmatrix} + \begin{bmatrix} 0 \\ x_5 \\ 0 \\ -\frac{1}{2} x_5 \\ x_5 \end{bmatrix} = x_3 \begin{bmatrix} 3 \\ 1 \\ 1 \\ 0 \\ 0 \end{bmatrix} + x_5 \begin{bmatrix} 0 \\ 1 \\ 0 \\ -\frac{1}{2} \\ 1 \end{bmatrix}$
"""

# ╔═╡ c77df6f8-feea-4f45-8a91-aa307332ec5f
md"## Problem Set 12"

# ╔═╡ dbe32a98-086b-485a-8999-22aed04b3d26
md"""
### PS12 #1

A problem with dot product: Let

$\vec{v}_1 = \begin{bmatrix} \frac{\sqrt{3}}{2\sqrt{2}} \\ \frac{1}{\sqrt{2}} \\ \frac{1}{2\sqrt{2}} \end{bmatrix}
,\quad \vec{v}_2 = \begin{bmatrix} -\frac{\sqrt{3}}{2\sqrt{2}} \\ \frac{1}{\sqrt{2}} \\ -\frac{1}{2\sqrt{2}} \end{bmatrix}
,\quad \vec{v}_3 = \begin{bmatrix} -\frac{1}{2} \\ 0 \\ \frac{\sqrt{3}}{2} \end{bmatrix}$

These vectors are unit vectors that are orthogonal to each other:
$\|\vec{v}_i\| = 1$ and $\vec{v}_i ⋅ \vec{v}_j = 0$ when $i ≠ j$.
You do not need to verify this, but can if you wish.
Suppose

$\vec{u} = a_1 \vec{v}_1 + a_2 \vec{v}_2 + a_3 \vec{v}_3$

with some numbers $a_1$, $a_2$, $a_3$.
Using only properties of the dot product, verify that $\vec{u} ⋅ \vec{v}_j = a_j$.
"""

# ╔═╡ b24cb37a-dfc2-4ca8-ab18-4e39a6054976
md"""
$\vec{u} ⋅ \vec{v}_1 = a_1 \|\vec{v}_1\| + 0 + 0 \implies \vec{u} ⋅ \vec{v}_1 = a_1$

$\vec{u} ⋅ \vec{v}_2 = 0 + a_2 \|\vec{v}_2\| + 0\implies \vec{u} ⋅ \vec{v}_2 = a_2$

$\vec{u} ⋅ \vec{v}_2 = 0 + 0 + a_3 \|\vec{v}_3\|\implies \vec{u} ⋅ \vec{v}_3 = a_3$

Hence

$\vec{u} ⋅ \vec{v}_j = a_j$
"""

# ╔═╡ 60449475-2834-423c-9620-bd1c001c6cca
md"""
### PS12 #2

With the setup of Problem 1, verify that $\vec{v}_1, \vec{v}_2, \vec{v}_3$ are linearly independent.
"""

# ╔═╡ f3afa449-4f62-4e4f-b57e-43eb0090a65b
let
	v1 = [sqrt(3) / 2sqrt(2); 1 / sqrt(2); 1 / 2sqrt(2)]
	v2 = [-sqrt(3) / 2sqrt(2); 1 / sqrt(2); -1 / 2sqrt(2)]
	v3 = [-1 / 2; 0; sqrt(3) / 2]
	
	A = [v1 v2 v3]
	rref(A)
end

# ╔═╡ d882dd0d-7ba8-40e6-9a56-410948ac1e7f
md"""
### PS12 #3

Let $\vec{w}_1, \vec{w}_2, \vec{w}_3$ be linearly independent vectors in $ℝ^n$, let

$A = \begin{bmatrix}
a_{11} & a_{12} & a_{13} \\
a_{21} & a_{22} & a_{23} \\
a_{31} & a_{32} & a_{33}
\end{bmatrix}$

be a matrix with rank 3.
Let

$\begin{align*}
\vec{v}_1 &= a_{11} \vec{w}_1 + a_{21} \vec{w}_2 + a_{31} \vec{w}_3 \\
\vec{v}_2 &= a_{12} \vec{w}_1 + a_{22} \vec{w}_2 + a_{32} \vec{w}_3 \\
\vec{v}_3 &= a_{13} \vec{w}_1 + a_{23} \vec{w}_2 + a_{33} \vec{w}_3
\end{align*}$

Argue that $\vec{v}_1$, $\vec{v}_2$, and $\vec{v}_3$ are linearly independent.
"""

# ╔═╡ 9c2cb2f9-6151-49e3-90db-7e6d9c95e6af
md"""
The $\vec{v}$ vectors are linearly independent because the solution of the system of equations will give rank 3, given that the $\vec{w}$ vectors are linearly independent and that $A$ has rank 3, hence there will be no free variables in the solution of $A\vec{v} = \vec{0}$.
"""

# ╔═╡ 79b287ae-bbcc-40bf-b667-b6214aa6ec30
md"""
### PS12 #4

Suppose the vector subspace $V$ of $ℝ^n$ has $\vec{v}_1, \vec{v}_2, …, \vec{v}_k$ as basis.
Let $\vec{u}$ be a vector in $ℝ^n$ that is not in $V$.
Why are the vectors

$\vec{v}_1, \vec{v}_2, …, \vec{v}_k, \vec{u}$

linearly independent?
"""

# ╔═╡ 03ab5746-802e-424c-97b0-a679290b0f3a
md"""
``\vec{v}_1, \vec{v}_2, …, \vec{v}_k`` are a basis of $V$ which is a vector subspace of $ℝ^n$.
Since $\vec{u}$ is in $ℝ^n$, but not in $V$, then the union of the basis of $V$ and $\vec{u}$ is an extension of the basis of $V$ in $ℝ^n$.
By definition of a basis, the given set of vectors is linearly independent.
"""

# ╔═╡ 1dd87c8f-790d-4cc6-8a2d-3118e7c08901
md"""
### PS12 #5

Let $V$ be a vector subspace of $ℝ^n$ of dimension $k$ with $k < n$.
Can $V$ be all of $ℝ^n$?
Explain.
"""

# ╔═╡ 02609b9a-c8f4-4f38-96e0-acef86136c81
md"""
No, $V$ has a smaller dimension than $ℝ^n$ so it has at most $k$ linearly independent vectors, thus it is impossible for $V$ to be all of $ℝ^n$.
"""

# ╔═╡ f26db46a-d1db-46f5-9e71-b0708834c434
md"## Problem Set #13"

# ╔═╡ e810e3e4-256f-481e-9061-8c1d67c28a7e
md"""
### PS13 #1

Let

$\vec{v} = \begin{bmatrix} \frac{\sqrt{3}}{2 \sqrt{2}} \\ \frac{1}{\sqrt{2}} \\ \frac{1}{2 \sqrt{2}} \end{bmatrix}$

and let $T : ℝ^3 → ℝ^3$ be the operator

$T(\vec{x}) = \text{proj}_{\vec{v}}(\vec{x}).$

Verify that $T$ is a linear operator.
"""

# ╔═╡ e782342f-b07c-4bec-8d9a-5a7a9effff3a
md"""
Preserves addition:

$\begin{align*}
\text{proj}_{\vec{v}}(\vec{x} + \vec{y}) &= \text{proj}_{\vec{v}}{\vec{x}} + \text{proj}_{\vec{v}}{\vec{y}} \\
\frac{(\vec{x} + \vec{y}) ⋅ \vec{v}}{\vec{v} ⋅ \vec{v}} \vec{v} &= \frac{\vec{x} ⋅ \vec{v}}{\vec{v} ⋅ \vec{v}} \vec{v} + \frac{\vec{y} ⋅ \vec{v}}{\vec{v} ⋅ \vec{v}} \vec{v} \\
\left(\frac{\vec{x} ⋅ \vec{v} + \vec{y} ⋅ \vec{v}}{\vec{v} ⋅ \vec{v}}\right) \vec{v} &= \left(\frac{\vec{x} ⋅ \vec{v} + \vec{y} ⋅ \vec{v}}{\vec{v} ⋅ \vec{v}}\right) \vec{v}
\end{align*}$

Preserves scalar multiplication:

$\begin{align*}
\text{proj}_{\vec{v}}(k\vec{x}) &= k\text{proj}_{\vec{v}}(\vec{x}) \\
\frac{k\vec{x} ⋅ \vec{v}}{\vec{v} ⋅ \vec{v}} \vec{v} &= k \left(\frac{\vec{x} ⋅ \vec{v}}{\vec{v} ⋅ \vec{v}} \vec{v}\right) \\
k \left(\frac{\vec{x} ⋅ \vec{v}}{\vec{v} ⋅ \vec{v}} \vec{v}\right) &= k \left(\frac{\vec{x} ⋅ \vec{v}}{\vec{v} ⋅ \vec{v}} \vec{v}\right)
\end{align*}$
"""

# ╔═╡ 4ac5a171-a82a-4a49-8ffd-43ba3bb71131
md"""
### PS13 #2

Let $\vec{v}_1$, $\vec{v}_2$ be a basis of $ℝ^2$, let $\vec{w}_1, \vec{w}_2, \vec{w}_3$ be one of $ℝ^3$.
Suppose $T : ℝ^2 → ℝ^3$ is given by the formula

$T(s_1 \vec{v}_1 + s_2 \vec{v}_2) = s_1 (6w_2 - 2w_3) + s_2 (3w_1 + 2w_2 - w_3).$

This is a linear transformation.

(a) Find $T(\vec{v}_1)$ and $T(\vec{v}_2)$.

(b) Find the matrix $A$ of $T$ with respect to the given bases.

(c) Find the rank of $A$.
"""

# ╔═╡ 31e5b3c4-1f85-412f-8c0b-eb351abb3eb9
md"""
**(a)**

$T(\vec{v}_1) = 6w_2 - 2w_3$

$T(\vec{v}_2) = 3w_1 + 2w_2 - w_3$
"""

# ╔═╡ 1b28eba2-7d10-41bd-a00b-a970369aa37a
md"""
**(b)**

$A = \begin{bmatrix}
0 & 1 \\
6 & 2 \\
-2 & -1
\end{bmatrix}$
"""

# ╔═╡ 98663ba6-3502-443b-9d9c-f32cddb9bd74
rref([0 1
	  6 2
	 -2 -1])

# ╔═╡ b6ae668c-aa4c-4ef0-82ab-2b1a9669997e
md"""
**(c)** The rank of $A$ is 2.
"""

# ╔═╡ d3c3260b-a38e-4e24-a1f6-8e42ca7b74ff
md"""
### PS13 #3

Let $T : ℝ^4 → ℝ^4$ be given by

$T{\left(\begin{bmatrix} a_0 \\ a_1 \\ a_2 \\ a_3 \end{bmatrix}\right)} = \begin{bmatrix} a_1 \\ 2a_2 \\ 3a_3 \\ 0 \end{bmatrix}.$

Write the matrix of $T$ with respect to the basis

$\vec{e}_1 = \begin{bmatrix} 1 \\ 0 \\ 0 \\ 0 \end{bmatrix},
\quad \vec{e}_2 = \begin{bmatrix} 0 \\ 1 \\ 0 \\ 0 \end{bmatrix},
\quad \vec{e}_3 = \begin{bmatrix} 0  \\ 0 \\ 1 \\ 0 \end{bmatrix},
\quad \vec{e}_4 = \begin{bmatrix} 0 \\ 0 \\ 0 \\ 1 \end{bmatrix}.$
"""

# ╔═╡ c5d083b9-0340-4a17-9b07-35638fbc5cef
md"""
$A = \begin{bmatrix}
1 & 0 & 0 & 0 \\
0 & 2 & 0 & 0 \\
0 & 0 & 3 & 0 \\
0 & 0 & 0 & 0
\end{bmatrix}$
"""

# ╔═╡ 0dc04569-51db-472b-835c-3d288c3d37ef
md"""
### PS13 #4

Let $P$ be the set of all polynomials in the variable $t$ of degree $≤ 3$:

$P = \{a_0 + a_1 t + a_2 t^2 + a_3 t^3\}.$

This vector space has as basis the polynomials $1, t, t^2, t^3$.
Let $D : P → P$ be the transformation $D(p) = p'$.
Find the matrix of $D$ with respect to this basis.
"""

# ╔═╡ ee7f6086-c337-47d8-b86f-a320f7596ee6
md"""
The basis can be represented as

$\vec{e}_1 = \begin{bmatrix} 1 \\ 0 \\ 0 \\ 0 \end{bmatrix},
\quad \vec{e}_2 = \begin{bmatrix} 0 \\ t \\ 0 \\ 0 \end{bmatrix},
\quad \vec{e}_3 = \begin{bmatrix} 0  \\ 0 \\ t^2 \\ 0 \end{bmatrix},
\quad \vec{e}_4 = \begin{bmatrix} 0 \\ 0 \\ 0 \\ t^3 \end{bmatrix}.$

$A = \begin{bmatrix}
0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 \\
0 & 0 & 2 & 0 \\
0 & 0 & 0 & 3
\end{bmatrix}$
"""

# ╔═╡ b659e17b-7ae7-4297-83e7-f29b4e8568ec
md"""
### PS13 #5

This problem is again about bases.
Let $\vec{w}_1, \vec{w}_2, \vec{w}_3$ be a basis of $ℝ^3$, let

$\begin{align*}
\vec{v}_1 &= \vec{w}_1 - \vec{w}_3 \\
\vec{v}_2 &= \vec{w}_2 + 2\vec{w}_3 \\
\vec{v}_3 &= 2\vec{w}_1 + \vec{w}_2 + \vec{w}_3
\end{align*}$

It so happens that then $\vec{v}_1, \vec{v}_2, \vec{v}_3$ are also linearly independent (you don't need to verify this).
Write $\vec{w}_1$ in terms of $\vec{v}_1, \vec{v}_2, \vec{v}_3$
"""

# ╔═╡ 2cd4c040-ed9e-48a8-bf88-83e00e1b1978
md"""
$\begin{align*}
\vec{v}_2 - \vec{v}_3 &= -2\vec{w}_1 + \vec{w}_3 \\
\vec{v}_1 + (\vec{v}_2 - \vec{v}_3) &= -\vec{w}_1 \\
\vec{w}_1 &= -\vec{v}_1 - \vec{v}_2 + \vec{v}_3
\end{align*}$
"""

# ╔═╡ 819b1e3e-77e9-45a6-83cf-60ee4c3af2c1
md"## Problem Set 14"

# ╔═╡ 796c2bf5-53c7-4b55-9584-bfcc4e09482a
md"""
### PS14 #1

Let $V = \text{span}\{e^t \cos{t} \sin{t}, e^t \cos^2{t}, e^t \sin^2{t}\}$, that is, $V$ consists of all functions

$v = a_1 e^t \cos{t} \sin{t} + a_2 e^t \cos^2{t} + a_3 e^t \sin^2{t}$

with $a_1, a_2, a_3 ∈ ℝ$.
Let $D : V → V$ be differentiation:
if $v ∈ V$, then $T(v) = v'$.
Let

$v_1 = e^t \cos{t} \sin{t}, \quad v_2 = e^t \cos^2{t}, \quad v_3 = e^t \sin^2{t}.$

These elements of $V$ are linearly independent (don't prove this) so they form a basis of $V$.
Find the matrix of $D$ with respect to this basis.
"""

# ╔═╡ 5c06aca8-1ec1-4b06-89e8-fd32a3be7e46
md"""
${v_1}' = e^t \cos{t} \sin{t} + e^t \cos^2{t} - e^t \sin^2{t} = v_1 + v_2 - v_3$

${v_2}' = e^t \cos^2{t} - 2e^t \cos{t} \sin{t} = -2v_1 + v_2$

${v_3}' = e^t \sin^2{t} + 2e^t \sin{t} \cos{t} = 2v_1 + v_3$

$D = \begin{bmatrix}
1 & -2 & 2 \\
1 & 1 & 0 \\
-1 & 0 & 1
\end{bmatrix}$
"""

# ╔═╡ 877a5b7c-83b2-45d7-af18-e54c45cc2bbe
md"## Problem Set 15"

# ╔═╡ f0ec76f1-4d0e-4d70-b1cc-a31efc6834aa
md"""
#### PS15 #1

Let $R_\theta : ℝ^2 → ℝ^2$ be the linear transformation given by the formula

$R_\theta \left(\begin{bmatrix} x_1 \\ x_2 \end{bmatrix}\right) = \begin{bmatrix} x_1 \cos{\theta} - x_2 \sin{\theta} \\ x_1 \sin{\theta} + x_2 \cos{\theta} \end{bmatrix}\;.$

Let $R_\psi$ be similarly defined.

(a) Evaluate

$R_\psi \left(R_\theta \left(\begin{bmatrix} x_1 \\ x_2 \end{bmatrix}\right)\right).$

(b) Use trigonometric identities to verify that

$R_\psi\left(R_\theta\left(\begin{bmatrix} x_1 \\ x_2 \end{bmatrix}\right)\right) = R_{\theta + \psi} \left(\begin{bmatrix} x_1 \\ x_2 \end{bmatrix}\right)$
"""

# ╔═╡ ed28d412-d601-4852-bb44-e93a029cb3c4
md"""
$\begin{align*}
R_{\psi}\left(R_\theta\left(\begin{bmatrix} x1 \\ x2 \end{bmatrix}\right)\right) &=  \begin{bmatrix} (x_1 \cos{\theta} - x_2 \sin{\theta}) \cos{\theta} - (x_1 \sin{\theta} -+ x_2 \cos{\theta}) \sin{\theta} \\
(x_1 \cos{\theta} - x_2 \sin{\theta}) \sin{\theta} + (x_1 \sin{\theta} -+ x_2 \cos{\theta}) \cos{\theta}
\end{bmatrix}
\end{align*}$
"""

# ╔═╡ 2eba100d-20c0-4b80-89a9-18f113f88cc3
md"## Problem Set 16"

# ╔═╡ 9daaf6fc-a3fa-4914-9f0e-2ce0f5df2004
md"""
### PS16 #1

Let

$A = \begin{bmatrix} -1 & 3 & 3 & 3 & -2 \\ 1 & -3 & -5 & -3 & 2 \\ -3 & 9 & 13 & 9 & -7 \end{bmatrix}$

(a) Find the rank and nullity of $A$.

(b) Find a basis for $\text{Col}(A)$ ($=$ the span of the columns of $A$)

(c) Find a basis for $\text{Null}(A)$ ($= \{\vec{x} : A\vec{x} = 0\}$).

(d) Find a basis for $\text{Row}(A)$ ($=$ the span of the rows of $A$).

(e) Find the rank and nullity of $A^t$.
"""

# ╔═╡ 6a4bb909-caf0-4b96-ba15-7c443499cbda
let
	A = [-1  3  3  3 -2
		  1 -3 -5 -3  2
		 -3  9 13  9 -7]
	
	:a => (rank(A), ndims(nullspace(A)))
end

# ╔═╡ 403b70ab-f14f-4725-8434-2b9e644955b5
let
	A = [-1  3  3  3 -2
		  1 -3 -5 -3  2
		 -3  9 13  9 -7]
	
	M = rref(A)
	
	basis = []
	
	for j ∈ 1:size(M)[2]
		if 1 ∈ M[:,j]
			push!(basis, A[:,j])
		end
	end
	
	:b => basis
end

# ╔═╡ 5cbb9223-33b0-406c-9391-2fa933440df7
let
	A = [-1  3  3  3 -2
		  1 -3 -5 -3  2
		 -3  9 13  9 -7]
	
	M = rref([A zeros(3)])
	
	:c => round.(Int, [-[M[:,2] M[:,4]][1:size(M)[1]-1,:]; I])
end

# ╔═╡ e845cbeb-37d6-464f-989c-f7adf6b8d6fb
let
	A = [-1  3  3  3 -2
		  1 -3 -5 -3  2
		 -3  9 13  9 -7]
	
	M = rref(A)
	
	basis = []
	
	for j ∈ 1:size(M)[1]
		if 1 ∈ M[j,:]
			push!(basis, round.(Int, M[j,:]))
		end
	end
	
	:d => basis
end

# ╔═╡ abeb3c03-29a2-44c1-b542-2af8eb35b8cc
let
	A = [-1  3  3  3 -2
		  1 -3 -5 -3  2
		 -3  9 13  9 -7]
	
	:e => (rank(A'), ndims(nullspace(A')))
	
	A'
end

# ╔═╡ ea442e30-106c-4f4b-9cc5-594dc9bd18a7
md"""
### PS16 #2

Suppose $A$ is an $m × n$ matrix.

(a) What is the biggest possible rank of $A$?

(b) Suppose that $A$ has rank $k$.
    What is the nullity of $A$?

(c) What is the biggest possible rank of $A^t$?

(d) Suppose that $A$ has rank $k$. Find the rank and nullity of $A^t$.
"""

# ╔═╡ 06ddaef2-dc36-4b11-b5a5-d78b72991ff4
md"""
(a) $n$

(b) $n - k$

(c) $n$

(d) $k, m - k$
"""

# ╔═╡ 54ddb1b0-b0a9-446c-a814-a551b0901bd1
md"""
### PS16 #3

Suppose $u_1, u_2, u_3$ are
"""

# ╔═╡ 4ede0ef8-8bc9-4a5f-820a-05e00b3bf82a
let
	A = [-2 -3 -2 -1
		 -1 -2 -3 2
		 1 -1 -3 1]
	
	round.(Int, rref(A))
end

# ╔═╡ 25d8aa91-7e21-48d6-871c-0f3aae09f07f
md"""
### PS16 #4
"""

# ╔═╡ cbc18af4-5a2c-4971-9939-5a757f90033c
md"""
### PS16 #5
"""

# ╔═╡ fb17c046-b1dc-46a3-9dd1-c7fe5792efd3
md"""
### PS16 #6
"""

# ╔═╡ 5b5ac2cc-ed80-4adf-9b02-016586ec5a68
md"""
### PS16 #7
"""

# ╔═╡ e9caa6fa-49bb-4816-9b78-84a7c43e93c0
md"""
### PS16 #8
"""

# ╔═╡ 54c05e91-35ef-44c4-b543-193ca2e98d1f
let
	A = [2 -3 1
		 1 0 -1]
	
	A'
end

# ╔═╡ ee649240-cd68-480c-b055-2ff6e9fd0d9a
md"""
### PS16 #9
"""

# ╔═╡ 516a0f4a-dcf5-4e92-8f91-5a1cb3336a5c
md"""
### PS16 #10
"""

# ╔═╡ 0451fa11-85cd-4e3c-b084-83e04a1276b3
let
	A = [2 -3 1
		 1 0 -1]
	
	det(A*A') == 0 ? "not invertible" : "invertible"
end

# ╔═╡ f8a4533c-98e3-4ecf-8514-a081fc51da7a
let
	A = [2 -3 1
		 1 0 -1]
	
	A' * 27 * inv(A*A')
end

# ╔═╡ ba4fd86b-2610-45ca-bb54-487b3029c982
md"""
### PS16 #11
"""

# ╔═╡ ce0942ab-074b-46e7-b4c2-e217022d68f7


# ╔═╡ b00d86b3-2d9a-4478-98bb-f09ca8d39d72
md"""
### PS16 #12
"""

# ╔═╡ d3374ce9-e308-4d93-9577-2b50e16805b7
let
	A = [-1 2
		  0 1
		  1 1]
	
	inv(A' * A) * 11
end

# ╔═╡ 583a73d4-5e64-4e6f-958e-74590ce21ce8
md"## Problem Set 17"

# ╔═╡ ff1cc355-2f53-43c2-860b-5617c206a2ec
md"""
### PS17 #1

Let

$A = \begin{bmatrix} 1 & -2 & 3 \\ 0 & 2 & 4 \\ 0 & -1 & 1 \end{bmatrix}\;.$

Find:

(a) The determinant of $A$.

(b) Each of the 9 minors of $A$, use them to make a matrix.

(c) Use the minors to write down the cofactor matrix.
    Call this matrix $C$.

(d) Let $B = C^t$ (the transpose of the cofactor matrix, which is called the adjugate of $A$).
    Find $AB$ and $BA$.
"""

# ╔═╡ 6427d602-c735-4f30-b41f-f102fcf7fa89
md"""
**(a)**

$\det(A) = 2 + 0 + 0 - 0 - (-4) - 0 = 2 + 4 = 6$

**(b)**

$\begin{bmatrix} 6 & 0 & 0 \\ 1 & 1 & -1 \\ -14 & 4 & 2 \end{bmatrix}$

**(c)**

$C = \begin{bmatrix} 6 & 0 & 0 \\ -1 & 1 & 1 \\ -14 & -4 & 2 \end{bmatrix}$

**(d)**

$B = \begin{bmatrix} 6 & -1 & -14 \\ 0 & 1 & -4 \\ 0 & 1 & 2 \end{bmatrix}$

$AB = \begin{bmatrix} 6 & 0 & 0 \\ 0 & 6 & 0 \\ 0 & 0 & 6 \end{bmatrix}$

$BA = \begin{bmatrix} 6 & 0 & 0 \\ 0 & 6 & 0 \\ 0 & 0 & 6 \end{bmatrix}$
"""

# ╔═╡ 3c5bb1e4-1d9f-4e52-97b2-5cc41ab2c106
let
	A = [1 -2 3
		 0  2 4
		 0 -1 1]
	C = Int.(cof(A))
	B = C'
	
	:a => det(A), :b => Int.(minor(A)), :c => Int.(cof(A)), :d => (A*B, B*A)
end

# ╔═╡ 1af14435-984c-49b2-b088-ff4505d10bd1
md"""
### PS17 #2

Let

$A = \begin{bmatrix} 1 & -2 & 0 \\ 0 & -2 & 1 \\ 1 & 0 & -2 \end{bmatrix}\;.$

Like in the previous problem, find:

(a) The determinant of $A$.

(b) Each of the 9 minors of $A$, use them to make a matrix.

(c) Use the minors to write down the cofactor matrix.
    Call this matrix $C$.

(d) Let $B = C^t$ (the transpose of the cofactor matrix, which is called the adjugate of $A$).
    Find $AB$ and $BA$.
"""

# ╔═╡ 68a2b285-fd90-4225-af15-10f54552c9f8
md"""
**(a)**

$\det(A) = 4 + (-2) + 0 - 0 - 0 - 0 = 4 - 2 = 2$

**(b)**

$\begin{bmatrix} 4 & -1 & 2 \\ 4 & -2 & 2 \\ -2 & 1 & -2 \end{bmatrix}$

**(c)**

$\begin{bmatrix} 4 & 1 & 2 \\ -4 & -2 & -2 \\ -2 & -1 & -2 \end{bmatrix}$

**(d)**

$\begin{bmatrix} 4 & -4 & -2 \\ 1 & -2 & -1 \\ 2 & -2 & -2 \end{bmatrix}$

$AB = \begin{bmatrix} 2 & 0 & 0 \\ 0 & 2 & 0 \\ 0 & 0 & 2 \end{bmatrix}$

$BA = \begin{bmatrix} 2 & 0 & 0 \\ 0 & 2 & 0 \\ 0 & 0 & 2 \end{bmatrix}$
"""

# ╔═╡ 63a22534-b547-4a50-ac80-197b279f4bd9
let
	A = [1 -2  0
		 0 -2  1
		 1  0 -2]
	
	C = Int.(cof(A))
	
	B = C'
	
	det(A), Int.(minor(A)), A*B, B*A
end

# ╔═╡ e31c86fe-713d-4301-90d4-d966d79f10e4
md"""
### PS17 #3

Let

$A = \begin{bmatrix} 2 & 1 & 0 \\ 1 & -1 & 4 \\ 1 & 2 & -4 \end{bmatrix} \;.$

Let $B$ be adjugate of $A$, find $AB$
"""

# ╔═╡ a1bbacbe-7683-4b0c-aac0-e60c2ad95e8f
let
	A = [2  1  0
		 1 -1  4
		 1  2 -4]
	
	B = Int.(cof(A)')
	
	A*B

	cof(A)
end

# ╔═╡ 6b6573e4-2f61-438d-9763-b3a4200ff790
md"""
### PS17 #4

Let

$A = \begin{bmatrix} -1 & 2 & 1 \\ 0 & -2 & 1 \\ 2 & -1 & 0 \end{bmatrix}, \quad B = \begin{bmatrix} -1 & 2 & 1 \\ x_1 & x_2 & x_3 \\ 2 & -1 & 0 \end{bmatrix}, \quad C = 
\begin{bmatrix} -1 & 2 & 1 \\ x_1 & x_2 - 2 & x_3 + 1 \\ 2 & -1 & 0 \end{bmatrix}$

Evaluate by hand the determinants of $A$, $B$, and $C$.
Did you get

$\det(C) = \det(A) + \det(B)?$
"""

# ╔═╡ 96fd1a78-be44-4e26-950b-2efb1b78973b
md"""
$\begin{align*}
\det(A) &= 0 \begin{vmatrix} 2 & 1 \\ -1 & 0 \end{vmatrix} + 2 \begin{vmatrix} -1 & 1 \\ 2 & 0 \end{vmatrix} + 1 \begin{vmatrix} -1 & 2 \\ 2 & -1 \end{vmatrix} \\
&= -0 - 2(-2) - (-3) \\
&= 7
\end{align*}$

$\begin{align*}
\det(B) &= -x_1 \begin{vmatrix} 2 & 1 \\ -1 & 0 \end{vmatrix} + x_2 \begin{vmatrix} -1 & 1 \\ 2 & 0 \end{vmatrix} - x_3 \begin{vmatrix} -1 & 2 \\ 2 & -1 \end{vmatrix} \\
&= -x_1(1) + x_2(-2) - x_3(-3) \\
&= -x_1 - 2x_2 + 3x_3
\end{align*}$

$\det(A) + \det(B) = -x_1 - 2x_2 + 3x_3 + 7$

$\begin{align*}
\det(C) &= -x_1 \begin{vmatrix} 2 & 1 \\ -1 & 0 \end{vmatrix} + (x_2 - 2) \begin{vmatrix} -1 & 1 \\ 2 & 0 \end{vmatrix} - (x_3 + 1) \begin{vmatrix} -1 & 2 \\ 2 & -1 \end{vmatrix} \\
&= -x_1(1) + (x_2 - 2)(-2) - (x_3 + 1)(-3) \\
&= -x_1 - 2x_2 + 4 - (-3x_3 - 3) \\
&= -x_1 - 2x_2 + 3x_3 + 7
\end{align*}$

$\det(C) = \det(A) + \det(B)$
"""

# ╔═╡ 2069b9b9-961e-4665-b2a4-5badf66a99d8
md"""
### PS17 #5

Let

$A = \begin{bmatrix} -1 & 2 & 1 \\ 0 & -2 & 1 \\ 2 & -1 & 0 \end{bmatrix}, \quad B = \begin{bmatrix} -1 & 2k & 1 \\ 0 & -2k & 1 \\ 2 & -4k & 0 \end{bmatrix}\;.$

Find $\det(A)$ and $\det(B)$.
"""

# ╔═╡ 18620a95-56db-4328-a3d4-b0343344d0b4
let
	A = [-1  2 1
		  0 -2 1
		  2 -1 0]
	
	det(A)
end

# ╔═╡ 130e838f-c6ac-4b54-9288-352e855de012
md"""
$\det(B) = 2k \det(A) = 14k$
"""

# ╔═╡ b0e4d322-8b8b-489d-bc85-d8d95bba024f
md"""
### PS17 #6

For the matrix $A$ in Problem 2, find $\det(A^t)$.
"""

# ╔═╡ ada1b995-c6e1-4f78-9a20-443adc40dc89
md"$\det(A) = \det(A^t) = 2$"

# ╔═╡ a0f249d1-2425-44b6-9454-784175d12657
md"## Problem Set 18"

# ╔═╡ c90450e4-2d02-435b-a7b3-5b9316f4d0a2
md"""
### PS18 #1

Let

$E_1 = \begin{bmatrix} 1 & 0 & 0 \\ -1 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix},
E_2 = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ -3 & 0 & 1 \end{bmatrix},
E_3 = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -1 & 1 \end{bmatrix},$
$E_4 = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & -1 \\ 0 & 0 & 1 \end{bmatrix},
E_5 = \begin{bmatrix} 1 & 0 & -2 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix},
E_6 = \begin{bmatrix} \frac{1}{2} & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix}$


(a) Find the determinant of each of these 6 matrices and their product.

(b) Evaluate $B = E_6 E_5 E_4 E_3 E_2 E_1$. Find $\det(B)$ directly.

(c) Let
$A = \begin{bmatrix} 2 & 0 & 2 \\ 2 & 1 & 3 \\ 6 & 1 & 8 \end{bmatrix}$.
    Find $\det(A)$ and $BA$.
"""

# ╔═╡ fe291037-01a5-4fba-978a-13c4af231aa7
md"""
**(a)**

$\det(E_1) = 1, \det(E_2) = 1, \det(E_3) = 1, \det(E_4) = 1, \det(E_5) = 1, \det(E_6) = \frac{1}{2}$
"""

# ╔═╡ ecfab479-361c-41a4-b304-433702da82dd
let
	E1 = [1 0 0
		  -1 1 0
		  0 0 1]
	E2 = [1 0 0
		  0 1 0
		  -3 0 1]
	E3 = [1 0 0
		  0 1 0
		  0 -1 1]
	E4 = [1 0 0
		  0 1 -1
		  0 0 1]
	E5 = [1 0 -2
		  0 1 0
		  0 0 1]
	E6 = [1//2 0 0
		  0 1 0
		  0 0 1]
	B = E6 * E5 * E4 * E3 * E2 * E1

	det.([E1, E2, E3, E4, E5, E6]), B, det(B)
end

# ╔═╡ fecee897-5961-4042-82d8-0a9facda5359
md"""
**(b)**

$B = E_6 E_5 E_4 E_3 E_2 E_1 = \begin{bmatrix} \frac{5}{2} & 1 & -1 \\ 1 & 2 & -1 \\ -2 & -1 & 1 \end{bmatrix}$

$\det(B) = \frac{1}{2}$
"""

# ╔═╡ 616f4cbd-cb7d-4042-b080-3a7e3f4a8dcf
md"""
**(c)**

$1 = \det(B) \det(A) \implies \det(A) = 2$

$BA = I$
"""

# ╔═╡ 708821d9-f22b-4765-88bb-72b88b57d87b
let
	A = [2 0 2
		 2 1 3
		 6 1 8]
	
	B = [5//2 1 -1
		 1 2 -1
		 -2 -1 1]
	
	det(A), B*A
end

# ╔═╡ 6f7c5fbf-b0d2-4912-afc6-65f64e2ebb7e
md"""
### PS18 #2

With the matrix $A$ and the $E_i$ of the previous problem, find each of

$E_1 A,\; E_2 E_1 A, …, E_6 ⋯ E_2 E_1 A.$
"""

# ╔═╡ 2f665e7a-daff-4320-b570-825f4e8fcc66
md"""
$\begin{align*}
E_1 A &= \begin{bmatrix} 2 & 0 & 2 \\ 0 & 1 & 1 \\ 6 & 1 & 8 \end{bmatrix} \\
E_2 E_1 A &= \begin{bmatrix} 2 & 0 & 2 \\ 0 & 1 & 1 \\ 0 & 1 & 2 \end{bmatrix} \\
E_3 E_2 E_1 A &= \begin{bmatrix} 2 & 0 & 2 \\ 0 & 1 & 1 \\ 0 & 0 & 1 \end{bmatrix} \\
E_4 E_3 E_2 E_1 A &= \begin{bmatrix} 2 & 0 & 2 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix} \\
E_5 E_4 E_3 E_2 E_1 A &= \begin{bmatrix} 2 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix} \\
E_6 E_5 E_4 E_3 E_2E_1 A &= \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix} \\
\end{align*}$
"""

# ╔═╡ 8ef3cd15-a03c-42f0-aa71-08b512fa2346
let
	E1 = [1 0 0
		  -1 1 0
		  0 0 1]
	E2 = [1 0 0
		  0 1 0
		  -3 0 1]
	E3 = [1 0 0
		  0 1 0
		  0 -1 1]
	E4 = [1 0 0
		  0 1 -1
		  0 0 1]
	E5 = [1 0 -2
		  0 1 0
		  0 0 1]
	E6 = [1//2 0 0
		  0 1 0
		  0 0 1]
	A = [2 0 2
		 2 1 3
		 6 1 8]
	E1 * A, E2 * E1 * A, E3 * E2 * E1 * A, E4 * E3 * E2 * E1 * A, E5 * E4 * E3 * E2 * E1 * A, Int.(E6 * E5 * E4 * E3 * E2 * E1 * A)
end

# ╔═╡ 654ca869-fd5f-468f-a19b-ab4f7bfd6733
md"""
### PS18 #3

Let 

$Q = \begin{bmatrix} 2 & 0 & 2 & 4 \\ 2 & 1 & 3 & 3 \\ 6 & 1 & 8 & 2 \end{bmatrix}$

Multiply $Q$ on the left by the matrix $B$ in 1b (that is, find $BQ$).
What does the last column of the result represent?
"""

# ╔═╡ d0767f8c-28a1-4c7d-a261-1c5950f39c2a
let
	B = [5//2 1 -1
		 1 2 -1
		 -2 -1 1]
	
	Q = [2 0 2 4
		 2 1 3 3
		 6 1 8 2]
	
	B * Q
end

# ╔═╡ 9783c6b9-a4cf-49c7-bdb5-b253de3a3f08
md"""
$BQ = \begin{bmatrix} 1 & 0 & 0 & 11 \\ 0 & 1 & 0 & 8 \\ 0 & 0 & 1 & -9 \end{bmatrix}$

The last column of the result represents the solution of $A\vec{x} = \vec{v}$ where $\vec{v} = \begin{bmatrix} 4 \\ 3 \\ 2 \end{bmatrix}$.
"""

# ╔═╡ f836b6ae-f4dd-4799-9198-d79e4a48e8d6
md"## Problem Set 19"

# ╔═╡ 45088ee2-ef86-481b-8ef0-0bf2a2482a38
md"""
### PS19 #1

The matrix of minors of

$A = \begin{bmatrix} 3 & -3 & 4 \\ -3 & -1 & 5 \\ -2 & 0 & -2 \end{bmatrix}$

is

$M = \begin{bmatrix} 2 & 16 & -2 \\ 6 & 2 & -6 \\ -11 & 27 & -12 \end{bmatrix}$

Use this information and that $\det(A) = 46$ to find the inverse of $A$.
"""

# ╔═╡ bee81924-5aa0-4ffd-afe1-0a52a70067af
md"""
$\text{cof}(A) = \begin{bmatrix} 2 & -16 & -2 \\ -6 & 2 & 6 \\ -11 & -27 & -12 \end{bmatrix}$

$\text{adj}(A) = \begin{bmatrix} 2 & -6 & -11 \\ -16 & 2 & -27 \\ -2 & 6 & -12 \end{bmatrix}$

$A^{-1} = \frac{1}{\det(A)} \text{adj}(A) = \frac{1}{46} \begin{bmatrix} 2 & -6 & -11 \\ -16 & 2 & -27 \\ -2 & 6 & -12 \end{bmatrix}$
"""

# ╔═╡ 40622c8c-1bda-41cd-8050-9347e003890d
md"""
### PS19 #2

The cofactor matrix of

$A = \begin{bmatrix} 1 & 0 & -2 \\ 5 & -5 & 0 \\ 0 & 3 & -2 \end{bmatrix}$

is mostly

$C = \begin{bmatrix} 10 & 10 & 15 \\ \fbox{\, ⠀ \,} & -2 & -3 \\ -10 & -10 & \fbox{\, ⠀ \,}\end{bmatrix}$

Find the missing entries, then multiply $C^t A$.

Based on your answer, what is $\det(A)$?
"""

# ╔═╡ 39166607-adb8-43ca-a300-63db130377fc
md"""
$C_{21} = (-1)^{2 + 1} \begin{vmatrix} 0 & -2 \\ 3 & -2 \end{vmatrix} = (-1) (6) = -6$

$C_{33} = (-1)^{3 + 3} \begin{vmatrix} 1 & 0 \\ 5 & -5 \end{vmatrix} = (1) (-5) = -5$

$C = \begin{bmatrix} 10 & 10 & 15 \\ -6 & -2 & -3 \\ -10 & -10 & -5 \end{bmatrix}$

$C^t A = \begin{bmatrix} 10 & -6 & -10 \\ 10 & -2 & -10 \\ 15 & -3 & -5 \end{bmatrix} \begin{bmatrix} 1 & 0 & -2 \\ 5 & -5 & 0 \\ 0 & 3 & -2 \end{bmatrix} = \begin{bmatrix} -20 & 0 & 0 \\ 0 & -20 & 0 \\ 0 & 0 & -20 \end{bmatrix}$

$\det(A) = -20$
"""

# ╔═╡ 4e845087-44d1-41f9-92aa-9e0298015993
let
	A = [1 0 -2
		 5 -5 0
		 0 3 -2]
	Int.(cof(A)), det(A)
end

# ╔═╡ 13c5e648-f3ea-4cde-8d7b-e9f714ae1d2a
md"""
### PS19 #3

Let

$A = \begin{bmatrix} 1 & 0 & -2 & 5 \\ -5 & 0 & 5 & 0 \\ -2 & 1 & 4 & 2 \\ 0 & -3 & -4 & 0 \end{bmatrix}.$

(a) Find the 4,3 minor of $A$.

(b) Except for the number $x$, the cofactor matrix of $A$ is

$C = \begin{bmatrix} -30 & 40 & x & 10 \\ -52 & 48 & -36 & -4 \\ 75 & -100 & 75 & 15 \\ 25 & -60 & 25 & 5 \end{bmatrix}$

With it we have

$C^t A = \begin{bmatrix} \fbox{\, ⠀ \,} & \fbox{\, ⠀ \,} & 0 & 0 \\ 0 & 80 & 0 & 0 \\ x + 30 & 0 & 20 - 2x & 5x + 150 \\ 0 & 0 & 0 & \fbox{\, ⠀ \,} \end{bmatrix}$

Find the missing numbers.


(c) Find $x$ so that the product $C^t A$ is a diagonal matrix.

(d) What is the meaning of the number 80?
"""

# ╔═╡ baf41a50-7f3e-4a07-aef6-561fb5c127d6
md"""
**(a)**

$M_{43} = \begin{vmatrix} 1 & 0 & 5 \\ -5 & 0 & 0 \\ -2 & 1 & 2 \end{vmatrix}$
"""

# ╔═╡ a2c6a412-6e24-4ae8-a652-c060c188f8b4
md"""
**(b)**

$C^t A = \begin{bmatrix} 80 & 0 & 0 & 0 \\ 0 & 80 & 0 & 0 \\ x + 30 & 0 & 20 - 2x & 5x + 150 \\ 0 & 0 & 0 & 80 \end{bmatrix}$
"""

# ╔═╡ f34e13d9-35ef-4d5d-98cb-e383d5ff9bc8
md"""
**(c)**

$x = -30$
"""

# ╔═╡ 2dc1ee27-399b-476a-9a29-dc8b9d0919fd
md"""
**(d)** 80 is the determinant of $A$.
"""

# ╔═╡ 3bc626df-b8e8-463c-ada4-e0cc139891aa
let
	A = [1 0 -2 5
		 -5 0 5 0
		 -2 1 4 2
		 0 -3 -4 0]
	
	round.(Int, cof(A)), det(A)
end

# ╔═╡ 88dbca07-05c4-4759-9cda-4be5c4ac0d30
md"""
### PS19 #4

Let

$A = \begin{bmatrix} 0 & -2 \\ 2 & 5 \end{bmatrix}\;,$

let $I$ be the $2 × 2$ identity matrix, and $\lambda$ be a variable.

(a) Compute $A - \lambda I$. This is a matrix in which some coefficients depend on the variable $\lambda$. Let $p(\lambda) = \det(A - \lambda I)$ (a polynomial in $\lambda$).

(b) $p(\lambda)$ is a polynomial. Find its roots, call them $\lambda_1$ and $\lambda_2$.

(c) Find a nonzero solution of each of the system of equations

$(A - \lambda_i I) x = 0, \quad i = 1, 2,$

where $x = \begin{bmatrix} x_1 \\ x_2 \end{bmatrix}$
"""

# ╔═╡ 8740405a-7a5c-4208-abf6-07dd00b703e4
md"""
**(a)**

$A - \lambda I = \begin{bmatrix} 0 & -2 \\ 2 & 5 \end{bmatrix} - \begin{bmatrix} \lambda & 0 \\ 0 & \lambda \end{bmatrix} = \begin{bmatrix} -\lambda & -2 \\ 2 & 5 - \lambda \end{bmatrix}$

$\begin{align*}
p(\lambda) &= \det(A - \lambda I) \\
&= (-\lambda)(5 - \lambda) - (-2)(2) \\
&= -5 \lambda + \lambda^2 + 4 \\
&= \lambda^2 - 5 \lambda + 4
\end{align*}$
"""

# ╔═╡ 5b52a6e4-38eb-488d-9b9e-c76ec24e13d2
md"""
**(b)**

$\begin{align*}
\lambda &= \frac{5 ± \sqrt{25 - 4(1)(4)}}{2(1)} \\
&= \frac{5 ± \sqrt{25 - 16}}{2} \\
&= \frac{5 ± \sqrt{9}}{2} \\
&= \frac{5 ± 3}{2}
\end{align*}$

$\lambda_1 = 4, \; \lambda_2 = 1$
"""

# ╔═╡ 4692a6e4-5808-40bf-b361-e8202cd97f6f
md"""
**(c)**

$A - \lambda_1 I = A - 4I = \begin{bmatrix} -4 & -2 \\ 2 & 1 \end{bmatrix}$

$\begin{align*}
(A - \lambda_1 I) 𝐱 = 0 &\implies \begin{bmatrix} -4 & -2 & 0 \\ 2 & 1 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & \frac{1}{2} & 0 \\ 0 & 0 & 0 \end{bmatrix} \\
\end{align*}$

$\begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} -\frac{1}{2} x_2 \\ x_2 \end{bmatrix} \implies \vec{x} = \begin{bmatrix} 1 \\ -2 \end{bmatrix} \text{ for } x_2 = -2$

$A - \lambda_2 I = A - I = \begin{bmatrix} -1 & -2 \\ 2 & 4 \end{bmatrix}$

$\begin{align*}
(A - \lambda_2 I) 𝐱 = 0 &\implies \begin{bmatrix} -1 & -2 & 0 \\ 2 & 4 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & 2 & 0 \\ 0 & 0 & 0 \end{bmatrix} \\
\end{align*}$

$\begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} -2 x_2 \\ x_2 \end{bmatrix} \implies \vec{x} = \begin{bmatrix} -2 \\ 1 \end{bmatrix} \text{ for } x_2 = 1$
"""

# ╔═╡ 565df250-b729-40d5-bba9-162307022d22
2.23607*eigvecs([0 -2; 2 5])

# ╔═╡ b3d4322e-001d-41ca-bd18-f66580699506
md"## Problem Set 20"

# ╔═╡ a8252a62-81a0-4739-b140-747714eaa880
md"""
### PS20 #1

Let

$A = \begin{bmatrix} 17 & -8 \\ 30 & -15 \end{bmatrix}$

Find:

(a) $p_A(\lambda)$, the characteristic polynomial of $A$.

(b) The roots of $p_A(\lambda)$, call them $\lambda_1$ and $\lambda_2$.

(c) For each of the two roots $\lambda_i$, a *nonzero* vector $v_i$ such that

$Av_i - \lambda_i v_i = 0.$

(d) Verify: $v_1$ and $v_2$ are linearly independent.
"""

# ╔═╡ 888a0de2-836d-4b9c-963b-8725dcd50ab4
md"""
**(a)**

$A - \lambda I = \begin{bmatrix} 17 - \lambda & -8 \\ 30 & -15 - \lambda \end{bmatrix}$

$\begin{align*}
\det(A - \lambda I) &= (17 - \lambda)(-15 - \lambda) - (-8)(30) \\
&= -255 - 17\lambda + 15\lambda + \lambda^2 + 240 \\
&= \lambda^2 - 2\lambda - 15
\end{align*}$

$p_A(\lambda) = \lambda^2 - 2\lambda - 15$
"""

# ╔═╡ 8bb3d6ae-30fd-4ad5-a1b4-600ddfad2cb6
md"""
**(b)**

$\begin{align*}
p_A(\lambda) &= (\lambda + 3)(\lambda - 5) \\
&\implies \lambda_1 = -3, \; \lambda_2 = 5
\end{align*}$
"""

# ╔═╡ de6a6773-d5e6-429c-9458-7ada3cb3f919
md"""
$\begin{align*}
(A + 3I) 𝐯_1 = 𝟎 &\implies \begin{bmatrix} 20 & -8 & 0 \\ 30 & -12 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 5 & -2 & 0 \\ 5 & -2 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & -\frac{2}{5} & 0 \\ 0 & 0 & 0 \end{bmatrix}
\end{align*}$

$\begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} \frac{2}{5} x_2 \\ x_2 \end{bmatrix} \implies 𝐯_1 = \begin{bmatrix} \frac{2}{5} \\ 1 \end{bmatrix}$

$\begin{align*}
(A - 5I) 𝐯_2 = 𝟎 &\implies \begin{bmatrix} 12 & -8 & 0 \\ 30 & -20 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 3 & -2 & 0 \\ 3 & -2 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & -\frac{2}{3} & 0 \\ 0 & 0 & 0 \end{bmatrix}
\end{align*}$

$\begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} \frac{2}{3} x_2 \\ x_2 \end{bmatrix} \implies 𝐯_2 = \begin{bmatrix} \frac{2}{3} \\ 1 \end{bmatrix}$
"""

# ╔═╡ dd9d5b05-15bc-4795-9d34-a2d23399daab
1/0.928477*eigvecs([17 -8; 30 -15])

# ╔═╡ b285f1a1-4c48-4a6d-a5eb-e013e7c5cae3
1/0.83205*eigvecs([17 -8; 30 -15])

# ╔═╡ 6616a5f6-2b20-442e-bcc2-e34d94e684a2
md"""
**(d)**

$a_1 𝐯_1 + a_2 𝐯_2 = 𝟎$

$a_1 \begin{bmatrix} \frac{2}{5} \\ 1 \end{bmatrix} + a_2 \begin{bmatrix} \frac{2}{3} \\ 1 \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \end{bmatrix}$

$\begin{align*}
\frac{2}{5} a_1 + \frac{2}{3} a_2 &= 0 \\
a_1 + a_2 &= 0
\end{align*}$

$\begin{align*}
\begin{bmatrix} \frac{2}{5} & \frac{2}{3} & 0 \\ 1 & 1 & 0 \end{bmatrix} &⇝ \begin{bmatrix} 1 & \frac{5}{3} & 0 \\ 1 & 1 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & \frac{5}{3} & 0 \\ 0 & -\frac{2}{3} & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{bmatrix} \\
\end{align*}$

$\implies a_1 = a_2 = 0$
"""

# ╔═╡ e4ef72e1-dd20-4537-a5f2-68e37ba0f26e
md"""
### PS20 #2

Do (a)-(d) of the previous problem for

$A = \begin{bmatrix} 1 & 2 & 2 \\ 1 & 1 & 0 \\ 1 & 0 & 1 \end{bmatrix}$

That is, find the eigenvalues of $A$, for each of these find an eigenvector, and verify that these are linearly independent.
"""

# ╔═╡ 0de75b5a-2216-45ac-a03a-3f9554e143cc
md"""
**(a)**

$A - \lambda I = \begin{bmatrix} 1 - \lambda & 2 & 2 \\ 1 & 1 - \lambda & 0 \\ 1 & 0 & 1 - \lambda\end{bmatrix}$

$\begin{align*}
\det(A - \lambda I) &= (1 - \lambda) \begin{vmatrix} 1 - \lambda & 0 \\ 0 & 1 - \lambda \end{vmatrix} - \begin{vmatrix} 2 & 2 \\ 0 & 1 - \lambda \end{vmatrix} + \begin{vmatrix} 2 & 2 \\ 1 - \lambda & 0 \end{vmatrix} \\
&= (1 - \lambda)(1 - \lambda)(1 - \lambda) - (2)(1 - \lambda) - (2)(1 - \lambda) \\
&= (1 - \lambda)(1 - 2 \lambda + \lambda^2) - 4(1 - \lambda) \\
&= (1 - 2 \lambda + \lambda^2 - \lambda + 2\lambda^2 - \lambda^3) - (4 - 4 \lambda) \\
&= -\lambda^3 + 3\lambda^2 - 3\lambda + 1 - 4 + 4 \lambda \\
p_A(\lambda) &= -\lambda^3 + 3\lambda^2 + \lambda - 3 \\
\end{align*}$
"""

# ╔═╡ 7b363010-4b09-4912-b099-0290b127148a
md"""
**(b)**

$\begin{array}{r|rrrr}
3 & -1 & 3 & 1 & -3 \\
&  & -3 & 0 & 3 \\
\hline
& -1 & 0 & 1 & 0
\end{array}$

$\begin{array}{r|rrrr}
1 & -1 & 3 & 1 & -3 \\
&  & -1 & 2 & 3 \\
\hline
& -1 & 2 & 3 & 0
\end{array}$

$\begin{array}{r|rrrr}
-1 & -1 & 3 & 1 & -3 \\
&  & 1 & -4 & 3 \\
\hline
& -1 & 4 & -3 & 0
\end{array}$

$\lambda_1 = -1,\; \lambda_2 = 1,\; \lambda_3 = 3$
"""

# ╔═╡ b84c139f-4919-4041-b742-ac41992d62b1
md"""
**(c)**

$\begin{align*}
(A - \lambda_1 I)𝐯_1 = (A + I)𝐯_1 = 𝟎 &\implies \begin{bmatrix} 2 & 2 & 2 & 0 \\ 1 & 2 & 0 & 0 \\ 1 & 0 & 2 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & 1 & 1 & 0 \\ 0 & 1 & -1 & 0 \\ 0 & 1 & -1 & 0\end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & 0 & 2 & 0 \\ 0 & 1 & -1 & 0 \\ 0 & 0 & 0 & 0 \end{bmatrix}
\end{align*}$

$𝐯_1 = \begin{bmatrix} -2 \\ 1 \\ 1 \end{bmatrix}$
"""

# ╔═╡ ccc2e160-bbe6-4bdc-83ef-37d217023059
md"""
$\begin{align*}
(A - \lambda_2 I)𝐯_1 = (A - I)𝐯_2 = 𝟎 &\implies \begin{bmatrix} 0 & 2 & 2 & 0 \\ 1 & 0 & 0 & 0 \\ 1 & 0 & 0 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & 0 & 0 & 0 \\ 0 & 2 & 2 & 0 \\ 0 & 0 & 0 & 0\end{bmatrix}
\end{align*}$

$𝐯_2 = \begin{bmatrix} 0 \\ -1 \\ 1 \end{bmatrix}$
"""

# ╔═╡ 85d46495-80a9-4715-a9ff-f15fdbb7d7f5
md"""
$\begin{align*}
(A - \lambda_3 I)𝐯_1 = (A - 3I)𝐯_2 = 𝟎 &\implies \begin{bmatrix} -2 & 2 & 2 & 0 \\ 1 & -2 & 0 & 0 \\ 1 & 0 & -2 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & -1 & -1 & 0 \\ 0 & -1 & 1 & 0 \\ 0 & 1 & -1 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & 0 & -2 & 0 \\ 0 & 1 & -1 & 0 \\ 0 & 0 & 0 & 0 \end{bmatrix}
\end{align*}$

$𝐯_3 = \begin{bmatrix} 2 \\ 1 \\ 1 \end{bmatrix}$
"""

# ╔═╡ 80d68dd0-d5e8-4fb0-bada-67a34f654036
md"""
**(d)**

$\begin{align*}
a_1 𝐯_1 + a_2 𝐯_2 + a_3 𝐯_3 = 𝟎 &\implies \begin{bmatrix} -2 & 0 & 2 & 0 \\ 1 & -1 & 1 & 0 \\ 1 & 1 & 1 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & 0 & -1 & 0 \\ 1 & -1 & 1 & 0 \\ 0 & 1 & 0 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & 0 & -1 & 0 \\ 0 & 1 & -2 & 0 \\ 0 & 0 & 1 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \end{bmatrix} \\
\end{align*}$
"""

# ╔═╡ 5f3d4b7e-78da-4111-b929-ddeca6e9e774
md"""
### PS20 #3

$A = \begin{bmatrix} 1 & -2 & 2 \\ 1 & 1 & 0 \\ -1 & 0 & 1 \end{bmatrix}.$

(a) Let $p_A(\lambda)$ be the characteristic polynomial of $A$.
One of its roots is 1.
Find the other roots.

(b) Find an eigenvector corrseponding to the eigenvalue 1.
"""

# ╔═╡ 907d5b75-0fe3-4a98-a528-e4c72323584d
md"""
**(a)**

$A - \lambda I = \begin{bmatrix} 1 - \lambda & -2 & 2 \\ 1 & 1 - \lambda & 0 \\ -1 & 0 & 1 - \lambda \end{bmatrix}$

$\begin{align*}
\det(A - \lambda I) &= -1 \begin{bmatrix} -2 & 2 \\ 1 - \lambda & 0\end{bmatrix} + (1 - \lambda) \begin{vmatrix} 1 - \lambda & -2 \\ 1 & 1 - \lambda \end{vmatrix} \\
&= 2 - 2 \lambda + (1 - \lambda)((1 - \lambda)^2 + 2) \\
&= 2 - 2 \lambda + (1 - \lambda)(3 - 2 \lambda + \lambda^2) \\
&= 2 - 2 \lambda + 3 - 2 \lambda + \lambda^2 - 3 \lambda + 2 \lambda^2 - \lambda^3 \\
&= -\lambda^3 + 3 \lambda^2 - 7 \lambda + 5 \\
&= -(\lambda^3 - 3\lambda^2 + 7 \lambda - 5) \\
&= -(\lambda - 1)(\lambda^2 + 2\lambda + 5)
\end{align*}$

$\begin{array}{c|cccc}
1 & 1 & -3 & 7 & -5 \\
& & 1 & -2 & -5 \\
\hline
& 1 & -2 & 5 & 0
\end{array}$

$\lambda = \frac{-2 ± \sqrt{4 - 4(5)}}{2} = \frac{-2 ± \sqrt{-16}}{2} = \frac{-2 ± 4i}{2} = -1 ± 2i$

$(A - I)𝐱 = 𝟎 \implies \begin{bmatrix} 0 & -2 & 2 \\ 1 & 0 & 0 \\ -1 & 0 & 0 \end{bmatrix} ⇝ \begin{bmatrix} 1 & -2 & 2 \\ 0 & 2 & -2 \\ 0 & 0 & 0 \end{bmatrix} ⇝ \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & -1 \\ 0 & 0 & 0 \end{bmatrix}$

$𝐱 = \begin{bmatrix} 0 \\ 1 \\ 1 \end{bmatrix}$
"""

# ╔═╡ 319b0a44-59ed-475a-ae6f-ac5a68b1a7ed
md"""
### PS20 #4

$A = \begin{bmatrix} 3 & 10 & 10 \\ 0 & -3 & -6 \\ 0 & -2 & 1 \end{bmatrix}$
"""

# ╔═╡ 5d0be0ee-8987-44d5-aa64-ff434a91acc4
md"""
**(a)**

$\begin{align*}
\det(A - \lambda) &= (3 - \lambda) \begin{vmatrix} -3 - \lambda & -6 \\ -2 & 1 - \lambda \end{vmatrix} \\
&= (3 - \lambda)((-3 - \lambda)(1 - \lambda) - 12) \\
&= (3 - \lambda)(-3 + 3 \lambda - \lambda + \lambda^2 - 12) \\
&= (3 - \lambda)(\lambda^2 + 2 \lambda - 15) \\
&= 3\lambda^2 + 6\lambda - 45 - \lambda^3 - 2\lambda^2 + 15\lambda \\
&= -\lambda^3 + \lambda^2 + 21\lambda - 45 \\
&= -(\lambda^3 - \lambda^2 - 21\lambda + 45) \\
&= -(\lambda + 5)(\lambda^2 - 6\lambda + 9) \\
&= -(\lambda + 5)(\lambda - 3)^2
\end{align*}$

$\begin{array}{c|cccc}
-5 & 1 & -1 & -21 & 45 \\
& & -5 & 30 & -45 \\
\hline
& 1 & -6 & 9 & 0
\end{array}$

"""

# ╔═╡ f8998419-f2db-49a6-975f-873896a3197a
md"""
$\text{Basis of } E_0 = \left\{\begin{bmatrix} 0 \\ -1 \\ 1 \end{bmatrix}, \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix}\right\}$
"""

# ╔═╡ 6eb690f1-664f-499d-85b9-e5078860464c
md"## Problem Set 21"

# ╔═╡ a43beff5-a615-4c88-9265-742c41ffa624
md"""
### PS21 #1

Let $Q$ be an $n × n$ matrix and $𝐯_1$, $𝐯_2$ be two nonzero vectors in $ℝ^n$ with the property that

$Q𝐯_1 = 0, \quad Q𝐯_2 ≠ 0$

Verify that $𝐯_1$ and $𝐯_2$ are linearly independent.
"""

# ╔═╡ da5aabe1-d5d1-4f93-8b6d-f9b31b0d8d89
md"""
**Solution:**
Suppose

$\tag{*} x_1 𝐯_1 + x_2 𝐯_2 = 0$

for some numbers $x_1, x_2$.
Then $Q(x_1 𝐯_1 + x_2 𝐯_2) = Q(0)$.
The left hand side of this equality is $x_1 Q 𝐯_1 + x_2 Q 𝐯_2$, the right hand side is 0:

$x_1 Q𝐯_1 + x_2 Q 𝐯_2 = 0$

Since $Q𝐯_1 = 0$, this reduces to $x_2 Q𝐯_2 = 0$.
Since $Q𝐯_2 ≠ 0$, it must be that $x_2 = 0$.
Replacing $x_2 = 0$ in $(*)$ we get $x_1 𝐯_1 = 0$.
Since $𝐯_1 ≠ 0$, $x_1 = 0$.
So, in $(*)$ $x_1$ and $x_2$ *must* both be zero.
This tells us $𝐯_1$ and $𝐯_2$ are linearly independent.
"""

# ╔═╡ 8956e768-4de6-4f69-b8f0-8e8396728e61
md"""
### PS21 #2

Suppose $A$ is an $n × n$ matrix and $𝐯_1$ is an eigenvector of $A$ with eigenvalue $\lambda_1$.
Suppose $𝐯_2$ is some nonzero vector that is not an eigenvector for the eigenvalue $\lambda_1$.
Verify that $𝐯_1$ and $𝐯_2$ are linearly independent.
"""

# ╔═╡ d3636101-679f-4895-92d4-4957585f9cee
md"""
**Solution:**

$\begin{align*}
(A - \lambda_1) 𝐯_1 &= 𝟎, \quad 𝐯_1 ≠ 𝟎 \\
(A - \lambda_2) 𝐯_2 &= 𝟎, \quad 𝐯_2 ≠ 𝟎
\end{align*}$

$(A - \lambda_1 I)(A - \lambda_2 I) = \begin{bmatrix} b_{11} & b_{12} \\ b_{21} & b_{22} \end{bmatrix}$

$(A - \lambda_1 I)(A - \lambda_2 I) \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} b_{11} x_1 & b_{12} x_2 \\ b_{21} x_1 & b_{22} x_3 \end{bmatrix}$

$\begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = y_1 𝐯_1 + y_2 𝐯_2$

$\begin{align*}
&(A - \lambda_1 I)(A - \lambda_2 I)(y_1 𝐯_1 + y_2 𝐯_2) \\
&= y_1(A - \lambda_1 I)(A - \lambda_2 I) 𝐯_1 + y_2(A - \lambda_1 I)(A - \lambda_2 I) 𝐯_2 \\
&= 𝟎
\end{align*}$

General fact:

$p_A(\lambda) = \lambda^n + c_{n - 1} \lambda^{n - 1} + ⋯ + c_1 \lambda^1 + c_0$

$p_A(A) = A^n + c_{n - 1} A^{n - 1} + ⋯ + c_1 A + c_0 I$
"""

# ╔═╡ 3584afba-39bc-4a15-afcc-0c53b22ae2ae
md"""
### PS21 #3

The characteristic polynomial of

$A = \begin{bmatrix} 5 & -3 & 3 \\ 8 & -6 & 3 \\ 8 & -8 & 5 \end{bmatrix}$

is $p_A(\lambda) = (\lambda - 2)(\lambda - 5)(\lambda + 3)$.
"""

# ╔═╡ 4eae4bdd-7103-44d3-8fda-4c3d8a665781
let
	A = [5 -3 3
		 8 -6 3
		 8 -8 5]
	
	λ1, λ2, λ3 = 2, 5, -3
	
	[round.(Int, rref(A - λ*I)) for λ ∈ (λ1, λ2, λ3)]
end

# ╔═╡ 6084e51c-caef-472a-a285-5aae4b2f11f1
md"""
$𝐯_1 = \begin{bmatrix} 1 \\ 1 \\ 0 \end{bmatrix} \quad 𝐯_2 = \begin{bmatrix} 1 \\ 1 \\ 1 \end{bmatrix} \quad 𝐯_3 = \begin{bmatrix} 0 \\ 1 \\ 1 \end{bmatrix}$
"""

# ╔═╡ 675b6c6f-7a0c-467a-99b3-347d810c9074
let
	round.(Int, rref([1 1 0; 1 1 1; 0 1 1]))
end

# ╔═╡ 923912b7-322f-4b2e-a3ec-46265f1c61c9
let
	A = [5 -3 3
		 8 -6 3
		 8 -8 5]
	
	T = [1 1 0; 1 1 1; 0 1 1]
	
	A * T
end

# ╔═╡ 758ea353-65ca-483a-9e44-d08e3789823b
md"## Problem Set 22"

# ╔═╡ 1bacfb53-3650-4272-9e2b-51f2fbf02c1a
md"""
### PS22 #1

The characteristic polynomial of

$A = \begin{bmatrix} -12 & 3 & 5 \\ -29 & 8 & 11 \\ -13 & 3 & 6 \end{bmatrix}$

has roots $-1, 1, 2$.
Find a matrix $P$ such that

$P^{-1}AP = \begin{bmatrix} -1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 2 \end{bmatrix}$
"""

# ╔═╡ 1c19cfdc-52c8-41a4-8347-eeabf853303e
md"""
$𝐯_1 = \begin{bmatrix} 1 \\ 2 \\ 1 \end{bmatrix} \quad 𝐯_2 = \begin{bmatrix} \frac{1}{2} \\ \frac{1}{2} \\ 1 \end{bmatrix} \quad 𝐯_3 = \begin{bmatrix} 1 \\ 3 \\ 1 \end{bmatrix}$

$P = \begin{bmatrix} 1 & \frac{1}{2} & 1 \\ 2 & \frac{1}{2} & 3 \\ 1 & 1 & 1 \end{bmatrix}$
"""

# ╔═╡ b58089bb-b178-4d5c-99cd-ea5fc254edf2
let
	A = [-12 3 5
		 -29 8 11
		 -13 3 6]
	
	(λ1, λ2, λ3) = (-1, 1, 2)
	[rref(M) for M ∈ [A - λ1*I, A - λ2*I, A - λ3*I]]
end

# ╔═╡ 72e3b4bb-5a09-44fa-84cb-3e442d599763
let
	A = [-12 3 5
		 -29 8 11
		 -13 3 6]
	
	P = [1 0.5 1
		 2 0.5 3
		 1 1 1]
	
	round.(Int, inv(P) * A * P)
end

# ╔═╡ 61608c42-ab66-4d31-85a0-5ad9cff5d904
md"""
### PS22 #2

Let $A$ be some matrix, let $b_0, b_1, b_2, c_1, c_2, c_3$ be some numbers.
Let

$P^{-1}B = b_0 I + b_1 A + b_2 A^2, \quad C = c_0 I + c_1 A + c_2 A^2.$

Verify that $BC = CB$.
"""

# ╔═╡ f94be157-8409-42e6-bd2a-4076e2e98f4f
md"""
### PS22 #3

Let

$N = \begin{bmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{bmatrix}$

Find $A^2$ and $A^3$.
"""

# ╔═╡ 98e27c9b-850c-41be-aac3-565bc8475687
md"""
$\det(\lambda I - N) = \begin{bmatrix} \lambda & -1 & 0 \\ 0 & \lambda & -1 \\ 0 & 0 & \lambda \end{bmatrix} = \lambda^3$

$N^2 = \begin{bmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{bmatrix} \begin{bmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{bmatrix} = \begin{bmatrix} 0 & 0 & 1 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{bmatrix}$

$N^3 = \begin{bmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{bmatrix} \begin{bmatrix} 0 & 0 & 1 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{bmatrix} = \begin{bmatrix} 0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{bmatrix}$

$N^k = 𝟎^3 \text{ for } k ≥ 3$
"""

# ╔═╡ 71a1cdbd-f12f-4ab6-83f2-b495dedaa953
md"## Problem Set 23"

# ╔═╡ 3618a2f0-c59c-4e07-8a1b-2d6f493612d6
md"""
### PS23 #1

Let

$A = \begin{bmatrix} 1 & 4 \\ -1 & 5 \end{bmatrix}.$

Find:

(a) The eigenvalue of $A$ (it only has one), call it $\lambda_1$.

(b) Verify that $(A - \lambda_1 I)^2 = 0$.

(c) Find an eigenvector of $A$, call it $v_1$, and a solution of $(A - \lambda_1 I) v_2 = v_1$.

(d) Verify that $v_1$ and $v_2$ are linearly independent.

(e) Let $T : ℝ^2 → ℝ^2$ be the linear transformation defined by $T(x) = Ax$.
Find the matrix of $T$ with respect to the basis $\{v_1, v_2\}$.
"""

# ╔═╡ 87e01065-9751-47ee-8144-d7aceb48cb9f
md"""
$\det(A - \lambda_1 I) = (1 - \lambda_1)(5 - \lambda_1) + 4 = 9 - 6\lambda_1 + {\lambda_1}^2$

${\lambda_1}^2 - 6\lambda_1 + 9 = (\lambda_1 - 3)^2 = 0 \implies \lambda_1 = 3$
"""

# ╔═╡ d862a4c7-1bd3-4d2f-a1c7-717904f8061b
md"""
## Problem Set 24

Let $A$ be some $n × n$ matrix, $\lambda_0$ one of its eigenvalues.
A generalized eigenvector of $A$ corresponding to $\lambda_0$ is a nonzero vector such that $(A - \lambda_0 I)^k 𝐯 = 𝟎$ for some $k$.
If the multiplicity of the root $\lambda_0$ in the characteristic polynomial of $A$ is $m$, then all generalized eigenvectors will appear as solutions of $(A - \lambda_0 I)^k 𝐱 = 𝟎$.
The generalized eigenspace of $A$ corresponding to that eigenvalue is $\text{Null}((A - \lambda_0 I)^m)$.
"""

# ╔═╡ f50a0569-4999-4ce1-9664-f9d9b6a26253
md"""
### PS24 #1

Let

$A = \begin{bmatrix} -2 & 0 & 6 \\ -6 & 4 & 7 \\ 0 & 0 & 4 \end{bmatrix}.$

Find the generalized eigenspaces of $A$ (write a basis for each of these spaces).
"""

# ╔═╡ c84211c9-a066-4120-9dd3-042af16f7044
md"""
$A - \lambda I = \begin{bmatrix} -2 - \lambda & 0 & 6 \\ -6 & 4 - \lambda & 7 \\ 0 & 0 & 4 - \lambda \end{bmatrix}$

$\begin{align*}
\det(A - \lambda) &= (4 - \lambda) \begin{vmatrix} -2 - \lambda & 0 \\ -6 & 4 - \lambda \end{vmatrix} \\
&= (4 - \lambda) (-2 - \lambda) (4 - \lambda) \\
&= (-2 - \lambda) (4 - \lambda)^2 \\
&\implies \lambda_1 = -2,\; \lambda_2 = 4
\end{align*}$
"""

# ╔═╡ c998e12b-dc4b-4e6e-a457-3b543dbc45cd
md"""
$\begin{align*}
(A + 2I) 𝐯_1 = 𝟎 &\implies \begin{bmatrix} 0 & 0 & 6 & 0 \\ -6 & 6 & 7 & 0 \\ 0 & 0 & 6 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 0 & 0 & 1 & 0 \\ 1 & -1 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{bmatrix} \\
&\implies 𝐯_1 = \begin{bmatrix} 1 \\ 1 \\ 0 \end{bmatrix}
\end{align*}$
"""

# ╔═╡ 377a28c9-3245-4c2b-b026-ac631fd06528
md"""
$\begin{align*}
(A - 4I)^2 𝐯_2 = 𝟎 &\implies \begin{bmatrix} 36 & 0 & -36 & 0 \\ 36 & 0 & -36 & 0 \\ 0 & 0 & 0 & 0 \end{bmatrix} \\
&⇝ \begin{bmatrix} 1 & 0 & -1 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{bmatrix} \\
&\implies \text{basis of } \text{Null}((A - 4I)^2) = \left\{\begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix}, \begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix}\right\}
\end{align*}$
"""

# ╔═╡ 08903063-f791-4533-a65e-74e6229a96d2
let
	A = [-2 0 6; -6 4 7; 0 0 4]
	nullspace(A + 2I), nullspace((A - 4I)^2)
end

# ╔═╡ a127a412-8301-4e9d-bf1c-0ec403bdaf21
md"""
### PS24 #2

Redo Problem 1 with

$\begin{bmatrix} 4 & 1 & 1 \\ -1 & 6 & 2 \\ 0 & 0 & 5 \end{bmatrix}.$
"""

# ╔═╡ 21d1467f-8723-4b2c-b029-150751d247ce
md"""
$\begin{align*}
\det(A - \lambda I) &= (5 - \lambda) \begin{vmatrix} 4 - \lambda & 1 \\ -1 & 6 - \lambda \end{vmatrix} \\
&= (5 - \lambda) \left[(4 - \lambda)(6 - \lambda) + 1\right] \\
&= (5 - \lambda) \left(25 - 10 \lambda + \lambda^2\right) \\
&= (5 - \lambda)^3 \implies \lambda = 5
\end{align*}$

$\begin{align*}
(A - 5I)^3 𝐯 = 𝟎 &\implies \begin{bmatrix} 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{bmatrix} \\
&\implies \text{Null}((A - 5I)^3) = \left\{\begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix},\begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix},\begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix}\right\}
\end{align*}$
"""

# ╔═╡ c765c755-cdae-4a18-88b2-a174bb1893ba
let
	A = [4 1 1; -1 6 2; 0 0 5]
	nullspace((A - 5I)^3)
end

# ╔═╡ 15d25922-8ae5-4ab4-a842-90aa1efab208
md"""
### PS24 #3

Let $A$ be some $n × n$ matrix, let $\lambda_0$ be a root of its characteristic polynomial, of multiplicity $m$.
Let $V$ be the generalized eigenspace of $A$ for the eigenvalue $\lambda_0$.
Verify the following:

$\text{If } v ∈ V, \text{ then } Av ∈ V$
"""

# ╔═╡ 834c91c6-08e6-409b-8de7-6dd1ec5c6d00
md"""
If $(A - \lambda_0 Id)v = 0$, then $(A - \lambda_0 I)(Av) = 0$.
"""

# ╔═╡ 71585849-3a8e-4f26-a476-589969fd0999
md"""
### PS24 #4

Newton's binomial formula for $(a + b)^k$ is $a^2 + 2ab + b^2$ if $k = 2$, $a^3 + 3a^2 b + 3ab^2 + b^3$ for $k = 3$ and so on.
Suppose $A$ and $B$ are square matrices:

a) Without assuming that $AB = BA$, expand:
  1. ``(A + B)^2``
  2. ``(A + B)^3``

b) Now assume that $AB = BA$. Expand and simplify:

  1. ``(A + B)^2``
  2. ``(A + B)^3``
"""

# ╔═╡ 5e011536-a63c-480a-a80f-e5154d2e9c78
md"""
**(a)**

$(A + B)^2 = (A + B)(A + B) = A^2 + AB + BA + B^2$

$(A + B)^3 = (A + B)(A^2 + AB + BA + B^2)$
$= A^3 + A^2 B + ABA + AB^2 + BA^2 + BAB + B^2 A + B^3$
"""

# ╔═╡ 13be3db5-3a54-43cf-bfac-7fc6ebba918d
md"""
**(b)**

$(A + B)^2 = (A + B)(A + B) = A^2 + 2AB + B^2$

$(A + B)^3 = (A + B)(A^2 + AB + BA + B^2) = A^3 + 3A^2 B + 3AB^2 + B^3$
"""

# ╔═╡ c2a6a003-4e2f-43d8-addf-1f3d6ba0a92f
md"""
### PS24 #5

Let

$N = \begin{bmatrix} 0 & 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 0 & 1 \\ 0 & 0 & 0 & 0 & 0 \end{bmatrix}.$

This matrix has the property that $N^5 = 0$.
Let $P = I - N$, let $Q = I + N + N^2 + N^3 + N^4$.
Verify that

$PQ = I, \quad QP = I.$
"""

# ╔═╡ a8eb7f21-aab2-4e47-8d27-cc5257a6b4fd
let
	N = [0 1 0 0 0
	     0 0 1 0 0
	     0 0 0 1 0
		 0 0 0 0 1
		 0 0 0 0 0]

	P = I - N

	Q = I + N + N^2 + N^3 + N^4

	P*Q, Q*P
end

# ╔═╡ b4944c17-c72d-421c-9009-0ad3fdec10db
md"## Problem Set 25"

# ╔═╡ 5b5e5890-bf08-4ab7-bf13-ac3aa27adce9
md"""
### PS25 #1

Let $V$ be a 2-dimensional vector space with inner product $⟨,⟩$.
Suppose that with a certain basis $u_1$, $u_2$ of $V$ we have

$⟨u_1, u_1⟩ = 2, \quad ⟨u_1, u_2⟩ = -1, \quad ⟨u_2, u_2⟩ = 4.$

Find:

(a) $⟨u_2, u_1⟩.$

(b) Use the values of the various inner products $⟨u_i, u_j⟩$ to construct the matrix

$B = \begin{bmatrix} ⟨u_1, u_1⟩ & ⟨u_1, u_2⟩ \\ ⟨u_2, u_1⟩ & ⟨u_2, u_2⟩ \end{bmatrix}$

Find the eigenvalues of $B$.

(c) For each eigenvalue, find an eigenvector.
These are vectors in $ℝ^2$; write them as

$a = \begin{bmatrix} a_1 \\ a_2 \end{bmatrix}, \quad b = \begin{bmatrix} b_1 \\ b_2 \end{bmatrix}$

Find $a ⋅ b$ (this is using the standard dot product of $ℝ^2$)

(d) Let $\alpha = a_1 u_1 + a_2 u_2$, let $\beta = b_1 u_1 + b_2 u_2$.
Find $⟨\alpha, \beta⟩$.

(e) Find $⟨\alpha, \alpha⟩$ and $⟨\beta, \beta⟩$
"""

# ╔═╡ 22ae14fb-1e25-4bcf-a49c-2bd9bab08a9d
let
	u1 = [0.5, -1.322877]
	u2 = [-2, 0]
 	dot(u2, u1)
end

# ╔═╡ b67c846e-acb6-4a75-8051-d722d88796cf
md"(a) $⟨u_2, u_1⟩ = -1$ from symmetry of inner products"

# ╔═╡ 29554616-3b16-4d85-b910-f8875184d98c
md"$B = \begin{bmatrix} 2 & -1 \\ -1 & 4 \end{bmatrix}$"

# ╔═╡ 69407621-152a-46b9-8259-898347efcc2f
let
	B = [ 2 -1
	     -1  4]
	eigen(B), dot(eigvecs(B)[:,1], eigvecs(B)[:,2])
end

# ╔═╡ 390d2d81-16fc-4f19-bbb3-c3623c1fe8e0
md"## Problem Set 27"

# ╔═╡ f27c2ade-291b-45b7-83e5-723a157d63f3
md"""
### PS27 #1

Let $A$ be some $3 × 3$ matrix, write $A_1, A_2, A_3$ for its *columns*.
Let

$C = \begin{bmatrix} 2 & 0 & 2 \\ 2 & 1 & 3 \\ 6 & 1 & 8 \end{bmatrix}$

Find $B = AC$.
"""

# ╔═╡ 27a55e85-10ec-4766-86b2-07d75541c55c
md"""
**Solution.**

$B = \begin{bmatrix} 2A_1 + 2A_2 + 6A_3 & A_2 + A_3 & 2A_1 + 3A_2 + 8A_3 \end{bmatrix}$
"""

# ╔═╡ 79262f6a-607c-4dfc-92ec-464ac5905a94
md"""
### PS27 #2

Let $A$ be some $3 × 3$ matrix, write $A_1, A_2, A_3$ for its *rows*.
Let

$Q = \begin{bmatrix} 1 & -2 & 0 \\ 0 & -2 & 1 \\ 1 & 0 & -2 \end{bmatrix}$

Find $B = QA$.
"""

# ╔═╡ c5f1d778-bde7-4856-8cff-5929840b49f1
md"""
**Solution.**

$B = \begin{bmatrix} A_1 - 2A_2 \\ -2A_2 + A_3 \\ A_1 - 2A_3 \end{bmatrix}$
"""

# ╔═╡ f02e8a80-84af-4749-8ba6-ab0b379821c6
md"""
### PS27 #3

Let $A$ be a $3 × 3$ matrix with columns $A_1, A_2, A_3$.
Let $B$ be the matrix with columns

$\begin{align*}
B_1 &= 3A_1 + A_2 + A_3 \\
B_2 &= -2A_1 - A_3, \\
B_3 &= 2A_1 + 2A_2 + 3A_3
\end{align*}$

Assuming $\det(A) = 7$, find $\det(B)$.
"""

# ╔═╡ 4e4cf97c-c98c-475e-9748-d68b2c67b808
let
	B1 = [3 1 1]
	B2 = [-2 0 -1]
	B3 = [2 2 3]
	C = [B1' B2' B3']
	7det(C)
end

# ╔═╡ 9bdcee0d-08a4-41b8-b286-5997907d1d05
md"""
### PS27 #4

Let $P$ be the space of all polynomials of degree $≤ 2$ in the variable $t$ with real coefficients:

$V = \{p_0 t^0 + p_1 t^1 + p_2 t^2 : p_0, p_1, p_2 \text{ are real numbers}\}$

Define $T : V → V$ by the formula

$T(p) = tp' + p'$

Find:

(a) The matrix of $T$ with respect to the basis $t^0, t^1, t^2$ of $V$.
Let $A$ denote this matrix.

(b) The eigenvalues of $A$ (all have multiplicity 1).

(c) An eigenvector for each of the eigenvalues.
Call these

$\alpha = \begin{bmatrix} a_0 \\ a_1 \\ a_2 \end{bmatrix}, \quad \beta = \begin{bmatrix} b_0 \\ b_1 \\ b_2 \end{bmatrix}, \quad \gamma = \begin{bmatrix} c_0 \\ c_1 \\ c_2 \end{bmatrix},$

Let $v_\alpha = a_0 t^0 + a_1 t^1 + a_2 t^2$, let $v_\beta$ and $v_\gamma$ be defined similarly using $\beta$ and $\gamma$.
Verify that these three elements of $V$ are linearly independent.

(d) Find the matrix of $T$ with respect to the basis $v_\alpha, v_\beta, v_\gamma$ of $V$.
"""

# ╔═╡ 0da3b8b5-b01d-417f-b462-d4372b91327b
md"""
**Solution.**

$\begin{align*}
T(t^0) &= 0t^0 + 0t^1 + 0t^2 \\
T(t^1) &= t^0 + 1t^1 + 0t^2 \\
T(t^2) &= 0t^0 + 2t^1 + 2t^2
\end{align*}$

$A = \begin{bmatrix} 0 & 1 & 0 \\ 0 & 1 & 2 \\ 0 & 0 & 2 \end{bmatrix}$

$\lambda I - A = \begin{bmatrix} \lambda & 1 & 0 \\ 0 & \lambda - 1 & 2 \\ 0 & 0 & \lambda - 2 \end{bmatrix}$

$\det(\lambda I - A) = \lambda (\lambda - 1)(\lambda - 2) \implies \lambda_1 = 0, \; \lambda_2 = 1, \; \lambda_3 = 2$

$\alpha = \begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix}, \quad \beta = \begin{bmatrix} 1 \\ 1 \\ 0 \end{bmatrix}, \quad \gamma = \begin{bmatrix} 1 \\ 2 \\ 1 \end{bmatrix}$
"""

# ╔═╡ d11d84b9-5a3c-42b4-a849-c74c4d30e46d
rref([1 1 1; 0 1 2; 1 0 1])

# ╔═╡ 5591399b-db5e-45fd-834a-f3f92b89eb03
md"""
### PS27 #5

Suppose $A$ is an $n × n$ matrix and that its characteristic polynomial has roots $\lambda_1, …, \lambda_k$ (here the $\lambda_j$ are all different).
Let the $\lambda_j$ all have multiplicity $m_j = 1$ except $\lambda_k$, which has multiplicity $m_k$.
Find out what $m_k$ must be, and also the dimension of the generalized eigenspace of $A$ corresponding to the eigenvalue $\lambda_k$.
"""

# ╔═╡ d2ccaef6-4b98-407d-989a-9002232e2d57
md"""
**Solution.**

$(A - \lambda_1)^{m_1} ⋯ (A - \lambda_k)^{m_k}$

$(A - \lambda_1) ⋯ (A - \lambda_{k - 1}) (A - \lambda_k)^{m_k}$

$n = k - 1 + m_k$
"""

# ╔═╡ 0e02b91d-43cd-49f6-9e7e-19f583481331
md"""
### PS27 #13

Let $A$ be a $3 × 3$ matrix with characteristic polynomial $p_A(\lambda) = (\lambda + 2)^2 (\lambda - 5)$. Let $v_1, v_2$ be a basis of the space of $\text{Null}(A + 2I)^2$ (the generalized eigenspace of $A$ for the eigenvalue $-2$), let $v_3$ be an eigenvector for the eigenvalue 5.

(a) Verify that $(A + 2I)^2 v_3 ≠ 0$.

(b) Use (a) to verify that $v_1, v_2, v_3$ are linearly independent.

(c) Give an argument as to why $(A + 2I)^2 (A - 5I) v = 0$.
"""

# ╔═╡ b06ef30d-0ec4-4e38-a931-19bf3dc093e3
md"""
**Solution.**

$(A + 2I)^2 (A - 5I)𝐯 = 𝟎?$

$(A + 2I)^2 (A - 5I)𝐯_2 = 𝟎$

$(A + 2I)^2 (A - 5I) 𝐯_2 = (A - 5I) (A + 2I)^2 𝐯_2 = 𝟎$

$𝐯_2 = \text{Null}(A + 2I)^2$
"""

# ╔═╡ e5c86d2a-5133-4b88-9818-62748974e749
let
	B = [2 -1; -1 4]
	(B - 7I)
end

# ╔═╡ 50eb9439-c478-4ece-85e3-2e49afc5c2be
md"## Problem Set 28"

# ╔═╡ c8d7cb79-3466-49e4-8117-6eaed1a4902f
md"""
### PS28 #1

Let $V$ be the vector space consisting of all polynomials of degree at most 2.
Let $v_0 = t^0, v_1 = t^1, v_2 = t^2$, and let

$w_0 = t^0 + t^2, \quad w_1 = -2t^0 - 2t^1, \quad w_2 = t^1 - 2t^2$

Verify:

(a) $[w_0, w_1, w_2]$ is obtained from $[v_0, v_1, v_2]$ by multiplying (and placing the coefficients on the left)

$[v_0, v_1, v_2] Q$

where

$Q = \begin{bmatrix} 1 & -2 & 0 \\ 0 & -2 & 1 \\ 1 & 0 & -2 \end{bmatrix}$

(b) Find

$[w_0, w_1, w_2] Q^{-1}.$

(c) Let $x_0, x_1, x_2$ denote unknowns.
Write $x_0 (t^0 + t^2) + x_1 (-2t^0 - 2t^1) + x_2 (t^1 - 2t^2)$ as a polynomial in $t$, let $c_0, c_1, c_2$ the coefficients (which are linear functions of $x_0, x_1, x_2$).
Find the matrix of the homogeneous system

$c_0 = 0, \quad c_1 = 0, \quad c_2 = 0.$

(d) Use $(*)$ directly to write formulas for $t^0, t^1, t^2$ in terms of $w_0, w_1, w_2$.
How are these formulas related to $Q^{-1}$?
"""

# ╔═╡ 5b31ed0e-bbc9-4901-9be7-10d38c9c8dc0
md"""
**(a)**

$\begin{align*}
[v_0, v_1, v_2] Q &= [t^0, t^1, t^2] \begin{bmatrix} 1 & -2 & 0 \\ 0 & -2 & 1 \\ 1 & 0 & -2 \end{bmatrix} \\
&= [t^0 + t^2, -2t^0 - 2t^1, t^1 - 2t^2] \\
&= [w_0, w_1, w_2]
\end{align*}$
"""

# ╔═╡ 75065d30-85eb-4ad9-9916-9c5ca36f9fc3
md"""
**(b)**

$\begin{align*}
[w_0, w_1, w_2] Q^{-1} &= [w_0, w_1, w_2] \begin{bmatrix} 2 & -2 & -1 \\ \frac{1}{2} & -1 & -\frac{1}{2} \\ 1 & -1 & -1 \end{bmatrix} \\
&= \begin{bmatrix} 2w_0 + \frac{1}{2} w_1 + w_2, -2w_0 - w_1 - w_2, -w_0 - \frac{1}{2} w_1 - w_2 \end{bmatrix}
\end{align*}$

$\begin{align*}
2w_0 + \frac{1}{2} w_1 + w_2 &= 2(t^0 + t^2) + \frac{1}{2} (-2t^0 - 2t^1) + (t^1 - 2t^2) \\
&= 2t^0 + 2t^2 - t^0 - t^1 + t^1 - 2t^2 \\
&= t^0
\end{align*}$

$\begin{align*}
-2w_0 - w_1 - w_2 &= -2(t^0 + t^2) - (-2t^0 - 2t^1) - (t^1 - 2t^2) \\
&= -2t^0 - 2t^2 + 2t^0 + 2t^1 - t^1 + 2t^2 \\
&= t^1
\end{align*}$

$\begin{align*}
-w_0 - \frac{1}{2} w_1 - w_2 &= -(t^0 + t^2) - \frac{1}{2} (-2t^0 - 2t^1) - (t^1 - 2t^2) \\
&= -t^0 - t^2 + t^0 + t^1 - t^1 + 2t^2 \\
&= t^2
\end{align*}$

$[w_0, w_1, w_2] Q^{-1} = [t^0, t^1, t^2]$
"""

# ╔═╡ 3dadd448-00c7-4ad7-8651-42a28c41eaf6
md"""
**(c)**

$\begin{align*}
x_0(t^0 + t^2) + x_1(-2t^0 - 2t^1) + x_2(t^1 - 2t^2) &= 0 \\
(x_0 t^0 + x_0 t^2) + (-2x_1 t^0 -2x_1 t^1) + (x_2 t^1 - 2x_2 t^2) &= 0 \\
(x_0 - 2x_1) t^0 + (-2x_1 + x_2) t^1 + (x_0 - 2x_2) t^2 &= 0
\end{align*}$

$\begin{align*}
x_0 - 2x_1 &= 0 \\
-2x_1 + x_2 &= 0 \\
x_0 - 2x_2 &= 0
\end{align*}\tag{*}$

$\begin{bmatrix}
1 & -2 & 0 \\
0 & -2 & 1 \\
1 & 0 & -2
\end{bmatrix}$
"""

# ╔═╡ 00c4a5ae-8c2b-43f5-a50a-3474998af970
md"""
**(d)**

$\begin{align*}
t^0 &= 2w_0 + \frac{1}{2} w_1 + w_2 \\
t^1 &= -2w_0 - w_1 - w_2 \\
t^2 &= -w_0 - \frac{1}{2} w_1 - w_2
\end{align*}$

They are related to $Q^{-1}$ because they show that $[w_0, w_1, w_2] Q^{-1}$ is a set of linear combinations of $w$.
"""

# ╔═╡ 04267575-996a-4c0a-b992-ea85eff2ca46
md"""
### PS28 #2

Let $V$ be the vector space in Problem 1.
We will use the basis $v_0, v_1, v_2$ in (a) and the basis $w_0, w_1, w_2$ in (b).
Define $T : V → V$ by

$T(f) = tf'(t) + f'(t) - 3f(t)$

(a) Find the matrix of $T$ using the basis $v_0, v_1, v_2$ of $V$, call it $A$, then find the characteristic polynomial of $A$.

(b) Find the matrix of $T$ using the basis $w_0, w_1, w_2$ of $V$, call it $B$, then find the characteristic polynomial of $B$.
"""

# ╔═╡ d935e71b-9357-49b0-b3bb-025dc8f8d7ce
md"""
**(a)**

$\begin{align*}
T(v_0) &= T(t^0) = t(0) + (0) - 3(t^0) = -3t^0 \\
T(v_1) &= T(t^1) = t(1) + (1) - 3(t^1) = t^0 - 2t^1 \\
T(v_2) &= T(t^2) = t(2t^1) + (2t^1) - 3(t^2) = 2t^1 - t^2
\end{align*}$

$A = \begin{bmatrix} -3 & 1 & 0 \\ 0 & -2 & 2 \\ 0 & 0 & -1 \end{bmatrix}$

$\begin{align*}p_A(\lambda) = \det(\lambda I - A) &= \begin{vmatrix} \lambda + 3 & 1 & 0 \\ 0 & \lambda + 2 & 2 \\ 0 & 0 & \lambda + 1 \end{vmatrix} \\
&= (\lambda + 3)(\lambda + 2)(\lambda + 1) \\
&= x^3 + 6x^2 + 11x + 6\end{align*}$
"""

# ╔═╡ 33023e36-c431-4a3a-8178-839a5ca2100d
md"""
**(b)**

$\begin{align*}
T(w_0) &= T(t^0 + t^2) = t(0 + 2t^1) + (0 + 2t^1) - 3(t^0 + t^2) = -3t^0 + 2t^1 - t^2 \\
T(w_1) &= T(-2t^0 - 2t^1) = t(0 - 2t^0) + (0 - 2t^0) - 3(-2t^0 - 2t^1) = 4t^0 + 4t^1 \\
T(w_2) &= T(t^1 - 2t^2) = t(t^0 - 4t^1) + (t^0 - 4t^1) - 3(t^1 - 2t^2) = t^0 - 6t^1 + 2t^2
\end{align*}$

$B = \begin{bmatrix} -3 & 4 & 1 \\ 2 & 4 & -6 \\ -1 & 0 & 2 \end{bmatrix}$

$\begin{align*}
p_B(\lambda) = \det(\lambda I - B) &= \begin{vmatrix} \lambda + 3 & 4 & 1 \\ 2 & \lambda - 4 & -6 \\ -1 & 0 & \lambda - 2 \end{vmatrix} \\
&= (-1) \begin{vmatrix} 4 & 1 \\ \lambda - 4 & -6 \end{vmatrix} + (\lambda - 2) \begin{vmatrix} \lambda + 3 & 4 \\ 2 & \lambda - 4 \end{vmatrix} \\
&= -[(4)(-6) - (1)(\lambda - 4)] + (\lambda - 2)[(\lambda + 3)(\lambda - 4) - (4)(2)] \\
&= -(-24 - \lambda + 4) + (\lambda - 2)(\lambda^2 - \lambda - 20) \\
&= (24 + \lambda - 4) + (\lambda^3 - \lambda^2 - 20 \lambda - 2\lambda^2 + 2\lambda + 40) \\
&= (20 + \lambda) + (\lambda^3 - 3 \lambda^2 - 18 \lambda + 40) \\
&= \lambda^3 - 3\lambda^2 - 17 \lambda + 60
\end{align*}$
"""

# ╔═╡ ac298ba7-8550-4b14-a4dd-3d0462370b4d
let
	@variables λ, t
	
	A = [-3 1 0; 0 -2 2; 0 0 -1]
	B = [-3 4 1; 2 4 -6; -1 0 2]
	det(λ*I - A), det(λ*I - B)
end

# ╔═╡ 06b44040-2288-4cc0-8ae9-b15ac15bfc42
md"""
### PS28 #3

Let $V$ be the span of the functions

$\begin{align*}
v_1 &= e^t \cos{t} \sin{t} \\
v_2 &= e^t \cos^2{t} \\
v_3 &= e^t \sin^2{t},
\end{align*}$

let $T : V → V$ be the operator $T(f) = f'$.

(a) Find the matrix of $T$ with respect to this basis of $V$, then the characteristic polynomial of that matrix.

(b) Now let

$\begin{align*}
w_1 &= e^t \sin(2t), \\
w_2 &= e^t, \\
w_3 &= e^t \cos(2t).
\end{align*}$

These functions form another basis of $V$.
Find the matrix of $T$ with respect to this basis, then its characteristic polynomial.
"""

# ╔═╡ a7120b38-5dd8-4fac-8a89-9bd24ecdce96
md"""
**(a)**

$\begin{align*}
T(v_1) = T(e^t \cos{t} \sin{t}) &= e^t \cos{t} \sin{t} + e^t \cos{(2t)} \\
&= e^t \cos{t} \sin{t} + e^t (\cos^2{t} - \sin^2{t}) \\
&= e^t \cos{t} \sin{t} + e^t \cos^2{t} - e^t \sin^2{t} \\
T(v_2) = T(e^t \cos^2{t}) &= -2e^t \cos{t} \sin{t} + e^t \cos^2{t} \\
T(v_3) = T(e^t \sin^2{t}) &= 2e^t \cos{t} \sin{t} + e^t \sin^2{t}
\end{align*}$

$A = \begin{bmatrix} 1 & -2 & 2 \\ 1 & 1 & 0 \\ -1 & 0 & 1 \end{bmatrix}$

$\begin{align*}
p_A(\lambda) = \det(\lambda I - A) &= \begin{vmatrix} \lambda - 1 & -2 & 2 \\ 1 & \lambda - 1 & 0 \\ -1 & 0 & \lambda - 1 \end{vmatrix} \\
&= (-1) \begin{vmatrix} -2 & 2 \\ \lambda - 1 & 0 \end{vmatrix} + (\lambda - 1) \begin{vmatrix} \lambda - 1 & -2 \\ 1 & \lambda - 1 \end{vmatrix} \\
&= -[0 - 2(\lambda - 1)] + (\lambda - 1) [(\lambda - 1)(\lambda - 1) + 2] \\
&= 2 \lambda - 2 + (\lambda - 1)(\lambda^2 - 2 \lambda + 3) \\
&= 2 \lambda - 2 + \lambda^3 - 2\lambda^2 + 3\lambda - \lambda^2 + 2\lambda - 3 \\
&= \lambda^3 - 3\lambda^2 + 7\lambda - 5
\end{align*}$
"""

# ╔═╡ f4268677-e115-4232-bd7f-6e3379c355c3
md"""
**(b)**

$\begin{align*}
T(w_1) &= T(2v_1) = 2T(v_1) = 2v_1 + 2v_2 - 2v_3 \\
T(w_2) &= T(v_2 + v_3) = T(v_2) + T(v_3) = v_2 + v_3 \\
T(w_3) &= T(v_2 - v_3) = T(v_2) - T(v_3) = -4v_1 + v_2 - v_3
\end{align*}$

$A = \begin{bmatrix} 2 & 0 & -4 \\ 2 & 1 & 1 \\ -2 & 1 & -1 \end{bmatrix}$

$\begin{align*}
p_A(\lambda) = \det(\lambda I - A) &= \begin{vmatrix} \lambda - 2 & 0 & -4 \\ 2 & \lambda - 1 & 1 \\ -2 & 1 & \lambda + 1 \end{vmatrix} \\
&= (\lambda - 2) \begin{vmatrix} \lambda - 1 & 1 \\ 1 & \lambda + 1 \end{vmatrix} + (-4) \begin{vmatrix} 2 & \lambda - 1 \\ -2 & 1 \end{vmatrix} \\
&= (\lambda - 2)[(\lambda - 1)(\lambda + 1) - 1] + (-4)[2 + 2(\lambda - 1)] \\
&= (\lambda - 2)(\lambda^2 - 2) - 4(2\lambda) \\
&= \lambda^3 - 2\lambda - 2\lambda^2 + 4 - 8\lambda \\
&= \lambda^3 - 2\lambda^2 - 10\lambda + 20
\end{align*}$
"""

# ╔═╡ 433d717b-2588-4f06-a8bf-dcbe91381c95
md"""
### PS28 #3 (Again, for practice)

$T : V → V \qquad T(f) = f'$
"""

# ╔═╡ da7dac05-c8e5-4792-98c1-ec133effef57
md"""
**(a)**

$\begin{align*}
T(v_1) &= T(e^t \cos{t} \sin{t}) = e^t \cos{t} \sin{t} + e^t (\cos^2{t} - \sin^2{t}) \\
&= e^t \cos{t} \sin{t} + e^t \cos^2{t} - e^t \sin^2{t} \\
&= v_1 + v_2 - v_3 \\
T(v_2) &= T(e^t \cos^2{t}) \\
&= e^t \cos^2{t} - 2e^t \cos{t} \sin{t} \\
&= -2v_1 + v_2 \\
T(v_3) &= T(e^t \sin^2{t}) \\
&= e^t \sin^2{t} + 2 e^t \cos{t} \sin{t} \\
&= 2v_1 + v_3
\end{align*}$

$A = \begin{bmatrix} 1 & -2 & 2 \\ 1 & 1 & 0 \\ -1 & 0 & 1 \end{bmatrix}$

$p_A(\lambda) = \begin{vmatrix} \lambda - 1 & -2 & 2 \\ 1 & \lambda - 1 & 0 \\ -1 & 0 & \lambda - 1 \end{vmatrix}$
"""

# ╔═╡ d577293d-a943-45d7-9d5b-e95a292ecb59
md"""
**(b)**

$\begin{align*}
T(w_1) &= T(2v_1) = 2T(v_1) \\
&= 2(v_1 + v_2 - v_3) \\
&= 2v_1 + 2v_2 - 2v_3 \\
T(w_2) &= T(v_2 + v_3) = T(v_2) + T(v_3) \\
&= (-2v_1 + v_2) + (2v_1 + v_3) \\
&= v_2 + v_3 \\
T(w_3) &= T(v_2 - v_3) = T(v_2) - T(v_3) \\
&= (-2v_1 + v_2) - (2v_1 + v_3) \\
&= -4v_1 + v_2 - v_3
\end{align*}$

$B = \begin{bmatrix} 2 & 0 & -4 \\ 2 & 1 & 1 \\ -2 & 1 & -1 \end{bmatrix}$

$p_B(\lambda) = \begin{vmatrix} \lambda - 2 & 0 & -4 \\ 2 & \lambda - 1 & 1 \\ -2 & 1 & \lambda + 1 \end{vmatrix}$
"""

# ╔═╡ c668ebdd-3e0c-4c98-a336-f8933c564565
let
	@variables t, x
	f(t) = exp(t) * cos(t)^2 + x^(2t)
	Symbolics.derivative(f(t), x)
end

# ╔═╡ 84d6fe1a-649c-4668-a5f6-cba0cf73d0d5
md"## Problem Set 29"

# ╔═╡ e7d69866-97fa-4914-a2db-31073f2bae12
md"""
### PS29 #1

For the system

$\begin{align*}
x_1 + 3x_2 + 2x_3 &= 1 \\
x_2 + 2x_3 &= 1 \\
x_1 + 2x_3 &= 4
\end{align*}$

find:

(a) The matrix of the system

(b) The augmented matrix of the system

(c) The solution.
"""

# ╔═╡ 48cd9ff9-ba8e-4f13-9c1f-6620db4d793a
md"""
**Solution (a).** The matrix of the system is

$\begin{bmatrix} 1 & 3 & 2 \\ 0 & 1 & 2 \\ 1 & 0 & 2 \end{bmatrix}$
"""

# ╔═╡ 56d69f1f-d55b-4d75-a291-7fb78eb92c2f
md"""
**Solution (b).** The augmented matrix of the system is

$\begin{bmatrix} 1 & 3 & 2 & 1 \\ 0 & 1 & 2 & 1 \\ 1 & 0 & 2 & 4 \end{bmatrix}$
"""

# ╔═╡ 5cb455fa-acb7-48a3-8c70-95a04c886c84
md"""
**Solution (c).**

$\begin{align*}
\begin{bmatrix} 1 & 3 & 2 & 1 \\ 0 & 1 & 2 & 1 \\ 1 & 0 & 2 & 4 \end{bmatrix}
&\rightsquigarrow \begin{bmatrix} 1 & 3 & 2 & 1 \\ 0 & 1 & 2 & 1 \\ 0 & 1 & 0 & -1 \end{bmatrix} \\
&\rightsquigarrow \begin{bmatrix} 1 & 3 & 2 & 1 \\ 0 & 1 & 2 & 1 \\ 0 & 0 & 1 & 1 \end{bmatrix} \\
&\rightsquigarrow \begin{bmatrix} 1 & 3 & 0 & -1 \\ 0 & 1 & 0 & -1 \\ 0 & 0 & 1 & 1 \end{bmatrix} \\
&\rightsquigarrow \begin{bmatrix} 1 & 0 & 0 & 2 \\ 0 & 1 & 0 & -1 \\ 0 & 0 & 1 & 1 \end{bmatrix}
\end{align*}$

$\begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} = \begin{bmatrix} 2 \\ -1 \\ 1 \end{bmatrix}$
"""

# ╔═╡ d2836e79-c66a-4eef-90ed-bf669774b7f6
md"""
### PS29 #2

Does the system

$\begin{align*}
-2x_1 + x_2 - 4x_3 &= -9 \\
x_1 - x_2 + 2x_3 &= 5 \\
x_1 + x_2 + 2x_3 &= 3
\end{align*}$

have infinitely many solutions?
"""

# ╔═╡ 5f41870b-5395-4642-bfc9-2a9aa7c47817
md"""
**Solution.**

$\begin{align*}
\begin{bmatrix} -2 & 1 & -4 & -9 \\ 1 & -1 & 2 & 5 \\ 1 & 1 & 2 & 3 \end{bmatrix}
&\rightsquigarrow \begin{bmatrix} -2 & 1 & -4 & -9 \\ 1 & -1 & 2 & 5 \\ 0 & 1 & 0 & -1 \end{bmatrix} \\
&\rightsquigarrow \begin{bmatrix} 1 & -1 & 2 & 5 \\ -2 & 1 & -4 & -9 \\ 0 & 1 & 0 & -1 \end{bmatrix} \\
&\rightsquigarrow \begin{bmatrix} 1 & -1 & 2 & 5 \\ 0 & -1 & 0 & 1 \\ 0 & 1 & 0 & -1 \end{bmatrix} \\
&\rightsquigarrow \begin{bmatrix} 1 & 0 & 2 & 4 \\ 0 & 1 & 0 & -1 \\ 0 & 0 & 0 & 0 \end{bmatrix} \\
\end{align*}$

Yes. There is a free variable ``x_3``.
"""

# ╔═╡ 2e160eac-4682-446c-a67c-d568718bac7b
md"""
### PS29 #3

Is the system

$\begin{align*}
-2x_1 + x_2 - 4x_3 &= -9 \\
x_1 - x_2 + 2x_3 &= 5 \\
x_1 + x_2 + 2x_3 &= 4
\end{align*}$

consistent?
"""

# ╔═╡ f6ed271c-8ace-432e-bb5a-5e04f4944629
md"""
**Solution.**

$\begin{align*}
\begin{bmatrix} -2 & 1 & -4 & -9 \\ 1 & -1 & 2 & 5 \\ 1 & 1 & 2 & 4 \end{bmatrix}
&\rightsquigarrow \begin{bmatrix} -2 & 1 & -4 & -9 \\ 1 & -1 & 2 & 5 \\ 0 & 2 & 0 & -1 \end{bmatrix} \\
&\rightsquigarrow \begin{bmatrix} 1 & -1 & 2 & 5 \\ -2 & 1 & -4 & -9 \\ 0 & 2 & 0 & -1 \end{bmatrix} \\
&\rightsquigarrow \begin{bmatrix} 1 & -1 & 2 & 5 \\ 0 & -1 & 0 & 1 \\ 0 & 2 & 0 & -1 \end{bmatrix} \\
&\rightsquigarrow \begin{bmatrix} 1 & 0 & 2 & 4 \\ 0 & 1 & 0 & -1 \\ 0 & 2 & 0 & -1 \end{bmatrix} \\
&\rightsquigarrow \begin{bmatrix} 1 & 0 & 2 & 4 \\ 0 & 1 & 0 & -1 \\ 0 & 0 & 0 & 1 \end{bmatrix} \\
\end{align*}$

No, it is not consistent---there is no solution.
"""

# ╔═╡ 55e273db-4389-4a7e-84c8-bdb788c20b23
md"""
### PS29 #4

Suppose that for a certain system ``Ax = b`` of three equations and three unknowns, the reduced echelon of the augmented matrix ``\begin{bmatrix} A ∣ b \end{bmatrix}`` happens to be

$\begin{bmatrix} 1 & 0 & -2 & 5 \\ 0 & 1 & 0 & -1 \\ 0 & 0 & 0 & -1 \end{bmatrix}$

Does the system ``Ax = b`` have a solution?
"""

# ╔═╡ c059e395-23cd-495c-98e7-1a820af5bf74
md"""
**Solution.** No.
"""

# ╔═╡ d37ae0d5-79ae-478d-b84a-be068ec27f20
md"""
### PS29 #5

Find the rank, nullity, and null space of the matrix

$A = \begin{bmatrix} -2 & 1 & -4 & -7 & 3 \\ 1 & -1 & 2 & 4 & -2 \\ 1 & 1 & 2 & 2 & 0 \\ 2 & 3 & 1 & 0 & -2 \end{bmatrix}$

given that its reduced echelon form is

$\begin{bmatrix} 1 & 0 & 0 & 1 & -3 \\ 0 & 1 & 0 & -1 & 1 \\ 0 & 0 & 1 & 1 & 1 \\ 0 & 0 & 0 & 0 & 0 \end{bmatrix}.$
"""

# ╔═╡ a3a9f6d6-e29d-4810-9d29-760af87dd4da
md"""
**Solution.**

$\text{rank}(A) = 3$

$\text{nullity}(A) = 2$

$\text{nullspace}(A) = \left\{\begin{bmatrix} -1 \\ 1 \\ -1 \\ 1 \\ 0 \end{bmatrix}, \begin{bmatrix} 3 \\ -1 \\ -1 \\ 0 \\ 1 \end{bmatrix}\right\}$
"""

# ╔═╡ b21cc3bd-64e3-4690-a54a-563b66f93403
let
	A = [-2 1 -4 -7 3; 1 -1 2 4 -2; 1 1 2 2 0; 2 3 1 0 -2]
	nullspace(A)
end

# ╔═╡ dfc4e6fe-694d-4669-b962-8e00d7fc53c3
md"""
### PS29 #6

Find bases for the column and row spaces of the matrix ``A`` in Problem 5.
"""

# ╔═╡ b102a4b9-d81f-4bcd-893f-e39deb0a83d9
md"""
**Solution.**

$\text{Basis of col}(A) = \left\{\begin{bmatrix} -2 \\ 1 \\ 1 \\ 2 \end{bmatrix}, \begin{bmatrix} 1 \\ -1 \\ 1 \\ 3 \end{bmatrix}, \begin{bmatrix} -4 \\ 2 \\ 2 \\ 1 \end{bmatrix}\right\}$

$\text{Basis of row}(A) = \left\{\begin{bmatrix} 1 & 0 & 0 & 1 & -3 \end{bmatrix}, \begin{bmatrix} 0 & 1 & 0 & -1 & 1 \end{bmatrix}, \begin{bmatrix} 0 & 0 & 1 & 1 & 1 \end{bmatrix}\right\}$
"""

# ╔═╡ 4bad5e75-862b-40cc-bdb8-aeb1014ca125
md"""
### PS29 #7

Let ``A`` be the matrix in Problem 5.
Suppose ``M`` is an invertible ``4 × 4`` matrix.
What is the relation between the reduced echelon forms of ``A`` and ``MA``?
(Equal or not equal, and why.)
"""

# ╔═╡ 9e05f711-ecb6-4923-a6eb-8ac42d128c2c
md"""
**Solution.** Equal.
The columns of ``MA`` are linear combination of the columns of ``A``, so the reduced echelon form of ``A`` and ``MA`` will be the same if ``M`` is invertible, but not if ``M`` is non-invertible because the reduced echelon form of a non-invertible matrix contains a row of zeros.
"""

# ╔═╡ 56fc50f7-88d0-46d6-abb3-de8a2cc58e12
md"""
### PS29 #8

Let

$P = \begin{bmatrix} 1 & 1 & 5 \\ 2 & -1 & 1 \end{bmatrix}, \quad Q = \begin{bmatrix} -2 & 1 \\ 1 & 4 \\ 1 & -1 \end{bmatrix}$

Find ``PQ``.
"""

# ╔═╡ b8ed2d91-ec62-4077-b7e2-fa37b16fdb9c
let
	P = [1 1 5; 2 -1 1]
	Q = [-2 1; 1 4; 1 -1]
	P*Q
end

# ╔═╡ 4de3766f-fc03-4a15-8e8c-4cea7eff3b35
md"""
### PS29 #9

The matrix ``PQ`` in Problem 8 is invertible.
Let ``R`` be its inverse.

(a) Is ``QR`` a right inverse of ``P``?

(b) If so, use it to find a solution of

$\begin{bmatrix} 1 & 1 & 5 \\ 2 & -1 & 1 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} = \begin{bmatrix} -1 \\ 1 \end{bmatrix}$

(c) Is this solution unique?
"""

# ╔═╡ 97e9c23b-b775-413b-960b-ac1c9e39f002
md"""
**(a)** Yes.

$P(QR) = (PQ)R = I$
"""

# ╔═╡ ce5e862f-2799-40ef-aa3d-f4fde538119a
md"""
**(b)**

$\begin{align*}
P \begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} &= \begin{bmatrix} -1 \\ 1 \end{bmatrix} \\
P(QR) \begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} &= (QR) \begin{bmatrix} -1 \\ 1 \end{bmatrix} \\
I \begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} &= QR \begin{bmatrix} -1 \\ 1 \end{bmatrix} \\
\begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} &= \begin{bmatrix} -\frac{5}{6} & -\frac{1}{3} \\ -\frac{13}{12} & -\frac{4}{3} \\ \frac{7}{12} & \frac{1}{3} \end{bmatrix} \begin{bmatrix} -1 \\ 1 \end{bmatrix} \\
&= \begin{bmatrix} \frac{5}{6} - \frac{1}{3} \\ \frac{13}{12} - \frac{4}{3} \\ -\frac{7}{12} + \frac{1}{3} \end{bmatrix} \\
&= \begin{bmatrix} \frac{1}{2} \\ -\frac{1}{4} \\ -\frac{1}{4} \end{bmatrix} \\
\end{align*}$
"""

# ╔═╡ 64f1c0b6-dd27-4549-a894-7f69c873a0c1
md"""
**(c)** No.

$\begin{bmatrix} 1 & 1 & 5 & -1 \\ 2 & -1 & 1 & 1 \end{bmatrix} \rightsquigarrow \begin{bmatrix} 1 & 0 & 2 & 0 \\ 0 & 1 & 3 & -1 \end{bmatrix}$

$\begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} = \begin{bmatrix} -2x_3 \\ -1 - 3x_3 \end{bmatrix} = \begin{bmatrix} 0 \\ -1 \end{bmatrix} + \begin{bmatrix} -2x_3 \\ -3x_3 \end{bmatrix} = \begin{bmatrix} 0 \\ -1 \end{bmatrix} + x_3 \begin{bmatrix} -2 \\ -3 \end{bmatrix}$
"""

# ╔═╡ d01b32ed-4d8d-4a69-b6d0-69cea9beb81b
md"""
### PS29 #10

Let

$K = \begin{bmatrix} 1 & 1 \\ 1 & -1 \\ 1 & 0 \end{bmatrix}$

(a) The matrix ``K^t K`` is invertible.

(b) Find a left inverse for ``K``.

(c) Suppose ``x = [x_1, x_2]`` solves ``Kx = 0``.
Use that ``K^t K`` (or (b)) is invertible to verify that ``x`` must be the ``0`` vector.
"""

# ╔═╡ 64239289-a2f6-4d9e-b67d-df540a439ef7
md"""
**Solution (a).**

$\begin{bmatrix} 1 & 1 & 1 \\ 1 & -1 & 0 \end{bmatrix} \begin{bmatrix} 1 & 1 \\ 1 & -1 \\ 1 & 0 \end{bmatrix} = \begin{bmatrix} 3 & 0 \\ 0 & 2 \end{bmatrix}$
"""

# ╔═╡ 3c1a2ea2-d8a9-4b7e-8360-842d1c053529
md"""
**Solution (b).**

$\begin{align*}
(K^t K)^{-1} (K^t K) &= I \\
((K^t K)^{-1} K^t) K &= I
\end{align*}$

$(K^t K)^{-1} K^t = \begin{bmatrix} 3 & 0 \\ 0 & 2 \end{bmatrix} \begin{bmatrix} 1 & 1 & 1 \\ 1 & -1 & 0 \end{bmatrix} = \begin{bmatrix} 3 & 3 & 3 \\ 2 & -2 & 0 \end{bmatrix}$
"""

# ╔═╡ b29a70f1-e6a7-49cc-8800-9ae4d0b8e0ab
md"""
**Solution (c).**

$\begin{align*}
((K^t K)^{-1} K^t) Kx &= ((K^t K)^{-1} K^t) 0 \\
x &= 0
\end{align*}$
"""

# ╔═╡ c8feede5-ecb0-4f54-a95b-013f773c05fe
md"""
### PS29 #11

The list of vectors

$v_1 = \begin{bmatrix} 1 \\ 3 \\ 2 \\ 0 \\ -2 \end{bmatrix}, \; v_2 = \begin{bmatrix} 0 \\ 7 \\ 3 \\ -5 \\ -8 \end{bmatrix}, \; v_3 = \begin{bmatrix} 2 \\ -1 \\ 1 \\ 5 \\ 4 \end{bmatrix}, \; v_4 = \begin{bmatrix} -2 \\ 0 \\ -1 \\ -4 \\ -4 \end{bmatrix}, v_5 = \begin{bmatrix} 0 \\ -1 \\ 0 \\ 1 \\ 0 \end{bmatrix}$

contains redundant vectors.
Remove enough of them so that those that remain form a basis of ``\text{span}\{v_1, v_2, v_3, v_4, v_5\}``.
"""

# ╔═╡ b657c4e2-4d1b-466e-bc51-5146eb03f3e6
md"""
**Solution.**
Consider the matrix $A = \begin{bmatrix} v_1 & v_2 & v_3 & v_4 & v_5 \end{bmatrix}$.
The basis of ``\text{col}(A) = \text{span}\{v_1, v_2, v_4\}``.
"""

# ╔═╡ ea1c0899-355e-49c9-a05d-28619faa978e
md"""
### PS29 #12

Suppose that the vectors ``v_1, v_2, v_3`` are linearly independent.
Let

$w_1 = v_1 + v_2 + 3v_3, \quad w_2 = v_1 + v_3, \quad w_3 = -2v_1 + v_2 + 2v_3$

Verify that ``w_1``, ``w_2`` and ``w_3`` are independent.
To do this, establish that if

$x_1 w_1 + x_2 w_2 + x_3 w_3 = 0,$

then ``x_1 = x_2 = x_3 = 0``.
"""

# ╔═╡ b2d46287-301b-4a9a-b4e7-e1ff9d0d5aa7
md"""
**Solution.**

$\begin{align*}
x_1 (v_1 + v_2 + 3v_3) + x_2 (v_1 + v_3) + x_3 (-2v_1 + v_2 + 2v_3) &= 0 \\
(x_1 + x_2 - 2x_3) v_1 + (x_1 + x_3) v_2 + (3x_1 + x_2 + 2x_3) v_3 &= 0
\end{align*}$

$\begin{align*}
x_1 + x_2 - 2x_3 &= 0 \\
x_1 + x_3 &= 0 \\
3x_1 + x_2 + 2x_3 &= 0
\end{align*}$

$\begin{bmatrix} 1 & 1 & -2 \\ 1 & 0 & 1 \\ 3 & 1 & 2 \end{bmatrix} \rightsquigarrow \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix}$

$\implies \begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \\ 0 \end{bmatrix}$
"""

# ╔═╡ 18170309-b97b-4ed3-ab69-3c2c640a0d9c
md"""
### PS29 #13

Let the ``v_i`` and ``w_j`` be as in Problem 12.
Let ``w_4 = v_1 - v_2 + v_3``.
Are the vectors ``w_1, w_2, w_3, w_4`` linearly independent?
"""

# ╔═╡ cec5575a-4c2d-4fd9-8a44-7b5efb564a22
md"""
**Solution.**

$\begin{bmatrix} 1 & 1 & -2 & 1 \\ 1 & 0 & 1 & -1 \\ 3 & 1 & 2 & 1 \end{bmatrix} \rightsquigarrow \begin{bmatrix} 1 & 0 & 0 & -2 \\ 0 & 1 & 0 & 5 \\ 0 & 0 & 1 & 1 \end{bmatrix}$

$\implies \begin{bmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{bmatrix} = x_4 \begin{bmatrix} 2 \\ -5 \\ -1 \\ 1 \end{bmatrix}$

$x_1 w_1 + x_2 w_2 + x_3 w_3 + x_4 w_4 = 0$

No, they are not independent because solutions of ``x`` are not strictly ``0``.
"""

# ╔═╡ f399418a-3ecd-4184-bde0-cdcc7734a442
md"""
### PS29 #14

The vectors ``v_1 = (1,0,1)``, ``v_2 = (0,1,1)``, ``v_3 = (1,-1,0)`` are linearly dependent.
Verify that any two of them are linearly independent.
"""

# ╔═╡ 9dcd9b99-9935-4ae1-8069-6080d27549d0
md"""
**Solution.**

$\begin{bmatrix} v_1 & v_2 \end{bmatrix} = \begin{bmatrix} 1 & 0 \\ 0 & 1 \\ 1 & 1 \end{bmatrix} \rightsquigarrow \begin{bmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{bmatrix}$

$\begin{align*}
v_3 &= x_1 v_1 + x_2 v_2 \\
\begin{bmatrix} 1 \\ -1 \\ 0 \end{bmatrix} &= x_1 \begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix} + x_2 \begin{bmatrix} 0 \\ 1 \\ 1 \end{bmatrix}
\end{align*}$
"""

# ╔═╡ a001e88c-c34a-4c8e-883c-0bbec2dbbc6b
let
	v1 = [1; 0; 1]
	v2 = [0; 1; 1]
	v3 = [1; -1; 0]
	rref([v1 v2]), rref([v1 v3]), rref([v2 v3])
end

# ╔═╡ 41978ff8-56d6-40ac-9d05-1493355a43e6
md"""
### PS29 #15

Let ``A`` be the matrix

$\begin{bmatrix} 1 & 1 & 3 & -1 \\ 2 & -1 & 0 & 4 \\ 3 & 0 & 3 & 4 \end{bmatrix},$

let

$E_1 = \begin{bmatrix} 1 & 0 & 0 \\ -2 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix}, \; E_2 = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ -3 & 0 & 1 \end{bmatrix}, \; E_3 = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -1 & 1 \end{bmatrix}$

$E_4 = \begin{bmatrix} 1 & 0 & 0 \\ 0 & -\frac{1}{3} & 0 \\ 0 & 0 & 1 \end{bmatrix}, \; E_5 = \begin{bmatrix} 1 & -1 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix}$

These elementary matrices yield the reduced echelon form of ``A``.

$E_5 E_4 E_3 E_2 E_1 A = \begin{bmatrix} 1 & 0 & 1 & 1 \\ 0 & 1 & 2 & -2 \\ 0 & 0 & 0 & 0 \end{bmatrix}$

Let ``Q = E_5 E_4 E_3 E_2 E_1``, find ``\det(Q)``.
"""

# ╔═╡ 35c96ee5-e003-4f64-979e-b4d466032a11
md"""
**Solution.**

$\det(Q) = -\frac{1}{3}$
"""

# ╔═╡ 61a03d8e-bf2c-48fa-b841-b81474fdafcb
md"""
### PS29 #16

Let ``B`` consist of the first three columns of the matrix ``A`` in Problem 15:

$B = \begin{bmatrix} 1 & 1 & 3 \\ 2 & -1 & 0 \\ 3 & 0 & 3 \end{bmatrix}$

The matrices ``E_j`` in Problem 15, whose product is

$Q = \begin{bmatrix} \frac{1}{3} & \frac{1}{3} & 0 \\ \frac{2}{3} & -\frac{1}{3} & 0 \\ -1 & -1 & 1 \end{bmatrix}$

were found taking into account only the first three columns of ``A``.
So they also work for ``B``:

$QB = \begin{bmatrix} 1 & 0 & 1 \\ 0 & 1 & 2 \\ 0 & 0 & 0 \end{bmatrix}$

Use this to find a solution of

$\begin{bmatrix} 1 & 1 & 3 \\ 2 & -1 & 0 \\ 3 & 0 & 3 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} = \begin{bmatrix} 0 \\ 3 \\ 3 \end{bmatrix}$

without starting over with row reduction; just use ``Q``.
"""

# ╔═╡ 1cde7920-ee1c-4877-b58d-a9dd41d590b3
md"""
**Solution.**

$\left[\begin{array}{c|c} QB & Q\begin{bmatrix} 0\\3\\3\end{bmatrix}\end{array}\right] = \begin{bmatrix} 1 & 0 & 1 & 1 \\ 0 & 1 & 2 & -1 \\ 0 & 0 & 0 & 0 \end{bmatrix}$

$\begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} = \begin{bmatrix} 1 - x_3 \\ -1 - 2x_3 \\ x_3 \end{bmatrix} = \begin{bmatrix} 1 \\ -1 \\ 0 \end{bmatrix} + x_3 \begin{bmatrix} -1 \\ -2 \\ 1 \end{bmatrix}$

$x = \begin{bmatrix} 1 \\ -1 \\ 0 \end{bmatrix}$
"""

# ╔═╡ 6bc1c224-ab48-4922-8c66-e20ababbda62
md"""
### PS29 #17

Find the cofactor matrix of

$R = \begin{bmatrix} 1 & 0 & 0 \\ -1 & 2 & 0 \\ -1 & 1 & 1 \end{bmatrix}.$
"""

# ╔═╡ 5ca73e2a-87c4-4d8b-831f-7dea471accf9
md"""
**Solution.**

$C = \begin{bmatrix} 2 & 1 & 1 \\ 0 & 1 & -1 \\ 0 & 0 & 2 \end{bmatrix}$

$C^t R = \begin{bmatrix} 2 & 0 & 0 \\ 1 & 1 & 0 \\ 1 & -1 & 2 \end{bmatrix} \begin{bmatrix} 1 & 0 & 0 \\ -1 & 2 & 0 \\ -1 & 1 & 1 \end{bmatrix} = \begin{bmatrix} 2 & 0 & 0 \\ 0 & 2 & 0 \\ 0 & 0 & 2 \end{bmatrix}$
"""

# ╔═╡ f5828a13-3c4d-49ac-8601-2f758b6bba84
md"""
### PS29 #18

Let ``A`` be a ``4 × 4`` matrix with columns ``A_1, A_2, A_3, A_4`` and determinant 3, let ``B`` be the matrix with columns

$B_1 = A_2, \; B_2 = 2A_1 - A_2, \; B_3 = A_1 + A_3 + A_4, \; B_4 = 3A_1 - 2A_4$

Find ``\det(B)``
"""

# ╔═╡ 6444d766-8635-4b50-bc7c-7223f0353656
md"""
**Solution.**

$B_1 \implies (0, 1, 0, 0)$

$C = \begin{bmatrix} 0 & 2 & 1 & 3 \\ 1 & -1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 1 & -2 \end{bmatrix}$

$\begin{align*}
B &= AC \\
\det(B) &= \det(AC) \\
&= \det(A) \det(C) \\
&= (3)(4) \\
&= 12
\end{align*}$
"""

# ╔═╡ 8a199215-6456-4cac-827b-2ea3c793b1f7
md"""
### PS29 #19

(a) Find a basis for each of the eigenspaces of ``Q = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ -2 & 0 & 3 \end{bmatrix}``.
Taken together, do they form a basis of ``ℝ^3``?

(b) Find a basis for each of the generalized eigenspaces of ``Q``.
"""

# ╔═╡ 3f95062c-986a-486c-bc19-31b6e3c19ece
md"""
**Solution (a).**

$\det(\lambda I - Q) = \begin{vmatrix} \lambda - 1 & 0 & 0 \\ 0 & \lambda - 1 & 0 \\ -2 & 0 & \lambda - 3 \end{vmatrix} = (\lambda - 1)^2 (\lambda - 3)$

$Q - I = \begin{bmatrix} 0 & 0 & 0 \\ 0 & 0 & 0 \\ -2 & 0 & 2 \end{bmatrix} \rightsquigarrow \begin{bmatrix} 1 & 0 & -1 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{bmatrix} \implies 𝐯 = \begin{bmatrix} x_3 \\ x_2 \\ x_3 \end{bmatrix} = x_2 \begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix} + x_3 \begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix}$

The basis for $\lambda = 1$ is

$\left\{\begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix}, \begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix}\right\}$

$Q - 3I = \begin{bmatrix} -2 & 0 & 0 \\ 0 & -2 & 0 \\ -2 & 0 & 0 \end{bmatrix} \rightsquigarrow \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{bmatrix} \implies 𝐯 = \begin{bmatrix} 0 \\ 0 \\ x_3 \end{bmatrix} = x_3 \begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix}$

The basis for $\lambda = 3$ is

$\left\{\begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix}\right\}$
"""

# ╔═╡ 15975b30-f690-4c03-ad2e-9c282417999f
md"""
**Solution (b).**

$(Q - I)^2 = \begin{bmatrix} 0 & 0 & 0 \\ 0 & 0 & 0 \\ -4 & 0 & 4 \end{bmatrix} \rightsquigarrow \begin{bmatrix} 1 & 0 & -1 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{bmatrix} \implies 𝐯 = x_2 \begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix} + x_3 \begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix}$

The basis of the generalized eigenspaces for $\lambda = 1$ and $\lambda = 3$ are

$\left\{\begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix}, \begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix}\right\}, \text{ and } \left\{\begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix}\right\}$

respectively.
"""

# ╔═╡ 8c68ed02-6e2a-45a6-87dc-2e5416cfbd05
md"""
### PS29 #20

(a) Find a basis for each of the eigenspaces of ``P = \begin{bmatrix} 1 & 1 & 0 \\ 0 & 1 & 0 \\ -2 & 1 & 3 \end{bmatrix}``.
Taken together, do they form a basis of ``ℝ^3``?

(b) Find a basis for each of generalized eigenspaces of ``P``.
"""

# ╔═╡ da18d4c7-5f2b-453d-84fe-3ad32605089e
md"""
**Solution (a).**

$\det(\lambda I - P) = \begin{vmatrix} \lambda - 1 & 1 & 0 \\ 0 & \lambda - 1 & 0 \\ -2 & 1 & \lambda - 3 \end{vmatrix} = (\lambda - 1)^2 (\lambda - 3)$

$P - I = \begin{bmatrix} 0 & 1 & 0 \\ 0 & 0 & 0 \\ -2 & 1 & 2 \end{bmatrix} \rightsquigarrow \begin{bmatrix} 1 & 0 & -1 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{bmatrix} \implies 𝐯 = \begin{bmatrix} x_3 \\ 0 \\ x_3 \end{bmatrix} = x_3 \begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix}$

$P - 3I = \begin{bmatrix} -2 & 1 & 0 \\ 0 & -2 & 0 \\ -2 & 1 & 0 \end{bmatrix} \rightsquigarrow \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{bmatrix} \implies 𝐯 = \begin{bmatrix} 0 \\ 0 \\ x_3 \end{bmatrix} = x_3 \begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix}$

No, they do not form a basis of ``ℝ^3`` because there are only two independent vectors.
"""

# ╔═╡ 9bf217de-c29e-4827-9c02-ffb916f2b2c9
rref([1 1 0; 0 1 0; -2 1 3] - I)^2

# ╔═╡ 663a0508-fb15-444d-bd42-9e4434e92016
md"""
**Solution (b).**

$(P - I.)^2 = \begin{bmatrix} 0 & 0 & 0 \\ 0 & 0 & 0 \\ -4 & 0 & 4 \end{bmatrix} \rightsquigarrow \begin{bmatrix} 1 & 0 & -1 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{bmatrix} \implies 𝐯 = x_2 \begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix} + x_3 \begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix}$

The basis of the generalized eigenspaces for $\lambda = 1$ and $\lambda = 3$ are

$\left\{\begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix}, \begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix}\right\}, \text{ and } \left\{\begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix}\right\}$

respectively.
"""

# ╔═╡ b9a063bd-3a71-4666-babe-9eca110188db
md"""
### PS29 #21

Let ``A`` be an ``n × n`` matrix, let ``v`` be an eigenvector of ``A`` with eigenvalue ``\lambda_0``.

(a) Let ``Q = 2A^2 - A + I``.
Verify that ``v`` is an eigenvector of ``Q`` with eigenvalue ``\mu = 2{\lambda_0}^2 - 2\lambda_0 + 1``. (Why is ``Qv = \mu v``?)

(b) Let ``P = A^3 - A``. Verify that ``PQ = QP``.

(c) Is ``({\lambda_0}^3 - \lambda_0)(2{\lambda_0}^2 - 2\lambda_0 + 1)`` an eigenvalue of ``PQ``? Yes or no, but explain.
"""

# ╔═╡ 9b6d402e-5c19-4c71-949a-8ed0a107f82c
md"""
**Solution (a).**

$\begin{align*}
Qv &= (2A^2 - A + I)v \\
&= 2A^2 v - Av + v \\
&= 2{\lambda_0}^2 v - 2{\lambda_0} v + 1 v \\
&= (2{\lambda_0}^2 - 2{\lambda_0} + 1) v \\
&= \mu v
\end{align*}$
"""

# ╔═╡ a4f6faf5-f770-4af6-884a-f53f49dd73d5
md"""
**Solution (b).**

$\begin{align*}
PQ &= (A^3 - A)(2A^2 - A + I) \\
&= 2A^5 - A^4 + A^3 - 2A^3 + A^2 - A \\
&= 2A^5 - A^4 - A^3 + A^2 - A \\
\end{align*}$

$\begin{align*}
QP &= (2A^2 - A + I)(A^3 - A) \\
&= 2A^5 - 2A^3 - A^4 + A^2 + A^3 - A \\
&= 2A^5 - A^4 - A^3 + A^2 - A
\end{align*}$

$PQ = QP$
"""

# ╔═╡ ac92fdd7-0b50-4bcc-8387-c76d74b70928
md"""
**Solution (c).**

$\begin{align*}
Qv &= (2{\lambda_0}^2 - 2\lambda_0 + 1) v \\
Pv &= ({\lambda_0}^3 - \lambda_0) v \\
PQv &= ({\lambda_0}^3 - \lambda_0)(2{\lambda_0}^2 - 2\lambda_0 + 1) v \\
\end{align*}$
"""

# ╔═╡ fdec34a2-462b-4e5b-af39-f76617f419ca
md"""
### PS29 #22

Let ``\lambda_1`` and ``\lambda_2`` be distinct eigenvalues of ``A``, with ``\lambda_2`` of multiplicity ``m_2``.
Let ``v`` be an eigenvector of ``A`` with eigenvalue ``\lambda_1``.
Can it be that ``v ∈ \text{Null}((A - \lambda_2)^{m_2})``?
"""

# ╔═╡ 55402cf2-623c-449a-9ab9-9ca057efa62a
md"""
**Solution.**

No, because the eigenvectors of ``A`` are linearly independent.

$(A - \lambda_2 I)^{m_2} v ≠ 0$
"""

# ╔═╡ f9e6cbe3-47bc-4f0a-a780-8000c8ad2a01
md"""
### PS29 #23
"""

# ╔═╡ ba60d592-8988-4709-b36e-bbdfd9d9e2a1
md"""
**Solution.**

$(A - \lambda_1 I)[(A - \lambda_1 I)^k v] = 0 \implies (A - \lambda_1 I) w = 0$
"""

# ╔═╡ c42fc532-f9d8-4c1f-bc65-3865fc0ec84e
md"""
### PS29 #24

Let ``A`` be an ``n × n`` matrix, let ``\lambda_1`` and ``\lambda_2`` be different eigenvalues of ``A``.
Let the multiplicites of ``\lambda_1`` and ``\lambda_2`` be ``m_1`` and ``m_2``, let ``\mathcal{E}_{\lambda_1}`` and ``\mathcal{E}_{\lambda_2}`` be the respective **generalized** eigenspaces.
Suppose ``v_1 ∈ \mathcal{E}_{\lambda_1}`` and ``v_2 ∈ \mathcal{E}_{\lambda_2}`` are nonzero.
Take advantage of Problem 22 to show:

> If ``v_1 ∈ \mathcal{E}_1`` and ``v_2 ∈ \mathcal{E}_2`` are nonzero vectors, then ``v_1`` and ``v_2`` are linearly independent.
"""

# ╔═╡ 468eb201-2ae0-4a35-9a43-86c6c7b45c43
md"""
### PS29 #25

Let ``V`` be the span of the functions

$v_1 = e^t \sin{2t}, \quad v_2 = e^t, \quad v_3 = e^t \cos{2t},$

let ``T : V → V`` be the operator ``T(f) = f'``.

(a) Find the matrix of ``T`` with respect to this basis of ``V``, then the characteristic polynomial of that matrix.

(b) Find the roots of the characteristic polynomial.
(One is real, the other two are complex.)
"""

# ╔═╡ 25e59fa9-e644-48ea-b31a-f548f29a3759
md"""
**Solution (a).**

$\begin{align*}
T(v_1) &= T(e^t \sin{2t}) = e^t \sin{2t} + 2e^t \cos{2t} = v_1 + 2v_3 \\
T(v_2) &= T(e^t) = e^t = v_2 \\
T(v_3) &= T(e^t \cos{2t}) = -2e^t \sin{2t} + e^t \cos{2t} = -2v_1 + v_3
\end{align*}$

$A = \begin{bmatrix} 1 & 0 & -2 \\ 0 & 1 & 0 \\ 2 & 0 & 1 \end{bmatrix}$

$\det(\lambda I - A) = \begin{vmatrix} \lambda - 1 & 0 & -2 \\ 0 & \lambda - 1 & 0 \\ 2 & 0 & \lambda - 1 \end{vmatrix} = (\lambda - 1)((\lambda - 1)^2 + 4) = (\lambda - 1)^3 + 4\lambda - 4$
"""

# ╔═╡ 559403e2-1119-4dea-97bd-3dc6dd8bcce8
let
	@variables λ
	A = [1 0 -2; 0 1 0; 2 0 1]
	expand(det(λ*I - A))
end

# ╔═╡ 908cedd3-0113-4c2a-9f56-aec2bf5ded66
md"""
**Solution (b).**

$\begin{array}{c|cccc}
1 & 1 & -3 & 7 & -5 \\
& & 1 & -2 & 5 \\
\hline
& 1 & -2 & 5 & 0
\end{array}$

$p_\lambda(A) = (x - 1)(x^2 - 2x + 5)$

$\lambda = \frac{2 ± \sqrt{4 - 4(5)}}{2} = \frac{2 ± \sqrt{-16}}{2} = \frac{2 ± 4i}{2} = 1 ± 2i$

$\lambda_1 = 1 \qquad \lambda_2 = 1 + 2i \qquad \lambda_3 = 1 - 2i$
"""

# ╔═╡ 019d3f9f-2f50-4aaf-8ba7-dc20cc7e1f4f
md"""
### PS29 #26

Let ``V`` be the span of the functions

$v_1 = e^t \sinh{2t}, \quad v_2 = e^t, \quad v_3 = e^t \cosh{2t},$

let ``T : V → V`` be the operator ``T(f) = f'``.

(a) Find the matrix of ``T`` with respect to this basis of ``V``, then the characteristic polynomial of that matrix.

(b) Find the roots of the characteristic polynomial. (All three are real.)
"""

# ╔═╡ 5920ba54-6e9e-495e-a0c8-dfc611e19b5b
md"""
**Solution (a).**

$\begin{align*}
T(v_1) &= T(e^t \sinh{2t}) = e^t \sinh{2t} + 2e^t \cosh{2t} = v_1 + 2v_3 \\
T(v_2) &= T(e^t) = e^t = v_2 \\
T(v_3) &= T(e^t \cosh{2t}) = 2e^t \sinh{2t} + e^t \cosh{2t} = 2v_1 + v_3
\end{align*}$

$A = \begin{bmatrix} 1 & 0 & 2 \\ 0 & 1 & 0 \\ 2 & 0 & 1 \end{bmatrix}$

$\det(\lambda I - A) = \begin{vmatrix} \lambda - 1 & 0 & 2 \\ 0 & \lambda - 1 & 0 \\ 2 & 0 & \lambda - 1 \end{vmatrix} = (\lambda - 1) ((\lambda - 1)^2 - 4)$
"""

# ╔═╡ 0f5552d5-9a9b-4dd8-a6ea-c000b4868abe
let
	@variables λ
	A = [1 0 2; 0 1 0; 2 0 1]
	expand(det(λ*I - A))
end

# ╔═╡ 30f40916-d34c-4903-b78f-35e73cc25d5d
md"""
$\begin{array}{c|cccc}
1 & 1 & -3 & -1 & 3 \\
& & 1 & -2 & -3 \\
\hline
& 1 & -2 & -3 & 0
\end{array}$

$p_A(\lambda) = (\lambda - 1)(\lambda^2 - 2\lambda - 3)$

$\lambda = \frac{2 ± \sqrt{4 - 4(-3)}}{2} = 1 ± 2$

$\lambda_1 = -1 \qquad \lambda_2 = 1 \qquad \lambda_3 = 3$
"""

# ╔═╡ a6a0dd74-3648-4dd7-b226-50aa35f8f75b
md"""
### PS29 #27

Let ``V`` be an inner product vector space of dimension 3.
In a certain basis ``v_1, v_2, v_3`` of ``V`` the inner product gives

$\begin{bmatrix}
⟨v_1,v_1⟩ & ⟨v_1,v_2⟩ & ⟨v_1,v_3⟩ \\
⟨v_2,v_1⟩ & ⟨v_2,v_2⟩ & ⟨v_2,v_3⟩ \\
⟨v_3,v_1⟩ & ⟨v_3,v_2⟩ & ⟨v_3,v_3⟩ \\
\end{bmatrix} = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 2 & 1 \\ 0 & 1 & 3 \end{bmatrix}.$

Find:

(a) ``\|v_2\|^2``,

(b) ``\|v_1 + v_2\|^2``,

(c) ``\|v_1 + v_2 - v_3\|^2``.
"""

# ╔═╡ eff0111d-48f8-40e3-a237-e70a9284fcd7
md"""
**Solution (a).**

$\|v_2\|^2 = ⟨v_2, v_2⟩ = 2$
"""

# ╔═╡ d6902368-c88b-445d-8dad-892f844f3bca
md"""
**Solution (b).**

$\begin{align*}
\|v_1 + v_2\|^2 &= ⟨v_1 + v_2, v_1 + v_2⟩ \\
&= ⟨v_1, v_1 + v_2⟩ + ⟨v_2, v_1 + v_2⟩ \\
&= ⟨v_1, v_1⟩ + ⟨v_1, v_2⟩ + ⟨v_2, v_1⟩ + ⟨v_2, v_2⟩ \\
&= 1 + 0 + 0 + 2 \\
&= 3
\end{align*}$
"""

# ╔═╡ df6c07ff-221a-42ff-88c7-4201c7b2e0da
md"""
$\begin{align*}
\|v_1 + v_2 - v_3\|^2 &= ⟨v_1 + v_2 - v_3, v_1 + v_2 - v_3⟩ \\
&= ⟨v_1, v_1 + v_2 - v_3⟩ + ⟨v_2, v_1 + v_2 - v_3⟩ - ⟨v_3, v_1 + v_2 - v_3⟩ \\
&= ⟨v_1, v_1⟩ + ⟨v_1, v_2⟩ - ⟨v_1, v_3⟩ + ⟨v_2, v_1⟩ + ⟨v_2, v_2⟩ - ⟨v_2, v_3⟩ \\
&\quad- ⟨v_3, v_1⟩ - ⟨v_3, v_2⟩ + ⟨v_3, v_3⟩ \\
&= 1 + 0 - 0 + 0 + 2 - 1 - 0 - 1 + 3 \\
&= 4
\end{align*}$
"""

# ╔═╡ 743ec978-23d2-40ec-aaee-705daf557c92
md"""
### PS29 #28

Let ``V``, the basis, and the inner product be as in Problem 27.
The subspace of ``V`` orthogonal to ``v_3`` consists of all vectors ``v = a_1 v_1 + a_2 v_2 + a_3 v_3`` such that

$⟨a_1 v_1 + a_2 v_2 + a_3 v_3, v_3⟩ = 0.$

This is a condition (a linear equation of rank 1) on the coefficients: one equation in three unknowns.
Find two independent solutions and use them to write a basis of ``\{v_3\}^\perp``.
"""

# ╔═╡ 9095f049-be4f-4d14-923c-4c0a41b06837
md"""
**Solution.**

$\begin{align*}
⟨a_1 v_1 + a_2 v_2 + a_3 v_3, v_3⟩ &= 0 \\
⟨a_1 v_1, v_3⟩ + ⟨a_2 v_2, v_3⟩ + ⟨a_3 v_3, v_3⟩ &= 0 \\
a_1 ⟨v_1, v_3⟩ + a_2 ⟨v_2, v_3⟩ + a_3 ⟨v_3, v_3⟩ &= 0 \\
a_2 + 3a_3 &= 0
\end{align*}$

$\begin{bmatrix} 0 & 1 & 3 \end{bmatrix}$

$\implies \begin{bmatrix} a_1 \\ a_2 \\ a_3 \end{bmatrix} = \begin{bmatrix} a_1 \\ -3a_3 \\ a_3 \end{bmatrix} = a_1 \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix} + a_3 \begin{bmatrix} 0 \\ -3 \\ 1 \end{bmatrix}$

Basis of ``\{v_3\}^\perp`` is ``\left\{v_1, -3v_2 + v_3\right\}``
"""

# ╔═╡ 15b1e91f-1a76-478e-b44d-b48cf19dd6d7
md"""
### PS29 #29

Let ``V`` be the set of polynomials in the variable ``t`` with real coefficients and degree ``≤ 2``.
For ``p, q ∈ V`` define

$⟨p, q⟩ = ∫_0^1 p(t) q(t) \;dt.$

This is an inner product.
Use Gram-Schmidt to orthogonalize the basis ``t^0, t^1, t^2`` of ``V``.
"""

# ╔═╡ 55a0390c-0243-411c-bead-74a6c112a8a8
md"""
**Solution.**

$\begin{align*}
w_0 &= t^0 = 1 \\
w_1 &= t - \frac{⟨w_0, t⟩}{⟨w_0, w_0⟩} w_0 = t - \frac{⟨1, t⟩}{⟨1, 1⟩} = t - \frac{1}{2} \\
w_2 &= t^2 - \frac{⟨w_0, t^2⟩}{⟨w_0, w_0⟩} w_0 - \frac{⟨w_1, t^2⟩}{⟨w_1, w_1⟩} w_1 \\
&= t^2 - \frac{⟨1, t^2⟩}{⟨1, 1⟩} - \frac{⟨t - 1/2, t^2⟩}{⟨t - 1/2, t - 1/2⟩} (t - 1/2) \\
&= t^2 - \frac{1}{3} + \frac{3}{12} (t - 1/2) \\
&= t^2 - \frac{1}{3} + \frac{1}{4} t - \frac{1}{8} \\
&= t^2 + \frac{1}{4} t - \frac{11}{24}
\end{align*}$

$⟨1, 1⟩ = ∫_0^1 1 \;dt = \left.t\right|_0^1 = 1$

$⟨1, t⟩ = ∫_0^1 t \;dt = \left.\frac{t}{2}\right|_0^1 = \frac{1}{2}$

$⟨1, t^2⟩ = ∫_0^1 t^2 \;dt = \left.\frac{t^3}{3}\right|_0^1 = \frac{1}{3}$

$⟨t - 1, t^2⟩ = ∫_0^1 (t - 1) t^2 \;dt = ∫_0^1 (t^3 - t^2) \;dt = \left.\frac{t^4}{4} - \frac{t^3}{3}\right|_0^1 = -\frac{1}{12}$

$⟨t - 1, t - 1⟩ = ∫_0^1 (t - 1)(t - 1) \;dt = ∫_0^1 (t^2 - 2t + 1) \;dt = \left.\frac{t^3}{3} - t^2 + t\right|_0^1 = \frac{1}{3}$
"""

# ╔═╡ 146403f6-0f16-4c28-b80c-77253eaaf778
md"""
### PS29 #30

Let

$A = \begin{bmatrix} 3 & 1 \\ 1 & 2 \end{bmatrix}$

(a) Verify (with the standard dot product of ``ℝ^2``) that

$\begin{bmatrix} x_1 \\ x_2 \end{bmatrix} ⋅ \left(A \begin{bmatrix} y_1 \\ y_2 \end{bmatrix}\right) = \left(A \begin{bmatrix} x_1 \\ x_2 \end{bmatrix}\right) ⋅ \begin{bmatrix} y_1 \\ y_2 \end{bmatrix}$

(b) Find the eigenvalues of ``A``, and for each of them (they are different) find an eigenvector.

(c) Let ``v_1`` and ``v_2`` be an eigenvector for different eigenvalues.
Find ``v_1 ⋅ v_2``.
"""

# ╔═╡ 195957b7-017e-46be-812a-2bdcadd408e5
md"""
**(a)**

$x ⋅ Ay = (A^T x) ⋅ y$
"""

# ╔═╡ 59d8d68c-8630-4a01-97a3-68541bae712a
md"""
**(b)**
"""

# ╔═╡ 773b603d-a895-4e37-ba17-b1fb3d0877ee
let
	@variables λ
	A = [3 1; 1 2]
	expand(det(A - λ*I))
end

# ╔═╡ 2f09b565-3e68-4985-aaa0-7d27f0f491ea
let
	a = 1
	b = -5
	c = 5

	(-b - sqrt(b^2 - 4*a*c)) / 2a, (-b + sqrt(b^2 - 4*a*c)) / 2a
end

# ╔═╡ 671293d3-d152-4635-9d49-55c17d84de79
md"""
**(c)**
"""

# ╔═╡ 90261640-dece-4136-8997-7a5560570d32
let
	A = [3 1; 1 2]
	rref(A - eigvals(A)[1]I), rref(A - eigvals(A)[2]I)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
RowEchelon = "af85af4c-bcd5-5d23-b03a-a909639aa875"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
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
git-tree-sha1 = "f885e7e7c124f8c92650d61b9477b9ac2ee607dd"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.1"

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

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "3cc368af3f110a767ac786560045dceddfc16758"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.3"

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

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

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
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

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
git-tree-sha1 = "def0718ddbabeb5476e51e5a43609bee889f285d"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.0"

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
# ╟─35bb1e85-5c56-48e5-a8ea-92ffb0ca6148
# ╟─0727e9d6-c906-410c-b7c3-b7f2922b3309
# ╟─2262a111-83be-4e51-8518-ef50b26fe450
# ╟─854c658c-c1ba-4dac-8982-17513b89f1df
# ╟─2468ff5b-b17b-41b5-983d-46d4229c2dc8
# ╟─db250a73-76e6-4288-945f-250e0a200df3
# ╟─965d09db-1fed-4706-9d7e-2b862a3ce44a
# ╟─7db67aff-3077-43d4-8177-6f5688b86eae
# ╟─0df73e47-63fb-46ba-ab1f-68e6f504cab3
# ╟─6b00b2de-35da-452b-953d-810f89e5fc3f
# ╟─6ac56c3d-8d4c-4a1c-a182-2c97e020ddc2
# ╟─4a056063-57d1-4596-b6c3-8892c6a78054
# ╟─0a904cfb-c4a9-48f8-8829-64a45fe362d9
# ╟─9c0518db-b2cd-4499-82bc-3735826373e7
# ╟─f2a2e166-6257-4efc-be59-636099c68740
# ╟─9c86875e-6070-4601-b24b-3da1a0e79136
# ╟─48b71ef8-d6ce-45ad-928a-2a2ca0e8ab5f
# ╟─8635836b-5e67-4f51-8501-881c358f08a4
# ╟─3f1cf813-9ea9-41fa-b5a5-a13343f3b20c
# ╟─9a1450a0-d356-42dd-806a-0bd6fc5e40de
# ╟─12dfce5c-39ce-46e3-878e-324ec03f73af
# ╟─05811194-b5d5-4cf8-af9d-37df311633cc
# ╠═20e520b8-bcdb-4e3e-aa22-21042e6ae15d
# ╠═d5bf5f56-3d4c-4afd-9c9d-b9f6a07ec145
# ╟─4f02f2d3-f72b-433f-ba23-6291ed773ddc
# ╟─59a6c778-3bc8-4462-abc7-1e5bdc1fd8c4
# ╟─b3d9f2fe-afe0-4e85-9216-354386bbdc51
# ╟─cb72d5be-289d-4e5a-aad1-360c716026b3
# ╠═59df2f01-4912-4da7-a349-35447eac202b
# ╟─1927374f-49d4-4d8c-97f3-e2ca80c945e3
# ╟─bd563257-69ea-41ce-bcc6-3c1cc62a0d8d
# ╠═e7dc3822-653a-4462-9e7a-1fde12cf638d
# ╟─f6bc8530-1822-44a7-b0ca-d39b1efd9225
# ╟─e42f0f08-fbbe-4547-a970-cc4204f447b1
# ╟─9f527f88-244a-4f60-8c45-fcf37b760656
# ╟─8b056c30-1e7c-4882-9a0c-f23427509f97
# ╟─0ebc87b1-5187-4106-9c95-6cfb75fc97a7
# ╠═07d9902a-af9f-474a-a015-56a7ce5b2476
# ╟─4b780acb-6fe8-446e-a0d8-293f520dd98c
# ╟─3c3a4efd-21cc-4d78-b86b-da0118a9eaff
# ╟─97ca197b-7d1f-4675-acd7-cd622946a24e
# ╠═53fdf8ed-9a96-41d8-81a8-7c5fa91bb00c
# ╟─12d3e4b5-b0df-4b90-baf7-3ec9a319d153
# ╟─133cf267-aa8d-46a6-a1f7-4d4a833fad75
# ╟─50a6df63-9753-46fc-bb20-d6b489108197
# ╟─a1409499-b628-4506-a5de-c3947f965433
# ╟─5bd95da9-b19b-4502-8f86-e81ec8af4acf
# ╟─904607da-67c6-4f45-a2bb-2e8d58514830
# ╟─509a9d05-fb66-49a8-945a-9b020627f745
# ╟─eb31d37b-ac26-4d04-8598-097dd16b332e
# ╟─a25c703f-537b-4cfb-ad60-559d027669a7
# ╟─080dd720-d75f-4108-aa22-cfe946675b41
# ╟─e9b8568e-a895-48fc-aa32-fabb03e574f2
# ╟─01f57e22-3edb-4a13-8d8d-24ebc30d4878
# ╟─4cb458b2-67a6-433c-8286-4b276f52a183
# ╟─6487ad98-683c-429b-a0d2-a370cc0b0c83
# ╟─7051900e-7863-45cd-b1dd-7066793b8360
# ╟─39204db3-380a-443f-b974-8bc652118ddf
# ╟─f976f042-3cc6-4f50-88c4-7dcd488d7cd5
# ╟─2bc90945-8e1f-47ca-a519-f8d93ddacbcf
# ╟─26dbb5a5-0005-40ea-bb88-963b4b995520
# ╟─7888fa3b-3074-4ccd-8f64-f7022e32aec7
# ╠═b5c21a12-d081-4686-ba67-6cd9256cc09e
# ╟─5135f07a-a843-4751-9c54-35712ab2ac04
# ╟─57bd2d47-b14c-4459-9c10-7d8d71705254
# ╠═4c41085c-3443-4716-bba2-7839ea337c40
# ╟─ea7c75d2-eb5e-4d51-8ee9-dee330304b6b
# ╠═fb5b591a-f32c-4bfe-a2cb-8183ba942a65
# ╟─8677382d-34c3-4de1-8894-6028380d5b02
# ╠═698abde3-1f79-4b28-9ffa-224972ba388a
# ╟─bd10b677-7b37-41ef-a4a0-3ad1017e5acd
# ╟─88d2cf92-46fa-4f9b-aab9-9a2fb8111d1c
# ╟─329ecc56-67a8-4ffb-8044-7dfca93a8116
# ╟─a29e9713-6550-43ee-95eb-eb35d46e51cc
# ╟─2415a79f-b72c-41d4-9784-0dec68ee4423
# ╠═f25e481c-3639-4980-947a-83128a148a81
# ╟─cf5ac8d4-d741-4946-a7a6-1a3bed5304d8
# ╟─ad6a34d8-9bd1-4105-8afd-f19e5328933a
# ╠═29860b91-212e-4b88-b199-af479a132c01
# ╟─b0927f96-ec7b-4a96-8a10-d69b55ec2c39
# ╟─590ba8f2-327c-43cf-8ea6-79bd2bccbbfe
# ╠═b5de1268-62a9-46c0-84b4-fb91160d729f
# ╟─6aab5a5e-2e68-4e1b-b4aa-daaf6d4482bf
# ╟─06cbfac7-ad60-4e22-a886-923ffebf04ad
# ╠═820b0411-9647-455d-badf-804886f89dc2
# ╠═26c08822-530b-4d78-9a27-381baa1f4bea
# ╟─901848a7-6f10-4730-91f8-aaeadb132b2a
# ╟─def19b3e-9ec8-4742-b378-415a3df7107d
# ╟─a96cd412-dccf-45bb-89cf-05b24774b9c2
# ╠═88dee7bb-53c3-4041-9712-ca66e6624501
# ╟─b8744e45-4c20-4fcc-911a-e5c57c38c2ef
# ╟─49a8acd9-0db5-4cb2-9cf5-0437483dbe71
# ╠═4c1e7119-db6d-413c-ad60-9fdfa0f58d18
# ╟─12693645-0282-46f3-9a4a-02c122f597fc
# ╟─16e43f06-127b-40d4-a26a-3e5a1e12a847
# ╠═37b94d4a-bd8b-412b-b20b-b4c8453b6d41
# ╟─b18830f5-dc83-4160-b87e-4212fd752175
# ╟─63db7ceb-b710-4830-acf0-61ede14c2fe5
# ╟─3eeb4786-1d9c-4cd3-8d2c-7910418f71bd
# ╟─77994a8b-d5e0-4673-899f-d16adefba744
# ╟─9b0bf121-c4a7-4dfc-a60e-29b8f616940b
# ╟─0dc779f0-3d01-4252-a5fd-3ba9d8d0ddd4
# ╟─8a798a6e-46b7-498e-be2f-253a45b13900
# ╟─b1a4e8f3-a5a0-478a-bd85-137e8cf4cd59
# ╟─32662135-aad8-47e1-bc33-976a3bf5445f
# ╟─eadcfe72-cf6c-4d97-855e-be9dac1a0b5b
# ╟─2775318a-9d4f-412a-9c5c-fbb698174070
# ╟─e985723e-144d-480a-bbd3-1eebe1331347
# ╠═a1ce5603-defc-4a95-9943-f38fbf40dfb5
# ╟─4455c750-19f2-4628-9201-60c18aecb291
# ╟─63020d78-c97a-4212-b729-87305143d597
# ╟─67aa6210-dd6b-4437-be6f-91e021eb75ef
# ╠═cf46c827-b872-42ac-b06b-84d93ce28970
# ╟─f1ec3d56-f091-4d3f-a992-afb994e0f2a7
# ╠═358f9c4a-ef4b-478b-b1ae-3f1f42d1da78
# ╟─75098ed4-b3b6-41df-b3ce-b4c92bf1f091
# ╠═253187a1-8942-4ce3-ba7c-c3e1b5275828
# ╟─f685369a-9768-496a-be6f-f006b70524c5
# ╟─1ff64173-6746-4029-8014-7afcaddb2e71
# ╠═29da4f6f-ce23-4092-8aaf-a6af74856b6b
# ╟─c5e60741-c889-41e9-86ae-f3d301f6bb21
# ╟─d9970466-8e57-4e98-86ac-114a9450df7b
# ╠═a234609f-b53d-492c-a34c-c264d1e14bab
# ╟─58c6e499-1ff5-42c6-9fb8-02708f56427b
# ╟─6a8cefbb-60d8-45d4-8d91-3f6896ed6ffa
# ╠═30993fa2-4fe8-4a65-b642-d53dbd0fcbc5
# ╟─ac6f7409-9feb-426e-b82a-27333ffaeecd
# ╟─d37e2898-05f0-4855-a81d-b06ab113c541
# ╠═b8365990-a6a1-433e-83bf-4950d88c61e6
# ╟─363a9d0f-4e24-4a48-ba34-463e50a5078b
# ╟─b3fb8b64-83eb-44b0-b721-63d6d28637ab
# ╟─1f2c250c-fd6c-435d-9548-99a11e8d354a
# ╟─cdc02ef6-fe47-45aa-a0f3-32a2f8a9f952
# ╟─c924378b-0d35-4a50-a804-6e03adc1bde6
# ╟─bb6f8b89-c199-4393-a7ee-ab5b06d6ad6a
# ╟─6e6fdde8-36eb-44a1-855e-21b676c6ece8
# ╟─065df513-1ffa-47eb-88ae-a69a834dcac5
# ╠═d070d235-c7fc-4082-b1f5-0d2cd4c59ea0
# ╟─4cfd479d-4c64-4351-9324-9db91c5a3178
# ╠═bf19bf1a-7801-4ffd-96fa-b6d5c281b55d
# ╟─52541443-187c-4ef3-a24e-3142d7ad9bbf
# ╟─c2532725-32ea-40e6-9d24-edd0c4b9bda2
# ╠═eef85d24-92ce-492d-87d5-40674a28d2ee
# ╟─a6ec0c87-245d-4548-8199-5438f3911433
# ╟─c77df6f8-feea-4f45-8a91-aa307332ec5f
# ╟─dbe32a98-086b-485a-8999-22aed04b3d26
# ╟─b24cb37a-dfc2-4ca8-ab18-4e39a6054976
# ╟─60449475-2834-423c-9620-bd1c001c6cca
# ╠═f3afa449-4f62-4e4f-b57e-43eb0090a65b
# ╟─d882dd0d-7ba8-40e6-9a56-410948ac1e7f
# ╟─9c2cb2f9-6151-49e3-90db-7e6d9c95e6af
# ╟─79b287ae-bbcc-40bf-b667-b6214aa6ec30
# ╟─03ab5746-802e-424c-97b0-a679290b0f3a
# ╟─1dd87c8f-790d-4cc6-8a2d-3118e7c08901
# ╟─02609b9a-c8f4-4f38-96e0-acef86136c81
# ╟─f26db46a-d1db-46f5-9e71-b0708834c434
# ╟─e810e3e4-256f-481e-9061-8c1d67c28a7e
# ╟─e782342f-b07c-4bec-8d9a-5a7a9effff3a
# ╟─4ac5a171-a82a-4a49-8ffd-43ba3bb71131
# ╟─31e5b3c4-1f85-412f-8c0b-eb351abb3eb9
# ╟─1b28eba2-7d10-41bd-a00b-a970369aa37a
# ╠═98663ba6-3502-443b-9d9c-f32cddb9bd74
# ╟─b6ae668c-aa4c-4ef0-82ab-2b1a9669997e
# ╟─d3c3260b-a38e-4e24-a1f6-8e42ca7b74ff
# ╟─c5d083b9-0340-4a17-9b07-35638fbc5cef
# ╟─0dc04569-51db-472b-835c-3d288c3d37ef
# ╟─ee7f6086-c337-47d8-b86f-a320f7596ee6
# ╟─b659e17b-7ae7-4297-83e7-f29b4e8568ec
# ╟─2cd4c040-ed9e-48a8-bf88-83e00e1b1978
# ╟─819b1e3e-77e9-45a6-83cf-60ee4c3af2c1
# ╟─796c2bf5-53c7-4b55-9584-bfcc4e09482a
# ╟─5c06aca8-1ec1-4b06-89e8-fd32a3be7e46
# ╟─877a5b7c-83b2-45d7-af18-e54c45cc2bbe
# ╟─f0ec76f1-4d0e-4d70-b1cc-a31efc6834aa
# ╟─ed28d412-d601-4852-bb44-e93a029cb3c4
# ╟─2eba100d-20c0-4b80-89a9-18f113f88cc3
# ╟─9daaf6fc-a3fa-4914-9f0e-2ce0f5df2004
# ╠═6a4bb909-caf0-4b96-ba15-7c443499cbda
# ╠═403b70ab-f14f-4725-8434-2b9e644955b5
# ╠═5cbb9223-33b0-406c-9391-2fa933440df7
# ╠═e845cbeb-37d6-464f-989c-f7adf6b8d6fb
# ╠═abeb3c03-29a2-44c1-b542-2af8eb35b8cc
# ╟─ea442e30-106c-4f4b-9cc5-594dc9bd18a7
# ╠═06ddaef2-dc36-4b11-b5a5-d78b72991ff4
# ╠═54ddb1b0-b0a9-446c-a814-a551b0901bd1
# ╠═4ede0ef8-8bc9-4a5f-820a-05e00b3bf82a
# ╠═25d8aa91-7e21-48d6-871c-0f3aae09f07f
# ╠═cbc18af4-5a2c-4971-9939-5a757f90033c
# ╠═fb17c046-b1dc-46a3-9dd1-c7fe5792efd3
# ╠═5b5ac2cc-ed80-4adf-9b02-016586ec5a68
# ╠═e9caa6fa-49bb-4816-9b78-84a7c43e93c0
# ╠═54c05e91-35ef-44c4-b543-193ca2e98d1f
# ╠═ee649240-cd68-480c-b055-2ff6e9fd0d9a
# ╠═516a0f4a-dcf5-4e92-8f91-5a1cb3336a5c
# ╠═0451fa11-85cd-4e3c-b084-83e04a1276b3
# ╠═f8a4533c-98e3-4ecf-8514-a081fc51da7a
# ╠═ba4fd86b-2610-45ca-bb54-487b3029c982
# ╠═ce0942ab-074b-46e7-b4c2-e217022d68f7
# ╠═b00d86b3-2d9a-4478-98bb-f09ca8d39d72
# ╠═d3374ce9-e308-4d93-9577-2b50e16805b7
# ╟─583a73d4-5e64-4e6f-958e-74590ce21ce8
# ╟─ff1cc355-2f53-43c2-860b-5617c206a2ec
# ╟─6427d602-c735-4f30-b41f-f102fcf7fa89
# ╠═3c5bb1e4-1d9f-4e52-97b2-5cc41ab2c106
# ╟─1af14435-984c-49b2-b088-ff4505d10bd1
# ╟─68a2b285-fd90-4225-af15-10f54552c9f8
# ╠═63a22534-b547-4a50-ac80-197b279f4bd9
# ╟─e31c86fe-713d-4301-90d4-d966d79f10e4
# ╠═a1bbacbe-7683-4b0c-aac0-e60c2ad95e8f
# ╟─6b6573e4-2f61-438d-9763-b3a4200ff790
# ╟─96fd1a78-be44-4e26-950b-2efb1b78973b
# ╟─2069b9b9-961e-4665-b2a4-5badf66a99d8
# ╠═18620a95-56db-4328-a3d4-b0343344d0b4
# ╟─130e838f-c6ac-4b54-9288-352e855de012
# ╟─b0e4d322-8b8b-489d-bc85-d8d95bba024f
# ╟─ada1b995-c6e1-4f78-9a20-443adc40dc89
# ╟─a0f249d1-2425-44b6-9454-784175d12657
# ╟─c90450e4-2d02-435b-a7b3-5b9316f4d0a2
# ╟─fe291037-01a5-4fba-978a-13c4af231aa7
# ╠═ecfab479-361c-41a4-b304-433702da82dd
# ╟─fecee897-5961-4042-82d8-0a9facda5359
# ╟─616f4cbd-cb7d-4042-b080-3a7e3f4a8dcf
# ╠═708821d9-f22b-4765-88bb-72b88b57d87b
# ╟─6f7c5fbf-b0d2-4912-afc6-65f64e2ebb7e
# ╟─2f665e7a-daff-4320-b570-825f4e8fcc66
# ╟─8ef3cd15-a03c-42f0-aa71-08b512fa2346
# ╟─654ca869-fd5f-468f-a19b-ab4f7bfd6733
# ╠═d0767f8c-28a1-4c7d-a261-1c5950f39c2a
# ╟─9783c6b9-a4cf-49c7-bdb5-b253de3a3f08
# ╟─f836b6ae-f4dd-4799-9198-d79e4a48e8d6
# ╟─45088ee2-ef86-481b-8ef0-0bf2a2482a38
# ╟─bee81924-5aa0-4ffd-afe1-0a52a70067af
# ╟─40622c8c-1bda-41cd-8050-9347e003890d
# ╟─39166607-adb8-43ca-a300-63db130377fc
# ╠═4e845087-44d1-41f9-92aa-9e0298015993
# ╟─13c5e648-f3ea-4cde-8d7b-e9f714ae1d2a
# ╟─baf41a50-7f3e-4a07-aef6-561fb5c127d6
# ╟─a2c6a412-6e24-4ae8-a652-c060c188f8b4
# ╟─f34e13d9-35ef-4d5d-98cb-e383d5ff9bc8
# ╟─2dc1ee27-399b-476a-9a29-dc8b9d0919fd
# ╠═3bc626df-b8e8-463c-ada4-e0cc139891aa
# ╟─88dbca07-05c4-4759-9cda-4be5c4ac0d30
# ╟─8740405a-7a5c-4208-abf6-07dd00b703e4
# ╟─5b52a6e4-38eb-488d-9b9e-c76ec24e13d2
# ╟─4692a6e4-5808-40bf-b361-e8202cd97f6f
# ╠═565df250-b729-40d5-bba9-162307022d22
# ╟─b3d4322e-001d-41ca-bd18-f66580699506
# ╟─a8252a62-81a0-4739-b140-747714eaa880
# ╟─888a0de2-836d-4b9c-963b-8725dcd50ab4
# ╟─8bb3d6ae-30fd-4ad5-a1b4-600ddfad2cb6
# ╟─de6a6773-d5e6-429c-9458-7ada3cb3f919
# ╠═dd9d5b05-15bc-4795-9d34-a2d23399daab
# ╠═b285f1a1-4c48-4a6d-a5eb-e013e7c5cae3
# ╟─6616a5f6-2b20-442e-bcc2-e34d94e684a2
# ╟─e4ef72e1-dd20-4537-a5f2-68e37ba0f26e
# ╟─0de75b5a-2216-45ac-a03a-3f9554e143cc
# ╟─7b363010-4b09-4912-b099-0290b127148a
# ╟─b84c139f-4919-4041-b742-ac41992d62b1
# ╟─ccc2e160-bbe6-4bdc-83ef-37d217023059
# ╟─85d46495-80a9-4715-a9ff-f15fdbb7d7f5
# ╟─80d68dd0-d5e8-4fb0-bada-67a34f654036
# ╟─5f3d4b7e-78da-4111-b929-ddeca6e9e774
# ╟─907d5b75-0fe3-4a98-a528-e4c72323584d
# ╟─319b0a44-59ed-475a-ae6f-ac5a68b1a7ed
# ╟─5d0be0ee-8987-44d5-aa64-ff434a91acc4
# ╟─f8998419-f2db-49a6-975f-873896a3197a
# ╟─6eb690f1-664f-499d-85b9-e5078860464c
# ╟─a43beff5-a615-4c88-9265-742c41ffa624
# ╟─da5aabe1-d5d1-4f93-8b6d-f9b31b0d8d89
# ╟─8956e768-4de6-4f69-b8f0-8e8396728e61
# ╟─d3636101-679f-4895-92d4-4957585f9cee
# ╟─3584afba-39bc-4a15-afcc-0c53b22ae2ae
# ╠═4eae4bdd-7103-44d3-8fda-4c3d8a665781
# ╟─6084e51c-caef-472a-a285-5aae4b2f11f1
# ╠═675b6c6f-7a0c-467a-99b3-347d810c9074
# ╠═923912b7-322f-4b2e-a3ec-46265f1c61c9
# ╟─758ea353-65ca-483a-9e44-d08e3789823b
# ╟─1bacfb53-3650-4272-9e2b-51f2fbf02c1a
# ╟─1c19cfdc-52c8-41a4-8347-eeabf853303e
# ╠═b58089bb-b178-4d5c-99cd-ea5fc254edf2
# ╠═72e3b4bb-5a09-44fa-84cb-3e442d599763
# ╟─61608c42-ab66-4d31-85a0-5ad9cff5d904
# ╟─f94be157-8409-42e6-bd2a-4076e2e98f4f
# ╟─98e27c9b-850c-41be-aac3-565bc8475687
# ╟─71a1cdbd-f12f-4ab6-83f2-b495dedaa953
# ╟─3618a2f0-c59c-4e07-8a1b-2d6f493612d6
# ╟─87e01065-9751-47ee-8144-d7aceb48cb9f
# ╟─d862a4c7-1bd3-4d2f-a1c7-717904f8061b
# ╟─f50a0569-4999-4ce1-9664-f9d9b6a26253
# ╟─c84211c9-a066-4120-9dd3-042af16f7044
# ╟─c998e12b-dc4b-4e6e-a457-3b543dbc45cd
# ╟─377a28c9-3245-4c2b-b026-ac631fd06528
# ╠═08903063-f791-4533-a65e-74e6229a96d2
# ╟─a127a412-8301-4e9d-bf1c-0ec403bdaf21
# ╟─21d1467f-8723-4b2c-b029-150751d247ce
# ╠═c765c755-cdae-4a18-88b2-a174bb1893ba
# ╟─15d25922-8ae5-4ab4-a842-90aa1efab208
# ╟─834c91c6-08e6-409b-8de7-6dd1ec5c6d00
# ╟─71585849-3a8e-4f26-a476-589969fd0999
# ╟─5e011536-a63c-480a-a80f-e5154d2e9c78
# ╟─13be3db5-3a54-43cf-bfac-7fc6ebba918d
# ╟─c2a6a003-4e2f-43d8-addf-1f3d6ba0a92f
# ╠═a8eb7f21-aab2-4e47-8d27-cc5257a6b4fd
# ╟─b4944c17-c72d-421c-9009-0ad3fdec10db
# ╟─5b5e5890-bf08-4ab7-bf13-ac3aa27adce9
# ╠═22ae14fb-1e25-4bcf-a49c-2bd9bab08a9d
# ╟─b67c846e-acb6-4a75-8051-d722d88796cf
# ╟─29554616-3b16-4d85-b910-f8875184d98c
# ╠═69407621-152a-46b9-8259-898347efcc2f
# ╟─390d2d81-16fc-4f19-bbb3-c3623c1fe8e0
# ╟─f27c2ade-291b-45b7-83e5-723a157d63f3
# ╟─27a55e85-10ec-4766-86b2-07d75541c55c
# ╟─79262f6a-607c-4dfc-92ec-464ac5905a94
# ╟─c5f1d778-bde7-4856-8cff-5929840b49f1
# ╟─f02e8a80-84af-4749-8ba6-ab0b379821c6
# ╠═4e4cf97c-c98c-475e-9748-d68b2c67b808
# ╟─9bdcee0d-08a4-41b8-b286-5997907d1d05
# ╟─0da3b8b5-b01d-417f-b462-d4372b91327b
# ╠═d11d84b9-5a3c-42b4-a849-c74c4d30e46d
# ╟─5591399b-db5e-45fd-834a-f3f92b89eb03
# ╟─d2ccaef6-4b98-407d-989a-9002232e2d57
# ╟─0e02b91d-43cd-49f6-9e7e-19f583481331
# ╟─b06ef30d-0ec4-4e38-a931-19bf3dc093e3
# ╠═e5c86d2a-5133-4b88-9818-62748974e749
# ╟─50eb9439-c478-4ece-85e3-2e49afc5c2be
# ╟─c8d7cb79-3466-49e4-8117-6eaed1a4902f
# ╟─5b31ed0e-bbc9-4901-9be7-10d38c9c8dc0
# ╟─75065d30-85eb-4ad9-9916-9c5ca36f9fc3
# ╟─3dadd448-00c7-4ad7-8651-42a28c41eaf6
# ╟─00c4a5ae-8c2b-43f5-a50a-3474998af970
# ╟─04267575-996a-4c0a-b992-ea85eff2ca46
# ╟─d935e71b-9357-49b0-b3bb-025dc8f8d7ce
# ╟─33023e36-c431-4a3a-8178-839a5ca2100d
# ╠═ac298ba7-8550-4b14-a4dd-3d0462370b4d
# ╟─06b44040-2288-4cc0-8ae9-b15ac15bfc42
# ╟─a7120b38-5dd8-4fac-8a89-9bd24ecdce96
# ╟─f4268677-e115-4232-bd7f-6e3379c355c3
# ╟─433d717b-2588-4f06-a8bf-dcbe91381c95
# ╟─da7dac05-c8e5-4792-98c1-ec133effef57
# ╟─d577293d-a943-45d7-9d5b-e95a292ecb59
# ╠═c668ebdd-3e0c-4c98-a336-f8933c564565
# ╟─84d6fe1a-649c-4668-a5f6-cba0cf73d0d5
# ╟─e7d69866-97fa-4914-a2db-31073f2bae12
# ╟─48cd9ff9-ba8e-4f13-9c1f-6620db4d793a
# ╟─56d69f1f-d55b-4d75-a291-7fb78eb92c2f
# ╟─5cb455fa-acb7-48a3-8c70-95a04c886c84
# ╟─d2836e79-c66a-4eef-90ed-bf669774b7f6
# ╟─5f41870b-5395-4642-bfc9-2a9aa7c47817
# ╟─2e160eac-4682-446c-a67c-d568718bac7b
# ╟─f6ed271c-8ace-432e-bb5a-5e04f4944629
# ╟─55e273db-4389-4a7e-84c8-bdb788c20b23
# ╟─c059e395-23cd-495c-98e7-1a820af5bf74
# ╟─d37ae0d5-79ae-478d-b84a-be068ec27f20
# ╟─a3a9f6d6-e29d-4810-9d29-760af87dd4da
# ╠═b21cc3bd-64e3-4690-a54a-563b66f93403
# ╟─dfc4e6fe-694d-4669-b962-8e00d7fc53c3
# ╟─b102a4b9-d81f-4bcd-893f-e39deb0a83d9
# ╟─4bad5e75-862b-40cc-bdb8-aeb1014ca125
# ╟─9e05f711-ecb6-4923-a6eb-8ac42d128c2c
# ╟─56fc50f7-88d0-46d6-abb3-de8a2cc58e12
# ╠═b8ed2d91-ec62-4077-b7e2-fa37b16fdb9c
# ╟─4de3766f-fc03-4a15-8e8c-4cea7eff3b35
# ╟─97e9c23b-b775-413b-960b-ac1c9e39f002
# ╟─ce5e862f-2799-40ef-aa3d-f4fde538119a
# ╟─64f1c0b6-dd27-4549-a894-7f69c873a0c1
# ╟─d01b32ed-4d8d-4a69-b6d0-69cea9beb81b
# ╟─64239289-a2f6-4d9e-b67d-df540a439ef7
# ╟─3c1a2ea2-d8a9-4b7e-8360-842d1c053529
# ╟─b29a70f1-e6a7-49cc-8800-9ae4d0b8e0ab
# ╟─c8feede5-ecb0-4f54-a95b-013f773c05fe
# ╟─b657c4e2-4d1b-466e-bc51-5146eb03f3e6
# ╟─ea1c0899-355e-49c9-a05d-28619faa978e
# ╟─b2d46287-301b-4a9a-b4e7-e1ff9d0d5aa7
# ╟─18170309-b97b-4ed3-ab69-3c2c640a0d9c
# ╟─cec5575a-4c2d-4fd9-8a44-7b5efb564a22
# ╟─f399418a-3ecd-4184-bde0-cdcc7734a442
# ╟─9dcd9b99-9935-4ae1-8069-6080d27549d0
# ╠═a001e88c-c34a-4c8e-883c-0bbec2dbbc6b
# ╟─41978ff8-56d6-40ac-9d05-1493355a43e6
# ╟─35c96ee5-e003-4f64-979e-b4d466032a11
# ╟─61a03d8e-bf2c-48fa-b841-b81474fdafcb
# ╟─1cde7920-ee1c-4877-b58d-a9dd41d590b3
# ╟─6bc1c224-ab48-4922-8c66-e20ababbda62
# ╟─5ca73e2a-87c4-4d8b-831f-7dea471accf9
# ╟─f5828a13-3c4d-49ac-8601-2f758b6bba84
# ╟─6444d766-8635-4b50-bc7c-7223f0353656
# ╟─8a199215-6456-4cac-827b-2ea3c793b1f7
# ╟─3f95062c-986a-486c-bc19-31b6e3c19ece
# ╟─15975b30-f690-4c03-ad2e-9c282417999f
# ╟─8c68ed02-6e2a-45a6-87dc-2e5416cfbd05
# ╟─da18d4c7-5f2b-453d-84fe-3ad32605089e
# ╠═9bf217de-c29e-4827-9c02-ffb916f2b2c9
# ╟─663a0508-fb15-444d-bd42-9e4434e92016
# ╟─b9a063bd-3a71-4666-babe-9eca110188db
# ╟─9b6d402e-5c19-4c71-949a-8ed0a107f82c
# ╟─a4f6faf5-f770-4af6-884a-f53f49dd73d5
# ╟─ac92fdd7-0b50-4bcc-8387-c76d74b70928
# ╟─fdec34a2-462b-4e5b-af39-f76617f419ca
# ╟─55402cf2-623c-449a-9ab9-9ca057efa62a
# ╟─f9e6cbe3-47bc-4f0a-a780-8000c8ad2a01
# ╟─ba60d592-8988-4709-b36e-bbdfd9d9e2a1
# ╟─c42fc532-f9d8-4c1f-bc65-3865fc0ec84e
# ╟─468eb201-2ae0-4a35-9a43-86c6c7b45c43
# ╟─25e59fa9-e644-48ea-b31a-f548f29a3759
# ╟─559403e2-1119-4dea-97bd-3dc6dd8bcce8
# ╟─908cedd3-0113-4c2a-9f56-aec2bf5ded66
# ╟─019d3f9f-2f50-4aaf-8ba7-dc20cc7e1f4f
# ╟─5920ba54-6e9e-495e-a0c8-dfc611e19b5b
# ╟─0f5552d5-9a9b-4dd8-a6ea-c000b4868abe
# ╟─30f40916-d34c-4903-b78f-35e73cc25d5d
# ╟─a6a0dd74-3648-4dd7-b226-50aa35f8f75b
# ╟─eff0111d-48f8-40e3-a237-e70a9284fcd7
# ╟─d6902368-c88b-445d-8dad-892f844f3bca
# ╟─df6c07ff-221a-42ff-88c7-4201c7b2e0da
# ╟─743ec978-23d2-40ec-aaee-705daf557c92
# ╟─9095f049-be4f-4d14-923c-4c0a41b06837
# ╟─15b1e91f-1a76-478e-b44d-b48cf19dd6d7
# ╟─55a0390c-0243-411c-bead-74a6c112a8a8
# ╟─146403f6-0f16-4c28-b80c-77253eaaf778
# ╟─195957b7-017e-46be-812a-2bdcadd408e5
# ╟─59d8d68c-8630-4a01-97a3-68541bae712a
# ╠═773b603d-a895-4e37-ba17-b1fb3d0877ee
# ╠═2f09b565-3e68-4985-aaa0-7d27f0f491ea
# ╟─671293d3-d152-4635-9d49-55c17d84de79
# ╠═90261640-dece-4136-8997-7a5560570d32
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
