### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 8c1fbb78-ec91-4027-ab81-740a097780fc
begin
	using Calculus
	using ForwardDiff
	using PlutoUI
	using Printf
end

# ╔═╡ 14954992-6e4b-11ef-1c5f-9184c4d14be8
A = [10.0 5.0; 2.0 9.0]

# ╔═╡ b2117c9f-2d47-4edf-a96a-132c5c323275
y = [6.0; 3.0;;]

# ╔═╡ ea96d72c-763e-4a00-acbc-e8686a87efcf
md"""
## Algoritmo de Jacobi

$x_k(i+1)=\frac{1}{\mathrm{~A}_{k k}}\left[y_k-\sum_{n=1}^{k-1} \mathrm{~A}_{k n} x_n(i)-\sum_{n=k+1}^N \mathrm{~A}_{k n} x_n(i)\right] k=1,2, \ldots, N$

The Jacobi method also can be written in the following matrix format:

$\begin{aligned}
&\mathbf{x}(i+1)=\mathbf{M x}(i)+\mathbf{D}^{-1} \mathbf{y}
\end{aligned}$

Com 

$\begin{aligned}
&\mathbf{M}=\mathbf{D}^{-1}(\mathbf{D}-\mathbf{A})
\end{aligned}$

E 

$\mathbf{D}=\left[\begin{array}{ccccc}
\mathrm{A}_{11} & 0 & 0 & \cdots & 0 \\
0 & \mathrm{~A}_{22} & 0 & \cdots & 0 \\
0 & \vdots & \vdots & & \vdots \\
\vdots & & & & 0 \\
0 & 0 & 0 & \cdots & \mathrm{A}_{N N}
\end{array}\right]$

"""

# ╔═╡ db8b3da7-1361-4f95-a0c4-c5a0f3729ff3
# Formulação do algoritmo de Jacobi sem formulação matricial
begin
	x = Float64[0.0; 0.0;;]
	x_ = copy(x)
	n = length(x)

	h = 0
	error = 1e6
	while (error >= 1e-4) && (h < 100)
		h = h + 1
		for k in 1:n
			a = 1.0 / A[k, k]
			b = y[k, 1]
			c = sum(Float64[A[k, i] * x_[i, 1] for i in 1:k-1])
			d = sum(Float64[A[k, i] * x_[i, 1] for i in k+1:n])
			
			x[k, 1] = a * (b - c - d)
		end
		
		error = maximum(abs.((x - x_) ./ x_))
		@printf "Erro na iteração %2i : %.4f. \n" h error 
		x_ = copy(x)
		
	end
	println("\nNumero de iterações: $h")

	println("Valor final de x: $x")
end

# ╔═╡ 4c55e5ef-d0bf-4ff5-af48-b8dc9fd3f5cd
A \ y

# ╔═╡ 677a125a-cf5f-49ce-8b09-94171a1e309e
md"""
## Algoritmo de Gauss-Seidel

$x_k(i+1)=\frac{1}{A_{k k}}\left[y_k-\sum_{n=1}^{k-1} \mathrm{~A}_{k n} x_n(i+1)-\sum_{n=k+1}^N \mathrm{~A}_{k n} x_n(i)\right]$

The Gauss-Sidel method also can be written in the following matrix format:

$\begin{aligned}
&\mathbf{x}(i+1)=\mathbf{M x}(i)+\mathbf{D}^{-1} \mathbf{y}
\end{aligned}$

Com 

$\begin{aligned}
&\mathbf{M}=\mathbf{D}^{-1}(\mathbf{D}-\mathbf{A})
\end{aligned}$

E

$\mathbf{D}=\left[\begin{array}{ccccc}
\mathrm{A}_{11} & 0 & 0 & \cdots & 0 \\
\mathrm{~A}_{21} & \mathrm{~A}_{22} & 0 & \cdots & 0 \\
\vdots & \vdots & & & \vdots \\
\mathrm{A}_{N 1} & \mathrm{~A}_{N 2} & \cdots & & \mathrm{A}_{N N}
\end{array}\right]$

"""

# ╔═╡ 3dec9587-a2aa-430d-a2bc-42bc26811eb2
# Formulação do algoritmo de Gaus-Seidel sem formulação matricial

let
	x = Float64[0.0; 0.0;;]
	x_ = copy(x)
	n = length(x)

	h = 0
	error = 1e6
	while (error >= 1e-4) && (h < 100)
		h = h + 1
		for k in 1:n
			a = 1.0 / A[k, k]
			b = y[k, 1]
			c = sum(Float64[A[k, i] * x[i, 1] for i in 1:k-1])
			d = sum(Float64[A[k, i] * x_[i, 1] for i in k+1:n])
			
			x[k, 1] = a * (b - c - d)
		end
		
		error = maximum(abs.((x - x_) ./ x_))
		@printf "Erro na iteração %2i : %.4f. \n" h error 
		x_ = copy(x)
		
	end
	println("\nNumero de iterações: $h")

	println("Valor final de x: $x")
end

# ╔═╡ bc894bd6-8de5-4392-9400-5bfafdc92200
md"""
## Exemplo de Não Convergência do Algoritmo de Gauss Seidel
"""

# ╔═╡ 6488b9a6-a040-4657-9aef-791553f6245e
let
	A = Float64[5 10; 9 2]
	x = Float64[0.0; 0.0;;]
	x_ = copy(x)
	n = length(x)

	h = 0
	error = 1e6
	while (error >= 1e-4) && (h < 20)
		h = h + 1
		for k in 1:n
			a = 1.0 / A[k, k]
			b = y[k, 1]
			c = sum(Float64[A[k, i] * x[i, 1] for i in 1:k-1])
			d = sum(Float64[A[k, i] * x_[i, 1] for i in k+1:n])
			
			x[k, 1] = a * (b - c - d)
		end
		
		error = maximum(abs.((x - x_) ./ x_))
		@printf "Erro na iteração %2i : %.4f. \n" h error 
		x_ = copy(x)
		
	end
	println("\nNumero de iterações: $h")

	println("Valor final de x: $x")
end

# ╔═╡ 8c3c6238-ceb6-4dc9-8e0a-87bceef5f785
md"""
## Método de Newton-Raphson
"""

# ╔═╡ b1de106c-0378-4ed4-865f-fde4a5cc531b
md"""
Exemplo de caso com função simples do tipo $f(x) = x^2$ e $y = 9$. Também considera-se $x_0 = 1.0$

Aqui utiliza-se a formulação clássica do método de Newton-Rhapson:

$\mathbf{x}(i+1)=\mathbf{x}(i)+\mathbf{J}^{-1}(i)\{\mathbf{y}-\mathbf{f}[\mathbf{x}(i)]\}$

Em que:

$\mathbf{J}(i)=\left.\frac{d \mathbf{f}}{d \mathbf{x}}\right|_{\mathbf{x}=\mathbf{x}(i)}=\left[\begin{array}{cccc}
        \frac{\partial f_{1}}{\partial x_{1}} & \frac{\partial f_{1}}{\partial x_{2}} & \cdots & \frac{\partial f_{1}}{\partial x_{N}} \\[0.1in]
        \frac{\partial f_{2}}{\partial x_{1}} & \frac{\partial f_{2}}{\partial x_{2}} & \cdots & \frac{\partial f_{2}}{\partial x_{N}} \\[0.1in]
        \vdots & \vdots & & \vdots \\[0.1in]
        \frac{\partial f_{N}}{\partial x_{1}} & \frac{\partial f_{N}}{\partial x_{2}} & \cdots & \frac{\partial f_{N}}{\partial x_{N}}
        \end{array}\right]_{\mathbf{x}=\mathbf{x}(i)}$

"""

# ╔═╡ 1012f4f9-e5c3-42f0-b2cb-a3c8d40b42ff
let
	f(x) = x^2
	J(x) = 2*x
	y = 9.0
	x = 1.0
	x_ = 1.0
	
	h = 0
	error = 1e6
	while error >= 1e-5 && h < 10
		h = h + 1

		x = x_ + 1.0 / J(x_) * (y - f(x_))
		
		error = abs((x - x_) / x_)
		x_ = x
		@printf "Erro na iteração %2i : %.3f. Valor de x = %.3f \n" h error x 
	end

	println("\nNumero de iterações: $h")

	println("Valor final de x: $x")	
end

# ╔═╡ a5eb2cde-3acd-436c-8adf-5255b3efdd2e
md"""
Exemplo de caso com sistema de equações não-lineares e com a formulação do método em quatro passos e sem a necessida de cálculo da matriz inversa:

$\mathbf{J}(i) \Delta \mathbf{x}(i)=\Delta \mathbf{y}(i)$

Em que:

$\Delta \mathbf{x}(i)=\mathbf{x}(i+1)-\mathbf{x}(i)$

e

$\Delta \mathbf{y}(i)=\mathbf{y}-\mathbf{f}[\mathbf{x}(\mathrm{i})]$


Then, during each iteration, the following four steps are completed:

STEP 1: Compute $\Delta \mathbf{y}(i)$ from $\Delta \mathbf{y}(i)=\mathbf{y}-\mathbf{f}[\mathbf{x}(\mathrm{i})]$ .

STEP 2: Compute $\mathbf{J}(i)$.

STEP 3: Solve $\mathbf{J}(i) \Delta \mathbf{x}(i)=\Delta \mathbf{y}(i)$ for $\Delta \mathbf{x}(i)$.

STEP 4: Compute $\mathbf{x}(i+1)$ from $\Delta \mathbf{x}(i)=\mathbf{x}(i+1)-\mathbf{x}(i)$.

"""

# ╔═╡ dd98c4ac-6741-4712-a30f-ec9ae4be7b97
let
	f(x) = [x[1] + x[2]; x[1]*x[2]]
	y = [15.0; 50.0;;]
	x = [4.0; 9.0]
	x_ = copy(x)

	h = 0
	error = 1e6
	while error >= 1e-5 && h < 10
		h = h + 1

		# step 1
		dy = y - f(x)

		# step 2
		J = ForwardDiff.jacobian(f, x)
		
		# step 3
		dx = J \dy

		# step 4
		x = x + dx

		# error calculation
		error = maximum(abs.((x - x_) ./ x_))
		x_ = copy(x)
		@printf "Erro na iteração %2i : %.3f. \n" h error 
	end
	println("\nNumero de iterações: $h")

	println("Valor final de x: $x")	
end

# ╔═╡ dcf48af9-5589-4cc4-9f0e-0d7a826a65fa
md"""
## Gauss-Seidel Algorithm applyed to Power-Flow Study


$V_k(i+1)=\frac{1}{Y_{k k}}\left[\frac{\mathrm{P}_k-j \mathrm{Q}_k}{V_k^*(i)}-\sum_{n=1}^{k-1} Y_{k n} V_n(i+1)-\sum_{n=k+1}^N Y_{k n} V_n(i)\right]$


For PV Bus:

$\mathrm{Q}_k=\mathrm{V}_k(i) \sum_{n=1}^N \mathrm{Y}_{k n} \mathrm{~V}_n(i) \sin \left[\delta_k(i)-\delta_n(i)-\theta_{k n}\right]$

"""

# ╔═╡ f4124102-d67c-45a8-a4bf-ae992f3da81a


# ╔═╡ 3ff63ec4-84f8-4466-b3e7-9f4a45b05178
md"""
## Newton-Raphson Applyed to Power-Flow Study


$\begin{aligned}
& y_k=\mathrm{P}_k=\mathrm{P}_k(\mathrm{x})=\mathrm{V}_k \sum_{n=1}^N \mathrm{Y}_{k n} \mathrm{~V}_n \cos \left(\delta_k-\delta_n-\theta_{k n}\right) \\
& y_{k+N}=\mathrm{Q}_{\mathrm{k}}=\mathrm{Q}_k(\mathrm{x})=\mathrm{V}_{\mathrm{k}} \sum_{n=1}^N \mathrm{Y}_{\mathrm{kn}} \mathrm{V}_n \sin \left(\delta_k-\delta_n-\theta_{k n}\right) \\
& k=2,3, \ldots, N
\end{aligned}$

The Jacobian matrix has the form:

$\mathbf{J}
=
\left[\begin{array}{c|c}
\mathbf{J1} & \mathbf{J2} \\
\hline
\mathbf{J3} & \mathbf{J4}
\end{array}\right]
=
\left[\begin{array}{ccc|ccc}
\frac{\partial \mathbf{P}_{2}}{\partial \delta_{2}} & \cdots & \frac{\partial \mathrm{P}_{2}}{\partial \delta_{N}} & \frac{\partial \mathbf{P}_{2}}{\partial \mathbf{V}_{2}} & \cdots & \frac{\partial \mathbf{P}_{2}}{\partial \mathbf{V}_{N}} \\
\vdots & & \vdots & \vdots & & \vdots \\
\frac{\partial \mathbf{P}_{N}}{\partial \delta_{2}} & \cdots & \frac{\partial \mathbf{P}_{N}}{\partial \delta_{N}} & \frac{\partial \mathbf{P}_{N}}{\partial \mathbf{V}_{2}} & \cdots & \frac{\partial \mathbf{P}_{N}}{\partial \mathbf{V}_{N}} \\[0.1in]
\hline
\frac{\partial \mathbf{Q}_{2}}{\partial \delta_{2}} & \cdots & \frac{\partial \mathbf{Q}_{2}}{\partial \delta_{N}} & \frac{\partial \mathbf{Q}_{2}}{\partial \mathbf{V}_{2}} & \cdots & \frac{\partial \mathbf{Q}_{2}}{\partial \mathbf{V}_{N}} \\
\vdots & & \vdots &  \vdots & & \vdots \\
\frac{\partial \mathbf{Q}_{N}}{\partial \delta_{2}} & \cdots & \frac{\partial \mathbf{Q}_{N}}{\partial \delta_{N}} & \frac{\partial \mathbf{Q}_{N}}{\partial \mathbf{V}_{2}} & \ldots & \frac{\partial \mathbf{Q}_{N}}{\partial \mathbf{V}_{N}}
\end{array}\right]$

### Newton-Raphson Algorithm

Starting with $\mathbf{x}(i)=\left[\begin{array}{c}\boldsymbol{\delta}(i) \\ \mathbf{V}(i)\end{array}\right]$ at the $i$ th iteration.

STEP 1: Compute:

$\Delta \mathbf{y}(i)=\left[\begin{array}{l}
\Delta \mathbf{P}(i) \\
\Delta \mathbf{Q}(i)
\end{array}\right]=\left[\begin{array}{l}
\mathbf{P}-\mathbf{P}[\mathbf{x}(i)] \\
\mathbf{Q}-\mathbf{Q}[\mathbf{x}(i)]
\end{array}\right]$

STEP 2: Calculate the Jacobian matrix.

STEP 3: Use Gauss elimination and back substitution to solve:

$\left[\begin{array}{l|l}
\mathbf{J} 1(i) & \mathbf{J} 2(i) \\
\hline \mathbf{J} 3(i) & \mathbf{J} 4(i)
\end{array}\right]\left[\begin{array}{c}
\Delta \boldsymbol{\delta}(i) \\
\Delta \mathbf{V}(i)
\end{array}\right]=\left[\begin{array}{l}
\Delta \mathbf{P}(i) \\
\Delta \mathbf{Q}(i)
\end{array}\right]$

STEP 4: Compute:
        
$\mathbf{x}(i+1)=\left[\begin{array}{l}
\boldsymbol{\delta}(i+1) \\
\mathbf{V}(i+1)
\end{array}\right]=\left[\begin{array}{l}
\boldsymbol{\delta}(i) \\
\mathbf{V}(i)
\end{array}\right]+\left[\begin{array}{c}
\Delta \boldsymbol{\delta}(i) \\
\Delta \mathbf{V}(i)
\end{array}\right]$

Starting with initial value $\mathbf{x}(0)$, the procedure continues until convergence is obtained or until the number of iterations exceeds a specified maximum.

Convergence criteria are often based on $\Delta \mathbf{y}(i)$ (called power mismatches) rather than on $\Delta \mathbf{x}(i)$ (phase angle and voltage magnitude mismatches).

"""

# ╔═╡ bd59a8ca-7ae2-4675-895d-3bb8f57de397
md"""
## Implementação de rede no sistema de arquivos MatPower 
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Calculus = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
Calculus = "~0.5.1"
ForwardDiff = "~0.10.36"
PlutoUI = "~0.7.60"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.5"
manifest_format = "2.0"
project_hash = "6fd6c462a0cef4f779019590ae0a780ab6701439"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

    [deps.ForwardDiff.weakdeps]
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eba4810d5e6a01f612b948c9fa94f905b49087b0"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.60"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╠═8c1fbb78-ec91-4027-ab81-740a097780fc
# ╠═14954992-6e4b-11ef-1c5f-9184c4d14be8
# ╠═b2117c9f-2d47-4edf-a96a-132c5c323275
# ╟─ea96d72c-763e-4a00-acbc-e8686a87efcf
# ╠═db8b3da7-1361-4f95-a0c4-c5a0f3729ff3
# ╠═4c55e5ef-d0bf-4ff5-af48-b8dc9fd3f5cd
# ╟─677a125a-cf5f-49ce-8b09-94171a1e309e
# ╠═3dec9587-a2aa-430d-a2bc-42bc26811eb2
# ╟─bc894bd6-8de5-4392-9400-5bfafdc92200
# ╠═6488b9a6-a040-4657-9aef-791553f6245e
# ╟─8c3c6238-ceb6-4dc9-8e0a-87bceef5f785
# ╟─b1de106c-0378-4ed4-865f-fde4a5cc531b
# ╠═1012f4f9-e5c3-42f0-b2cb-a3c8d40b42ff
# ╟─a5eb2cde-3acd-436c-8adf-5255b3efdd2e
# ╠═dd98c4ac-6741-4712-a30f-ec9ae4be7b97
# ╟─dcf48af9-5589-4cc4-9f0e-0d7a826a65fa
# ╠═f4124102-d67c-45a8-a4bf-ae992f3da81a
# ╟─3ff63ec4-84f8-4466-b3e7-9f4a45b05178
# ╠═bd59a8ca-7ae2-4675-895d-3bb8f57de397
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
