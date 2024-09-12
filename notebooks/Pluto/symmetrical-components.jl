### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 230b7342-9a11-410d-80ea-928bb13e0161
begin
	using LinearAlgebra
	using Printf
	using PlutoUI
end

# ╔═╡ bf724ff4-c45b-4ae3-8e08-6fb670c924fd
PlutoUI.TableOfContents(include_definitions=true)

# ╔═╡ 414f4e30-5a5a-11ef-0515-3535f8040e0e
p(m, a) = m * cis(deg2rad(a))

# ╔═╡ 8fea84bc-560f-4d0d-baf5-8098b1a174d5
function dv(v)
	for i in v
		m = abs(i)
		a = rad2deg(angle(i))
		@printf "%.2f ∠%.2fº\n" m a
	end
end

# ╔═╡ 54289664-b5a4-4bfe-9c4a-a1787e0ae7a6
md"""
## Example 2.4: Balanced $\Delta$ and Y loads
"""

# ╔═╡ 8bcc270a-5202-4e72-b679-a85b612eb887
Eab = p(480.0, 0.0)

# ╔═╡ 0c78ef74-8226-4152-aaae-7c90adb9f52d
Zd = p(30.0, 40.0)

# ╔═╡ bf2e1db2-1aab-4575-92c5-59dcc82d84d4
Zl = p(1.0, 85.0)

# ╔═╡ a0f5b366-aeb4-40dd-b4e6-43ffd59409d4
Zy = Zd / 3.0

# ╔═╡ 1e5db9b3-889a-42f8-9e33-a232d4aced8d
Ean = Eab * p(1.0/√3, -30.0)

# ╔═╡ e39920fa-11a7-406d-aefc-bdeb8062228e
Ia = Ean / (Zl + Zy)

# ╔═╡ 1cb9038e-4112-42e0-853c-84647505be06
abs(Ia), angle(Ia) * 180.0 / π

# ╔═╡ ce6088c0-1a21-4aa2-887a-a6a7b99d0ee3
md"""
As correntes das fases B e C são facilmente obtidas ao atrasar e adiantar em 120º o ângulo da corrente $I_a$
"""

# ╔═╡ 684aca68-c092-48c3-a6c2-e14142dddeb4
md"""
As correntes no delta são dadas por:

$I_{ab} = \frac{|I_{a}|}{\sqrt{3}} \angle (\varphi_{a} + 30^\circ)$
"""

# ╔═╡ 7d850b1c-c8c8-4622-996d-183de9bd4aef
Iab = Ia * p(1.0/√3, 30.0)

# ╔═╡ 681b0883-b4d7-4e97-937a-b31431b40382
abs(Iab), angle(Iab) * 180.0/π

# ╔═╡ 95aa15ed-2ec7-405d-a385-a63543003ad9
md"""
As tensões nos terminais da carga serão:

$V_{ab} = Z_{\Delta} \cdot I_{ab}$

"""

# ╔═╡ 2b11f000-afac-4155-9b4c-999cf816ab55
begin
	Eab_ = Zd * Iab
	dv(Eab_)
end

# ╔═╡ 799f8560-93ce-4926-9dc1-d22c8b62a168
md"""
## Example 8.5: Currents in sequence networks 
"""

# ╔═╡ 25d46568-5226-49db-ac20-742fb85f6258


# ╔═╡ 48d49310-0fca-4c19-9920-128ce58e83bf
md"""
Considera-se que o gerador do exemplo anterior, possui ligação trifásica em Y com neutro aterrado por impedância.

O valor da impedância de aterramento do neutro e de sequência zero, positiva e negativa são:
"""

# ╔═╡ 5a6a2c25-8cd4-4a4d-b093-94b5faba6c23
Zgn = 10.0im

# ╔═╡ d3535419-7d1c-460d-a5c3-d7511df9e100
Zg0 = 1.0im

# ╔═╡ b76cecd9-3728-4e20-b481-eed3edbbea6f
Zg1 = 15.0im

# ╔═╡ 669282f1-febc-46b2-a482-c2169e03cec5
Zg2 = 3.0im

# ╔═╡ 2722e4cf-e71b-4eee-baca-83cffcccabe5


# ╔═╡ af41f30e-3d96-4795-b5df-902bb77187ae
md"""
Para a condição de fonte e carga balanceadas e linha completamente transposta, só haverá corrente no circuito de sequência positiva. 
"""

# ╔═╡ b8f65c7a-30cb-4827-911a-85fefc2477c5
md"""
Definição das tensões de sequência nos terminais do gerador:
"""

# ╔═╡ 97c10bda-2f89-49aa-9c0b-aaf2c53c6171
V1 = Ean;

# ╔═╡ e9e3b772-3209-4a6d-a1fd-9e9e20cac785
V0 = V2 = 0.0;

# ╔═╡ 62a0bd46-7a98-4bc5-8530-88143a5c0f84


# ╔═╡ aabec488-1a08-46bc-bd1d-a270daaa6906
md"""
Definição das impedâncias de sequência na linha:
"""

# ╔═╡ 0c494ca5-9c32-4e87-8931-c0e73318717d
Zl0 = Zl1 = Zl2 = Zl

# ╔═╡ 214eba83-b629-4793-9064-fe9c4d872bf6


# ╔═╡ 35e1da14-a9c7-4ed4-9811-c479ff3c8412
md"""
Definição das impedâncias de sequência da carga:
"""

# ╔═╡ 788d8617-9de2-4b8c-ba62-9ecb28629e32
Zlo1 = Zlo2 = Zy

# ╔═╡ 14cfd883-1729-4a70-91bf-760f64b3bec8
Zlo0 = Zy + Inf

# ╔═╡ 95720c9c-f134-4cfe-9e47-3f93776f0574


# ╔═╡ f80d911a-edea-42cb-a7fb-231be5396c26
md"""

Como o circuito é totalmente balanceado, não haverá correntes de sequência zero e negativa, somente de sequência positiva.

Cálculo da corrente de sequência positiva: 
"""

# ╔═╡ 73d70622-5efa-4086-ad0a-7f952519ad76
begin
	I1 = V1 / (Zl1 + Zlo1)
	dv(I1)
end

# ╔═╡ b727fa22-8237-4234-9fd5-4139a921eb63


# ╔═╡ 24d24ed5-26b6-41f0-8705-ce5b9210eed4
md"""
## Example 8.6: Solving unbalanced three-phase networks using sequence components
"""

# ╔═╡ a4c50eab-282c-4aa9-b81e-531ad5563be3
md"""
Agora considere que o seguinte vetor de tensões represente tensões desbalanceadas nos terminais do gerador:
"""

# ╔═╡ 860b78c4-5a3f-4fe0-b51e-b86c3b057d1a
md"""
$\left[
\begin{array}{c}
V_{ag} \\
V_{bg} \\
V_{cg}
\end{array}
\right]
=
\left[
\begin{array}{c}
277,0 \angle{0^\circ} \\
260,0 \angle{-120^\circ}\\
295,0 \angle{+115^\circ}
\end{array}
\right]
~Volts$
"""

# ╔═╡ 4b19aeea-0a77-49e9-b656-d6cf34011cc6
Vabc = [p(277.0, 0.0); p(260.0, -120.0); p(295.0, 115.0);;]

# ╔═╡ b021632e-02c7-4e8a-a5a3-5564d821ff44
a = p(1.0, 120.0)

# ╔═╡ 912c093f-9537-4a48-a919-3ab9924e4ded
A = [1.0 1.0 1.0;
	 1.0 a^2 a;
	 1.0 a   a^2]

# ╔═╡ d6ef13ae-d766-4514-b7ed-e87029eb8815
md"""
O vetor de tensões de sequência nos terminais do gerador será dado por:
"""

# ╔═╡ 89c5861c-f2ee-48fc-b7b9-ed855d957785
V012 = inv(A) * Vabc;

# ╔═╡ 497a4515-ea91-4dee-ab63-63e027e031e7
md"""
Tensões de sequência:
"""

# ╔═╡ beff506f-47da-4515-8cef-009c8c60eaa5
dv(V012)

# ╔═╡ 35a01d8c-31d7-4cb6-abce-a3f692da0eec
md"""
Cálculo da corrente de sequência zero:
"""

# ╔═╡ feaac5e0-1735-4188-9c59-ea96cdc1a956
I0 = V012[1] / (Zl0 + Zlo0)

# ╔═╡ 11fd5eaf-427d-48b5-9dcc-01b3bc1237de
md"""
Cálculo da corrente de sequência positiva:
"""

# ╔═╡ 7418c0bd-9483-485a-a93c-458a5176bdae
I1_ = V012[2] / (Zl1 + Zlo1);

# ╔═╡ a3226aa6-cf7f-40e4-bf9f-a510233dfef1
dv(I1_)

# ╔═╡ bbba001a-74b7-4a58-9631-7af43ead80c1
md"""
Cálculo da corrente de sequência negativa:
"""

# ╔═╡ 66df8921-cb84-4ba0-b2b9-5600b3919efe
I2_ = V012[3] / (Zl2 + Zlo2);

# ╔═╡ d951e7ef-5ff6-4bfd-a458-edd091d0f6c1
dv(I2_)

# ╔═╡ 7986d5ef-dae8-47ce-8a72-fe6d77407a24


# ╔═╡ c2f5d6d0-4df6-4382-9f00-b7adb6ae6f6b
md"""
Cálculo do vetor de correntes de fase, a partir do vetor de correntes de sequência:
"""

# ╔═╡ a3a4c781-25b4-428a-8062-809753affdb7
Iabc = A * [I0; I1_; I2_;;];

# ╔═╡ 5082bca4-4ff1-4857-bb08-84496cdba3e0
md"""
Correntes de carga:
"""

# ╔═╡ 98fa8cb4-b15d-4a75-8ecb-4c606af182e6
dv(Iabc)

# ╔═╡ 92f99aa2-2fbf-4f05-9a70-ef0ff6527dd4


# ╔═╡ af20cb9f-b8b3-4564-89fd-89d88c6ffb9a
md"""
## Example 8.7: Solving unbalanced three-phase networks with transformers using per-unit sequence components
"""

# ╔═╡ b4c2d45f-0140-4076-a4b5-0062d1641e19
md"""
A 75-kVA, 480-volt D/208-volt Y transformer with a solidly grounded neutral
is connected between the source and line of Example 8.6. 

The transformer leakage reactance is Xeq = 0.10 per unit;

winding resistances and exciting current are neglected.

Using the transformer ratings as base quantities, draw the per-unit sequence networks and calculate the phase a source current Ia.
"""

# ╔═╡ 483f6dbf-d1fa-4b4e-9450-ada2b80f063e
md"""
Impedância do transformador em pu:
"""

# ╔═╡ 0eef1ba0-1149-48a9-bf9e-c954f78b4c4e
Xtpu = 0.1im

# ╔═╡ 0e963201-d654-4ad8-aa7b-a60de9f58147
md"""
Potência de base monofásica:
"""

# ╔═╡ 85c9d4ba-f124-4653-9c1a-ac38bbab9b89
Sbase1 = 75e3/3

# ╔═╡ 071e445b-273a-4178-9ead-25921c4142ea
md"""
Tensão de base fase-neutro no primário:
"""

# ╔═╡ 44d00482-41c6-4f7f-ac7b-1abe42b3675e
Vbase1 = 480.0 / √3

# ╔═╡ 79ff612a-4fee-4c16-9688-7568335d3256
md"""
Tensão de base fase-neutro no secundário:
"""

# ╔═╡ 35352a29-9199-40cb-bf2b-647ac5f8d7bf
Vbase2 = 208.0 / √3

# ╔═╡ 41c1321f-8076-49e1-93e6-339884bb8621
Zbase2 = Vbase2^2 / Sbase1

# ╔═╡ 3936d41e-97e0-43f3-a4ad-5031fbe6ad8e
begin
	V012pu = V012 / Vbase1;
	dv(V012pu)
end

# ╔═╡ 8220f9c1-e769-4c30-8698-986e757544fc
begin
	Zl0pu = Zl1pu = Zl2pu = Zl / Zbase2;
	abs(Zl1pu), angle(Zl1pu) * 180.0 / π
end

# ╔═╡ 7e7e0601-e0b3-482a-934d-e9bf1c38b073
begin
	Zlo1pu = Zlo2pu = Zy / Zbase2;
	abs(Zlo1pu), angle(Zlo1pu) * 180.0 / π
end

# ╔═╡ f8ae6105-def2-4e6a-b666-44b5c754462c
I0t = 0.0

# ╔═╡ f1d80d3e-af82-48eb-9411-5445d0afba30
begin
	I1t = V012pu[2] / (Xtpu + Zl1pu + Zlo1pu);
	abs(I1t), angle(I1t) * 180.0 / π
end

# ╔═╡ 8c9d0e19-bcfe-4b19-a0ea-6e6ddf7ddfca
begin
	I2t = V012pu[3] / (Xtpu + Zl2pu + Zlo2pu);
	abs(I2t), angle(I2t) * 180.0 / π
end

# ╔═╡ f5d50cf9-9697-4615-848a-e2b4e103e177
begin
	Iabct = A * [I0t; I1t; I2t];
	dv(Iabct)
end

# ╔═╡ 5b690a46-d5ca-4ebc-869a-9820442d7a60
md"""
Para o primário teremos a corrente de base dada por: 
"""

# ╔═╡ 77ad25d8-a4b1-4809-b948-11812895fac6
Ibaset = 75e3 / (√3 * 480.0)

# ╔═╡ ba0077e3-992f-4b2b-8b19-a691227be0aa
md"""
Correntes de fase no primário do transformador, para a condição de carga analisada: 
"""

# ╔═╡ 6926741a-bccc-4e6c-84a7-b643f77a3c02
dv(Iabct * Ibaset)

# ╔═╡ f6e41482-80e7-4215-b346-50c0abfff829


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
PlutoUI = "~0.7.59"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "be6bb64aaa5b6d21578592978fee0b673444b366"

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

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

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

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

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

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

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
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

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
version = "5.8.0+1"

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
# ╠═230b7342-9a11-410d-80ea-928bb13e0161
# ╠═bf724ff4-c45b-4ae3-8e08-6fb670c924fd
# ╠═414f4e30-5a5a-11ef-0515-3535f8040e0e
# ╠═8fea84bc-560f-4d0d-baf5-8098b1a174d5
# ╟─54289664-b5a4-4bfe-9c4a-a1787e0ae7a6
# ╠═8bcc270a-5202-4e72-b679-a85b612eb887
# ╠═0c78ef74-8226-4152-aaae-7c90adb9f52d
# ╠═bf2e1db2-1aab-4575-92c5-59dcc82d84d4
# ╠═a0f5b366-aeb4-40dd-b4e6-43ffd59409d4
# ╠═1e5db9b3-889a-42f8-9e33-a232d4aced8d
# ╠═e39920fa-11a7-406d-aefc-bdeb8062228e
# ╠═1cb9038e-4112-42e0-853c-84647505be06
# ╟─ce6088c0-1a21-4aa2-887a-a6a7b99d0ee3
# ╟─684aca68-c092-48c3-a6c2-e14142dddeb4
# ╠═7d850b1c-c8c8-4622-996d-183de9bd4aef
# ╠═681b0883-b4d7-4e97-937a-b31431b40382
# ╟─95aa15ed-2ec7-405d-a385-a63543003ad9
# ╠═2b11f000-afac-4155-9b4c-999cf816ab55
# ╟─799f8560-93ce-4926-9dc1-d22c8b62a168
# ╟─25d46568-5226-49db-ac20-742fb85f6258
# ╟─48d49310-0fca-4c19-9920-128ce58e83bf
# ╠═5a6a2c25-8cd4-4a4d-b093-94b5faba6c23
# ╠═d3535419-7d1c-460d-a5c3-d7511df9e100
# ╠═b76cecd9-3728-4e20-b481-eed3edbbea6f
# ╠═669282f1-febc-46b2-a482-c2169e03cec5
# ╟─2722e4cf-e71b-4eee-baca-83cffcccabe5
# ╟─af41f30e-3d96-4795-b5df-902bb77187ae
# ╟─b8f65c7a-30cb-4827-911a-85fefc2477c5
# ╠═97c10bda-2f89-49aa-9c0b-aaf2c53c6171
# ╠═e9e3b772-3209-4a6d-a1fd-9e9e20cac785
# ╟─62a0bd46-7a98-4bc5-8530-88143a5c0f84
# ╟─aabec488-1a08-46bc-bd1d-a270daaa6906
# ╠═0c494ca5-9c32-4e87-8931-c0e73318717d
# ╟─214eba83-b629-4793-9064-fe9c4d872bf6
# ╟─35e1da14-a9c7-4ed4-9811-c479ff3c8412
# ╠═788d8617-9de2-4b8c-ba62-9ecb28629e32
# ╠═14cfd883-1729-4a70-91bf-760f64b3bec8
# ╟─95720c9c-f134-4cfe-9e47-3f93776f0574
# ╟─f80d911a-edea-42cb-a7fb-231be5396c26
# ╠═73d70622-5efa-4086-ad0a-7f952519ad76
# ╟─b727fa22-8237-4234-9fd5-4139a921eb63
# ╟─24d24ed5-26b6-41f0-8705-ce5b9210eed4
# ╟─a4c50eab-282c-4aa9-b81e-531ad5563be3
# ╟─860b78c4-5a3f-4fe0-b51e-b86c3b057d1a
# ╠═4b19aeea-0a77-49e9-b656-d6cf34011cc6
# ╠═b021632e-02c7-4e8a-a5a3-5564d821ff44
# ╠═912c093f-9537-4a48-a919-3ab9924e4ded
# ╟─d6ef13ae-d766-4514-b7ed-e87029eb8815
# ╠═89c5861c-f2ee-48fc-b7b9-ed855d957785
# ╟─497a4515-ea91-4dee-ab63-63e027e031e7
# ╠═beff506f-47da-4515-8cef-009c8c60eaa5
# ╟─35a01d8c-31d7-4cb6-abce-a3f692da0eec
# ╠═feaac5e0-1735-4188-9c59-ea96cdc1a956
# ╟─11fd5eaf-427d-48b5-9dcc-01b3bc1237de
# ╠═7418c0bd-9483-485a-a93c-458a5176bdae
# ╠═a3226aa6-cf7f-40e4-bf9f-a510233dfef1
# ╟─bbba001a-74b7-4a58-9631-7af43ead80c1
# ╠═66df8921-cb84-4ba0-b2b9-5600b3919efe
# ╠═d951e7ef-5ff6-4bfd-a458-edd091d0f6c1
# ╟─7986d5ef-dae8-47ce-8a72-fe6d77407a24
# ╟─c2f5d6d0-4df6-4382-9f00-b7adb6ae6f6b
# ╠═a3a4c781-25b4-428a-8062-809753affdb7
# ╟─5082bca4-4ff1-4857-bb08-84496cdba3e0
# ╠═98fa8cb4-b15d-4a75-8ecb-4c606af182e6
# ╟─92f99aa2-2fbf-4f05-9a70-ef0ff6527dd4
# ╟─af20cb9f-b8b3-4564-89fd-89d88c6ffb9a
# ╟─b4c2d45f-0140-4076-a4b5-0062d1641e19
# ╟─483f6dbf-d1fa-4b4e-9450-ada2b80f063e
# ╠═0eef1ba0-1149-48a9-bf9e-c954f78b4c4e
# ╟─0e963201-d654-4ad8-aa7b-a60de9f58147
# ╠═85c9d4ba-f124-4653-9c1a-ac38bbab9b89
# ╟─071e445b-273a-4178-9ead-25921c4142ea
# ╠═44d00482-41c6-4f7f-ac7b-1abe42b3675e
# ╟─79ff612a-4fee-4c16-9688-7568335d3256
# ╠═35352a29-9199-40cb-bf2b-647ac5f8d7bf
# ╠═41c1321f-8076-49e1-93e6-339884bb8621
# ╠═3936d41e-97e0-43f3-a4ad-5031fbe6ad8e
# ╠═8220f9c1-e769-4c30-8698-986e757544fc
# ╠═7e7e0601-e0b3-482a-934d-e9bf1c38b073
# ╠═f8ae6105-def2-4e6a-b666-44b5c754462c
# ╠═f1d80d3e-af82-48eb-9411-5445d0afba30
# ╠═8c9d0e19-bcfe-4b19-a0ea-6e6ddf7ddfca
# ╠═f5d50cf9-9697-4615-848a-e2b4e103e177
# ╟─5b690a46-d5ca-4ebc-869a-9820442d7a60
# ╠═77ad25d8-a4b1-4809-b948-11812895fac6
# ╟─ba0077e3-992f-4b2b-8b19-a691227be0aa
# ╠═6926741a-bccc-4e6c-84a7-b643f77a3c02
# ╠═f6e41482-80e7-4215-b346-50c0abfff829
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
