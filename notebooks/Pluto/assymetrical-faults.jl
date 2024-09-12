### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 84181e0c-5e4e-11ef-12c9-b9603ac2cbdb
begin
	using PlutoUI
	using LinearAlgebra
end

# ╔═╡ e08927f0-cf6e-4d76-9794-df9c224f6972
md"""
# Chapter 9: Unsymmetrical Faults
"""

# ╔═╡ d88c87b7-b654-4ecf-ab5e-9bef6bbec900
PlutoUI.TableOfContents()

# ╔═╡ 300b0cd1-b493-488a-806e-9b449a40e625
md"""
## Equacionamento utilizado para cálculo de faltas assimétricas
"""

# ╔═╡ e705c595-b289-4b60-891f-fb1183e615fd
md"""
Definição da Matriz de transformação A e de sua inversa:

Conforme definido da teoria das componentes simétricas, a matriz de transformação entre o domínio de fase e o domínio das componentes simétricas pode ser dada pela seguinte operação linear:

$\left[\begin{array}{l}
V_a \\
V_b \\
V_c
\end{array}\right]=\left[\begin{array}{ccc}
\mathrm{1} & \mathrm{1} & \mathrm{1} \\
\mathrm{1} & a^2 & a \\
\mathrm{1} & a & a^2
\end{array}\right]\left[\begin{array}{c}
V_a^{(\mathrm{o})} \\
V_a^{(\mathrm{I})} \\
V_a^{(2)}
\end{array}\right]$

Em que:

$A=\left[\begin{array}{ccc}
\mathrm{1} & \mathrm{1} & \mathrm{1} \\
\mathrm{1} & a^2 & a \\
\mathrm{1} & a & a^2
\end{array}\right]$

"""

# ╔═╡ 76553e35-2f5a-46a7-9db1-c72376c0177f
md"""
Tensões de sequência nos terminais da falta:

$\left[\begin{array}{l}
V_0 \\
V_1 \\
V_2
\end{array}\right]=\left[\begin{array}{c}
0 \\
V_{\mathrm{F}} \\
0
\end{array}\right]-\left[\begin{array}{ccc}
Z_0 & 0 & 0 \\
0 & Z_1 & 0 \\
0 & 0 & Z_2
\end{array}\right]\left[\begin{array}{l}
I_0 \\
I_1 \\
I_2
\end{array}\right]$

"""

# ╔═╡ e4c4eaf0-c5c8-4f7d-8f5c-44d10fe9134b
md"""
### Equações para faltas monofásicas
"""

# ╔═╡ e7cea086-700f-4061-91c8-5363e3d9faed
md"""
Correntes de sequência:

$I_0=I_1=I_2=\frac{V_{\mathrm{F}}}{Z_0+\mathrm{Z}_1+\mathrm{Z}_2+\left(3 Z_{\mathrm{F}}\right)}$

Corrente na fase atingida:

$I_a=I_0+I_1+I_2=3 I_1=\frac{3 V_{\mathrm{F}}}{Z_0+\mathrm{Z}_1+\mathrm{Z}_2+\left(3 Z_{\mathrm{F}}\right)}$

"""

# ╔═╡ 51c3e84d-8df3-4fd8-a1f4-a9513c2209c3
md"""
### Equações para faltas bifásicas

Corrente de sequência:

$I_1=-I_2=\frac{V_{\mathrm{F}}}{\left(Z_1+Z_2+Z_{\mathrm{F}}\right)}$

$I_0=0$

Corrente na fase atingida:

$\begin{aligned}
I_b & =I_0+a^2 I_1+a I_2=\left(a^2-a\right) I_1 \\
& =-j \sqrt{3} I_1=\frac{-j \sqrt{3} V_{\mathrm{F}}}{\left(Z_1+Z_2+Z_{\mathrm{F}}\right)}
\end{aligned}$

$I_a = I_0 + I_1 + I_2 = 0$

$I_c = - I_b$
"""

# ╔═╡ 5f20c381-3a1e-4709-9fcd-dfdcab9c3363
md"""
### Equações para faltas bifásicas à terra

Correntes de sequência:

$\begin{aligned}
& I_1=\frac{V_{\mathrm{F}}}{Z_1+\left[Z_2 / /\left(Z_0+3 Z_{\mathrm{F}}\right)\right]}=\frac{V_{\mathrm{F}}}{Z_1+\left[\frac{Z_2\left(Z_0+3 Z_{\mathrm{F}}\right)}{Z_2+Z_0+3 Z_{\mathrm{F}}}\right]} \\
& I_2=\left(-I_1\right)\left(\frac{Z_0+3 Z_{\mathrm{F}}}{Z_0+3 Z_{\mathrm{F}}+Z_2}\right) \\
& I_0=\left(-I_1\right)\left(\frac{Z_2}{Z_0+3 Z_{\mathrm{F}}+Z_2}\right)
\end{aligned}$

Correntes de fase:

$\left[\begin{array}{l}
I_a \\
I_b \\
I_c
\end{array}\right]=\left[\begin{array}{ccc}
\mathrm{1} & \mathrm{1} & \mathrm{1} \\
\mathrm{1} & a^2 & a \\
\mathrm{1} & a & a^2
\end{array}\right]\left[\begin{array}{c}
I_a^{(\mathrm{o})} \\
I_a^{(\mathrm{I})} \\
I_a^{(2)}
\end{array}\right]$

"""

# ╔═╡ 12971eca-f21e-4c3a-917e-ec0af879da60
p(m, a) = abs(m) * cosd(a) + 1im * abs(m) * sind(a)

# ╔═╡ 40a32e16-8787-4b98-909c-45c93b7672d7
function dv(v)
	r = zeros(3, 2)
	for (i, j) in zip(v, 1:3)
		r[j, 1] = abs(i)
		r[j, 2] = rad2deg(angle(i))
		# append!(r, (abs(i), rad2deg(angle(i))))
	end
	return r
end

# ╔═╡ 360e119e-87ae-43f4-8f54-4ad21166002e
a = p(1.0, 120.0)

# ╔═╡ 984284fd-67c8-4f3f-8882-f52914321867
A = [1.0 1.0 1.0;
	 1.0 a^2 a;
	 1.0 a a^2]

# ╔═╡ 893d06a2-e16f-4e89-b3aa-9605f95956d9
invA = inv(A)

# ╔═╡ feeea186-3682-43f1-80eb-22d640f458a8
md"""
## Example 9.1: Power-system sequence networks and their Thévenin equivalents

A single-line diagram of the power system considered in Example 7.3 is shown in Figure 9.3, where negative- and zero-sequence reactances are also given. 

The neutrals of the generator and D–Y transformers are solidly grounded. 

The motor neutral is grounded through a reactance Xn = 0.05 per unit on the motor base.

!!! note "Items"
	- (a) Draw the per-unit zero-, positive-, and negative-sequence networks on a 100-MVA, 13.8-kV base in the zone of the generator.
	- (b) Reduce the sequence networks to their Thévenin equivalents, as viewed from bus 2. 

Prefault voltage is VF = 1.05 /_ 0º per unit.

!!! important "Consideration"
	The prefault load current and D–Y transformer phase shift are neglected.

![Single Line Diagram of this example](https://lh3.googleusercontent.com/d/1nn3V7zmCJE9K9hft3PrJoNnNmsSn3xgO)

"""

# ╔═╡ edbe576f-a0cb-48f2-b48c-b39a6b4ae902
md"""
### Resolução
"""

# ╔═╡ 380fdaeb-2aa1-427b-a557-3fcbd1ae48e4
md"""
Assim como realizado no exemplo 7.4 a única impedância que precisa ser convertida para pu é a impedância da linha de transmissão.

Considerando então uma potência de base de 100MVA e uma tensão de base idêntica à tensão nominal do sistema, teremos: 
"""

# ╔═╡ d9782a87-db49-4110-b1d6-74bd71a78628
Zbasel = 138.0^2 / 100.0

# ╔═╡ 7735ff4d-a08e-468c-b781-89d5135b150b
Xline1 = Xline2 = 20.0

# ╔═╡ 689093f4-1ebf-4915-b4b4-f5721efd6c6e
Xlinepu1 = Xlinepu2 = Xline1 / Zbasel

# ╔═╡ 41a555ee-9b09-430a-8c00-de073805b55d
md"""
Também a impedância de sequência zero precisa ser transformada para p.u:
"""

# ╔═╡ 6f2d7d97-743a-4052-a26a-b45a9be08f78
Xline0 = 60.0

# ╔═╡ e9bfb3db-4544-4b50-a8d4-b84d3dc7f0ee
Xlinepu0 = Xline0 / Zbasel

# ╔═╡ 1546e93c-4f2d-4e4b-9a2d-71d3dfcefb24
md"""
Podemos agora definir todas as impedâncias de sequência do problema
"""

# ╔═╡ b8c4bdc8-1605-48f3-8ec9-ddaa998eaf24
md"""
Impedâncias de sequência do gerador:
"""

# ╔═╡ 193994a5-509f-442a-af92-8a2e4de54ddd
begin
	Xg1 = 0.15
	Xg2 = 0.17
	Xg0 = 0.05
	Xgn = 0.0
end;

# ╔═╡ 2f70870c-c553-4ddc-a481-173647d0d71e
md"""
Impedâncias de sequência do transformador 1 e do transformador 2:
"""

# ╔═╡ efddf058-44f2-4ca4-8ef6-6b0192aa88e8
begin
	Xt0 = Xt1 = Xt2 = 0.1
end;

# ╔═╡ 380f0924-91a1-4703-9e70-eb62f7ac8c65
md"""
Impedâncias de sequência da linha de transmissão:
"""

# ╔═╡ bb9541aa-7a12-4b19-ae80-06e532af8d08
begin
	Xl1 = Xl2 = Xlinepu1
	Xl0 = Xlinepu0
end;

# ╔═╡ 88287684-16d5-4c8e-a553-90a9fe24a7e9
md"""
Impedância de sequência do motor síncrono:
"""

# ╔═╡ 58e4de6d-9b23-41a9-a39d-132cd217443c
begin
	Xm1 = 0.2
	Xm2 = 0.21
	Xm0 = 0.1
	Xmn = 0.05
end;

# ╔═╡ 618ccc1e-8006-403b-981c-5b39d42f7fb3
md"""
Em relação à barra 2 do sistema, podemos calcular a impeância de Thevennin equivalente para cada um dos circuitos de sequência:
"""

# ╔═╡ 09ca922b-d8e1-44c9-ac56-143ad66c6a1d
Zth2_0 = 1im * (Xm0 + 3 * Xmn)

# ╔═╡ ea84f5cd-7ae6-4c39-9084-62854c655dff
begin
	_1 = Xg1 + 2 * Xt1 + Xl1
	_2 = Xm1
	Zth2_1 = 1im * (_1 * _2 / (_1 + _2))
end

# ╔═╡ f0da4182-4506-4423-910d-5fdb4e070a3c
begin
	_3 = Xg2 + 2 * Xt2 + Xl2
	_4 = Xm2
	Zth2_2 = 1im * (_3 * _4 / (_3 + _4))
end

# ╔═╡ bbe49832-b653-4e5d-becc-ddce8cc32b35


# ╔═╡ 2200f673-3d81-4d72-b987-23087ded3851
md"""
## Example 9.2: Three-phase short-circuit calculations using sequence networks

!!! important "Question"
	Calculate the per-unit subtransient fault currents in phases a, b, and c for a bolted three-phase-to-ground short circuit at bus 2 in Example 9.1.
"""

# ╔═╡ 925adf5e-bdf8-4377-bd75-2aca5c8c3ad6
md"""
### Resolução
"""

# ╔═╡ 84b51bc4-557f-43ef-8c8c-2c0e83997e49
md"""
Como para o curto-circuito trifásico só há a presença de correntes de sequência positiva, então a determinação do cálculo de suas correntes resume-se à razão entre a tensão de pré-falta e a impedância de sequência:
"""

# ╔═╡ 9ea8e62f-f8dc-4e6d-8330-661f48901ac8
Icc_3ph_0 = 0.0

# ╔═╡ ccda01a4-c3e5-45e9-89d5-09a76009c4db
Icc_3ph_1 = 1.05 / Zth2_1

# ╔═╡ a08d7997-57d0-4c9b-a587-1c0059dae14a
Icc_3ph_2 = 0.0

# ╔═╡ 7f0861bc-13b3-49fc-beea-5773635aa2e2
Icc_3ph_012 = [Icc_3ph_0; Icc_3ph_1; Icc_3ph_2;;]

# ╔═╡ 5d07a2a4-9199-4165-a338-f0bfd419ef4f
md"""
Calculando a corrente de curto-circuito para cada uma das fases:
"""

# ╔═╡ c53172dd-ed14-4014-8a26-628e266909d9
Iabc_cc_3ph = A * Icc_3ph_012

# ╔═╡ 406fa646-dbaf-43e9-b085-7f469e213865
dv(Iabc_cc_3ph) # exibe os resultados em módulo e ângulo

# ╔═╡ 17e5fa14-9932-4674-98e9-410372e406b1
md"""
As tensões de sequência nos terminais da falta podem ser obtidas fazendo: 
"""

# ╔═╡ ce0556fd-400b-4edf-95f6-e4134a07b69f
V012_3ph = [0.0; 1.05; 0.0] - diagm([Zth2_0; Zth2_1; Zth2_2]) * Icc_3ph_012

# ╔═╡ c98c4c36-b80e-434e-bbb7-3345e418c354
Vabc_3ph = A * V012_3ph

# ╔═╡ 54ac7a03-15fd-4aba-9acc-8ae6949e6639


# ╔═╡ 6fb52297-c521-44fc-9568-10b6fe883cb3
md"""
## Example 9.3: Single line-to-ground short-circuit calculations using sequence networks

!!! important "Question"
	Calculate the subtransient fault current in per-unit and in kA for a bolted single line-to-ground short circuit from phase a to ground at bus 2 in Example 9.1. 

	Also calculate the per-unit line-to-ground voltages at faulted bus 2.
"""

# ╔═╡ c31ca637-7e19-41f6-b694-a0a6fcec193b
md"""
### Resolução
"""

# ╔═╡ 44c97ab1-9f6c-4f52-93bb-c73451cce2ed
md"""
Utilizando o equacionamento que resulta da aplicação de componentes simétricas para o caso de faltas monofásicas, teremos:
"""

# ╔═╡ b27c134f-6736-465f-b473-5060f2099cf9
Zf_1ph = 0.0 # impedância de falta

# ╔═╡ c72b4480-7912-4fb6-aef6-399c07489e51
Icc_1ph_1 = Icc_1ph_2 = Icc_1ph_0 = 1.05 / (Zth2_0 + Zth2_1 + Zth2_2 + Zf_1ph)

# ╔═╡ c4de09dd-abac-4606-b702-8b795c12075d
Icc_1ph_012 = [Icc_1ph_0; Icc_1ph_1; Icc_1ph_2;;]

# ╔═╡ bf4eeeb3-7d07-4644-aa3f-3a47cce827c4
md"""
Convertendo para o domínio de fase:
"""

# ╔═╡ e35721b0-2f1e-4205-a66f-2a048d4e7edb
Icc_1ph_abc = A * Icc_1ph_012

# ╔═╡ 44613696-6087-4c8b-8065-a90ea846c81b
dv(Icc_1ph_abc)

# ╔═╡ 0ac51b1c-df6b-43bd-8d30-e535d3b91218
Ibase = (100.0 / 3.0) / (13.8 / √3) * 1e3

# ╔═╡ 168b3cd6-25d8-4ee9-9ab7-824ba3f66777
dv(Icc_1ph_abc * Ibase)

# ╔═╡ 7a893665-21a8-4a10-bec2-81db48e4cf0b
md"""
Cálculo das tensões de sequência nos terminais da falta:
"""

# ╔═╡ c668ce3e-c6eb-4d9c-a779-ecc11ce44b73
V012_1ph = [0.0; 1.05; 0.0] - diagm([Zth2_0; Zth2_1; Zth2_2]) * Icc_1ph_012

# ╔═╡ 4edbadc7-8ec1-4068-b36c-74739e40816f
begin
	Vabc_1ph = A * V012_1ph
	dv(Vabc_1ph)
end

# ╔═╡ eba51237-c21c-476d-8ca7-2848ce012aac


# ╔═╡ 607e74ad-ef8b-4c5c-9e62-da9129593fac
md"""
## Example 9.4:  Line-to-line short-circuit calculations using sequence networks

!!! important "Question"
	Calculate the subtransient fault current in per-unit and in kA for a bolted line-to- line fault from phase b to c at bus 2 in Example 9.1.
"""

# ╔═╡ b660051c-ea82-4c62-ae44-2452de3197da
md"""
### Resolução
"""

# ╔═╡ 9c220de1-3a66-4b94-93af-596c8fa6c2f3
md"""
Utilizando o equacionamento que resulta da aplicação de componentes simétricas para o caso de faltas bifásicas, teremos:
"""

# ╔═╡ efe8408e-c652-4c22-b68d-d90dd3066213
Zf_2ph = 0.0 #impedância para a falta bifásica 

# ╔═╡ b30cb39e-dd05-4b46-b1df-f64ab36fa1fd
md"""
Por meio da montagem do circuito é possível obter:
"""

# ╔═╡ 3ab747ea-ab51-4280-adbb-8932e2ed3b9f
begin
	Icc_2ph_1 = 1.05 / (Zth2_1 + Zth2_2 + Zf_2ph)
	Icc_2ph_2 = - Icc_2ph_1
	Icc_2ph_0 = 0.0
end

# ╔═╡ ebcdd089-c1f8-4a75-8c18-f47901e8ee2b
Icc_2ph_012 = [Icc_2ph_0; Icc_2ph_1; Icc_2ph_2;;]

# ╔═╡ cc8b842d-eace-423d-b6dd-0286e6c08527
md"""
Convertendo para o domínio de fase:
"""

# ╔═╡ 5e40efa2-0ded-4259-b4bd-2da7e9ff0f32
Icc_2ph_abc = A * Icc_2ph_012

# ╔═╡ 7b6a5ab0-b90c-4818-a4f2-a214b4f2c25c
dv(Icc_2ph_abc)

# ╔═╡ bd7f495f-43f4-48cf-83ae-32059462bd77
dv(Icc_2ph_abc * Ibase)

# ╔═╡ 5baeaa55-31ef-4d1e-8073-28dc761d4429
md"""
Cálculo das tensões de sequência nos terminais da falta:
"""

# ╔═╡ 8cb55546-4282-4d7d-a1db-b7ed55451709
V012_2ph = [0.0; 1.05; 0.0] - diagm([Zth2_0; Zth2_1; Zth2_2]) * Icc_2ph_012

# ╔═╡ 7a97c09a-da58-4ec8-8551-67fa1bafcb5d
begin
	Vabc_2ph = A * V012_2ph
	dv(Vabc_2ph)
end

# ╔═╡ 7d1e1286-bea2-4ef8-9937-10c0495127a6


# ╔═╡ cb0173e8-d981-4e0c-ab3d-d71c5cfd18bc
md"""
## Example 9.5: Double line-to-ground short-circuit calculations using sequence networks

!!! important "Question"
	Calculate:
	
	- (a) the subtransient fault current in each phase, 
	- (b) neutral fault current, and
	- (c) contributions to the fault current from the motor and from the transmission line, for a bolted double line-to-ground fault from phase b to c to ground at bus 2 in Example 9.1. Neglect the D–Y transformer phase shifts.
"""

# ╔═╡ 059ad344-23d4-4efd-b81b-d801d9a1fe4f
md"""
### Resolução
"""

# ╔═╡ 9743cc24-f97a-458f-95d0-a0f8d0618fda
md"""
Utilizando o equacionamento que resulta da aplicação de componentes simétricas para o caso de faltas bifásicas à terra, teremos:
"""

# ╔═╡ d225fe5c-6a09-4404-ac9e-241b653fcf6b
Zf_2phg = 0.0

# ╔═╡ c03a2907-9bfb-4c63-b4d2-d0e9d89df250
begin
	Icc_2phg_1 = 1.05 / (Zth2_1 + (Zth2_2 * (Zth2_0 + 3.0 * Zf_2phg)) / (Zth2_2 + Zth2_0 + 3.0 * Zf_2phg));
	Icc_2phg_2 = - Icc_2phg_1 * (Zth2_0 + 3.0 * Zf_2phg) / (Zth2_2 + Zth2_0 + 3.0 * Zf_2phg);
	Icc_2phg_0 = - Icc_2phg_1 * Zth2_2 / (Zth2_2 + Zth2_0 + 3.0 * Zf_2phg);
end

# ╔═╡ 99317f2d-2fe9-48da-b064-190d2a8166ca
Icc_2phg_012 = [Icc_2phg_0; Icc_2phg_1; Icc_2phg_2;;]

# ╔═╡ 2a7d0c6a-e59a-46c9-9c22-f90379cbc46b
Icc_2phg_abc = A * Icc_2phg_012

# ╔═╡ babc68db-8b16-4e7e-b7ee-da21b048ce10
dv(Icc_2phg_abc)

# ╔═╡ e64fde16-c53a-4e6b-b33e-b9907e042daa
dv(Icc_2phg_abc * Ibase)

# ╔═╡ f9ec36bf-6083-4a9d-9241-16ff2e6fbcd7
md"""
Corrente de falta no condutor neutro:
"""

# ╔═╡ ccb6dbdd-d6bd-4f78-9445-c5b67f406b3c
Icc_2phg_n = Icc_2phg_abc[2] + Icc_2phg_abc[3]

# ╔═╡ 9d54f2fd-83f3-469d-96d0-5b4ece42981e
abs(Icc_2phg_n)

# ╔═╡ 683d1289-36da-4dab-a88f-ce5f8d632dc3
abs(Icc_2phg_n * Ibase)

# ╔═╡ f0ac233f-c3d7-48de-89ef-b78444e43514
md"""
Cálculo das tensões de sequência nos terminais da falta:
"""

# ╔═╡ e43b4d19-efe4-4ce5-b525-4eadece7bbcb
V012_2phg = [0.0; 1.05; 0.0] - diagm([Zth2_0; Zth2_1; Zth2_2]) * Icc_2phg_012

# ╔═╡ 135b4bf9-4ca4-495d-ba40-a9e6056d500b
begin
	Vabc_2phg = A * V012_2phg
	dv(Vabc_2phg)
end

# ╔═╡ 41fe4d65-0883-4aef-8a27-fd2089de94b3


# ╔═╡ aeb2f2dc-b9d6-4727-a537-41149a904664
md"""
## Example 9.7:  Single line-to-ground short-circuit calculations using Zbus 0, Zbus 1, and Zbus 2

!!! important "Question"
	Faults at buses 1 and 2 for the three-phase power system given in Example 9.1 are of interest. The prefault voltage is 1.05 per unit. Prefault load current is neglected.

	- (a) Determine the per-unit zero-, positive-, and negative-sequence bus impedance matrices. 
	- (b) Find the subtransient fault current in per-unit for a bolted single line-to-ground fault current from phase a to ground at bus 1, and 
	- (c) at bus 2.

	Find the per-unit line-to-ground voltages at:	
	
	- (d) bus 1, and
	- (e) bus 2 during the single line-to-ground fault at bus 1.
"""

# ╔═╡ c882d846-a5c1-4d70-8351-ea27c6eb0deb
md"""
Montagem das matrizes de admitância nodal de sequência:
"""

# ╔═╡ 595e55f7-7cf7-48ba-b614-30f868f7e953
Ybus_0 = -1.0im * [20.0 0.0; 
				   0.0  4.0] 

# ╔═╡ 6008c5d3-bd85-4428-9fb8-8ff25d4d4078
Zbus_0 = inv(Ybus_0)

# ╔═╡ 9ad1de87-79ca-4d67-84d8-b0e420172310
md"""
Note that the transformer leakage reactances and the zero-sequence transmission-line reactance have no effect on Zbus 0.

The transformer D connections block the flow of zero-sequence current from the transformers to bus 1 and 2.
"""

# ╔═╡ 095b0fcc-407a-4be9-bb1e-405535ada979
Ybus_1 = -1.0im * [9.9454 -3.2787; -3.2787 8.2787]

# ╔═╡ 9f0dae14-b6df-4b98-b6fd-3720f54127d8
Zbus_1 = inv(Ybus_1)

# ╔═╡ 3171bad9-9151-4cd9-920e-a1063de2352f
Ybus_2 = -1.0im * [9.1611 -3.2787; -3.2787 8.0406]

# ╔═╡ d4fb0942-2fe2-4850-9172-38e84e576ebe
Zbus_2 = inv(Ybus_2)

# ╔═╡ bc9596d7-3cb4-4725-ad3f-d45225b0d541
md"""
### Calculando o curto-circuito monofásico na barra 1
"""

# ╔═╡ 7f7a2c85-042e-4729-8bc6-9e003e5a764f
Icc_1ph_1_ = Icc_1ph_2_ = Icc_1ph_0_ = 1.05 / (Zbus_1[1, 1] + Zbus_2[1, 1] + Zbus_0[1, 1])

# ╔═╡ 6dee2f88-9e60-4c04-a423-3487e00e8c30
Icc_1ph_012_ = [Icc_1ph_1_; Icc_1ph_2_;  Icc_1ph_0_;;]

# ╔═╡ 716b20f7-0207-4bb7-be6e-e4e5805fb9ce
Icc_1ph_abc_ = A * Icc_1ph_012_

# ╔═╡ 7543bf17-8094-4989-ad6e-8c51cbb561dc
dv(Icc_1ph_abc_)

# ╔═╡ 29c20476-24f8-4a98-80f4-c04f3a6157ac
md"""
Para o curto-circuito na barra 2
"""

# ╔═╡ 17dc49e9-30cd-4bd1-aeee-a9d21189a115


# ╔═╡ d20470f8-db45-4eb4-9483-5cdbd189f259
Icc_1ph_1_k = Icc_1ph_2_k = Icc_1ph_0_k = 1.05 / (Zbus_1[2, 2] + Zbus_2[2, 2] + Zbus_0[2, 2])

# ╔═╡ 7cc8c2dd-e296-499f-9a3e-11412cedae1f
Icc_1ph_012_k = [Icc_1ph_1_k; Icc_1ph_2_k;  Icc_1ph_0_k;;]

# ╔═╡ e9f03b27-dbc9-447c-b9e9-5ce1fa8fff41
Icc_1ph_abc_k = A * Icc_1ph_012_k

# ╔═╡ b6f1eebd-cf8c-45fa-bb56-654b594d45b1
dv(Icc_1ph_abc_k)

# ╔═╡ fd733333-f88c-4772-ad32-1214ad51a0b4
7200 / 120

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.59"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.5"
manifest_format = "2.0"
project_hash = "e41793cbd1124ea5d05573eb874098b20a960e1d"

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
# ╠═e08927f0-cf6e-4d76-9794-df9c224f6972
# ╠═84181e0c-5e4e-11ef-12c9-b9603ac2cbdb
# ╠═d88c87b7-b654-4ecf-ab5e-9bef6bbec900
# ╠═300b0cd1-b493-488a-806e-9b449a40e625
# ╟─e705c595-b289-4b60-891f-fb1183e615fd
# ╟─76553e35-2f5a-46a7-9db1-c72376c0177f
# ╟─e4c4eaf0-c5c8-4f7d-8f5c-44d10fe9134b
# ╟─e7cea086-700f-4061-91c8-5363e3d9faed
# ╟─51c3e84d-8df3-4fd8-a1f4-a9513c2209c3
# ╟─5f20c381-3a1e-4709-9fcd-dfdcab9c3363
# ╠═12971eca-f21e-4c3a-917e-ec0af879da60
# ╠═40a32e16-8787-4b98-909c-45c93b7672d7
# ╠═360e119e-87ae-43f4-8f54-4ad21166002e
# ╠═984284fd-67c8-4f3f-8882-f52914321867
# ╠═893d06a2-e16f-4e89-b3aa-9605f95956d9
# ╟─feeea186-3682-43f1-80eb-22d640f458a8
# ╟─edbe576f-a0cb-48f2-b48c-b39a6b4ae902
# ╟─380fdaeb-2aa1-427b-a557-3fcbd1ae48e4
# ╠═d9782a87-db49-4110-b1d6-74bd71a78628
# ╠═7735ff4d-a08e-468c-b781-89d5135b150b
# ╠═689093f4-1ebf-4915-b4b4-f5721efd6c6e
# ╟─41a555ee-9b09-430a-8c00-de073805b55d
# ╠═6f2d7d97-743a-4052-a26a-b45a9be08f78
# ╠═e9bfb3db-4544-4b50-a8d4-b84d3dc7f0ee
# ╟─1546e93c-4f2d-4e4b-9a2d-71d3dfcefb24
# ╟─b8c4bdc8-1605-48f3-8ec9-ddaa998eaf24
# ╠═193994a5-509f-442a-af92-8a2e4de54ddd
# ╟─2f70870c-c553-4ddc-a481-173647d0d71e
# ╠═efddf058-44f2-4ca4-8ef6-6b0192aa88e8
# ╟─380f0924-91a1-4703-9e70-eb62f7ac8c65
# ╠═bb9541aa-7a12-4b19-ae80-06e532af8d08
# ╟─88287684-16d5-4c8e-a553-90a9fe24a7e9
# ╠═58e4de6d-9b23-41a9-a39d-132cd217443c
# ╟─618ccc1e-8006-403b-981c-5b39d42f7fb3
# ╠═09ca922b-d8e1-44c9-ac56-143ad66c6a1d
# ╠═ea84f5cd-7ae6-4c39-9084-62854c655dff
# ╠═f0da4182-4506-4423-910d-5fdb4e070a3c
# ╟─bbe49832-b653-4e5d-becc-ddce8cc32b35
# ╟─2200f673-3d81-4d72-b987-23087ded3851
# ╟─925adf5e-bdf8-4377-bd75-2aca5c8c3ad6
# ╟─84b51bc4-557f-43ef-8c8c-2c0e83997e49
# ╠═9ea8e62f-f8dc-4e6d-8330-661f48901ac8
# ╠═ccda01a4-c3e5-45e9-89d5-09a76009c4db
# ╠═a08d7997-57d0-4c9b-a587-1c0059dae14a
# ╠═7f0861bc-13b3-49fc-beea-5773635aa2e2
# ╟─5d07a2a4-9199-4165-a338-f0bfd419ef4f
# ╠═c53172dd-ed14-4014-8a26-628e266909d9
# ╠═406fa646-dbaf-43e9-b085-7f469e213865
# ╟─17e5fa14-9932-4674-98e9-410372e406b1
# ╠═ce0556fd-400b-4edf-95f6-e4134a07b69f
# ╠═c98c4c36-b80e-434e-bbb7-3345e418c354
# ╟─54ac7a03-15fd-4aba-9acc-8ae6949e6639
# ╟─6fb52297-c521-44fc-9568-10b6fe883cb3
# ╟─c31ca637-7e19-41f6-b694-a0a6fcec193b
# ╟─44c97ab1-9f6c-4f52-93bb-c73451cce2ed
# ╠═b27c134f-6736-465f-b473-5060f2099cf9
# ╠═c72b4480-7912-4fb6-aef6-399c07489e51
# ╠═c4de09dd-abac-4606-b702-8b795c12075d
# ╟─bf4eeeb3-7d07-4644-aa3f-3a47cce827c4
# ╠═e35721b0-2f1e-4205-a66f-2a048d4e7edb
# ╠═44613696-6087-4c8b-8065-a90ea846c81b
# ╠═0ac51b1c-df6b-43bd-8d30-e535d3b91218
# ╠═168b3cd6-25d8-4ee9-9ab7-824ba3f66777
# ╟─7a893665-21a8-4a10-bec2-81db48e4cf0b
# ╠═c668ce3e-c6eb-4d9c-a779-ecc11ce44b73
# ╠═4edbadc7-8ec1-4068-b36c-74739e40816f
# ╟─eba51237-c21c-476d-8ca7-2848ce012aac
# ╟─607e74ad-ef8b-4c5c-9e62-da9129593fac
# ╟─b660051c-ea82-4c62-ae44-2452de3197da
# ╟─9c220de1-3a66-4b94-93af-596c8fa6c2f3
# ╠═efe8408e-c652-4c22-b68d-d90dd3066213
# ╟─b30cb39e-dd05-4b46-b1df-f64ab36fa1fd
# ╠═3ab747ea-ab51-4280-adbb-8932e2ed3b9f
# ╠═ebcdd089-c1f8-4a75-8c18-f47901e8ee2b
# ╟─cc8b842d-eace-423d-b6dd-0286e6c08527
# ╠═5e40efa2-0ded-4259-b4bd-2da7e9ff0f32
# ╠═7b6a5ab0-b90c-4818-a4f2-a214b4f2c25c
# ╠═bd7f495f-43f4-48cf-83ae-32059462bd77
# ╟─5baeaa55-31ef-4d1e-8073-28dc761d4429
# ╠═8cb55546-4282-4d7d-a1db-b7ed55451709
# ╠═7a97c09a-da58-4ec8-8551-67fa1bafcb5d
# ╟─7d1e1286-bea2-4ef8-9937-10c0495127a6
# ╟─cb0173e8-d981-4e0c-ab3d-d71c5cfd18bc
# ╟─059ad344-23d4-4efd-b81b-d801d9a1fe4f
# ╟─9743cc24-f97a-458f-95d0-a0f8d0618fda
# ╠═d225fe5c-6a09-4404-ac9e-241b653fcf6b
# ╠═c03a2907-9bfb-4c63-b4d2-d0e9d89df250
# ╠═99317f2d-2fe9-48da-b064-190d2a8166ca
# ╠═2a7d0c6a-e59a-46c9-9c22-f90379cbc46b
# ╠═babc68db-8b16-4e7e-b7ee-da21b048ce10
# ╠═e64fde16-c53a-4e6b-b33e-b9907e042daa
# ╟─f9ec36bf-6083-4a9d-9241-16ff2e6fbcd7
# ╠═ccb6dbdd-d6bd-4f78-9445-c5b67f406b3c
# ╠═9d54f2fd-83f3-469d-96d0-5b4ece42981e
# ╠═683d1289-36da-4dab-a88f-ce5f8d632dc3
# ╟─f0ac233f-c3d7-48de-89ef-b78444e43514
# ╠═e43b4d19-efe4-4ce5-b525-4eadece7bbcb
# ╠═135b4bf9-4ca4-495d-ba40-a9e6056d500b
# ╟─41fe4d65-0883-4aef-8a27-fd2089de94b3
# ╟─aeb2f2dc-b9d6-4727-a537-41149a904664
# ╟─c882d846-a5c1-4d70-8351-ea27c6eb0deb
# ╠═595e55f7-7cf7-48ba-b614-30f868f7e953
# ╠═6008c5d3-bd85-4428-9fb8-8ff25d4d4078
# ╟─9ad1de87-79ca-4d67-84d8-b0e420172310
# ╠═095b0fcc-407a-4be9-bb1e-405535ada979
# ╠═9f0dae14-b6df-4b98-b6fd-3720f54127d8
# ╠═3171bad9-9151-4cd9-920e-a1063de2352f
# ╠═d4fb0942-2fe2-4850-9172-38e84e576ebe
# ╟─bc9596d7-3cb4-4725-ad3f-d45225b0d541
# ╠═7f7a2c85-042e-4729-8bc6-9e003e5a764f
# ╠═6dee2f88-9e60-4c04-a423-3487e00e8c30
# ╠═716b20f7-0207-4bb7-be6e-e4e5805fb9ce
# ╠═7543bf17-8094-4989-ad6e-8c51cbb561dc
# ╟─29c20476-24f8-4a98-80f4-c04f3a6157ac
# ╠═17dc49e9-30cd-4bd1-aeee-a9d21189a115
# ╠═d20470f8-db45-4eb4-9483-5cdbd189f259
# ╠═7cc8c2dd-e296-499f-9a3e-11412cedae1f
# ╠═e9f03b27-dbc9-447c-b9e9-5ce1fa8fff41
# ╠═b6f1eebd-cf8c-45fa-bb56-654b594d45b1
# ╠═fd733333-f88c-4772-ad32-1214ad51a0b4
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
