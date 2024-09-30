### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 9b2a16fa-5782-4202-a60c-e8e73e893ee5
begin
	using LinearAlgebra
	LA = LinearAlgebra
end

# ╔═╡ e78573f7-7a77-4b3b-9fa1-04f516bccd31
"""
Function to extract the main and secondary diagonals of a band matrix in form of vectors.

M ::AbstractArray;
bandnum ::Int64 = -1

	# Examples
    ```jldoctest
    julia> M = [1 -1 0 0 0; 2 1 -1 0 0; 0 2 1 -1 0; 0 0 2 1 -1; 0 0 0 2 1]
	julia> bandMatrixTrim(M)
	julia> Any[[2, 2, 2, 2], [1, 1, 1, 1, 1], [-1, -1, -1, -1]]
    ```
"""
function bandMatrixTrim(M::AbstractArray;bandnum::Int=3)
	diag_vector::Array = []
	s = bandnum÷2
	for k in collect(-s:+s)
			push!(diag_vector,diag(M,k))
	end
	
	return diag_vector
end
	

# ╔═╡ 4f5acabc-5687-4ec2-bb0b-7fcb20cac011
begin
	A = LA.Tridiagonal([2, 2, 2, 0],[1, 2, 3, 4, 5], [3, 3, 3, 3])
	M = [1 -1 0 0 0; 2 1 -1 0 0; 0 2 1 -1 0; 0 0 2 1 -1; 0 0 0 2 1]
	println(bandMatrixTrim(M))
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "ac1187e548c6ab173ac57d4e72da1620216bce54"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"
"""

# ╔═╡ Cell order:
# ╠═9b2a16fa-5782-4202-a60c-e8e73e893ee5
# ╠═e78573f7-7a77-4b3b-9fa1-04f516bccd31
# ╠═4f5acabc-5687-4ec2-bb0b-7fcb20cac011
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
