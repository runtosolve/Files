### A Pluto.jl notebook ###
# v0.18.0

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

# ╔═╡ f9d65a36-1708-4a60-b91d-c23fdef2788b
begin
	import Pkg
	Pkg.add(url="https://github.com/runtosolve/CUFSM.jl.git")
	Pkg.add(url="https://github.com/runtosolve/Geometry.jl.git")
	Pkg.add(url="https://github.com/runtosolve/InternalForces.jl.git")
	Pkg.add(url="https://github.com/runtosolve/S100AISI.jl.git")
	Pkg.add(url="https://github.com/runtosolve/ScrewConnections.git")
	Pkg.add(url="https://github.com/runtosolve/SectionProperties.jl.git")
	Pkg.add(url="https://github.com/runtosolve/ThinWalledBeam.jl.git")
	Pkg.add(url="https://github.com/runtosolve/ThinWalledBeamColumn.jl.git")
	Pkg.add(url="https://github.com/runtosolve/PurlinLine.jl.git")
	Pkg.add(url="https://github.com/runtosolve/TopHatDesigner.jl.git")
	Pkg.add("ImageMagick")
	
	Pkg.activate(".")
	Pkg.instantiate()



	using TopHatDesigner, PlutoUI, Images, Dates, CSV, DataFrames, UrlDownload

	purlin_data = DataFrame(urldownload("https://raw.githubusercontent.com/runtosolve/Files/main/TH/database/Purlins.csv"))

	top_hat_data = DataFrame(urldownload("https://raw.githubusercontent.com/runtosolve/Files/main/TH/database/TopHats.csv"))

	existing_deck_data = DataFrame(urldownload("https://raw.githubusercontent.com/runtosolve/Files/main/TH/database/Existing_Deck.csv"))

	new_deck_data = DataFrame(urldownload("https://raw.githubusercontent.com/runtosolve/Files/main/TH/database/New_Deck.csv"))

	logo = urldownload("https://www.tophatframing.com/assets/images/tophat-horz.png")

end;

# ╔═╡ 83c92777-b410-4a95-9ba7-54d725e2ed9b
logo

# ╔═╡ c2353bb3-ee8c-4d55-9447-470427c22b06
@bind project_details TextField((30,5); default="Project details")

# ╔═╡ 96f90537-0b4d-4d48-927b-01492e3789ef
@bind report_date DateField(default=today())

# ╔═╡ 99299f0c-30ee-4807-a7a2-d4509b4680ab
md" ## Build existing roof system."

# ╔═╡ 1d9b00aa-7f6b-4f7e-9da8-e1f0b1ace647
md"""
Purlin sizes $(@bind purlin_type_1 Select(["Z8x2.5 054"; purlin_data[:, 1]]))
$(@bind purlin_type_2 Select(["none"; purlin_data[:, 1]]))
"""

# ╔═╡ 5d180e53-27aa-4bb2-9ab1-81cc2737ab3b
purlin_spans = (25.0)  #ft

# ╔═╡ 111ad395-5261-41e3-bd71-0a08ebe97119
purlin_size_span_assignment = (1)

# ╔═╡ 710a97bb-6cd1-457c-a352-23428408de55
purlin_laps = ()

# ╔═╡ dde7f4c2-2212-4244-a78b-8fe12b6c8d0e
purlin_spacing = 5.0  #ft

# ╔═╡ 1a1727db-c828-4317-8ccd-2be491ee48c0
frame_flange_width = 10.0  #in

# ╔═╡ 9eebfba0-7913-40fd-bda5-b1d1178c741a
roof_slope = 1/12

# ╔═╡ e20e6735-ae8a-4ae0-99bb-f563a602afbc
md" Existing roof deck type $(@bind existing_deck_type Select(existing_deck_data[:, 1]))"

# ╔═╡ dd6bbcfe-b781-4bf5-8485-1c4b25380ebc
md" ## Calculate existing roof system strength."

# ╔═╡ 45651a24-ebc3-4ab4-b7a9-ea5a1ab7fa7f
purlin_line, purlin_line_uplift = UI.existing_roof_UI_mapper(purlin_spans, purlin_laps, purlin_spacing, roof_slope, purlin_data, existing_deck_type, existing_deck_data, frame_flange_width, purlin_type_1, purlin_type_2, purlin_size_span_assignment);

# ╔═╡ 32fbe54d-b709-4f80-bf10-913abcd63e41
UI.plot_purlin_geometry(purlin_line.inputs.cross_section_dimensions[1][2], purlin_line.cross_section_data[1].node_geometry[:,1], purlin_line.cross_section_data[1].node_geometry[:,2], roof_slope)

# ╔═╡ d2246077-5cf8-4f00-81af-9e5922bec619
md"**Existing roof system downward (gravity) strength = $(round(purlin_line.applied_pressure*1000*144, digits=1)) psf**"

# ╔═╡ c66ea01f-f1ea-4071-a490-730d02db486a
md"**Existing roof system uplift strength = $(round(purlin_line_uplift.applied_pressure*1000*144, digits=1)) psf**"

# ╔═╡ ec8da739-a36b-4683-98a0-34572d403660
purlin_line

# ╔═╡ da993059-338a-42c2-b1c5-34b1b4434967
purlin_line_uplift

# ╔═╡ f67cf06d-2fe6-401f-ab71-1ecea9aa373f
md" ## Add TopHat and new roof."

# ╔═╡ 6cebdfda-7d1c-4df6-9e64-b6afad1c7f6d
md" TopHat size $(@bind top_hat_type Select(top_hat_data[:, 1]))"
		

# ╔═╡ bfac357a-4b13-471e-84e8-edf4272065ba
md" New roof deck type $(@bind new_deck_type Select(new_deck_data[:, 1]))"

# ╔═╡ edb7ac4d-3b87-4c73-bb3a-97711bc9c0dd
top_hat_purlin_line, top_hat_purlin_line_uplift = UI.retrofit_UI_mapper(purlin_line, top_hat_data, top_hat_type, existing_deck_type, existing_deck_data, new_deck_type, new_deck_data);

# ╔═╡ 444cb829-2b66-48a8-a052-b16f5f5529b8
UI.plot_top_hat_purlin_geometry(top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[1][1], top_hat_purlin_line.purlin_cross_section_data[1].node_geometry[:,1], top_hat_purlin_line.purlin_cross_section_data[1].node_geometry[:,2], roof_slope, top_hat_purlin_line.top_hat_purlin_cross_section_data[1].node_geometry[:,1], top_hat_purlin_line.top_hat_purlin_cross_section_data[1].node_geometry[:,2], top_hat_purlin_line.top_hat_cross_section_data[1].node_geometry[:,1], top_hat_purlin_line.top_hat_cross_section_data[1].node_geometry[:,2])

# ╔═╡ 680fdccf-76fc-451e-8caa-bdb6ff359818
UI.plot_net_section_top_hat_purlin_geometry(purlin_line.inputs.cross_section_dimensions[1][2], top_hat_purlin_line.purlin_cross_section_data[1].node_geometry[:,1], top_hat_purlin_line.purlin_cross_section_data[1].node_geometry[:,2], roof_slope, top_hat_purlin_line.top_hat_purlin_net_cross_section_data[1].node_geometry[:,1], top_hat_purlin_line.top_hat_purlin_net_cross_section_data[1].node_geometry[:,2])

# ╔═╡ 4d47c72c-52a9-4caa-9516-446296e73443
md" ## Calculate retrofitted roof system strength."

# ╔═╡ 8fe61b54-e9e1-4aa6-be50-b0812d86a1cf
md"**Retrofitted roof system downward (gravity) strength = $(round(top_hat_purlin_line.applied_pressure*1000*144, digits=1)) psf**"

# ╔═╡ af2fc099-ded8-430f-bd28-a82af438a280
md"**Retrofitted roof system uplift strength = $(round(top_hat_purlin_line_uplift.applied_pressure*1000*144, digits=1)) psf**"

# ╔═╡ e73dde74-872b-42fa-8159-42f354468861
top_hat_purlin_line

# ╔═╡ 1444f7b3-bcf4-4fbc-8126-d069345d8389
top_hat_purlin_line_uplift

# ╔═╡ Cell order:
# ╟─f9d65a36-1708-4a60-b91d-c23fdef2788b
# ╟─83c92777-b410-4a95-9ba7-54d725e2ed9b
# ╟─c2353bb3-ee8c-4d55-9447-470427c22b06
# ╟─96f90537-0b4d-4d48-927b-01492e3789ef
# ╟─99299f0c-30ee-4807-a7a2-d4509b4680ab
# ╟─1d9b00aa-7f6b-4f7e-9da8-e1f0b1ace647
# ╠═5d180e53-27aa-4bb2-9ab1-81cc2737ab3b
# ╠═111ad395-5261-41e3-bd71-0a08ebe97119
# ╠═710a97bb-6cd1-457c-a352-23428408de55
# ╠═dde7f4c2-2212-4244-a78b-8fe12b6c8d0e
# ╠═1a1727db-c828-4317-8ccd-2be491ee48c0
# ╠═9eebfba0-7913-40fd-bda5-b1d1178c741a
# ╟─e20e6735-ae8a-4ae0-99bb-f563a602afbc
# ╟─32fbe54d-b709-4f80-bf10-913abcd63e41
# ╟─dd6bbcfe-b781-4bf5-8485-1c4b25380ebc
# ╟─45651a24-ebc3-4ab4-b7a9-ea5a1ab7fa7f
# ╟─d2246077-5cf8-4f00-81af-9e5922bec619
# ╟─c66ea01f-f1ea-4071-a490-730d02db486a
# ╠═ec8da739-a36b-4683-98a0-34572d403660
# ╠═da993059-338a-42c2-b1c5-34b1b4434967
# ╟─f67cf06d-2fe6-401f-ab71-1ecea9aa373f
# ╟─6cebdfda-7d1c-4df6-9e64-b6afad1c7f6d
# ╟─bfac357a-4b13-471e-84e8-edf4272065ba
# ╟─edb7ac4d-3b87-4c73-bb3a-97711bc9c0dd
# ╟─444cb829-2b66-48a8-a052-b16f5f5529b8
# ╟─680fdccf-76fc-451e-8caa-bdb6ff359818
# ╟─4d47c72c-52a9-4caa-9516-446296e73443
# ╟─8fe61b54-e9e1-4aa6-be50-b0812d86a1cf
# ╟─af2fc099-ded8-430f-bd28-a82af438a280
# ╠═e73dde74-872b-42fa-8159-42f354468861
# ╠═1444f7b3-bcf4-4fbc-8126-d069345d8389
