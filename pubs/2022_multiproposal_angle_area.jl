using HDF5, Printf, Plots, Plots.Measures

# TODO: Add a comment to the top describing which plots in the paper these generate.
# TODO: Replace 18 and 29 magic numbers with named constants.

include("../src/fourierBasis.jl");
include("../src/computeRadii.jl");
include("../src/plotSave.jl");

scenario = "svsector";

# Run setup to get radius squash (uses default :/ ).
include("../scenarios/$(scenario)/setup.jl");

dir = "/projects/SIAllocation/stokes/$(scenario)";

const n_chains_per_case = 2
files = ["$(scenario)_0$(i).h5" for i ∈ 18:n_chains_per_case:29];

# Build a list of beta values to color the figure.
betas = [h5read(dir * "/" * file, "final_mcmc/beta") for file in files];
unique!(betas);

# Angles and fourier basis with which we'll take the inner product of.
const angles_stride = pi / 180;
angles = 0:angles_stride:2*pi;
fourier_basis = exp.(im * angles[1:end-1]);

# Start the plots.
pphi = plot(xlab = "Samples", ylab = "Mean Angle (Radians)");
parea = plot(xlab = "Samples", ylab = "Mean Area");

for file in files
	f = h5open(dir * "/" * file)
	local mcmc = read(f, "mcmc")
	beta = read(f, "final_mcmc/beta")
	samples = read(f, "samples")
	close(f)

	rho = sqrt(1 - beta^2)
	r = computeRadii(samples, angles)
	area = angles_stride * sum(0.5 * r[:, 1:end-1] .^ 2, dims = 2)[:]
	fftr = r[:, 1:end-1] * fourier_basis
	phi = atan.(imag.(fftr), real.(fftr))

	# Running average.
	mphi = cumsum(phi) ./ (1:length(phi))
	marea = cumsum(area) ./ (1:length(area))

	# Running stdev.
	sphi = sqrt.(cumsum((phi .- mean(phi)) .^ 2) ./ (1:length(phi)))
	sarea = sqrt.(cumsum((area .- mean(area)) .^ 2) ./ (1:length(area)))

	# Plot attributes.
	# MCMC type determines label and line style.
	mcmc_sp = split(mcmc, "|")
	if mcmc_sp[1] == "pcn"
		label = "pCN" # "Vanilla pCN";
		ls = :dash
	else
		label = "$(mcmc_sp[3]) Proposals"
		ls = :solid
	end

	label *= @sprintf(", \$\\rho\$=%5.3f", rho)

	# HACK: Only include unique labels -- use empty string if label is already used.
	labels_so_far = [p.plotattributes[:label] for p in pphi.series_list]
	label = (label in labels_so_far) ? "" : label

	# Set color based on `beta`.
	col = findall(==(beta), betas)[1]

	plot!(pphi, 1:length(mphi), mphi, c = col, ls = ls, lab = label)
	plot!(parea, 1:length(marea), marea, c = col, ls = ls, lab = label)
end

plot!(pphi, leg = :bottomright, ylims = (-1.00, 1.00));
plot!(parea, leg = :bottomright, ylims = (3.65, 3.85));

plotSave(pphi, "pubs/2022_multiproposal_angle");
plotSave(parea, "pubs/2022_multiproposal_area");
