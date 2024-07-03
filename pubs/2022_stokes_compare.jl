using HDF5, Printf, Plots, Plots.Measures

## TODO: Add a comment to the top describing which plots in the paper these generate.

include("../src/fourierBasis.jl");
include("../src/computeRadii.jl");
include("../src/plotSave.jl");

scenario = "vortsensor";

# Run setup to get radius squash (uses default :/ ).
include("../scenarios/$(scenario)/setup.jl");

dir = "/projects/SIAllocation/stokes/$(scenario)";
files = ["$(scenario)_0$(i).h5" for i ∈ 17:22];

# Angles and fourier basis with which we'll take the inner product of.
const angles_stride = pi / 180;
angles = 0:angles_stride:2*pi;
fourier_basis = exp.(im * angles[1:end-1]);

# Lists.
nsamp = h5read(dir * "/" * files[1], "sampComplete"); # Assume the same for all files.
areas = zeros(nsamp, 0);
phis = zeros(nsamp, 0);

pphi = plot(xlab = "Samples", ylab = "Mean Angle (Radians)");
parea = plot(xlab = "Samples", ylab = "Mean Area");
for file in files
	f = h5open(dir * "/" * file)
	samples = read(f, "samples")
	close(f)

	r = computeRadii(samples, angles)
	area = angles_stride * sum(0.5 * r[:, 1:end-1] .^ 2, dims = 2)[:]
	fftr = r[:, 1:end-1] * fourier_basis
	phi = atan.(imag.(fftr), real.(fftr))

	# Running average.
	mphi  = cumsum(phi) ./ (1:length(phi))
	marea = cumsum(area) ./ (1:length(area))

	# Running stdev.
	sphi  = sqrt.(cumsum((phi .- mean(phi)) .^ 2) ./ (1:length(phi)))
	sarea = sqrt.(cumsum((area .- mean(area)) .^ 2) ./ (1:length(area)))

	# Plot attributes.
	label = "Chain $i"

	plot!(pphi, 1:length(mphi), mphi, lab = label)
	plot!(parea, 1:length(marea), marea, lab = label)

	global areas = hcat(areas, area)
	global phis  = hcat(phis, phi)
end

plot!(pphi, leg = :bottomright);
plot!(parea, leg = :bottomright, ylims = (0, 5));

plotSave(pphi, "pubs/2022_stokes_compare_angle");
plotSave(parea, "pubs/2022_stokes_compare_area");

# See Gelman BDA3, pg284.
function statGelmanRubin(psi)
	n, m    = size(psi)
	psibar  = mean(psi)
	psibarj = mean(psi, dims = 1)[:]

	B = n / (m - 1) * norm(psibar .- psibarj)^2

	s = zeros(m)
	for j ∈ 1:m
		for i ∈ 1:n
			s[j] += 1 / (n - 1) * (psi[i, j] .- psibarj[j])^2
		end
	end
	W = mean(s)

	varp = (W * (n - 1) + B) / n
	Rhat = sqrt(varp / W)

	return Rhat
end

#print some stats
grArea = statGelmanRubin(areas);
grPhi = statGelmanRubin(phis);
println("Gelman-Rubin for area is : $grArea");
println("Gelman-Rubin for angle is: $grPhi");

for k ∈ 1:8
	comp = zeros(nsamp, length(files))
	for j ∈ 1:length(files)
		comp[:, j] = h5read(dir * "/" * files[j], "samples")[:, k]
	end
	grComp = statGelmanRubin(comp)
	println("Gelman-Rubin for component $k is : $grComp")
	if grComp > 1.1
		println("means for component $k are:")
		display(mean(comp, dims = 1))
	end
end
