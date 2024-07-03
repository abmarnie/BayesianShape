using HDF5, Printf, Plots, Plots.Measures

## TODO: Add a comment to the top describing which plots in the paper these generate.

include("../src/fourierBasis.jl");
include("../src/computeRadii.jl");
include("../src/plotSave.jl");

# Labels.
r_min_lab = "\$r_{min}\$"; # Min inner boundary.
r_max_lab = "\$r_{max}\$"; # Max inner boundary.
r_out_lab = "\$\\Gamma^{\\rm outer}\$ (radius \$R\$)"; # Outer boundary.
boundary_lab = "\$\\Gamma_b^{\\rm inner}\$ (radius \$c(b_0+b)\$)";  # Inner boundary.

r_out_lab_rad = "\$R\$";             # Outer boundary.
boundary_lab_rad = "\$c(b_0+b)\$";   # Inner boundary.
unsquashed_lab = "\$b_0+b\$";        # Unsquashed boundary.

# Get `a0`, and run setup to get radius squash (uses default :/ ).
include("../scenarios/svsector/setup.jl");

file = "pubs/2022_radius_diagram.h5";

if isfile(file)
	println("Reading from: $file")
	f = h5open(file, "r")
	r = read(f, "r")
	b = read(f, "b")
	angles = read(f, "th")
	r_min = read(f, "rMin")
	r_max = read(f, "rMax")
	r_out = read(f, "rOut")
	close(f)
else
	# Radii.
	r_out = 2.0

	# Angles.
	const angles_stride = pi / 180
	angles = 0:angles_stride:2*pi

	fb = fourierBasis(unkDim, angles)

	r = zeros(2 * unkDim)
	b = zeros(2 * unkDim)
	samples = zeros(2 * unkDim)

	stop_iter = false
	count = 1

	while !stop_iter
		global stop_iter, count, r, b, samples # Must be global, because of Julia jank.

		println("Attempt $count")

		# Draw randomly from the prior.
		samples = rand(mcmcP.prior)

		r = computeRadii(samples, angles)

		# Compute unsquashed radii.
		b = a0 .+ fb' * samples

		pct_at_boundary = (sum(r .== r_max) + sum(r .== r_min)) / length(r)

		stop_iter = (pct_at_boundary < 0.25) && (count <= 50)
		count += 1
	end
end

# Plot (radial).
sz = 500;
p = plot(proj = :polar, size = (sz, sz), margin = 5Plots.mm);
plot!(p, angles, [r_min], ls = :dash, lab = r_min_lab);
plot!(p, angles, r, lab = boundary_lab);
plot!(p, angles, [r_max], ls = :dash, lab = r_max_lab);
plot!(p, angles, [r_out], c = :black, lab = r_out_lab);
plotSave(p, "pubs/2022_radius_diagram");

# Plot (rectangular).
# The lc is set to manually match above.
th_deg = angles * 180 / pi;
sz = 500;
p = plot(size = (sz, sz), margin = 5Plots.mm, xlab = "Angle, \$\\phi\$ (Degrees)", ylab = "Radius", xticks = 0:90:360);
hline!(p, [r_min], lc = 1, ls = :dash, lab = r_min_lab);
plot!(p, th_deg, b, lc = 4, lab = unsquashed_lab);
plot!(p, th_deg, r, lc = 2, lab = boundary_lab_rad);
hline!(p, [r_max], lc = 3, ls = :dash, lab = r_max_lab);
hline!(p, [r_out], c = :black, lab = r_out_lab_rad);
plot!(legend = :bottomright);
plotSave(p, "pubs/2022_radius_diagram_rect");

if !isfile(file)
	f = h5open(file, "w")
	write(f, "r", r)
	write(f, "b", b)
	write(f, "samples", samples)
	write(f, "th", collect(angles))
	write(f, "rMin", r_min)
	write(f, "rMax", r_max)
	write(f, "rOut", r_out)
	close(f)
	println("Wrote data to: $file")
end
