using PyPlot, HDF5
include("RiemannSolver.jl")

N = 2000
xmin = -0.5
xmax = 0.5
x = collect(LinRange(xmin, xmax, N))
save_figs = true #if this is false the figures will use show(), otherwise they will be saved locally in 'figs' dir

#create directory for figures
if save_figs && !isdir("figs")
	mkdir("figs")
end

struct TestCase
	name::String
	t::Float64
	left::HydroStatus
	right::HydroStatus
end

function calculate_and_plot(test::TestCase, h5group)
	name = test.name

	println("Time for the test $name:")

	#dummy run to compile all code
	sample_riemann(x, test.t, test.left, test.right)

	@time profiles = sample_riemann(x, test.t, test.left, test.right)
	density = [status.rho for status in profiles]
	velocity = [status.u for status in profiles]
	pressure = [status.p for status in profiles]
	energy = pressure ./ density ./ (test.left.gamma - 1.)
#ASSUMES CV = 1.0
	cv = 1.0
	specific_entropy = cv * log.(pressure ./ (density.^test.left.gamma))
	entropy_density = density .* specific_entropy
	temperature = energy / cv



	write(create_dataset(h5group, "density", Float64, (N,)), density)
	write(create_dataset(h5group, "velocity", Float64, (N,)), velocity)
	write(create_dataset(h5group, "pressure", Float64, (N,)), pressure)
	write(create_dataset(h5group, "energy", Float64, (N,)), energy)
	write(create_dataset(h5group, "entropy_density", Float64, (N,)), entropy_density)
	write(create_dataset(h5group, "specific_entropy", Float64, (N,)), specific_entropy)
	write(create_dataset(h5group, "temperature", Float64, (N,)), temperature)

	figure()
	suptitle(name)
	subplot(231)
	plot(x, density)
	title("density")
	xlabel("x")
	grid(true)
	subplot(232)
	plot(x, velocity)
	title("velocity")
	xlabel("x")
	grid(true)
	subplot(233)
	plot(x, pressure)
	title("pressure")
	xlabel("x")
	grid(true)
	subplot(234)
	plot(x, energy)
	title("energy")
	xlabel("x")
	grid(true)
	subplot(235)
	plot(x, specific_entropy)
	title("specific entropy")
	xlabel("x")
	grid(true)
	subplot(236)
	plot(x, entropy_density)
	title("entropy density")
	xlabel("x")
	grid(true)
	tight_layout()

	if !save_figs
		show()
	else
		savefig("figs/$name.png",format="png",dpi=200)
	end
end

gamma = 7. / 5.


fid = h5open("RP.h5", "w")

#Sod Shock Tube
create_group(fid, "RP1")
create_group(fid["RP1"], "0.0")
create_group(fid["RP1"], "0.1")
create_group(fid["RP1"], "0.2")
create_group(fid["RP1"], "0.25")
attributes(fid["RP1"])["altname"] = "sod"

RP1 = TestCase("RP1-0.0",  0.0, HydroStatus(1., 0., 1., gamma), HydroStatus(0.125, 0., 0.1, gamma))
calculate_and_plot(RP1, fid["RP1/0.0"])
RP1 = TestCase("RP1-0.1",  0.1, HydroStatus(1., 0., 1., gamma), HydroStatus(0.125, 0., 0.1, gamma))
calculate_and_plot(RP1, fid["RP1/0.1"])
RP1 = TestCase("RP1-0.2",  0.2, HydroStatus(1., 0., 1., gamma), HydroStatus(0.125, 0., 0.1, gamma))
calculate_and_plot(RP1, fid["RP1/0.2"])
RP1 = TestCase("RP1-0.25",  0.25, HydroStatus(1., 0., 1., gamma), HydroStatus(0.125, 0., 0.1, gamma))
calculate_and_plot(RP1, fid["RP1/0.25"])

#Enfield 123
create_group(fid, "RP3")
create_group(fid["RP3"], "0.0")
create_group(fid["RP3"], "0.05")
create_group(fid["RP3"], "0.10")
create_group(fid["RP3"], "0.15")
attributes(fid["RP3"])["altname"] = "vaccumexpansion"

RP3 = TestCase("RP3-0.0",  0.0, HydroStatus(1., -2., 0.4, gamma), HydroStatus(1., 2., 0.4, gamma))
calculate_and_plot(RP3, fid["RP3/0.0"])
RP3 = TestCase("RP3-0.05",  0.05, HydroStatus(1., -2., 0.4, gamma), HydroStatus(1., 2., 0.4, gamma))
calculate_and_plot(RP3, fid["RP3/0.05"])
RP3 = TestCase("RP3-0.1",  0.10, HydroStatus(1., -2., 0.4, gamma), HydroStatus(1., 2., 0.4, gamma))
calculate_and_plot(RP3, fid["RP3/0.10"])
RP3 = TestCase("RP3-0.15",  0.15, HydroStatus(1., -2., 0.4, gamma), HydroStatus(1., 2., 0.4, gamma))
calculate_and_plot(RP3, fid["RP3/0.15"])

create_group(fid, "RP2")
create_group(fid["RP2"], "0.0")
create_group(fid["RP2"], "0.015")
create_group(fid["RP2"], "0.025")
create_group(fid["RP2"], "0.035")
attributes(fid["RP2"])["altname"] = "ToroTest5"

RP2 =  TestCase("RP2-0.0",  0.0, HydroStatus(5.99924, 19.5975, 460.894, gamma), HydroStatus(5.99242, -6.19633, 46.0950, gamma))
calculate_and_plot(RP2, fid["RP2/0.0"])
RP2 =  TestCase("RP2-0.015",  0.015, HydroStatus(5.99924, 19.5975, 460.894, gamma), HydroStatus(5.99242, -6.19633, 46.0950, gamma))
calculate_and_plot(RP2, fid["RP2/0.015"])
RP2 =  TestCase("RP2-0.025",  0.025, HydroStatus(5.99924, 19.5975, 460.894, gamma), HydroStatus(5.99242, -6.19633, 46.0950, gamma))
calculate_and_plot(RP2, fid["RP2/0.025"])
RP2 =  TestCase("RP2-0.035",  0.035, HydroStatus(5.99924, 19.5975, 460.894, gamma), HydroStatus(5.99242, -6.19633, 46.0950, gamma))
calculate_and_plot(RP2, fid["RP2/0.035"])

#ADD A LOT MORE HERE

#case 3:
create_group(fid, "ToroTest3")
create_group(fid["ToroTest3"], "0.012")
toro3 =  TestCase("ToroTest3",  0.012, HydroStatus(1., 0., 1000., gamma),HydroStatus(1., 0., 0.1, gamma))
calculate_and_plot(toro3, fid["ToroTest3/0.012"])

#case 4:
create_group(fid, "ToroTest4")
create_group(fid["ToroTest4"], "0.035")
toro4 =  TestCase("ToroTest4",  0.035, HydroStatus(1., 0., 0.01, gamma),HydroStatus(1., 0., 100., gamma))
calculate_and_plot(toro4, fid["ToroTest4/0.035"])

#case 6:
create_group(fid, "VaccumExpansionLeft")
create_group(fid["VaccumExpansionLeft"], "0.75")
vaccum_left =  TestCase("VaccumExpansionLeft",  0.75, HydroStatus(0., 0., 0., gamma), HydroStatus(1., 0., 1., gamma))
calculate_and_plot(vaccum_left, fid["VaccumExpansionLeft/0.75"])

#case 7:
create_group(fid, "VaccumExpansionRight")
create_group(fid["VaccumExpansionRight"], "0.75")
vaccum_right =  TestCase("VaccumExpansionRight",  0.75, HydroStatus(1., 0., 1., gamma), HydroStatus(0., 0., 0., gamma))
calculate_and_plot(vaccum_right, fid["VaccumExpansionRight/0.75"])

#case 8:
create_group(fid, "RCVCR")
create_group(fid["RCVCR"], "0.75")
rcvcr =  TestCase("RCVCR",  0.75, HydroStatus(1., -4., 0.4, gamma), HydroStatus(1., 4., 0.4, gamma))
calculate_and_plot(rcvcr, fid["RCVCR/0.75"])

#case 9:
create_group(fid, "ModifiedSod")
create_group(fid["ModifiedSod"], "0.2")
modified_sod = TestCase("ModifiedSod",  0.2, HydroStatus(1., 0.75, 1., gamma), HydroStatus(0.125, 0., 0.1, gamma))
calculate_and_plot(modified_sod, fid["ModifiedSod/0.2"])

#case 9:
create_group(fid, "StreamCollision")
create_group(fid["StreamCollision"], "0.8")
stream = TestCase("StreamCollision",  0.8, HydroStatus(1., 2., 0.1, gamma), HydroStatus(1., -2., 0.1, gamma))
calculate_and_plot(stream, fid["StreamCollision/0.8"])

#case 10:
create_group(fid, "LeBlanc")
create_group(fid["LeBlanc"], "0.5")
leblanc = TestCase("LeBlanc",  0.5, HydroStatus(1., 0., (2. / 3.)*1.e-1, gamma), HydroStatus(1.e-3, 0., (2. / 3.)*1.e-10, gamma))
calculate_and_plot(leblanc, fid["LeBlanc/0.5"])

#case 11:
create_group(fid, "PeakProblem")
create_group(fid["PeakProblem"], "3.9e-3")
peak = TestCase("PeakProblem",  3.9e-3, HydroStatus(0.1261192, 8.9047029, 782.92899, gamma), HydroStatus(6.591493, 2.2654207, 3.1544874, gamma))
calculate_and_plot(peak, fid["PeakProblem/3.9e-3"])

#case 12:
create_group(fid, "SlowShock")
create_group(fid["SlowShock"], "2.")
slow = TestCase("SlowShock",  2., HydroStatus(3.857143, -0.810631, 10.33333, gamma), HydroStatus(1.0, -3.44, 1.0, gamma))
calculate_and_plot(slow, fid["SlowShock/2."])

#case 13:
create_group(fid, "StationaryContact")
create_group(fid["StationaryContact"], "0.012")
stationary = TestCase("StationaryContact",  0.012, HydroStatus(1.0, -19.59745, 1.e3, gamma), HydroStatus(1.0, -19.59745, 1.e-2, gamma))
calculate_and_plot(stationary, fid["StationaryContact/0.012"])
