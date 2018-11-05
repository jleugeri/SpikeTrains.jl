using SpikeTrains, Plots

rates = [1.0, 2.0, 3.0]
durations = [0.0, 50.0, 100.0]
durations2 = [0.0, 25.0, 50.0, 75.0, 100.0]
correlations = [0.0, 0.25, 0.5]
stimulus = cat([1 2 2; 1 2 3], [3 3 3; 3 3 3], dims=3)
stimulus2 = cat([true true false; true false false], [false false true; false true true],[false true true; true true false], [true true true; true true true], dims=3)

spikes1 = draw_uncorrelated_spikes((0,100), rates)

spikes2 = draw_correlated_spikes((0,100), rates, 0.5)

spikes3 = draw_correlated_spikes((0,100), [1.0, 1.0], 0.9)

spikes4 = correlation_code(durations, stimulus, rates, correlations)

spikes5 = merge(spikes4, dims=(3,))

spikes6 = draw_coincident_spikes((0,100), 1, ones(Int64, 5))

spikes7 = coincidence_code(durations2, stimulus2, 0.05, 1.0)

spikes8 = merge(spikes7, dims=(3,))

p1=plot_spike_trains(spikes1;legend=false)
p2=plot_spike_trains(spikes2;legend=false)
p3=plot_spike_trains(spikes3;legend=false)
p4=plot_spike_trains(spikes5;legend=false)
p5=plot_spike_trains(spikes6;legend=false, xlims=(0,100))
p6=plot_spike_trains(spikes8;legend=false)

pp = plot(p1,p2,p3,p4,p5,p6, layout=grid(6,1), size=(500,1000))
savefig(pp, "test.svg")
