using SpikeTrains, Plots

rates = [1.0, 2.0, 3.0]
durations = [0.0, 50.0, 100.0]
correlations = [0.0, 0.25, 0.5]
stimulus = cat([1 2 2; 1 2 3], [3 3 3; 3 3 3], dims=3)

spikes1 = draw_uncorrelated_spikes((0,100), rates)

spikes2 = draw_correlated_spikes((0,100), rates, 0.5)

spikes3 = draw_correlated_spikes((0,100), [1.0, 1.0], 0.9)

spikes4 = correlation_code(durations, stimulus, rates, correlations)

spikes5 = merge(spikes4, dims=(3,))


p1=plot_spike_trains(spikes1)
p2=plot_spike_trains(spikes2)
p3=plot_spike_trains(spikes3)
p4=plot_spike_trains(spikes5)

pp = plot(p1,p2,p3,p4, layout=grid(4,1))
savefig(pp, "test.svg")
