module SpikeTrains

export SpikeTrain, draw_uncorrelated_spikes, draw_correlated_spikes, length, iterate, convert, vcat, merge, make_exponentialShift, correlation_code, plot_spike_trains

using Distributions, Plots#, PlotRecipes


struct SpikeTrain
    times::Vector{Float64}
    SpikeTrain(times=Float64[]) = new(sort(times))
end

Base.length(s::SpikeTrain) = length(s.times)
Base.iterate(s::SpikeTrain, args...) = iterate(s.times, args...)
Base.convert(SpikeTrain, v) = SpikeTrain(v)
Base.vcat(s::SpikeTrain...) = SpikeTrain(vcat(getfield.(s,:times)...))
Base.merge(s::Array{SpikeTrain}; dims=Base.OneTo(ndims(s))) = dropdims(mapslices(x->vcat(x...), s, dims=dims), dims=dims)


function draw_uncorrelated_spikes(trange, rates; sorted=true)
    duration = trange[2]-trange[1]
    source = Vector{SpikeTrain}(undef, length(rates))
    # draw uncorrelated spike trains with high rate
    for (i,r) ∈ enumerate(rates)
        num_spikes = rand(Poisson(duration*r))
        spikes = rand(num_spikes).*duration.+trange[1]
        if sorted
            sort!(spikes)
        end
        source[i] = SpikeTrain(spikes)
    end

    return source
end

make_exponentialShift(τ) = t -> t+rand(Exponential(τ))

"""
    draw_correlated_spikes(trange, rates, c, τ=1.0)

Draws spike trains with given `rates` within the given time interval `trange` and
pairwise correlation coefficients c. Exponential noise is added to result in an
exponential cross correlation with time constant `τ`. The algorithm implemented
is the offline version of the mixture process (4.6.1) by Brette 2008 [1].

[1](http://romainbrette.fr/WordPress3/wp-content/uploads/2014/06/Brette2008NC.pdf)
"""
function draw_correlated_spikes(trange, rates, c, shift = make_exponentialShift(1.0); sorted=true)
    if c≈0.0
        return draw_uncorrelated_spikes(trange, rates; sorted=sorted)
    end

    v = 1/c * mean(rates.^2)/sum(rates)
    p = (c/mean(rates.^2)) .* (rates * rates')

    @assert all(0 .<= p .<= 1) "P not all valid probabilities."

    source = draw_uncorrelated_spikes(trange, fill(v, length(rates)), sorted=false)
    target = Vector{SpikeTrain}(undef, length(rates))

    # draw correlated spike trains
    for i ∈ eachindex(target)
        t = target[i] = SpikeTrain()
        for (k,s) ∈ enumerate(source)
            num_spikes = rand(Binomial(length(s.times), p[i,k]))
            append!(t.times, sample(s.times, num_spikes; replace=false))
        end

        # shift each value randomly
        t.times .= shift.(t.times)

        if sorted
            sort!(t.times)
        end
    end

    return target
end


function plot_spike_trains(spiketrains::Array{SpikeTrain}, colors=fill(:auto, length(spiketrains)), plt=plot())
    for (i,(spiketrain, color)) ∈ enumerate(zip(spiketrains, colors))
        plot!(plt, ([spiketrain.times spiketrain.times fill(NaN, length(spiketrain.times))])'[:], i .- 0.5 .+ ([zeros(length(spiketrain.times)) ones(length(spiketrain.times)) fill(NaN, length(spiketrain.times))])'[:], lc=color)
    end
    return plt
end

function correlation_code(trange, stimulus, rates, correlations)
    @assert length(trange)-1 == size(stimulus)[end] "trange must include n+1 steps for stimuli with last dimension n"
    spiketrains = Array{SpikeTrain}(undef, size(stimulus))
    classes = unique(stimulus)

    for c ∈ classes
        for i ∈ Base.OneTo(length(trange)-1)
            idx = findall(s->c==s, view(stimulus, fill(Colon(), ndims(stimulus)-1)...,i))
            spiketrains[idx,i] = draw_correlated_spikes(trange[i:i+1], fill(rates[c], length(idx)), correlations[c])
        end
    end

    return spiketrains
end

end
