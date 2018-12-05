module SpikeTrains

export SpikeTrain, draw_uncorrelated_spikes, draw_correlated_spikes, draw_coincident_spikes, length, iterate, convert, vcat, merge, make_exponentialShift, correlation_code, coincidence_code, plot_spike_trains

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

function has_min_distance(times, min_distance)
    keep = ones(Bool, length(times))
    last_time = -Inf
    for (i,t) ∈ enumerate(times)
        if t < last_time
            keep[i] = false
        else
            last_time = t+min_distance
        end
    end
    return keep
end

ensure_min_distance(times, min_distance) = times[has_min_distance(times, min_distance)]
    
"""
    draw_uncorrelated_spikes(trange, rates; sorted=true)
    
Draw uncorrelated spike trains in the interval `trange` with the respective `rates`.
"""
function draw_uncorrelated_spikes(trange, rates)
    duration = trange[2]-trange[1]
    source = Vector{SpikeTrain}(undef, length(rates))
    # draw uncorrelated spike trains with high rate
    for (i,r) ∈ enumerate(rates)
        num_spikes = rand(Poisson(duration*r))
        spikes = rand(num_spikes).*duration.+trange[1]
        sort!(spikes)

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
function draw_correlated_spikes(trange, rates, c, shift = make_exponentialShift(1.0); min_master_distance=-0)
    return if c≈0.0
        draw_uncorrelated_spikes(trange, rates)
    else
        v = 1/c * mean(rates.^2)/sum(rates)
        p = (c/mean(rates.^2)) .* (rates * rates')

        @assert all(0 .<= p .<= 1) "P not all valid probabilities."

        source = draw_uncorrelated_spikes(trange, fill(v, length(rates)))
        target = Vector{SpikeTrain}(undef, length(rates))

        # ensure minimum distance between all source spike trains
        if min_master_distance > 0
            n = sum(length, source)
            times = Vector{Float64}(undef, n)
            sources = Vector{Int}(undef, n)
            i = 1
            for (s,spikes) ∈ enumerate(source)
                inew = i+length(spikes)
                times[i:inew-1] = spikes.times
                sources[i:inew-1] .= s
                i=inew 
            end
            
            # sort times and corresponding sources
            idx=sortperm(times)
            times=times[idx]
            sources=sources[idx]
            
            # determine which spikes to keep ...
            keep=has_min_distance(times, min_master_distance)
            
            # ... and keep only those spikes for the respective source
            for (s,spikes) ∈ enumerate(source)
                source[s] = SpikeTrain(times[(sources.==s) .& keep])
            end
        end

        # draw correlated spike trains
        for i ∈ eachindex(target)
            t = target[i] = SpikeTrain()
            for (k,s) ∈ enumerate(source)
                num_spikes = rand(Binomial(length(s.times), p[i,k]))
                append!(t.times, sample(s.times, num_spikes; replace=false))
            end

            # shift each value randomly
            t.times .= shift.(t.times)
            
            sort!(t.times)
        end
        target
    end
end

"""
    draw_coincident_spikes(trange, num_spikes, p=ones(Float64,1), shift = make_exponentialShift(1.0); sorted=true)

Draws `num_spikes` Poisson spikes in the interval `trange`, each of which appears 
in each spike-train `i` with the respective probability `p[i]`.
"""
function draw_coincident_spikes(trange, num_spikes, p=ones(Float64,1), shift = make_exponentialShift(1.0); sorted=true)
    times = rand(num_spikes).*(trange[2]-trange[1]).+trange[1]
    
    if sorted
        sort!(times)
    end
    
    target = Vector{SpikeTrain}(undef, length(p))
    for i ∈ eachindex(target)
        n = rand(Binomial(num_spikes, p[i]))
        target[i] = SpikeTrain(sample(times, n; replace=false))
    end
    
    return target
end

"""
    correlation_code(trange, stimulus, rates, correlations; kwargs...)

Each slice `stimulus[:,...,:,i]` contains an array indicating the class that the
spike-train with corresponding index represents in the time interval `trange[i:i+1]`.
Each class `c` results in spike-trains with rate `rates[c]` and mutual pairwise
correlation `correlations[c]`.
"""
function correlation_code(trange, stimulus, rates, correlations; kwargs...)
    @assert length(trange)-1 == size(stimulus)[end] "trange must include n+1 steps for stimuli with last dimension n"
    spiketrains = Array{SpikeTrain}(undef, size(stimulus))
    classes = unique(stimulus)

    for c ∈ classes
        for i ∈ Base.OneTo(length(trange)-1)
            idx = findall(s->c==s, view(stimulus, fill(Colon(), ndims(stimulus)-1)...,i))
            spiketrains[idx,i] = draw_correlated_spikes(trange[i:i+1], fill(rates[c], length(idx)), correlations[c]; kwargs...)
        end
    end

    return spiketrains
end

"""
    coincidence_code(trange, stimulus, background_rate, p; kwargs...)
    
Each slice `stimulus[:,...,:,i]` contains a boolean array corresponding to 
whether or not a spike train with corresponding index participates or not
in the coincident spiking during the time interval `trange[i:i+1]`. 
A single spike is drawn within each interval, that appears with probability `p`
in each of the spike  trains for which the corresponding entry in `stimulus` is `true`.
In addition to the signal spike for each interval, uncorrelated poisson spike-trains,
drawn with a given `background_rate`, are added as background noise.
"""
function coincidence_code(trange, stimulus, background_rate, p; kwargs...)
    @assert length(trange)-1 == size(stimulus)[end] "trange must include n+1 steps for stimuli with last dimension n"
    spiketrains = reshape(draw_uncorrelated_spikes((trange[1],trange[end]), fill(background_rate, length(stimulus))), size(stimulus))
    
    for i ∈ Base.OneTo(length(trange)-1)
        idx = findall(view(stimulus, fill(Colon(), ndims(stimulus)-1)...,i))
        nd = ndims(idx)+1
        spiketrains[idx,i] = merge(cat(spiketrains[idx,i], draw_coincident_spikes(trange[i:i+1], 1, fill(p, length(idx)); kwargs...), dims=nd),dims=nd)
    end
    
    return spiketrains
end

"""
    plot_spike_trains(spiketrains::Array{SpikeTrain}, colors=fill(:auto, length(spiketrains)), plt=plot(); kwargs...)
    
Draw an array of `spiketrains` with given `colors` within a single given or new plot `plt`.
"""
function plot_spike_trains(spiketrains::Array{SpikeTrain}, colors=fill(:auto, length(spiketrains)), plt=plot(); kwargs...)
    for (i,(spiketrain, color)) ∈ enumerate(zip(spiketrains, colors))
        plot!(plt, ([spiketrain.times spiketrain.times fill(NaN, length(spiketrain.times))])'[:], i .- 0.5 .+ ([zeros(length(spiketrain.times)) ones(length(spiketrain.times)) fill(NaN, length(spiketrain.times))])'[:], lc=color; kwargs...)
    end
    return plt
end

end
