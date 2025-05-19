## Program to run PPSeq on a section of data
# Other Imports
import DelimitedFiles: readdlm
import Random
import StatsBase: quantile
using CSV, DataFrames, Dates
import JSON


# INPUTS:
using ArgParse
s = ArgParseSettings()
@add_arg_table! s begin
    "--data-directory"
        default = "./data/preparedData/" #where is data saved
    "--data-filename"
        default = "defaultData" #data name (don't add the ".txt")
    "--number-of-sequence-types"
        default = "5" # how many sequence types to try fit 
    "--PPSeq_file"
        default = "./PPSeq.jl/src/PPSeq.jl" #the PPSeq.jl file to use 
    "--num-threads"
        default = "1"
    "--results-directory"
        default = "./data/resultsData/"
    "--results-name"
        default = nothing
    "--slurm-array-task-id"
        default = -1
    "--sacred-directory"
        default = "./data/resultsData/"
    "--sacred-data-file"
        default = nothing
end


parsed_args = parse_args(ARGS, s)

PPSeq_file = parsed_args["PPSeq_file"]

# Import PPSeq
include(PPSeq_file)
include("./save_results.jl")
const seq = PPSeq

Data_filename = string(parsed_args["data-directory"],parsed_args["data-filename"],".txt")
num_of_seq_types = parse(Int64, parsed_args["number-of-sequence-types"])
num_threads = parse(Int64,parsed_args["num-threads"])
if parsed_args["results-name"] == nothing
    parsed_args["results-name"] = parsed_args["data-filename"]
end 
#get number of neurons and maxtime from the data JSON file 
json_filename = string(parsed_args["data-directory"],"params_" , parsed_args["data-filename"],".json")
json_params = JSON.parsefile(json_filename)
num_neurons = json_params["number_of_neurons"]
avg_firing_rate = json_params["average_firing_rate"]

timeslicelist = json_params["time_span"]
global max_time = 0
for timeslice in timeslicelist
    	global max_time
  	max_time = max_time + timeslice[2] - timeslice[1]
end
max_time = float(max_time)


# Load spikes.
spikes = seq.Spike[]
file_name = Data_filename
for (n, t) in eachrow(readdlm(file_name, '\t', Float64, '\n'))
    push!(spikes, seq.Spike(Int(n), t))
end

# If there are some sacred sequences extract the relevant parameters here
if parsed_args["sacred-data-file"] != nothing
    print("Using a sacred data file to fix some sequences")
    sacred_data_file = parsed_args["sacred-directory"]*parsed_args["sacred-data-file"]*"/neuron_response.csv"

    # Extract these parameters from the data files
    neuron_responses = CSV.read(sacred_data_file, DataFrame)

    sacred_sequences = ncol(neuron_responses)÷3
    print(typeof(sacred_sequences))

    if mod(nrow(neuron_responses), num_neurons) != 0
        error("Sacred Sequences have different number of neurons to input data!")
    end

    sacred_neuron_responses = neuron_responses[nrow(neuron_responses)-num_neurons+1:nrow(neuron_responses),:]
else
    sacred_sequences = 0
end

# For our grid search, set up parameters to sweep and which values for this slurm_id
slurm_id = parse(Int64, parsed_args["slurm-array-task-id"])
# Define parameter lists
num_of_seqs = [6]
seq_event_factor = [1]
seq_event_conc_param = [1]
event_amp_fudge = [0.3]
event_amp_std_multiplier = [10]
neu_resp_conc = [0.6]
neuron_width_prior = [0.5]
neuron_with_scale = [1.0]
neur_offset_scale_ks = [0.5]
mean_bckgrnd_rate_multiplier = [0.3]
var_bckgrnd_rate_multiplier = [0.3]
bkgd_spikes_conc = [0.3]
max_sequence_length = [60]

# Generate all combinations
base_combinations = collect(Base.Iterators.product(
    num_of_seqs,
    seq_event_factor,
    seq_event_conc_param,
    event_amp_fudge,
    event_amp_std_multiplier,
    neu_resp_conc,
    neuron_width_prior,
    neuron_with_scale,
    neur_offset_scale_ks,
    mean_bckgrnd_rate_multiplier,
    var_bckgrnd_rate_multiplier,
    bkgd_spikes_conc,
    max_sequence_length
))

# Repeat each combination, convert only the tail to Float64
param_grid = []
for combo in base_combinations
    first = combo[1]                      # Int
    rest = map(x -> Float64(x), combo[2:end])  # Float64
    new_combo = (first, rest...)          # Combine
    append!(param_grid, fill(new_combo, n_repeats))
end

config = Dict(
    # Model hyperparameters
    :num_sequence_types         => param_grid[slurm_id + 1][1],
    :seq_type_conc_param        => param_grid[slurm_id + 1][3],
    :seq_event_rate             => 0.2 * param_grid[slurm_id + 1][2] * param_grid[slurm_id + 1][1],


    :mean_event_amplitude       => num_neurons * avg_firing_rate * param_grid[slurm_id + 1][4],
    :var_event_amplitude        => (num_neurons * avg_firing_rate * param_grid[slurm_id + 1][4]) *
                                    param_grid[slurm_id + 1][5],

    :neuron_response_conc_param => param_grid[slurm_id + 1][6],
    :neuron_offset_pseudo_obs   => param_grid[slurm_id + 1][9],
    :neuron_width_pseudo_obs    => param_grid[slurm_id + 1][8],
    :neuron_width_prior         => param_grid[slurm_id + 1][7],

    :num_warp_values            => 1,
    :max_warp                   => 30.0,
    :warp_variance              => 1.0,
    :warp_type                  => 1,

    :mean_bkgd_spike_rate       => num_neurons * avg_firing_rate * param_grid[slurm_id + 1][10],
    :var_bkgd_spike_rate        => (num_neurons * avg_firing_rate * param_grid[slurm_id + 1][11])^2,
    :bkgd_spikes_conc_param     => param_grid[slurm_id + 1][12],
    :max_sequence_length        => param_grid[slurm_id + 1][13],  

    # MCMC Sampling parameters.
    :num_threads => num_threads,
    :num_anneals => 20, # Never set to 1! Makes it break (presumably 0 would as well?)
    :samples_per_anneal => 100,
    :max_temperature => 40.0,
    :samples_after_anneal => 1000,
    :split_merge_moves_during_anneal => 10,
    :split_merge_moves_after_anneal => 10,
    :split_merge_window => 1.0,

    :save_every_after_anneal => 10,
    :save_every_during_anneal => 10,

    # Masking specific parameters
    :are_we_masking => 1, # Binary var, 1 = masking
    :mask_lengths => 5, # In seconds
    :percent_masked => 10, # percentage = number between 0 and 100 (not 0 and 1)
    :num_spike_resamples_per_anneal => 20, # How many resamplings of masked spikes per anneal
    :num_spike_resamples => 100, # How many times to resample masked spikes after annealing
    :samples_per_resample => 10, # How many times the unmasked spikes are sampled for each sampling

    # For training on some data and resampling on another set
    :sacred_sequences => sacred_sequences
);


# 2022 update, should now work!
if config[:are_we_masking] == 1
    println("WARNING: MASKING IS ON, YOU SURE YOU WANT TO MASK?")
end

# Then train the PPSeq Model
# Initialize all spikes to background process.
init_assignments = fill(-1, length(spikes))

# Construct model struct (PPSeq instance).
model = seq.construct_model(config, max_time, num_neurons)

# Run Gibbs sampling with an initial annealing period.
if config[:are_we_masking] == 1
    # First create a random set of masks

    config[:save_every_after_anneal] = min(config[:samples_per_resample], config[:save_every_after_anneal])
    config[:save_every_during_anneal] = min(config[:samples_per_resample], config[:save_every_during_anneal])

    masks = seq.create_random_mask(
        num_neurons,
        max_time,
        config[:mask_lengths] + 0.000001,
        config[:percent_masked]
    )

    masks = seq.merge_contiguous_masks(masks, num_neurons)

    # Then run the easy sampler using tßhese masks
    results = seq.easy_sample_masked!(model, spikes, masks, init_assignments, config)
    results_directory = parsed_args["results-directory"]*parsed_args["results-name"]*"_"*Dates.format(Dates.now(),"ddmmyyy_HHMM")*"/"
	masked_spikes, unmasked_spikes = seq.split_spikes_by_mask(spikes, masks)
	@time save_results_masked(results, config, results_directory, masked_spikes, unmasked_spikes)
    mkdir(results_directory*"trainingData/")
    cp(parsed_args["data-directory"]*parsed_args["data-filename"]*".txt", results_directory*"trainingData/"*parsed_args["data-filename"]*".txt")
    cp(parsed_args["data-directory"]*"params_"*parsed_args["data-filename"]*".json", results_directory*"trainingData/"*"params_"*parsed_args["data-filename"]*".json")

    # Save the results
    #@time save_results_masked(results, config, "./data/resultsData/n_seq"*string(num_of_seq_types)*"_"*Dates.format(Dates.now(),"HHMM"))
    # save_results_masked(results, config, "results/masked_n_seq"*string(num_of_seq_types))
else
    # Check if there are any saved sequences that shouldn't be sampled
    if sacred_sequences != 0
        # Then make a new model that has these fixed parameters
        print("Sactifying model")
        model = seq.sanctify_model(model, Matrix(sacred_neuron_responses), config)
    end

    @time results = seq.easy_sample!(model, spikes, init_assignments, config);

    # Save the results
    # save_results(results, config, "../Simple"*Dates.format(Dates.now(),"HHMM"))
    results_directory = parsed_args["results-directory"]*parsed_args["results-name"]*"_"*Dates.format(Dates.now(),"ddmmyyy_HHMM")*"/"
    save_results(results, config, results_directory)
    mkdir(results_directory*"trainingData/")
    cp(parsed_args["data-directory"]*parsed_args["data-filename"]*".txt", results_directory*"trainingData/"*parsed_args["data-filename"]*".txt")
    cp(parsed_args["data-directory"]*"params_"*parsed_args["data-filename"]*".json", results_directory*"trainingData/"*"params_"*parsed_args["data-filename"]*".json")
end
