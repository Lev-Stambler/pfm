[Hx, Hz, N_Qubits, N_Stabilizers] = get_hypergraph('tmp/hypergraph.json');
Hz_t = transpose(Hz);
H_transpose_z_projection = Hz_t * inv(Hz * Hz_t) * Hz;
format shortg

[p, n_trials, cutoff, file_save_block_error] = get_params('tmp/params.json');

simulate_n_rounds(file_save_block_error, p, n_trials, cutoff, N_Qubits, Hx, Hz, N_Stabilizers, H_transpose_z_projection)

% Compare with success rate of non problalistic
function success_probability = simulate_n_rounds(file_save_block_err, p, n, cutoff, n_qubits, Hx, Hz, n_stabilizers, H_transpose_z_projection)
 numb_successes = 0;
 for r = 1:n
   fprintf("\nRunning round %d\n\n", r);
   numb_successes = numb_successes + simulate_round_filter(p, cutoff, n_qubits, Hx, Hz, n_stabilizers, H_transpose_z_projection);
 end
 block_error_rate = 1 - numb_successes / n;
 fprintf("\n Block error rate: %d\n", block_error_rate);
 fileID = fopen(file_save_block_err, 'w');
 fwrite(fileID, string(block_error_rate));
 fclose(fileID);
end

% Currently only runs one loop iteration of the small set flip error
function success = simulate_round_filter(p, cutoff, n_qubits, Hx, Hz, n_stabilizers, H_transpose_z_projection)
  error = get_random_error(p, n_qubits);
  current_correction = sparse(1:n_qubits, repelem(1, n_qubits), repelem(0, n_qubits), n_qubits, 1);
  syndrome = matmul(Hx, error);
  current_syndrome = syndrome;
  success = 1;
  
  n_checked_over = 0;
  n_unchecked_over = 0;
  
  % this is used for testing if only including an odd number of
  % flips changes the outcome
  % Setting to 1 will check both even and odd number of flips.
  % 2 will check only odd
  skip_by = 1;
  
  % Run until the syndrome is 0 and thus the error is corrected
  while nnz(current_syndrome) > 0
      fprintf('Number of non zero syndromes: %d\n', nnz(current_syndrome));
      best_correction = repelem(0, n_qubits);
      best_syndrome = repelem(0, n_stabilizers);
      
      % Get only the qubits which could have maybe errored
      qubit_counts_for_syndrome = transpose(Hx) * current_syndrome; 
      numb_overlapping_potential_qubits = Hz * qubit_counts_for_syndrome;
      [sorted_overlap, generators_sorted_incr_idxs] = sort(numb_overlapping_potential_qubits);
      
      % Find the "important" qubits
      generators_sorted = flip(generators_sorted_incr_idxs);
      
      % Keeps track of the maximum value to be optimized (i.e. the decrease
      % in syndrome size per qubit flipped
      target = -Inf;
      % Only range over the top 30 "important" stabilizers
      if cutoff == -1
        to_range_over = 1:n_stabilizers;
      else
        to_range_over = 1:cutoff;
      end
      for overall_idx = to_range_over
        r = generators_sorted(overall_idx);
        [a, generators] = find(Hz(r, :));
        
        % Skip a generator which cannot possibly decrease the syndrome
        if numb_overlapping_potential_qubits(r) < 1 
          %disp("NO OVERLAP");
          n_unchecked_over = n_unchecked_over + 1;
          continue;
        end
        
        n_checked_over = n_checked_over + 1;
        
        POWER_SET_MAX_SIZE = 100000;
        % Get the power set of generators
        generators_power = power_set(generators, skip_by, POWER_SET_MAX_SIZE);
        %size(generators_power)
         % Go over every set in the powerset
         for power_set_size_div_skip_by = 1:min(POWER_SET_MAX_SIZE, floor(size(generators, 2) / skip_by))
           power_set_size = power_set_size_div_skip_by * skip_by;
           for i = 1:size(generators_power{power_set_size}, 1)
               gRaw = generators_power{power_set_size}(i, :);
               g = sparse(gRaw(1, :), repelem(1, power_set_size), repelem(1, power_set_size), n_qubits, 1);
               % Get the syndrome
               syndrome_x = matmul(Hx, g);
               v = (nnz(current_syndrome) - nnz(xor(current_syndrome, syndrome_x))) / power_set_size;
               % Update the target if necessary
               if v > target
                %  fprintf("V: %f, idx: %d, T: %d\n", v, overall_idx, numb_overlapping_potential_qubits(r));
                 for qubitInd = generators
                  %  fprintf("With qubit values: %d\n", qubit_counts_for_syndrome(qubitInd));
                 end
                 target  = v;
                 best_correction = g;
                 best_syndrome = syndrome_x;
               end
           end
         end
      end
      % there is no generator g such that its application decreases the
      % number of nnz syndrome entries
      if target <= 0
        success = 0;
        fprintf("Could not correct the error");
        return
      end
      current_syndrome = xor(current_syndrome, best_syndrome);
      current_correction = xor(current_correction, best_correction);
  end
  correct = is_correction_correct(error, current_correction, H_transpose_z_projection);
  if correct
    success = 1;
    fprintf("Found correction with size: %d", nnz(current_correction));
  else
    fprintf("Found correction with size: %d but not in row space", nnz(current_correction));
  end
end

function res = matmul(A, B)
  res = mod(A * B, 2);
end

function random_error = get_random_error(p, n_qubits)
    n =  n_qubits; m = 1;
    random_error = rand(n,m) < p;
end

% Check if the generated correction - error is in the row space of opposing
% H
function correct = is_correction_correct (random_error, correction, opposing_Ht_projection)
  size(opposing_Ht_projection);
  difference = random_error - correction;
  correct = isequal(opposing_Ht_projection * difference, difference);
end

function [p_error, n_trials, cutoff, file_save_block_error] = get_params(json_file)
  fid = fopen(json_file);
  raw = fread(fid, inf);
  str = char(raw');
  fclose(fid);
  jsonData = jsondecode(str);
  p_error = jsonData.p;
  n_trials = jsonData.n_trials;
  cutoff = jsonData.cutoff;
  file_save_block_error = jsonData.file_out;
end

function [Hx, Hz, N_Qubits, N_Stabilizers] = get_hypergraph(json_file)
    fid = fopen(json_file);
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    jsonData = jsondecode(str);
    
    N_Qubits = jsonData.N_Qubits;
    N_Stabilizers = jsonData.N_Stabilizers;
    Hx = sparse(jsonData.Ix, jsonData.Jx, jsonData.Vx, jsonData.N_Stabilizers, jsonData.N_Qubits);
    Hz = sparse(jsonData.Iz, jsonData.Jz, jsonData.Vz, jsonData.N_Stabilizers, jsonData.N_Qubits);
end

function p = power_set(set, skip_by, max_size)
    p = {};
    for nn = 1:min(max_size, floor(size(set, 2) / skip_by))
        % all the combinations taken nn items at a time
        p{nn * skip_by} = combnk(set, nn * skip_by);
    end

end
