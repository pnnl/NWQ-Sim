from CrossLS import CrossLS
from typing import Tuple, List
import pandas as pd
from utils import save_dataframe, convert_to_per_round
from datetime import datetime
from utils import dem_to_check_matrices
import stim
import numpy as np
import multiprocessing as mp
# 


# %% [markdown]
# ### 1. Regular BPOSD

# %%

def gpu_decode_worker(dem, chkDenseForNV, decoder_params, fail_target, print_interval,
 max_shots, sim_num, fail_num, lock, gpu_idx, postselect, m):
    import os
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_idx) 
    import cudaq_qec as qec
    # initialize the decoder
    nvdec = qec.get_decoder('nv-qldpc-decoder', chkDenseForNV, **decoder_params)
    dem_sampler: stim.CompiledDemSampler = dem.compile_sampler(seed=os.getpid())
    batch_size = 10000 

    while True:
        with lock:
            if fail_num.value >= fail_target or sim_num.value >= max_shots:
                break
        current_errors = 0
        det_data, obs_data, _ = dem_sampler.sample(shots=batch_size, return_errors=False, bit_packed=False)
        if postselect:
            postselect_indices = np.all(det_data[:, -m:] == 0, axis=1)
            det_data = det_data[postselect_indices, :]
            obs_data = obs_data[postselect_indices, :]

        results = nvdec.decode_batch(det_data)
        decoded = np.array([e.result for e in results])
        # update the statistics
        for i in range(decoded.shape[0]):
            ans = (obs @ decoded[i] + obs_data[i]) % 2 
            current_errors += ans.any()
        with lock:
            fail_num.value += current_errors
            sim_num.value += batch_size
            post_counter.value += det_data.shape[0] if postselect else batch_size
            if sim_num.value % print_interval ==0:
                ler = fail_num.value/post_counter.value
                error_bar = 1.96 * np.sqrt(ler*(1-ler)/post_counter.value)
                print(f"error: {fail_num.value}, shots: {sim_num.value:,}, LER: {ler:.3e} ± {error_bar:.3e}, "
                 f"selected: {post_counter.value:,}, post-selction rate: {post_counter.value/sim_num.value:.5f}")


if __name__ == "__main__":
    tasks = []
    noise = 0.001
    noise_profile = [1e-6, noise, noise, noise]
    postselect = True
    # All experiments
    PQRM_para_list = [(1,2,4), (1,3,5), (1,4,6)]
    for PQRM_para in PQRM_para_list:
        if PQRM_para == (1,2,4):
            d_surf_list = [3,5,7]
        if PQRM_para in [(1,3,5), (1,4,6)]:
            d_surf_list = [7]
        for d_surf in d_surf_list:
            rounds = d_surf
            for PQRM_state in ["Z", "X"]:
                circuit = CrossLS(PQRM_para, d_surf, rounds, PQRM_state, noise_profile, punctured=True)
                tasks.append((circuit, {'d_surf': d_surf, 'PQRM_para': PQRM_para, 'rounds': rounds, 'PQRM_state': PQRM_state, 
                                    'p1': noise_profile[0], 'p2': noise_profile[1], 'pM': noise_profile[2], 'pR': noise_profile[3],
                                    'postselect': postselect, 'm': PQRM_para[2]}))

    # replace sinter with homemade simulator
    max_shots = 1e8
    max_errors = 50
    print_interval = 100000
    gpu_batch = 1000 # increase this if you have more GPU memory
    gpu_num = 4
    process_per_gpu = 2 # number of processes per GPU
    collected_stats = []
    for task in tasks:
        print("start simulation task")
        print(task[1])
        metadata = task[1]
        circ = task[0]
        current_shots = 0
        current_errors = 0
        # initialize the sampler
        dem = circ.detector_error_model()
        pcm, obs, priors, _ = dem_to_check_matrices(dem, return_col_dict=True)
        # decoder parameters
        decoder_opts = dict() 
        decoder_opts['error_rate_vec'] = priors
        decoder_opts['max_iterations'] = 1000 # num of bp iterations
        decoder_opts['use_osd'] = True
        decoder_opts['osd_method'] = 3 # using OSD-CS
        decoder_opts['osd_order'] = 10 # OSD depth
        decoder_opts['bp_batch_size'] = gpu_batch
        decoder_opts['use_sparsity']=True
        chkDenseForNV = np.array(pcm.todense(order='C'))
        # initialize the multiprocessing manager
        manager = mp.Manager()
        fail_counter = manager.Value('i', 0)  # Counter for failures
        sim_counter = manager.Value('i', 0)  # Counter for total simulations
        post_counter = manager.Value('i', 0)  # Counter for post-selected shots
        lock = manager.Lock()
        start_time = datetime.now()
        # start the multiprocessing pool
        processes = []
        for i in range(gpu_num*process_per_gpu):
            p_ = mp.Process(target=gpu_decode_worker, args=(dem, chkDenseForNV, decoder_opts, 
            max_errors, print_interval, max_shots, sim_counter, fail_counter, lock, i%gpu_num, postselect, metadata['m']))
            p_.start()
            processes.append(p_)

        # Wait for all processes to finish
        for p_ in processes:
            p_.join()
        end_time = datetime.now()
        current_shots = sim_counter.value
        post_shots = post_counter.value
        current_errors = fail_counter.value
        ler = current_errors/post_shots
        error_bar = 1.96 * np.sqrt(ler*(1-ler)/post_shots) # 95% confidence interval
        print(f"error: {current_errors}, shots: {current_shots:,}, LER: {ler:.3e} ± {error_bar:.3e},"
         f"selected: {post_shots:,}, post-selction rate: {post_shots/current_shots:.5f}")
        
        stat = {
            'decoder': 'nv-qldpc-decoder',
            'shots': current_shots,
            'errors': current_errors,
            'seconds': (end_time - start_time).total_seconds(),
            'json_metadata': metadata,
            'post_shots': post_shots,
        }
        collected_stats.append(stat)
        print("simulation task done\n\n")
    # Convert collected_stats into Dataframe
    records = []
    for stat in collected_stats:
        record = {
            'decoder': stat['decoder'],
            'shots': stat['shots'],
            'errors': stat['errors'],
            'seconds': stat['seconds'],
            'post_shots': stat['post_shots'],
            'LER': convert_to_per_round(stat['errors'] / stat['post_shots'], stat['json_metadata']['rounds']) if stat['shots'] > 0 else None,
            **stat['json_metadata'],
        }
        records.append(record)
    df_sf = pd.DataFrame(records)

    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    filename = f"surface_sim_{timestamp}.csv"
    save_dataframe(df_sf, "data/results", filename)
