# %%
import stim
from planar_PQRM import PlanarSurfaceCode
from CrossLS import CrossLS
import pymatching
import sinter
from typing import Tuple, List
from itertools import product
import copy
import pandas as pd
import numpy as np
from bposdd import BPOSD
# from stimbposd import BPOSD
from stimbposd import SinterDecoder_BPOSD, sinter_decoders


# %%
# Compute the logical error rate (LER) per round
def convert_to_per_round(total_LER: float, rounds: int) -> float:
    per_round_LER = 1/2 * (1 - np.abs(1-2*total_LER)**(1/rounds))
    return per_round_LER

# %% [markdown]
# ### 1. Regular BPOSD

# %%
if __name__ == "__main__":
    tasks = []
    noise = 0.001
    noise_profile = [1e-6, noise, noise, noise]
    PQRM_para_list = [(1,3,5)]
    for PQRM_para in PQRM_para_list:
        if PQRM_para == (1,2,4):
            d_surf_list = [3,5,7]
        if PQRM_para == (1,3,5):
            d_surf_list = [7]
        for d_surf in d_surf_list:
            rounds = d_surf
            for PQRM_state in ["Z"]:
                circuit = CrossLS(PQRM_para, d_surf, rounds, PQRM_state, noise_profile, punctured=True)
                tasks.append(sinter.Task(circuit=circuit, 
                    json_metadata={'d_surf': d_surf, 'PQRM_para': PQRM_para, 'rounds': rounds, 'PQRM_state': PQRM_state, 
                                    'p1': noise_profile[0], 'p2': noise_profile[1], 'pM': noise_profile[2], 'pR': noise_profile[3]}))

    collected_stats: List[sinter.TaskStats] = sinter.collect(
        num_workers=48,
        tasks=tasks,
        decoders=['bposd'],
        max_shots=1000000, # change to large number for full scale experiment
        max_errors=15, # Whichever first hits the threshold stops the simulation
        custom_decoders=sinter_decoders(),
        # custom_decoders={'bposd': BPOSD(max_iter = 1000, bp_method = 'ms', osd_order = 10, osd_method='osd_cs')},
        print_progress=True,
    )

    # Convert collected_stats into Dataframe
    records = []
    for stat in collected_stats:
        record = {
            'decoder': stat.decoder,
            'shots': stat.shots,
            'errors': stat.errors,
            'seconds': stat.seconds,
            'LER': convert_to_per_round(stat.errors / stat.shots, 1) if stat.shots > 0 else None, # don't average
            'rounds': stat.json_metadata['rounds'],
            'd_surf': stat.json_metadata['d_surf'],
            'PQRM_para': stat.json_metadata['PQRM_para'],
            'PQRM_state': stat.json_metadata['PQRM_state'],
            'p1': stat.json_metadata['p1'], 'p2': stat.json_metadata['p2'], 'pM': stat.json_metadata['pM'], 'pR': stat.json_metadata['pR']
        }
        records.append(record)
    df_sf = pd.DataFrame(records)
    df_sf.to_csv("Cross_sim_731.csv", index=False) # If you wanna save the experiment data to csv