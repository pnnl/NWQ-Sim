from CrossLS import CrossLS
import sinter
from typing import Tuple, List
import pandas as pd
from stimbposd import sinter_decoders
from utils import save_dataframe, convert_to_per_round
from datetime import datetime


# 

# %% [markdown]
# ### 1. Regular BPOSD

# %%
if __name__ == "__main__":
    tasks = []
    # All experiments
    PQRM_para_list = [(1,2,4), (1,3,5), (1,4,6)]
    noise_list = [1e-3,5e-4]
    for noise in noise_list:
        noise_profile = [1e-6, noise, noise, noise]
        for PQRM_para in PQRM_para_list:
            if PQRM_para == (1,2,4):
                d_surf_list = [3,5,7]
            if PQRM_para in [(1,3,5), (1,4,6)]:
                d_surf_list = [7]
            for d_surf in d_surf_list:
                rounds = d_surf
                for PQRM_state in ["Z", "X"]:
                    circuit = CrossLS(PQRM_para, d_surf, rounds, PQRM_state, noise_profile, punctured=True)
                    tasks.append(sinter.Task(circuit=circuit, 
                        json_metadata={'d_surf': d_surf, 'PQRM_para': PQRM_para, 'rounds': rounds, 'PQRM_state': PQRM_state, 
                                        'p1': noise_profile[0], 'p2': noise_profile[1], 'pM': noise_profile[2], 'pR': noise_profile[3]}))

    collected_stats: List[sinter.TaskStats] = sinter.collect(
        num_workers=6, # change if needed
        tasks=tasks,
        decoders=['bposd'],
        max_shots=100, # change to large number for full scale experiment
        max_errors=10, # Whichever first hits the threshold stops the simulation
        custom_decoders=sinter_decoders(),
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

    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    filename = f"surface_sim_{timestamp}.csv"
    save_dataframe(df_sf, "data/results", filename)