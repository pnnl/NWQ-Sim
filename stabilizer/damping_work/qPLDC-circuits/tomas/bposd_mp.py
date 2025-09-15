import numpy as np

_DECODER = None
_LOGICALS = None

def worker_init(check_matrix, priors, logicals_matrix, decoder_kwargs):
    global _DECODER, _LOGICALS
    from ldpc import bposd_decoder
    _DECODER = bposd_decoder(
        check_matrix,
        channel_probs=priors,
        **decoder_kwargs
    )
    _LOGICALS = logicals_matrix

def count_errors(det_chunk, obs_chunk):
    errors = 0
    for shot, obs in zip(det_chunk, obs_chunk):
        phys = _DECODER.decode(shot)
        log_pred = (_LOGICALS @ phys) % 2
        if np.any(log_pred != obs):
            errors += 1
    return errors