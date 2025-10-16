import numpy as np
from typing import Dict, Any

# Global variables for worker processes
_DECODER = None
_LOGICALS = None

def worker_init(decoder, logicals_matrix: np.ndarray):
    """
    Initializer for each worker process in the pool.
    Sets up global variables for the decoder and logicals matrix.
    """
    global _DECODER, _LOGICALS
    _DECODER = decoder
    _LOGICALS = logicals_matrix

def count_errors(det_chunk: np.ndarray, obs_chunk: np.ndarray) -> int:
    """
    Worker function to decode a chunk of detection events and count logical errors.
    """
    if _DECODER is None or _LOGICALS is None:
        raise RuntimeError("Worker process not initialized. Call worker_init.")

    if det_chunk.shape[0] == 0:
        return 0

    # Bit-pack detection events for the compiled decoder API
    det_u8 = det_chunk.astype(np.uint8, copy=False)
    bit_packed_dets = np.packbits(det_u8, axis=1, bitorder="little")
    bit_packed_dets = np.ascontiguousarray(bit_packed_dets)

    # Decode
    pred_obs_b8 = _DECODER.decode_shots_bit_packed(
        bit_packed_detection_event_data=bit_packed_dets
    )
    
    # Unpack predictions
    pred_obs = np.unpackbits(pred_obs_b8, bitorder="little", axis=1)[:, : _LOGICALS.shape[0]]

    # Count errors
    errors = int(np.any(pred_obs != obs_chunk.astype(np.uint8), axis=1).sum())
    
    return errors
