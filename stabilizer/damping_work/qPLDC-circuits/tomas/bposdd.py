import pathlib
import numpy as np
from scipy.sparse import csc_matrix
from sinter import CompiledDecoder, Decoder
from stim import DetectorErrorModel

from beliefmatching import detector_error_model_to_check_matrices
from ldpc import BpOsdDecoder
'''This file is copied from Josh Viszlai's gbstim repo'''
class CompiledBPOSD(CompiledDecoder):

    def __init__(self, check_matrices, decoder):
        self.check_matrices = check_matrices
        self.decoder = decoder

    def decode_shots_bit_packed(self, 
                                bit_packed_detection_event_data: np.ndarray
                               ) -> np.ndarray:
        obs_flip_data = []
        for shot_data in bit_packed_detection_event_data:
            unpacked_data = np.unpackbits(shot_data, bitorder='little', count=self.check_matrices.check_matrix.shape[0])
            pred_errors = self.decoder.decode(unpacked_data)
            obs_pred = (self.check_matrices.observables_matrix @ pred_errors) % 2
            obs_flip_data.append(np.packbits(obs_pred, bitorder='little'))

        return np.array(obs_flip_data)

class BPOSD(Decoder):

    def __init__(self, **kwargs):
        self.decoder_kwargs = kwargs

    def compile_decoder_for_dem(self, 
                                dem: DetectorErrorModel
                               ) -> CompiledDecoder:
        check_matrices = detector_error_model_to_check_matrices(dem, allow_undecomposed_hyperedges=True)
        decoder = BpOsdDecoder(check_matrices.check_matrix, channel_probs=check_matrices.priors, **self.decoder_kwargs)
        return CompiledBPOSD(check_matrices, decoder)

    def decode_via_files(self, 
                         *, 
                         num_shots: int, 
                         num_dets: int, 
                         num_obs: int, 
                         dem_path: pathlib.Path, 
                         dets_b8_in_path: pathlib.Path, 
                         obs_predictions_b8_out_path: pathlib.Path, 
                         tmp_dir: pathlib.Path
                        ) -> None:
        raise NotImplementedError()

class CompiledCSSWindowBPOSD(CompiledDecoder):
    
    def __init__(self, num_dets, num_obs, z_decoder_schedule, x_decoder_schedule):
        self.num_dets = num_dets
        self.num_obs = num_obs
        self.z_decoder_schedule = z_decoder_schedule
        self.x_decoder_schedule = x_decoder_schedule
    
    def decode_shots_bit_packed(self, bit_packed_detection_event_data: np.ndarray) -> np.ndarray:
        obs_flip_data = []
        for shot_data in bit_packed_detection_event_data:
            obs_pred = np.zeros(self.num_obs, dtype=np.uint8)
            unpacked_data = np.unpackbits(shot_data, bitorder='little', count=self.num_dets)
            for t, (z_decoder, det_map, error_map, commit_until) in self.z_decoder_schedule.items():
                window = np.array([unpacked_data[det] for det, _ in det_map.items()])
                pred_errors = z_decoder.decode(window)
                for i, err in enumerate(pred_errors):
                    t, dets, obs = error_map[i]
                    if t < commit_until:
                        unpacked_data[dets] ^= err
                        obs_pred[obs] ^= err
            for t, (x_decoder, det_map, error_map, commit_until) in self.x_decoder_schedule.items():
                window = np.array([unpacked_data[det] for det, _ in det_map.items()])
                pred_errors = x_decoder.decode(window)
                for i, err in enumerate(pred_errors):
                    t, dets, obs = error_map[err]
                    if t < commit_until:
                        unpacked_data[dets] ^= err
                        obs_pred[obs] ^= err
            obs_flip_data.append(np.packbits(obs_pred, bitorder='little'))
        
        return np.array(obs_flip_data)
            

class CSSWindowBPOSD(Decoder):

    def __init__(self, window_size, stride, **kwargs):
        self.window_size = window_size
        self.stride = stride
        self.decoder_kwargs = kwargs
    
    def dicts_to_matrix(self, dicts, t0):
        det_list = []
        det_map = {}
        error_map = {}
        priors = []
        row_idx = []
        col_idx = []
        i = 0
        for j, d in enumerate(dicts):
            for dets, (obs, prior) in d.items():
                error_map[i] = (t0 + j, list(dets), obs)
                priors += [prior]
                for det in dets:
                    if not det in det_list:
                        det_map[det] = len(det_list)
                        det_list += [det]
                    row_idx += [det_map[det]]
                    col_idx += [i]
                i += 1
        
        return csc_matrix((np.ones(len(row_idx)), (row_idx, col_idx)), shape=(len(det_map.keys()), len(error_map.keys()))), priors, det_map, error_map

    def sparse_dem_to_checks(self, dem, w, stride):
        x_errors = {} 
        z_errors = {}
        det_coords = dem.get_detector_coordinates()
        for instr in dem.flattened():
            if instr.type == 'error':
                x_dets = []
                z_dets = []
                obs = []
                for targ in instr.targets_copy():
                    if targ.is_relative_detector_id():
                        if det_coords[targ.val][0] % 2:
                            x_dets += [targ.val]
                        else:
                            z_dets += [targ.val]
                    elif targ.is_logical_observable_id():
                        obs += [targ.val]
                prior = instr.args_copy()[0]
                if len(z_dets):
                    z_t = int(min(map(lambda det_id: det_coords[det_id][-1], z_dets)))
                    z_sig = tuple(sorted(z_dets))
                    obs_sig = sorted(obs)
                    if not z_t in x_errors:
                        x_errors[z_t] = {}
                    if z_sig in x_errors[z_t]:
                        x_errors[z_t][z_sig][1] += prior # Accumulate prior
                    else:
                        x_errors[z_t][z_sig] = [obs_sig, prior] # Assume Z obs tracking
                if len(x_dets):
                    x_t = int(min(map(lambda det_id: det_coords[det_id][-1], x_dets)))
                    x_sig = tuple(sorted(x_dets))
                    if not x_t in z_errors:
                        z_errors[x_t] = {}
                    if x_sig in z_errors[x_t]:
                        z_errors[x_t][x_sig][1] += prior # Accumulate prior
                    else:
                        z_errors[x_t][x_sig] = [[], prior]
        z_window_checks = {}
        x_window_checks = {}
        z_ub = max(x_errors.keys())
        x_ub = max(z_errors.keys())
        for t in range(min(x_errors.keys()), z_ub + 1, stride):
            w_size = w if t + w <= z_ub + 1 else z_ub + 1 - t
            z_window_checks[t] = self.dicts_to_matrix([x_errors[t + i] for i in range(w_size)], t)
        for t in range(min(z_errors.keys()), x_ub + 1, stride):
            w_size = w if t + w <= x_ub + 1 else x_ub + 1 - t
            x_window_checks[t] = self.dicts_to_matrix([z_errors[t + i] for i in range(w_size)], t)
        return z_window_checks, x_window_checks
                
                    
    def compile_decoder_for_dem(self, dem: DetectorErrorModel) -> CompiledDecoder:
        z_window_checks, x_window_checks = self.sparse_dem_to_checks(dem, self.window_size, self.stride)
        z_decoder_schedule = {}
        x_decoder_schedule = {}
        for t, (check_mat, priors, det_map, error_map) in z_window_checks.items():
            decoder = BpOsdDecoder(check_mat, channel_probs=priors, **self.decoder_kwargs)
            z_decoder_schedule[t] = (decoder, det_map, error_map, t + self.stride)
        for t, (check_mat, priors, det_map, error_map) in x_window_checks.items():
            decoder = BpOsdDecoder(check_mat, channel_probs=priors, **self.decoder_kwargs)
            x_decoder_schedule[t] = (decoder, det_map, error_map, t + self.stride)
        return CompiledCSSWindowBPOSD(dem.num_detectors, dem.num_observables, z_decoder_schedule, x_decoder_schedule)
        

    def decode_via_files(self, 
                         *, 
                         num_shots: int, 
                         num_dets: int, 
                         num_obs: int, 
                         dem_path: pathlib.Path, 
                         dets_b8_in_path: pathlib.Path, 
                         obs_predictions_b8_out_path: pathlib.Path, 
                         tmp_dir: pathlib.Path
                        ) -> None:
        raise NotImplementedError()