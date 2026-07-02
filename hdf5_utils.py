"""Generic nested-dict <-> HDF5 serialization, replacing the JSON outputs.

Nested dicts become HDF5 groups; scalars, strings, lists, and numpy arrays
become datasets. Mirrors the shape produced by the alpha calculation and
secondary-quantity fit functions without hardcoding their keys.
"""
import h5py
import numpy as np


def save_dict_to_hdf5(path, data):
    with h5py.File(path, 'w') as f:
        _write_group(f, data)


def _write_group(group, data):
    for key, value in data.items():
        key = str(key)
        if isinstance(value, dict):
            _write_group(group.create_group(key), value)
        elif isinstance(value, (list, tuple)) and any(isinstance(v, dict) for v in value):
            # e.g. alphas[galaxy_type][region] = one result dict per z-bin
            list_group = group.create_group(key)
            list_group.attrs['is_list'] = True
            for i, item in enumerate(value):
                _write_group(list_group.create_group(str(i)), item)
        elif value is None:
            group.create_dataset(key, data=h5py.Empty('f'))
        else:
            group.create_dataset(key, data=value)


def load_dict_from_hdf5(path):
    with h5py.File(path, 'r') as f:
        return _read_group(f)


def _read_group(group):
    if group.attrs.get('is_list', False):
        return [_read_group(group[str(i)]) for i in range(len(group))]
    data = {}
    for key, item in group.items():
        if isinstance(item, h5py.Group):
            data[key] = _read_group(item)
        elif item.shape is None:
            data[key] = None
        else:
            value = item[()]
            if isinstance(value, bytes):
                value = value.decode('utf-8')
            data[key] = value
    return data
