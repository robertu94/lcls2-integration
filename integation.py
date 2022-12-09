#!/usr/bin/env python
import numpy as np
import random
import libpressio
from pprint import pprint

MAX_PEAKS = 4
MAX_X, MAX_Y = 1024, 1024

n_peaks = np.uint16(MAX_PEAKS)
rows = np.zeros((n_peaks,), np.uint16)
cols = np.zeros((n_peaks,), np.uint16)
peaks = np.zeros((n_peaks, 3), np.uint16)

# make synthetic data
for i in range(n_peaks):
    rows[i] = random.randint(0, MAX_X - 1)
    cols[i] = random.randint(0, MAX_Y - 1)
    peaks[i, 0] = rows[i]
    peaks[i, 1] = cols[i]
    peaks[i, 2] = 0 # event_idx

data = np.random.random((1, MAX_X, MAX_Y)).astype(np.float32)
print(data)
print(peaks)


#### formulate compression

comp = libpressio.PressioCompressor.from_config({
        "compressor_id": "pressio",
        "early_config": {
            "pressio": {
                "pressio:compressor": "roibin",
                "roibin": {
                    "roibin:metric": "composite",
                    "roibin:background": "binning",
                    "roibin:roi": "fpzip",
                    "background": {
                        "binning:compressor": "pressio",
                        "pressio": {
                            "pressio:compressor": "sz3"
                        }
                    },
                    "composite": {
                        "composite:plugins": ["size", "time"]
                    }
                }
            }
        },
        "compressor_config": {
            "pressio": {
                "roibin": {
                    "roibin:roi_size": [8, 8, 0],
                    "roibin:centers": peaks,
                    "roibin:nthreads": 1,
                    "roi": {
                        "fpzip:prec": 0
                    },
                    "background": {
                        "binning:shape": [2, 2, 1],
                        "binning:nthreads": 4,
                        "pressio": {
                            "pressio:abs": 90.0
                        }
                    },
                }
            }
        },
        "name": "pressio"
    })
pprint(comp.get_config())

output = data.copy()
compressed = comp.encode(data)
comp.decode(compressed, output)
pprint(comp.get_metrics())

