#!/usr/bin/env python
import numpy as np
from psana import DataSource
from itertools import islice
from psana.dgramedit import DgramEdit, AlgDef, DetectorDef
from psana.psexp import TransitionId
import libpressio

rdetector = 'libpressio'
ifname = "/cds/home/r/robertu/test.xtc2"
ds = DataSource(files=ifname)

for r, run in enumerate(ds.runs()):
    max_npeaks = 2048
    det = run.Detector('libpressio')
    row_array = np.zeros((max_npeaks,) , dtype=np.uint16)
    col_array = np.zeros((max_npeaks,) , dtype=np.uint16)
    for i, evt in enumerate(run.events()):
        det.raw.calib(evt)
        print(f"read evt={i}, run={r}")
  
