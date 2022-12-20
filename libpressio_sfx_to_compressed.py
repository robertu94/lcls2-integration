#!/usr/bin/env python
import numpy as np
from psana import DataSource
from itertools import islice
from psana.dgramedit import DgramEdit, AlgDef, DetectorDef  # type: ignore
from psana.psexp import TransitionId
import libpressio


def make_compressor(n_peaks, rows, cols):
    peaks = np.zeros((n_peaks, 3), np.uint64)
    for i in range(n_peaks):
        peaks[i, 0] = rows[i]
        peaks[i, 1] = cols[i]
        peaks[i, 2] = 0  # event_idx== 0 because 1 event at a time
    comp = libpressio.PressioCompressor.from_config(
        {
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
                            "pressio": {"pressio:compressor": "sz3"},
                        },
                        "composite": {"composite:plugins": ["size", "time"]},
                    },
                }
            },
            "compressor_config": {
                "pressio": {
                    "roibin": {
                        "roibin:roi_size": [8, 8, 0],
                        "roibin:centers": peaks,
                        "roibin:nthreads": 1,
                        "roi": {"fpzip:prec": 0},
                        "background": {
                            "binning:shape": [2, 2, 1],
                            "binning:nthreads": 4,
                            "pressio": {"pressio:abs": 90.0},
                        },
                    }
                }
            },
            "name": "pressio",
        }
    )
    return comp


memsize = 164000000
outbuf = bytearray(memsize)


def save_dgramedit(dg_edit, outbuf, outfile):
    """Save dgram edit to output buffer and write to file"""
    dg_edit.save(outbuf)
    outfile.write(outbuf[: dg_edit.size])


sfx_detector = "sfx"
rdetector = "Rayonix"
odetector = "libpressio"
nodeId = 1
namesId = {
    rdetector: 0,
    odetector: 0,
    "runinfo": 1,
    "scan": 2,
}


ifname = ("/cds/home/s/smarches/git/xtc1to2/examples/converted/"
          "mfxp22820_13b.xtc2")

# TODO I am only doing 100 events so $HOME is ok
ofname = "/cds/home/r/robertu/test.xtc2"

with open(ofname, "wb") as xtc2file:
    ds = DataSource(files=ifname)

    # Configure
    config = DgramEdit(transition_id=TransitionId.Configure)
    alg = AlgDef("raw", 1, 2, 3)
    det_def = DetectorDef(odetector, 'libpressio', "detnum1234")
    datadef = {
        "compressed": (np.uint8, 1),
        "npeaks": (np.uint16, 0),
        "row": (np.uint16, 1),
        "col": (np.uint16, 1),
        "shape": (np.uint16, 1),
    }
    det_out = config.Detector(det_def, alg, datadef, nodeId=nodeId,
                          namesId=namesId[odetector])

    runinfo_alg = AlgDef("runinfo", 0, 0, 1)
    runinfo_det = DetectorDef("runinfo", "runinfo", "")
    runinfodef = {
        "expt": (str, 1),
        "runnum": (np.uint32, 0),  # "mask": (np.uint8, 3),
    }
    runinfo_out = config.Detector(
        runinfo_det, runinfo_alg, runinfodef, nodeId=nodeId,
        namesId=namesId["runinfo"]
    )

    scan_alg = AlgDef("raw", 2, 0, 0)
    scan_det = DetectorDef("scan", "scan", "detnum1234")
    scandef = {
        "pixel_position": (np.float32, 4),
        "pixel_index_map": (np.int16, 4),
        "mask": (np.uint16, 3),
        "pf_dict": (str, 1),
    }
    scan_out = config.Detector(
        scan_det, scan_alg, scandef, nodeId=nodeId, namesId=namesId["scan"]
    )
    config_timestamp = 0  # TODO set this from the ds
    config.updatetimestamp(config_timestamp)
    save_dgramedit(config, outbuf, xtc2file)

    for r, run in enumerate(islice(ds.runs(), 1)):
        max_npeaks = 2048
        nevts = 10


        det = run.Detector(rdetector)
        row_array = np.zeros((max_npeaks,), dtype=np.uint16)
        col_array = np.zeros((max_npeaks,), dtype=np.uint16)
        for i, evt in enumerate(run.events()):
            if i == 0:
                # BeginRun
                beginrun = DgramEdit(transition_id=TransitionId.BeginRun, config=config, ts=config_timestamp + 1)
                runinfo_out.runinfo.expt = run.expt
                runinfo_out.runinfo.runnum = run.runnum
                beginrun.adddata(runinfo_out.runinfo)
                scan_out.raw.pixel_position = run.Detector('pixel_position')(evt)
                scan_out.raw.pixel_index_map = run.Detector('pixel_index_map')(evt)
                m = run.Detector('mask')(evt)
                print(m.dtype, m.shape)
                scan_out.raw.mask = m
                scan_out.raw.pf_dict = run.Detector('pf_dict')(evt)
                beginrun.adddata(scan_out.raw)
                save_dgramedit(beginrun, outbuf, xtc2file)

                # BeginStep
                beginstep = DgramEdit(transition_id=TransitionId.BeginStep, config=config, ts=config_timestamp + 2)
                save_dgramedit(beginstep, outbuf, xtc2file)

                # Enable
                enable = DgramEdit(transition_id=TransitionId.Enable, config=config, ts=config_timestamp + 3)
                save_dgramedit(enable, outbuf, xtc2file)
                current_timestamp = run.timestamp

            npeaks = det.raw.npeaks(evt)
            row_array[0:npeaks] = det.raw.row(evt)
            col_array[0:npeaks] = det.raw.col(evt)
            data = det.raw.calib(evt)

            comp = make_compressor(npeaks, row_array, col_array)
            compressed = comp.encode(data)
            cr = (data.size * data.itemsize) / len(compressed)
            print(f"doing run={r} event {i}, CR={cr}, data={data.shape}")

            # L1Accept
            d0 = DgramEdit(transition_id=TransitionId.L1Accept, config=config, ts=evt.timestamp)
            det_out.raw.row = row_array
            det_out.raw.col = col_array
            det_out.raw.npeaks = npeaks
            det_out.raw.compressed = compressed
            det_out.raw.shape = np.array(data.shape).astype(np.uint16)
            d0.adddata(det_out.raw)
            save_dgramedit(d0, outbuf, xtc2file)
            current_timestamp = evt.timestamp

        disable = DgramEdit(transition_id=TransitionId.Disable, config=config, ts=current_timestamp + 1)
        save_dgramedit(disable, outbuf, xtc2file)
        current_timestamp = config_timestamp + 3
        endstep = DgramEdit(transition_id=TransitionId.EndStep, config=config, ts=current_timestamp + 2)
        save_dgramedit(endstep, outbuf, xtc2file)
        endrun = DgramEdit(transition_id=TransitionId.EndRun, config=config, ts=current_timestamp + 3)
        save_dgramedit(endrun, outbuf, xtc2file)
