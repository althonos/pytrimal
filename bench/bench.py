import argparse
import glob
import json
import os
import time
import statistics
import sys
import warnings
from itertools import islice

sys.path.insert(0, os.path.realpath(os.path.join(__file__, "..", "..")))

import numpy
import pandas
import rich.progress
from pytrimal import Alignment, AutomaticTrimmer

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--runs", default=3, type=int)
parser.add_argument("-d", "--data", required=True)
parser.add_argument("-o", "--output", required=True)
args = parser.parse_args()

with rich.progress.Progress(transient=True) as progress:
    warnings.showwarning = lambda msg, c, f, l, file=None, line=None: progress.print(msg)

    results = dict(results=[])
    examples = glob.glob(os.path.join(args.data, "example.*.AA.*.fasta"))

    example = Alignment.load(os.path.join(args.data, "example.014.AA.EggNOG.COG0591.fasta"))
    subsets = numpy.linspace(10, len(example.sequences), 10, dtype=numpy.int32)

    task = progress.add_task(total=len(subsets), description="Working...")
    for subset in progress.track(subsets, task_id=task):
        alignment = Alignment(example.names[:subset], example.sequences[:subset])
        task2 = progress.add_task(total=2, description="Backends...")
        for backend in progress.track(("sse", "generic", None), task_id=task2):
            times = []
            task3 = progress.add_task(total=args.runs, description="Repeat...")
            for run in progress.track(range(args.runs), task_id=task3):
                trimmer = AutomaticTrimmer(backend=backend, method="strict")
                t1 = time.time()
                trimmed = trimmer.trim(alignment)
                t2 = time.time()
                times.append(t2 - t1)
            progress.remove_task(task3)
            results["results"].append({
                "backend": backend,
                "residues": len(alignment.residues),
                "sequences": len(alignment.sequences),
                "times": times,
                "mean": statistics.mean(times),
                "stddev": statistics.stdev(times),
                "median": statistics.median(times),
                "min": min(times),
                "max": max(times),
            })
        progress.remove_task(task_id=task2)
    progress.remove_task(task_id=task)

with open(args.output, "w") as f:
    json.dump(results, f, sort_keys=True, indent=4)
