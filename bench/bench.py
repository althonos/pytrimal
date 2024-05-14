import argparse
import glob
import json
import os
import platform
import time
import statistics
import sys
import warnings
import typing
from itertools import islice
from typing import Dict, Callable, List

import numpy
import pandas
import rich.progress

sys.path.insert(0, os.path.realpath(os.path.join(__file__, "..", "..")))

from pytrimal import (
    _trimal,
    Alignment,
    BaseTrimmer,
    AutomaticTrimmer,
    ManualTrimmer,
    OverlapTrimmer,
)

if typing.TYPE_CHECKING:
    from pytrimal._trimal import COMPUTE_PLATFORM

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--runs", default=3, type=int)
parser.add_argument("-d", "--data", required=True)
parser.add_argument("-o", "--output", required=True)
args = parser.parse_args()


PLATFORMS = [None]
if _trimal._SSE2_RUNTIME_SUPPORT:
    PLATFORMS.append("sse2")
if _trimal._AVX2_RUNTIME_SUPPORT:
    PLATFORMS.append("avx2")
if _trimal._NEON_RUNTIME_SUPPORT:
    PLATFORMS.append("neon")

STATISTIC: Dict[str, Callable[["COMPUTE_PLATFORM"], BaseTrimmer]] = {
    "Gaps": lambda platform: ManualTrimmer(gap_threshold=0.5, platform=platform),
    "Similarity": lambda platform: ManualTrimmer(
        similarity_threshold=0.5, platform=platform
    ),
    "Overlap": lambda platform: OverlapTrimmer(
        sequence_overlap=60, residue_overlap=0.5, platform=platform
    ),
}

with rich.progress.Progress(transient=True) as progress:
    rich.print("[bold green]Benchmarking[/]: [cyan italic]{}[/]".format(" ".join(map(str, PLATFORMS))))

    warnings.showwarning = lambda msg, c, f, l, file=None, line=None: None

    results: Dict[str, List[Dict[str, object]]] = dict(results=[])
    examples = glob.glob(os.path.join(args.data, "example.*.AA.*.fasta"))

    example = Alignment.load(
        os.path.join(args.data, "example.014.AA.EggNOG.COG0591.fasta")
    )
    subsets = numpy.linspace(10, len(example.sequences), 10, dtype=numpy.int32)

    task0 = progress.add_task(total=len(STATISTIC), description="Statistic (...)")
    for statistic, get_trimmer in progress.track(STATISTIC.items(), task_id=task0):
        progress.update(task_id=task0, description=f"Statistic ({statistic})")
        task1 = progress.add_task(total=len(PLATFORMS), description=" Platform (...)")
        for backend in progress.track(PLATFORMS, task_id=task1):
            progress.update(task_id=task1, description=f" Platform ({backend})")
            trimmer = get_trimmer(backend)
            task2 = progress.add_task(
                total=len(subsets), description=f"  Subset (0/{len(example.sequences)})"
            )
            for subset in progress.track(subsets, task_id=task2):
                progress.update(
                    task_id=task2,
                    description=f"  Subset ({subset}/{len(example.sequences)})",
                )
                times = []
                alignment = Alignment(
                    example.names[:subset], example.sequences[:subset]
                )
                task3 = progress.add_task(
                    total=args.runs, description=f"  Repeat (0/{args.runs})"
                )
                for run in progress.track(range(args.runs), task_id=task3):
                    progress.update(
                        task_id=task3, description=f"  Repeat ({run}/{args.runs})"
                    )
                    t1 = time.time()
                    trimmed = trimmer.trim(alignment)
                    t2 = time.time()
                    times.append(t2 - t1)
                progress.remove_task(task3)
                results["results"].append(
                    {
                        "backend": backend,
                        "residues": len(alignment.residues),
                        "sequences": len(alignment.sequences),
                        "times": times,
                        "mean": statistics.mean(times),
                        "stddev": statistics.stdev(times),
                        "median": statistics.median(times),
                        "min": min(times),
                        "max": max(times),
                        "statistic": statistic,
                        "trimmer": repr(trimmer),
                    }
                )
            progress.remove_task(task_id=task2)
        progress.remove_task(task_id=task1)
    progress.remove_task(task_id=task0)

with open(args.output, "w") as f:
    json.dump(results, f, sort_keys=True, indent=4)
