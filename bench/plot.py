import argparse
import itertools
import json
import os
import re
import math

import numpy
import matplotlib.pyplot as plt
import scipy.stats
from numpy.polynomial import Polynomial
from palettable.cartocolors.qualitative import Bold_9


COLORS = dict(zip(
    ["None", "Generic", "MMX", "SSE2", "AVX2", "NEON"], 
    Bold_9.hex_colors
))
ORDER = {x:i for i,x in enumerate(COLORS)}

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True)
parser.add_argument("-o", "--output")
parser.add_argument("-s", "--show", action="store_true")
args = parser.parse_args()

with open(args.input) as f:
    data = json.load(f)
for result in data["results"]:
    if result["backend"] is None:
        result["backend"] = "None"
    elif result["backend"] == "generic":
        result["backend"] = "Generic"
    else:
        result["backend"] = result["backend"].upper()

plt.figure(1, figsize=(18, 6))

data["results"].sort(key=lambda r: r["statistic"])
for i, (statistic, statistic_group) in enumerate(
    itertools.groupby(data["results"], key=lambda r: r["statistic"])
):
    plt.subplot(1, 3, i + 1)
    plt.title(statistic)

    statistic_group = sorted(
        statistic_group, key=lambda r: (ORDER[r["backend"]], r["sequences"])
    )
    for backend, group in itertools.groupby(statistic_group, key=lambda r: r["backend"]):
        group = list(group)
        X = numpy.array([r["sequences"] for r in group])
        Y = numpy.array([r["mean"] for r in group])

        if statistic == "Gaps":
            reg = scipy.stats.linregress(X, Y)
            plt.scatter(
                X,
                Y,
                marker="+",
                color=COLORS[backend],
                label=f"{backend} (R²={reg.rvalue**2:.3f})",
            )
            plt.plot(
                [0, max(X)],
                [reg.intercept, reg.slope * max(X) + reg.intercept],
                color=COLORS[backend],
                linestyle="--",
                marker="",
            )
        else:
            p = Polynomial.fit(X, Y, 2)
            pX = numpy.linspace(0, max(r["sequences"] for r in group), 1000)
            # reg = scipy.stats.linregress(X, Y)
            # plt.plot([ 0, max(X) ], [ reg.intercept, reg.slope*max(X) + reg.intercept ], color=COLORS[backend], linestyle="--", marker="")
            # ci = [1.96 * r["stddev"] / math.sqrt(len(r["times"])) for r in group]
            plt.scatter(X, Y, marker="+", color=COLORS[backend], label=f"{backend}")
            plt.plot(pX, p(pX), color=COLORS[backend], linestyle="--")

    plt.legend()
    plt.xlabel("Number of sequences")
    plt.ylabel("Time (s)")

plt.tight_layout()
output = args.output or args.input.replace(".json", ".svg")
plt.savefig(output, transparent=True)
if args.show:
    plt.show()
