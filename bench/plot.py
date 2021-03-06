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
from palettable.cartocolors.qualitative import Bold_4


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

plt.figure(1, figsize=(6, 6))

plt.subplot(1, 1, 1)
data["results"].sort(key=lambda r: (r["backend"], r["sequences"]))
for color, (backend, group) in zip(
    Bold_4.hex_colors, itertools.groupby(data["results"], key=lambda r: r["backend"])
):
    group = list(group)
    X = numpy.array([r["sequences"] for r in group])
    Y = numpy.array([r["mean"] for r in group])

    p = Polynomial.fit(X, Y, 2)
    pX = numpy.linspace(0, max(r["sequences"] for r in group), 1000)
    # reg = scipy.stats.linregress(X, Y)
    # plt.plot([ 0, max(X) ], [ reg.intercept, reg.slope*max(X) + reg.intercept ], color=color, linestyle="--", marker="")
    # ci = [1.96 * r["stddev"] / math.sqrt(len(r["times"])) for r in group]
    plt.scatter(X, Y, marker="+", color=color, label=f"{backend}")
    plt.plot(pX, p(pX), color=color, linestyle="--")

plt.legend()
plt.xlabel("Number of sequences")
plt.ylabel("Time (s)")

# plt.subplot(1, 2, 2)
# data["results"].sort(key=lambda r: (r["backend"], r["residues"]))
# for color, (backend, group) in zip(
#     Bold_4.hex_colors, itertools.groupby(data["results"], key=lambda r: r["backend"])
# ):
#     group = list(group)
#     X = numpy.array([r["residues"] for r in group])
#     Y = numpy.array([r["mean"] for r in group])
#     p = Polynomial.fit(X, Y, 2)
#     # reg = scipy.stats.linregress(X, Y)
#     # plt.plot([ 0, max(X) ], [ reg.intercept, reg.slope*max(X) + reg.intercept ], color=color, linestyle="--", marker="")
#     # ci = [1.96 * r["stddev"] / math.sqrt(len(r["times"])) for r in group]
#     plt.scatter(X, Y, marker="+", color=color, label=f"{backend}")
#     plt.plot(X, p(X), color=color, linestyle="--")
#
# plt.legend()
# plt.xlabel("Number of residues")
# plt.ylabel("Time (s)")

plt.tight_layout()
output = args.output or args.input.replace(".json", ".svg")
plt.savefig(output, transparent=True)
if args.show:
    plt.show()
