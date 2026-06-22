#!/usr/bin/env python3
"""Assemble constraint_explorer.ipynb from _constraint_explorer_src.py cell markers.

Cells split on lines starting with "# %%". "# %% [markdown]" -> markdown cell;
its body has the leading "# " comment prefix stripped.
"""
import json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "_constraint_explorer_src.py")
OUT = os.path.join(HERE, "constraint_explorer.ipynb")


def parse_cells(text):
    lines = text.splitlines()
    cells = []
    cur_type = None
    cur = []
    # skip leading header block before first marker
    started = False
    for ln in lines:
        m = re.match(r"^# %%(.*)$", ln)
        if m:
            if started:
                cells.append((cur_type, cur))
            started = True
            cur_type = "markdown" if "[markdown]" in m.group(1) else "code"
            cur = []
        elif started:
            cur.append(ln)
    if started:
        cells.append((cur_type, cur))
    return cells


def to_md(lines):
    out = []
    for ln in lines:
        if ln.startswith("# "):
            out.append(ln[2:])
        elif ln == "#":
            out.append("")
        else:
            out.append(ln)
    # trim leading/trailing blanks
    while out and out[0] == "":
        out.pop(0)
    while out and out[-1] == "":
        out.pop()
    return out


def to_code(lines):
    while lines and lines[0] == "":
        lines.pop(0)
    while lines and lines[-1] == "":
        lines.pop()
    return lines


def mk_source(lines):
    # nbformat wants list of strings each ending in newline except last
    if not lines:
        return []
    return [l + "\n" for l in lines[:-1]] + [lines[-1]]


def main():
    with open(SRC) as f:
        text = f.read()
    cells = parse_cells(text)
    nb_cells = []
    for ctype, body in cells:
        if ctype == "markdown":
            src = mk_source(to_md(body))
            nb_cells.append({"cell_type": "markdown", "metadata": {}, "source": src})
        else:
            src = mk_source(to_code(body))
            if not src:
                continue
            nb_cells.append({"cell_type": "code", "metadata": {}, "execution_count": None,
                             "outputs": [], "source": src})
    nb = {
        "cells": nb_cells,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3"},
        },
        "nbformat": 4, "nbformat_minor": 5,
    }
    with open(OUT, "w") as f:
        json.dump(nb, f, indent=1)
    print(f"wrote {OUT} with {len(nb_cells)} cells")


if __name__ == "__main__":
    main()
