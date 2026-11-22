"""Read the plotted curves of Figures 3-5 from the paper's vector figures.

One-off provenance tool. Figures 3-5 of the paper (sig2_lomihi.eps, sig2n.eps,
sig2n_example2.eps in the paper's LaTeX source; MATLAB R2014b output) are
vector figures: every curve is a polyline whose points are printed in the file.

Usage:  python digitize_published_figures.py <figure.eps> <output.csv>
Output: CSV with one row per polyline point: Line (polyline number), Colour
(RGB), Width, Dash (dash pattern), Point (index in the polyline), X, Y (inner
coordinates of the figure, y pointing down).

The axis tick labels are drawn as glyph outlines, so the data coordinates are
not calibrated here: verify/verify_against_paper.m fits the log-log axis
scales of each panel (two parameters per axis) to the published-run summary
and reports the residuals in pixels.
"""
import csv
import re
import sys

POINT = re.compile(r'^(-?[\d.]+) (-?[\d.]+) ([ML])$')


def polylines(lines):
    """Yield (colour, width, dash, points) for every path drawn with M/L."""
    colour, width, dash, points = '', '', '', []
    for line in lines:
        text = line.strip()
        match = POINT.match(text)
        if match:
            if match.group(3) == 'M' and points:
                yield colour, width, dash, points
                points = []
            points.append((float(match.group(1)), float(match.group(2))))
            continue
        if text in ('S', 'GR') and points:
            yield colour, width, dash, points
            points = []
        if text.endswith(' RC'):
            colour = text[:-3]
        elif text.endswith(' LW'):
            width = text[:-3]
        elif text.endswith(' setdash') or text.endswith(' d'):
            dash = text.rsplit(' ', 1)[0]
        elif text == 'GR':
            dash = ''


def main(eps_file, csv_file):
    with open(eps_file, encoding='latin-1') as handle:
        lines = handle.read().splitlines()
    rows = []
    for number, (colour, width, dash, points) in enumerate(polylines(lines), start=1):
        for index, (x, y) in enumerate(points, start=1):
            rows.append([number, colour, width, dash, index, x, y])
    with open(csv_file, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['Line', 'Colour', 'Width', 'Dash', 'Point', 'X', 'Y'])
        writer.writerows(rows)
    print(f'{len(rows)} points written to {csv_file}')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
