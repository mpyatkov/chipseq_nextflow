#!/usr/bin/python3

import sys
import random

debug_path = "/projectnb/wax-dk/max/G223_H3K27ac/bin/bigwig_for_hub.csv"

def header(hub_name):
    s =f"""
hub {hub_name} hub
shortLabel {hub_name}
longLabel {hub_name}
useOneFile on
email mpyatkov@bu.edu

genome mm9

track {hub_name}
type bigWig
container multiWig
shortLabel {hub_name}
longLabel {hub_name}
visibility full
aggregate none
autoScale group
showSubtrackColorOnUi on
maxHeightPixels 128:64:8
priority 1
"""
    return s.strip()

def render_track(hub_name, dpath, short_name, long_name, color):
    return f"""
            track {short_name}
            bigDataUrl {dpath}
            shortLabel {short_name}
            longLabel {long_name}
            windowingFunction mean
            parent {hub_name}
            graphTypeDefault bar
            type bigWig
            color {color}"""

if __name__ == "__main__":
    hub_name = sys.argv[1]   
    track_file = sys.argv[2]

    # hub_name="TEST_HUB1"
    # track_file= "/projectnb/wax-dk/max/G223_H3K27ac/bin/bigwig_for_hub.csv"
    with open(track_file, "r") as f:
        hdr = header(hub_name)
        result = [hdr]
        track_lines = list(map(lambda x: x.strip(),f.readlines()))
        for track_line in track_lines:
            dpath, shortname, longname, color = track_line.split(";")
            result.append(render_track(hub_name, dpath, shortname, longname, color))

        final_hub = "\n".join(result)
        with open("all_bigwig_group_autoscale_hub.txt", "w") as out:
            out.write(final_hub)

