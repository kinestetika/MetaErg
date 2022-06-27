import re
from metaerg.run_and_read import subsystems_data

CUES = {}
current_subsystem = None
for line in subsystems_data.subsystem_data().split('\n'):
    line = line.strip()
    if line.startswith("#") or not len(line):
        continue
    elif line.startswith(">"):
        current_subsystem = line[1:]
    elif current_subsystem is not None:
        CUES[line] = current_subsystem


def match(descriptions) -> str:
    for d in descriptions:
        for cue, subsystem in CUES.items():
            if len(d) > len(cue) + 20:
                continue
            match = re.search(r'\b' + cue + r'\b', d)
            if match and match.start() < 10:
                return subsystem
    return ''
