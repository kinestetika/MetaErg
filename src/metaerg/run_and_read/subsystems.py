class SubSystem:
    def __init__(self, name:str):
        self.name = name
        self.targets = []
        self.hits = {}

    def get_cues_as_hash(self) -> {}:
        return {target: self.name for target in self.targets}


def match_subsystem(descr, subsystem_cues):
    for cue, subsystem in subsystem_cues:
        if len(descr) > len(cue) + 20:
            return None, None
        match = re.search(r'\b' + cue + r'\b', descr)
        if match and match.start() < 10:
            return cue, subsystem
    return None, None



def init_subsystems() -> {}:
    subsystem_hash = {}
    current_subsystem = None
    for line in '''>[secondary-metabolites]
            # (populated by antismash)
            >[hydrocarbon degradation]
            # (populated by cant-hyd)
            >[ribosome]
            ribosomal protein L1
            ribosomal protein L2
            ribosomal protein L3
            ribosomal protein L4
            ribosomal protein L5
            ribosomal protein L6
            ribosomal protein L9
            ribosomal protein L11
            ribosomal protein L11 methyltransferase
            ribosomal protein L13
            ribosomal protein L14
            ribosomal protein L15
            ribosomal protein L16
            ribosomal protein L18
            ribosomal protein L19
            ribosomal protein L20
            ribosomal protein L21
            ribosomal protein L22
            ribosomal protein L23
            ribosomal protein L24
            ribosomal protein L25
            ribosomal protein L27
            ribosomal protein L29
            ribosomal protein L31
            ribosomal protein L34
            ribosomal protein L35
            ribosomal protein S1
            ribosomal protein S2
            ribosomal protein S3
            ribosomal protein S4
            ribosomal protein S5
            ribosomal protein S6
            ribosomal protein S7
            ribosomal protein S8
            ribosomal protein S9
            ribosomal protein S10
            ribosomal protein S11
            ribosomal protein S12
            ribosomal protein S13
            ribosomal protein S14
            ribosomal protein S15
            ribosomal protein S16
            ribosomal protein S17
            ribosomal protein S18
            ribosomal protein S19
            ribosomal protein S20
            16S ribosomal RNA
            23S ribosomal RNA
            5S ribosomal RNA
            ribosomal protein PSRP-3
            small ribosomal subunit biogenesis GTPase RsgA
            ribosomal protein S12 methylthiotransferase RimO
            ribosomal protein S6--L-glutamate ligase
            ribosome-binding factor RbfA
            '''.split('\n'):
        line = line.strip()
        if line.startswith("#") or not len(line):
            continue
        if line.startswith(">"):
            current_subsystem = SubSystem(line[1:])
            subsystem_hash[current_subsystem.name] = current_subsystem
            continue
        current_subsystem.targets.append(line)
    return subsystem_hash
