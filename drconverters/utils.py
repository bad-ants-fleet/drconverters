import os
from glob import glob
import numpy as np


def parse_vienna_stdin(stdin, chars='ACGUNTacgunt'):
    """Parse name and sequence from file with fasta format.

    Only one input-sequence is allowed at a time.

    Args:
      stdin (list): Input to parse, ususally :obj:`sys.stdin`
      chars (string, optional): Allowed characters in a sequence.

    Returns:
      str, str: name and sequence.
    """
    name = None
    seq = ''
    for line in stdin:
        line = line.strip()
        if line[0] == '>':
            assert name is None, "Only single-sequence fasta format supported!"
            assert seq == '', "Only single-sequence fasta format supported!"
            name = line.split()[0][1:]
        else:
            warn = set([x for x in line if x not in chars])
            if len(warn):
                raise SystemExit(f"Unsupported character(s) in RNA: {warn}")
            seq += line
    return name, seq

def get_drf_output_times(seqlen, t1, t8, t_lin, t_log):
    times = np.array([0])
    for nuc in range(seqlen):
        ntime = np.linspace(times[-1], times[-1] + t1, t_lin + 1)
        times = np.concatenate([times, ntime[1:]])
    ntime = np.logspace(np.log10(times[-1]), np.log10(times[-1] + t8), t_log + 1)
    return np.concatenate([times, ntime[1:]])

def combine_drfs(drffiles, oname, seqlen, times, use_counts = False, get_kp8 = False):
    #
    # Collect data from all drf output files.
    #
    cdict = {t: dict() for t in range(len(times))} # Counts
    edict = {t: dict() for t in range(len(times))} # Energy
    idict, nid = dict(), 0 # Identity
    nfiles, nsim = 0, 0
    for data in glob(drffiles):
        nfiles += 1
        with open(data) as dat:
            t = 0
            for i, line in enumerate(dat):
                if i ==0:
                    continue
                _, time, occ, ss, en = line.split()
                # NOTE: If the line below breaks, then probably because of old
                # data that was generated using a different t-lin and/or t-log.
                assert np.isclose(float(time), times[t])
                cdict[t][ss] = cdict[t].get(ss, 0) + 1
                edict[t][ss] = int(round(float(en)*100))
                future = '.' * (seqlen - len(ss))
                if ss+future not in idict:
                    idict[ss+future] = nid
                    nid += 1
                t += 1
                if t == len(times):
                    t = 0
                    nsim += 1
    print(f'[collecting data:] Parsed {nsim} simulations from {nfiles} files.')
    #
    # Write the final vector into a separate file for potential further analysis
    #
    if get_kp8:
        st = len(times)-1
        with open(f'{oname}.kp8', 'w') as df:
            for s in sorted(edict[st], key = lambda x: edict[st][x]):
                df.write(f'{s} {cdict[st][s]:>5d} {edict[st][s]/100:6.2f}\n')
    #
    # Transform the count dict into an occupancy dict!
    #
    if use_counts:
        odict = cdict
    else:
        odict = {t: dict() for t in range(len(times))} # Occupancy
        for t in sorted(cdict):
            for ss in cdict[t]:
                odict[t][ss] = cdict[t][ss]/nsim
    #
    # Write *.drf output file.
    #
    if os.path.exists(oname):
        print(f"[WARNING:] Overwriting existing file: {oname}")
    with open(oname, 'w') as df:
        df.write(f"id time occupancy structure energy\n")
        for t in sorted(odict):
            time = times[t]
            for (ss, en) in sorted(edict[t].items(), key = lambda x: x[1]):
                if ss not in odict[t]:
                    continue
                occu = odict[t][ss]
                future = '.' * (seqlen - len(ss))
                ni = idict[ss+future]
                if use_counts:
                    df.write(f'{ni:5d} {times[t]:03.3f} {occu:5d} {ss} {en/100:6.2f}\n')
                else:
                    df.write(f'{ni:5d} {times[t]:03.3f} {occu:03.4f} {ss} {en/100:6.2f}\n')
 
