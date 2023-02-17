#!/usr/bin/env python
#
# DrKinefold: Produce the DrForna *.drf file format from Kinefold simulations.
# Written by Stefan Badelt (stefan.badelt@univie.ac.at)
#
import os
import sys
import math
import string
import argparse
import numpy as np
import subprocess as sub
from random import randint
from glob import glob

import RNA
from .utils import (parse_vienna_stdin, 
                   get_drf_output_times, 
                   combine_drfs)

def parse_kinefold_structure(line1, line2):
    """ For example:

     A U[C G G G G G]C U C U[G U U G]G U U[C U C C C G^C A A C]G C U A C C
     - - -5- - - - - - - - - -6- - - - - - -5' - - - - -6' - - - - - - - -
     =>
     . . ( ( ( ( ( ( . . . . ( ( ( ( . . . ) ) ) ) ) ) ) ) ) ) . . . . . .
    """
    hdict = {}
    havail = [True for _ in range(26+1)]
    sseq, sstr = '', ''
    lx, lc, c = '', '', '.'
    for j, (x, y) in enumerate(zip(line1, line2)):
        if x == y:
            assert x == ' '
            continue
        if x in ('ACGU'):
            # Append the previous x, c to strings.
            sseq += lx
            sstr += lc
            lx = x
            lc = c
        elif x == '[' or x == '^':
            # We need to resolve the structure character
            c = '|'
        elif x == ']':
            # An unpaired strech will follow
            c = '.'
        if c == '|' and y not in (' ', '-', "'"):
            # y = the full ID of the current helix.
            jj = j
            while line2[jj+1] not in (' ', '-', "'"):
                y += line2[jj+1] 
                jj += 1
    
            if y not in hdict:
                myid = havail.index(True) # First occurrence of True in ID list
                hdict[y] = [myid, j, None]
                hid = hdict[y][0]
                c = '(' if hid == 0 else string.ascii_letters[hid - 1 + 26]
                havail[hid] = False
            else:
                hdict[y][2] = j
                hid = hdict[y][0]
                c = ')' if hid == 0 else string.ascii_letters[hid - 1 ]
                havail[hid] = True
            lc = c
    sseq += lx
    sstr += lc
    return sseq, sstr

def rnm_to_drf(rnmfile, drffile, times, t_ext):
    """Translates Kinefold *.rnm file to DrForna *.drf file.
    """ 
    idc = 0
    with open(rnmfile, 'r') as rnm, open(drffile, 'w') as drf, open(rnmfile + '.log', 'w') as log:
        drf.write(f"id time occupancy structure energy\n")
        t, delay = 0, None
        lsstr = '.'
        for i, line in enumerate(rnm, 1):
            line = line.rstrip()
            if i == 1:
                assert line[0] == '<'
                _, name = line.split()
                continue
            if i == 2:
                seq = line
                continue
            if i % 2: # process sequence info
                subseq = line
            else:
                subseq, info = subseq.split('|')
                substr = line.split('H')[0]
                en, eunit, _, _, ms, tunit = info.split()[0:6]
                assert eunit == 'kcal/mol'
                assert tunit == 'ms,'
                sseq, sstr = parse_kinefold_structure(subseq.rstrip(), substr.rstrip())
                sstr = RNA.db_from_ptable(RNA.ptable(sstr, RNA.BRACKETS_ANY))
                stime = float(ms) * 10**(-3)
                if delay is None:
                    delay = len(sstr) * t_ext
                stime += delay
                log.write(f'# {float(ms) * 10**(-3):13.9f} {stime:13.9f} {sstr} {float(en):6.2f}\n')
                while t < len(times) and times[t] <= stime:
                    tlen = round(times[t]/t_ext) if np.isclose(round(times[t]/t_ext, 9), 
                                                               times[t]/t_ext) else math.ceil(times[t]/t_ext)
                    if tlen == 0:
                        tlen = 1
                    elif tlen > len(seq):
                        tlen = len(seq)
                    if len(lsstr) < tlen:
                        lsstr += '.' * (tlen-len(lsstr))
                    elif len(lsstr) > tlen:
                        # If this happens, check that only unpaired nucleotides
                        # are removed at the end... otherwise it's a problem!
                        raise ValueError('This case has never been observed before!')
                        lsstr = lsstr[:tlen]
                    assert len(lsstr) == tlen
                    drf.write(f'{idc:>5d} {times[t]:13.9f} 1 {lsstr} {float(en):6.2f}\n')
                    t += 1
                lsstr = sstr
                idc += 1
        while t < len(times):
            drf.write(f'{idc:>5d} {times[t]:13.9f} 1 {lsstr} {float(en):6.2f}\n')
            t += 1
    return seq, name

def parse_drkinefold_args(parser):
    parser.add_argument("--name", default = '', metavar = '<str>',
            help = """Name your output files, name the header of your plots, etc.
            this option overwrites the fasta-header.""")

    parser.add_argument("--tmpdir", default = 'drkinefold', action = 'store', metavar = '<str>',
            help = """Specify path for storing Kinefold output files.""")

    parser.add_argument("-p", "--processes", type = int, default = 0,
            help="Number of individual Kinefold system calls. By default, only existing data is processed.")

    parser.add_argument("--t-ext", type = float, default = 0.02, metavar = '<flt>',
            help = """Time per nucleotide extension (the inverse of the transcription rate)
            [s/nt].""")

    parser.add_argument("--t-end", type = float, default = 30, metavar = '<flt>',
            help = "Post-transcriptional simulation time [s].")

    parser.add_argument("--t-lin", type = int, default = 10, metavar = '<int>',
            help = """Evenly space output *--t-lin* times during transcription on a linear time scale.""")

    parser.add_argument("--t-log", type = int, default = 30, metavar = '<int>',
            help = """Evenly space output *--t-log* times after transcription on a logarithmic time scale.""")
    return

def get_kinefold_input(name, i, seq, t_ext, t_end):
    wdir = os.getcwd()
    return f"""\
{randint(1, 10000)}	# random seed
{wdir}/{name}.w
{wdir}/{name}.w
{wdir}/{name}.{i:03d}.rnm
{wdir}/{name}.w
{wdir}/{name}.w
{wdir}/{name}.w
{wdir}/{name}.dat
0		# 0=RNA ; 1=DNA
6.3460741	# helix minimum free energy in kcal/mol: 6.3460741=10kT
10000000	# NA
{(len(seq)*t_ext + t_end) * 10**3:.0f}      # folding time requested in msec
1		# pseudoknots   1=yes 0=no 
0		# entanglements	1=yes 0=no
2 {t_ext*10**3}	# simulation type: 1=renaturation; 2 20=cotrans. @ 20msec/nt 
		# add T i j k or F i j k options here 
"""

def main():
    """Translate Kinefold cotranscriptional folding output to DrForna input format.
    """
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'DrKinefold: Produce DrForna input from Kinefold simulation output.')
    parse_drkinefold_args(parser)
    args = parser.parse_args()

    if args.processes and not os.path.exists('kinefold_long_static'):
        raise SystemExit(f'Kinfold executable "kinefold_long_static" not found.')

    #
    # Read Input & Update Arguments
    #
    name, seq = parse_vienna_stdin(sys.stdin)
    if args.name:
        name = args.name
    print(f'>{name}\n{seq}')
    # NOTE: Adding ext/end, as they are necessary to adjust Kinefold simulations ...
    #name = f'{name}_ext-{args.t_ext}_end-{args.t_end}'
    
    #
    # Prepare the output times in the *.drf file format.
    #
    times = get_drf_output_times(len(seq), args.t_ext, args.t_end, args.t_lin, args.t_log)

    #
    # Put everything in one directory, update the file ID in case there are
    # existing simulations.
    #
    fid = 1 # Set initial file ID according to what can already be found in tmpdir.
    if os.path.exists(args.tmpdir):
        for data in glob(f'{args.tmpdir}/{name}*.rnm'):
            ndata = data.split('/')[-1]
            *pre, nfid, suf = ndata.split('.')
            fid = max(fid, int(nfid)+1)
    else:
        os.mkdir(args.tmpdir)
    
    #
    # The *.dat file seems to be required input?!?
    #
    if args.processes:
        datf = os.path.join(args.tmpdir, f'{name}.dat')
        if os.path.exists(datf):
            print(f"WARNING: Overwriting existing file: {datf}")
        with open(datf, 'w') as dat:
            dat.write(f'< {name}\n')
            dat.write(f'{seq}\n')

    #
    # Do --processes separate simulations.
    #
    for i in range(fid, args.processes+fid):
        print(f'[in progress:] Calling Kinefold #{i}.')
        infile = os.path.join(args.tmpdir, f'{name}.{i:03d}.in')
        with open(infile, 'w') as k:
            k.write(get_kinefold_input(f'{args.tmpdir}/{name}', i, seq, args.t_ext, args.t_end))
        kcall = ['./kinefold_long_static', infile, '-noprint']
        sub.run(kcall, capture_output = True) 
    if args.processes: # clean up 
        os.remove(f'{args.tmpdir}/{name}.w')
        os.remove(f'{args.tmpdir}/{name}.i')

    #
    # Convert all rnmfiles to drffiles using the current time vector.
    #
    for rnmfile in glob(f'{args.tmpdir}/{name}.*.rnm'):
        drffile = rnmfile[:-3]+'drf'
        kseq, kname = rnm_to_drf(rnmfile, drffile, times, args.t_ext)
        assert kseq == seq and kname == name

    #
    # Combine all drf files from individual simulations to one lage output file.
    #
    combine_drfs(f'{args.tmpdir}/{name}*.drf', f'{name}.drf', len(seq), times, use_counts = False)
    return

if __name__ == '__main__':
    main()

