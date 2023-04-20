#!/usr/bin/env python
#
# DrKinfold: Produce the DrForna *.drf file format from Kinfold simulations.
# Written by Stefan Badelt (stefan.badelt@univie.ac.at)
#
import os
import sys
import glob
import argparse
import numpy as np
from subprocess import Popen, PIPE
from multiprocessing import Pool

from . import __version__
from .utils import (parse_vienna_stdin, 
                    get_drf_output_times, 
                    combine_drfs)


def syscall_kinfold(name, seq,
                    start = None,
                    stop = None,
                    fpt  = False,
                    rect = False,
                    time = 5000,
                    num = 1,
                    ratemodel = 'Metropolis',
                    moves = 'single-base-pair',
                    noLP = False,
                    logML = False,
                    dangle = 2,
                    temp = 37,
                    params = None,
                    glen = 1,
                    grow = None,
                    # Output
                    erange = 20, # Kinfold Default.
                    lmin = False,
                    silent = False,
                    force = False):
    """Perform a system-call of the program ``Kinfold``.

    The print the results into a file and return the respective filename. This
    wrapper will return the output files, including ``STDIN`` and ``STDERR``.

    Args:
      ...

    Returns:
      ...
    """
    klfile = name + '.log'
    if force and os.path.exists(klfile):
        os.remove(klfile)

    kinput = seq + '\n'
    syscall = ['Kinfold']
    syscall.extend(['--num', str(int(num))])
    syscall.extend(['--time', str(time)])
    syscall.extend(['--log', name])
    syscall.extend(['--cut', str(erange)])
    #syscall.extend(['--seed', "62159=58010=26254"])
    if params:
        syscall.extend(["--Par", params])
    if dangle != 2:
        syscall.extend(['--dangle', str(dangle)])
    if temp != 37:
        syscall.extend(['-T', str(temp)])
    if ratemodel == 'Kawasaki':
        pass
    elif ratemodel == 'Metropolis':
        syscall.extend(['--met'])
    else:
        raise NotImplementedError('unknown rate model')

    if lmin:
        syscall.extend(['--lmin'])

    if not fpt: 
        # NOTE: fpt switches first passage time off (!!!!)
        syscall.extend(['--fpt'])

    if rect:
        syscall.extend(['--rect'])

    if not logML:
        # NOTE: logML switches logarithmic multiloop evaluation off (!!!!)
        syscall.extend(['--logML'])

    if silent:
        syscall.extend(['--silent'])

    if noLP:
        syscall.extend(['--noLP'])

    if moves == 'single-base-pair':
        syscall.extend(['--noShift'])
    elif moves == 'shift':
        pass

    if grow:
        syscall.extend(['--glen', str(glen)])
        syscall.extend(['--grow', str(grow)])

    if start:
        assert isinstance(start, str)
        syscall.extend(['--start'])
        kinput += start + '\n'

    if stop:
        syscall.extend(['--stop'])
        if isinstance(stop, str):
            kinput += stop + '\n'
        else:
            assert isinstance(stop, list)
            kinput += '\n'.join(stop)
    return kinput, syscall

def sub_kinfold(*kargs, **kwargs):
    name, seq = kargs
    kinput, kcall = syscall_kinfold(*kargs, **kwargs)
    print('[in progress:] ' + ' '.join(kcall))
    kefile = name + '.err'
    with open(kefile, 'w') as ehandle:
        with Popen(kcall, stdin = PIPE, stdout = PIPE, 
                   bufsize = 1, 
                   universal_newlines = True, 
                   stderr = ehandle) as proc:
            proc.stdin.write(kinput)
            proc.stdin.close()
            for line in proc.stdout:
                yield line
    return

def run_kinfold(times, basename, seq, num, atupernuc, atupersec, totkftime, temperature, params):
    idc = 0
    with open(f'{basename}.drf', 'w') as drf:
        drf.write(f"id time occupancy structure energy\n")
        t, nsim = 0, 0
        for line in sub_kinfold(basename, seq, num = num, glen = 1, temp = temperature,
                                params = params, grow = atupernuc, time = totkftime, 
                                erange = 999999):
            [ss, en, st] = line.split()[0:3]
            stime = float(st)
            # Add all drf output times until the give time step
            while t < len(times) and times[t]*atupersec <= stime:
                drf.write(f'{idc:>5d} {times[t]:13.9f} 1 {ss} {float(en):6.2f}\n')
                t += 1
            if len(line.split()) == 4:
                if t < len(times):
                    assert np.isclose(times[t]*atupersec, stime)
                    drf.write(f'{idc:>5d} {times[t]:13.9f} 1 {ss} {float(en):6.2f}\n')
                    t += 1
                assert t == len(times)
                t = 0
                nsim += 1
                print(f'[status update:] Done with simulation {nsim} in {basename}.drf. ', end = '\r')
            idc += 1
    print(f'[Done:] Kinfold call for {basename} finished after {nsim} simulations. ')

def parse_drkinfold_args(parser):
    parser.add_argument('--version', action = 'version', 
            version = '%(prog)s ' + __version__)

    parser.add_argument("--name", default = '', metavar = '<str>',
            help = """Name your output files, name the header of your plots, etc.
            this option overwrites the fasta-header.""")

    parser.add_argument("--tmpdir", default = 'drkinfold', action = 'store', metavar = '<str>',
            help = """Specify path for storing Kinfold output files.""")

    parser.add_argument("-p", "--processes", type = int, default = 0,
            help="Number of individual Kinfold system calls. By default, only existing data is processed.")

    parser.add_argument("-c", "--cpus", type = int, default = None,
            help="Maximal number of cpus used for threading.")

    parser.add_argument("-n", "--num", type = int, default = 1,
            help="Number of simulations per Kinfold call.")

    parser.add_argument("--k0", type = float, default = 1e5, metavar = '<flt>',
            help = """Arrhenius rate constant. Adjust to relate free energy
            changes to experimentally determined folding time [atu/s].""")

    parser.add_argument("--t-ext", type = float, default = 0.02, metavar = '<flt>',
            help = """Time per nucleotide extension (the inverse of the transcription rate)
            [s/nt].""")

    parser.add_argument("--t-end", type = float, default = 30, metavar = '<flt>',
            help = "Post-transcriptional simulation time [s].")

    parser.add_argument("--t-lin", type = int, default = 10, metavar = '<int>',
            help = """Evenly space output *--t-lin* times during transcription on a linear time scale.""")

    parser.add_argument("--t-log", type = int, default = 30, metavar = '<int>',
            help = """Evenly space output *--t-log* times after transcription on a logarithmic time scale.""")

    parser.add_argument("-T", "--temp", type = float, default = 37.0, 
        metavar = '<flt>',
        help = 'Rescale energy parameters to a temperature of temp C.')

    parser.add_argument("-P", "--paramFile", action = "store", default = None,
        metavar = '<str>',
        help = """Read energy parameters from a parameter file, instead of 
        using the default ViennaRNA parameter set.""")

    return

def main():
    """Call Kinfold for co-transcriptional folding and provide *.drf output format.
    """
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'DrKinfold: Cotranscriptional folding using Kinfold.')
    parse_drkinfold_args(parser)
    args = parser.parse_args()

    # Read Input & Update Arguments
    name, seq = parse_vienna_stdin(sys.stdin)
    if args.name:
        name = args.name
    print(f'>{name}\n{seq}')

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
        for data in glob.glob(f'{args.tmpdir}/{name}.*.drf'):
            ndata = data.split('/')[-1]
            *pre, nfid, suf = ndata.split('.')
            fid = max(fid, int(nfid)+1)
    else:
        os.mkdir(args.tmpdir)

    #
    # Do all the Kinfold calculations.
    #
    if args.processes:
        # Conversion factors between seconds and Kinfold's internal time units.
        atupersec = args.k0
        atupernuc = atupersec * args.t_ext
        totkftime = atupernuc * len(seq) + atupersec * args.t_end
        with Pool(processes = args.cpus) as q:
            multiple_results = [q.apply_async(run_kinfold, 
                (times, f'{args.tmpdir}/{name}.{fid+x:03d}', seq, 
                 args.num, atupernuc, atupersec, totkftime, args.temp, args.paramFile)) for x in range(args.processes)]
            [res.get() for res in multiple_results]

    #
    # Combine all drf files from individual simulations to one lage output file.
    #
    combine_drfs(f'{args.tmpdir}/{name}*.drf', f'{name}.drf', len(seq), times, use_counts = False)

if __name__ == '__main__':
    main()

