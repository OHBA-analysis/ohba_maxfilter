import os
import sys
import argparse
import subprocess


def _add_headpos(cmd, args):
    """stimates and stores head position parameters but does not transform data"""
    if ('headpos' in args) and args['headpos']:
        hp_name = outname.replace('.fif','_headpos.log')
        hp_name = os.path.join(args['outdir'], hp_name)
        cmd += ' -hp {0}'.format(hp_name)
    return cmd

def _add_movecomp(cmd, args):
    """estimates head movements and transforms data to reference head position in continuous raw data"""
    # Add option for Automatic head-movement compensation
    if ('movecomp' in args) and args['movecomp']:
        cmd += ' -movecomp'
    return cmd

def _add_movecompinter(cmd, args):
    """estimates head movements and transforms data to reference head position in continuous raw data"""
    # Add option for Automatic head-movement compensation
    if ('movecompinter' in args) and args['movecompinter']:
        cmd += ' -movecomp inter'
    return cmd

def _add_hpie(cmd, args):
    """sets the error limit for hpi coil fitting (def 5 mm)"""
    if ('hpie' in args) and (args['hpie'] is not None):
        cmd += ' -hpie {0}'.format(args['hpie'])
    return cmd

def _add_hpig(cmd, args):
    """sets the g-value limit for hpi coil fitting (def 0.98)"""
    if ('hpig' in args) and (args['hpig'] is not None):
        cmd += ' -hpig {0}'.format(args['hpig'])
    return cmd

def _add_autobad(cmd, args):
    """sets automated bad channel detection on: scan the whole raw data file, off: no autobad"""
    if ('autobad' in args) and args['autobad']:
        cmd += ' -autobad on'
    else:
        cmd += ' -autobad off'
    return cmd

def _add_bad(cmd, args):
    """sets the list of static bad channels (logical chnos, e.g.: 0323 1042 2631)"""
    if ('bads' in args) and args['bads'] is not None:
        cmd += ' -bad {0}'.format(args['bads'])
    return cmd

def _add_tsss(cmd, args):
    """Add all tsss related args"""
    # Add option for Temporal Extension of maxfilter
    if ('tsss' in args) and args['tsss']:
        cmd += ' -st {0} -corr {1}'.format(args['st'], args['corr'])
    return cmd

def _add_trans(cmd, args):
    """transforms the data into head position in <fiff_file>"""
    if ('trans' in args) and args['trans'] is not None:
        cmd += ' -trans {0} -force'.format(args['trans'])
    return cmd

def _add_inorder(cmd, args):
    """sets the order of the inside expansion"""
    if ('inorder' in args) and args['inorder'] is not None:
        cmd += ' -in {0}'.format(args['inorder'])
    return cmd

def _add_outorder(cmd, args):
    """sets the order of the outside expansion"""
    if ('outorder' in args) and args['outorder'] is not None:
        cmd += ' -out {0}'.format(args['outorder'])
    return cmd

def _add_ctc(cmd, args):
    """uses the cross-talk matrix in <ctcfile>"""
    if ('ctc' in args) and args['ctc']:
        cmd += ' -ctc {0}'.format(args['ctc'])
    return cmd

def _add_cal(cmd, args):
    """uses the fine-calibration in <calfile>"""
    if ('cal' in args) and args['cal']:
        cmd += ' -cal {0}'.format(args['cal'])
    return cmd

def _add_scanner(cmd, args):
    if ('scanner' in args) is False:
        return cmd
    if args['scanner'] == 'VectorView':
        # Just use defaults
        pass
    elif args['scanner'] == 'VectorView2':
        args['cal'] = '/net/aton/meg_pool/neuromag/databases/sss/sss_cal_3026_171220.dat'
        cmd = _add_cal(cmd, args)
        args['ctc'] = '/net/aton/meg_pool/neuromag/databases/ctc/ct_sparse.fif'
        cmd = _add_ctc(cmd, args)
    elif args['scanner'] == 'Neo':
        args['cal'] = '/net/aton/meg_pool/data/TriuxNeo/system/sss/sss_cal.dat'
        cmd = _add_cal(cmd, args)
        args['ctc'] = '/net/aton/meg_pool/data/TriuxNeo/system/ctc/ct_sparse.fif'
        cmd = _add_ctc(cmd, args)

    return cmd


def run_maxfilter(infif, outfif, args, logfile_tag=''):

    #print(args)
    #if ('maxpath' in args) and ('maxpath' is not None):
    #    maxpath = args['maxpath']
    #else:
    #    maxpath = '/neuro/bin/util/maxfilter-2.2'

    basecmd = '{maxpath} -f {infif} -o {outfif}'

    # --------------
    # Format Maxfilter options

    if ('tsss' in args) and args['tsss']:
        outfif = outfif.replace('.fif','tsss.fif')
    else:
        outfif = outfif.replace('.fif','sss.fif')

    # Create base command
    cmd = basecmd.format(maxpath=args['maxpath'], infif=fif, outfif=outfif)

    cmd = _add_headpos(cmd, args)
    cmd = _add_movecomp(cmd, args)
    cmd = _add_movecompinter(cmd, args)
    cmd = _add_hpie(cmd, args)
    cmd = _add_hpig(cmd, args)
    cmd = _add_autobad(cmd, args)
    cmd = _add_bad(cmd, args)
    cmd = _add_tsss(cmd, args)
    cmd = _add_trans(cmd, args)
    cmd = _add_inorder(cmd, args)
    cmd = _add_outorder(cmd, args)
    if ('scanner' in args) and args['scanner'] is not None:
        cmd = _add_scanner(cmd, args)
    else:
        cmd = _add_ctc(cmd, args)
        cmd = _add_cal(cmd, args)

    # Add verbose and logfile
    stdlog = outname.replace('.fif','{0}.log'.format(logfile_tag))
    stdlog = os.path.join(args['outdir'], stdlog)
    errlog = outname.replace('.fif','{0}_err.log'.format(logfile_tag))
    errlog = os.path.join(args['outdir'], errlog)

    # Set tee to capture both stdout and stderr into separate files
    # https://stackoverflow.com/a/692407
    cmd += ' -v > >(tee -a {stdlog}) 2> >(tee -a {errlog} >&2)'.format(stdlog=stdlog, errlog=errlog)

    # --------------
    # Run Maxfilter

    if args['dryrun']:
        # Dry-run just prints the command
        print(cmd)
        print('\n')
    else:
        print(cmd)
        print('\n')
        # Call maxfilter in a subprocess
        os.system("bash -c '{}'".format(cmd))

    return outfif, stdlog

# -------------------------------------------------

def run_multistage_maxfilter(infif, outbase, args):
    """
    https://imaging.mrc-cbu.cam.ac.uk/meg/Maxfilter
    https://imaging.mrc-cbu.cam.ac.uk/meg/maxbugs

    Camb advice - don't use trans with movecomp
                - don't use autobad with headpos or movecomp
                - don't use autobad with st
    """

    # --------------------------------------
    # Stage 1 - Find Bad Channels

    outfif = outbase.format('autobad_.fif')
    outlog = outbase.format('autobad.log')

    if os.path.exists(outfif):
        os.remove(outfif)

    # Fixed Args
    stage1_args = {'autobad': True}
    # User args
    for key in ['inorder', 'outorder', 'hpie', 'hpig', 'maxpath',
                'scanner', 'ctc', 'cal', 'dryrun', 'overwrite', 'outdir']:
        if key in args:
            stage1_args[key] = args[key]

    outfif, outlog = run_maxfilter(infif, outfif, stage1_args, '_autobad')

    if args['dryrun'] is False:
        # Read in bad channels from logfile
        with open(outlog,'r') as f:
            txt = f.readlines()

        for ii in range(len(txt)):
            if txt[ii][:19] == 'Static bad channels':
                bads = txt[ii].split(': ')[1].split(' ')
                bads = ' '.join([b.strip('\n') for b in bads])
                break
            else:
                bads = None
    else:
        bads = None

    # --------------------------------------
    # Stage 2 - Signal Space Separation

    outfif = outbase.format('.fif')
    outlog = outbase.format('.log')

    if os.path.exists(outfif):
        os.remove(outfif)

    # Fixed Args
    stage2_args = {'autobad': None, 'bads': bads}
    # User args
    for key in ['tsss', 'st', 'corr' , 'inorder', 'outorder', 'maxpath',
                'scanner', 'ctc', 'cal', 'dryrun','overwrite', 'hpig', 'hpie',
                'movecomp', 'movecompinter', 'headpos', 'outdir']:
        if key in args:
            stage2_args[key] = args[key]

    outfif, outlog = run_maxfilter(infif, outfif, stage2_args, '_tsss')

    # --------------------------------------
    # Stage 3 - Translate to reference file
    if ('trans' in args) and args['trans'] is not None:

        infif = outfif  # input is output from previous stage
        outfif = outbase.format('trans.fif')
        outlog = outbase.format('trans.log')

        if os.path.exists(outfif):
            os.remove(outfif)

        # Fixed Args
        stage3_args = {'autobad': None}
        # User args
        for key in ['scanner','ctc','cal', 'dryrun','overwrite','trans','outdir']:
            if key in args:
                stage3_args[key] = args[key]

        outfif, outlog = run_maxfilter(infif, outfif, stage3_args, '_trans')


# -------------------------------------------------

parser = argparse.ArgumentParser(description='Batch run Maxfilter on some fif files.')
parser.add_argument('files', type=str,
                    help='plain text file containing full paths to files to be processed')
parser.add_argument('outdir', type=str,
                    help='Path to output directory to save data in')

parser.add_argument('--maxpath', type=str, default='/neuro/bin/util/maxfilter-2.2',
                     help='Path to maxfilter command to use')

parser.add_argument('--mode', type=str, default='standard',
                    help="Running mode for maxfilter. Either 'standard' or 'multistage'")

parser.add_argument('--headpos', action='store_true',
                    help='Output additional head movement parameter file')
parser.add_argument('--movecomp', action='store_true',
                    help='Apply movement compensation')
parser.add_argument('--movecompinter', action='store_true',
                    help='Apply movement compensation on data with intermittent HPI')
parser.add_argument('--autobad', action='store_true',
                    help='Apply automatic bad channel detection')
parser.add_argument('--trans', type=str, default=None,
                    help='Transforms the data to the head position in defined file')

parser.add_argument('--tsss', action='store_true',
                    help='Apply temporal extension of maxfilter')
parser.add_argument('--st', type=float, default=10,
                    help='Data buffer length for TSSS processing')
parser.add_argument('--corr', type=float, default=0.98,
                    help='Subspace correlation limit for TSSS processing')

parser.add_argument('--inorder', type=int, default=None,
                    help='Set the order of the inside expansion')
parser.add_argument('--outorder', type=int, default=None,
                    help='Set the order of the outside expansion')

parser.add_argument('--hpie', type=int, default=None,
                    help="sets the error limit for hpi coil fitting (def 5 mm)")
parser.add_argument('--hpig', type=float, default=None,
                    help="ets the g-value limit (goodness-of-fit) for hpi coil fitting (def 0.98))")

parser.add_argument('--scanner', type=str, default=None,
                    help="Set CTC and Cal for the OHBA scanner the dataset was collected with (VectorView, VectorView2 or Neo). This overrides the --ctc and --cal options.")
parser.add_argument('--ctc', type=str, default=None,
                    help='Specify cross-talk calibration file')
parser.add_argument('--cal', type=str, default=None,
                    help='Specify fine-calibration file')


parser.add_argument('--overwrite', action='store_true',
                    help="Overwrite previous output files if they're in the way")
parser.add_argument('--dryrun', action='store_true',
                    help="Don't actually run anything, just print commands that would have been run")

args = parser.parse_args()
print(args)

# -------------------------------------------------

print('\n\nOHBA-Maxfilter\n\n')

with open(args.files, 'r') as f:
    infifs = f.readlines()
infifs = [fif.strip('\n') for fif in infifs]

good_fifs = [1 for ii in range(len(infifs))]
for idx, fif in enumerate(infifs):
    if os.path.isfile(fif) == False:
        good_fifs[idx] = 0
        print('File not found: {0}'.format(fif))

if args.trans is not None:
    if os.path.isfile(args.trans) is False:
        sys.exit('Trans file not found ({0})'.format(args.trans))

print('Processing {0} files'.format(sum(good_fifs)))
print('Outputs saving to: {0}\n\n'.format(args.outdir))

# -------------------------------------------------

for idx, fif in enumerate(infifs):

    # --------------
    # Format input and output files and run some checks
    print('Processing run {0}/{1} : {2}'.format(idx+1,len(infifs), fif))

    # Skip run if we couldn't find the input file on disk
    if good_fifs[idx] == 0:
        print('Input file not found, skipping run ({0})'.format(fif))
        continue

    # Make an output name : myscan.fif -> myscan_tsss.fif
    outname = os.path.split(fif)[1]

    # Outputfile is output dir + output name
    outfif = os.path.join(args.outdir, outname)

    # Skip run if the output exists and we don't want to overwrite
    if os.path.isfile(outfif) and (args.overwrite is False):
        print('Existing output found, skipping run ({0})'.format(fif))
        continue

    # Delete previous output if it output exists and we do want to overwrite
    if os.path.isfile(outfif) and args.overwrite:
        print('Deleting previous output: {0}'.format(outfif))
        os.remove(outfif)

    if args.mode == 'standard':
        outfif, outlog = run_maxfilter(infifs[idx], outfif, vars(args))
    elif args.mode == 'multistage':
        outbase = outfif[:-4] + '_{0}'
        run_multistage_maxfilter(infifs[idx], outbase, vars(args))

print('\nProcessing complete. OHBA-and-out.\n')
