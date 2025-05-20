import os, sys
import optparse
from subprocess import Popen, PIPE, STDOUT

# Global variable for options
opt = None
args = None

def parseOptions():
    global opt, args

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    parser.add_option('', '--obsName', dest='OBSNAME', type='string', default='costhetaZ1',
                      help='Observable name (e.g. inclusive, pT4l, eta4l, etc.)')
    parser.add_option('', '--obsBins', dest='OBSBINS', type='string',
                      default='|-1.0|-0.75|-0.50|-0.25|0.0|0.25|0.50|0.75|1.0|',
                      help='Bin boundaries separated by "|", e.g. "|0|50|100|"')
    parser.add_option('', '--year', dest='YEAR', type='string', default='2022',
                      help='Year: 2016, 2017, 2018, or Full')
    parser.add_option('', '--unblind', action='store_true', dest='UNBLIND', default=False,
                      help='Use real data (unblind)')
    parser.add_option('', '--fitOnly', action='store_true', dest='FITONLY', default=False,
                      help='Run only the fit step')

    (opt, args) = parser.parse_args()

    if not opt.OBSBINS and opt.OBSNAME != 'inclusive':
        parser.error('Bin boundaries not specified for differential measurement. Exiting...')
        sys.exit(1)

def processCmd(cmd, quiet=0):
    print(f">>> Running command: {cmd}")
    output = ''
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT, universal_newlines=True)
    for line in iter(p.stdout.readline, ''):
        output += line
        if not quiet:
            print(line, end='')
    p.stdout.close()
    if p.wait() != 0:
        print(f"!!! Command failed with status {p.returncode}: {cmd}")
        raise RuntimeError(f"{cmd!r} failed, exit status: {p.returncode}")
    return output

def safe_chdir(path):
    try:
        os.chdir(path)
        print(f">>> Changed directory to: {path}")
    except FileNotFoundError:
        print(f"!!! Directory not found: {path}")
        sys.exit(1)

def pipeline():
    print(">>> pipeline() is running")

    observableBins = {0: (opt.OBSBINS.split("|")[1:-1]), 1: ['0', 'inf']}[opt.OBSBINS == 'inclusive']
    obsName = opt.OBSNAME
    obsBins = opt.OBSBINS
    year = opt.YEAR
    unblind = opt.UNBLIND
    fitOnly = opt.FITONLY

    print('============================================================')
    print('=== RUNNING THE WORKFLOW FOR THE FOLLOWING CONFIGURATION ===')
    print(f'=== obsName: {obsName}, year: {year}, bins: {obsBins}, unblind: {unblind} ===')
    print('============================================================')

    if not fitOnly:
        safe_chdir('./coefficients')
        processCmd(f'python3 -u RunCoefficients.py --obsName "{obsName}" --obsBins "{obsBins}" --year "{year}"')

        safe_chdir('../templates')
        processCmd(f'python3 -u RunTemplates.py --obsName "{obsName}" --obsBins "{obsBins}" --year "{year}"')

        safe_chdir('../fit')

    else:
        safe_chdir('./fit')

    processCmd(f'python3 -u addConstrainedModel.py --obsName "{obsName}" --year "{year}"')

    runFidCmd = f'python3 -u RunFiducialXS.py --obsName "{obsName}" --obsBins "{obsBins}" --year "{year}"'
    if unblind:
        runFidCmd += ' --unblind'
    processCmd(runFidCmd)

    processCmd(f'python3 -u impacts.py --obsName "{obsName}" --year "{year}"')

    print('=== PLOTTING ===')
    safe_chdir('../LHScans')
    llScanCmd = f'python3 -u plot_LLScan.py --obsName "{obsName}" --obsBins "{obsBins}" --year "{year}"'
    if unblind:
        llScanCmd += ' --unblind'
    processCmd(llScanCmd)

    safe_chdir('../coefficients')
    processCmd(f'python3 -u RunPlotCoefficients.py --obsName "{obsName}" --obsBins "{obsBins}" --year "{year}"')

    safe_chdir('../templates')
    processCmd(f'python3 -u plot_templates.py --obsName "{obsName}" --year "{year}"')

    safe_chdir('..')
    processCmd(f'./copy_to_www.sh {obsName} {year}')

if __name__ == "__main__":
    parseOptions()
    pipeline()
