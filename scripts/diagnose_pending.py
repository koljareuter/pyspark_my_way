"""One-line diagnosis per unfinished OB: stage that failed + error."""
import argparse
import configparser
import glob
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument('config', nargs='?', default='scripts/kmos_server_new.cfg')
parser.add_argument('--log', default='logs/reduce_new.log')
parser.add_argument('--verbose', action='store_true', help='show all matching log lines too')
args = parser.parse_args()

cp = configparser.ConfigParser()
cp.read(args.config)
procdir = cp['pathing']['procdir'].strip()
subdirs = [s.strip() for s in re.split(r',|\n', cp['pathing']['subdirs']) if s.strip()]

pending = [ob for ob in subdirs
           if not glob.glob(os.path.join(procdir, ob, 'science', '*_rv*.fits'))]

# parse log: keep only the LAST run section per OB (pipeline may have run twice)
log_sections = {}
current = None
if os.path.exists(args.log):
    with open(args.log) as f:
        for line in f:
            m = re.search(r'Reducing (KMOS\.\S+)', line)
            if m:
                current = m.group(1).rstrip()
                log_sections[current] = []   # reset → keeps last run
            if current:
                log_sections[current].append(line.rstrip())

CALIB_KEYS = ['MASTER_DARK', 'MASTER_FLAT', 'XCAL', 'YCAL', 'DET_IMG_WAVE']

print(f"{'OB':<42}  {'calib':^6}  {'sci':^4}  {'failed at / note'}")
print("-" * 100)

for ob in pending:
    calib_dir = os.path.join(procdir, ob, 'calib')
    sci_dir   = os.path.join(procdir, ob, 'science')

    if os.path.isdir(calib_dir):
        have = os.listdir(calib_dir)
        calib_ok = sum(1 for k in CALIB_KEYS if any(k in f for f in have))
        calib_str = f"{calib_ok}/{len(CALIB_KEYS)}"
    else:
        calib_str = "none"

    if os.path.isdir(sci_dir):
        sci_n = len(os.listdir(sci_dir))
        sci_str = str(sci_n)
    else:
        sci_str = "-"

    lines = log_sections.get(ob, [])

    # find the fatal error from the LAST run
    fatal = ""
    for i, l in enumerate(lines):
        if '!!!! An error occured during' in l:
            m = re.search(r'during (\S+)', l)
            fatal = m.group(1) if m else "unknown"
            # grab the first non-!!!! context line after it
            break

    # if no fatal error, check if still running or some other note
    if not fatal:
        if not lines:
            fatal = "(no log entry)"
        elif any('Finished reducing' in l for l in lines):
            fatal = "completed? (check rv files)"
        elif any('sky subtraction' in l for l in lines):
            fatal = "sky_subtraction"
        else:
            last = next((l for l in reversed(lines) if '[INFO' in l or '[ERROR' in l), "")
            fatal = f"last: {last[-60:]}" if last else "(no error found)"

    print(f"{ob:<42}  {calib_str:^6}  {sci_str:^4}  {fatal}")

    if args.verbose:
        interesting = [l for l in lines if any(k in l for k in
                       ['!!!! An error', 'sky subtraction', 'File not found',
                        'Finished reducing', 'Finished calib', 'Finished recon'])]
        for l in interesting:
            print(f"    {l}")
