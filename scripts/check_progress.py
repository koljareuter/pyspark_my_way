import argparse
import configparser
import glob
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument('config', nargs='?', default='scripts/kmos_server_new.cfg')
parser.add_argument('--list', action='store_true', help='print every unfinished OB')
args = parser.parse_args()

cp = configparser.ConfigParser()
cp.read(args.config)

rawdir  = cp['pathing']['rawdir'].strip()
procdir = cp['pathing']['procdir'].strip()
subdirs = [s.strip() for s in re.split(r',|\n', cp['pathing']['subdirs']) if s.strip()]

done    = [ob for ob in subdirs if glob.glob(os.path.join(procdir, ob, 'science', '*_rv*.fits'))]
pending = [ob for ob in subdirs if ob not in done]

print(f"rawdir : {rawdir}")
print(f"procdir: {procdir}")
print(f"OBs in config : {len(subdirs)}")
print(f"OBs with rv.* : {len(done)}")
print(f"Progress       : {100 * len(done) / len(subdirs):.1f}%")

if args.list:
    print(f"\nNot yet reduced ({len(pending)}):")
    for ob in pending:
        print(f"  {ob}")
