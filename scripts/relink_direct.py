"""Replace two-level kmos_ob_sorted symlinks with direct links.

kmos_ob_sorted/OB/X.fits -> input_all/X.fits -> input3/X.fits (fails in kmos_flat CPL)
becomes
kmos_ob_sorted/OB/X.fits -> input3/X.fits (single hop, bypasses input_all)
"""
import os

OB_DIR = '/home/kreuter/kmos_ob_sorted'

updated = 0
already_direct = 0
broken = 0

for ob in sorted(os.listdir(OB_DIR)):
    ob_path = os.path.join(OB_DIR, ob)
    if not os.path.isdir(ob_path):
        continue
    for fname in sorted(os.listdir(ob_path)):
        fpath = os.path.join(ob_path, fname)
        if not os.path.islink(fpath):
            continue
        current_target = os.readlink(fpath)
        real_target = os.path.realpath(fpath)
        if not os.path.exists(real_target):
            broken += 1
            continue
        if current_target == real_target:
            already_direct += 1
            continue
        os.unlink(fpath)
        os.symlink(real_target, fpath)
        updated += 1

print(f'Updated to direct: {updated}')
print(f'Already direct:    {already_direct}')
print(f'Broken (skipped):  {broken}')
