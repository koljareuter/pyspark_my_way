import argparse
import datetime
import json
import os
import re
import shutil

KMOS_TIMESTAMP_RE = re.compile(r"(\d{4}-\d{2}-\d{2})T(\d{2}):(\d{2}):(\d{2})")


def parse_kmos_timestamp(name):
    """Extract the UT timestamp embedded in a KMOS dataset_name/filename,
    e.g. 'KMOS.2026-01-15T01:00:33.152' -> datetime(2026, 1, 15, 1, 0, 33)."""
    match = KMOS_TIMESTAMP_RE.search(name)
    date, hour, minute, second = match.groups()
    year, month, day = (int(x) for x in date.split("-"))
    return datetime.datetime(year, month, day, int(hour), int(minute), int(second))


def night_of(timestamp):
    """The observing night a UT timestamp belongs to.

    This is just the UT calendar date of the dataset's anchor timestamp: on
    the server, /local_storage/EMPOWER_raw_data/cosmos/2026-01-12/ holds
    exposures timestamped 05:17-12:37 UT -- i.e. the existing convention
    already spans past UT noon without splitting into two folders, so a
    noon cutoff would disagree with it. Grouping itself doesn't depend on
    this label at all (files are grouped via the JSON's own dataset tree,
    not by per-file date), so this only decides the folder name -- and it's
    picked to match cosmos/<date> exactly."""
    return timestamp.date()


def collect_dataset_files(node, files=None):
    """Recursively collect (path, category) pairs referenced anywhere in a
    dataset's job tree (main/associated input files, and the same for every
    nested main/associated input job)."""
    if files is None:
        files = []
    for entry in node.get("main_input_files", []) + node.get("associated_input_files", []):
        files.append((entry["name"], entry["category"]))
    for job in node.get("main_input_jobs", []) + node.get("associated_input_jobs", []):
        collect_dataset_files(job["dataset"], files)
    return files


def is_static(path):
    """Static calibration reference files live under a 'static' subfolder and
    are shared across every night, so they shouldn't be sorted per dataset."""
    return "static" in os.path.normpath(path).split(os.sep)


def group_datasets(data):
    """Map each dataset's own name to its observing night and the
    deduplicated list of file paths that belong together (the same
    calibration frame is often referenced by several sub-jobs of one
    dataset).

    Grouping key is the dataset_name itself, not a computed night -- this
    adapts to whatever granularity the JSON provides. kmos_nights.json has
    one dataset per *night* (bundling several OBs internally), so this
    produces one folder per night, same as before. kmos_ob_structure.json
    has one dataset per *OB*, so this produces one folder per OB. Keying by
    a computed night would silently collapse multiple same-night OBs onto
    one folder, overwriting all but the last."""
    groups = {}
    for dataset in data["datasets"]:
        name = dataset["dataset_name"]
        night = night_of(parse_kmos_timestamp(name))
        files = collect_dataset_files(dataset["dataset"])
        seen = dict.fromkeys(path for path, category in files if not is_static(path))
        groups[name] = {"dataset_name": name, "night": night, "files": list(seen)}
    return groups


def resolve_source_path(path):
    """Return the real, existing path for a raw-data JSON path.

    kmos_nights.json records raw files under .../EMPOWER_raw_data/input_all/,
    but on the server the actual files live under
    .../EMPOWER_raw_data/cosmos/<night>/ instead (confirmed via `ls` on
    empower.hq.eso.org) -- input_all may be a stale or partial pool. Fall
    back to the cosmos/<night>/ layout, derived from the timestamp in the
    filename itself, before giving up on a file."""
    if os.path.isfile(path):
        return path
    marker = os.path.join("EMPOWER_raw_data", "input_all")
    if marker not in path:
        return path
    basename = os.path.basename(path)
    match = KMOS_TIMESTAMP_RE.search(basename)
    if match is None:
        return path
    night = parse_kmos_timestamp(basename).date().isoformat()
    candidate = path.replace(marker, os.path.join("EMPOWER_raw_data", "cosmos", night))
    return candidate if os.path.isfile(candidate) else path


def check_dest_is_safe(dest_root, groups):
    """Refuse to sort into a destination that sits inside the raw-data tree
    itself -- sorting must never place anything among (or overwrite) the
    original files it links to."""
    dest_real = os.path.realpath(dest_root)
    for group in groups.values():
        for path in group["files"]:
            source_dir = os.path.realpath(os.path.dirname(resolve_source_path(path)))
            if source_dir == dest_real or source_dir.startswith(dest_real + os.sep):
                raise SystemExit(f"Refusing to sort into {dest_root!r}: it contains raw source data ({source_dir}).")


def sort_files(groups, dest_root, use_copy=False):
    """Link (or, with use_copy, copy) each group's files into
    dest_root/<dataset_name>/. Existing files/links are left untouched --
    nothing is ever deleted or overwritten, and the source files are never
    touched (a symlink only reads the source path, it never modifies it)."""
    check_dest_is_safe(dest_root, groups)
    for group in groups.values():
        dest_dir = os.path.join(dest_root, group["dataset_name"])
        os.makedirs(dest_dir, exist_ok=True)
        for path in group["files"]:
            real_path = resolve_source_path(path)
            if not os.path.isfile(real_path):
                print("  !! missing: " + path)
                continue
            dest_path = os.path.join(dest_dir, os.path.basename(real_path))
            if os.path.lexists(dest_path):
                continue
            if use_copy:
                shutil.copy2(real_path, dest_path)
            else:
                os.symlink(os.path.abspath(real_path), dest_path)


def write_filelist(groups, filelist_path):
    """Write the real (non-static) source paths for the selected nights to
    a plain text file, one per line -- for use with e.g.
    `rsync -avR --files-from=<filelist> user@host:/ localdir/` to pull down
    just those files for local testing."""
    with open(filelist_path, "w") as f:
        for group in groups.values():
            for path in group["files"]:
                f.write(path + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Group KMOS raw files by dataset (one night in kmos_nights.json, one OB in kmos_ob_structure.json).")
    parser.add_argument("--json", default=os.path.join(os.path.dirname(__file__), "kmos_nights.json"))
    parser.add_argument("--dest", help="If given, link (or copy) each selected dataset's files into <dest>/<dataset_name>/ instead of only printing.")
    parser.add_argument("--copy", action="store_true", help="Copy files into --dest instead of symlinking them (symlinking is the default).")
    parser.add_argument("--filelist", help="Write the real source paths for the selected datasets to this file, one per line "
                                            "(for `rsync -avR --files-from=<filelist> user@host:/ localdir/`).")
    parser.add_argument("--limit", type=int, default=3, help="Only act on the first N datasets, chronologically (0 = all). Ignored if --nights is given. Default: 3.")
    parser.add_argument("--nights", help="Comma-separated list of nights to restrict to (selects every dataset whose own night falls on one of these dates), "
                                          "e.g. 2026-01-12,2026-01-13,2026-01-14.")
    args = parser.parse_args()

    with open(args.json, "r") as file:
        data = json.load(file)

    groups = group_datasets(data)
    all_names = sorted(groups)

    if args.nights:
        wanted = {datetime.date.fromisoformat(s.strip()) for s in args.nights.split(",")}
        selected = [n for n in all_names if groups[n]["night"] in wanted]
        found_nights = {groups[n]["night"] for n in selected}
        missing = wanted - found_nights
        if missing:
            print("!! no dataset found for night(s): " + ", ".join(sorted(m.isoformat() for m in missing)))
    elif args.limit:
        selected = all_names[: args.limit]
    else:
        selected = all_names
    selected_groups = {name: groups[name] for name in selected}

    for name in selected:
        group = groups[name]
        print(f"\n{group['dataset_name']}  (night {group['night']}, {len(group['files'])} files)")
        for path in group["files"]:
            print('  ... ' + path)
    if len(selected) < len(all_names):
        print(f"\n... ({len(all_names) - len(selected)} more datasets not selected, use --limit 0 to include all)")

    if args.filelist:
        write_filelist(selected_groups, args.filelist)
        n_files = sum(len(g["files"]) for g in selected_groups.values())
        print(f"\nWrote {n_files} source path(s) for {len(selected)} dataset(s) to {args.filelist}")

    if args.dest:
        sort_files(selected_groups, args.dest, use_copy=args.copy)
