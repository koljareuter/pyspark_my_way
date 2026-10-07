"""
Bidirectional lookup between KMOS group_id (synthetic) and the catalog
arm_target_id (= ESO PRO REFLEX SUFFIX = ESO OCS ARMx NAME in FITS headers).

Typical usage in a notebook, after _ext_gid_map and _Z_TABLE are built::

    from pyspark.utils.targetutils import TargetLookup
    targets = TargetLookup.from_gid_map(_ext_gid_map, "scripts/target_redshifts.csv")

    targets[90003]          # group_id  → info dict
    targets["101639"]       # arm_target_id (str) → info dict
    targets[101639]         # arm_target_id (int)  → info dict
    targets.tid(90003)      # → '101639'
    targets.gid("101639")   # → 90003  (or list if tid appears in multiple fields)
    print(targets[90003])   # pretty summary line
"""

import csv
import pathlib


class TargetLookup:
    """Translate between group_id and arm_target_id / field / redshift metadata."""

    def __init__(self):
        self._by_gid: dict[int, dict] = {}
        self._by_tid: dict[str, list[dict]] = {}

    # ── construction ────────────────────────────────────────────────────────────

    @classmethod
    def from_gid_map(cls, gid_map, csv_path=None):
        """Build from the in-memory _ext_gid_map dict and optionally a redshift CSV.

        Parameters
        ----------
        gid_map : dict
            ``{(field, arm_target_id): group_id}`` — the ``_ext_gid_map`` built
            in the sc-load cell.
        csv_path : str or Path, optional
            Path to ``target_redshifts.csv``.  When supplied, each record gains
            ``ra_deg``, ``dec_deg``, ``z``, ``z_type``, ``source_catalog``,
            ``source_id``, ``quality`` fields.
        """
        lut = cls()
        # build z-table from CSV
        zt: dict[tuple, dict] = {}
        if csv_path is not None:
            with open(csv_path, newline="", encoding="utf-8") as fh:
                for row in csv.DictReader(fh):
                    key = (row["field"], row["arm_target_id"])
                    zt[key] = {
                        "ra_deg":         float(row["ra_deg"])  if row["ra_deg"]  else None,
                        "dec_deg":        float(row["dec_deg"]) if row["dec_deg"] else None,
                        "z":              float(row["redshift"]) if row["redshift"] else None,
                        "z_type":         row["redshift_type"],
                        "source_catalog": row["source_catalog"],
                        "source_id":      row["source_id"],
                        "quality":        row["quality"],
                    }

        for (field, tid), gid in gid_map.items():
            tid_str = str(tid)
            rec = {
                "group_id": gid,
                "field":    field,
                "tid":      tid_str,
            }
            rec.update(zt.get((field, tid_str), {}))
            lut._by_gid[gid] = rec
            lut._by_tid.setdefault(tid_str, []).append(rec)

        return lut

    @classmethod
    def from_fits_dir(cls, directory, csv_path=None):
        """Build by scanning a directory tree for DATACUBE_COMBINED_*.fits files
        and reading ``ESO PRO REFLEX SUFFIX`` from each primary header.  Useful
        when the in-memory gid_map is not available (e.g. in a separate script).
        Assigns group_ids starting at 90000 in discovery order.
        """
        import glob, re
        from astropy.io import fits as _fits

        OLD_RE = re.compile(r"DATACUBE_COMBINED_(.+)_(IZ|YJ|HK|H|K)\.fits$")
        NEW_RE = re.compile(r"DATACUBE_COMBINED_(.+)\.fits$")

        lut = cls()
        gid_counter = [90000]
        seen: dict[tuple, int] = {}

        for fp in sorted(glob.glob(
                str(pathlib.Path(directory) / "**" / "DATACUBE_COMBINED_*.fits"),
                recursive=True)):
            fn = pathlib.Path(fp).name
            mo = OLD_RE.match(fn) or NEW_RE.match(fn)
            if not mo:
                continue
            with _fits.open(fp, memmap=False) as hf:
                hdr0 = hf[0].header
                suffix = str(hdr0.get("HIERARCH ESO PRO REFLEX SUFFIX", "")).strip()
                field  = str(hdr0.get("HIERARCH ESO OBS TARG NAME", "")).strip() or mo.group(1).rsplit("_", 1)[0]
                if not suffix:
                    # fall back: one OCS ARMx NAME keyword in DATA ext
                    dh = hf["DATA"].header
                    for k in dh.keys():
                        if "OCS" in k and "NAME" in k:
                            suffix = str(dh[k]).strip()
                            break
            if not suffix:
                continue
            key = (field, suffix)
            if key not in seen:
                seen[key] = gid_counter[0]
                gid_counter[0] += 1
                gid = seen[key]
                rec = {"group_id": gid, "field": field, "tid": suffix}
                lut._by_gid[gid] = rec
                lut._by_tid.setdefault(suffix, []).append(rec)

        if csv_path is not None:
            # enrich with redshift metadata
            zt = {}
            with open(csv_path, newline="", encoding="utf-8") as fh:
                for row in csv.DictReader(fh):
                    zt[(row["field"], row["arm_target_id"])] = row
            for rec in lut._by_gid.values():
                r = zt.get((rec["field"], rec["tid"]))
                if r:
                    rec.update({
                        "ra_deg":         float(r["ra_deg"])   if r["ra_deg"]  else None,
                        "dec_deg":        float(r["dec_deg"])  if r["dec_deg"] else None,
                        "z":              float(r["redshift"]) if r["redshift"] else None,
                        "z_type":         r["redshift_type"],
                        "source_catalog": r["source_catalog"],
                        "source_id":      r["source_id"],
                        "quality":        r["quality"],
                    })
        return lut

    # ── query interface ──────────────────────────────────────────────────────────

    def __getitem__(self, key):
        """Look up by group_id (int ≥90000) or arm_target_id (str or int <90000).

        Returns a dict.  If a tid matches records in multiple fields, returns a
        list of dicts; use ``.gid(tid, field=...)`` to disambiguate.
        """
        if isinstance(key, int) and key >= 90000:
            return self._by_gid[key]
        tid_str = str(key)
        hits = self._by_tid.get(tid_str)
        if hits is None:
            raise KeyError(key)
        return hits[0] if len(hits) == 1 else hits

    def tid(self, gid: int) -> str:
        """Return the arm_target_id string for a group_id."""
        return self._by_gid[gid]["tid"]

    def gid(self, tid, field: str = None) -> int:
        """Return the group_id for an arm_target_id.

        If the same tid appears in multiple fields, supply ``field`` to
        disambiguate; otherwise the first match is returned.
        """
        tid_str = str(tid)
        hits = self._by_tid.get(tid_str, [])
        if not hits:
            raise KeyError(tid_str)
        if field is not None:
            for h in hits:
                if h["field"] == field:
                    return h["group_id"]
            raise KeyError(f"tid={tid_str!r} not found in field={field!r}")
        return hits[0]["group_id"]

    def info(self, key) -> str:
        """One-line human-readable summary for a group_id or arm_target_id."""
        rec = self[key]
        if isinstance(rec, list):
            return "\n".join(self._fmt(r) for r in rec)
        return self._fmt(rec)

    @staticmethod
    def _fmt(rec: dict) -> str:
        z_str = f"z={rec['z']:.4f} [{rec.get('z_type','')}]" if rec.get("z") else "z=?"
        cat   = rec.get("source_catalog", "")
        return (
            f"gid={rec['group_id']}  tid={rec['tid']}  "
            f"field={rec['field']}  {z_str}"
            + (f"  cat={cat}" if cat else "")
        )

    def __repr__(self) -> str:
        return f"TargetLookup({len(self._by_gid)} entries)"

    def __len__(self) -> int:
        return len(self._by_gid)

    def __contains__(self, key) -> bool:
        if isinstance(key, int) and key >= 90000:
            return key in self._by_gid
        return str(key) in self._by_tid

    def summary_table(self):
        """Return a list of dicts sorted by group_id — handy for DataFrame construction."""
        return [self._by_gid[g] for g in sorted(self._by_gid)]
