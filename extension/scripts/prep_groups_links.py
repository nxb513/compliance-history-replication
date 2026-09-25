"""Build plant groups (chronic water violators and siblings) and FRS program links for the primary sample.
Recovered from the inline commands used on 2026-09-24/25 and saved as a script for reproducibility (writes to a temp path when --check is given)."""
import sys, duckdb, pandas as pd
ROOT = "."
check = "--check" in sys.argv
py = pd.read_parquet(f"{ROOT}/data/interim/one_bad_plant_panel.parquet"); pr = py[py.nplant >= 3]
pl = pr.groupby("REG").agg(pnorm=("pnorm", "first"), st=("STATE_CODE", "first"), major=("is_major", "max"),
                           yrs_bad=("any_e90", "sum"), water_rate=("any_e90", "mean"), nyr=("yr", "nunique"))
pl["chronic"] = pl.yrs_bad >= 6
pl["sib_chronic"] = (pl.groupby("pnorm").chronic.transform("sum") - pl.chronic.astype(int)) > 0
con = duckdb.connect(); con.register("pl", pl.reset_index()[["REG"]])
links = con.execute(f"""SELECT f.PGM_SYS_ACRNM, f.PGM_SYS_ID, f.REGISTRY_ID FROM read_csv('{ROOT}/data/interim/frs/FRS_PROGRAM_LINKS.csv', all_varchar=true) f
       JOIN pl ON f.REGISTRY_ID = pl.REG WHERE f.PGM_SYS_ACRNM IN ('AIR','AIRS/AFS','RCRAINFO','ICIS','NPDES','TRIS','EIS')""").df()
if check:
    old_pl = pd.read_parquet(f"{ROOT}/extension/plant_groups_primary.parquet"); old_l = pd.read_parquet(f"{ROOT}/extension/frs_links_primary.parquet")
    same_pl = old_pl.sort_index().equals(pl.sort_index())
    key = ["PGM_SYS_ACRNM", "PGM_SYS_ID", "REGISTRY_ID"]
    same_l = old_l.sort_values(key).reset_index(drop=True).equals(links.sort_values(key).reset_index(drop=True))
    print("plant groups identical:", same_pl, "| FRS links identical:", same_l, "| rows", len(pl), len(links))
else:
    pl.to_parquet(f"{ROOT}/extension/plant_groups_primary.parquet"); links.to_parquet(f"{ROOT}/extension/frs_links_primary.parquet")
    print(pl.shape, int(pl.chronic.sum()), int(pl.sib_chronic.sum()), links.groupby("PGM_SYS_ACRNM").REGISTRY_ID.nunique().to_dict())
