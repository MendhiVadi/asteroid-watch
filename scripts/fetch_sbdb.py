"""Fetch the full NEO catalogue from JPL SBDB Query API -> data/raw/sbdb_neo.json"""
import json, sys, urllib.request, urllib.parse
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "raw" / "sbdb_neo.json"
FIELDS = "spkid,full_name,pdes,name,a,e,i,om,w,ma,epoch,per,moid,H,diameter,albedo,pha,neo"
URL = "https://ssd-api.jpl.nasa.gov/sbdb_query.api?" + urllib.parse.urlencode(
    {"fields": FIELDS, "sb-kind": "a", "sb-group": "neo", "full-prec": "true"})

def fetch(tries=6):
    import time
    last = None
    for k in range(tries):
        try:
            with urllib.request.urlopen(URL, timeout=180) as r:
                return r.read()
        except Exception as e:          # IncompleteRead etc.: retry
            last = e
            print(f"attempt {k + 1} failed: {e!r}", flush=True)
            time.sleep(3 * (k + 1))
    raise last


def main():
    try:
        raw = fetch()
    except Exception as e:
        print("NETWORK_FAILURE:", repr(e)); sys.exit(2)
    d = json.loads(raw)
    OUT.write_bytes(raw)
    print("saved", OUT, len(raw), "bytes; count =", d.get("count"), "rows =", len(d["data"]))
    print("fields:", d["fields"]); print("signature:", d.get("signature"))

if __name__ == "__main__":
    main()
