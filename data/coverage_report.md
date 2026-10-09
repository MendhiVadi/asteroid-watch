# NEO coverage report

Source catalogue: JPL SBDB Query API (`sb-kind=a`, `sb-group=neo`), fetched 2026-10-09; raw: `data/raw/sbdb_neo.json`.

## Headline

- SBDB total known NEO asteroids: **42,670**
- CSV (`neo_v2.csv`) rows: 90,836; unique (id,name) objects: **27,423**
- CSV objects matched to an SBDB NEO: **27,379** (99.8%) -> 27,181 distinct SBDB objects (198 SBDB objects appear under >1 CSV id, i.e. duplicated/renumbered designations)
- CSV objects NOT found in SBDB NEO list: **44** (36 confirmed by SBDB as no longer NEO-class, 8 not resolvable)
- SBDB NEOs missing from CSV: **15,489** (36.3% of all NEOs); of these 422 are flagged PHA
- **The CSV does NOT contain all asteroids**: it covers 63.7% of today's NEO catalogue (it is a ~2014-2022 close-approach set; later discoveries and objects without a listed approach are absent).

## Matching method (CSV id -> SBDB)

CSV ids are old-style SPK ids (2,000,000+number) while SBDB uses 20,000,000+number, so direct spkid matching only hits some. Order tried:

- designation/pdes: 17,557
- spkid: 6,259
- number/pdes: 2,727
- designation/full_name: 799
- unmatched: 44
- SBDB alias lookup: 37

## Element table (`data/neo_elements.json`)

- SBDB NEOs in: 42,670
- Dropped: 0 
- Kept: **42,670** (27,181 in CSV, 15,489 not in CSV)
- PHA flag null in SBDB (treated as false): 134
- Diameter source: sbdb=1,245, H_albedo0.14=41,419, none=6
- Units: a in AU, angles in degrees, epoch as JD (TDB), per in days, moid in AU, diameter_km in km.

## CSV objects unmatched (first 30)

- 54250323 (2022 CL12) (SBDB neo flag: False, SBDB: (2022 CL12))
- 3076637 (2001 DC77) (SBDB neo flag: False, SBDB: (2001 DC77))
- 54249331 (2022 DD5) (SBDB neo flag: False, SBDB: (2022 DD5))
- 3137774 (2002 EH1) (SBDB neo flag: False, SBDB: (2002 EH1))
- 54054450 (2020 SO) (SBDB neo flag: None, SBDB: None)
- 2518426 518426 (2002 ON4) (SBDB neo flag: None, SBDB: None)
- 54253551 (2022 EQ4) (SBDB neo flag: False, SBDB: (2022 EQ4))
- 54138703 (2021 HS) (SBDB neo flag: True, SBDB: P/2021 HS (PANSTARRS))
- 3297355 (2005 UK5) (SBDB neo flag: False, SBDB: (2005 UK5))
- 3291226 (2005 SR9) (SBDB neo flag: False, SBDB: (2005 SR9))
- 3321533 (2006 CV9) (SBDB neo flag: False, SBDB: 849638 (2006 CV9))
- 3358226 (2006 VB3) (SBDB neo flag: False, SBDB: (2006 VB3))
- 3367152 (2007 BB50) (SBDB neo flag: False, SBDB: (2007 BB50))
- 3465674 (2009 PN) (SBDB neo flag: False, SBDB: (2009 PN))
- 2217683 217683 (1999 RP36) (SBDB neo flag: None, SBDB: None)
- 3515378 (2010 GJ23) (SBDB neo flag: False, SBDB: (2010 GJ23))
- 3512280 (2010 FP) (SBDB neo flag: False, SBDB: (2010 FP))
- 3020949 (1999 HW1) (SBDB neo flag: False, SBDB: (1999 HW1))
- 3720880 (2015 JF11) (SBDB neo flag: False, SBDB: (2015 JF11))
- 3588058 (2011 WK2) (SBDB neo flag: False, SBDB: (2011 WK2))
- 3602598 (2012 GM11) (SBDB neo flag: False, SBDB: (2012 GM11))
- 54106000 (2021 AD7) (SBDB neo flag: False, SBDB: (2021 AD7))
- 2038091 38091 (1999 JT3) (SBDB neo flag: None, SBDB: None)
- 3612840 (2012 UW27) (SBDB neo flag: False, SBDB: (2012 UW27))
- 3824112 (2018 HT3) (SBDB neo flag: True, SBDB: 463P/NEOWISE (2018 HT3))
- 3838007 (2019 AY14) (SBDB neo flag: False, SBDB: (2019 AY14))
- 3754315 (2016 LF51) (SBDB neo flag: False, SBDB: (2016 LF51))
- 3785710 (2017 SD33) (SBDB neo flag: False, SBDB: (2017 SL18 = 2017 SD33))
- 3791301 (2017 XV2) (SBDB neo flag: False, SBDB: (2017 XV2))
- 3802101 (2018 FH5) (SBDB neo flag: False, SBDB: (2018 FH5))
