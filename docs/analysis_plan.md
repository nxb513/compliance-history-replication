# Phase 4 Extension: Pre-specified Analysis Plan
_2026-09-24. Written BEFORE running any extension analysis, so that results cannot steer which analyses are reported. Results of an earlier water-only analysis are not modified. Every analysis below is descriptive; no causal claims. All analyses use the primary sample (924 private firms with at least 3 plants; 7,215 plants; 55,560 plant-years; 2010-2025) unless stated, and the same firm-clustered bootstrap (resample parent firms). All results will be reported whether or not they support the paper's thesis._

## Motivation
An earlier water-only analysis shows that persistent effluent noncompliance is mainly plant-specific. Reviewers will ask: (a) is the plant component really a plant attribute, or a regulator/state effect, a permit-stringency effect, or an artifact of counting violations? (b) what kind of plant-specific problem is it? (c) how does it relate to pollution volumes and to other environmental programs? The extensions below use data already on disk that the paper has not used.

## E1. Technical versus procedural obligations
- **Data:** QNCR counts of effluent violations (E90, technical) and DMR reporting violations (D80/D90, procedural), already in the panel (`e90`, `rep`).
- **Question:** Is the organizational level the same for technical compliance (meeting effluent limits) and procedural compliance (filing reports on time)?
- **Expectation from the two pathways:** procedural compliance depends on administrative systems that firms often standardize, so its firm share should be larger than for effluent violations. The opposite result is also informative.
- **Test:** nested partition and REML for reporting outcomes; difference in D (plant minus firm) between effluent and reporting outcomes, firm-clustered bootstrap.

## E2. Regulator (state) versus firm versus plant
- **Data:** plant state (`STATE_CODE`), industry (NAICS), firm.
- **Question:** Does the plant component partly reflect the state regulator (states administer NPDES) or the industry?
- **Tests:** (a) four-level nested partition firm > firm-within-state > plant > year; the firm-within-state level captures a firm's operations under one state regulator. (b) Sibling concordance for pairs of plants of the same firm that are in the same state versus different states, and in the same versus different NAICS-4 industry. (c) Share of plant-level variance explained by state and NAICS-4 fixed effects.
- **Reading:** a small firm-within-state component and small same-state premium mean the plant component is not a regulator effect.

## E3. Pollution volume versus legal noncompliance (same plants, same medium)
- **Data:** TRI "5.3 - WATER" (on-site releases to water) and "ON-SITE RELEASE TOTAL", pounds, reliable years 2010 and 2017-2024, summed per plant-year.
- **Question:** Do water releases and water violations have the same organizational structure? Are chronic violators the same plants as heavy dischargers?
- **Tests:** nested partition and REML of log(1 + water releases) for plants in the primary sample that report to TRI; compare firm share with that of violations for the same plants; plant-level correlation of long-run releases and long-run violations; overlap of the top decile of each.
- **Expectation (Galli Robertson and Collins 2019):** pollution volumes carry a larger parent-level component than violations. Either result is reported.

## E4. Pollutant-level anatomy of chronic noncompliance
- **Data:** `NPDES_EFF_VIOLATIONS.csv` (15.6 GB), filtered to E90 effluent violations 2010-2025 at the NPDES permits of primary-sample plants.
- **Questions and tests:** (a) concentration of a plant's violations in its top pollutant parameter; (b) persistence at the pollutant level (is the same parameter violated year after year?); (c) do sibling plants violate the same parameters (Jaccard overlap of violated-parameter sets) more than plants of different firms in the same NAICS-4 industry and state? (d) severity measured by exceedance magnitude (EXCEEDENCE_PCT), decomposed into plant and firm components. Item (d) addresses the unresolved severity question with a better measure than counts of "serious" violations.
- **Reading:** concentration in one or two parameters that persist, with little sibling overlap beyond industry, points to plant-specific treatment or process problems.

## E5. What explains the plant component?
- **Data:** permit limits (`NPDES_LIMITS.csv`: number of distinct enforceable limited parameters and outfalls per permit), permit components (stormwater, pretreatment), design flow and original permit issue date (facility master), major status, NAICS, state, TRI reporting.
- **Test:** regress each plant's long-run noncompliance (mean any-E90, and REML plant effect) on observables, with and without firm fixed effects; report the share of between-plant variance explained, especially within firms.
- **Reading:** a large unexplained share means the plant component reflects unobserved operational factors; a large permit-complexity share means regulatory burden drives it.

## E6. Cross-program enforcement
- **Data:** ICIS-FEC federal enforcement cases (`CASE_FACILITIES`, `CASE_PROGRAMS`) 2010-2025, programs classified as CWA, CAA, RCRA, other.
- **Question:** Is chronic CWA noncompliance a plant trait that extends to other programs, or a firm trait?
- **Test:** probability of a CAA or RCRA federal case at a chronically CWA-noncompliant plant versus at its sibling plants versus at other plants; firm-clustered intervals.
- **Limitation noted in advance:** federal cases only; state enforcement and non-enforced violations are not observed. Full cross-media compliance histories would require an additional EPA download (ECHO Exporter / ICIS-Air / RCRAInfo), which requires the author's approval.

## Reporting rules
- Report every pre-specified test. Label any additional analysis as exploratory.
- Descriptive language only ("associated with", "accounts for").
- Outputs: `extension/tables/`, `extension/figures/`, `extension/extension_report.md`.

---
# Addendum (2026-09-25): cross-program analyses with newly downloaded data
_Written after the author approved the download and BEFORE opening the new files. Data: ICIS-AIR_downloads.zip (Clean Air Act), rcra_downloads.zip (hazardous waste), echo_exporter.zip (facility-level multi-program summary), EPA ECHO. Same primary sample, methods, bootstrap, and reporting rules as above._

## E7. Organizational structure of noncompliance in other programs
- **Measures:** for primary-sample plants regulated under the Clean Air Act and/or RCRA, annual indicators 2010-2025 of any air violation (ICIS-Air violation history) and any RCRA violation (RCRAInfo violations), plus high priority violations (air) where available.
- **Test:** nested partition and REML (firm / plant / year) for each program; paired firm bootstrap comparing D with water noncompliance on the same plants.
- **Expectation:** if chronic noncompliance is mainly a plant attribute, air and RCRA noncompliance are also plant-dominated.

## E8. Cross-program covariance: plant or firm?
- **Test:** plants' long-run violation rates in water, air, and RCRA. Compare (a) the correlation between two programs at the same plant, (b) the correlation between program A at a plant and program B at its sibling plants (same state and different state), and (c) the same cross-program correlation for unrelated plants in the same state.
- **Expectation:** (a) exceeds (b), and (b) is close to (c), meaning the cross-program propensity sits in the plant.

## E9. Replication of E6 with full compliance histories
- **Test:** rates of air violations, high priority air violations, and RCRA violations for chronic water violators, their non-chronic siblings, and plants in firms without a chronic plant; linear probability models with own and sibling chronic status, with and without firm fixed effects; firm-clustered errors.
- **Reading:** unlike E6 (federal cases only), this uses violations recorded by state and federal regulators.

## E10. Cross-program regulatory attention
- **Test:** air compliance evaluations and RCRA evaluations per plant, 2010-2025, by the same three groups. Descriptive; higher attention can itself raise detected violations, which is noted when reading E7 to E9.

## E11 (exploratory, labeled as such). Community context within firms
- **Data:** demographic fields in the ECHO Exporter, if present (for example percent minority or low-income population near the facility).
- **Test:** within-firm comparison (firm fixed effects) of the community characteristics of a firm's chronically noncompliant plants versus its other plants.
- **Reading:** descriptive association only; reported regardless of direction.

### Measurement note added after inspecting file structure (before any outcome was computed)
ICIS-Air records federally reportable violations (FRV) and high priority violations (HPV). Reporting is required mainly for major and synthetic-minor sources, so air violations at minor sources are under-recorded. Air analyses are therefore reported (i) for all primary-sample plants with an ICIS-Air ID and (ii) for plants whose air pollutant class is major or synthetic minor. Air and RCRA outcomes are measured on the same plant-years as the water panel, so programs are compared on identical plant-years. Violation year = earliest FRV determination date (FRV) or HPV day-zero date (HPV); RCRA violation year = date the violation was determined.

### RCRA measurement note (added after inspecting file structure, before computing outcomes)
RCRA violations are recorded when regulators inspect, so they measure detected violations and depend on inspection effort; this is noted in reading E7 to E10 for RCRA. Plants are linked to RCRAInfo handler IDs through FRS (program acronym RCRAINFO). Outcomes per plant-year: any violation determined in the year; any month flagged as significant noncomplier (SNC) in RCRA_VIOSNC_HISTORY; number of evaluations. Results are reported for all linked plants and for large quantity generators or operating treatment, storage, and disposal facilities, which face the most regular inspection.

### E11 measurement note (after inspecting ECHO Exporter columns, before computing outcomes)
The ECHO Exporter has one demographic field that is populated for the sample: FAC_PERCENT_MINORITY, which ECHO defines as percent people of color within 3 miles of the facility (ACS 5-year summary; ECHO search-results help). FAC_POP_DEN is empty for all sample plants. Contextual flags used as secondary E11 outcomes: FAC_NAA_FLAG (facility in an air nonattainment area) and FAC_IMP_WATER_FLG (impaired-water flag; populated as Y or blank). Income and poverty fields exist only in the separate ECHO Facility Demographic download (ECHO_DEMOGRAPHICS.csv), which has not been downloaded.

### E11 extension with the ECHO demographic download (written after the author approved the download on 2026-09-25 and BEFORE opening the file; exploratory, like E11)
- **Data:** `echo_demographics.zip` (ECHO_DEMOGRAPHICS.csv; ACS 2019-2023; 1, 3, 5 mile radii; keyed by REGISTRY_ID). Field definitions confirmed from the ECHO Exporter data dictionary (`echo_exporter_columns_7-16-2025_0.xlsx`): FAC_IMP_WATER_FLG = facility discharges into a water identified as impaired (category 4 or 5 in WATERS; the dictionary lists values 4/5, the 2026 file codes Y or blank).
- **Outcomes (3-mile radius is primary, matching FAC_PERCENT_MINORITY; 1 and 5 miles reported as sensitivity):** share low-income (income below twice the poverty level), share below the poverty level, share people of color (from the same file).
- **Tests:** same as E11: group means (chronic, sibling, other firms); pooled own/sibling LPM; own chronic with firm fixed effects; with firm and state fixed effects; long-run water violation rate with firm fixed effects; firm-clustered SE. Reported whatever the direction.

### Exploratory checks added AFTER seeing the E7 to E10 results (labeled exploratory in all reporting)
- X-a. Detection check for E9: repeat the own/sibling LPMs controlling for the plant's own evaluation intensity in the program (log of 1 + evaluations per year), because E10 shows that chronic water violators are evaluated more often.
- X-b. Worst-plant coincidence: in firms with at least 3 plants observed in both water and another program, how often the plant with the highest water violation rate is also a plant with the highest violation rate in the other program, compared with random assignment within firm (1,000 within-firm permutations of the other program's rates).
- X-c. Enforcement-history descriptives from the ECHO Exporter (5-year windows): formal actions, penalties, inspections, federal cases across all programs for chronic plants, their siblings, and plants in firms without a chronic plant.

---
# Addendum 2 (2026-09-25): further analyses with data already on disk (E12 to E14)
_Written BEFORE running any of these analyses. No new download. Same primary sample, reporting rules, firm-clustered inference, descriptive language. REML fits use the best log-likelihood across lbfgs, Powell and Nelder-Mead (see the REML optimizer check)._

## E12. Pollution volumes versus violations in the air medium
- **Data:** TRI "5.1 - FUGITIVE AIR" plus "5.2 - STACK AIR" (pounds; grams converted), 2010 and 2017-2024, summed per plant-year; air panel (`air_panel.parquet`).
- **Question:** Does the E3 contrast (volumes carry a larger firm component than violations) also hold for air?
- **Tests:** as in E3, on plants in both the air panel and TRI (at least 2 TRI years; firms with at least 2 such plants): nested partition and REML of log(1 + TRI air releases) and of any air violation on the same plant-years; firm bootstrap of the firm-share difference; plant-level Spearman correlation of long-run air releases and long-run air violation rate; top-decile overlap.
- **Expectation:** air releases have a larger firm share than air violations. Either result is reported.

## E13. Stable trait or shared bad years? Decomposing cross-media co-occurrence
- **Data:** water panel with air and RCRA panels (same plant-years).
- **Tests:** (a) split the plant-year covariance between water and each other program into a between-plant part (plant means) and a within-plant part (deviations from plant and year means), and report the within-plant correlation of annual deviations; (b) lead-lag linear probability models with plant and year fixed effects: water violation in t+1 on air (RCRA) violation in t, and the reverse; firm-clustered SE.
- **Reading:** a small within-plant correlation and weak lead-lag relations mean the cross-media link is a stable plant attribute rather than shared temporary shocks; a sizable within-plant correlation means plants also have bad years across media at the same time.

## E14. Does formal enforcement at one plant coincide with changes at its sister plants? (descriptive event study)
- **Data:** dated formal NPDES enforcement actions (`data/interim/formal_enf.parquet`, SETTLEMENT_ENTERED_DATE year, EPA and state), mapped to plants through NPDES IDs; water panel.
- **Event:** the first year in which any plant of a firm receives a formal NPDES action, for firms with no formal action at any plant in the previous 3 years; event years 2013 to 2022 (3 years before and after are observed).
- **Units:** the firm's other plants (siblings of the enforced plant). Controls: plants of firms with no formal action from 3 years before to 3 years after that event year (clean controls). Stacked by event year; window -3 to +3; reference year -1.
- **Specification:** any effluent violation on event-time-by-treated indicators, with plant-by-stack and year-by-stack fixed effects; firm-clustered SE. Also reported: the enforced plant's own path; siblings in the same state versus other states (pre-specified because E2 and E8 show a local firm signal).
- **Pre-trend check:** joint test that the lead coefficients (-3, -2) are zero. If it fails, the post-event coefficients are not interpreted.
- **Reading (stated in advance):** enforcement timing responds to violations, and violations at a firm's plants in the same state are correlated, so any change at siblings is an association, not a deterrence effect. The analysis is reported whatever it shows.

### Robustness checks added during the full audit (2026-09-25; labeled "robustness" in reporting)
- R1. Firm-bootstrap confidence intervals (200 draws, firms resampled; pairs weighted by resampled firm counts) for the E8 correlations and for two contrasts: same plant minus same-state sibling, and different-state sibling minus unrelated same-state plant.
- R2. E9 with alternative chronic definitions: effluent violations in at least 4 years, at least 8 years, and top decile of the long-run water violation rate.
- R3. Linkage coverage: water outcomes of primary plants with and without an ICIS-Air or RCRAInfo link.

---
# Addendum 3 (2026-09-25): B1 to B5, written BEFORE running (author approved running them)
_Data on disk only. Same primary sample, firm-clustered inference, descriptive language. All results reported._

## B1. Information value for targeting: whose history predicts a plant's future violations?
- **Split:** history window 2010-2017; outcome window 2018-2025. Plants observed at least 3 years in each window, in firms with at least 2 such plants in the program sample.
- **Outcomes (plant level, 2018-2025):** any air violation; any high priority air violation; any RCRA violation; any RCRA significant noncomplier month; any effluent (water) violation.
- **Information sets (2010-2017):** (1) the plant's own water violation rate; (2) the plant's own history in the outcome program (air or RCRA violation rate; for water, the same as 1); (3) mean water violation rate of the plant's same-state siblings; (4) of its other-state siblings; (5) of all its siblings (firm history, leave the plant out); (6) the same sibling and firm means of the program history.
- **Metrics:** AUC of each single information set (ranking plants; no fitting, so the evaluation is out of sample in time); top-decile precision (share of the top 10% by the predictor that violate in 2018-2025) against the base rate. Combinations (own program history + own water; + firm history): logistic models fitted on a random half of firms and evaluated on the other half, 100 splits, mean test AUC. Firm-bootstrap (200 draws) CIs for the AUC differences own water minus firm water and own water minus same-state siblings' water.
- **Reading:** if own history beats sibling and firm history, plant-level records carry most of the targeting information; incremental AUC from firm history shows whether firm records add anything.

## B2. Detection source
- Air violations split by the agency recorded in ICIS-Air (state or local vs U.S. EPA); RCRA violations by VIOL_DETERMINED_BY_AGENCY (S = state, E = EPA). Re-estimate the E9 own/sibling LPM (with and without firm FE) for "ever a state- or local-detected violation" and "ever an EPA-detected violation", 2010-2025.

## B3. Heterogeneity of the cross-media gap
- E9 LPM with interactions of own chronic status with: NPDES major; NAICS-3 sector group from FRS (manufacturing 31-33, utilities 221, mining and extraction 21, other or unknown); firm size (number of plants in the primary sample: 3-5, 6-15, 16 or more). Report group-specific own and sibling coefficients.

## B4. Stability over time
- Separately for 2010-2017 and 2018-2025: chronic defined within the period (effluent violations in at least 3 of the period's years); E9 own/sibling LPM with outcomes measured in the same period; E8 correlations (same plant, same-state sibling, other-state sibling, unrelated same state) of within-period long-run rates.

## B5. Multiple testing
- Holm and Benjamini-Hochberg adjusted p-values within two families: (F1) all own-chronic coefficients in E9, R2, B2 (pooled specs); (F2) all sibling coefficients in the same tables. Plus the E14 post-event joint tests as a third family.

---
# Addendum 4 (2026-09-25): Extension, written BEFORE running (author: extend and tighten; planning note redacted)
_Data already on disk (ICIS-Air formal actions and stack tests extracted from the approved ICIS-Air download; ECHO Exporter coordinates; FRS NAICS; NPDES formal actions). Same primary sample and inference rules. Each test below is reported whatever it shows._

## E15. What is the "local" firm component? Distance versus state border (sibling pairs)
- **Idea:** E2, E8 and R1 show that sister plants in the same state are much more alike than sister plants in different states. Two explanations predict different patterns: (a) a shared state regulator (then the same-state premium should appear for unrelated plants too, and should hold at any distance within a state); (b) shared local management, staff or resources (then similarity should fall with physical distance whether or not the plants are in the same state, and nearby cross-border siblings should look like nearby same-state siblings).
- **Data:** plant coordinates (ECHO Exporter FAC_LAT, FAC_LONG); outcome = plant's long-run any-effluent-violation rate 2010-2025, standardized; secondary outcome = the same rate residualized on plant attributes and NAICS-3 (not on state), as in E5b.
- **Pairs:** all sibling pairs; unrelated pairs = random sample of pairs of plants of different firms (sampled to cover the distance range; 1,000,000 draws with oversampling of pairs within 400 km, weights restore population shares).
- **Cells:** sibling vs unrelated; same state vs different state; distance bins 0-50, 50-100, 100-200, 200-400, 400-800, over 800 km.
- **Statistics:** pairwise correlation in each cell; pair-level regression of the product of standardized outcomes on sibling, same state, sibling x same state, distance-bin dummies, sibling x distance-bin dummies, and same NAICS-3; firm bootstrap (200 draws; unrelated-pair weights = product of the two firms' draw counts).
- **Pre-specified contrasts:** (T1) state border at fixed distance: same-state minus cross-state correlation for pairs within 200 km, separately for siblings and unrelated plants; (T2) firm premium: sibling minus unrelated, by state status and distance; (T3) distance decay of the sibling premium: within 100 km vs beyond 400 km.
- **Reading:** (a) supported if the border contrast is similar for siblings and unrelated plants and the sibling premium does not depend on distance; (b) supported if the sibling premium falls with distance and nearby cross-border siblings resemble nearby same-state siblings; a sibling premium only within the same state at every distance points to the firm's relationship with a particular state regulator.
- **Secondary:** the same analysis for the water-air cross-program correlation.

## E16. Targeting: budget-constrained allocation, and whether regulators already use cross-program information
- **E16a (simulation):** B1 set-up (history 2010-2017, outcomes 2018-2025, 100 random half-splits of firms). For inspection budgets of 5%, 10% and 20% of plants, select plants by predicted risk from each information set (NPDES major only; own program history; + own water history; + sister plants' histories; sister plants' histories only) and report the share of all 2018-2025 violators captured (recall) and the violation rate among selected plants (precision), for air violations, high priority air violations, RCRA violations and RCRA significant noncompliers.
- **E16b (actual allocation):** plant-level regressions of 2018-2025 evaluations per year (air FCE/PCE; RCRA evaluations) on the plant's 2010-2017 own-program violation rate, its 2010-2017 water violation rate, sister plants' histories, NPDES major and program class, firm-clustered SE; compared with the same regression for 2018-2025 violations. If water history predicts future air/RCRA violations but not future air/RCRA evaluations, cross-program information is not used in allocating evaluations.

## E17. Intra-firm enforcement spillovers: direct comparison with Li & Lyon (2026)
- **Design as in Li & Lyon:** facility-year panel 2010-2025 (air panel; and water panel for NPDES), outcome = any violation in the program in year t; regressors = indicators for a formal action with a penalty in t-1 at the plant itself and at sister plants in four categories (same/different industry x same/different state; industry = shared FRS NAICS-6 code; robustness NAICS-4); plant and year fixed effects; firm-clustered SE. Air penalties from ICIS-AIR_FORMAL_ACTIONS (PENALTY_AMOUNT > 0); water penalties from NPDES formal actions (federal or state penalty > 0).
- **Also:** the same regression with any formal action (with or without penalty); and sibling penalties in the other program (a water penalty at a sibling and air violations at the plant).
- **Reading:** Li & Lyon report +0.017 to +0.023 for same-industry-same-state siblings (air, 2005-2017). A similar estimate here would replicate them; a null would indicate that the result depends on period, sample or definitions.

## E18. Compliance measured independently of inspector discretion: air stack tests
- **Outcome:** any failed stack test in the year (ICIS-AIR_STACK_TESTS, status Fail vs Pass; pending and incomplete excluded), on air-panel plant-years with at least one completed test.
- **Tests:** E7 variance partition (nested and REML) and E9 own/sibling LPM (ever failed a stack test 2010-2025, among plants with tests; controls include log number of tests).
- **Reading:** stack tests are scheduled measurements with a pass/fail result; if chronic water violators also fail stack tests more often, the cross-media pattern is not an artifact of where inspectors look.

## E19. Firm-wide settlements (feasibility first)
- Count ICIS-FEC cases that name two or more primary-sample plants of the same firm (2010-2025). If at least 30 such cases exist, run a stacked event study of water violations at covered plants and at uncovered sister plants versus clean controls, as in E14. Otherwise report the count and stop.

---
# Addendum 5 (2026-09-25): regulatory regime changes inside the sample period, written BEFORE running
_Author's instruction: check the effective dates of legal texts and amendments inside the data period. Two amendments change how violations are recorded: (1) the revised CAA High Priority Violation Policy (Aug. 25, 2014; implemented from Oct. 1, 2014, i.e., FY2015; criteria reduced from 10 to 6) replaced the 1998 HPV Policy; (2) the NPDES Electronic Reporting Rule (80 FR, Oct. 22, 2015; Phase 1 electronic DMRs by Dec. 21, 2016) changed DMR data flows into ICIS-NPDES. Compliance monitoring strategies were also revised (CAA CMS 2014 and Oct. 2016; NPDES CMS July 2014 replacing 2007; RCRA CMS 2015 and Dec. 2021), but recommended baseline frequencies remained program-specific throughout._

## E20. Results within regulatory regimes
- **Air regimes:** 2010-2014 (1998 HPV Policy) vs 2015-2025 (2014 HPV Policy).
- **Water regimes:** 2010-2016 (before electronic DMR reporting) vs 2017-2025 (after Phase 1).
- **Within each regime:** chronic water violator = effluent violations in at least 40% of the plant's observed years in the regime (at least 2 observed years); sibling = another chronic plant in the firm.
- **Tests in each regime:** (a) E9 own/sibling LPM for ever an air violation and ever an HPV (air regimes), ever an RCRA violation (water regimes), controls major and program class, firm-clustered SE; (b) water REML firm/plant shares (best of three optimizers) in each water regime; (c) same-state and different-state sibling correlations of within-regime water violation rates.
- **Reading:** the conclusions hold across regimes if the own-plant gap is large and the sibling coefficient near zero in each regime, and the plant share exceeds the firm share in each water regime.
- *Note added after E20 was run:* the description "criteria reduced from 10 to 6" came from a search summary and is not verified; verified: the revised HPV Policy took effect Oct. 1, 2014 and replaced the 1998 policy's violation matrix with six HPV categories. The regime split (2010-2014 vs 2015-2025) is unaffected.

---
# Addendum 6 (2026-09-25): a causal test using the Florida enforcement shock, written BEFORE running
_Shock: after EPA's 2012 review of the Florida DEP, penalties for Clean Air Act high priority violations in Florida rose sharply, with little change for other plants and no change in neighbouring states (Blundell 2020, JEEM 101, 102288; control states Alabama, South Carolina, Georgia; effects persist into 2017). Window 2010-2017; 2012 = transition year; reference year 2011; post = 2013-2017. Pre-trend check: the 2010 coefficient._

## E22a. Validation (replicates Blundell in this sample)
Air violations at air-regulated primary-sample plants in Florida vs Alabama, South Carolina, Georgia; plant and year fixed effects; event-study coefficients by year and a pooled post (2013-2017) coefficient; firm-clustered SE. Inference check with few treated states: placebo distribution from re-running the same regression with each other state (at least 30 air-regulated plants) as the "treated" state against the same controls; report FL's rank.

## E22b. Cross-program effect at the same plant (main test 1)
Same sample and design, outcome = any effluent violation (self-reported DMR data, so not affected by what air inspectors detect). If the air enforcement shock lowers water violations at the same Florida plants, compliance responds as a plant-level, cross-program decision.

## E22c. Transmission through the firm (main test 2)
Plants outside Florida in the water panel; treated = the plant's firm had an air-regulated Florida plant in 2010-2012; controls = plants outside Florida of other firms; plant fixed effects and state-by-year fixed effects (plants in the same state and year are compared); outcome = any effluent violation; secondary outcome = any air violation (air-panel plants outside Florida). Placebo: plants outside FL, AL, SC and GA whose firm has a plant in AL, SC or GA but none in Florida, treated as if exposed. If out-of-state sister plants do not respond, the firm is not the channel through which the Florida shock travels.

## Reporting
All coefficients reported whatever their sign; language: effects of the Florida enforcement change, conditional on the identifying assumption of parallel trends, supported or not by the pre-period coefficient and the placebos.

---
# Addendum 7 (2026-09-25): natural-disaster shocks as exogenous plant-level shocks (E23), written BEFORE opening the storm data
_Author approved downloading NOAA NCEI Storm Events detail files 2009-2025 (17 csv.gz, about 193 MB). Aim: one credible causal piece. Floods are not chosen by plants or firms, so a major flood in a plant's county is an exogenous shock to that plant._

## Exposure
- Events: EVENT_TYPE "Flood" or "Flash Flood" recorded by county (CZ_TYPE = "C"), matched to the plant's county (ECHO Exporter FAC_DERIVED_STCTY_FIPS). Zone-based events (hurricane, storm surge, recorded by NWS zone) are not used because a zone-to-county crosswalk would be needed; this is stated as a limitation.
- County-year property damage = sum of DAMAGE_PROPERTY (K, M, B suffixes converted to dollars, nominal).
- Major flood: county-year property damage of at least $10 million (sensitivity: $2 million and $50 million).
- Event for a plant: the first year e in 2011-2022 with a major flood in its county and no major flood in its county in e-3 to e-1.

## Groups (stacked by event year e, window e-3 to e+3, reference e-1)
- Hit plants: plants with an event in year e.
- Unaffected sister plants: plants of a firm with a hit plant in year e, located in a different state from every hit plant of that firm in year e, and whose own county had no flood damage of $1 million or more in e-3 to e+3.
- Controls: plants of firms with no hit plant in e-3 to e+3, whose county had no flood damage of $1 million or more in e-3 to e+3.

## Outcomes and tests
- Hit plants (first stage and cross-program test): any effluent violation (self-reported DMRs; the first-stage check), any air violation, any RCRA violation.
- Unaffected sister plants (firm-channel test): the same three outcomes.
- Specification: event-time x group indicators, plant-by-stack and year-by-stack fixed effects, firm-clustered SE; pre-trend test (leads -3, -2 jointly zero); post effect = average of e to e+2.
- Reading, stated in advance: if hit plants' water violations rise (the shock bites) and their air and RCRA violations also rise, compliance capacity is shared across programs within the plant. If unaffected sister plants do not change, the firm does not transmit (or reallocate) compliance problems across plants; a rise would indicate reallocation (the Li & Lyon mechanism). If the first stage fails, the sibling results are not interpreted. After disasters, regulators may suspend or relax inspections, which affects inspection-detected air and RCRA violations at hit plants; this is why the water (DMR) outcome is the first-stage check and why the sibling test (plants in other states, normal regulation) is the main causal test.


---
# Addendum 8 (2026-09-25): E21, is the corporate-history presumption applied in penalty practice? Written BEFORE running
_Motivation: the penalty policies (GM-22 1984; CAA Stationary Source Civil Penalty Policy 1991; RCRA Civil Penalty Policy 2003; CWA Settlement Penalty Policy 1995) instruct staff to count prior violations at other facilities of the same corporation as history of noncompliance, and the statutes list the violator's history as a penalty factor. E21 asks whether assessed penalties in the data move with the history of the firm's other plants, holding the plant's own history fixed. This compares the policy on paper with the policy in practice. Descriptive association only._

## Sample
- Formal enforcement actions with a settlement or action date in 2013-2025 at primary-sample plants (firms with at least 3 plants), by program: water (NPDES formal actions, federal plus state/local penalty), air (ICIS-Air formal actions, PENALTY_AMOUNT), RCRA (formal enforcement types: numeric part of the RCRAInfo type code 200-399 or 500-699, i.e. initial and final administrative orders and civil judicial actions; informal 100-199 and referrals 400-499 excluded). RCRA penalty of a unit = sum of final monetary penalties (FMP_AMOUNT) of its formal actions; only if none of them has an FMP, the sum of proposed penalties (PMP_AMOUNT), to avoid counting the proposed and the final amount of one case twice.
- Unit: plant x program x year with at least one formal action (amounts summed within the unit). Start year 2013 so that every unit has a full three-year history window inside 2010-2025.

## Variables (history window t-3 to t-1, t = action year)
- Own, same program: share of the three years with a violation of that program at the plant (water: any effluent violation; air: any FRV/HPV; RCRA: any violation).
- Own, other programs: indicator that the plant had a violation in any of the other two programs in the window (0 for programs the plant is not in; an indicator for being in at least one other program is included).
- Sister plants, same state: share of the firm's other plants in the same state with a violation in any program in the window; plus an indicator for having any same-state sister plant (share set to 0 when none).
- Sister plants, other states: same for the firm's plants in other states.
- Controls: year fixed effects, state fixed effects, lead agency (EPA vs state), judicial vs administrative action, log number of plants of the firm, program class (NPDES major; air major/synthetic minor; RCRA LQG or TSDF).

## Outcomes and tests
- (a) Any penalty > 0 (linear probability model), all units.
- (b) log(total penalty) among units with a penalty > 0.
- By program (water, air, RCRA) and pooled with program fixed effects. Firm-clustered SE.
- Coefficients of interest: the two sister-plant shares; for comparison, own same-program history and own other-program history.
- Multiple testing: Holm and Benjamini-Hochberg over the 12 sister-plant coefficients (3 programs x 2 outcomes x 2 sister groups).

## Reading, stated in advance
- Sister-plant history positive and material (relative to the own-history coefficient) and surviving the multiple-testing adjustment: the corporate-history presumption is applied in practice. Combined with L4 (sister histories carry little predictive information), penalties would then load on information with little value for predicting future violations.
- Sister-plant coefficients near zero: the presumption on paper does not show up in assessed penalties; penalties follow the plant's own record. This is a gap between the written policy and practice, and practice would then be aligned with where the information is.
- Either result is reported. Caveats stated in advance: only settled formal actions are observed (selection on reaching a formal action); state penalty amounts in ICIS and RCRAInfo may be incomplete; amounts reflect negotiated settlements and ability to pay, not only history.

---
# Addendum 9 (2026-09-25): R4, parent-linkage robustness of the plant/firm decomposition. Written BEFORE running
_Reason: parent firms come from self-reported TRI parent names, fixed over time. Misassigned plants would bias the firm share of variance toward zero, which is the paper's central quantity. Readiness check before writing._
- Sample: the "DUNS-confirmed" subset defined in the original design (docs/one_bad_plant_formal.md, sample C): firms with at least 3 plants whose plants all carry the same TRI parent D&B number (427 firms, 2,445 plants in the water panel). No new cut-offs are introduced.
- Models: the same REML variance decomposition (best of three optimizers, `reml_best`) as E7 for water (any effluent violation), air (any violation) and RCRA (any violation), on the DUNS-confirmed subset; primary-sample estimates shown alongside.
- Reading, stated in advance: if the firm share stays well below the plant share in all three programs, L1 is robust to linkage error. If the firm share rises materially (for example to the level of the plant share), the paper reports it and states the firm component as a lower bound.
