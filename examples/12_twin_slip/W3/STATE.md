# W3 milestone / closeout

2026-10-29, branch `codex/w3-source-contact`; W1/W2 preserved; no main merge.

- A0 PASS: five configuration/geometry contracts, ten native response checks including three W1 native pin contracts, upstream 12 unit + 2 FR cases. Frozen A0 interface files unchanged.
- A1 PASS in effective FCC benchmark: 240/80 nm source, 200 MPa pure HI RSS, 16 um nonperiodic cell, three independently growing detached free components at 1.082511 ns. Parent retained; pin/Burgers error zero. Native 20 nm first emission agrees within 2.22 percent; 3-emission spatial convergence is not claimed.
- B0 PASS at 300 MPa: existing 60 nm loops expand without walls in HI/HII. The 200 MPa shrink/annihilation counterexample is retained.
- Minimal B1 PASS: 16 same-initial-array cases to common 0.1 ns; all integration/collision/topology/remesh contact states admissible; native tangent/release tests; four temporal refinements and two actual reversed-load releases passed.
- Second wall reduces cumulative sweep relative to lower-only by 59.17/51.80 percent (HI/HII loop) and 67.10/60.54 percent (prepared threader). Two-wall nonadditivity is small and not claimed as qualified material strengthening.
- Raw incorrect counters, narrow-domain trajectories, empty-network failures and two wall-time-limited 0.25 ns pilots remain in evidence. Formal verdict uses only documented accepted cases.

Unresolved / NOT_RUN: VCoNi friction/mobility/interface calibration; actual CTB transmission/absorption; source operation with TB; both crystallographic sides of one TB; periodic multi-lamella contact and FFT/integrator transfer; full local plastic tensor field; four-thickness collective strengthening/localization; full spatial/topological convergence of small nonadditivity.

All native dynamic runners stopped at closeout. Runtime memory registration unavailable; no agent-attributed memory writes performed.
