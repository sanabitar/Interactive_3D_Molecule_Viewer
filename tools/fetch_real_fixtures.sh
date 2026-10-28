#!/bin/sh
# Run where network access exists: downloads real fixtures listed in tests/fixtures/manifest.json (RCSB data are CC0 1.0).
cd "$(dirname "$0")/../tests/fixtures" && mkdir -p real && for id in 1CRN 4HHB 1UBQ 1HSG 2HHB; do curl -fsSL -o real/$id.pdb https://files.rcsb.org/download/$id.pdb; curl -fsSL -o real/$id.cif https://files.rcsb.org/download/$id.cif; done; ls -l real
