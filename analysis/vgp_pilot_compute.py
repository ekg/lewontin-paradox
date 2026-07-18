#!/usr/bin/env python3
"""Scientific compute element for an already-authorized repaired VGP row.

This entrypoint is intentionally fail-closed.  It records the exact SweepGA
and IMPG contracts and refuses to manufacture a success sentinel until the
full production implementation is independently reviewed.  The repaired gate
currently is NO_GO, so the login-node runner cannot submit this dormant worker.
Keeping this boundary explicit prevents a later GO from silently falling back
to an older WFMASH-first or demographic workflow.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from analysis.tier3_common import Tier3ValidationError, sha256_file


SWEEPGA_CONTRACT = ["sweepga", "H1.fa", "H2.fa", "--output-file", "whole.1to1.paf", "--num-mappings", "1:1", "--scaffold-jump", "0"]
IMPG_CONTRACTS = {
    "index": ["impg", "index", "-a", "whole.1to1.paf", "-i", "whole.impg"],
    "partition": ["impg", "partition", "-a", "whole.1to1.paf", "-i", "whole.impg", "-o", "bed"],
    "query": ["impg", "query", "-a", "whole.1to1.paf", "-i", "whole.impg", "-b", "native-focus.bed", "-o", "vcf:poa"],
}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-ledger", type=Path, required=True)
    parser.add_argument("--staged-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--gate-sha256", required=True)
    parser.add_argument("--authorization-tuple", required=True)
    args = parser.parse_args(argv)
    if not args.input_ledger.is_file() or not args.staged_root.is_dir():
        raise Tier3ValidationError("authorized staged inputs are unavailable")
    args.output.mkdir(parents=True, exist_ok=True)
    contract = {
        "run_id": args.run_id,
        "gate_sha256": args.gate_sha256,
        "authorization_tuple_digest": args.authorization_tuple,
        "input_ledger_sha256": sha256_file(args.input_ledger),
        "sweepga": SWEEPGA_CONTRACT,
        "impg": IMPG_CONTRACTS,
        "annotation_targets": "original exact H1-native annotation only",
        "denominators": ["callable_bases", "callable_fraction", "queryable_gene_count", "queryable_gene_bases", "target_gene_total", "target_base_total"],
        "forbidden": ["compute-node download", "population interpretation of H1/H2", "PSMC", "MSMC2", "SMC++", "demographic inference"],
    }
    (args.output / "refused_unreviewed_compute.json").write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    raise Tier3ValidationError("authorized scientific compute implementation has not been independently reviewed; refusing without success sentinel")


if __name__ == "__main__":
    raise SystemExit(main())
