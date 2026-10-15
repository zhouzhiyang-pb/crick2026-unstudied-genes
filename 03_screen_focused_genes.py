#!/usr/bin/env python3
"""
Step 3: screen OMA and PANTHER members against UniProt, and output one
candidate table per gene.

Keep rule:
  A member is kept if it is reviewed (Swiss-Prot) OR it carries experimental
  evidence (ECO:0000269 or ECO:0000314).

Tier rule (from strongest to weakest):
  star = has experimental evidence (ECO:0000269 or ECO:0000314)
  ring = has author-statement evidence (ECO:0000303), no experimental code
  triangle = has curator-inference or similarity evidence
             (ECO:0000305 or ECO:0000250), no stronger code
  dash  = reviewed but none of the above codes

For every kept member we also collect the PubMed identifiers that back its
annotations, so the user knows which papers to read.

Output (two files per gene, inside the gene's folder):
  <GENE>/candidate_table.csv    -- which proteins are worth reading about
  <GENE>/annotations_table.csv  -- what each protein does (function, catalytic
                                   activity, cofactor, subunit) and any annotated
                                   residues, so the user can skim instead of
                                   reading every paper
Each protein appears once. If a protein is listed by both OMA and PANTHER,
its "source" column is "OMA & PANTHER".

OMA is queried only through its largest HOG (the top level), because every
sub-HOG is a subset of its parent HOG, so the largest one already contains
all members.

Evidence classification reference:
  https://www.uniprot.org/help/evidences

Run with:  python3 03_screen_focused_genes.py
"""

import csv
import glob
import json
import re
import subprocess
import time

# The two genes under study. Each has a folder of the same name.
GENES = ["SPBC115.03", "SPAPB1A10.12c"]

# UniProt REST base URL.
UNIPROT = "https://rest.uniprot.org/uniprotkb"

# A UniProt accession looks like 6-10 alphanumerics (no underscore).
ACCESSION_RE = re.compile(r"^[A-Z0-9]{6,10}$")

# A UniProt entry name looks like GENE_ORG where ORG is a 5-letter
# organism code (e.g. YEG9_SCHPO, ALO1_YEAST, A0A5P8YFB9_YERPE).
# Other ids with underscores (RefSeq XP_..., model-organism ids like
# MGT_11527) are not in UniProt and are skipped.
ENTRY_NAME_RE = re.compile(r"^[A-Za-z0-9]+_[A-Z]{5}$")

# Evidence codes grouped by tier, strongest first.
EXPERIMENTAL_CODES = ["ECO:0000269", "ECO:0000314"]
AUTHOR_STATEMENT_CODE = "ECO:0000303"
INFERENCE_CODES = ["ECO:0000305", "ECO:0000250"]


def fetch(url):
    """GET a URL with curl, retrying a few times. Returns the body or None."""
    for attempt in range(4):
        result = subprocess.run(
            ["curl", "-sL", "-w", "\n%{http_code}", "--max-time", "60", url],
            capture_output=True,
            text=True,
        )
        output = result.stdout
        if output.endswith("200"):
            return output[:-4]  # strip the HTTP code we appended
        time.sleep(2)
    return None


def collect_oma_members(fasta_path):
    """
    Parse an OMA fasta.txt. Each header looks like:
      >OMA_ID | UniProtAC_or_entryname | sub-HOG | [species]
    Return a list of (id, id2, subhog) tuples.
    """
    members = []
    for line in open(fasta_path):
        if not line.startswith(">"):
            continue
        parts = [p.strip() for p in line[1:].strip().split("|")]
        if len(parts) >= 3:
            members.append((parts[0], parts[1], parts[2]))
    return members


def collect_panther_members(fasta_path):
    """
    Parse a PANTHER cluster.fasta. Headers look like:
      >PREFIX|...|UniProtKB=ACC
    Return a list of (label, accession) tuples.
    """
    members = []
    for line in open(fasta_path):
        if not line.startswith(">"):
            continue
        match = re.search(r"UniProtKB=([A-Za-z0-9_]+)", line)
        if match:
            members.append((line[1:].strip().split("|")[0], match.group(1)))
    return members


def lookup_accessions_batch(accessions):
    """Batch-lookup accessions via the accessions endpoint (chunks of 25)."""
    records = {}
    for i in range(0, len(accessions), 25):
        chunk = accessions[i : i + 25]
        # Fetch the full record (no "fields" filter) so the ECO evidence codes are
        # included; the default field subset drops the comments/features that carry
        # the experimental evidence (ECO:0000269 / ECO:0000314).
        url = (
            f"{UNIPROT}/accessions?accessions="
            + ",".join(chunk)
            + "&format=json"
        )
        body = fetch(url)
        if not body:
            continue
        try:
            for record in json.loads(body).get("results", []):
                records[record["primaryAccession"]] = record
        except json.JSONDecodeError:
            pass
        time.sleep(0.4)
    return records


def lookup_entry_name(entry_name):
    """Look up a single entry name. Returns the record or None."""
    body = fetch(f"{UNIPROT}/{entry_name}.json")
    if not body:
        return None
    try:
        record = json.loads(body)
    except json.JSONDecodeError:
        return None
    if "primaryAccession" in record:
        return record
    return None


def extract_evidence(record):
    """
    Walk a UniProt record and collect every (evidence_code, source, source_id) triple.
    UniProt stores them as dicts like
      {"evidenceCode": "ECO:0000269", "source": "PubMed", "id": "10477285"}.
    Return a list of (code, source, id) tuples.
    """
    triples = []

    def walk(node):
        if isinstance(node, dict):
            if "evidenceCode" in node:
                triples.append((
                    node["evidenceCode"],
                    node.get("source", ""),
                    node.get("id", ""),
                ))
            for value in node.values():
                walk(value)
        elif isinstance(node, list):
            for item in node:
                walk(item)

    walk(record)
    return triples


def tier_of(codes):
    """Map a set of evidence codes to a tier symbol (star/ring/triangle/dash)."""
    if any(code in codes for code in EXPERIMENTAL_CODES):
        return "star"
    if AUTHOR_STATEMENT_CODE in codes:
        return "ring"
    if any(code in codes for code in INFERENCE_CODES):
        return "triangle"
    return "dash"


def count_sequences(fasta_path):
    """Count the number of FASTA headers (sequences) in a file."""
    return sum(1 for line in open(fasta_path) if line.startswith(">"))


def comment_text(record, comment_type):
    """Return the joined text of the first comment of the given type, or ""."""
    for c in record.get("comments", []):
        if c.get("commentType") == comment_type:
            return " ".join(t.get("value", "") for t in c.get("texts", []))
    return ""


def catalytic_activity(record):
    """Return all CATALYTIC ACTIVITY reaction names, joined by '; '."""
    reactions = []
    for c in record.get("comments", []):
        if c.get("commentType") == "CATALYTIC ACTIVITY":
            name = c.get("reaction", {}).get("name", "")
            if name:
                reactions.append(name)
    return "; ".join(reactions)


def cofactor_names(record):
    """Return all COFACTOR molecule names, joined by '; '."""
    names = []
    for c in record.get("comments", []):
        if c.get("commentType") == "COFACTOR":
            for cf in c.get("cofactors", []):
                if cf.get("name"):
                    names.append(cf["name"])
    return "; ".join(names)


def important_residues(record):
    """
    Return a string of position-annotated residues from the feature table:
    active sites, binding sites, metal binding, and modified residues.
    These are not annotated for every protein, so the string is often empty.
    """
    keep = {"Active site", "Binding site", "Metal binding", "Modified residue", "Site"}
    parts = []
    for f in record.get("features", []):
        if f.get("type") in keep:
            pos = f.get("location", {}).get("start", {}).get("value", "")
            desc = f.get("description", "") or ""
            parts.append("{} {} {}".format(f["type"], pos, desc).strip())
    return "; ".join(parts)


def main():
    # --- 1. collect members, one entry per unique UniProt identifier ---
    # members[gene][identifier] = set of sources ("OMA", "PANTHER") that listed it.
    members = {gene: {} for gene in GENES}

    for gene in GENES:
        # OMA: use only the largest HOG (the top level). A sub-HOG is a subset of
        # its parent HOG, so the largest one already contains every member.
        oma_fastas = glob.glob(f"{gene}/OMA/*/fasta.txt")
        if oma_fastas:
            top_fasta = max(oma_fastas, key=count_sequences)
            for omaid, id2, subhog in collect_oma_members(top_fasta):
                if ACCESSION_RE.match(id2) or ENTRY_NAME_RE.match(id2):
                    members[gene].setdefault(id2, set()).add("OMA")
                # else: a non-UniProt id (RefSeq / model-organism) -- skip.

        # PANTHER: one family per gene.
        for fasta in glob.glob(f"{gene}/PANTHER/*/cluster.fasta"):
            for label, ac in collect_panther_members(fasta):
                if ACCESSION_RE.match(ac) or ENTRY_NAME_RE.match(ac):
                    members[gene].setdefault(ac, set()).add("PANTHER")
                # else: not a UniProt id -- skip.

    # --- 2. look up every identifier ---
    # Accessions are fetched in one batch. Entry names are fetched one by one
    # and then re-keyed under their resolved accession.
    accessions = sorted({
        ident for gene in GENES for ident in members[gene] if ACCESSION_RE.match(ident)
    })
    records = lookup_accessions_batch(accessions)

    for gene in GENES:
        for ident in list(members[gene]):
            if not ENTRY_NAME_RE.match(ident):
                continue
            record = lookup_entry_name(ident)
            if record:
                accession = record["primaryAccession"]
                records[accession] = record
                # Move the entry name's sources under its resolved accession.
                members[gene].setdefault(accession, set()).update(members[gene].pop(ident))

    # --- 3. apply the keep + tier rules and write one table per gene ---
    header = ["source", "uniprot_ac", "reviewed",
              "tier", "evidence_codes", "pubmed_ids", "protein_name"]
    annotation_header = ["uniprot_ac", "protein_name", "tier",
                         "function", "catalytic_activity", "cofactor", "subunit",
                         "residues", "pubmed_ids"]

    for gene in GENES:
        rows = []
        annotation_rows = []
        for identifier, sources in sorted(members[gene].items()):
            record = records.get(identifier)
            if record is None:
                continue

            reviewed = record.get("entryType", "").startswith("UniProtKB reviewed")
            evidence = extract_evidence(record)
            codes = {code for code, _, _ in evidence}
            has_experimental = any(code in codes for code in EXPERIMENTAL_CODES)

            # Keep only reviewed members or members with experimental evidence.
            if not reviewed and not has_experimental:
                continue

            tier = tier_of(codes)
            # PubMed IDs live in the "id" field of evidence whose source is "PubMed".
            pubmed_ids = sorted({
                source_id for _, source, source_id in evidence if source == "PubMed"
            })
            name = (
                record.get("proteinDescription", {})
                .get("recommendedName", {})
                .get("fullName", {})
                .get("value", "")
            )

            rows.append([
                " & ".join(sorted(sources)),
                record["primaryAccession"],
                "reviewed" if reviewed else "unreviewed",
                tier,
                ",".join(sorted(codes)),
                ",".join(pubmed_ids),
                name,
            ])

            # Extra curated annotations, so the user can see what each protein
            # does and which residues matter without opening every paper.
            annotation_rows.append([
                record["primaryAccession"],
                name,
                tier,
                comment_text(record, "FUNCTION"),
                catalytic_activity(record),
                cofactor_names(record),
                comment_text(record, "SUBUNIT"),
                important_residues(record),
                ",".join(pubmed_ids),
            ])

        out = f"{gene}/candidate_table.csv"
        with open(out, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(header)
            writer.writerows(rows)

        ann_out = f"{gene}/annotations_table.csv"
        with open(ann_out, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(annotation_header)
            writer.writerows(annotation_rows)

        print(f"{gene}: kept {len(rows)} -> {out} and {ann_out}")


if __name__ == "__main__":
    main()
