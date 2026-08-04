# Mouse Dmd exon 44 / C3H / C2C12 pipeline

This command-line script selects the annotated full-length muscle **Dmd/Dp427m** transcript in GRCm39, numbers exons in transcriptional (5′→3′) order, extracts exon 44 and introns 43/44, maps the exon-43-to-exon-45 locus to C3H, and compares transcript-oriented sequences. Mouse Dmd is normally on the minus strand; all reported genomic coordinates are **1-based inclusive**, while FASTA sequences are transcript-oriented and clearly labeled.

## Install

Python 3.10+ and Biopython are required. `minimap2` is preferred when available. For closely related mouse assemblies, the script can fall back to exact 41-nt anchors sampled across the locus; it refuses the mapping unless at least 10 anchors support a consistent chromosome-X location.

A ready-to-use Windows environment is included in `.venv`. Activate it in PowerShell with:

```powershell
& ".\.venv\Scripts\Activate.ps1"
python .\dmd_exon44_pipeline.py --help
```

```bash
python -m pip install biopython rapidfuzz
```

## Run discovery/comparison

```bash
python dmd_exon44_pipeline.py \
  --grcm39-fasta GRCm39.genome.fa \
  --grcm39-gtf GRCm39.annotation.gtf \
  --c3h-fasta C3H.genome.fa \
  --c3h-gbff C3H.sequence_and_annotation.gbff \
  --outdir dmd_results
```

Review `report.json` before relying on the result. It records the selected transcript and why it was chosen, every Dmd transcript seen, the C3H mapping target/strand/MAPQ, GBFF Dmd features, coordinates, sequences, and differences. If the annotation does not explicitly call a transcript Dp427m, the script selects the most full-length protein-coding Dmd transcript and says that it used the fallback. To remove ambiguity, rerun with `--transcript-id ENSMUST...`.

Outputs include:

- `report.json`: complete machine-readable audit report
- `regions.fasta`: exon 44 and introns 43/44 for both assemblies, transcript-oriented
- `variants.tsv`: substitutions and indels by region
- `candidate_spcas9_guides.tsv`: preliminary 20-nt SpCas9 candidates with NGG PAMs in the exon-44 flank (default ±500 bp)

Candidate guides are sequence candidates only. They are not checked for whole-genome off-targets, efficacy, repeats, polymorphic PAMs, or experimental suitability.

## Add C2C12 Sanger reads

Pass cleaned single-record FASTA files or ABI/AB1 chromatograms:

```bash
python dmd_exon44_pipeline.py [the four genome options above] \
  --sanger sample1.ab1 sample2.fasta \
  --outdir dmd_results_with_sanger
```

For AB1, the report includes mean PHRED quality, positions below Q20, ambiguous base calls, and positions where the second-highest raw trace channel is at least 33% of the primary peak. Every substitution/indel includes local read and target coordinates, PHRED quality, and a nearby overlapping-peak flag. Classification is based on local alignment against both assembly loci. ABI channel scaling and heterozygous/indel mixtures still require trace-aware manual review and replicate/strand confirmation.

Candidate-guide overlap can be determined directly from the reported target coordinate and `candidate_spcas9_guides.tsv`: the guide interval, PAM (last 3 nt), and PAM-proximal seed (last 10 spacer nt) are all explicit. Primer overlap requires the actual primer sequences or binding coordinates; those cannot be inferred safely before primers are supplied.

## Important validation notes

- C2C12 is historically derived from C3H muscle, but cell-line drift and subclonal variation make direct PCR/Sanger confirmation essential.
- The script uses transcript order, not increasing genomic coordinate, to define exon and intron numbers on the minus strand.
- Exact exon numbering can differ among annotations. The report flags disagreement between supplied `exon_number` values and computed transcript order.
- A high-confidence C3H mapping should be unique, span the locus, and have a strong MAPQ. Inspect it before guide or primer design.
