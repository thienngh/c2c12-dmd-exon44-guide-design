#!/usr/bin/env python3
"""Locate/compare mouse Dmd Dp427m exon 44 and adjacent introns.

Requires Biopython. C3H fallback mapping additionally requires minimap2 on PATH.
Coordinates in TSV/JSON are 1-based inclusive; BED is 0-based half-open.
"""
from __future__ import annotations

import argparse, csv, gzip, json, re, shutil, subprocess, sys, tempfile
import statistics
from collections import defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional

from Bio import Align, SeqIO
from Bio.Seq import Seq
from rapidfuzz.distance import Levenshtein


def opn(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def attrs_gtf(text):
    out = {}
    for item in text.strip().strip(";").split(";"):
        item = item.strip()
        if not item: continue
        if " " in item:
            k, v = item.split(None, 1); out[k] = v.strip().strip('"')
        elif "=" in item:
            k, v = item.split("=", 1); out[k] = v
    return out


@dataclass
class Exon:
    chrom: str; start: int; end: int; strand: str; transcript_id: str
    annotated_number: Optional[int] = None; order_number: Optional[int] = None


@dataclass
class Region:
    name: str; chrom: str; start: int; end: int; strand: str; sequence: str

    def oriented(self):
        return self.sequence if self.strand == "+" else str(Seq(self.sequence).reverse_complement())


def parse_gtf(path):
    tx = defaultdict(list); meta = defaultdict(dict)
    with opn(path) as fh:
        for line in fh:
            if not line or line.startswith("#"): continue
            p = line.rstrip().split("\t")
            if len(p) != 9: continue
            chrom, _, typ, start, end, _, strand, _, raw = p
            a = attrs_gtf(raw); tid = a.get("transcript_id")
            is_dmd = a.get("gene_name", "").lower() == "dmd" or a.get("gene_id", "").lower() == "dmd"
            if not tid or not is_dmd: continue
            meta[tid].update({k:a.get(k) for k in ("gene_id","gene_name","transcript_name","transcript_biotype","transcript_type") if a.get(k)})
            if typ == "exon":
                n = a.get("exon_number")
                tx[tid].append(Exon(chrom,int(start),int(end),strand,tid,int(n) if n and n.isdigit() else None))
    for tid, exons in tx.items():
        exons.sort(key=lambda e:e.start, reverse=exons[0].strand == "-")
        for i,e in enumerate(exons,1): e.order_number=i
    return tx, meta


def select_dp427m(tx, meta, requested=None):
    if requested:
        if requested not in tx: raise SystemExit(f"Transcript {requested!r} was not found among Dmd transcripts")
        return requested, "explicit --transcript-id"
    def score(tid):
        m = " ".join(str(v) for v in meta[tid].values()).lower()
        named = 1 if re.search(r"dp\s*427\s*m|dp427m|dystrophin.*muscle", m) else 0
        coding = 1 if "protein_coding" in m else 0
        span = sum(e.end-e.start+1 for e in tx[tid])
        return named, coding, len(tx[tid]), span
    tid = max(tx, key=score)
    reason = "annotation name identifies Dp427m" if score(tid)[0] else "fallback: most full-length protein-coding Dmd transcript by exon count/spliced length"
    return tid, reason


def fasta_index(path, cache):
    return SeqIO.index_db(str(cache), [str(path)], "fasta")


def fasta_key(index, chrom):
    if chrom in index: return chrom
    aliases = [chrom.removeprefix("chr"), "chr"+chrom]
    for k in index:
        if k.split()[0] in aliases or k.split(".")[0] == chrom.split(".")[0]: return k
    raise KeyError(f"Sequence {chrom!r} not found in FASTA")


def fetch(index, chrom, start, end):
    key=fasta_key(index,chrom); return str(index[key].seq[start-1:end]).upper(), key


def regions_from_exons(exons, index):
    bynum={e.order_number:e for e in exons}
    if not all(n in bynum for n in (43,44,45)): raise SystemExit("Selected transcript has fewer than 45 exons")
    e43,e44,e45=(bynum[n] for n in (43,44,45)); strand=e44.strand
    if strand == "+": i43=(e43.end+1,e44.start-1); i44=(e44.end+1,e45.start-1)
    else: i43=(e44.end+1,e43.start-1); i44=(e45.end+1,e44.start-1)
    specs=[("exon44",e44.start,e44.end),("intron43",*i43),("intron44",*i44)]
    out={}
    for name,s,e in specs:
        seq,key=fetch(index,e44.chrom,s,e); out[name]=Region(name,key,s,e,strand,seq)
    locus_s=min(e43.start,e44.start,e45.start); locus_e=max(e43.end,e44.end,e45.end)
    seq,key=fetch(index,e44.chrom,locus_s,locus_e); out["locus43_45"]=Region("locus43_45",key,locus_s,locus_e,strand,seq)
    return out


def scan_gbff(path):
    hits=[]
    with opn(path,"rt") as fh:
        for rec in SeqIO.parse(fh,"genbank"):
            for f in rec.features:
                q=f.qualifiers; gene=" ".join(q.get("gene",[])+q.get("gene_synonym",[]))
                if gene.lower() == "dmd" or re.search(r"\bdmd\b",gene,re.I):
                    hits.append({"record":rec.id,"type":f.type,"start":int(f.location.start)+1,"end":int(f.location.end),"strand":"+" if f.location.strand==1 else "-", "qualifiers":q})
    return hits


def run_minimap(ref_fa, query_seq, work):
    if not shutil.which("minimap2"):
        return run_anchor_mapping(ref_fa, query_seq)
    q=work/"grcm39_locus.fa"; q.write_text(">GRCm39_Dmd_ex43_45\n"+query_seq+"\n")
    cp=subprocess.run(["minimap2","-x","asm5","--secondary=no",str(ref_fa),str(q)],capture_output=True,text=True,check=True)
    rows=[x.split("\t") for x in cp.stdout.splitlines() if x.strip()]
    if not rows: raise SystemExit("No C3H mapping found for GRCm39 Dmd exon43-45 locus")
    r=max(rows,key=lambda x:int(x[11])); return {"qstart":int(r[2]),"qend":int(r[3]),"strand":r[4],"target":r[5],"tstart":int(r[7]),"tend":int(r[8]),"mapq":int(r[11])}


def run_anchor_mapping(ref_fa, query_seq, k=41, step=500):
    """Map a closely related mouse locus using consistent exact anchors."""
    anchors={}
    # Twenty-five distributed anchors keep the fallback fast on chromosome-scale DNA.
    step=max(step,max(1,(len(query_seq)-k)//24))
    for pos in range(0,max(1,len(query_seq)-k+1),step):
        s=query_seq[pos:pos+k]
        if len(s)==k and set(s)<=set("ACGT") and max(s.count(b) for b in "ACGT")<int(k*.8):
            anchors.setdefault(s,[]).append(pos)
    if len(anchors)<10: raise SystemExit("Too few usable reference anchors for C3H fallback mapping")
    rcanchors={}
    for s,poss in anchors.items(): rcanchors.setdefault(str(Seq(s).reverse_complement()),[]).extend(poss)
    candidates=[]
    with opn(ref_fa) as fh:
        for rec in SeqIO.parse(fh,"fasta"):
            desc=rec.description.lower()
            if "chromosome: x" not in desc and "chromosome x" not in desc and rec.id!="OW971861.1": continue
            target_seq=str(rec.seq).upper()
            for strand,lookup in (("+",anchors),("-",rcanchors)):
                starts=[]
                for seed,ref_positions in lookup.items():
                    hit=target_seq.find(seed)
                    # Ignore non-unique anchors; unique exact seeds give auditable placement.
                    if hit<0 or target_seq.find(seed,hit+1)>=0: continue
                    for rp in ref_positions:
                        starts.append(hit-rp if strand=="+" else hit-len(query_seq)+rp+k)
                if starts:
                    bins=defaultdict(list)
                    for x in starts: bins[round(x/1000)].append(x)
                    cluster=max(bins.values(),key=len)
                    candidates.append((len(cluster),strand,int(statistics.median(cluster)),rec.id,len(starts)))
            break
    if not candidates: raise SystemExit("No exact-anchor mapping found on C3H chromosome X")
    support,strand,start,target_id,total=max(candidates)
    if support<10: raise SystemExit(f"C3H anchor mapping was not reliable: only {support} consistent anchors")
    start=max(0,start); end=start+len(query_seq)
    return {"qstart":0,"qend":len(query_seq),"strand":strand,"target":target_id,"tstart":start,"tend":end,
            "mapq":None,"method":"exact_anchor_fallback","consistent_anchor_support":support,"all_anchor_hits":total,"anchor_k":k,"anchor_step":step}


def global_alignment(a,b):
    al=Align.PairwiseAligner(); al.mode="global"; al.match_score=2; al.mismatch_score=-3; al.open_gap_score=-7; al.extend_gap_score=-1
    return al.align(a,b)[0]


def variants(ref, alt):
    out=[]
    for tag,rs,re,as_,ae in Levenshtein.opcodes(ref,alt):
        if tag=="equal": continue
        if tag=="delete": out.append({"type":"deletion","ref_pos":rs+1,"ref":ref[rs:re],"alt":"-"})
        elif tag=="insert": out.append({"type":"insertion","ref_pos":rs+1,"ref":"-","alt":alt[as_:ae]})
        elif (re-rs)==(ae-as_):
            for i,(x,y) in enumerate(zip(ref[rs:re],alt[as_:ae])):
                if x!=y: out.append({"type":"substitution","ref_pos":rs+i+1,"ref":x,"alt":y})
        else: out.append({"type":"replacement","ref_pos":rs+1,"ref":ref[rs:re],"alt":alt[as_:ae]})
    return out


def project_regions(ref_locus, alt_oriented, mapping, regions):
    """Project reference offsets through an optimized Levenshtein alignment."""
    ops=list(Levenshtein.opcodes(ref_locus.sequence,alt_oriented))
    def project(p):
        for tag,rs,re,as_,ae in ops:
            if rs<=p<re:
                if tag=="delete": return as_
                return as_+min(p-rs,max(0,ae-as_-1))
        return len(alt_oriented)-1
    out={}
    for name,r in regions.items():
        if name=="locus43_45": continue
        lo=r.start-ref_locus.start; hi=r.end-ref_locus.start
        a,b=project(lo),project(hi); a,b=min(a,b),max(a,b); seq=alt_oriented[a:b+1]
        if mapping["strand"]=="+": gs=mapping["tstart"]+a+1; ge=mapping["tstart"]+b+1
        else: gs=mapping["tend"]-b; ge=mapping["tend"]-a
        genomic=seq if mapping["strand"]=="+" else str(Seq(seq).reverse_complement())
        out[name]=Region(name,mapping["target"],gs,ge, "+" if (r.strand==mapping["strand"]) else "-",genomic)
    return out


def guides(flank_seq, genomic_start, strand, pam="NGG"):
    ans=[]
    for orientation,s in [("+",flank_seq),("-",str(Seq(flank_seq).reverse_complement()))]:
        for m in re.finditer(r"(?=([ACGT]{20}[ACGT]GG))",s):
            spacer=m.group(1)[:20]; pos=m.start()
            if orientation=="+": gs=genomic_start+pos; ge=gs+22
            else: ge=genomic_start+len(s)-pos-1; gs=ge-22
            ans.append({"spacer":spacer,"pam":m.group(1)[20:],"guide_strand":orientation,"start":gs,"end":ge,"seed":spacer[-10:]})
    return ans


def read_sanger(path):
    fmt="abi" if str(path).lower().endswith(".ab1") else "fasta"
    rec=SeqIO.read(path,fmt); q=rec.letter_annotations.get("phred_quality",[])
    return str(rec.seq).upper(), q, rec


def trace_peak_metrics(rec):
    """Return second/highest channel ratio at each ABI base call when raw tags exist."""
    raw=rec.annotations.get("abif_raw",{}); loc=raw.get("PLOC2") or raw.get("PLOC1")
    order=raw.get("FWO_1",b"GATC"); order=order.decode(errors="ignore") if isinstance(order,bytes) else str(order)
    channels=[raw.get(f"DATA{i}") for i in range(9,13)]
    if not loc or any(x is None for x in channels): return []
    ans=[]
    for p in loc:
        vals=sorted((int(ch[p]),order[i] if i<len(order) else "?") for i,ch in enumerate(channels))
        hi=max(vals[-1][0],1); ans.append({"primary":vals[-1][1],"secondary":vals[-2][1],"secondary_primary_ratio":round(vals[-2][0]/hi,3)})
    return ans


def local_result(target, read):
    al=Align.PairwiseAligner(); al.mode="local"; al.match_score=2; al.mismatch_score=-3; al.open_gap_score=-6; al.extend_gap_score=-1
    a=al.align(target,read)[0]; vs=[]; rb=int(a.aligned[0][0][0]); ab=int(a.aligned[1][0][0]); matches=0; compared=0
    for (rs,re),(qs,qe) in zip(a.aligned[0],a.aligned[1]):
        rs,re,qs,qe=map(int,(rs,re,qs,qe))
        if rs>rb: vs.append({"type":"deletion","target_pos":rb+1,"read_pos":ab+1,"ref":target[rb:rs],"alt":"-"})
        if qs>ab: vs.append({"type":"insertion","target_pos":rb+1,"read_pos":ab+1,"ref":"-","alt":read[ab:qs]})
        for i,(x,y) in enumerate(zip(target[rs:re],read[qs:qe])):
            compared+=1
            if x==y: matches+=1
            else: vs.append({"type":"substitution","target_pos":rs+i+1,"read_pos":qs+i+1,"ref":x,"alt":y})
        rb,ab=re,qe
    return {"score":a.score,"target_aligned_start":int(a.aligned[0][0][0])+1,"target_aligned_end":int(a.aligned[0][-1][1]),
            "read_aligned_start":int(a.aligned[1][0][0])+1,"read_aligned_end":int(a.aligned[1][-1][1]),
            "identity":matches/compared if compared else 0,"variants":vs}


def classify_sanger(seq, refs):
    scored=[]
    for label,target in refs.items():
        for orient,s in [("forward",seq),("reverse_complement",str(Seq(seq).reverse_complement()))]:
            x=local_result(target,s); scored.append((x["score"],label,orient,x))
    scored.sort(key=lambda x:x[0],reverse=True); best=scored[0]
    tied=sorted({x[1] for x in scored if x[0]==best[0]})
    return {"best_score":best[0],"matches":tied,"orientation":best[2],"classification":"both" if len(tied)>1 else tied[0],"alignment":best[3]}


def main():
    p=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--grcm39-fasta",required=True); p.add_argument("--grcm39-gtf",required=True)
    p.add_argument("--c3h-fasta",required=True); p.add_argument("--c3h-gbff")
    p.add_argument("--transcript-id"); p.add_argument("--flank",type=int,default=500)
    p.add_argument("--sanger",nargs="*",default=[]); p.add_argument("--outdir",default="dmd_exon44_results")
    a=p.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    tx,meta=parse_gtf(a.grcm39_gtf)
    if not tx: raise SystemExit("No Dmd transcripts found in GTF (expected gene_name Dmd)")
    tid,reason=select_dp427m(tx,meta,a.transcript_id)
    refidx=fasta_index(a.grcm39_fasta,out/"grcm39.seqio.idx"); c3idx=fasta_index(a.c3h_fasta,out/"c3h.seqio.idx")
    rr=regions_from_exons(tx[tid],refidx); locus=rr["locus43_45"]
    gbhits=scan_gbff(a.c3h_gbff) if a.c3h_gbff else []
    with tempfile.TemporaryDirectory() as td:
        mp=run_minimap(a.c3h_fasta,locus.sequence,Path(td))
    altseq,_=fetch(c3idx,mp["target"],mp["tstart"]+1,mp["tend"])
    if mp["strand"]=="-": altseq=str(Seq(altseq).reverse_complement())
    cr=project_regions(locus,altseq,mp,rr)
    report={"transcript_id":tid,"selection_reason":reason,"transcript_metadata":meta[tid],"strand":tx[tid][0].strand,
            "coordinate_system":"1-based inclusive","grcm39":{},"c3h":{},"c3h_mapping":mp,"c3h_gbff_dmd_features":gbhits,"comparisons":{},"other_dmd_transcripts":[]}
    for name,r in rr.items():
        if name!="locus43_45": report["grcm39"][name]={**asdict(r),"transcript_oriented_sequence":r.oriented()}
    for name,r in cr.items(): report["c3h"][name]={**asdict(r),"transcript_oriented_sequence":r.oriented()}
    for name in ("exon44","intron43","intron44"):
        if name in cr: report["comparisons"][name]=variants(rr[name].oriented(),cr[name].oriented())
    for other,ex in sorted(tx.items()):
        report["other_dmd_transcripts"].append({"transcript_id":other,"exon_count":len(ex),"metadata":meta[other],"selected":other==tid,
            "annotated_vs_order_number_disagree":any(e.annotated_number and e.annotated_number!=e.order_number for e in ex)})
    e44=rr["exon44"]; fs=max(1,e44.start-a.flank); fe=e44.end+a.flank; fseq,_=fetch(refidx,e44.chrom,fs,fe)
    report["grcm39_exon44_flank"]={"chrom":e44.chrom,"start":fs,"end":fe,"sequence":fseq}
    report["candidate_spcas9_guides"]=guides(fseq,fs,e44.strand)
    if a.sanger:
        refs={"GRCm39":locus.sequence,"C3H":altseq}
        report["sanger"]={}
        for f in a.sanger:
            seq,q,rec=read_sanger(f); x=classify_sanger(seq,refs)
            peaks=trace_peak_metrics(rec); oriented_seq=seq if x["orientation"]=="forward" else str(Seq(seq).reverse_complement())
            for v in x["alignment"]["variants"]:
                rp=v["read_pos"]-1
                # In reverse orientation, translate the position back to original base-call order.
                original_rp=rp if x["orientation"]=="forward" else len(seq)-rp-1
                v["phred"]=q[original_rp] if q and 0<=original_rp<len(q) else None
                v["ambiguous_peak"]=bool(peaks and 0<=original_rp<len(peaks) and peaks[original_rp]["secondary_primary_ratio"]>=0.33)
                v["feature_overlaps"]=[]
                if x["classification"] in ("GRCm39","both"):
                    gp=locus.start+v["target_pos"]-1; v["grcm39_genomic_position"]=gp
                    for g in report["candidate_spcas9_guides"]:
                        if not (g["start"]<=gp<=g["end"]): continue
                        if g["guide_strand"]=="+":
                            pam=(g["end"]-2,g["end"]); seed=(g["start"]+10,g["start"]+19)
                        else:
                            pam=(g["start"],g["start"]+2); seed=(g["start"]+3,g["start"]+12)
                        part="PAM" if pam[0]<=gp<=pam[1] else ("seed" if seed[0]<=gp<=seed[1] else "spacer")
                        v["feature_overlaps"].append({"spacer":g["spacer"],"part":part,"guide_start":g["start"],"guide_end":g["end"]})
            x.update({"length":len(seq),"mean_phred":sum(q)/len(q) if q else None,"low_quality_positions":[i+1 for i,v in enumerate(q) if v<20],
                      "ambiguous_positions":[i+1 for i,b in enumerate(seq) if b not in "ACGT"],
                      "overlapping_peak_positions":[i+1 for i,v in enumerate(peaks) if v["secondary_primary_ratio"]>=0.33],
                      "overlapping_peak_threshold":"secondary/primary >= 0.33", "oriented_sequence":oriented_seq})
            report["sanger"][str(f)]=x
    (out/"report.json").write_text(json.dumps(report,indent=2))
    with open(out/"regions.fasta","w") as fh:
        for strain,d in (("GRCm39",rr),("C3H",cr)):
            for name,r in d.items():
                if name=="locus43_45": continue
                fh.write(f">{strain}|{name}|{r.chrom}:{r.start}-{r.end}({r.strand})|transcript_oriented\n{r.oriented()}\n")
    with open(out/"variants.tsv","w",newline="") as fh:
        w=csv.writer(fh,delimiter="\t"); w.writerow(["region","type","region_position","reference","C3H"])
        for reg,vs in report["comparisons"].items():
            for v in vs:w.writerow([reg,v["type"],v["ref_pos"],v["ref"],v["alt"]])
    with open(out/"candidate_spcas9_guides.tsv","w",newline="") as fh:
        keys=["spacer","pam","guide_strand","start","end","seed"]; w=csv.DictWriter(fh,keys,delimiter="\t"); w.writeheader(); w.writerows(report["candidate_spcas9_guides"])
    print(f"Selected {tid}: {reason}"); print(f"Results: {out.resolve()}")

if __name__=="__main__": main()
