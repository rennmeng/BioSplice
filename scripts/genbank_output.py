#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GenBank record generation.

Two GenBank records are written per event:
  - EXON-*.gb : mRNA sequence spliced from the reference transcript exons
                (exon/CDS features plus event-region annotations)
  - Tran-*.gb : genomic sequence spanned by the reference transcript
                (exon/CDS features plus event-region annotations)
"""

import logging
import os
import re
from datetime import datetime

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio import SeqIO

from .events import _extract_event_regions
from .matching import _exon_in_ref
from .sequences import get_exon_sequence, safe_get_seq

logger = logging.getLogger(__name__)


def sanitize_filename(s: str) -> str:
    return re.sub(r'[^\w\-\.]', '_', str(s).strip())


def generate_exon_gb(ref_tid, exon_bounds, cds_ph_list, genome, chrom, strand,
                     gene_symbol, etype, row, output_dir, job_id, event_id=""):
    if not exon_bounds or not ref_tid:
        return None
    try:
        is_minus = (strand == '-')
        n_exons = len(exon_bounds)

        mrna_seq = get_exon_sequence(exon_bounds, genome, chrom, strand)
        if not mrna_seq:
            return None
        seq_len = len(mrna_seq)
        feat_strand = 1

        def _mirror(pos):
            return seq_len - pos if is_minus else pos

        safe_gene = sanitize_filename(gene_symbol)
        id_prefix = f"{job_id}_" if job_id else ""
        locus_name = f"{id_prefix}{event_id}_EXON-{etype}-{safe_gene}"

        record = SeqRecord(Seq(mrna_seq), id=ref_tid, name=locus_name[:16],
                           description=f"{gene_symbol} {etype} mRNA (ID:{event_id})")
        record.annotations["molecule_type"] = "mRNA"
        record.annotations["date"] = datetime.now().strftime("%d-%b-%Y").upper()

        cum = 0
        for idx, (s, e) in enumerate(exon_bounds):
            elen = e - s
            rel_s, rel_e = cum, cum + elen
            if is_minus:
                rel_s, rel_e = _mirror(cum + elen), _mirror(cum)
            tx_num = (n_exons - idx) if is_minus else (idx + 1)
            record.features.append(SeqFeature(
                FeatureLocation(rel_s, rel_e, strand=feat_strand),
                type="exon", qualifiers={"number": str(tx_num)}))
            cum += elen

        if cds_ph_list:
            cds_sorted = sorted(cds_ph_list, key=lambda x: x[0])
            cds_rel_parts = []
            for cs, ce, _ in cds_sorted:
                c = 0
                for es, ee in exon_bounds:
                    elen = ee - es
                    if cs >= es and ce <= ee:
                        off = cs - es
                        rs, re_ = c + off, c + off + (ce - cs)
                        if is_minus:
                            rs, re_ = _mirror(c + off + (ce - cs)), _mirror(c + off)
                        cds_rel_parts.append((rs, re_))
                        break
                    elif cs < ee and ce > es:
                        os_ = max(cs, es)
                        oe = min(ce, ee)
                        off = os_ - es
                        rs, re_ = c + off, c + off + (oe - os_)
                        if is_minus:
                            rs, re_ = _mirror(c + off + (oe - os_)), _mirror(c + off)
                        cds_rel_parts.append((rs, re_))
                    c += elen

            if cds_rel_parts:
                cds_rel_parts.sort(key=lambda x: x[0])
                locs = [FeatureLocation(s, e, strand=feat_strand) for s, e in cds_rel_parts]
                cds_loc = locs[0] if len(locs) == 1 else CompoundLocation(locs)

                cds_seq = "".join(mrna_seq[s:e] for s, e in cds_rel_parts)
                translated = str(Seq(cds_seq).translate(to_stop=False)) if cds_seq else ""
                sp = translated.find('*')
                if sp >= 0:
                    translated = translated[:sp + 1]

                quals = {
                    "gene": gene_symbol,
                    "transcript_id": ref_tid,
                    "codon_start": "1",
                    "product": f"{gene_symbol} protein"
                }
                if translated:
                    quals["translation"] = translated

                record.features.append(SeqFeature(cds_loc, type="CDS", qualifiers=quals))
        else:
            record.annotations["comment"] = f"Non-coding RNA transcript for {gene_symbol}"

        try:
            anchors, diff_regions = _extract_event_regions(row, etype)
            event_labels = []
            if etype == "SE":
                if len(anchors) >= 2:
                    event_labels.append(("Upstream_Anchor_Exon", anchors[0]))
                    event_labels.append(("Downstream_Anchor_Exon", anchors[1]))
                for d in diff_regions:
                    lbl = "Skipped_Exon_IN_WT" if _exon_in_ref(d, exon_bounds) else "Skipped_Exon_ALT_FORM"
                    event_labels.append((lbl, d))
            elif etype == "RI":
                if len(anchors) >= 2:
                    event_labels.append(("Upstream_Exon", anchors[0]))
                    event_labels.append(("Downstream_Exon", anchors[1]))
                for d in diff_regions:
                    event_labels.append(("Retained_Intron_REGION", d))
            elif etype == "MXE":
                if len(anchors) >= 2:
                    event_labels.append(("Upstream_Exon", anchors[0]))
                    event_labels.append(("Downstream_Exon", anchors[1]))
                for i, d in enumerate(diff_regions):
                    in_wt = _exon_in_ref(d, exon_bounds)
                    lbl = f"Exon{i + 1}_Candidate_IN_WT" if in_wt else f"Exon{i + 1}_Candidate_ALT_FORM"
                    event_labels.append((lbl, d))
            elif etype in ("A5SS", "A3SS"):
                short = (int(float(row['shortES'])), int(float(row['shortEE'])))
                long_ = (int(float(row['longExonStart_0base'])), int(float(row['longExonEnd'])))
                flank = (int(float(row['flankingES'])), int(float(row['flankingEE'])))

                short_in_wt = _exon_in_ref(short, exon_bounds)
                long_in_wt = _exon_in_ref(long_, exon_bounds)
                s_label = "Short_Exon_IN_WT" if short_in_wt else "Short_Exon_ALT_FORM"
                l_label = "Long_Exon_IN_WT" if long_in_wt else "Long_Exon_ALT_FORM"
                event_labels.append((s_label, short))
                event_labels.append((l_label, long_))
                event_labels.append(("Flanking_Exon", flank))
                for d in diff_regions:
                    event_labels.append(("Alt_SS_DIFF_Region", d))

            for label, (gs, ge) in event_labels:
                rel_parts = []
                c = 0
                found = False
                for es, ee in exon_bounds:
                    elen = ee - es
                    os_ = max(gs, es)
                    oe = min(ge, ee)
                    if os_ < oe:
                        off = os_ - es
                        rs, re_ = c + off, c + off + (oe - os_)
                        if is_minus:
                            rs, re_ = _mirror(re_), _mirror(rs)
                        rel_parts.append((rs, re_))
                        found = True
                    c += elen
                if found and rel_parts:
                    rel_parts.sort(key=lambda x: x[0])
                    merged = [rel_parts[0]]
                    for s, e in rel_parts[1:]:
                        if s <= merged[-1][1]:
                            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
                        else:
                            merged.append((s, e))
                    locs = [FeatureLocation(s, e, strand=feat_strand) for s, e in merged]
                    feat_loc = locs[0] if len(locs) == 1 else CompoundLocation(locs)
                    record.features.append(SeqFeature(
                        feat_loc, type="misc_feature",
                        qualifiers={"label": label, "note": f"event={etype};genomic={chrom}:{gs}-{ge}"}))
                else:
                    existing_note = record.annotations.get("comment", "")
                    add_note = f"{label}: {chrom}:{gs}-{ge} (not in this WT transcript)"
                    record.annotations["comment"] = f"{existing_note}; {add_note}" if existing_note else add_note
        except Exception as e:
            logger.debug(f"EXON-GB misc_feature annotation failed: {e}")

        gb_path = os.path.join(output_dir, f"{locus_name}.gb")
        SeqIO.write(record, gb_path, "genbank")
        return gb_path

    except Exception as e:
        logger.warning(f"EXON-GB failed {ref_tid}: {e}")
        return None


def generate_tran_gb(ref_tid, full_tx_exons, cds_ph_list, genome, chrom, strand,
                     gene_symbol, etype, row, output_dir, job_id, event_id=""):
    if not ref_tid:
        return None
    try:
        all_coords = []
        if full_tx_exons:
            for s, e in full_tx_exons:
                all_coords.extend([s, e])
        if cds_ph_list:
            for cs, ce, _ in cds_ph_list:
                all_coords.extend([cs, ce])
        if not all_coords:
            return None

        gene_start = min(all_coords)
        gene_end = max(all_coords)
        genomic_seq = safe_get_seq(genome, chrom, gene_start, gene_end)
        if not genomic_seq:
            return None

        is_minus = (strand == '-')
        display_seq = str(Seq(genomic_seq).reverse_complement()) if is_minus else genomic_seq
        seq_len = len(display_seq)
        feat_strand = 1

        def _mirror(pos):
            return seq_len - pos

        safe_gene = sanitize_filename(gene_symbol)
        id_prefix = f"{job_id}_" if job_id else ""
        locus_name = f"{id_prefix}{event_id}_Tran-{etype}-{safe_gene}"

        record = SeqRecord(Seq(display_seq), id=ref_tid, name=locus_name[:16],
                           description=f"{gene_symbol} {etype} genomic (ID:{event_id})")
        record.annotations["molecule_type"] = "DNA"
        record.annotations["date"] = datetime.now().strftime("%d-%b-%Y").upper()

        if full_tx_exons:
            n_exons = len(full_tx_exons)
            for idx, (s, e) in enumerate(full_tx_exons):
                cs = max(gene_start, min(gene_end, s))
                ce = max(gene_start, min(gene_end, e))
                if cs >= ce:
                    continue
                rel_s, rel_e = cs - gene_start, ce - gene_start
                if is_minus:
                    rel_s, rel_e = _mirror(rel_e), _mirror(rel_s)
                tx_num = (n_exons - idx) if is_minus else (idx + 1)
                record.features.append(SeqFeature(
                    FeatureLocation(rel_s, rel_e, strand=feat_strand),
                    type="exon", qualifiers={"number": str(tx_num)}))

        if cds_ph_list:
            raw_parts = []
            for gs, ge, _ in sorted(cds_ph_list, key=lambda x: x[0]):
                cs = max(gene_start, min(gene_end, gs))
                ce = max(gene_start, min(gene_end, ge))
                if cs >= ce:
                    continue
                rel_s, rel_e = cs - gene_start, ce - gene_start
                if is_minus:
                    rel_s, rel_e = _mirror(rel_e), _mirror(rel_s)
                raw_parts.append((rel_s, rel_e))
            raw_parts.sort(key=lambda x: x[0])
            merged = []
            for s, e in raw_parts:
                if merged and s <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], e))
                else:
                    merged.append((s, e))
            if merged:
                locs = [FeatureLocation(s, e, strand=feat_strand) for s, e in merged]
                cds_loc = locs[0] if len(locs) == 1 else CompoundLocation(locs)
                record.features.append(SeqFeature(
                    cds_loc, type="CDS",
                    qualifiers={"gene": gene_symbol, "transcript_id": ref_tid,
                                "codon_start": "1", "product": f"{gene_symbol} protein"}))

        try:
            anchors, diff_regions = _extract_event_regions(row, etype)
            event_labels = []
            if etype == "SE":
                if len(anchors) >= 2:
                    event_labels.append(("Upstream_Anchor", anchors[0]))
                    event_labels.append(("Downstream_Anchor", anchors[1]))
                for d in diff_regions:
                    event_labels.append(("Skipped_Exon_DIFF", d))
            elif etype == "RI":
                if len(anchors) >= 2:
                    event_labels.append(("Upstream_Exon", anchors[0]))
                    event_labels.append(("Downstream_Exon", anchors[1]))
                for d in diff_regions:
                    event_labels.append(("Retained_Intron_DIFF", d))
            elif etype == "MXE":
                event_labels = []
                if len(anchors) >= 2:
                    event_labels.append(("Upstream_Exon", anchors[0]))
                    event_labels.append(("Downstream_Exon", anchors[1]))
                for i, d in enumerate(diff_regions):
                    in_wt = _exon_in_ref(d, full_tx_exons)
                    lbl = f"Exon{i + 1}_Candidate_IN_WT" if in_wt else f"Exon{i + 1}_Candidate_ALT_FORM"
                    event_labels.append((lbl, d))
            elif etype in ("A5SS", "A3SS"):
                short = (int(float(row['shortES'])), int(float(row['shortEE'])))
                long_ = (int(float(row['longExonStart_0base'])), int(float(row['longExonEnd'])))
                flank = (int(float(row['flankingES'])), int(float(row['flankingEE'])))

                short_in_wt = _exon_in_ref(short, full_tx_exons)
                long_in_wt = _exon_in_ref(long_, full_tx_exons)
                s_label = "Short_Exon_IN_WT" if short_in_wt else "Short_Exon_ALT_FORM"
                l_label = "Long_Exon_IN_WT" if long_in_wt else "Long_Exon_ALT_FORM"
                event_labels = [(s_label, short), (l_label, long_), ("Flanking_Exon", flank)]
                for d in diff_regions:
                    event_labels.append(("Alt_SS_DIFF", d))
            for label, (gs, ge) in event_labels:
                cs = max(gene_start, min(gene_end, gs))
                ce = max(gene_start, min(gene_end, ge))
                if cs < ce:
                    rel_s, rel_e = cs - gene_start, ce - gene_start
                    if is_minus:
                        rel_s, rel_e = _mirror(rel_e), _mirror(rel_s)
                    record.features.append(SeqFeature(
                        FeatureLocation(rel_s, rel_e, strand=feat_strand),
                        type="misc_feature",
                        qualifiers={"label": label, "note": f"event={etype}"}))
        except Exception as e:
            logger.debug(f"Tran-GB misc_feature annotation failed: {e}")

        gb_path = os.path.join(output_dir, f"{locus_name}.gb")
        SeqIO.write(record, gb_path, "genbank")
        return gb_path

    except Exception as e:
        logger.warning(f"Tran-GB failed {ref_tid}: {e}")
        return None
