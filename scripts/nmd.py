#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""PTC genomic mapping and NMD sensitivity classification."""


def _map_ptc_to_genomic(cds_ph_list, strand, ptc_nt_pos):
    if not cds_ph_list or ptc_nt_pos <= 0:
        return None

    if strand == '+':
        ordered_cds = sorted(cds_ph_list, key=lambda x: x[0])
    else:
        ordered_cds = sorted(cds_ph_list, key=lambda x: x[0], reverse=True)

    cum = 0
    for cs, ce, ph in ordered_cds:
        seg_len = ce - cs
        if cum + seg_len >= ptc_nt_pos:
            offset = ptc_nt_pos - cum - 1
            if strand == '+':
                return cs + offset
            else:
                return ce - offset - 1
        cum += seg_len

    return None


def _classify_nmd_by_ptc_position(ptc_genomic_pos, full_tx_exons, strand, start_codon_pos):
    if not full_tx_exons:
        return "NMD_Target"

    sorted_exons = sorted(full_tx_exons, key=lambda x: x[0])
    n_exons = len(sorted_exons)

    if n_exons == 1:
        return "NMD_Escape_3"

    ptc_exon_idx = -1
    for i, (es, ee) in enumerate(sorted_exons):
        if es <= ptc_genomic_pos < ee:
            ptc_exon_idx = i
            break

    if ptc_exon_idx == -1:
        return "NMD_Target"

    if start_codon_pos is not None:
        if strand == '+':
            dist_from_start = ptc_genomic_pos - start_codon_pos
        else:
            dist_from_start = start_codon_pos - ptc_genomic_pos
        if dist_from_start <= 150:
            return "NMD_Escape_1"

    if ptc_exon_idx == n_exons - 1:
        return "NMD_Escape_3"

    if ptc_exon_idx == n_exons - 2:
        if strand == '+':
            junction_pos = sorted_exons[-2][1]
            dist_to_downstream_ej = junction_pos - ptc_genomic_pos
        else:
            junction_pos = sorted_exons[-2][0]
            dist_to_downstream_ej = ptc_genomic_pos - junction_pos
        if dist_to_downstream_ej <= 50:
            return "NMD_Escape_3"

    if strand == '+':
        downstream_junction_pos = sorted_exons[ptc_exon_idx][1]
        dist_to_downstream_ej = downstream_junction_pos - ptc_genomic_pos
    else:
        downstream_junction_pos = sorted_exons[ptc_exon_idx][0]
        dist_to_downstream_ej = ptc_genomic_pos - downstream_junction_pos

    if dist_to_downstream_ej > 400:
        return "NMD_Escape_2"

    return "NMD_Target"


def judge_nmd_status(wt_prot, alt_prot, is_frameshift,
                     alt_cds_ph_truncated, full_tx_exons, strand, start_codon_pos=None):
    if not wt_prot or not alt_prot:
        return "NO_PROTEIN"

    if not is_frameshift:
        return "Express"

    alt_star_pos = alt_prot.find('*')
    if alt_star_pos == -1:
        return "NO_STOP_CODON"

    ptc_nt_pos = (alt_star_pos + 1) * 3
    ptc_genomic_pos = _map_ptc_to_genomic(alt_cds_ph_truncated, strand, ptc_nt_pos)

    if ptc_genomic_pos is None:
        return "NMD_Target"

    return _classify_nmd_by_ptc_position(ptc_genomic_pos, full_tx_exons, strand, start_codon_pos)


def _get_start_codon_position(ref_cds_ph_trunc, strand):
    if not ref_cds_ph_trunc:
        return None
    if strand == '+':
        return ref_cds_ph_trunc[0][0]
    else:
        return ref_cds_ph_trunc[-1][1]
