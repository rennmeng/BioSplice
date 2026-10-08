#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Derive the alternative-isoform CDS segments from the reference CDS.

Both input and output are lists of half-open segments
[(cds_start_0, cds_end_0, phase), ...].
"""


def compute_alt_from_ref(ref_cds_ph, row, etype, strand):
    ref_sorted = sorted(ref_cds_ph, key=lambda x: x[0])

    if etype == "RI":
        up_end = int(float(row['upstreamEE']))
        dn_start = int(float(row['downstreamES']))
        if up_end >= dn_start:
            return ref_sorted
        new_cds, inserted = [], False
        for cs, ce, ph in ref_sorted:
            if not inserted and cs >= dn_start:
                new_cds.append((up_end, dn_start, 0))
                inserted = True
            if not inserted and cs < up_end < ce:
                new_cds.append((cs, up_end, ph))
                new_cds.append((up_end, dn_start, 0))
                inserted = True
                new_cds.append((dn_start, ce, 0))
                continue
            new_cds.append((cs, ce, ph))
        if not inserted:
            new_cds.append((up_end, dn_start, 0))
        return sorted(new_cds, key=lambda x: x[0])

    elif etype == "SE":
        skip_s = int(float(row['exonStart_0base']))
        skip_e = int(float(row['exonEnd']))
        new_cds = []
        for cs, ce, ph in ref_sorted:
            os_ = max(cs, skip_s)
            oe = min(ce, skip_e)
            if os_ >= oe:
                new_cds.append((cs, ce, ph))
            else:
                if cs < os_:
                    new_cds.append((cs, os_, ph))
                if ce > oe:
                    new_cds.append((oe, ce, ph))
        return sorted(new_cds, key=lambda x: x[0]) if new_cds else ref_sorted

    elif etype == "MXE":
        # ==================== MXE handling v2 (exclusive-exon fix) ====================
        e1_s = int(float(row['1stExonStart_0base']))
        e1_e = int(float(row['1stExonEnd']))
        e2_s = int(float(row['2ndExonStart_0base']))
        e2_e = int(float(row['2ndExonEnd']))

        def _ov_len(rs, re_):
            tot = 0
            for cs, ce, _ in ref_sorted:
                o_s, o_e = max(cs, rs), min(ce, re_)
                if o_s < o_e:
                    tot += (o_e - o_s)
            return tot

        ov1 = _ov_len(e1_s, e1_e)
        ov2 = _ov_len(e2_s, e2_e)

        if ov1 == 0 and ov2 == 0:
            if strand == '+':
                before_len = sum(ce - cs for cs, ce, _ in ref_sorted if ce <= e2_s)
            else:
                before_len = sum(ce - cs for cs, ce, _ in ref_sorted if cs >= e2_e)
            new_cds = list(ref_sorted)
            new_cds.append((e2_s, e2_e, before_len % 3))
            return sorted(new_cds, key=lambda x: x[0])

        if ov1 >= ov2:
            ref_s, ref_e = e1_s, e1_e
            alt_s, alt_e = e2_s, e2_e
        else:
            ref_s, ref_e = e2_s, e2_e
            alt_s, alt_e = e1_s, e1_e

        kept = []
        for cs, ce, ph in ref_sorted:
            o_s, o_e = max(cs, ref_s), min(ce, ref_e)
            if o_s >= o_e:
                kept.append((cs, ce, ph))
                continue
            if cs < o_s:
                kept.append((cs, o_s, ph))
            if ce > o_e:
                kept.append((o_e, ce, ph))

        has_alt_in_ref = any(max(cs, alt_s) < min(ce, alt_e) for cs, ce, _ in kept)

        if not has_alt_in_ref:
            if strand == '+':
                before_len = sum(ce - cs for cs, ce, _ in kept if ce <= alt_s)
            else:
                before_len = sum(ce - cs for cs, ce, _ in kept if cs >= alt_e)
            kept.append((alt_s, alt_e, before_len % 3))

        kept.sort(key=lambda x: x[0])

        merged = []
        for s, e, ph in kept:
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e), merged[-1][2])
            else:
                merged.append((s, e, ph))
        return merged

    elif etype in ("A5SS", "A3SS"):
        short_s = int(float(row['shortES']))
        short_e = int(float(row['shortEE']))
        long_s = int(float(row['longExonStart_0base']))
        long_e = int(float(row['longExonEnd']))

        def _overlap_len(rs, re_):
            tot = 0
            for cs, ce, _ in ref_sorted:
                os_ = max(cs, rs)
                oe = min(ce, re_)
                if os_ < oe:
                    tot += (oe - os_)
            return tot

        ov_short = _overlap_len(short_s, short_e)
        ov_long = _overlap_len(long_s, long_e)

        if ov_short == 0 and ov_long == 0:
            return ref_sorted

        ref_uses_short = (ov_short >= ov_long)
        if ref_uses_short:
            ref_es, ref_ee = short_s, short_e
            alt_es, alt_ee = long_s, long_e
        else:
            ref_es, ref_ee = long_s, long_e
            alt_es, alt_ee = short_s, short_e

        variable_is_start = (etype == "A5SS")
        TOL = 2

        new_cds = []
        handled = False

        for cs, ce, ph in ref_sorted:
            os_ = max(cs, ref_es)
            oe = min(ce, ref_ee)

            if os_ >= oe:
                new_cds.append((cs, ce, ph))
                continue

            if handled:
                new_cds.append((cs, ce, ph))
                continue

            if cs < ref_es:
                new_cds.append((cs, ref_es, ph))
            if ce > ref_ee:
                new_cds.append((ref_ee, ce, ph))

            if variable_is_start:
                if os_ <= ref_es + TOL:
                    new_s = alt_es
                else:
                    offset_from_fixed_end = ref_ee - os_
                    new_s = alt_ee - offset_from_fixed_end

                offset_end_from_fixed = ref_ee - oe
                new_e = alt_ee - offset_end_from_fixed
            else:
                if oe >= ref_ee - TOL:
                    new_e = alt_ee
                else:
                    offset_end_from_fixed = oe - ref_es
                    new_e = alt_es + offset_end_from_fixed

                offset_from_fixed = os_ - ref_es
                new_s = alt_es + offset_from_fixed

            lo, hi = min(alt_es, alt_ee), max(alt_es, alt_ee)
            new_s = max(new_s, lo)
            new_e = min(new_e, hi)

            if new_s < new_e:
                new_cds.append((new_s, new_e, ph))

            handled = True

        if not new_cds:
            return ref_sorted

        new_cds = sorted(new_cds, key=lambda x: x[0])
        merged = []
        for s, e, ph in new_cds:
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e), merged[-1][2])
            else:
                merged.append((s, e, ph))
        return merged

    return ref_sorted
