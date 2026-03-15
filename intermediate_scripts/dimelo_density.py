#!/usr/bin/env python

import argparse
import csv
import itertools
import multiprocessing as mp
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO


Region = Tuple[int, int]
RegionDict = Dict[str, List[Region]]


def bed_to_region_dict(bed_path: str) -> RegionDict:
    """
    Read a BED file and return regions per chromosome.

    Parameters
    ----------
    bed_path : str
        Path to BED file with at least 3 columns: chrom, start, end.

    Returns
    -------
    dict
        {chrom: [(start, end), ...]}
    """
    regions: RegionDict = {}
    with open(bed_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            regions.setdefault(chrom, []).append((start, end))
    return regions


def load_reference_lengths(ref_fasta_path: str) -> Dict[str, int]:
    """
    Load reference FASTA and return chromosome lengths.

    Parameters
    ----------
    ref_fasta_path : str
        Path to reference genome FASTA.

    Returns
    -------
    dict
        {chrom: length}
    """
    lengths: Dict[str, int] = {}
    with open(ref_fasta_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            lengths[record.id] = len(record.seq)
    return lengths


def mod_subset_producing_step(
    mod_no_dash: np.ndarray,
    alignment_dash: str,
    target_start_no_dash: int,
    target_end_no_dash: int,
) -> np.ndarray:
    """
    Subset mod_no_dash to the segment corresponding to [target_start_no_dash, target_end_no_dash)
    in the *no-dash* coordinate system, while alignment_dash still contains '-' for indels.
    """
    # mask: True for non-dash positions
    mask = [char != "-" for char in alignment_dash]

    # cumulative counts only for True entries
    cumulative_counts = list(itertools.accumulate(mask))

    # map dashed positions to no-dash indices (or '-' for dashes)
    indexes = [
        count - 1 if is_non_dash else "-"
        for count, is_non_dash in zip(cumulative_counts, mask)
    ]

    target_start_dash = indexes.index(target_start_no_dash)
    try:
        target_end_dash = indexes.index(target_end_no_dash)
    except ValueError:
        target_end_dash = len(indexes) - 1

    # Get dashed alignment for pre-subset and subset
    alignment_dash_sequence_pre_subset = alignment_dash[0:target_start_dash]
    alignment_dash_sequence_subset = alignment_dash[target_start_dash:target_end_dash]

    # Remove dashes to get no-dash sequences
    alignment_no_dash_sequence_pre_subset = alignment_dash_sequence_pre_subset.replace(
        "-", ""
    )
    alignment_no_dash_sequence_subset = alignment_dash_sequence_subset.replace("-", "")

    subset_no_dash_start = len(alignment_no_dash_sequence_pre_subset)
    subset_no_dash_end = subset_no_dash_start + len(alignment_no_dash_sequence_subset)

    # Subset mod array
    mod_subset = mod_no_dash[subset_no_dash_start:subset_no_dash_end]
    return mod_subset


def tile_regions(region_dict: RegionDict, window_size: int) -> RegionDict:
    """
    Split each (start, end) region into fixed-size windows.

    For each region [start, end), generate windows:
        [start, start + window),
        [start + window, start + 2*window),
        ...
    with the last window clipped to <= end.

    Parameters
    ----------
    region_dict : dict
        {chrom: [(start, end), ...]}
    window_size : int
        Window size in bp (e.g. 1000).

    Returns
    -------
    dict
        {chrom: [(win_start, win_end), ...]} tiled windows.
    """
    tiled: RegionDict = {}
    for chrom, regions in region_dict.items():
        tiled_regions: List[Region] = []
        for start, end in regions:
            pos = start
            while pos < end:
                win_end = min(pos + window_size, end)
                tiled_regions.append((pos, win_end))
                pos += window_size
        if tiled_regions:
            tiled[chrom] = tiled_regions
    return tiled


def _compute_for_chromosome(
    args: Tuple[str, List[Region], str, str, float]
) -> List[Tuple[str, str, float, float]]:
    """
    Internal worker: compute region densities for a single chromosome.

    Returns a list of rows:
        (chrom, "[start, end]", density, coverage)
    """
    chrom, regions, bam_path, mod_tag, filtering_val = args
    bamfile = pysam.AlignmentFile(bam_path, "rb")

    rows: List[Tuple[str, str, float, float]] = []

    for region_start_index, region_end_index in regions:
        region_density_vals: List[float] = []
        region_base = 0

        for read in bamfile.fetch(chrom, region_start_index, region_end_index):
            read_start_position = read.reference_start
            read_end_position = read.reference_end

            # Sequence with indels explicitly represented as '-'
            sequence = read.get_aligned_pairs(matches_only=False, with_seq=True)

            read_sequence_insertion_included = ""
            genomic_alignment_sequence_deletion_mistach_included = ""

            for ref_pos, read_pos, base in sequence:
                if ref_pos is None:
                    # insertion
                    read_sequence_insertion_included += "-"
                elif read_pos is None:
                    # deletion or mismatch on reference
                    genomic_alignment_sequence_deletion_mistach_included += "-"
                else:
                    base = base if base is not None else "N"
                    read_sequence_insertion_included += base
                    genomic_alignment_sequence_deletion_mistach_included += base

            read_sequence_insertion_included = read_sequence_insertion_included.upper()
            genomic_alignment_sequence_deletion_mistach_included = (
                genomic_alignment_sequence_deletion_mistach_included.upper()
            )

            genomic_alignment_mask = np.array(
                [char != "-" for char in genomic_alignment_sequence_deletion_mistach_included]
            )

            insertions = read_sequence_insertion_included.count("-")
            no_insertion_no_deletion_sequence_length = len(
                read_sequence_insertion_included
            )

            # Build mod_score array (same length as alignment with deletions)
            try:
                mod = read.modified_bases_forward
            except AttributeError:
                # No modified base info
                continue

            mod_score = np.zeros(
                len(genomic_alignment_sequence_deletion_mistach_included),
                dtype=float,
            )

            try:
                if mod_tag == "A":
                    # 6mA on A
                    for indices, values in mod[("A", 0, "a")]:
                        mod_score[indices] = values
                elif mod_tag == "CG":
                    # 5mC on C (CpG context)
                    for indices, values in mod[("C", 0, "m")]:
                        mod_score[indices] = values
                else:
                    raise ValueError(f"Unsupported mod_tag: {mod_tag}")
            except KeyError:
                # read has no entries for this mod
                continue

            # Remove deletions
            mod_score = mod_score[genomic_alignment_mask]

            # Reverse orientation if read is reverse
            if read.is_reverse:
                mod_score = mod_score[::-1]

            # Apply filtering threshold
            mod_score = np.where(mod_score < filtering_val, 0, mod_score)

            # Decide overlap between region and read
            if (region_end_index - region_start_index) > (
                read_end_position - read_start_position
            ):
                # region longer than read
                if (region_end_index >= read_end_position) and (
                    region_start_index <= read_start_position
                ):
                    # read fully inside region
                    mod_start = 0
                    mod_end = no_insertion_no_deletion_sequence_length
                elif (region_end_index < read_end_position) and (
                    region_start_index > read_start_position
                ):
                    # read covers later part of region
                    mod_start = 0
                    mod_end = (
                        no_insertion_no_deletion_sequence_length
                        - read_end_position
                        - region_end_index
                    )
                elif (region_end_index > read_end_position) and (
                    region_start_index > read_start_position
                ):
                    # region covers starting part of read
                    mod_start = region_start_index - read_start_position
                    mod_end = no_insertion_no_deletion_sequence_length
                else:
                    # no valid overlap scenario
                    continue
            else:
                # read longer than or equal to region
                if (read_start_position <= region_start_index) and (
                    read_end_position >= region_end_index
                ):
                    # region inside read
                    mod_start = region_start_index - read_start_position
                    mod_end = region_end_index - read_start_position
                elif (read_end_position < region_end_index) and (
                    read_end_position > region_start_index
                ):
                    # region covers end of read
                    mod_start = region_start_index - read_start_position
                    mod_end = no_insertion_no_deletion_sequence_length
                elif (read_start_position > region_start_index) and (
                    read_start_position < region_end_index
                ):
                    # region covers beginning of read
                    mod_start = 0
                    mod_end = region_end_index - read_start_position
                else:
                    continue

            # Guard against impossible positions
            if (region_start_index - read_start_position) > (
                no_insertion_no_deletion_sequence_length - insertions
            ):
                continue

            try:
                trimmed_mod_score = mod_subset_producing_step(
                    mod_score, read_sequence_insertion_included, mod_start, mod_end
                )
            except ValueError:
                # mis-mapping in index search
                continue

            # Number of bases in this read that contribute to region coverage
            region_base += (mod_end - mod_start)

            # Remove zeros (non-modified or below threshold)
            mod_no_zeros = trimmed_mod_score[trimmed_mod_score != 0]
            m_mod_tag = len(mod_no_zeros)

            # Count motif occurrences in this read segment
            segment_seq = read_sequence_insertion_included[mod_start:mod_end]
            motif = "A" if mod_tag == "A" else "CG"
            total_mod_tag = segment_seq.count(motif)

            try:
                read_density = m_mod_tag / total_mod_tag
            except ZeroDivisionError:
                read_density = 0.0

            region_density_vals.append(read_density)

        # Aggregate per-region
        if region_density_vals:
            region_density_average = float(sum(region_density_vals) / len(region_density_vals))
        else:
            region_density_average = 0.0

        region_length = max(1, region_end_index - region_start_index)
        coverage_fraction = region_base / region_length

        # Column 2 format: [start, end]
        coord_str = f"[{region_start_index}, {region_end_index}]"
        rows.append((chrom, coord_str, region_density_average, coverage_fraction))

    bamfile.close()
    return rows


def run_region_density(
    bam_path: str,
    bed_path: str,
    ref_fasta_path: str,
    mod_tag: str,
    filtering_val: float,
    out_csv_path: str,
    threads: int = 1,
    window_size: int = 1000,
) -> pd.DataFrame:
    """
    High-level entry point.

    Output file columns (no header):
        chrom    [start, end]    density    coverage
    """
    # 1) Load regions and reference lengths
    regions = bed_to_region_dict(bed_path)
    ref_lengths = load_reference_lengths(ref_fasta_path)

    # 2) Sanity-filter regions: valid lengths, non-empty
    filtered_regions: RegionDict = {}
    for chrom, reg_list in regions.items():
        if chrom not in ref_lengths:
            # Skip chromosomes not present in reference
            continue
        chrom_len = ref_lengths[chrom]
        valid: List[Region] = []
        for start, end in reg_list:
            if end <= start:
                continue
            if start < 0:
                start = 0
            if end > chrom_len:
                end = chrom_len
            if end > start:
                valid.append((start, end))
        if valid:
            filtered_regions[chrom] = valid

    # 3) Tile into fixed-size windows (e.g. 1000 bp)
    tiled_regions = tile_regions(filtered_regions, window_size=window_size)

    tasks = [
        (chrom, reg_list, bam_path, mod_tag, filtering_val)
        for chrom, reg_list in tiled_regions.items()
    ]

    # 4) Compute densities (optionally in parallel)
    if threads > 1 and len(tasks) > 1:
        with mp.Pool(processes=threads) as pool:
            results = pool.map(_compute_for_chromosome, tasks)
        rows = list(itertools.chain.from_iterable(results))
    else:
        rows: List[Tuple[str, str, float, float]] = []
        for t in tasks:
            rows.extend(_compute_for_chromosome(t))

    # 5) Build DataFrame and write TSV (no header)
    df = pd.DataFrame(rows, columns=["chrom", "coord", "density", "coverage"])
    df.to_csv(out_csv_path, sep="\t", header=False, index=False, quoting=csv.QUOTE_MINIMAL)
    return df


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute region-level modified base density from a DiMeLo-style BAM and BED, "
                    "tiling regions into fixed-size windows (default 1000 bp)."
    )
    parser.add_argument(
        "--bam",
        required=True,
        help="Input BAM file with modified base tags.",
    )
    parser.add_argument(
        "--bed",
        required=True,
        help="BED file with regions of interest (chrom, start, end).",
    )
    parser.add_argument(
        "--ref",
        required=True,
        help="Reference genome FASTA file.",
    )
    parser.add_argument(
        "--mod-tag",
        required=True,
        choices=["A", "CG"],
        help="Modification to analyze: 'A' for 6mA, 'CG' for CpG 5mC.",
    )
    parser.add_argument(
        "--threshold",
        required=True,
        type=float,
        help="Filtering threshold for mod score; values below this are set to 0.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV file path. Columns: chrom, [start, end], density, coverage (no header).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of worker processes to use (default: 1).",
    )
    parser.add_argument(
        "--window-size",
        type=int,
        default=1000,
        help="Window size in bp for tiling regions (default: 1000).",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    run_region_density(
        bam_path=args.bam,
        bed_path=args.bed,
        ref_fasta_path=args.ref,
        mod_tag=args.mod_tag,
        filtering_val=args.threshold,
        out_csv_path=args.output,
        threads=args.threads,
        window_size=args.window_size,
    )


if __name__ == "__main__":
    main()
