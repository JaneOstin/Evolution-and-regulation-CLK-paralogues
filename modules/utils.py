import gzip
import pandas as pd
import numpy as np

from intervaltree import Interval, IntervalTree
from collections import Counter


iupac_codes = {
    frozenset("A"): "A",
    frozenset("C"): "C",
    frozenset("G"): "G",
    frozenset("T"): "T",
    frozenset("AC"): "M",
    frozenset("AG"): "R",
    frozenset("AT"): "W",
    frozenset("CG"): "S",
    frozenset("CT"): "Y",
    frozenset("GT"): "K",
    frozenset("ACG"): "V",
    frozenset("ACT"): "H",
    frozenset("AGT"): "D",
    frozenset("CGT"): "B",
    frozenset("ACG"): "N"
}


iupac_rna_codes = {
    frozenset("A"): "A",
    frozenset("C"): "C",
    frozenset("G"): "G",
    frozenset("U"): "U",
    frozenset("AC"): "M",
    frozenset("AG"): "R",
    frozenset("AU"): "W",
    frozenset("CG"): "S",
    frozenset("CU"): "Y",
    frozenset("GU"): "K",
    frozenset("ACG"): "V",
    frozenset("ACU"): "H",
    frozenset("AGU"): "D",
    frozenset("CGU"): "B",
    frozenset("ACGU"): "N"
}


def make_junction_table(meta, j_path, colnames):
    '''Function for combining files into one table'''
    junctions = pd.DataFrame(columns=colnames + ['file', 'chx'])
    
    for f0 in meta['files']:
        '''Read files in a row, add column names, sample category and merge everything into one table.'''
        file_path = j_path / f"{f0}.gz"
        if file_path.exists():
            f=gzip.open(file_path,'rb')
            ff=pd.read_csv(f, sep = '\t', header = None, names = colnames)
            ff['file'] = f0
            ff['experiment'] = meta.loc[meta['files'] == f0, 'group'].values[0]
            junctions = pd.concat([junctions, ff], ignore_index=True)
    
    junctions = junctions.drop(columns=['staggered_count', 'entropy'])
    
    if junctions.columns[0] == 'junction_id':
        junctions = junctions.groupby([junctions.columns[0], 'experiment']).agg({ # Aggregate by site and type of experiment - stack reads by replicate
            'total_count': 'sum',
            'annotation_status': 'max'
        }).reset_index()
        junctions[['chrom', 'start', 'end', 'strand']] = junctions[junctions.columns[0]].str.split('_', expand=True) # Split the name of the reads in the first column
    else:
        junctions = junctions.groupby([junctions.columns[0], 'experiment']).agg({
            'total_count': 'sum'
        }).reset_index()
        junctions[['chrom', 'site', 'strand']] = junctions[junctions.columns[0]].str.split('_', expand=True) # Split the name of the reads in the first column
    
    return junctions


def split_rows(row):
    '''A function that checks that S6 has only the end or the beginning of a reid, or both versions of reids '''
    rows = []
    if pd.notna(row['site_id']):
        new_row = row.copy()
        new_row['total_count_site'] = new_row['total_count_start']
        rows.append(new_row)
    if pd.notna(row['site_id_end']):
        new_row = row.copy()
        new_row['site_id'] = new_row['site_id_end']
        new_row['total_count_site'] = new_row['total_count_end']
        new_row['site'] = new_row['site_end']
        rows.append(new_row)
    if pd.notna(row['site_id']) != True and pd.notna(row['site_id_end']) != True:
        rows.append(row.copy())
    return rows


def get_cassete_exons_from_ipsa(df):
    '''Function for finding cassette exons: searches for intersection of intervals within a split reid'''
    # Turn a DataFrame into a list of tuples
    intervals = []
    for row in df.itertuples():
        intervals.append((row.Index, row.start, row.end, row.chrom, row.experiment))
    
    # Build IntervalTree
    itree = IntervalTree()
    for (id_val, start, end, chrom, experiment) in intervals:
        itree.add(Interval(int(start), int(end), data=(int(id_val), chrom, experiment)))
    
    # For each segment, find who lies inside it
    results = {}
    for (id_val, start, end, chrom, experiment) in intervals:
        contained_intervals = itree.envelop(int(start), int(end))
        contained_ids = [iv.data[0] for iv in contained_intervals if iv.data[1] == chrom and iv.data[2] == experiment]
        contained_ids.remove(int(id_val))
        if len(contained_ids) > 1:
            results[id_val] = contained_ids
    
    cassete_exon = []
    for big_interval, intervals in results.items():
        row_big = df.loc[big_interval]
        exclusion_count = row_big.total_count
        chromosome = row_big.chrom
        strand_ = row_big['strand']
        exp = row_big.experiment
        start_pos = row_big.start
        end_pos = row_big.end
        start_exon = []
        end_exon = []
        i1 = []
        i2 = []

        # Collect end/beginning of cassette exon and read counts
        for idx in intervals:
            row = df.loc[idx]
            start, end, counts = row.start, row.end, row.total_count
            if start == start_pos:
                start_exon.append(end)
                i1.append(counts)
            elif end == end_pos:
                end_exon.append(start)
                i2.append(counts)
        
        # Pair up all the variations
        cassete_exon.extend(
            {"start": start_pos, "end": end_pos, "start_exon": s, "end_exon": e, "i1": c1, "i2": c2, "exclusion_count": exclusion_count, "chrom": chromosome, "strand": strand_, "experiment": exp}
            for s, c1 in zip(start_exon, i1)
            for e, c2 in zip(end_exon, i2)
        )

    cassete_exon = pd.DataFrame(cassete_exon)
    cassete_exon["psi"] = (2*cassete_exon["exclusion_count"]) / (cassete_exon["i1"] + cassete_exon["i2"] + 2*cassete_exon["exclusion_count"]) # Count PSI
    return pd.DataFrame(cassete_exon)


def make_track(ce, colorx = "103,137,33"):
    '''Function for a track at UCSC genome browser'''
    ce = ce.copy()

    # Calculation  score == dPSI
    ce["score"] = np.round(ce["-chx_psi"] - ce["+chx_psi"], 2)
    
    # Filling in BED fields
    ce["chromStart"] = ce["start"] - 2
    ce["chromEnd"] = ce["end"] + 2
    ce["thickStart"] = ce["start"] - 2
    ce["thickEnd"] = ce["end"] + 2
    ce["blockCount"] = 3
    ce["itemRgb"] = colorx
    ce["blockSizes"] = ce.apply(lambda row: f"2,{str(int(row['end_exon'] - row['start_exon']))},2", axis=1)
    ce["blockStarts"] = ce.apply(lambda row: f"0,{str(int(row['start_exon'] - row['start']))},{str(int(row['end'] - row['start'] + 2))}", axis=1)
    
    # Formation of the name field
    ce["name"] = ce.apply(lambda row: (
        f"dpsi={row['score']},psi_chx={np.round(row['+chx_psi'], 2)},psi_-chx={np.round(row['-chx_psi'], 2)},"
        f"excl_chx={row['+chx_exclusion_count']},excl_-chx={row['-chx_exclusion_count']},"
        f"i1_chx={row['+chx_i1']},i1_ctrl={row['-chx_i1']},"
        f"i2_chx={row['+chx_i2']},i2_ctrl={row['-chx_i2']}"
    ), axis=1)

    # Selecting the right columns
    ce_bed = ce[["chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"]]
    return ce_bed


def select_exon(row, intervals):
    '''Selects an exon label for a given row based on genomic intervals'''
    if pd.isna(row['exon_number']):
        for key, value in intervals.items():
            if key != row['gene']:
                pass
            else:
                for intron_number, intervals in value.items():
                    if row['Start'] >= intervals[0] and row['End'] <= intervals[1]:
                        return str(intron_number) + '_poison'
    else:
        return row['exon_number']


def conservation_score(column):
    '''Calculates the degree of conservativity of the alignment column'''
    bases = [b for b in column if b != '-']
    if not bases:
        return 0.0
    most_common = Counter(bases).most_common(1)[0][1]
    return most_common / len(bases)


def is_fully_covered(column):
    '''Check if we have any gaps'''
    return '-' not in column


def find_fully_conserved_windows(alignment, window=10, threshold=0.9):
    '''Finds sections in a multiple alignment that satisfy the condition'''
    alignment_length = alignment.get_alignment_length()
    conserved_windows = []

    for i in range(alignment_length - window + 1):
        window_cols = [alignment[:, i + j] for j in range(window)]

        scores = [conservation_score(col) for col in window_cols]
        no_gaps = all(is_fully_covered(col) for col in window_cols)

        if (sum(scores) / window) >= threshold and no_gaps:
            conserved_windows.append((i + 1, i + window))

    return conserved_windows


def get_consensus_motif(sequences):
    '''Function for obtaining the sequence motif '''
    motif = ""
    for pos in zip(*sequences):
        nucs = frozenset(pos)
        symbol = iupac_codes.get(nucs, "N")
        motif += symbol
    return motif


def get_rna_consensus_motif(sequences):
    '''Function for obtaining the sequence RNA motif '''
    motif = ""
    for pos in zip(*sequences):
        nucs = frozenset(pos)
        symbol = iupac_rna_codes.get(nucs, "N")
        motif += symbol
    return motif
