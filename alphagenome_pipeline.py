# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "alphagenome",
#   "pandas",
#   "pyarrow",
#   "requests",
#   "matplotlib",
#   "gtfparse",
# ]
# ///

import importlib.util
import argparse
import os
import subprocess
import sys
import re
import warnings
# Silence harmless warnings raised inside third-party libraries that only confuse
# non-technical users: pandas FutureWarnings (AlphaGenome's groupby calls) and the
# matplotlib "identical low and high ylims" UserWarning (flat/constant plot tracks).
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', message='Attempting to set identical low and high ylims.*')
warnings.filterwarnings('ignore', message='Transforming to str index.*')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt





if sys.version_info < (3, 10):
    print(f"ERROR: Python 3.10 or higher is required.\nYou are running Python {sys.version_info.major}.{sys.version_info.minor}.")

    print("This script requires Python 3.10 or higher to run. Please upgrade your Python version and try again. https://www.python.org/downloads/")
    sys.exit(1)

genome = None
dna_client = None
variant_scorers = None
plot_components = None
junction_data = None
gene_annotation = None



def ensure_packages():
    
    packages = ['alphagenome', 'requests', 'pandas']
    for package in packages:
        if importlib.util.find_spec(package) is None:
            print(f"{package} is not installed, installing now")
            try:
                subprocess.check_call([sys.executable, '-m', 'pip', 'install', package])
                print(f"{package} has been installed")
            except subprocess.CalledProcessError:
                print(f"ERROR: Failed to install {package}")
                print(f"Please install {package} manually using pip and try again.")
                sys.exit(1)

CANONICAL_COLUMNS = {
    'id':               'ID',
    'chrom':            'CHROM',
    'pos':              'POS',
    'ref':              'REF',
    'alt':              'ALT',
    'gene':             'Gene',
    'strandedness':     'Strandedness',
    'deletion_length':  'deletion_length',
    'window_bp':        'WINDOW_BP',
    'centre_pos':       'CENTRE_POS',
    'cdna':             'cDNA',
    'hgvs':             'HGVS',
    'phenotype':        'Phenotype',
}

def load_variants(path):
    import pandas as pd
    if path.endswith('.tsv') or path.endswith('.csv'):
        try:
            # encoding='utf-8-sig' strips the BOM that Excel-for-Windows prepends to UTF-8 files,
            # which would otherwise mangle the first column name (e.g. 'ID' -> '﻿ID').
            # keep_default_na=False prevents empty REF/ALT cells from becoming NaN, which breaks deletion handling.
            if path.endswith('.tsv'):
                df = pd.read_csv(path, sep='\t', keep_default_na=False, encoding='utf-8-sig')
            else:
                # sep=None + engine='python' auto-detects comma vs semicolon (European Excel exports)
                df = pd.read_csv(path, sep=None, engine='python', keep_default_na=False, encoding='utf-8-sig')

            # normalise column names to canonical casing (case-insensitive input)
            df.columns = [CANONICAL_COLUMNS.get(c.strip().lower(), c.strip()) for c in df.columns]

            # Rows shorter than the header (ragged rows, common in Excel exports) get their
            # missing trailing cells filled with NaN. Normalise those to empty strings so blank-row
            # detection and downstream .get() lookups behave consistently.
            df = df.fillna('')

            # check if required columns are present
            required_columns = {'ID', 'CHROM', 'POS', 'Gene', 'REF', 'ALT'}
            missing_columns = required_columns - set(df.columns)
            if missing_columns:
                print(f"ERROR: Missing required columns: {', '.join(missing_columns)}")
                print("Please ensure the input file contains the following columns: ID, CHROM, POS, Gene, REF, ALT")
                sys.exit(1)

            # strip whitespace from all string cells (common when editing in Excel)
            for col in df.columns:
                if df[col].dtype == object:
                    df[col] = df[col].str.strip()

            # drop fully-blank rows (Excel often appends many empty trailing rows)
            blank_mask = df.apply(lambda r: all(str(x).strip() == '' for x in r), axis=1)
            if blank_mask.any():
                print(f"NOTE: Ignoring {int(blank_mask.sum())} blank row(s) in the input file.")
                df = df[~blank_mask].reset_index(drop=True)

            # empty file / header only
            if len(df) == 0:
                print("ERROR: The input file contains no variant rows (only a header or completely empty).")
                print("Please add at least one variant row and try again.")
                sys.exit(1)

            # cast required columns to correct types after renaming
            for col, typ in {'ID': str, 'CHROM': str, 'Gene': str, 'REF': str, 'ALT': str}.items():
                if col in df.columns:
                    df[col] = df[col].astype(typ)
            try:
                df['POS'] = df['POS'].astype(int)
            except ValueError:
                print("ERROR: The POS column contains values that are not whole numbers.")
                print("Please ensure every POS cell is a genomic position (e.g. 11022418) with no letters, commas, or decimals.")
                sys.exit(1)

            # duplicate IDs would cause the wrong row to be looked up during plotting
            dup_ids = df['ID'][df['ID'].duplicated()].unique()
            if len(dup_ids) > 0:
                print(f"ERROR: Duplicate variant IDs found: {', '.join(dup_ids)}")
                print("Each variant must have a unique value in the ID column.")
                sys.exit(1)

            # empty IDs break per-variant output folders and lookups
            if (df['ID'].str.len() == 0).any():
                print("ERROR: One or more rows have an empty ID.")
                print("Every variant must have a unique, non-empty value in the ID column.")
                sys.exit(1)

            # Strandedness is optional (the gene annotation is authoritative). If the column is
            # present, any value that is filled in must be '+' or '-'; blank cells are allowed
            # and simply mean "use the annotation".
            if 'Strandedness' in df.columns:
                df['Strandedness'] = df['Strandedness'].astype(str)
                bad_strand = df[~df['Strandedness'].isin(['+', '-', ''])]
                if len(bad_strand) > 0:
                    ids = ', '.join(bad_strand['ID'].astype(str).tolist())
                    print(f"ERROR: Invalid Strandedness for variant(s): {ids}")
                    print("The Strandedness column may only contain '+' or '-' (or be left blank). You can also remove the column entirely.")
                    sys.exit(1)

            # normalise chromosome format: accept '1', 'Chr1', 'CHR1' etc. -> 'chr1'
            def _normalise_chrom(c):
                c = str(c).strip()
                if not c.lower().startswith('chr'):
                    return f'chr{c}'
                return 'chr' + c[3:]   # lowercase the prefix, keep the rest (handles chrX, chrM)
            original = df['CHROM'].copy()
            df['CHROM'] = df['CHROM'].apply(_normalise_chrom)
            changed = original[original != df['CHROM']]
            if len(changed) > 0:
                print(f"NOTE: Reformatted {len(changed)} chromosome value(s) to the 'chr...' style required by AlphaGenome.")

            return df

        except SystemExit:
            raise
        except Exception as e:
            print(f"ERROR: Failed to load the input file: {e}")
            print("Please ensure the file is in the correct format and try again.")
            sys.exit(1)
    else:
        print("ERROR: Unsupported file format. Please provide a .tsv or .csv file.")
        sys.exit(1)

def load_ontology_terms(path, organism='human'):
    default_ontologies = {
        'human': [
            'UBERON:0000305',   # neural tissue
            'CL:0000540',       # neural cells
            'UBERON:0000955',   # brain tissue
            'CL:0000100',       # motor neuron
            'UBERON:0001870',   # frontal cortex
            'UBERON:0001871',   # temporal lobe
            'CL:0000188',       # skeletal muscle
        ],
        'mouse': [
            'CL:0002608',   # hippocampal neuron
            'CL:0002609',   # neuron of cerebral cortex
            'CL:0002610',   # raphe nuclei neuron
            'CL:0002611',   # neuron of the dorsal spinal cord
            'CL:0002612',   # neuron of the ventral spinal cord
            'CL:0002613',   # striatum neuron
            'CL:0002614',   # neuron of the substantia nigra
        ],
    }

    if path is None:
        return default_ontologies.get(organism, default_ontologies['human'])
    else:
        try:
            with open(path, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            print(f"ERROR: Ontology terms file not found: {path}")
            sys.exit(1)

        ONTOLOGY_TERMS = []
        invalid = []   # (line_number, raw_term) for every line that isn't a valid term
        for line_no, raw in enumerate(lines, start=1):
            term = raw.strip()
            if not term:
                continue  # skip blank lines
            if re.match(r'^(UBERON|CL):\d+$', term):
                ONTOLOGY_TERMS.append(term)
            else:
                invalid.append((line_no, term))

        if invalid:
            print(f"ERROR: Invalid ontology terms in {path}:")
            for line_no, term in invalid:
                print(f"  line {line_no}: '{term}'")
            print("Each term must look like 'UBERON:0000000' or 'CL:0000000'.")
            sys.exit(1)

        if not ONTOLOGY_TERMS:
            print(f"ERROR: Ontology terms file is empty: {path}")
            print("Please provide a file with valid ontology terms or leave the argument out to use the default set.")
            sys.exit(1)

        return ONTOLOGY_TERMS

ENSEMBL_SPECIES = {'human': 'human', 'mouse': 'mus_musculus'}

def fetch_ref_seq(chrom, start, end, organism='human'):
    import requests
    """Fetch reference sequence from Ensembl REST API (1-based, inclusive coords)."""
    species = ENSEMBL_SPECIES.get(organism, 'human')
    chrom_num = chrom.replace('chr', '')
    url = f"https://rest.ensembl.org/sequence/region/{species}/{chrom_num}:{start}..{end}:1"
    r = requests.get(url, headers={"Content-Type": "text/plain"}, timeout=15)
    r.raise_for_status()
    return r.text.strip()

def build_variant(row, organism='human'):
    """Build a genome.Variant from a row of variants_df.

    - SNPs / anchored indels: REF and ALT come straight from the input file.
    - Deletions (empty REF and empty ALT): REF = fetched from Ensembl using deletion_length, ALT = ''.
    - Insertions (empty REF, ALT given): REF = '', ALT = the inserted bases.
    """
    if not row.REF and not row.ALT:
        # deletion specified only by length -> fetch the deleted reference bases
        if not row.get('deletion_length', ''):
            print(f"ERROR: Deletion variant with ID {row.ID} is missing deletion_length column, which is required to fetch the reference sequence from Ensembl.")
            print("Please add a deletion_length column to your input file with the length of each deletion variant, or provide the reference bases directly in the REF column.")
            sys.exit(1)
        print(f'  Fetching reference sequence from Ensembl...')
        ref = fetch_ref_seq(row.CHROM, row.POS, row.POS + int(row.get('deletion_length', '')) - 1, organism)
        alt = ''
    elif not row.REF:
        # insertion: no reference bases, ALT holds the inserted sequence
        ref = ''
        alt = row.ALT
    else:
        ref = row.REF
        alt = row.ALT

    return genome.Variant(
        chromosome=str(row.CHROM),
        position=int(row.POS),
        reference_bases=ref,
        alternate_bases=alt,
        name=row.ID,
    )

def resolve_gene_symbol(gtf, gene_symbol):
    if gene_symbol in gtf['gene_name'].values:
        return gene_symbol
    matches = gtf[gtf['gene_name'].str.upper() == gene_symbol.upper()]['gene_name'].unique()
    if len(matches) == 1:
        return matches[0]
    elif len(matches) == 0:
        raise ValueError(f"Gene '{gene_symbol}' not found in gene annotations.")
    else:
        raise ValueError(f"Ambiguous gene symbol '{gene_symbol}': matches {list(matches)}")

def gtf_to_alphagenome_schema(gtf_df):
    """Convert a gtfparse-parsed GTF into the schema AlphaGenome expects.

    gtfparse returns the 8 standard GTF fields in lowercase ('seqname', 'feature',
    'start', ...) with 1-based coordinates. AlphaGenome's hosted feather files use
    capitalised names ('Chromosome', 'Feature', 'Start', ...) with 0-based 'Start'.
    """
    rename = {
        'seqname': 'Chromosome', 'seqid': 'Chromosome', 'chromosome': 'Chromosome',
        'source': 'Source',
        'feature': 'Feature',
        'start': 'Start',
        'end': 'End',
        'score': 'Score',
        'strand': 'Strand',
        'frame': 'Frame',
    }
    gtf_df = gtf_df.rename(columns={c: rename[c.lower()] for c in gtf_df.columns if c.lower() in rename})

    missing = {'Chromosome', 'Feature', 'Start', 'End', 'Strand'} - set(gtf_df.columns)
    if missing:
        print(f"ERROR: Converted GTF is missing required columns: {', '.join(sorted(missing))}")
        print("The file does not look like a standard GENCODE/Ensembl GTF.")
        sys.exit(1)

    # GTF is 1-based inclusive; AlphaGenome's feather stores 0-based Start.
    gtf_df['Start'] = gtf_df['Start'].astype('int64') - 1
    gtf_df['End'] = gtf_df['End'].astype('int64')
    return gtf_df

def setup_alphagenome(api_key, organism, gtf_path=None):
    global genome, dna_client, variant_scorers, plot_components, junction_data, gene_annotation

    print('Loading software libraries (this can take up to a minute the first time)...', flush=True)
    from alphagenome.data import gene_annotation, genome, junction_data
    from alphagenome.data import transcript as transcript_utils
    from alphagenome.models import dna_client, variant_scorers
    from alphagenome.visualization import plot_components
    import pandas as pd

    print(f'Loading AlphaGenome model...')
    dna_model = dna_client.create(api_key)

    GTF_URLS = {
        'human': 'https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather',
        'mouse': 'https://storage.googleapis.com/alphagenome/reference/gencode/mm10/gencode.vM23.annotation.gtf.gz.feather',
    }

    if gtf_path:
        if not gtf_path.endswith('.feather'):
            try:
                import gtfparse
            except ImportError:
                print('Installing gtfparse (one-time)...')
                subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'gtfparse'])
                import gtfparse
            feather_path = re.sub(r'\.gtf(\.gz)?$', '.feather', gtf_path)
            if feather_path == gtf_path:
                feather_path = gtf_path + '.feather'
            print(f'Reading GTF file {gtf_path} (this can take a few minutes for a full genome)...')
            gtf_df = gtfparse.read_gtf(gtf_path)
            if hasattr(gtf_df, 'to_pandas'):
                gtf_df = gtf_df.to_pandas()
            print('Converting to AlphaGenome format...')
            gtf_df = gtf_to_alphagenome_schema(gtf_df)
            print(f'Saving converted annotations to {feather_path} (so future runs can skip the conversion)...')
            gtf_df.to_feather(feather_path)
            print('Gene annotations ready.')
            gtf = gtf_df
        else:
            print(f'Loading gene annotations from {gtf_path}...')
            gtf = pd.read_feather(gtf_path)
    else:
        if organism not in GTF_URLS:
            print(f"ERROR: Unsupported organism: {organism}. Supported organisms are: human and mouse.")
            sys.exit(1)

        cache_path = os.path.join(os.path.expanduser('~'), '.alphagenome_cache', f'gtf_{organism}.feather')
        if os.path.exists(cache_path):
            print('Gene annotations loaded from cache.')
            gtf = pd.read_feather(cache_path)
        else:
            print('Downloading gene annotations (~100 MB, this happens only once and may take a minute)...')
            gtf = pd.read_feather(GTF_URLS[organism])
            os.makedirs(os.path.dirname(cache_path), exist_ok=True)
            gtf.to_feather(cache_path)
            print('Gene annotations downloaded and cached for future runs.')



    print('Preparing transcript annotations (this can take a minute on a full genome)...')
    gtf_protein_coding = gene_annotation.filter_protein_coding(gtf)
    if organism == 'human':
        gtf_mane = gene_annotation.filter_to_mane_select_transcript(gtf_protein_coding)
    else:
        gtf_mane = gene_annotation.filter_to_longest_transcript(gtf_protein_coding)

    transcript_extractor         = transcript_utils.TranscriptExtractor(gtf_mane)
    longest_transcript_extractor = transcript_utils.TranscriptExtractor(
        gene_annotation.filter_to_longest_transcript(gtf_mane)
    )

    return dna_model, gtf, gtf_mane, transcript_extractor, longest_transcript_extractor

def plot_output_type(output, variant, plot_interval, longest_transcripts, track_attr, title, has_biosample=True, save_dir=None):
    ref = getattr(output.reference, track_attr).filter_to_positive_strand()
    alt = getattr(output.alternate, track_attr).filter_to_positive_strand()

    transcript_panel = plot_components.TranscriptAnnotation(longest_transcripts)
    variant_marker   = plot_components.VariantAnnotation([variant])
    ylabel = '{biosample_name} ({strand})\n{name}' if has_biosample else '{name} ({strand})'

    REF_ALT_COLORS    = {'REF': 'blue', 'ALT': 'red'}
    _ = plot_components.plot(
        [transcript_panel, plot_components.OverlaidTracks(tdata={'REF': ref, 'ALT': alt}, colors=REF_ALT_COLORS, ylabel_template=ylabel)],
        annotations=[variant_marker], interval=plot_interval,
        title=f'{title}: REF vs ALT\n{variant.name}',
    )
    if save_dir:
        plt.savefig(os.path.join(save_dir, f'{track_attr}_overlay_{variant.name}.png'), bbox_inches='tight', dpi=150)
        plt.close()

    _ = plot_components.plot(
        [transcript_panel, plot_components.Tracks(tdata=alt - ref, ylabel_template=ylabel, filled=True)],
        annotations=[variant_marker], interval=plot_interval,
        title=f'{title}: ALT - REF\n{variant.name}',
    )
    if save_dir:
        plt.savefig(os.path.join(save_dir, f'{track_attr}_diff_{variant.name}.png'), bbox_inches='tight', dpi=150)
        plt.close()


def main():
    import pandas as pd
    ensure_packages()
    parser = argparse.ArgumentParser(
        description="Run AlphaGenome analysis", 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("Optional arguments")

    required.add_argument("--api-key", required=True, help="Your AlphaGenome API key, which can be obtained from https://deepmind.google.com/science/alphagenome/")

    required.add_argument("--input-file", required=False, help="Path to the input file containing genomic data, file format should be .tsv or .csv")

    optional.add_argument("--ontology-terms", required=False, help="Path to the file with ontology terms")

    optional.add_argument("--output-dir", required=False, help="Directory where the results will be saved, defaults to results_alphagenome which will be created if it doesn't exist in the current directory", default="./results_alphagenome")

    optional.add_argument("--organism", required=False, help="The organism for which the analysis will be run, defaults to human. At the moment, only human and mouse are supported", default="human")

    optional.add_argument("--gtf", required=False, help="Path to the GTF file containing gene annotation data")

    optional.add_argument("--list-ontologies", action="store_true", help="List all available ontology terms for the given organism (only human and mouse available) and write them to a TSV file, then exit")

    args = parser.parse_args()

    if args.list_ontologies:
        from alphagenome.models import dna_client as _dna_client
        import pandas as pd
        print(f'Loading AlphaGenome model...')
        _model = _dna_client.create(args.api_key)
        organism_enum = _dna_client.Organism.MUS_MUSCULUS if args.organism == 'mouse' else _dna_client.Organism.HOMO_SAPIENS
        print(f'Fetching ontology terms for {args.organism}...')
        meta = _model.output_metadata(organism_enum).concatenate()
        unique = (
            meta[['ontology_curie', 'biosample_name', 'biosample_type']]
            .drop_duplicates('ontology_curie')
            .sort_values('ontology_curie')
            .reset_index(drop=True)
        )
        out_path = f'{args.output_dir}/available_ontologies_{args.organism}.tsv'
        unique.to_csv(out_path, sep='\t', index=False)
        print(f'Wrote {len(unique)} ontology terms to {out_path}')
        sys.exit(0)


    if not args.input_file:
        parser.error("--input-file is required")

    variants_df = load_variants(args.input_file)

    PREDICTION_LENGTH = 2**20
    ONTOLOGY_TERMS = load_ontology_terms(args.ontology_terms, args.organism)

    dna_model, gtf, gtf_mane, transcript_extractor, longest_transcript_extractor = setup_alphagenome(args.api_key, args.organism, args.gtf)

    REQUESTED_OUTPUTS = {
    dna_client.OutputType.RNA_SEQ,
    dna_client.OutputType.CAGE,
    dna_client.OutputType.ATAC,
    dna_client.OutputType.DNASE,
    dna_client.OutputType.CHIP_HISTONE,
    dna_client.OutputType.CHIP_TF,
    dna_client.OutputType.PROCAP,
    dna_client.OutputType.SPLICE_SITES,
    dna_client.OutputType.SPLICE_SITE_USAGE,
    dna_client.OutputType.SPLICE_JUNCTIONS,
    }



    all_scores           = []   # accumulates scores_df from every variant
    outputs              = {}   # stores predict_variant output keyed by variant ID (for plotting later)
    variants             = {}   # stores genome.Variant objects keyed by ID
    prediction_intervals = {}   # stores prediction interval keyed by variant ID (for plot clamping)
    gene_strands         = {}   # stores authoritative gene strand (from GTF) keyed by variant ID

    ORGANISM = dna_client.Organism.MUS_MUSCULUS if args.organism == 'mouse' else dna_client.Organism.HOMO_SAPIENS
    # Use only the scorers that support this organism (e.g. polyadenylation is human-only),
    # otherwise score_variant fails on a single unsupported scorer.
    SCORERS = variant_scorers.get_recommended_scorers(ORGANISM.value)

    for _, row in variants_df.iterrows():
        print(f'\n=== {row.ID} ({row.get("cDNA", "")}, {row.get("HGVS", "")}) ===')

        try:
            gene_symbol = resolve_gene_symbol(gtf, row.Gene)
            gene_interval = gene_annotation.get_gene_interval(gtf, gene_symbol=gene_symbol)
            prediction_interval = gene_interval.resize(PREDICTION_LENGTH)

            # The gene annotation is the source of truth for strand. If the user supplied a
            # Strandedness value that disagrees, warn; always use the annotated strand for
            # plotting so a wrong value can't silently corrupt the splice-junction tracks.
            true_strand = gene_interval.strand
            user_strand = row.get('Strandedness', '')
            if user_strand and user_strand != true_strand:
                print(f"  WARNING: Strandedness for {row.ID} says '{user_strand}', but gene "
                      f"'{gene_symbol}' is on the '{true_strand}' strand. Using '{true_strand}'.")
            gene_strands[row.ID] = true_strand

            v = build_variant(row, args.organism)
            variants[row.ID] = v
            print(f'  Variant: {v}')

            # Make sure the variant actually falls inside the gene's prediction window before
            # calling the API. The most common cause of a mismatch is using variant coordinates
            # from the wrong organism (e.g. human chr1 coordinates with --organism mouse).
            if v.chromosome != prediction_interval.chromosome:
                print(f"  WARNING: skipping {row.ID} — variant is on {v.chromosome} but gene "
                      f"'{gene_symbol}' is on {prediction_interval.chromosome}. "
                      f"Check that your coordinates match --organism {args.organism}.")
                continue
            if not (prediction_interval.start <= v.position <= prediction_interval.end):
                print(f"  WARNING: skipping {row.ID} — position {v.position} is outside gene "
                      f"'{gene_symbol}' ({prediction_interval.chromosome}:{prediction_interval.start}-{prediction_interval.end}). "
                      f"Check that POS and Gene refer to the same location and organism.")
                continue

            flags = [name for name, val in [
                ('snv',        v.is_snv),
                ('deletion',   v.is_deletion),
                ('insertion',  v.is_insertion),
                ('indel',      v.is_indel),
                ('frameshift', v.is_frameshift),
                ('structural', v.is_structural),
            ] if val]
            print(f'  [{", ".join(flags) if flags else "unknown"}]')

            # --- predict ---
            print('  Running predictions (this may take a minute)...')

            output = dna_model.predict_variant(
                interval=prediction_interval,
                variant=v,
                requested_outputs=REQUESTED_OUTPUTS,
                ontology_terms=ONTOLOGY_TERMS,
                organism=ORGANISM,
            )
            outputs[row.ID] = output
            prediction_intervals[row.ID] = prediction_interval
            print('  predict_variant: done')

            # --- score ---
            print('  Scoring variant...')
            scores = dna_model.score_variant(
                interval=prediction_interval,
                variant=v,
                variant_scorers=SCORERS,
                organism=ORGANISM,
            )

            # tidy_scores formats the raw AnnData output into one score per row and handles
            # variant-centric vs gene-centric scorers, gene/track strand matching, missing
            # quantile scores, and empty results (returns None) — all robustly.
            scores_df = variant_scorers.tidy_scores(scores)

            if scores_df is None or scores_df.empty:
                print(f'  WARNING: no scores were produced for {row.ID}; skipping score file (predictions and plots are still saved).')
            else:
                # Add the variant's input metadata so the score file is self-contained.
                scores_df['input_gene']   = row.Gene
                scores_df['CHROM']        = row.CHROM
                scores_df['POS']          = row.POS
                scores_df['cDNA']         = row.get('cDNA', 'N/A')
                scores_df['HGVS']         = row.get('HGVS', 'N/A')
                scores_df['Phenotype']    = row.get('Phenotype', 'N/A')
                scores_df['is_snv']       = v.is_snv
                scores_df['is_deletion']  = v.is_deletion
                scores_df['is_insertion'] = v.is_insertion
                scores_df['is_indel']     = v.is_indel
                scores_df['is_frameshift'] = v.is_frameshift
                scores_df['is_structural'] = v.is_structural

                # --- save ---
                safe_id = re.sub(r'[^\w\-.]', '_', str(v.name))
                out_dir = f'{args.output_dir}/{safe_id}'
                os.makedirs(out_dir, exist_ok=True)
                out_path = f'{out_dir}/scores_{safe_id}.tsv'
                scores_df.to_csv(out_path, sep='\t', index=False)
                print(f'  Wrote {out_path}  ({scores_df.shape[0]} rows)')

                all_scores.append(scores_df)

        except Exception as e:
            err = str(e).lower()
            if any(word in err for word in ['permission', 'unauthorized', 'unauthenticated', 'api key', '403', '401']):
                print(f"ERROR: API authentication failed — check your API key.\n  {e}")
                sys.exit(1)
            print(f'  WARNING: {row.ID} failed — {e}')
            continue

    if all_scores:
        combined_df = pd.concat(all_scores, ignore_index=True)
        combined_path = f'{args.output_dir}/all_variants_scores.tsv'
        combined_df.to_csv(combined_path, sep='\t', index=False)
        print(f'Combined: {combined_df.shape}  -> {combined_path}')

    if not outputs:
        print("ERROR: No variants were successfully processed. Check your API key and input file.")
        sys.exit(1)

    OUTPUT_DIR = args.output_dir

    for variant_id, output in outputs.items():
        print(f'\n=== plotting {variant_id} ===')
        row     = variants_df[variants_df['ID'] == variant_id].iloc[0]
        variant = variants[variant_id]

        safe_id = re.sub(r'[^\w\-.]', '_', str(variant_id))
        plots_dir = os.path.join(OUTPUT_DIR, safe_id, 'plots')
        os.makedirs(plots_dir, exist_ok=True)

        WINDOW_BP = int(row.get('WINDOW_BP', '') or 24000)
        CENTRE_POS = int(row.get('CENTRE_POS', '') or variant.position)

        # use the annotated strand from the GTF (validated during prediction), not the user's input
        GENE_STRAND = gene_strands[variant_id]

        # build plot interval centered on variant, clamped to prediction interval
        pred_iv = prediction_intervals[variant_id]
        raw_start = CENTRE_POS - WINDOW_BP // 2
        raw_end   = CENTRE_POS + WINDOW_BP // 2
        if raw_start < pred_iv.start or raw_end > pred_iv.end or variant.chromosome != pred_iv.chromosome:
            print(f"  WARNING: plot interval chr{variant.chromosome}:{raw_start}-{raw_end} falls outside "
                  f"the prediction region. Check CENTRE_POS and WINDOW_BP values. "
                  f"Falling back to gene-centered window.")
            raw_start = pred_iv.start + (pred_iv.end - pred_iv.start) // 2 - WINDOW_BP // 2
            raw_end   = raw_start + WINDOW_BP
            raw_start = max(raw_start, pred_iv.start)
            raw_end   = min(raw_end,   pred_iv.end)
        v_plot_interval = genome.Interval(
            chromosome=variant.chromosome,
            start=raw_start,
            end=raw_end,
        )
        v_longest_transcripts = longest_transcript_extractor.extract(v_plot_interval)
        v_all_transcripts     = transcript_extractor.extract(v_plot_interval)

        # per-track plots
        plot_output_type(output, variant, v_plot_interval, v_longest_transcripts, 'rna_seq',           'RNA-seq',          save_dir=plots_dir)
        plot_output_type(output, variant, v_plot_interval, v_longest_transcripts, 'cage',              'CAGE',             save_dir=plots_dir)
        plot_output_type(output, variant, v_plot_interval, v_longest_transcripts, 'atac',              'ATAC',             save_dir=plots_dir)
        plot_output_type(output, variant, v_plot_interval, v_longest_transcripts, 'dnase',             'DNase',            save_dir=plots_dir)
        plot_output_type(output, variant, v_plot_interval, v_longest_transcripts, 'chip_histone',      'ChIP-histone',     save_dir=plots_dir)
        plot_output_type(output, variant, v_plot_interval, v_longest_transcripts, 'chip_tf',           'ChIP-TF',          save_dir=plots_dir)
        plot_output_type(output, variant, v_plot_interval, v_longest_transcripts, 'splice_sites',      'Splice sites',     has_biosample=False, save_dir=plots_dir)
        plot_output_type(output, variant, v_plot_interval, v_longest_transcripts, 'splice_site_usage', 'Splice site usage', save_dir=plots_dir)

        # splice junctions overlay
        ref_junc = output.reference.splice_junctions.filter_to_strand(GENE_STRAND)
        alt_junc = output.alternate.splice_junctions.filter_to_strand(GENE_STRAND)
        _ = plot_components.plot(
            [
                plot_components.TranscriptAnnotation(v_all_transcripts),
                plot_components.Sashimi(ref_junc, ylabel_template='REF {biosample_name} ({strand})\n{name}'),
                plot_components.Sashimi(alt_junc, ylabel_template='ALT {biosample_name} ({strand})\n{name}'),
            ],
            interval=v_plot_interval,
            annotations=[plot_components.VariantAnnotation([variant])],
            title=f'Splice junctions: REF vs ALT\n{variant_id}',
        )
        plt.savefig(os.path.join(plots_dir, f'splice_junctions_overlay_{variant_id}.png'), bbox_inches='tight', dpi=150)
        plt.close()

        # splice junctions diff
        alt_full = output.alternate.splice_junctions
        ref_full = output.reference.splice_junctions
        diff = junction_data.JunctionData(
            junctions=alt_full.junctions,
            values=alt_full.values - ref_full.values,
            metadata=alt_full.metadata,
            interval=alt_full.interval,
        )
        _ = plot_components.plot(
            [
                plot_components.TranscriptAnnotation(v_all_transcripts),
                plot_components.Sashimi(
                    junction_track=diff.filter_to_strand(GENE_STRAND),
                    ylabel_template='{biosample_name} ({strand})\n{name}',
                    filter_threshold=0.01,
                ),
            ],
            annotations=[plot_components.VariantAnnotation([variant])],
            interval=v_plot_interval,
            title=f'Splice junctions: ALT - REF\n{variant_id}',
        )
        plt.savefig(os.path.join(plots_dir, f'splice_junctions_diff_{variant_id}.png'), bbox_inches='tight', dpi=150)
        plt.close()

        print(f'  saved {len(os.listdir(plots_dir))} plots -> {plots_dir}')

    print(f'\nAll done! Results saved to {OUTPUT_DIR}')


if __name__ == "__main__":
    main()