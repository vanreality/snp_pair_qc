#!/usr/bin/env python3
"""
SNP VAF Correlation Heatmap Plotter

This script generates a correlation heatmap comparing variant allele frequencies (VAF)
between multiple pileup samples. It reads a samplesheet containing cfDNA and tissue
pileup files, calculates pairwise Pearson correlations between all samples, and
visualizes the results as a heatmap.

The tool performs the following steps:
1. Reads pileup files and calculates VAF for each sample
2. Filters SNP sites by minimum depth threshold
3. For each pair of pileup files, computes Pearson correlation on SNPs common to that pair only
4. Builds an N x N correlation matrix and generates a heatmap visualization
"""

import sys
import gzip
from pathlib import Path
from typing import Dict, Tuple, List
from itertools import combinations

import click
import polars as pl
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn
from rich.table import Table
from rich.panel import Panel

# Initialize rich console for output formatting
console = Console()

# Set matplotlib backend to avoid display issues
plt.switch_backend('Agg')


def read_samplesheet(samplesheet_file: Path) -> pl.DataFrame:
    """
    Read and validate the input samplesheet CSV file.

    The samplesheet must contain three columns:
    - sample: Sample identifier (used for labels)
    - cfDNA_pileup: Path to cfDNA pileup file
    - pileup: Path to tissue pileup file

    Args:
        samplesheet_file (Path): Path to the samplesheet CSV file

    Returns:
        pl.DataFrame: Validated samplesheet with required columns

    Raises:
        FileNotFoundError: If samplesheet file doesn't exist
        ValueError: If required columns are missing or file is empty
    """
    if not samplesheet_file.exists():
        raise FileNotFoundError(f"Samplesheet file not found: {samplesheet_file}")

    console.print(f"[cyan]●[/cyan] Reading samplesheet: {samplesheet_file}")

    try:
        df = pl.read_csv(samplesheet_file)

        # Validate required columns
        required_columns = ['sample', 'cfDNA_pileup', 'pileup']
        missing_cols = set(required_columns) - set(df.columns)
        if missing_cols:
            raise ValueError(
                f"Samplesheet missing required columns: {missing_cols}. "
                f"Required columns: {required_columns}"
            )

        if df.is_empty():
            raise ValueError("Samplesheet is empty")

        # Check for duplicate sample names
        dup_mask = df.filter(pl.col('sample').is_duplicated())
        if not dup_mask.is_empty():
            duplicate_samples = dup_mask['sample'].unique().to_list()
            raise ValueError(f"Duplicate sample names found: {duplicate_samples}")

        console.print(f"[green]✓[/green] Loaded {len(df)} samples from samplesheet")
        return df

    except pl.exceptions.NoDataError:
        raise ValueError("Samplesheet file is empty")
    except Exception as e:
        console.print(f"[red]Error reading samplesheet: {e}[/red]")
        raise


def read_pileup_file(pileup_file: Path, sample_name: str,
                     min_depth: int = 10) -> pl.DataFrame:
    """
    Read and parse a pileup file into a polars DataFrame with VAF calculation.

    Handles both gzipped (.gz) and uncompressed TSV files. The file is expected to
    have columns: chr, pos, ref, alt, af, cfDNA_ref_reads, cfDNA_alt_reads, current_depth.
    Filters SNP sites by minimum depth and calculates variant allele frequency (VAF).

    Args:
        pileup_file (Path): Path to the pileup file
        sample_name (str): Sample identifier for display and tracking
        min_depth (int): Minimum sequencing depth threshold (default: 10)

    Returns:
        pl.DataFrame: DataFrame with filtered pileup data and calculated VAF

    Raises:
        FileNotFoundError: If pileup file doesn't exist
        ValueError: If file format is invalid or required columns are missing
    """
    if not pileup_file.exists():
        raise FileNotFoundError(f"Pileup file not found for {sample_name}: {pileup_file}")

    try:
        if pileup_file.suffix == '.gz':
            open_func = gzip.open
            mode = 'rt'
        else:
            open_func = open
            mode = 'r'

        expected_columns = ['chr', 'pos', 'ref', 'alt', 'af',
                           'cfDNA_ref_reads', 'cfDNA_alt_reads', 'current_depth']

        with open_func(pileup_file, mode) as f:
            df = pl.read_csv(
                f,
                separator='\t',
                schema_overrides={
                    'chr': pl.Utf8,
                    'pos': pl.Int64,
                    'ref': pl.Utf8,
                    'alt': pl.Utf8,
                    'af': pl.Float64,
                    'cfDNA_ref_reads': pl.Int64,
                    'cfDNA_alt_reads': pl.Int64,
                    'current_depth': pl.Int64,
                }
            )

        missing_cols = set(expected_columns) - set(df.columns)
        if missing_cols:
            raise ValueError(
                f"Pileup file for {sample_name} missing required columns: {missing_cols}"
            )

        if df.is_empty():
            raise ValueError(f"Pileup file for {sample_name} is empty")

        initial_count = len(df)
        df = df.filter(pl.col('current_depth') >= min_depth)

        if df.is_empty():
            raise ValueError(
                f"All {initial_count:,} SNP sites for {sample_name} filtered out by "
                f"min_depth={min_depth} threshold"
            )

        df = df.with_columns(
            (pl.when(pl.col('current_depth') > 0)
             .then(pl.col('cfDNA_alt_reads') / pl.col('current_depth'))
             .otherwise(0.0)
             .alias('vaf'))
        )
        df = df.with_columns(
            (pl.col('chr').cast(pl.Utf8) + ':' +
             pl.col('pos').cast(pl.Utf8) + ':' +
             pl.col('ref') + ':' + pl.col('alt') + ':' +
             pl.col('af').cast(pl.Utf8))
            .alias('snp_id')
        )
        return df

    except pl.exceptions.NoDataError:
        console.print(f"[red]Error: Pileup file for {sample_name} is empty[/red]")
        raise ValueError(f"Pileup file for {sample_name} is empty")
    except Exception as e:
        console.print(f"[red]Error reading pileup file for {sample_name}: {e}[/red]")
        raise


def load_all_pileups(samplesheet: pl.DataFrame, min_depth: int,
                     progress: Progress) -> Dict[str, Dict[str, pl.DataFrame]]:
    """
    Load all pileup files from the samplesheet and calculate VAF for each.

    Reads both cfDNA and tissue pileup files for all samples, applying depth
    filtering and VAF calculation using the same strategy across all samples.

    Args:
        samplesheet (pl.DataFrame): Samplesheet containing sample info and file paths
        min_depth (int): Minimum depth threshold for filtering
        progress (Progress): Rich progress bar instance for tracking

    Returns:
        Dict[str, Dict[str, pl.DataFrame]]: Nested dict:
            {sample_name: {'cfDNA': df_cfDNA, 'tissue': df_tissue}}
    """
    pileup_data = {}
    total_files = len(samplesheet) * 2
    task = progress.add_task("[cyan]Loading pileup files...", total=total_files)

    for row in samplesheet.iter_rows(named=True):
        sample = row['sample']
        pileup_data[sample] = {}
        progress.update(task, description=f"[cyan]Loading {sample} cfDNA pileup...")
        cfDNA_file = Path(row['cfDNA_pileup'])
        pileup_data[sample]['cfDNA'] = read_pileup_file(
            cfDNA_file, f"{sample}_cfDNA", min_depth
        )
        progress.advance(task)
        progress.update(task, description=f"[cyan]Loading {sample} tissue pileup...")
        tissue_file = Path(row['pileup'])
        pileup_data[sample]['tissue'] = read_pileup_file(
            tissue_file, f"{sample}_tissue", min_depth
        )
        progress.advance(task)

    progress.update(task, description="[green]✓ All pileup files loaded")
    return pileup_data


def _pileup_id_to_vaf_df(
    pileup_data: Dict[str, Dict[str, pl.DataFrame]], samples: List[str]
) -> Dict[str, pl.DataFrame]:
    """Build map pileup_id -> DataFrame with columns (snp_id, vaf)."""
    out = {}
    for sample in samples:
        out[f"{sample}_cfDNA"] = pileup_data[sample]['cfDNA'].select(['snp_id', 'vaf'])
        out[f"{sample}_tissue"] = pileup_data[sample]['tissue'].select(['snp_id', 'vaf'])
    return out


def calculate_correlation_matrix(
    pileup_data: Dict[str, Dict[str, pl.DataFrame]],
    samples: List[str],
    progress: Progress,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Calculate pairwise Pearson correlations between all pileup files.

    For each pair of pileup files, inner-joins on snp_id to get SNPs common to
    that pair only, then computes Pearson correlation on the two VAF vectors.
    Returns N x N correlation and p-value matrices plus list of pileup_ids.

    Args:
        pileup_data: Nested dict {sample: {'cfDNA': df, 'tissue': df}}
        samples: List of sample names
        progress: Rich progress bar instance

    Returns:
        Tuple of (corr_matrix 2D array, pval_matrix 2D array, pileup_ids list)
    """
    console.print("[yellow]●[/yellow] Calculating pairwise correlations (per-pair common SNPs)...")

    pileup_ids = []
    for sample in samples:
        pileup_ids.append(f"{sample}_cfDNA")
        pileup_ids.append(f"{sample}_tissue")
    id_to_df = _pileup_id_to_vaf_df(pileup_data, samples)
    n = len(pileup_ids)

    corr_matrix = np.zeros((n, n))
    pval_matrix = np.zeros((n, n))
    for i in range(n):
        corr_matrix[i, i] = 1.0
        pval_matrix[i, i] = 0.0

    total_pairs = (n * (n - 1)) // 2
    task = progress.add_task("[cyan]Computing correlations...", total=total_pairs)

    for i in range(n):
        for j in range(i + 1, n):
            id1, id2 = pileup_ids[i], pileup_ids[j]
            left = id_to_df[id1].rename({'vaf': 'vaf1'})
            right = id_to_df[id2].rename({'vaf': 'vaf2'})
            joined = left.join(right, on='snp_id', how='inner')
            if joined.is_empty() or len(joined) < 2:
                corr, pval = 0.0, 1.0
            else:
                vaf1 = joined['vaf1'].to_numpy()
                vaf2 = joined['vaf2'].to_numpy()
                if vaf1.std() == 0 or vaf2.std() == 0:
                    console.print(
                        f"[yellow]Warning: Zero variance in {id1} or {id2} "
                        f"(n_common={len(joined)}), setting correlation to 0[/yellow]"
                    )
                    corr, pval = 0.0, 1.0
                else:
                    try:
                        corr, pval = pearsonr(vaf1, vaf2)
                    except Exception as e:
                        console.print(
                            f"[red]Error calculating correlation between {id1} and {id2}: {e}[/red]"
                        )
                        raise
            corr_matrix[i, j] = corr_matrix[j, i] = corr
            pval_matrix[i, j] = pval_matrix[j, i] = pval
            progress.advance(task)

    progress.update(task, description="[green]✓ Correlation matrix computed")
    return corr_matrix, pval_matrix, pileup_ids


def plot_heatmap(
    corr_matrix: np.ndarray,
    pileup_ids: List[str],
    samples: List[str],
    output_prefix: str,
    title: str = 'SNP VAF Correlation Heatmap\n(Tissue vs cfDNA)',
    xtitle: str = 'cfDNA Samples',
    ytitle: str = 'Tissue Samples',
) -> Path:
    """
    Generate and save correlation heatmap as PDF file.

    Creates a heatmap with tissue samples on y-axis and cfDNA samples on x-axis,
    with configurable title and axis labels.

    Args:
        corr_matrix: N x N correlation matrix (numpy)
        pileup_ids: List of pileup identifiers in matrix order
        samples: List of sample names for tick labels
        output_prefix: Output file prefix for saving PDF
        title: Plot title
        xtitle: X-axis label
        ytitle: Y-axis label

    Returns:
        Path to the saved PDF file.
    """
    console.print("[yellow]●[/yellow] Generating correlation heatmap...")

    try:
        cfDNA_cols = [f"{s}_cfDNA" for s in samples]
        tissue_rows = [f"{s}_tissue" for s in samples]
        idx_row = [pileup_ids.index(r) for r in tissue_rows]
        idx_col = [pileup_ids.index(c) for c in cfDNA_cols]
        heatmap_data = corr_matrix[np.ix_(idx_row, idx_col)]

        fig_height = max(8, len(samples) * 0.5)
        fig_width = max(10, len(samples) * 0.5 + 2)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        sns.heatmap(
            heatmap_data,
            annot=True,
            fmt='.3f',
            cmap='RdYlGn',
            center=0.7,
            vmin=0.0,
            vmax=1.0,
            cbar_kws={'label': 'Pearson Correlation'},
            linewidths=0.5,
            linecolor='gray',
            square=True,
            ax=ax
        )

        ax.set_xlabel(xtitle, fontsize=12, fontweight='bold')
        ax.set_ylabel(ytitle, fontsize=12, fontweight='bold')
        ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
        ax.set_xticklabels(samples, rotation=45, ha='right')
        ax.set_yticklabels(samples, rotation=0)
        plt.tight_layout()

        output_file = Path(f"{output_prefix}.pdf")
        plt.savefig(output_file, format='pdf', dpi=300, bbox_inches='tight')
        plt.close()
        console.print(f"[green]✓[/green] Heatmap saved to: {output_file}")
        return output_file

    except Exception as e:
        console.print(f"[red]Error generating heatmap: {e}[/red]")
        raise


def display_summary_statistics(
    corr_matrix: np.ndarray, pileup_ids: List[str], samples: List[str]
) -> None:
    """
    Display summary statistics of correlation values.

    Calculates and displays statistics for tissue-cfDNA correlation pairs:
    mean, median, min, max, std, and counts for high/QC-pass ranges.

    Args:
        corr_matrix: N x N correlation matrix (numpy)
        pileup_ids: List of pileup identifiers in matrix order
        samples: List of sample names
    """
    console.print("\n[bold yellow]Correlation Summary Statistics[/bold yellow]")

    tissue_cfDNA_corrs = []
    for sample in samples:
        tissue_idx = pileup_ids.index(f"{sample}_tissue")
        cfDNA_idx = pileup_ids.index(f"{sample}_cfDNA")
        tissue_cfDNA_corrs.append(corr_matrix[tissue_idx, cfDNA_idx])
    corr_array = np.array(tissue_cfDNA_corrs)
    
    stats_table = Table(
        title="Tissue-cfDNA Correlation Statistics",
        show_header=True,
        header_style="bold cyan"
    )
    stats_table.add_column("Statistic", style="cyan")
    stats_table.add_column("Value", style="white", justify="right")
    
    stats_table.add_row("Number of Pairs", f"{len(corr_array)}")
    stats_table.add_row("Mean", f"{corr_array.mean():.4f}")
    stats_table.add_row("Median", f"{np.median(corr_array):.4f}")
    stats_table.add_row("Std Dev", f"{corr_array.std():.4f}")
    stats_table.add_row("Min", f"{corr_array.min():.4f}")
    stats_table.add_row("Max", f"{corr_array.max():.4f}")
    stats_table.add_row(
        "High Corr (>= 0.8)",
        f"{(corr_array >= 0.8).sum()} ({(corr_array >= 0.8).sum() / len(corr_array) * 100:.1f}%)"
    )
    stats_table.add_row(
        "QC Pass (0.6-0.8)",
        f"{((corr_array >= 0.6) & (corr_array <= 0.8)).sum()} "
        f"({((corr_array >= 0.6) & (corr_array <= 0.8)).sum() / len(corr_array) * 100:.1f}%)"
    )
    
    console.print(stats_table)
    console.print()


@click.command()
@click.option(
    '--input-samplesheet',
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help='Path to input samplesheet CSV with columns: sample, cfDNA_pileup, pileup'
)
@click.option(
    '--output-prefix',
    required=True,
    type=str,
    help='Output file prefix (will create {prefix}.pdf)'
)
@click.option(
    '--min-depth',
    default=10,
    type=int,
    help='Minimum sequencing depth threshold for filtering (default: 10)'
)
@click.option(
    '--title',
    default='SNP VAF Correlation Heatmap\n(Tissue vs cfDNA)',
    type=str,
    help='Plot title (default: SNP VAF Correlation Heatmap with Tissue vs cfDNA)'
)
@click.option(
    '--xtitle',
    default='cfDNA Samples',
    type=str,
    help='X-axis label (default: cfDNA Samples)'
)
@click.option(
    '--ytitle',
    default='Tissue Samples',
    type=str,
    help='Y-axis label (default: Tissue Samples)'
)
def main(
    input_samplesheet: Path,
    output_prefix: str,
    min_depth: int,
    title: str,
    xtitle: str,
    ytitle: str,
) -> None:
    """
    Generate SNP VAF correlation heatmap from multiple pileup samples.
    
    This tool reads a samplesheet containing cfDNA and tissue pileup files for
    multiple samples, calculates variant allele frequencies (VAF), computes pairwise
    Pearson correlations, and generates a heatmap visualization showing the correlation
    matrix with tissue samples on the y-axis and cfDNA samples on the x-axis.
    
    Input Samplesheet Format:
    
        The CSV file must contain three columns:
        - sample: Sample identifier (used for heatmap labels)
        - cfDNA_pileup: Path to cfDNA pileup file (TSV or TSV.GZ)
        - pileup: Path to tissue pileup file (TSV or TSV.GZ)
    
    Output:
    
        - {prefix}.pdf: Correlation heatmap visualization
    
    Exit Codes:
    
        0: Success
        1: Error during processing
    
    Examples:
    
        # Basic usage
        python snp_correlation_heatmap.py \\
            --input-samplesheet samples.csv \\
            --output-prefix correlation_results
        
        # With custom depth threshold
        python snp_correlation_heatmap.py \\
            --input-samplesheet samples.csv \\
            --output-prefix correlation_results \\
            --min-depth 20
    """
    console.print("\n[bold blue]SNP VAF Correlation Heatmap Plotter[/bold blue]")
    console.print("=" * 70 + "\n")
    
    # Display input parameters
    params_table = Table(
        title="Input Parameters",
        show_header=True,
        header_style="bold magenta"
    )
    params_table.add_column("Parameter", style="cyan", no_wrap=True)
    params_table.add_column("Value", style="white")
    
    params_table.add_row("Input Samplesheet", str(input_samplesheet))
    params_table.add_row("Output Prefix", output_prefix)
    params_table.add_row("Min Depth", str(min_depth))
    
    console.print(params_table)
    console.print()
    
    try:
        console.print("[bold yellow]Step 1: Reading samplesheet[/bold yellow]")
        samplesheet = read_samplesheet(input_samplesheet)
        samples = samplesheet['sample'].to_list()
        console.print()

        console.print("[bold yellow]Step 2: Loading pileup files[/bold yellow]")
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            console=console
        ) as progress:
            pileup_data = load_all_pileups(samplesheet, min_depth, progress)
        console.print()

        load_table = Table(
            title="Loaded Pileup Statistics",
            show_header=True,
            header_style="bold green"
        )
        load_table.add_column("Sample", style="cyan")
        load_table.add_column("cfDNA SNPs", style="white", justify="right")
        load_table.add_column("Tissue SNPs", style="white", justify="right")
        for sample in samples:
            cfDNA_count = len(pileup_data[sample]['cfDNA'])
            tissue_count = len(pileup_data[sample]['tissue'])
            load_table.add_row(sample, f"{cfDNA_count:,}", f"{tissue_count:,}")
        console.print(load_table)
        console.print()

        console.print("[bold yellow]Step 3: Calculating correlation matrix[/bold yellow]")
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            console=console
        ) as progress:
            corr_matrix, pval_matrix, pileup_ids = calculate_correlation_matrix(
                pileup_data, samples, progress
            )
        console.print()

        display_summary_statistics(corr_matrix, pileup_ids, samples)

        console.print("[bold yellow]Step 4: Generating heatmap[/bold yellow]")
        output_file = plot_heatmap(
            corr_matrix, pileup_ids, samples, output_prefix,
            title=title, xtitle=xtitle, ytitle=ytitle
        )
        console.print()

        console.print(Panel(
            f"[green]✓ Analysis completed successfully![/green]\n\n"
            f"Pair-wise correlations computed on per-pair common SNPs.\n"
            f"Output file: {output_file}",
            title="Success",
            style="green",
            expand=False
        ))
        console.print()
        sys.exit(0)
        
    except Exception as e:
        console.print(f"\n[bold red]✗ Error during analysis:[/bold red] {e}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()

