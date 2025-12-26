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
3. Merges all samples on common SNP sites
4. Computes pairwise Pearson correlations to create an N x N correlation matrix
5. Generates a heatmap visualization with sample labels
"""

import sys
import gzip
from pathlib import Path
from typing import Dict, Tuple, List
from datetime import datetime
from itertools import combinations

import click
import pandas as pd
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


def read_samplesheet(samplesheet_file: Path) -> pd.DataFrame:
    """
    Read and validate the input samplesheet CSV file.
    
    The samplesheet must contain three columns:
    - sample: Sample identifier (used for labels)
    - cfDNA_pileup: Path to cfDNA pileup file
    - pileup: Path to tissue pileup file
    
    Args:
        samplesheet_file (Path): Path to the samplesheet CSV file
        
    Returns:
        pd.DataFrame: Validated samplesheet with required columns
        
    Raises:
        FileNotFoundError: If samplesheet file doesn't exist
        ValueError: If required columns are missing or file is empty
    """
    if not samplesheet_file.exists():
        raise FileNotFoundError(f"Samplesheet file not found: {samplesheet_file}")
    
    console.print(f"[cyan]●[/cyan] Reading samplesheet: {samplesheet_file}")
    
    try:
        df = pd.read_csv(samplesheet_file)
        
        # Validate required columns
        required_columns = ['sample', 'cfDNA_pileup', 'pileup']
        missing_cols = set(required_columns) - set(df.columns)
        
        if missing_cols:
            raise ValueError(
                f"Samplesheet missing required columns: {missing_cols}. "
                f"Required columns: {required_columns}"
            )
        
        if df.empty:
            raise ValueError("Samplesheet is empty")
        
        # Check for duplicate sample names
        duplicates = df['sample'].duplicated()
        if duplicates.any():
            duplicate_samples = df.loc[duplicates, 'sample'].tolist()
            raise ValueError(f"Duplicate sample names found: {duplicate_samples}")
        
        console.print(f"[green]✓[/green] Loaded {len(df)} samples from samplesheet")
        
        return df
        
    except pd.errors.EmptyDataError:
        raise ValueError("Samplesheet file is empty")
    except Exception as e:
        console.print(f"[red]Error reading samplesheet: {e}[/red]")
        raise


def read_pileup_file(pileup_file: Path, sample_name: str, 
                     min_depth: int = 10) -> pd.DataFrame:
    """
    Read and parse a pileup file into a pandas DataFrame with VAF calculation.
    
    Handles both gzipped (.gz) and uncompressed TSV files. The file is expected to
    have columns: chr, pos, ref, alt, af, cfDNA_ref_reads, cfDNA_alt_reads, current_depth.
    Filters SNP sites by minimum depth and calculates variant allele frequency (VAF).
    
    Args:
        pileup_file (Path): Path to the pileup file
        sample_name (str): Sample identifier for display and tracking
        min_depth (int): Minimum sequencing depth threshold (default: 10)
        
    Returns:
        pd.DataFrame: DataFrame with filtered pileup data and calculated VAF
        
    Raises:
        FileNotFoundError: If pileup file doesn't exist
        ValueError: If file format is invalid or required columns are missing
        pd.errors.EmptyDataError: If file is empty or all sites filtered out
    """
    if not pileup_file.exists():
        raise FileNotFoundError(f"Pileup file not found for {sample_name}: {pileup_file}")
    
    try:
        # Determine if file is gzipped based on extension
        if pileup_file.suffix == '.gz':
            open_func = gzip.open
            mode = 'rt'
        else:
            open_func = open
            mode = 'r'
        
        # Read pileup file with specified columns
        expected_columns = ['chr', 'pos', 'ref', 'alt', 'af', 
                          'cfDNA_ref_reads', 'cfDNA_alt_reads', 'current_depth']
        
        with open_func(pileup_file, mode) as f:
            df = pd.read_csv(
                f,
                sep='\t',
                dtype={
                    'chr': str,
                    'pos': int,
                    'ref': str,
                    'alt': str,
                    'af': float,
                    'cfDNA_ref_reads': int,
                    'cfDNA_alt_reads': int,
                    'current_depth': int
                }
            )
        
        # Validate required columns are present
        missing_cols = set(expected_columns) - set(df.columns)
        if missing_cols:
            raise ValueError(
                f"Pileup file for {sample_name} missing required columns: {missing_cols}"
            )
        
        if df.empty:
            raise pd.errors.EmptyDataError(f"Pileup file for {sample_name} is empty")
        
        initial_count = len(df)
        
        # Filter by minimum depth threshold
        df = df[df['current_depth'] >= min_depth].copy()
        
        if df.empty:
            raise ValueError(
                f"All {initial_count:,} SNP sites for {sample_name} filtered out by "
                f"min_depth={min_depth} threshold"
            )
        
        # Calculate VAF (Variant Allele Frequency)
        # Handle division by zero: set VAF to 0 where depth is 0 (should not happen after filtering)
        df['vaf'] = np.where(
            df['current_depth'] > 0,
            df['cfDNA_alt_reads'] / df['current_depth'],
            0.0
        )
        
        # Create unique SNP identifier for merging
        df['snp_id'] = (
            df['chr'].astype(str) + ':' + 
            df['pos'].astype(str) + ':' + 
            df['ref'] + ':' + 
            df['alt'] + ':' + 
            df['af'].astype(str)
        )
        
        filtered_count = len(df)
        removed_count = initial_count - filtered_count
        
        return df
        
    except pd.errors.EmptyDataError as e:
        console.print(f"[red]Error: Pileup file for {sample_name} is empty[/red]")
        raise
    except Exception as e:
        console.print(f"[red]Error reading pileup file for {sample_name}: {e}[/red]")
        raise


def load_all_pileups(samplesheet: pd.DataFrame, min_depth: int,
                     progress: Progress) -> Dict[str, Dict[str, pd.DataFrame]]:
    """
    Load all pileup files from the samplesheet and calculate VAF for each.
    
    Reads both cfDNA and tissue pileup files for all samples, applying depth
    filtering and VAF calculation using the same strategy across all samples.
    
    Args:
        samplesheet (pd.DataFrame): Samplesheet containing sample info and file paths
        min_depth (int): Minimum depth threshold for filtering
        progress (Progress): Rich progress bar instance for tracking
        
    Returns:
        Dict[str, Dict[str, pd.DataFrame]]: Nested dictionary structure:
            {sample_name: {'cfDNA': df_cfDNA, 'tissue': df_tissue}}
        
    Raises:
        FileNotFoundError: If any pileup file is not found
        ValueError: If any pileup file fails validation or filtering
    """
    pileup_data = {}
    
    # Calculate total tasks for progress tracking
    total_files = len(samplesheet) * 2  # cfDNA + tissue for each sample
    task = progress.add_task(
        "[cyan]Loading pileup files...",
        total=total_files
    )
    
    for idx, row in samplesheet.iterrows():
        sample = row['sample']
        pileup_data[sample] = {}
        
        # Load cfDNA pileup
        progress.update(
            task, 
            description=f"[cyan]Loading {sample} cfDNA pileup..."
        )
        cfDNA_file = Path(row['cfDNA_pileup'])
        pileup_data[sample]['cfDNA'] = read_pileup_file(
            cfDNA_file, f"{sample}_cfDNA", min_depth
        )
        progress.advance(task)
        
        # Load tissue pileup
        progress.update(
            task,
            description=f"[cyan]Loading {sample} tissue pileup..."
        )
        tissue_file = Path(row['pileup'])
        pileup_data[sample]['tissue'] = read_pileup_file(
            tissue_file, f"{sample}_tissue", min_depth
        )
        progress.advance(task)
    
    progress.update(task, description="[green]✓ All pileup files loaded")
    
    return pileup_data


def merge_all_samples(pileup_data: Dict[str, Dict[str, pd.DataFrame]]) -> pd.DataFrame:
    """
    Merge all pileup data across samples to find common SNP sites.
    
    Performs successive inner joins on SNP ID (chr:pos:ref:alt:af) to identify
    SNP sites present in all samples. Each sample contributes both cfDNA and
    tissue VAF values to the final merged dataset.
    
    Args:
        pileup_data (Dict[str, Dict[str, pd.DataFrame]]): Nested dictionary with
            pileup data for all samples
            
    Returns:
        pd.DataFrame: Merged DataFrame with VAF columns for all samples at common SNPs
            Columns: snp_id, chr, pos, ref, alt, af, 
                    {sample1}_cfDNA_vaf, {sample1}_tissue_vaf, ...
        
    Raises:
        ValueError: If no common SNP sites found across all samples
    """
    console.print("[yellow]●[/yellow] Merging all samples on common SNP sites...")
    
    # Start with the first sample as base
    samples = list(pileup_data.keys())
    first_sample = samples[0]
    
    # Initialize merged dataframe with first sample's cfDNA data
    merged_df = pileup_data[first_sample]['cfDNA'][['snp_id', 'chr', 'pos', 'ref', 'alt', 'af', 'vaf']].copy()
    merged_df = merged_df.rename(columns={'vaf': f'{first_sample}_cfDNA_vaf'})
    
    # Merge first sample's tissue data
    tissue_df = pileup_data[first_sample]['tissue'][['snp_id', 'vaf']].copy()
    tissue_df = tissue_df.rename(columns={'vaf': f'{first_sample}_tissue_vaf'})
    merged_df = merged_df.merge(tissue_df, on='snp_id', how='inner')
    
    # Successively merge remaining samples
    for sample in samples[1:]:
        # Merge cfDNA data
        cfDNA_df = pileup_data[sample]['cfDNA'][['snp_id', 'vaf']].copy()
        cfDNA_df = cfDNA_df.rename(columns={'vaf': f'{sample}_cfDNA_vaf'})
        merged_df = merged_df.merge(cfDNA_df, on='snp_id', how='inner')
        
        # Merge tissue data
        tissue_df = pileup_data[sample]['tissue'][['snp_id', 'vaf']].copy()
        tissue_df = tissue_df.rename(columns={'vaf': f'{sample}_tissue_vaf'})
        merged_df = merged_df.merge(tissue_df, on='snp_id', how='inner')
        
        if merged_df.empty:
            raise ValueError(
                f"No common SNP sites after merging sample {sample}. "
                f"Consider lowering the min_depth threshold."
            )
    
    if merged_df.empty:
        raise ValueError(
            "No common SNP sites found across all samples. "
            "Consider lowering the min_depth threshold."
        )
    
    console.print(
        f"[green]✓[/green] Found {len(merged_df):,} common SNP sites across all samples"
    )
    
    return merged_df


def calculate_correlation_matrix(merged_df: pd.DataFrame, 
                                 samples: List[str],
                                 progress: Progress) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate pairwise Pearson correlations between all pileup files.
    
    Computes correlations for all pairwise combinations of pileup files
    (both cfDNA and tissue) across all samples, creating an N x N correlation
    matrix where N = number_of_samples * 2.
    
    Args:
        merged_df (pd.DataFrame): Merged DataFrame with VAF columns for all samples
        samples (List[str]): List of sample names
        progress (Progress): Rich progress bar instance for tracking
        
    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: 
            - Correlation matrix (N x N DataFrame with correlation coefficients)
            - P-value matrix (N x N DataFrame with p-values)
        
    Raises:
        ValueError: If correlation calculation fails for any pair
        RuntimeError: If matrix construction fails
    """
    console.print("[yellow]●[/yellow] Calculating pairwise correlations...")
    
    # Build list of all pileup identifiers (sample_type format)
    pileup_ids = []
    for sample in samples:
        pileup_ids.append(f"{sample}_cfDNA")
        pileup_ids.append(f"{sample}_tissue")
    
    n = len(pileup_ids)
    
    # Initialize correlation and p-value matrices
    corr_matrix = pd.DataFrame(
        np.zeros((n, n)),
        index=pileup_ids,
        columns=pileup_ids
    )
    pval_matrix = pd.DataFrame(
        np.zeros((n, n)),
        index=pileup_ids,
        columns=pileup_ids
    )
    
    # Set diagonal to 1.0 (perfect correlation with self)
    for i in range(n):
        corr_matrix.iloc[i, i] = 1.0
        pval_matrix.iloc[i, i] = 0.0
    
    # Calculate total number of pairs for progress tracking
    total_pairs = (n * (n - 1)) // 2
    task = progress.add_task(
        "[cyan]Computing correlations...",
        total=total_pairs
    )
    
    # Calculate pairwise correlations
    for i in range(n):
        for j in range(i + 1, n):
            id1 = pileup_ids[i]
            id2 = pileup_ids[j]
            
            # Extract VAF columns
            vaf1 = merged_df[f"{id1}_vaf"].values
            vaf2 = merged_df[f"{id2}_vaf"].values
            
            # Check for variance
            if vaf1.std() == 0 or vaf2.std() == 0:
                console.print(
                    f"[yellow]Warning: Zero variance in {id1} or {id2}, "
                    f"setting correlation to 0[/yellow]"
                )
                corr = 0.0
                pval = 1.0
            else:
                try:
                    corr, pval = pearsonr(vaf1, vaf2)
                except Exception as e:
                    console.print(
                        f"[red]Error calculating correlation between {id1} and {id2}: {e}[/red]"
                    )
                    raise
            
            # Fill both symmetric positions
            corr_matrix.loc[id1, id2] = corr
            corr_matrix.loc[id2, id1] = corr
            pval_matrix.loc[id1, id2] = pval
            pval_matrix.loc[id2, id1] = pval
            
            progress.advance(task)
    
    progress.update(task, description="[green]✓ Correlation matrix computed")
    
    return corr_matrix, pval_matrix


def plot_heatmap(corr_matrix: pd.DataFrame, samples: List[str],
                 output_prefix: str) -> Path:
    """
    Generate and save correlation heatmap as PDF file.
    
    Creates a clustered heatmap visualization with:
    - X-axis: cfDNA pileup samples
    - Y-axis: tissue pileup samples
    - Sample labels for both axes
    - Color scale representing correlation coefficients
    - Annotations showing correlation values
    
    Args:
        corr_matrix (pd.DataFrame): N x N correlation matrix
        samples (List[str]): List of sample names for labeling
        output_prefix (str): Output file prefix for saving PDF
        
    Returns:
        Path: Path to the saved PDF file
        
    Raises:
        IOError: If PDF file cannot be written
        ValueError: If plotting fails
    """
    console.print("[yellow]●[/yellow] Generating correlation heatmap...")
    
    try:
        # Reorder matrix to show cfDNA on x-axis and tissue on y-axis
        # X-axis: cfDNA samples (columns)
        cfDNA_cols = [f"{sample}_cfDNA" for sample in samples]
        # Y-axis: tissue samples (rows)
        tissue_rows = [f"{sample}_tissue" for sample in samples]
        
        # Extract submatrix for visualization
        heatmap_data = corr_matrix.loc[tissue_rows, cfDNA_cols]
        
        # Create figure with appropriate size
        fig_height = max(8, len(samples) * 0.5)
        fig_width = max(10, len(samples) * 0.5 + 2)
        
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        
        # Create heatmap with seaborn
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
        
        # Set labels
        ax.set_xlabel('cfDNA Samples', fontsize=12, fontweight='bold')
        ax.set_ylabel('Tissue Samples', fontsize=12, fontweight='bold')
        ax.set_title(
            'SNP VAF Correlation Heatmap\n(Tissue vs cfDNA)',
            fontsize=14,
            fontweight='bold',
            pad=20
        )
        
        # Format tick labels to show only sample names
        ax.set_xticklabels(samples, rotation=45, ha='right')
        ax.set_yticklabels(samples, rotation=0)
        
        # Adjust layout to prevent label cutoff
        plt.tight_layout()
        
        # Save figure
        output_file = Path(f"{output_prefix}.pdf")
        plt.savefig(output_file, format='pdf', dpi=300, bbox_inches='tight')
        plt.close()
        
        console.print(f"[green]✓[/green] Heatmap saved to: {output_file}")
        
        return output_file
        
    except Exception as e:
        console.print(f"[red]Error generating heatmap: {e}[/red]")
        raise


def display_summary_statistics(corr_matrix: pd.DataFrame, samples: List[str]) -> None:
    """
    Display summary statistics of correlation values.
    
    Calculates and displays statistics for tissue-cfDNA correlation pairs:
    - Mean, median, min, max correlation values
    - Standard deviation
    - Count of high correlations (>= 0.8)
    
    Args:
        corr_matrix (pd.DataFrame): N x N correlation matrix
        samples (List[str]): List of sample names
        
    Returns:
        None
    """
    console.print("\n[bold yellow]Correlation Summary Statistics[/bold yellow]")
    
    # Extract tissue-cfDNA correlations (diagonal in the tissue vs cfDNA submatrix)
    tissue_cfDNA_corrs = []
    for sample in samples:
        tissue_idx = f"{sample}_tissue"
        cfDNA_idx = f"{sample}_cfDNA"
        corr_value = corr_matrix.loc[tissue_idx, cfDNA_idx]
        tissue_cfDNA_corrs.append(corr_value)
    
    # Calculate statistics
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
def main(input_samplesheet: Path, output_prefix: str, min_depth: int) -> None:
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
        # Step 1: Read samplesheet
        console.print("[bold yellow]Step 1: Reading samplesheet[/bold yellow]")
        samplesheet = read_samplesheet(input_samplesheet)
        samples = samplesheet['sample'].tolist()
        console.print()
        
        # Step 2: Load all pileup files
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
        
        # Display loading statistics
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
            load_table.add_row(
                sample,
                f"{cfDNA_count:,}",
                f"{tissue_count:,}"
            )
        
        console.print(load_table)
        console.print()
        
        # Step 3: Merge all samples
        console.print("[bold yellow]Step 3: Merging samples on common SNPs[/bold yellow]")
        merged_df = merge_all_samples(pileup_data)
        console.print()
        
        # Step 4: Calculate correlation matrix
        console.print("[bold yellow]Step 4: Calculating correlation matrix[/bold yellow]")
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            console=console
        ) as progress:
            corr_matrix, pval_matrix = calculate_correlation_matrix(
                merged_df, samples, progress
            )
        console.print()
        
        # Step 5: Display summary statistics
        display_summary_statistics(corr_matrix, samples)
        
        # Step 6: Generate heatmap
        console.print("[bold yellow]Step 5: Generating heatmap[/bold yellow]")
        output_file = plot_heatmap(corr_matrix, samples, output_prefix)
        console.print()
        
        # Success message
        console.print(Panel(
            f"[green]✓ Analysis completed successfully![/green]\n\n"
            f"Common SNPs: {len(merged_df):,}\n"
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

