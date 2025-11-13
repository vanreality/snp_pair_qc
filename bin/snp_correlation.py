#!/usr/bin/env python3
"""
SNP VAF Correlation Calculator

This script calculates the Pearson correlation between tissue and cfDNA variant allele 
frequencies (VAF) at matching SNP sites. It performs quality control checks to determine 
if the correlation meets expected thresholds for sample matching.

The tool reads pileup data from tissue and cfDNA samples, filters by depth, merges on
common SNP sites, and computes correlation statistics to assess sample concordance.
"""

import sys
import gzip
from pathlib import Path
from typing import Tuple, Optional
from datetime import datetime

import click
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.table import Table
from rich.panel import Panel

# Initialize rich console for output formatting
console = Console()


def read_pileup_file(pileup_file: Path, progress: Progress, task_id: int,
                     sample_type: str) -> pd.DataFrame:
    """
    Read and parse a pileup file into a pandas DataFrame.
    
    Handles both gzipped (.gz) and uncompressed TSV files. The file is expected to
    have columns: chr, pos, ref, alt, af, cfDNA_ref_reads, cfDNA_alt_reads, current_depth.
    
    Args:
        pileup_file (Path): Path to the pileup file
        progress (Progress): Rich progress bar instance
        task_id (int): Task ID for progress tracking
        sample_type (str): Type of sample (e.g., "tissue" or "cfDNA") for display
        
    Returns:
        pd.DataFrame: DataFrame with pileup data including calculated VAF
        
    Raises:
        FileNotFoundError: If pileup file doesn't exist
        ValueError: If file format is invalid or required columns are missing
        pd.errors.EmptyDataError: If file is empty
    """
    if not pileup_file.exists():
        raise FileNotFoundError(f"{sample_type} pileup file not found: {pileup_file}")
    
    progress.update(task_id, description=f"Reading {sample_type} pileup file...")
    
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
                f"{sample_type} pileup file missing required columns: {missing_cols}"
            )
        
        if df.empty:
            raise pd.errors.EmptyDataError(f"{sample_type} pileup file is empty")
        
        # Calculate VAF (Variant Allele Frequency)
        # Handle division by zero: set VAF to 0 where depth is 0
        df['vaf'] = np.where(
            df['current_depth'] > 0,
            df['cfDNA_alt_reads'] / df['current_depth'],
            0.0
        )
        
        progress.update(task_id, advance=50)
        console.print(
            f"[green]✓[/green] Read {len(df):,} SNP sites from {sample_type} pileup"
        )
        
        return df
        
    except pd.errors.EmptyDataError as e:
        console.print(f"[red]Error: {sample_type} pileup file is empty[/red]")
        raise
    except Exception as e:
        console.print(f"[red]Error reading {sample_type} pileup file: {e}[/red]")
        raise


def filter_by_depth(df: pd.DataFrame, min_depth: int, sample_type: str) -> pd.DataFrame:
    """
    Filter pileup DataFrame by minimum depth threshold.
    
    Removes SNP sites with coverage below the specified minimum depth threshold.
    
    Args:
        df (pd.DataFrame): Input pileup DataFrame
        min_depth (int): Minimum depth threshold
        sample_type (str): Type of sample (e.g., "tissue" or "cfDNA") for display
        
    Returns:
        pd.DataFrame: Filtered DataFrame
        
    Raises:
        ValueError: If all sites are filtered out
    """
    initial_count = len(df)
    
    # Filter by minimum depth
    filtered_df = df[df['current_depth'] >= min_depth].copy()
    
    filtered_count = len(filtered_df)
    removed_count = initial_count - filtered_count
    
    if filtered_count == 0:
        raise ValueError(
            f"All {initial_count:,} {sample_type} SNP sites filtered out by "
            f"min_depth={min_depth} threshold"
        )
    
    console.print(
        f"[blue]●[/blue] {sample_type}: {filtered_count:,} sites remaining after "
        f"depth filtering (removed {removed_count:,} sites)"
    )
    
    return filtered_df


def merge_pileup_data(tissue_df: pd.DataFrame, cfDNA_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge tissue and cfDNA pileup data on common SNP sites.
    
    Performs an inner join on ['chr', 'pos', 'ref', 'alt', 'af'] to find SNP sites
    present in both samples. Suffixes '_tissue' and '_cfDNA' are added to distinguish
    columns from each sample.
    
    Args:
        tissue_df (pd.DataFrame): Filtered tissue pileup data
        cfDNA_df (pd.DataFrame): Filtered cfDNA pileup data
        
    Returns:
        pd.DataFrame: Merged DataFrame with matching SNP sites from both samples
        
    Raises:
        ValueError: If no common SNP sites are found
    """
    merge_keys = ['chr', 'pos', 'ref', 'alt', 'af']
    
    # Perform inner join to find common SNP sites
    merged_df = pd.merge(
        tissue_df,
        cfDNA_df,
        on=merge_keys,
        how='inner',
        suffixes=('_tissue', '_cfDNA')
    )
    
    if merged_df.empty:
        raise ValueError(
            "No common SNP sites found between tissue and cfDNA samples after merging"
        )
    
    console.print(
        f"[green]✓[/green] Merged data: {len(merged_df):,} common SNP sites"
    )
    
    return merged_df


def calculate_correlation(merged_df: pd.DataFrame) -> Tuple[float, float]:
    """
    Calculate Pearson correlation between tissue and cfDNA VAF values.
    
    Computes the Pearson correlation coefficient and p-value for the VAF values
    from tissue and cfDNA samples at matched SNP sites.
    
    Args:
        merged_df (pd.DataFrame): Merged DataFrame with VAF values from both samples
        
    Returns:
        Tuple[float, float]: Pearson correlation coefficient and p-value
        
    Raises:
        ValueError: If insufficient data points for correlation calculation
        RuntimeError: If correlation calculation fails
    """
    if len(merged_df) < 2:
        raise ValueError(
            f"Insufficient data points for correlation calculation: {len(merged_df)} sites"
        )
    
    try:
        # Extract VAF values for correlation
        tissue_vaf = merged_df['vaf_tissue'].values
        cfDNA_vaf = merged_df['vaf_cfDNA'].values
        
        # Check for variance in data
        if tissue_vaf.std() == 0 or cfDNA_vaf.std() == 0:
            raise ValueError("One or both VAF datasets have zero variance")
        
        # Calculate Pearson correlation
        correlation, p_value = pearsonr(tissue_vaf, cfDNA_vaf)
        
        console.print(
            f"[cyan]●[/cyan] Pearson correlation: {correlation:.4f} (p-value: {p_value:.2e})"
        )
        
        return correlation, p_value
        
    except Exception as e:
        console.print(f"[red]Error calculating correlation: {e}[/red]")
        raise RuntimeError(f"Correlation calculation failed: {e}")


def assess_qc_status(merged_snp_count: int, correlation: float,
                    min_snp_count: int) -> Tuple[str, str]:
    """
    Assess QC status based on merged SNP count and correlation value.
    
    QC passes if:
    - Merged SNP count >= min_snp_count AND
    - Correlation is in the range [0.6, 0.8]
    
    Args:
        merged_snp_count (int): Number of merged SNP sites
        correlation (float): Pearson correlation coefficient
        min_snp_count (int): Minimum required SNP count
        
    Returns:
        Tuple[str, str]: QC status ("PASS" or "FAIL") and reason message
    """
    reasons = []
    
    # Check SNP count threshold
    if merged_snp_count < min_snp_count:
        reasons.append(
            f"merged SNP count ({merged_snp_count:,}) < min_snp_count ({min_snp_count:,})"
        )
    
    # Check correlation range
    if not (0.6 <= correlation <= 0.8):
        reasons.append(
            f"correlation ({correlation:.4f}) not in range [0.6, 0.8]"
        )
    
    if reasons:
        status = "FAIL"
        reason = "; ".join(reasons)
    else:
        status = "PASS"
        reason = "All QC criteria met"
    
    return status, reason


def save_report(output_prefix: str, tissue_file: Path, cfDNA_file: Path,
                tissue_filtered_count: int, cfDNA_filtered_count: int,
                merged_snp_count: int, correlation: float, p_value: float,
                qc_status: str, qc_reason: str, min_depth: int,
                min_snp_count: int) -> Path:
    """
    Save correlation analysis report to a text file.
    
    Creates a detailed report including input parameters, filtering statistics,
    correlation results, and QC assessment.
    
    Args:
        output_prefix (str): Output file prefix
        tissue_file (Path): Path to tissue pileup file
        cfDNA_file (Path): Path to cfDNA pileup file
        tissue_filtered_count (int): Number of tissue SNPs after filtering
        cfDNA_filtered_count (int): Number of cfDNA SNPs after filtering
        merged_snp_count (int): Number of merged SNP sites
        correlation (float): Pearson correlation coefficient
        p_value (float): Correlation p-value
        qc_status (str): QC status ("PASS" or "FAIL")
        qc_reason (str): Reason for QC status
        min_depth (int): Minimum depth threshold used
        min_snp_count (int): Minimum SNP count threshold used
        
    Returns:
        Path: Path to the saved report file
        
    Raises:
        IOError: If report file cannot be written
    """
    report_file = Path(f"{output_prefix}_report.txt")
    
    try:
        with open(report_file, 'w') as f:
            # Write header
            f.write("=" * 80 + "\n")
            f.write("SNP VAF Correlation Analysis Report\n")
            f.write("=" * 80 + "\n\n")
            
            # Write timestamp
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            f.write(f"Analysis Date: {timestamp}\n\n")
            
            # Write input parameters
            f.write("Input Parameters:\n")
            f.write("-" * 80 + "\n")
            f.write(f"  Tissue Pileup:     {tissue_file}\n")
            f.write(f"  cfDNA Pileup:      {cfDNA_file}\n")
            f.write(f"  Min Depth:         {min_depth}\n")
            f.write(f"  Min SNP Count:     {min_snp_count:,}\n")
            f.write(f"  Output Prefix:     {output_prefix}\n\n")
            
            # Write filtering statistics
            f.write("Filtering Statistics:\n")
            f.write("-" * 80 + "\n")
            f.write(f"  Tissue SNPs (after depth >= {min_depth} filter):  "
                   f"{tissue_filtered_count:,}\n")
            f.write(f"  cfDNA SNPs (after depth >= {min_depth} filter):   "
                   f"{cfDNA_filtered_count:,}\n")
            f.write(f"  Merged SNP sites (common):                         "
                   f"{merged_snp_count:,}\n\n")
            
            # Write correlation results
            f.write("Correlation Analysis:\n")
            f.write("-" * 80 + "\n")
            f.write(f"  Pearson Correlation Coefficient:  {correlation:.6f}\n")
            f.write(f"  P-value:                          {p_value:.6e}\n\n")
            
            # Write QC assessment
            f.write("Quality Control Assessment:\n")
            f.write("-" * 80 + "\n")
            f.write(f"  QC Status:    {qc_status}\n")
            f.write(f"  QC Criteria:  Merged SNP count >= {min_snp_count:,} AND "
                   f"correlation in [0.6, 0.8]\n")
            f.write(f"  Reason:       {qc_reason}\n\n")
            
            # Write footer
            f.write("=" * 80 + "\n")
            f.write("End of Report\n")
            f.write("=" * 80 + "\n")
        
        console.print(f"[green]✓[/green] Report saved to: {report_file}")
        return report_file
        
    except IOError as e:
        console.print(f"[red]Error writing report file: {e}[/red]")
        raise


@click.command()
@click.option(
    '--tissue-pileup',
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help='Path to tissue pileup file (TSV or TSV.GZ format)'
)
@click.option(
    '--cfDNA-pileup',
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help='Path to cfDNA pileup file (TSV or TSV.GZ format)'
)
@click.option(
    '--min-snp-count',
    default=300,
    type=int,
    help='Minimum number of merged SNPs required for QC pass (default: 3000)'
)
@click.option(
    '--min-depth',
    default=10,
    type=int,
    help='Minimum sequencing depth threshold for filtering (default: 30)'
)
@click.option(
    '--output',
    required=True,
    type=str,
    help='Output file prefix (will create {prefix}_report.txt)'
)
def main(tissue_pileup: Path, cfdna_pileup: Path, min_snp_count: int,
         min_depth: int, output: str) -> None:
    """
    Calculate Pearson correlation between tissue and cfDNA SNP VAF values.
    
    This tool compares variant allele frequencies (VAF) between tissue and cfDNA samples
    to assess sample concordance. It filters pileup data by depth, merges on common SNP
    sites, calculates correlation, and performs quality control checks.
    
    QC Pass Criteria:
      - Merged SNP count >= min_snp_count
      - Pearson correlation in range [0.6, 0.8]
    
    Exit Codes:
      0: Success (QC Pass or Fail with valid results)
      1: Error during processing
    """
    console.print("\n[bold blue]SNP VAF Correlation Calculator[/bold blue]")
    console.print("=" * 70 + "\n")
    
    # Display input parameters
    params_table = Table(
        title="Input Parameters",
        show_header=True,
        header_style="bold magenta"
    )
    params_table.add_column("Parameter", style="cyan", no_wrap=True)
    params_table.add_column("Value", style="white")
    
    params_table.add_row("Tissue Pileup", str(tissue_pileup))
    params_table.add_row("cfDNA Pileup", str(cfdna_pileup))
    params_table.add_row("Min Depth", str(min_depth))
    params_table.add_row("Min SNP Count", f"{min_snp_count:,}")
    params_table.add_row("Output Prefix", output)
    
    console.print(params_table)
    console.print()
    
    try:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console
        ) as progress:
            # Read tissue pileup
            tissue_task = progress.add_task("Reading tissue pileup...", total=100)
            tissue_df = read_pileup_file(
                tissue_pileup, progress, tissue_task, "tissue"
            )
            progress.update(tissue_task, completed=100)
            
            # Read cfDNA pileup
            cfdna_task = progress.add_task("Reading cfDNA pileup...", total=100)
            cfDNA_df = read_pileup_file(
                cfdna_pileup, progress, cfdna_task, "cfDNA"
            )
            progress.update(cfdna_task, completed=100)
        
        console.print()
        
        # Filter by depth
        console.print(f"[bold yellow]Filtering by depth (>= {min_depth})...[/bold yellow]")
        tissue_filtered = filter_by_depth(tissue_df, min_depth, "tissue")
        cfDNA_filtered = filter_by_depth(cfDNA_df, min_depth, "cfDNA")
        
        console.print()
        
        # Merge pileup data
        console.print("[bold yellow]Merging pileup data...[/bold yellow]")
        merged_df = merge_pileup_data(tissue_filtered, cfDNA_filtered)
        
        console.print()
        
        # Calculate correlation
        console.print("[bold yellow]Calculating correlation...[/bold yellow]")
        correlation, p_value = calculate_correlation(merged_df)
        
        console.print()
        
        # Assess QC status
        qc_status, qc_reason = assess_qc_status(
            len(merged_df), correlation, min_snp_count
        )
        
        # Display QC results
        qc_color = "green" if qc_status == "PASS" else "red"
        qc_symbol = "✓" if qc_status == "PASS" else "✗"
        
        results_table = Table(
            title="Analysis Results",
            show_header=True,
            header_style="bold green"
        )
        results_table.add_column("Metric", style="cyan", no_wrap=True)
        results_table.add_column("Value", style="white", justify="right")
        
        results_table.add_row(
            "Tissue SNPs (filtered)",
            f"{len(tissue_filtered):,}"
        )
        results_table.add_row(
            "cfDNA SNPs (filtered)",
            f"{len(cfDNA_filtered):,}"
        )
        results_table.add_row(
            "Merged SNPs (common)",
            f"{len(merged_df):,}"
        )
        results_table.add_row("", "")  # Spacer
        results_table.add_row(
            "Pearson Correlation",
            f"{correlation:.6f}"
        )
        results_table.add_row(
            "P-value",
            f"{p_value:.6e}"
        )
        results_table.add_row("", "")  # Spacer
        results_table.add_row(
            f"[bold]QC Status[/bold]",
            f"[{qc_color}]{qc_symbol} {qc_status}[/{qc_color}]"
        )
        
        console.print(results_table)
        console.print()
        
        # Display QC reason in a panel
        panel_style = "green" if qc_status == "PASS" else "red"
        console.print(Panel(
            qc_reason,
            title="QC Assessment",
            style=panel_style,
            expand=False
        ))
        console.print()
        
        # Save report
        console.print("[bold yellow]Saving report...[/bold yellow]")
        report_file = save_report(
            output, tissue_pileup, cfdna_pileup,
            len(tissue_filtered), len(cfDNA_filtered), len(merged_df),
            correlation, p_value, qc_status, qc_reason,
            min_depth, min_snp_count
        )
        
        console.print(
            f"\n[bold green]{qc_symbol} Analysis completed successfully![/bold green]"
        )
        console.print(f"Report file: [cyan]{report_file}[/cyan]\n")
        
        # Exit with success (0) regardless of QC pass/fail
        # The QC status is informational only
        sys.exit(0)
        
    except Exception as e:
        console.print(f"\n[bold red]✗ Error during analysis:[/bold red] {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

