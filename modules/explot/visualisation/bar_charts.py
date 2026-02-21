#!/usr/bin/env python3
"""
Bar chart visualisation for SCEPTR ExPlot.

This module provides functions to create bar charts for visualising
functional or cellular category distributions across expression tiers.

Author: James McCabe
Module: SCEPTR ExPlot
"""

import os
import logging
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Any, Optional

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('sceptr.explot.visualisation.bar_charts')

def prepare_bar_data(
    tier_data: Dict[str, Any],
    method: str = 'multi'
) -> Tuple[List[str], List[float]]:
    """
    Prepare category percentages for bar chart display.
    
    Args:
        tier_data: Dictionary containing category counts and totals
        method: Percentage calculation method:
            - 'multi': Use sum of all category counts as denominator
            - 'single': Use total genes as denominator
        
    Returns:
        tuple: (categories, percentages)
    """
    counts = tier_data['counts']
    
    if method == 'multi':
        # Calculate total based on sum of category counts, which may exceed total genes due to multi-assignment
        total = sum(counts.values())
    else:
        # Calculate based on total genes
        total = tier_data.get('categorised', 0) + tier_data.get('uncategorised', 0)
    
    if total == 0:
        logger.warning("No genes categorised in this tier, chart will be empty")
        return list(counts.keys()), [0] * len(counts)
    
    # Calculate percentages
    percentages = []
    categories = []
    
    for category, count in counts.items():
        categories.append(category)
        percentages.append(round(count / total * 100, 1))
    
    return categories, percentages

def create_bar_charts(
    results: Dict[str, Dict[str, Any]],
    tier_names: List[str],
    categories: List[str],
    output_prefix: str,
    chart_type: str = "BP_MF",
    color_scheme: Optional[List[str]] = None,
    percentage_method: str = 'multi',
    sort_by: str = 'value',
    fig_size: Tuple[int, int] = (15, 8)
) -> List[str]:
    """
    Create improved bar charts with better readability and aesthetics.
    
    Args:
        results: Dictionary containing analysis results for each tier
        tier_names: List of tier names to include in the chart
        categories: List of category names to include in the chart
        output_prefix: Prefix for output file
        chart_type: Type of chart (BP_MF for functional or CC for cellular)
        color_scheme: Optional list of colors for bars
        percentage_method: Method to calculate percentages ('multi' or 'single')
        sort_by: How to sort bars: 'value' (highest to lowest) or 'name' (alphabetical)
        fig_size: Figure size as (width, height)
        
    Returns:
        list: Paths to the saved charts
    """
    logger.info("Creating bar charts...")
    
    # Define a better color palette if not provided
    if color_scheme is None:
        color_scheme = ['#2196f3', '#4caf50', '#ff9800', '#f44336', '#9c27b0', '#3f51b5', 
                      '#009688', '#ffeb3b', '#795548', '#607d8b']
    
    chart_paths = []
    
    for tier_name in tier_names:
        tier_data = results[tier_name]
        
        # Calculate percentages
        cats, percentages = prepare_bar_data(tier_data, method=percentage_method)
        
        # Create a mapping of categories to percentages
        cat_to_pct = {cat: pct for cat, pct in zip(cats, percentages)}
        
        # Get percentages in the order of the input categories
        ordered_cats = []
        ordered_pcts = []
        for cat in categories:
            if cat in cat_to_pct:
                ordered_cats.append(cat)
                ordered_pcts.append(cat_to_pct[cat])
        
        # Sort categories by percentage if requested
        if sort_by == 'value':
            sorted_data = sorted(zip(ordered_cats, ordered_pcts), key=lambda x: x[1], reverse=True)
            ordered_cats = [x[0] for x in sorted_data]
            ordered_pcts = [x[1] for x in sorted_data]
        elif sort_by == 'name':
            sorted_data = sorted(zip(ordered_cats, ordered_pcts), key=lambda x: x[0])
            ordered_cats = [x[0] for x in sorted_data]
            ordered_pcts = [x[1] for x in sorted_data]
        
        # Create figure with white background
        plt.figure(figsize=fig_size, facecolor='white')
        
        # Create bars with custom colors
        bars = plt.bar(
            ordered_cats, 
            ordered_pcts, 
            color=[color_scheme[i % len(color_scheme)] for i in range(len(ordered_cats))]
        )
        
        # Add a subtle grid for readability
        plt.grid(axis='y', linestyle='--', alpha=0.3)
        
        # Set chart properties with improved styling
        plt.xticks(rotation=45, ha='right', fontsize=10)
        
        if chart_type == "BP_MF":
            title = f"Functional Category Distribution - {tier_name}"
        else:
            title = f"Cellular Component Distribution - {tier_name}"
            
        # Add percentage calculation method to title
        if percentage_method == 'multi':
            title += " (% of all assignments)"
        else:
            title += " (% of all genes)"
            
        plt.title(title, fontsize=16, pad=20, fontweight='bold')
        plt.ylabel("Percentage (%)", fontsize=12, fontweight='bold')
        
        # Set y-axis limits with a little buffer at top
        max_pct = max(ordered_pcts) if ordered_pcts else 0
        plt.ylim(0, min(100, max_pct * 1.2))
        
        # Add value labels with better positioning
        for i, v in enumerate(ordered_pcts):
            if v > 0:  # Only add labels for non-zero values
                plt.text(i, v + (max_pct * 0.02), f"{v:.1f}%", 
                        ha='center', fontsize=9, fontweight='bold')
        
        # Customise background and remove top and right spines
        ax = plt.gca()
        ax.set_facecolor('#f8f9fa')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Add enrichment indicators if available
        if 'enrichment' in tier_data:
            # Add asterisks for significantly enriched categories
            for i, cat in enumerate(ordered_cats):
                if cat in tier_data['enrichment'] and tier_data['enrichment'][cat].get('significant', False):
                    plt.text(i, ordered_pcts[i] + (max_pct * 0.05), '*', 
                            ha='center', fontsize=16, fontweight='bold', color='red')
            
            # Add a note about significance
            plt.figtext(0.5, 0.01, "* Significantly enriched (adj. p < 0.05)", 
                       ha='center', fontsize=9, style='italic')
        
        # Save the chart with the output prefix
        plt.tight_layout()
        chart_path = f"{output_prefix}_{chart_type}_{tier_name}_bar.png"
        plt.savefig(chart_path, dpi=300, bbox_inches='tight', facecolor='white')
        svg_path = chart_path.replace('.png', '.svg')
        plt.savefig(svg_path, bbox_inches='tight', facecolor='white', format='svg')
        logger.info(f"Saved bar chart: {chart_path}, {svg_path}")
        plt.close()
        
        chart_paths.append(chart_path)
    
    return chart_paths

def create_grouped_bar_chart(
    results: Dict[str, Dict[str, Any]],
    tier_names: List[str],
    categories: List[str],
    output_prefix: str,
    chart_type: str = "BP_MF",
    color_scheme: Optional[List[str]] = None,
    percentage_method: str = 'multi',
    sort_by: str = 'value',
    fig_size: Tuple[int, int] = (18, 10)
) -> str:
    """
    Create a grouped bar chart comparing categories across tiers.
    
    Args:
        results: Dictionary containing analysis results for each tier
        tier_names: List of tier names to include in the chart
        categories: List of category names to include in the chart
        output_prefix: Prefix for output file
        chart_type: Type of chart (BP_MF for functional or CC for cellular)
        color_scheme: Optional list of colors for different tiers
        percentage_method: Method to calculate percentages ('multi' or 'single')
        sort_by: How to sort bars: 'value' (highest to lowest) or 'name' (alphabetical)
        fig_size: Figure size as (width, height)
        
    Returns:
        str: Path to the saved chart
    """
    logger.info("Creating grouped bar chart...")
    
    # Define a better color palette if not provided
    if color_scheme is None:
        color_scheme = ['#2196f3', '#4caf50', '#ff9800', '#f44336', '#9c27b0', '#3f51b5']
    
    # Calculate percentages for each tier
    tier_percentages = {}
    for tier_name in tier_names:
        tier_data = results[tier_name]
        cats, percentages = prepare_bar_data(tier_data, method=percentage_method)
        tier_percentages[tier_name] = {cat: pct for cat, pct in zip(cats, percentages)}
    
    # Ensure all categories are present in all tiers
    for tier_name in tier_names:
        for category in categories:
            if category not in tier_percentages[tier_name]:
                tier_percentages[tier_name][category] = 0
    
    # Sort categories if requested
    if sort_by == 'value':
        # Sort by the average percentage across tiers
        avg_percentages = {}
        for category in categories:
            avg_percentages[category] = sum(tier_percentages[tier][category] for tier in tier_names) / len(tier_names)
        sorted_categories = sorted(categories, key=lambda x: avg_percentages[x], reverse=True)
    elif sort_by == 'name':
        sorted_categories = sorted(categories)
    else:
        sorted_categories = categories
    
    # Create figure with white background
    plt.figure(figsize=fig_size, facecolor='white')
    
    # Set up bar positions
    bar_width = 0.8 / len(tier_names)
    index = np.arange(len(sorted_categories))
    
    # Create grouped bars
    for i, tier_name in enumerate(tier_names):
        tier_data = results[tier_name]
        percentages = [tier_percentages[tier_name][cat] for cat in sorted_categories]
        
        # Plot bars for this tier
        plt.bar(
            index + i * bar_width - (len(tier_names) - 1) * bar_width / 2, 
            percentages,
            bar_width,
            color=color_scheme[i % len(color_scheme)],
            label=tier_name
        )
    
    # Add a subtle grid for readability
    plt.grid(axis='y', linestyle='--', alpha=0.3)
    
    # Set chart properties with improved styling
    plt.xticks(index, sorted_categories, rotation=45, ha='right', fontsize=10)
    
    if chart_type == "BP_MF":
        title = "Functional Category Distribution Across Expression Tiers"
    else:
        title = "Cellular Component Distribution Across Expression Tiers"
        
    # Add percentage calculation method to title
    if percentage_method == 'multi':
        title += " (% of all assignments)"
    else:
        title += " (% of all genes)"
        
    plt.title(title, fontsize=16, pad=20, fontweight='bold')
    plt.ylabel("Percentage (%)", fontsize=12, fontweight='bold')
    plt.xlabel("Category", fontsize=12, fontweight='bold')
    
    # Set y-axis limits with a little buffer at top
    max_pct = max(
        max(tier_percentages[tier_name][cat] for cat in sorted_categories)
        for tier_name in tier_names
    )
    plt.ylim(0, min(100, max_pct * 1.2))
    
    # Customise background and remove top and right spines
    ax = plt.gca()
    ax.set_facecolor('#f8f9fa')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add legend
    plt.legend(loc='upper right', frameon=True, fontsize=10)
    
    # Add enrichment indicators if available
    enrichment_found = False
    for tier_name in tier_names:
        if 'enrichment' in results[tier_name]:
            enrichment_found = True
            break
    
    if enrichment_found:
        plt.figtext(0.5, 0.01, "* Categories with * are significantly enriched in at least one tier (adj. p < 0.05)", 
                   ha='center', fontsize=9, style='italic')
    
    # Save the chart
    plt.tight_layout()
    chart_path = f"{output_prefix}_{chart_type}_grouped_bar.png"
    plt.savefig(chart_path, dpi=300, bbox_inches='tight', facecolor='white')
    svg_path = chart_path.replace('.png', '.svg')
    plt.savefig(svg_path, bbox_inches='tight', facecolor='white', format='svg')
    logger.info(f"Saved grouped bar chart: {chart_path}, {svg_path}")
    plt.close()
    
    return chart_path

def create_enrichment_bar_chart(
    results: Dict[str, Dict[str, Any]],
    tier_name: str,
    output_prefix: str,
    chart_type: str = "BP_MF",
    p_value_threshold: float = 0.05,
    use_adjusted_p: bool = True,
    fig_size: Tuple[int, int] = (12, 8)
) -> str:
    """
    Create a bar chart showing fold enrichment for categories.
    
    Args:
        results: Dictionary containing analysis results for each tier
        tier_name: Name of the tier to visualise
        output_prefix: Prefix for output file
        chart_type: Type of chart (BP_MF for functional or CC for cellular)
        p_value_threshold: P-value threshold for significance
        use_adjusted_p: Whether to use adjusted p-values
        fig_size: Figure size as (width, height)
        
    Returns:
        str: Path to the saved chart
    """
    logger.info(f"Creating enrichment bar chart for {tier_name}...")
    
    # Get enrichment data for this tier
    tier_data = results[tier_name]
    if 'enrichment' not in tier_data:
        logger.error(f"No enrichment data found for {tier_name}")
        return None
    
    enrichment_data = tier_data['enrichment']
    
    # Extract fold changes and p-values
    categories = []
    fold_changes = []
    p_values = []
    significant = []
    
    for category, stats in enrichment_data.items():
        if stats['fold_change'] == 0:
            continue  # Skip categories with no fold change
            
        categories.append(category)
        fold_changes.append(stats['fold_change'])
        
        # Get p-value (adjusted or raw)
        if use_adjusted_p and 'adjusted_p_value' in stats:
            p_val = stats['adjusted_p_value']
        else:
            p_val = stats['p_value']
            
        p_values.append(p_val)
        significant.append(p_val < p_value_threshold)
    
    # Sort categories by fold change
    sorted_data = sorted(zip(categories, fold_changes, p_values, significant), 
                        key=lambda x: x[1], reverse=True)
    
    if not sorted_data:
        logger.warning(f"No enrichment data to plot for {tier_name}")
        return None
        
    categories = [x[0] for x in sorted_data]
    fold_changes = [x[1] for x in sorted_data]
    p_values = [x[2] for x in sorted_data]
    significant = [x[3] for x in sorted_data]
    
    # Create figure with white background
    plt.figure(figsize=fig_size, facecolor='white')
    
    # Create bars with color based on significance
    bar_colors = ['#e53935' if sig else '#2196f3' for sig in significant]
    bars = plt.bar(categories, fold_changes, color=bar_colors)
    
    # Add a subtle grid for readability
    plt.grid(axis='y', linestyle='--', alpha=0.3)
    
    # Add horizontal line at fold change = 1 (no enrichment)
    plt.axhline(y=1, color='gray', linestyle='--', alpha=0.7)
    
    # Set chart properties with improved styling
    plt.xticks(rotation=45, ha='right', fontsize=10)
    
    if chart_type == "BP_MF":
        title = f"Functional Category Enrichment - {tier_name}"
    else:
        title = f"Cellular Component Enrichment - {tier_name}"
    
    p_value_type = "Adjusted p-value" if use_adjusted_p else "p-value"
    plt.title(title, fontsize=16, pad=20, fontweight='bold')
    plt.ylabel("Fold Enrichment", fontsize=12, fontweight='bold')
    
    # Add value labels and significance markers
    for i, (v, sig) in enumerate(zip(fold_changes, significant)):
        label_text = f"{v:.1f}x"
        if sig:
            label_text += "*"
            
        plt.text(i, v + 0.1, label_text, ha='center', fontsize=9, fontweight='bold')
    
    # Customise background and remove top and right spines
    ax = plt.gca()
    ax.set_facecolor('#f8f9fa')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add a note about significance
    plt.figtext(0.5, 0.01, f"* Significant enrichment ({p_value_type} < {p_value_threshold})", 
               ha='center', fontsize=9, style='italic')
    
    # Save the chart
    plt.tight_layout()
    chart_path = f"{output_prefix}_{chart_type}_{tier_name}_enrichment.png"
    plt.savefig(chart_path, dpi=300, bbox_inches='tight', facecolor='white')
    svg_path = chart_path.replace('.png', '.svg')
    plt.savefig(svg_path, bbox_inches='tight', facecolor='white', format='svg')
    logger.info(f"Saved enrichment bar chart: {chart_path}, {svg_path}")
    plt.close()
    
    return chart_path

def create_multi_tier_enrichment_bar_chart(
    results: Dict[str, Dict[str, Any]],
    tier_names: List[str],
    output_prefix: str,
    chart_type: str = "BP_MF",
    p_value_threshold: float = 0.05,
    use_adjusted_p: bool = True,
    fig_size: Tuple[int, int] = (14, 7)
) -> Optional[str]:
    """
    Create a grouped bar chart showing fold enrichment across all tiers.

    Categories on x-axis, one bar per tier per category, color-coded by tier.
    Bars with FDR < 0.05 get full opacity; non-significant get alpha=0.4.

    Args:
        results: Dictionary containing analysis results for each tier
        tier_names: List of tier names (e.g. ['top_50', 'top_100', 'top_250', 'top_500'])
        output_prefix: Prefix for output file
        chart_type: Type of chart (BP_MF or CC)
        p_value_threshold: P-value threshold for significance
        use_adjusted_p: Whether to use adjusted p-values
        fig_size: Figure size as (width, height)

    Returns:
        str: Path to the saved chart, or None if no data
    """
    logger.info("Creating multi-tier enrichment bar chart...")

    # Tier colors matching radar chart palette
    tier_colors = {
        'top_50': '#1f77b4',    # Blue
        'top_100': '#ff7f0e',   # Orange
        'top_250': '#2ca02c',   # Green
        'top_500': '#d62728'    # Red
    }

    # Collect fold-change data across all tiers
    all_categories = set()
    tier_data_map = {}

    for tier_name in tier_names:
        tier_data = results.get(tier_name, {})
        if 'enrichment' not in tier_data:
            logger.warning(f"No enrichment data for {tier_name}, skipping")
            continue

        enrichment = tier_data['enrichment']
        tier_entry = {}
        for category, stats in enrichment.items():
            fc = stats.get('fold_change', 0)
            if fc == 0:
                continue
            if use_adjusted_p and 'adjusted_p_value' in stats:
                p_val = stats['adjusted_p_value']
            else:
                p_val = stats.get('p_value', 1.0)
            tier_entry[category] = {
                'fold_change': fc,
                'significant': p_val < p_value_threshold
            }
            all_categories.add(category)
        tier_data_map[tier_name] = tier_entry

    if not all_categories or not tier_data_map:
        logger.warning("No enrichment data to plot for multi-tier chart")
        return None

    # Sort categories by fold-change at the first tier (descending)
    first_tier = tier_names[0] if tier_names[0] in tier_data_map else list(tier_data_map.keys())[0]
    sorted_categories = sorted(
        all_categories,
        key=lambda c: tier_data_map.get(first_tier, {}).get(c, {}).get('fold_change', 0),
        reverse=True
    )

    # Set up figure
    plt.style.use('default')
    plt.rcParams.update({
        'font.size': 11,
        'axes.linewidth': 1.0,
        'figure.facecolor': 'white',
        'axes.facecolor': 'white'
    })

    fig, ax = plt.subplots(figsize=fig_size, facecolor='white')

    active_tiers = [t for t in tier_names if t in tier_data_map]
    n_tiers = len(active_tiers)
    n_cats = len(sorted_categories)
    bar_width = 0.8 / n_tiers
    index = np.arange(n_cats)

    # Plot grouped bars
    for i, tier_name in enumerate(active_tiers):
        tier_entry = tier_data_map[tier_name]
        color = tier_colors.get(tier_name, f'C{i}')

        fold_changes = []
        significances = []
        for cat in sorted_categories:
            entry = tier_entry.get(cat, {'fold_change': 0, 'significant': False})
            fold_changes.append(entry['fold_change'])
            significances.append(entry['significant'])

        x_positions = index + i * bar_width - (n_tiers - 1) * bar_width / 2
        alphas = [1.0 if sig else 0.4 for sig in significances]

        # Plot each bar individually to control alpha
        for j, (x, fc, alpha, sig) in enumerate(zip(x_positions, fold_changes, alphas, significances)):
            ax.bar(x, fc, bar_width, color=color, alpha=alpha,
                   label=tier_name.replace('_', ' ').title() if j == 0 else None,
                   edgecolor='white', linewidth=0.5)

            # Add fold-change label above bar
            if fc > 0:
                label_text = f"{fc:.1f}x"
                if sig:
                    label_text += "*"
                ax.text(x, fc + 0.05, label_text, ha='center', va='bottom',
                        fontsize=7, fontweight='bold' if sig else 'normal',
                        color='#333333')

    # Reference line at fold-change = 1.0
    ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.7, linewidth=1.0, zorder=0)

    # Styling
    ax.set_xticks(index)
    ax.set_xticklabels(sorted_categories, rotation=45, ha='right', fontsize=10)

    if chart_type == "BP_MF":
        title = "Functional Category Enrichment Across Expression Tiers"
    else:
        title = "Cellular Component Enrichment Across Expression Tiers"

    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
    ax.set_ylabel("Fold Enrichment", fontsize=12, fontweight='bold')

    # Y-axis: buffer above max
    all_fc = [tier_data_map[t].get(c, {}).get('fold_change', 0)
              for t in active_tiers for c in sorted_categories]
    max_fc = max(all_fc) if all_fc else 2.0
    ax.set_ylim(0, max_fc * 1.25)

    # Grid and spines
    ax.grid(axis='y', linestyle='--', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Legend
    ax.legend(loc='upper right', frameon=True, fontsize=10,
              facecolor='white', edgecolor='#cccccc')

    # Significance note
    p_label = "FDR" if use_adjusted_p else "p-value"
    fig.text(0.5, 0.01,
             f"* Significant enrichment ({p_label} < {p_value_threshold}); "
             f"faded bars = not significant",
             ha='center', fontsize=9, style='italic', color='#666666')

    # Save
    plt.tight_layout(rect=[0, 0.03, 1, 1])
    chart_path = f"{output_prefix}_{chart_type}_multi_tier_enrichment.png"
    svg_path = chart_path.replace('.png', '.svg')

    plt.savefig(chart_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(svg_path, bbox_inches='tight', facecolor='white', format='svg')
    logger.info(f"Saved multi-tier enrichment bar chart: {chart_path}, {svg_path}")
    plt.close()

    return chart_path


def create_multi_category_overlap_chart(
    results: Dict[str, Dict[str, Any]],
    tier_name: str,
    output_prefix: str,
    chart_type: str = "BP_MF",
    max_categories: int = 10,
    fig_size: Tuple[int, int] = (14, 8)
) -> str:
    """
    Create a bar chart showing the overlap between categories.
    
    Args:
        results: Dictionary containing analysis results for each tier
        tier_name: Name of the tier to visualise
        output_prefix: Prefix for output file
        chart_type: Type of chart (BP_MF for functional or CC for cellular)
        max_categories: Maximum number of categories to show
        fig_size: Figure size as (width, height)
        
    Returns:
        str: Path to the saved chart
    """
    # This is a placeholder for a future UpSet-like plot
    # that shows the overlaps between categories
    # For now, let's create a simple bar chart of category counts
    
    logger.info(f"Creating multi-category overlap chart for {tier_name}...")
    
    # Get data for this tier
    tier_data = results[tier_name]
    if 'genes_in_category' not in tier_data:
        logger.error(f"No category assignment data found for {tier_name}")
        return None
    
    category_counts = tier_data['counts']
    
    # Sort categories by count
    sorted_data = sorted(category_counts.items(), key=lambda x: x[1], reverse=True)
    
    # Limit to max_categories
    if len(sorted_data) > max_categories:
        sorted_data = sorted_data[:max_categories]
        
    categories = [x[0] for x in sorted_data]
    counts = [x[1] for x in sorted_data]
    
    # Create figure with white background
    plt.figure(figsize=fig_size, facecolor='white')
    
    # Create bars
    bars = plt.bar(categories, counts, color='#2196f3')
    
    # Add a subtle grid for readability
    plt.grid(axis='y', linestyle='--', alpha=0.3)
    
    # Set chart properties with improved styling
    plt.xticks(rotation=45, ha='right', fontsize=10)
    
    title = f"Category Gene Counts - {tier_name}"
    plt.title(title, fontsize=16, pad=20, fontweight='bold')
    plt.ylabel("Number of Genes", fontsize=12, fontweight='bold')
    
    # Add value labels
    for i, v in enumerate(counts):
        plt.text(i, v + 0.5, str(v), ha='center', fontsize=9, fontweight='bold')
    
    # Customise background and remove top and right spines
    ax = plt.gca()
    ax.set_facecolor('#f8f9fa')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add a note about the multi-category nature
    multi_category_genes = tier_data.get('multi_category_genes', 0)
    total_categorised = tier_data.get('categorised', 0)
    
    if total_categorised > 0:
        multi_category_pct = (multi_category_genes / total_categorised) * 100
        plt.figtext(0.5, 0.01, 
                   f"Note: {multi_category_genes} genes ({multi_category_pct:.1f}%) are assigned to multiple categories", 
                   ha='center', fontsize=9, style='italic')
    
    # Save the chart
    plt.tight_layout()
    chart_path = f"{output_prefix}_{chart_type}_{tier_name}_category_counts.png"
    plt.savefig(chart_path, dpi=300, bbox_inches='tight', facecolor='white')
    svg_path = chart_path.replace('.png', '.svg')
    plt.savefig(svg_path, bbox_inches='tight', facecolor='white', format='svg')
    logger.info(f"Saved category counts chart: {chart_path}, {svg_path}")
    plt.close()
    
    return chart_path
