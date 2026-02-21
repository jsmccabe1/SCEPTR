#!/usr/bin/env python3
"""
Enhanced radar chart visualisation for SCEPTR ExPlot.

This module provides functions to create publication-quality radar (spider) charts 
for visualising functional or cellular category distributions across expression tiers.

Author: James McCabe
Module: SCEPTR ExPlot
"""

import os
import logging
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import font_manager
from typing import Dict, List, Tuple, Any, Optional

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('sceptr.explot.visualisation.radar_charts')

def prepare_radar_data(
    tier_data: Dict[str, Any],
    method: str = 'multi'
) -> Tuple[List[str], List[float]]:
    """
    Prepare category percentages for radar chart display.
    
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

def create_publication_radar_chart(
    results: Dict[str, Dict[str, Any]],
    tier_names: List[str],
    categories: List[str],
    output_prefix: str,
    chart_type: str = "BP_MF",
    fig_size: Tuple[int, int] = (10, 10),
    dpi: int = 300
) -> Tuple[str, str]:
    """
    Create a publication-quality radar chart matching the target aesthetic.
    
    Args:
        results: Dictionary containing analysis results for each tier
        tier_names: List of tier names to include in the chart
        categories: List of category names to include in the chart
        output_prefix: Prefix for output file
        chart_type: Type of chart (BP_MF for functional or CC for cellular)
        fig_size: Figure size as (width, height)
        dpi: Resolution for output files
        
    Returns:
        tuple: (png_path, svg_path)
    """
    logger.info("Creating publication-quality radar chart...")
    
    # Set up clean, publication-ready style
    plt.style.use('default')
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 11,
        'axes.linewidth': 1.0,
        'axes.edgecolor': '#000000',
        'grid.linewidth': 0.8,
        'grid.alpha': 0.3,
        'figure.facecolor': 'white',
        'axes.facecolor': 'white'
    })
    
    num_cats = len(categories)
    
    # Set up angles for radar chart
    angles = np.linspace(0, 2 * np.pi, num_cats, endpoint=False).tolist()
    angles += angles[:1]  # Close the loop
    
    # Create figure with clean white background
    fig = plt.figure(figsize=fig_size, facecolor='white', dpi=dpi)
    ax = fig.add_subplot(111, polar=True, facecolor='white')
    
    # Publication-quality color scheme - distinct, colorblind-friendly colors
    tier_colors = {
        'top_50': '#1f77b4',    # Blue
        'top_100': '#ff7f0e',   # Orange  
        'top_250': '#2ca02c',   # Green
        'top_500': '#d62728'    # Red
    }
    
    # Calculate max percentage for consistent scaling
    max_percentage = 0
    for tier_name in tier_names:
        tier_data = results[tier_name]
        _, percentages = prepare_radar_data(tier_data, method='multi')
        max_percentage = max(max_percentage, max(percentages) if percentages else 0)
    
    # Round up to nearest 20 for clean scaling
    max_percentage = max(20, ((max_percentage // 20) + 1) * 20)
    max_percentage = min(100, max_percentage)  # Cap at 100%
    
    # Store data and plot each tier
    for i, tier_name in enumerate(tier_names):
        tier_data = results[tier_name]
        cats, percentages = prepare_radar_data(tier_data, method='multi')
        
        # Create a mapping from category to percentage
        cat_to_pct = {cat: pct for cat, pct in zip(cats, percentages)}
        
        # Reorder percentages to match the input categories order
        ordered_percentages = [cat_to_pct.get(cat, 0) for cat in categories]
        
        # Close the loop for plotting
        percentages_closed = ordered_percentages + [ordered_percentages[0]]
        
        # Get color for this tier
        color = tier_colors.get(tier_name, plt.cm.tab10(i))
        
        # Plot line with clean styling
        line = ax.plot(angles, percentages_closed, 
                      linewidth=2.5, 
                      color=color, 
                      marker='o', 
                      markersize=6,
                      markerfacecolor=color,
                      markeredgecolor='white',
                      markeredgewidth=1,
                      label=tier_name.replace('_', ' ').title(),
                      zorder=10-i)  # Higher tiers on top
        
        # Add subtle semi-transparent fill
        ax.fill(angles, percentages_closed, 
               color=color, 
               alpha=0.15, 
               zorder=i)
    
    # Clean up the chart appearance
    ax.set_theta_offset(np.pi / 2)  # Start from top
    ax.set_theta_direction(-1)  # Clockwise
    
    # Set category labels with clean formatting
    ax.set_xticks(angles[:-1])
    
    # Format category names for better display
    formatted_categories = []
    for cat in categories:
        # Split long category names for better display
        if len(cat) > 20:
            words = cat.split()
            mid = len(words) // 2
            cat = ' '.join(words[:mid]) + '\n' + ' '.join(words[mid:])
        formatted_categories.append(cat)
    
    ax.set_xticklabels(formatted_categories, 
                       fontsize=10, 
                       fontweight='normal',
                       color='#333333')
    
    # Set radial ticks - clean, minimal grid
    yticks = list(range(0, int(max_percentage) + 20, 20))
    ax.set_yticks(yticks)
    ax.set_yticklabels([f'{y}%' for y in yticks], 
                       fontsize=9, 
                       color='#666666',
                       alpha=0.8)
    ax.set_rlabel_position(0)  # Position y-labels at 0 degrees
    
    # Set limits
    ax.set_ylim(0, max_percentage)
    
    # Clean grid styling
    ax.grid(True, 
           linestyle='-', 
           alpha=0.3, 
           linewidth=0.8,
           color='#666666')
    ax.set_facecolor('white')
    
    # Add title with publication formatting
    if chart_type == "BP_MF":
        title = "Gene Expression by Biological Process and Molecular Function"
    else:
        title = "Gene Expression by Cellular Component"
    
    ax.set_title(title, 
                fontsize=14, 
                fontweight='bold', 
                pad=20, 
                color='#000000')
    
    # Clean legend styling - positioned outside plot area
    legend = ax.legend(loc='center left', 
                      bbox_to_anchor=(1.15, 0.5),
                      frameon=True, 
                      fontsize=11,
                      facecolor='white',
                      edgecolor='#cccccc',
                      framealpha=1.0,
                      columnspacing=0.8,
                      handletextpad=0.5)
    
    # Style legend
    legend.get_frame().set_linewidth(1)
    
    # Add percentage method note at bottom
    plt.figtext(0.5, 0.02, 
               "Percentages based on sum of all category assignments", 
               ha='center', 
               fontsize=8, 
               color='#666666',
               style='italic')
    
    # Adjust layout for clean appearance
    plt.tight_layout()
    
    # Save both PNG and SVG versions
    png_path = f"{output_prefix}_{chart_type}_radar.png"
    svg_path = f"{output_prefix}_{chart_type}_radar.svg"
    
    # Save PNG with high quality
    plt.savefig(png_path, 
               dpi=dpi, 
               bbox_inches='tight', 
               facecolor='white',
               edgecolor='none',
               format='png')
    
    # Save SVG for vector editing
    plt.savefig(svg_path,
               bbox_inches='tight',
               facecolor='white',
               edgecolor='none', 
               format='svg')
    
    logger.info(f"Saved radar chart to: {png_path}")
    logger.info(f"Saved radar chart (SVG) to: {svg_path}")
    
    plt.close()
    
    return png_path, svg_path

def create_individual_tier_charts_publication(
    results: Dict[str, Dict[str, Any]],
    tier_names: List[str],
    categories: List[str],
    output_prefix: str,
    chart_type: str = "BP_MF"
) -> List[Tuple[str, str]]:
    """
    Create individual publication-quality radar charts for each tier.
    
    Args:
        results: Dictionary containing analysis results for each tier
        tier_names: List of tier names to include in the chart
        categories: List of category names to include in the chart
        output_prefix: Prefix for output file
        chart_type: Type of chart (BP_MF for functional or CC for cellular)
        
    Returns:
        list: List of tuples (png_path, svg_path) for each tier chart
    """
    chart_paths = []
    
    for tier_name in tier_names:
        # Create a results dict with just this one tier
        single_tier_results = {tier_name: results[tier_name]}
        
        png_path, svg_path = create_publication_radar_chart(
            single_tier_results, 
            [tier_name], 
            categories, 
            f"{output_prefix}_{tier_name}", 
            chart_type
        )
        chart_paths.append((png_path, svg_path))
    
    return chart_paths

def create_multi_panel_figure(
    results_dict: Dict[str, Dict[str, Dict[str, Any]]],
    categories: List[str],
    output_prefix: str,
    chart_type: str = "BP_MF",
    fig_size: Tuple[int, int] = (15, 5)
) -> Tuple[str, str]:
    """
    Create a multi-panel radar chart figure with multiple organisms/datasets.
    Similar to the target publication figure with multiple radar charts side by side.
    
    Args:
        results_dict: Dictionary mapping organism names to results dictionaries
        categories: List of category names
        output_prefix: Prefix for output files
        chart_type: Type of chart
        fig_size: Figure size
        
    Returns:
        tuple: (png_path, svg_path)
    """
    logger.info("Creating multi-panel publication figure...")
    
    # Set up publication style
    plt.style.use('default')
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 10,
        'axes.linewidth': 1.0,
        'figure.facecolor': 'white'
    })
    
    n_panels = len(results_dict)
    fig, axes = plt.subplots(1, n_panels, figsize=fig_size, 
                            subplot_kw=dict(projection='polar'),
                            facecolor='white')
    
    if n_panels == 1:
        axes = [axes]  # Make it iterable
    
    # Color scheme for tiers
    tier_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    for panel_idx, (organism_name, results) in enumerate(results_dict.items()):
        ax = axes[panel_idx]
        
        # Set up angles
        num_cats = len(categories)
        angles = np.linspace(0, 2 * np.pi, num_cats, endpoint=False).tolist()
        angles += angles[:1]
        
        # Calculate max for this organism
        max_percentage = 0
        tier_names = list(results.keys())
        
        for tier_name in tier_names:
            _, percentages = prepare_radar_data(results[tier_name], method='multi')
            max_percentage = max(max_percentage, max(percentages) if percentages else 0)
        
        max_percentage = max(20, ((max_percentage // 20) + 1) * 20)
        max_percentage = min(100, max_percentage)
        
        # Plot each tier
        for i, tier_name in enumerate(tier_names):
            tier_data = results[tier_name]
            cats, percentages = prepare_radar_data(tier_data, method='multi')
            
            # Reorder percentages
            cat_to_pct = {cat: pct for cat, pct in zip(cats, percentages)}
            ordered_percentages = [cat_to_pct.get(cat, 0) for cat in categories]
            percentages_closed = ordered_percentages + [ordered_percentages[0]]
            
            color = tier_colors[i % len(tier_colors)]
            
            # Plot with clean styling
            ax.plot(angles, percentages_closed,
                   linewidth=2,
                   color=color,
                   marker='o',
                   markersize=4,
                   alpha=0.8,
                   label=tier_name.replace('_', ' ').title())
            
            ax.fill(angles, percentages_closed,
                   color=color,
                   alpha=0.1)
        
        # Format this subplot
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        
        # Category labels
        ax.set_xticks(angles[:-1])
        formatted_categories = [cat.replace(' & ', '\n& ') if ' & ' in cat else cat 
                               for cat in categories]
        ax.set_xticklabels(formatted_categories, fontsize=8, color='#333333')
        
        # Radial labels
        yticks = list(range(0, int(max_percentage) + 20, 20))
        ax.set_yticks(yticks)
        ax.set_yticklabels([f'{y}%' if y > 0 else '' for y in yticks], 
                          fontsize=8, color='#666666', alpha=0.7)
        ax.set_ylim(0, max_percentage)
        
        # Grid
        ax.grid(True, linestyle='-', alpha=0.3, linewidth=0.5)
        ax.set_facecolor('white')
        
        # Organism title
        ax.set_title(f"{organism_name}", fontsize=12, fontweight='bold', pad=15)
        
        # Legend only on first panel
        if panel_idx == 0:
            legend = ax.legend(loc='center left', bbox_to_anchor=(1.1, 0.5),
                             frameon=True, fontsize=9, facecolor='white',
                             edgecolor='#cccccc')
    
    # Overall title
    if chart_type == "BP_MF":
        main_title = "Parasite Gene Expression by Biological Process and Molecular Function"
    else:
        main_title = "Parasite Gene Expression by Cellular Component"
    
    fig.suptitle(main_title, fontsize=16, fontweight='bold', y=0.95)
    
    plt.tight_layout()
    
    # Save files
    png_path = f"{output_prefix}_{chart_type}_multi_panel.png"
    svg_path = f"{output_prefix}_{chart_type}_multi_panel.svg"
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(svg_path, bbox_inches='tight', facecolor='white', format='svg')
    
    logger.info(f"Saved multi-panel figure: {png_path}")
    logger.info(f"Saved multi-panel figure (SVG): {svg_path}")
    
    plt.close()
    
    return png_path, svg_path

# Update the main create_radar_chart function to use the new publication version
def create_radar_chart(
    results: Dict[str, Dict[str, Any]],
    tier_names: List[str],
    categories: List[str],
    output_prefix: str,
    chart_type: str = "BP_MF",
    min_percentage: float = 5.0,
    color_scheme: Optional[List[str]] = None,
    normalise: bool = False,
    percentage_method: str = 'multi',
    fig_size: Tuple[int, int] = (12, 10)
) -> str:
    """
    Main radar chart creation function - now uses publication quality version.
    
    This is kept for backwards compatibility but now calls the publication version.
    """
    png_path, svg_path = create_publication_radar_chart(
        results, tier_names, categories, output_prefix, chart_type, fig_size
    )
    return png_path  # Return PNG path for backwards compatibility
