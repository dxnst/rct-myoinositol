#!/usr/bin/env python3
"""
Script to visualize the improvements made to the forest plot
"""

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

def display_improvements():
    """
    Display the improved forest plot with annotations showing the improvements
    """
    # Load the improved forest plot
    img = mpimg.imread('ejemplo_forest_plots/forest_plot_combinado_bmi_kg_m2.png')
    
    fig, ax = plt.subplots(figsize=(18, 12))
    ax.imshow(img)
    ax.axis('off')
    
    # Add title
    ax.set_title('Enhanced Forest Plot with Improvements', fontsize=20, fontweight='bold', pad=30)
    
    # Add improvement annotations
    improvements = [
        "✓ Absolute Hedges' g values (|g|) for clarity",
        "✓ Better vertical spacing between obesity groups", 
        "✓ Clear group headers: 'Peso normal' and 'Sobrepeso/Obesidad'",
        "✓ Preserved signs in change values (Δ)",
        "✓ Improved text positioning to avoid overlap",
        "✓ Enhanced readability with better group separation"
    ]
    
    # Add text box with improvements
    improvement_text = '\n'.join(improvements)
    ax.text(0.02, 0.98, improvement_text, transform=ax.transAxes, 
            verticalalignment='top', fontsize=12, 
            bbox=dict(boxstyle="round,pad=0.5", facecolor='lightgreen', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('forest_plot_improvements_preview.png', dpi=200, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    display_improvements()