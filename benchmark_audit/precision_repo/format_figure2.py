import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    os.makedirs("./images", exist_ok=True)
    
    scores_path = "./Evaluation/XAI/occlusion_scores.pkl"
    if not os.path.exists(scores_path):
        raise FileNotFoundError(f"Occlusion scores not found at {scores_path}")
        
    with open(scores_path, "rb") as f:
        data = pickle.load(f)
        
    # Rank 1 entry (Influenza LPRRSGAAGA)
    entry = data[0]
    peptide = entry["peptide"]
    baseline_prob = entry["baseline_prob"]
    hla_allele = entry["hla_allele"]
    
    print(f"Loading Rank 1: Peptide={peptide}, Prob={baseline_prob:.4f}")
    
    # We want strictly the peptide sequence characters and scores
    # Peptide length is 10
    pep_chars = entry["pephla_chars"][:10]
    scores_w1 = entry["pephla_scores_w1"][:10]
    scores_w3 = entry["pephla_scores_w3"][:10]
    
    matrix = np.array([scores_w1, scores_w3])
    
    # Set up styling
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    
    fig, ax = plt.subplots(figsize=(3.5, 2.2))
    
    # Draw heatmap
    sns.heatmap(
        matrix,
        xticklabels=pep_chars,
        yticklabels=["W=1", "W=3"],
        cmap="Reds",
        ax=ax,
        cbar_kws={"label": r"Sensitivity ($\Delta$P)", "shrink": 0.8},
        vmin=min(0.0, matrix.min()),
        vmax=max(0.1, matrix.max()),
        linewidths=0.5,
        linecolor='grey'
    )
    
    # Make labels large and readable for 3.5-inch column scaling
    ax.set_xticklabels(pep_chars, rotation=0, fontsize=12, fontweight='bold', fontfamily='monospace')
    ax.set_yticklabels(["W=1", "W=3"], rotation=0, fontsize=10)
    
    # Title
    ax.set_title(f"Influenza Peptide: {peptide}\n(Prob: {baseline_prob:.3f}, HLA: {hla_allele})", 
                 fontsize=9, fontweight='bold', pad=10)
    
    plt.tight_layout()
    
    # Save files
    png_path = "./images/Figure2_OcclusionMap.png"
    pdf_path = "./images/Figure2_OcclusionMap.pdf"
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight')
    plt.close()
    
    print(f"Successfully generated Figure 2 at {png_path} and {pdf_path}")

if __name__ == "__main__":
    main()
