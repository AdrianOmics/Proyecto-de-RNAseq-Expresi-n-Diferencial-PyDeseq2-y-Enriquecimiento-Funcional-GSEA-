#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Firmas DPN 48h y GSEA sobre ranking DPN 24h (log2FC)
----------------------------------------------------

1) A partir de los resultados DESeq2 (DPN vs Control, 48h):
   - Generar dos firmas:
     - DPN48_UP  : 100 genes más sobreexpresados (log2FC > 0)
     - DPN48_DOWN: 100 genes más infraexpresados (log2FC < 0)

2) A partir de los resultados DESeq2 (DPN vs Control, 24h):
   - Construir un ranking de genes por log2FoldChange.

3) Ejecutar GSEA preranked (gseapy.prerank):
   - gene_sets = {'DPN48_UP': [...], 'DPN48_DOWN': [...]}
   - ranking = resultados 24h (serie ordenada por log2FoldChange)

4) Guardar:
   - tablas de resultados GSEA
   - enrichment plots
   - dotplot conjunto (NES vs Term) para ambas firmas
   - dotplot de enriquecimiento individual de genes
"""

from pathlib import Path
import pandas as pd
import numpy as np
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

sns.set_theme(style="whitegrid", palette="deep")

# =============================================================================
# 1. CONFIGURACIÓN DE RUTAS Y PARÁMETROS
# =============================================================================

working_dir = Path.home()

# Directorio donde guardaste los resultados DESeq2 (ajusta si es distinto)
DEG_dir = working_dir / "Documentos/Transcriptomica/Trabajo_final_RNAseq/Deseq2"

# Ficheros de resultados DESeq2 (AJUSTAR NOMBRES SI SON OTROS)
res_48h_path = DEG_dir / "pydeseq2_parathyroid_DEG/deseq2_results_DPN_vs_Control_48h.csv"
res_24h_path = DEG_dir / "pydeseq2_parathyroid_DEG/deseq2_results_DPN_vs_Control_24h.csv"

# Carpeta de salida para firmas y GSEA
output_dir = working_dir / "Documentos/Transcriptomica/Trabajo_final_RNAseq/GSEApr"
output_dir.mkdir(parents=True, exist_ok=True)

# Número de genes por firma
N_TOP = 100

# Parámetros de GSEA preranked
N_PERMUTATIONS = 1000
MIN_SIZE = 5
MAX_SIZE = 500
SEED = 7

# =============================================================================
# 2. GENERAR FIRMAS A PARTIR DE DPN vs CONTROL 48h
# =============================================================================

def generate_dpn48_signatures(results_48h_path: Path, n_top: int = 100):
    """
    A partir de la tabla DESeq2 (DPN_vs_Control_48h), generar:
      - top N_UP genes (log2FC > 0)
      - top N_DOWN genes (log2FC < 0)
    
    Criterio:
      - Filtrar por log2FoldChange > 0 ó < 0
      - Ordenar por padj ascendente y, dentro, por log2FC en la dirección adecuada
    """
    if not results_48h_path.exists():
        raise FileNotFoundError(f"No se encontró el archivo de resultados 48h: {results_48h_path}")
    
    logger.info(f"Cargando resultados 48h desde: {results_48h_path}")
    res48 = pd.read_csv(results_48h_path, index_col=0)
    
    if "log2FoldChange" not in res48.columns or "padj" not in res48.columns:
        raise ValueError("El fichero 48h debe contener las columnas 'log2FoldChange' y 'padj'.")
    
    res48.index.name = "gene"
    
    # Limpiar NaN
    res48_clean = res48.dropna(subset=["log2FoldChange", "padj"]).copy()
    
    # --- UP: genes más expresados en DPN (log2FC > 0) ---
    up = res48_clean[res48_clean["log2FoldChange"] > 0].copy()
    up_sorted = up.sort_values(
        by=["padj", "log2FoldChange"],
        ascending=[True, False],
    )
    top_up = up_sorted.head(n_top)
    
    # --- DOWN: genes menos expresados en DPN (log2FC < 0) ---
    down = res48_clean[res48_clean["log2FoldChange"] < 0].copy()
    # Orden: padj ascendente, log2FC ascendente (más negativo primero)
    down_sorted = down.sort_values(
        by=["padj", "log2FoldChange"],
        ascending=[True, True],
    )
    top_down = down_sorted.head(n_top)
    
    logger.info(f"Nº genes UP disponibles 48h: {up.shape[0]} | seleccionados: {top_up.shape[0]}")
    logger.info(f"Nº genes DOWN disponibles 48h: {down.shape[0]} | seleccionados: {top_down.shape[0]}")
    
    # Guardar tablas de firmas
    cols_to_save = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    cols_to_save = [c for c in cols_to_save if c in res48_clean.columns]
    
    up_csv = output_dir / f"firmas_DPN48h_up_{n_top}.csv"
    down_csv = output_dir / f"firmas_DPN48h_down_{n_top}.csv"
    
    top_up[cols_to_save].to_csv(up_csv)
    top_down[cols_to_save].to_csv(down_csv)
    
    logger.info(f"Tabla firma UP 48h guardada en: {up_csv}")
    logger.info(f"Tabla firma DOWN 48h guardada en: {down_csv}")
    
    # Guardar listas de genes (útil para otros análisis)
    up_genes_txt = output_dir / f"firmas_DPN48h_up_{n_top}_genes.txt"
    down_genes_txt = output_dir / f"firmas_DPN48h_down_{n_top}_genes.txt"
    
    top_up.index.to_series().to_csv(up_genes_txt, index=False, header=False)
    top_down.index.to_series().to_csv(down_genes_txt, index=False, header=False)
    
    logger.info(f"Lista de genes UP 48h guardada en: {up_genes_txt}")
    logger.info(f"Lista de genes DOWN 48h guardada en: {down_genes_txt}")
    
    up_gene_list = top_up.index.tolist()
    down_gene_list = top_down.index.tolist()
    
    return up_gene_list, down_gene_list, top_up, top_down

# =============================================================================
# 3. CONSTRUIR RANKING DPN vs CONTROL 24h POR log2FoldChange
# =============================================================================

def build_dpn24_ranking(results_24h_path: Path):
    """
    Construye el ranking de genes usando SIEMPRE log2FoldChange.
    """
    if not results_24h_path.exists():
        raise FileNotFoundError(f"No se encontró el archivo de resultados 24h: {results_24h_path}")
    
    logger.info(f"Cargando resultados 24h desde: {results_24h_path}")
    res24 = pd.read_csv(results_24h_path, index_col=0)
    res24.index.name = "gene"
    
    if "log2FoldChange" not in res24.columns:
        raise ValueError("El fichero 24h debe contener la columna 'log2FoldChange' para usarla como ranking.")
    
    logger.info("Usando la columna 'log2FoldChange' como métrica de ranking para GSEA.")
    
    rnk = res24[["log2FoldChange"]].dropna().copy()
    rnk = rnk.rename(columns={"log2FoldChange": "score"})
    
    # Genes más sobreexpresados en DPN24h arriba del ranking
    rnk = rnk.sort_values("score", ascending=False)
    
    # DataFrame con columnas: gene, score (formato aceptado por gseapy.prerank)
    rnk_df = rnk.reset_index()
    
    logger.info(f"Ranking 24h construido con {rnk_df.shape[0]} genes.")
    
    return rnk_df, res24

# =============================================================================
# 4. GSEA PRERANKED + DOTPLOTS MODIFICADOS
# =============================================================================

def run_prerank_gsea_with_signatures(
    rnk_df: pd.DataFrame,
    up_genes: list,
    down_genes: list,
    res24_full: pd.DataFrame,
    outdir: Path,
):
    """
    Ejecuta GSEA preranked con gseapy.prerank usando como conjuntos de genes:
      - 'DPN48_UP'  : up_genes
      - 'DPN48_DOWN': down_genes
    
    y genera:
      - tabla de resultados res2d
      - enrichment plots
      - dotplot CONJUNTO (NES vs Term) para ambas firmas
      - dotplot de enriquecimiento individual de genes
    """
    
    gene_sets = {
        "DPN48_UP": up_genes,
        "DPN48_DOWN": down_genes,
    }
    
    logger.info("Iniciando GSEA preranked con firmas DPN48 sobre ranking DPN24 (log2FC)...")
    
    gs_res = gp.prerank(
        rnk=rnk_df,
        gene_sets=gene_sets,
        processes=4,
        permutation_num=N_PERMUTATIONS,
        min_size=MIN_SIZE,
        max_size=MAX_SIZE,
        seed=SEED,
        outdir=str(outdir),
        format="png",
    )
    
    # Tabla de resultados
    res_gsea = gs_res.res2d.copy()
    res_table_path = outdir / "GSEA_DPN48_signatures_on_DPN24_ranking_log2FC.csv"
    res_gsea.to_csv(res_table_path)
    logger.info(f"Tabla de resultados GSEA guardada en: {res_table_path}")
    
    # Enrichment plots para las dos firmas
    terms_to_plot = ["DPN48_UP", "DPN48_DOWN"]
    enrichment_plot_path = outdir / "GSEA_DPN48_signatures_on_DPN24_enrichment_plots_log2FC.png"
    logger.info(f"Generando enrichment plots para: {terms_to_plot}")
    gs_res.plot(terms=terms_to_plot, ofname=str(enrichment_plot_path))
    logger.info(f"Enrichment plots guardados en: {enrichment_plot_path}")
    
    # -------------------------------------------------------------------------
    # MODIFICACIÓN 1: Dotplot CONJUNTO para DPN48_UP y DPN48_DOWN
    # -------------------------------------------------------------------------
    
    # Aseguramos que tenemos las columnas con mayúsculas
    colmap = {}
    if "es" in res_gsea.columns and "ES" not in res_gsea.columns:
        colmap["es"] = "ES"
    if "nes" in res_gsea.columns and "NES" not in res_gsea.columns:
        colmap["nes"] = "NES"
    if "pval" in res_gsea.columns and "NOM p-val" not in res_gsea.columns:
        colmap["pval"] = "NOM p-val"
    if "fdr" in res_gsea.columns and "FDR q-val" not in res_gsea.columns:
        colmap["fdr"] = "FDR q-val"
    
    if colmap:
        res_gsea = res_gsea.rename(columns=colmap)
    
    # Si no existe 'Term', usamos el índice como nombre de conjunto
    if "Term" not in res_gsea.columns:
        res_gsea = res_gsea.reset_index().rename(columns={"index": "Term"})
    
    # Asegurar que FDR q-val es numérica
    res_gsea["FDR q-val"] = pd.to_numeric(res_gsea["FDR q-val"], errors="coerce")
    
    df_plot = res_gsea.copy()
    df_plot = df_plot.sort_values("NES")  # ordenar por NES
    
    logger.info("Generando dotplot CONJUNTO para DPN48_UP y DPN48_DOWN...")
    
    # Calcular -log10(FDR) de forma segura
    fdr_values = df_plot["FDR q-val"].values
    fdr_values = np.where(fdr_values <= 0, 1e-10, fdr_values)  # reemplazar <=0 por valor muy pequeño
    log_fdr = -np.log10(fdr_values + 1e-10)
    
    plt.figure(figsize=(8, 4))
    scatter = plt.scatter(
        df_plot["NES"],
        df_plot["Term"],
        s=log_fdr * 100,  # tamaño proporcional a -log10(FDR)
        c=log_fdr,
        cmap="viridis",
        edgecolors="black",
        linewidth=0.8,
        alpha=0.8,
    )
    
    plt.axvline(0, color="gray", linestyle="--", linewidth=1)
    plt.xlabel("NES", fontsize=12, fontweight='bold')
    plt.ylabel("Firma", fontsize=12, fontweight='bold')
    plt.title("Dotplot GSEA: Firmas DPN48 en ranking DPN24", fontsize=13, fontweight="bold")
    plt.colorbar(scatter, label="-log10(FDR q-val)")
    plt.tight_layout()
    
    combined_path = outdir / "dotplot_GSEA_DPN48_COMBINED_log2FC.png"
    plt.savefig(combined_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Dotplot conjunto guardado en: {combined_path}")
    
    # -------------------------------------------------------------------------
    # MODIFICACIÓN 2: Dotplot de Enriquecimiento Individual de Genes
    # -------------------------------------------------------------------------
    
    logger.info("Generando dotplot de enriquecimiento individual de genes...")
    
    # Preparar datos para mostrar genes individuales
    detailed_plot_data = []
    
    for idx, row in df_plot.iterrows():
        term = row['Term']
        nes = row['NES']
        fdr = row['FDR q-val']
        
        # Obtener genes que están en el leading edge
        leading_genes = row['Lead_genes'].split(';') if isinstance(row.get('Lead_genes'), str) else []
        
        # Para cada gen de la firma original, verificar si está en leading edge
        original_signature = gene_sets[term]
        
        # Obtener log2FC de cada gen desde res24
        for gene in original_signature[:20]:  # Mostrar top 20 para claridad
            in_leading = gene in leading_genes
            
            # Obtener log2FC del gen en DPN24h si está disponible
            gene_log2fc = res24_full.loc[gene, 'log2FoldChange'] if gene in res24_full.index else np.nan
            
            detailed_plot_data.append({
                'Firma': term.replace('DPN48_', ''),
                'Gen': gene,
                'NES': nes,
                'FDR_q_val': fdr,
                'In_Leading_Edge': in_leading,
                'Point_size': 120 if in_leading else 60,
                'log2FC_24h': gene_log2fc
            })
    
    detailed_df = pd.DataFrame(detailed_plot_data)
    
    # Crear figura para dotplot detallado
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Definir colores para cada firma
    colors = {'UP': '#e74c3c', 'DOWN': '#3498db'}
    
    # Scatter plot con colores según firma y tamaño según leading edge
    for firma in detailed_df['Firma'].unique():
        firma_data = detailed_df[detailed_df['Firma'] == firma]
        scatter = ax.scatter(
            firma_data['NES'],
            firma_data['Gen'],
            s=firma_data['Point_size'],
            alpha=0.7,
            label=f'Firma {firma}',
            color=colors.get(firma, 'gray'),
            edgecolors='black',
            linewidth=0.5
        )
    
    # Añadir línea vertical en NES=0
    ax.axvline(x=0, color='gray', linestyle='--', linewidth=1.5, alpha=0.5)
    
    # Etiquetas
    ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Genes (Top 20 de cada firma)', fontsize=14, fontweight='bold')
    ax.set_title('Enriquecimiento Individual de Genes de Firmas DPN 48h en DPN 24h',
                 fontsize=14, fontweight='bold', pad=20)
    
    # Leyenda mejorada
    handles, labels = ax.get_legend_handles_labels()
    from matplotlib.lines import Line2D
    
    # Añadir información sobre tamaño de puntos
    legend_elements = handles + [
        Line2D([0], [0], marker='o', color='w', label='Leading Edge',
               markerfacecolor='gray', markersize=12, markeredgecolor='black'),
        Line2D([0], [0], marker='o', color='w', label='No Leading Edge',
               markerfacecolor='gray', markersize=8, markeredgecolor='black')
    ]
    
    ax.legend(handles=legend_elements, loc='best', frameon=True, fontsize=10)
    
    plt.tight_layout()
    detailed_path = outdir / 'dotplot_enriquecimiento_individual_genes_DPN.png'
    plt.savefig(detailed_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.success(f"Dotplot de enriquecimiento individual guardado en: {detailed_path}")
    
    return res_gsea

# =============================================================================
# 5. PIPELINE COMPLETO
# =============================================================================

def main():
    # 1) Firmas DPN 48h
    up48_genes, down48_genes, top_up_df, top_down_df = generate_dpn48_signatures(
        results_48h_path=res_48h_path,
        n_top=N_TOP,
    )
    
    # 2) Ranking DPN 24h (por log2FoldChange)
    rnk_24_df, res24_full = build_dpn24_ranking(results_24h_path=res_24h_path)
    
    # 3) GSEA preranked + dotplots modificados
    res_gsea = run_prerank_gsea_with_signatures(
        rnk_df=rnk_24_df,
        up_genes=up48_genes,
        down_genes=down48_genes,
        res24_full=res24_full,
        outdir=output_dir,
    )
    
    # 4) Resumen final
    print("\n=== RESUMEN GSEA: Firmas DPN48 en ranking DPN24 (log2FC) ===")
    print(res_gsea[["ES", "NES", "NOM p-val", "FDR q-val"]])
    print(f"\nAnálisis completado. Todos los resultados están en:\n {output_dir}")

if __name__ == "__main__":
    main()
