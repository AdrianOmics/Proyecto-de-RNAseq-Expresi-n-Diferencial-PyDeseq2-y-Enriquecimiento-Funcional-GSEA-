#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Análisis de expresión diferencial (24h y 48h) con PyDESeq2 usando AnnData

Datos de entrada:
- Tabla de cuentas crudas (genes x muestras):
    ;GSM913873;GSM913874;GSM913875;...
    ENSG00000000003;673;916;391;...
    ...

- Tabla de metadatos:
    ;patient;agent;time
    GSM913873;1;Control;24h
    GSM913874;1;Control;48h
    GSM913875;1;DPN;24h
    ...

Objetivo:
- Crear AnnData (muestras x genes) + metadatos (obs)
- Para cada tiempo (24h, 48h):
    - Diseño pareado: ~ patient + agent
    - Contrastes: DPN vs Control, OHT vs Control
    - Tablas completas y filtradas de DE
    - Volcano plots
"""

# =============================================================================
# 0. IMPORTAR PAQUETES
# =============================================================================

from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

sns.set_theme(style="whitegrid", palette="deep")


# =============================================================================
# 1. CONFIGURACIÓN
# =============================================================================

BASE_DIR = Path("/home/adrian/Documentos/Transcriptomica/Trabajo_final_RNAseq/Deseq2")  # Ajusta si lo necesitas

# Rutas de entrada (AJUSTA LOS NOMBRES A TUS FICHEROS REALES)
COUNTS_PATH = BASE_DIR / "rawcounts.csv"
METADATA_PATH = BASE_DIR / "metadata.csv"

COUNTS_SEP = ";"    
METADATA_SEP = ";"   

# Parámetros DESeq2
N_CPUS = 8
DESIGN_FORMULA = "~ patient + agent"  # diseño pareado por paciente
CONTRAST_VAR = "agent"
REF_LEVEL = "Control"
TEST_LEVELS = ["DPN", "OHT"]         # comparaciones de interés

# Filtros de significancia
FDR_THRESHOLD = 0.05
LFC_THRESHOLD = 1.0
MIN_COUNTS_TOTAL = 10                # filtro de genes con muy baja expresión

OUTPUT_DIR = BASE_DIR / "pydeseq2_parathyroid_DEG"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

inference = DefaultInference(n_cpus=N_CPUS)

# =============================================================================
# 2. CARGA Y CREACIÓN DEL OBJETO AnnData COMPLETO
# =============================================================================

print("\n" + "=" * 80)
print("BLOQUE 1: CARGA DE DATOS Y CREACIÓN DE AnnData")
print("=" * 80)

# -------------------------------------------------------------------------
# 2.1. Cargar matriz de cuentas
# -------------------------------------------------------------------------
print("\n[1.1] Cargando matriz de cuentas crudas...")

# Estructura esperada:
#   primera columna: gene_id (ENSG...)
#   columnas siguientes: GSM913873, GSM913874, ...
counts_raw = pd.read_csv(COUNTS_PATH, sep=COUNTS_SEP, index_col=0)

print(f"Shape cuentas (genes x muestras): {counts_raw.shape}")
print("Primeras columnas (muestras):", counts_raw.columns[:6].tolist())
print("Primeros genes:", counts_raw.index[:5].tolist())

# -------------------------------------------------------------------------
# 2.2. Cargar metadatos
# -------------------------------------------------------------------------
print("\n[1.2] Cargando metadatos...")

metadata = pd.read_csv(METADATA_PATH, sep=METADATA_SEP, index_col=0)

print(f"Shape metadatos (muestras x variables): {metadata.shape}")
print("Columnas de metadatos:", metadata.columns.tolist())
print("Ejemplo de metadatos:")
print(metadata.head())

# -------------------------------------------------------------------------
# 2.3. Alinear muestras entre cuentas y metadatos
# -------------------------------------------------------------------------
print("\n[1.3] Alineando IDs de muestras entre cuentas y metadatos...")

samples_in_both = metadata.index.intersection(counts_raw.columns)

if len(samples_in_both) == 0:
    raise ValueError(
        "No hay solapamiento entre IDs de muestras en cuentas y metadatos.\n"
        "Revisa que las columnas de la matriz de cuentas sean los GSM y "
        "que el índice de metadatos sean los mismos IDs."
    )

# Reordenar y filtrar para quedarnos con el mismo conjunto de muestras
metadata = metadata.loc[samples_in_both].copy()
counts_raw = counts_raw[samples_in_both].copy()

print(f"Nº de muestras comunes: {len(samples_in_both)}")
print("IDs comunes (primeros):", samples_in_both[:6].tolist())

# -------------------------------------------------------------------------
# 2.4. Construir AnnData completo (todas las muestras, ambos tiempos)
# -------------------------------------------------------------------------
print("\n[1.4] Creando objeto AnnData global (todas las muestras)...")

# AnnData espera: filas = muestras, columnas = genes
counts_matrix = counts_raw.T  # ahora: filas = muestras, columnas = genes

# Preparar obs (metadatos de muestras)
# Aseguramos tipos adecuados
metadata["patient"] = metadata["patient"].astype(str)
metadata["agent"] = metadata["agent"].astype(str)
metadata["time"] = metadata["time"].astype(str)

# Metadatos de genes (var)
var_df = pd.DataFrame(
    {"gene_id": counts_raw.index.tolist()},
    index=counts_raw.index.tolist(),
)

# Crear AnnData
adata_full = ad.AnnData(
    X=counts_matrix.values.astype("int32"),
    obs=metadata,
    var=var_df,
)

print("Resumen AnnData global:")
print(f"  Nº muestras (obs): {adata_full.n_obs}")
print(f"  Nº genes (var):    {adata_full.n_vars}")
print(f"  Tipo de X:         {adata_full.X.dtype}")

print("\nPrimeras filas de obs (metadatos):")
print(adata_full.obs.head())

print("\nPrimeras filas de var (genes):")
print(adata_full.var.head())

print("\n✅ BLOQUE 1 COMPLETADO: AnnData global creado correctamente.")


# =============================================================================
# 3. FUNCIÓN AUXILIAR: subsetting por tiempo, DESeq2 y Volcano plot
# =============================================================================

def run_deseq_for_time(
    adata_full,
    time_label,
    design_formula=DESIGN_FORMULA,
    contrast_var=CONTRAST_VAR,
    ref_level=REF_LEVEL,
    test_levels=TEST_LEVELS,
    fdr_threshold=FDR_THRESHOLD,
    lfc_threshold=LFC_THRESHOLD,
    min_counts_total=MIN_COUNTS_TOTAL,
    output_dir=OUTPUT_DIR,
):
    """
    Toma el AnnData completo y:
      - filtra por tiempo (24h o 48h)
      - aplica filtro de baja expresión
      - ajusta DESeq2 con diseño ~ patient + agent
      - ejecuta contrastes (DPN vs Control, OHT vs Control)
      - genera tablas y Volcano plots

    Devuelve:
      - dict con results_df y de_genes por contraste
    """
    print("\n" + "=" * 80)
    print(f"ANÁLISIS PARA TIEMPO: {time_label}")
    print("=" * 80)

    # ---------------------------------------------------------------------
    # 3.1. Subconjunto por tiempo
    # ---------------------------------------------------------------------
    print(f"\n[2.1 - {time_label}] Filtrando muestras por tiempo...")

    mask_time = adata_full.obs["time"] == time_label
    if mask_time.sum() == 0:
        raise ValueError(f"No se encuentran muestras con time == '{time_label}'")

    adata_t = adata_full[mask_time, :].copy()

    print(f"Muestras en {time_label}: {adata_t.n_obs}")
    print("Distribución de 'agent' en este tiempo:")
    print(adata_t.obs["agent"].value_counts())

    # ---------------------------------------------------------------------
    # 3.2. Preparar factores y filtro de genes de baja expresión
    # ---------------------------------------------------------------------
    print(f"\n[2.2 - {time_label}] Preparando factores y filtrando genes...")

    # patient como string/categórica
    adata_t.obs["patient"] = adata_t.obs["patient"].astype(str)

    # agent como categórica con Control como referencia
    adata_t.obs["agent"] = pd.Categorical(
        adata_t.obs["agent"],
        categories=["Control", "DPN", "OHT"],
    )

    # Filtro de genes por cuentas totales
    gene_sums = np.asarray(adata_t.X.sum(axis=0)).ravel()
    genes_to_keep = gene_sums >= min_counts_total

    n_genes_before = adata_t.n_vars
    adata_t = adata_t[:, genes_to_keep].copy()
    n_genes_after = adata_t.n_vars

    print(f"Genes antes del filtro: {n_genes_before}")
    print(f"Genes después del filtro (suma >= {min_counts_total}): {n_genes_after}")

    # ---------------------------------------------------------------------
    # 3.3. Crear DeseqDataSet y ejecutar DESeq2
    # ---------------------------------------------------------------------
    print(f"\n[2.3 - {time_label}] Creando DeseqDataSet y ejecutando DESeq2...")

    dds = DeseqDataSet(
        adata=adata_t,
        design=design_formula,  # "~ patient + agent"
        refit_cooks=True,
        inference=inference,
    )

    dds.deseq2()

    print("DESeq2 completado.")
    print("Dimensiones internas del dds:", dds.shape)
    print("Matriz de diseño (columnas):", dds.obsm["design_matrix"].columns.tolist())

    # ---------------------------------------------------------------------
    # 3.4. Función interna para un contraste concreto
    # ---------------------------------------------------------------------
    def run_contrast(
        dds,
        test_level,
        ref_level,
        time_label,
        contrast_var="agent",
        fdr_threshold=0.05,
        lfc_threshold=1.0,
        output_dir=output_dir,
    ):
        """
        Ejecuta un contraste (p.ej. DPN vs Control) para un tiempo concreto.
        """

        print("\n" + "-" * 70)
        print(f"[{time_label}] Contraste: {test_level} vs {ref_level} (variable: {contrast_var})")
        print("-" * 70)

        ds = DeseqStats(
            dds,
            contrast=[contrast_var, test_level, ref_level],
            inference=inference,
        )

        ds.summary()

        results_df = ds.results_df.copy()
        results_df.index.name = "gene"

        # Intentar LFC shrinkage
        coeff_name = f"{contrast_var}[T.{test_level}]"
        design_cols = list(ds.design_matrix.columns)

        if coeff_name in design_cols:
            print(f"Aplicando LFC shrinkage sobre coeficiente: {coeff_name}")
            ds.lfc_shrink(coeff=coeff_name)
            results_df = ds.results_df.copy()
        else:
            print(f"AVISO: coeficiente '{coeff_name}' no encontrado en la matriz de diseño.")
            print("       Se continúa sin shrinkage de LFC.")

        # Filtrado de genes significativos
        sig_mask = (
            (results_df["padj"] < fdr_threshold)
            & (results_df["log2FoldChange"].abs() >= lfc_threshold)
        )

        de_genes = results_df[sig_mask].copy()
        up_genes = de_genes[de_genes["log2FoldChange"] > 0].copy()
        down_genes = de_genes[de_genes["log2FoldChange"] < 0].copy()

        print("\nResumen de genes diferencialmente expresados:")
        print(f"  Genes analizados: {len(results_df)}")
        print(f"  Genes DE (padj < {fdr_threshold}, |log2FC| >= {lfc_threshold}): {len(de_genes)}")
        print(f"    - Up-regulados en {test_level}:   {len(up_genes)}")
        print(f"    - Down-regulados en {test_level}: {len(down_genes)}")

        # Guardar tablas
        tag = f"{test_level}_vs_{ref_level}_{time_label}"

        full_path = output_dir / f"deseq2_results_{tag}.csv"
        de_path = output_dir / f"deseq2_DEgenes_{tag}_FDR{fdr_threshold}_LFC{lfc_threshold}.csv"

        results_df.to_csv(full_path)
        de_genes.to_csv(de_path)

        print("\nFicheros guardados:")
        print(f"  - Tabla completa: {full_path}")
        print(f"  - Tabla filtrada: {de_path}")

        # Volcano plot
        print("\nGenerando Volcano plot...")

        volcano_df = results_df.replace([np.inf, -np.inf], np.nan).dropna(
            subset=["padj", "log2FoldChange"]
        )

        # Evitar infinito si padj=0
        if (volcano_df["padj"] == 0).any():
            min_nonzero = volcano_df["padj"][volcano_df["padj"] > 0].min()
            volcano_df.loc[volcano_df["padj"] == 0, "padj"] = min_nonzero / 10

        x = volcano_df["log2FoldChange"]
        y = -np.log10(volcano_df["padj"])

        sig = (
            (volcano_df["padj"] < fdr_threshold)
            & (volcano_df["log2FoldChange"].abs() >= lfc_threshold)
        )

        plt.figure(figsize=(8, 6))
        plt.scatter(
            x[~sig], y[~sig],
            c="lightgray", s=20, alpha=0.6,
            label="No significativo",
        )
        plt.scatter(
            x[sig], y[sig],
            c="red", s=25, alpha=0.8,
            label=f"DE (padj<{fdr_threshold}, |log2FC|>={lfc_threshold})",
        )

        # Líneas de referencia
        plt.axhline(-np.log10(fdr_threshold), color="blue", ls="--", lw=1, alpha=0.7,
                    label=f"FDR={fdr_threshold}")
        plt.axvline(lfc_threshold, color="green", ls="--", lw=1, alpha=0.7,
                    label=f"log2FC=±{lfc_threshold}")
        plt.axvline(-lfc_threshold, color="green", ls="--", lw=1, alpha=0.7)

        plt.xlabel("log2 Fold Change")
        plt.ylabel("-log10(p-adj)")
        plt.title(f"Volcano plot ({time_label}): {test_level} vs {ref_level}")
        plt.legend()
        plt.tight_layout()

        volcano_path = output_dir / f"volcano_{tag}.png"
        plt.savefig(volcano_path, dpi=300)
        plt.close()

        print(f"Volcano plot guardado en: {volcano_path}")

        return results_df, de_genes

    # ---------------------------------------------------------------------
    # 3.5. Ejecutar todos los contrastes para este tiempo
    # ---------------------------------------------------------------------
    results_time = {}
    de_time = {}

    for test in test_levels:
        res_df, de_df = run_contrast(
            dds=dds,
            test_level=test,
            ref_level=ref_level,
            time_label=time_label,
            contrast_var=contrast_var,
            fdr_threshold=fdr_threshold,
            lfc_threshold=lfc_threshold,
            output_dir=output_dir,
        )
        key = f"{test}_vs_{ref_level}_{time_label}"
        results_time[key] = res_df
        de_time[key] = de_df

    print(f"\nRESUMEN GLOBAL ({time_label}):")
    for comp, de_df in de_time.items():
        print(f"  {comp}: {len(de_df)} genes DE (padj<{fdr_threshold}, |log2FC|>={lfc_threshold})")

    return results_time, de_time


# =============================================================================
# 4. EJECUTAR ANÁLISIS PARA 24h Y 48h
# =============================================================================

if __name__ == "__main__":
    # 24h
    results_24h, de_24h = run_deseq_for_time(
        adata_full=adata_full,
        time_label="24h",
        design_formula=DESIGN_FORMULA,
        contrast_var=CONTRAST_VAR,
        ref_level=REF_LEVEL,
        test_levels=TEST_LEVELS,
        fdr_threshold=FDR_THRESHOLD,
        lfc_threshold=LFC_THRESHOLD,
        min_counts_total=MIN_COUNTS_TOTAL,
        output_dir=OUTPUT_DIR,
    )

    # 48h
    results_48h, de_48h = run_deseq_for_time(
        adata_full=adata_full,
        time_label="48h",
        design_formula=DESIGN_FORMULA,
        contrast_var=CONTRAST_VAR,
        ref_level=REF_LEVEL,
        test_levels=TEST_LEVELS,
        fdr_threshold=FDR_THRESHOLD,
        lfc_threshold=LFC_THRESHOLD,
        min_counts_total=MIN_COUNTS_TOTAL,
        output_dir=OUTPUT_DIR,
    )

    print("\n=== ANÁLISIS COMPLETO FINALIZADO ===")
    print(f"Resultados (tablas + volcano plots) en: {OUTPUT_DIR.resolve()}")

