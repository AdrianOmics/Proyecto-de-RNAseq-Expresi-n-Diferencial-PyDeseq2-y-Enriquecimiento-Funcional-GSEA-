#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""

Created on Mon Jan 26 11:25:58 2026

Work Flow RNAseq (BASH Y PYTHON)
Parte 1: Análisis desde Ficheros FASTQ (6 puntos)

@author: Adrián Aguilera Martín.

"""

### [ANÁLISIS] SECCIÓN 1: PRE-PROCESAMIENTO Y CONTROL DE CALIDAD (BASH)
# Esta sección está escrita en sintaxis Bash. En un entorno de producción, 
# se ejecutaría en una terminal o script .sh separado.

# Comienzo realizando un FASTQC de las muestras.
# El comando 'fastqc' genera un informe HTML para visualizar la calidad base por base, 
# contenido de GC, niveles de duplicación y presencia de adaptadores.
# fastqc *.fastq

# Tras visualizar la calidad de las muestras, realizo trimming con FASTP.

"""

--qualified_quality_phred 28: Una base se considera "buena" si tiene Phred ≥ 28. Como tu calidad es > 30, es permisivo.

--unqualified_percent_limit 10: Descarta reads si > 10% de sus bases están por debajo de Q20. Esto elimina reads muy degradados sin ser agresivo.

--length_required 40: Descarta reads < 40 pb después de trimming (tus reads son 50 pb, así que es suave).

--cut_front y --cut_tail: Elimina bases de baja calidad desde los extremos 5' y 3' del read usando una ventana deslizante.

--cut_window_size 4 y --cut_mean_quality 20: Ventana de 4 bases; si la calidad media en la ventana cae por debajo de Q20, recorta. Esto es suave (para comparación: trimming agresivo usaría Q25-Q28).

--html y --json: Genera reportes de calidad antes/después.

-j 4: Usa 4 threads de procesamiento (acelera el análisis).

"""

# Array con los nombres de las muestras
# Se definen los archivos de entrada brutos (raw data).
samples=("SRR3480371.chr21.fastq" "SRR3480372.chr21.fastq" "SRR3480373.chr21.fastq" "SRR3480374.chr21.fastq")

# Crear directorios de salida
# -p asegura que no de error si la carpeta ya existe.
mkdir -p trimmed_results
mkdir -p reports_fastp

# Bucle para procesar cada muestra
# Itera sobre cada elemento del array 'samples'.
for sample in "${samples[@]}"; do
# Extraer nombre sin extensión
# ${sample%.fastq} es una expansión de parámetros de Bash que elimina el sufijo ".fastq".
    sample_name="${sample%.fastq}"
    echo "Trimming: $sample_name"
# Comando fastp para single-end con parámetros específicos
# NOTA TÉCNICA: Al no especificar --detect_adapter_for_pe (pues es single-end), 
# fastp intentará detectar adaptadores automáticamente o por solapamiento si fuera PE. 
# Para SE, busca adaptadores comunes de Illumina si no se especifican.
    fastp -i "$sample" \
        -o "trimmed_results/${sample_name}.trimmed.fastq" \
        --qualified_quality_phred 28 \
        --unqualified_percent_limit 10 \
        --length_required 40 \
        --cut_front \
        --cut_tail \
        --cut_window_size 4 \
        --cut_mean_quality 20 \
        --html "reports_fastp/${sample_name}.html" \
        --json "reports_fastp/${sample_name}.json" \
        -j 4

    echo "Completado: $sample_name"
done

echo "Todas las muestras han sido procesadas"

### [ANÁLISIS] SECCIÓN 2: CUANTIFICACIÓN CON SALMON (BASH)
# Salmon es un "pseudoalineador". No alinea base a base contra el genoma completo,
# sino que mapea lecturas a transcritos para estimar abundancia. Es mucho más rápido que STAR/HISAT2.

# Crear índice para Salmon (Bash)

# Parámetros:

# - -t chr21.fa: FASTA de transcritos. 
#   IMPORTANTE: Debe ser un FASTA de *transcriptoma* (cDNA), no de genoma completo (DNA), 
#   a menos que uses modo "decoy" (recomendado para evitar mapeos espurios).

# - -i salmon_index: directorio del índice de Salmon.

# - -k 31: tamaño de k-mer. 
#   NOTA: 31 es el valor por defecto. Para reads muy cortos (<50pb), un k menor (ej. 21-25) 
#   puede mejorar la sensibilidad, pero para >=50pb, k=31 es óptimo para especificidad.

# - --gencode: indica formato Gencode de los IDs.
#   Esto ayuda a Salmon a manejar los IDs compuestos tipo "ENSTXXX.Y|GENEXXX|..." limpiando los nombres.

# - -p 8: hilos.

#-------------------------------------------------------------------------------------------------
# Antes de hacer el índice con salmon, debo obtener el transcriptoma para el chr21.
mkdir -p ref
cd ref

# cDNA de Homo sapiens GRCh38 (ajusta el "release" si hace falta)
wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# Uso awk para quedarme únicamente con las entradas cuyo header contenga chromosome:GRCh38:21: (patrón típico de Ensembl).
zcat Homo_sapiens.GRCh38.cdna.all.fa.gz \
  | awk 'BEGIN{keep=0}
         /^>/{
           keep=0
           if ($0 ~ /chromosome:GRCh38:21:/) keep=1
         }
         { if (keep) print }' \
  | gzip > Homo_sapiens.GRCh38.cdna.chr21.fa.gz

# Cambio el nombre del archivo final
mv Homo_sapiens.GRCh38.cdna.chr21.fa.gz chr21.fa.gz

# Hacer un index para salmon (si el .fa esta en el formato adecuado) (Bash).

if [ ! -f salmon_index/hash.bin ]; then
    echo "Creando índice de Salmon..."
    salmon index \
        -t ./chr21.fa.gz \
        -i salmon_index \
        -k 31 \
        --gencode \
        -p 8
fi

#Pseudoalineamiento con Salmon (Bash)

# ---------------------------------------------------------------------

# Propósito:

# - Cuantificar la expresión a nivel de transcrito usando Salmon

# en modo quasi-mapping (single-end).

# - Generar quant.sf (y luego comprimirlo) para cada muestra.

# -i salmon_index: índice del transcriptoma.
# * -l U: librería single-end, unstranded. 
#   NOTA: En el 
# * comando abajo usas '-l A' (automático), lo cual es mejor que 'U' si desconoces la preparación.
#   Salmon inferirá si es stranded (forward/reverse) o unstranded basándose en los primeros miles de reads.

# * -r FASTQ: archivo FASTQ de entrada.

# * -p 8: hilos.

# * --validateMappings: 
#   CRÍTICO. Activa una fase extra de validación geométrica del alineamiento. 
#   Mejora la precisión (reduce falsos positivos) a costa de un poco de tiempo. Muy recomendado.

# * --seqBias: Corrige el sesgo de secuencia aleatoria (hexámeros) en el inicio de los reads (común en RNA-seq).
# * --gcBias: Corrige el sesgo por contenido de GC, que puede afectar la eficiencia de secuenciación de ciertos fragmentos.

# * -o salmon_results/SRR...: carpeta de salida.

# - Verifica que quant.sf se genere y lo comprime con gzip.

#-------------------------------------------------------------------------------------------------

# Hacer un pseudoalineamiento con salmon. De esta manera se generará ya el archivo comprimido para usarlo posteriormente con Pymportx (Bash).

mkdir -p salmon_results
WD="/home/adrian/Documentos/Transcriptomica/Trabajo_Final_RNAseq"
TRIM_DIR="${WD}/Seq_Raw/trimmed_results"

for FASTQ in SRR3480371 SRR3480372 SRR3480373 SRR3480374; do
    echo "Procesando: $FASTQ"

    salmon quant \
        -i ref/salmon_index \
        -l U \
        -r "${TRIM_DIR}/${FASTQ}.chr21.trimmed.fastq" \
        -p 8 \
        --validateMappings \
        --seqBias \
        --gcBias \
        -o "salmon_results/${FASTQ}_salmon_quant_sample"
# Verificar si se generó el archivo antes de comprimir
# Salmon genera un archivo 'quant.sf' que contiene las abundancias (TPM, NumReads).
# Muchos parsers esperan 'quant.sf', pero comprimirlo ahorra espacio. Pymportx suele manejar .gz transparentemente.

    if [ -f "salmon_results/${FASTQ}_salmon_quant_sample/quant.sf" ]; then
        echo "Éxito: quant.sf generado para $FASTQ"
# gzip -f fuerza la compresión y reemplaza el original
        gzip -f "salmon_results/${FASTQ}_salmon_quant_sample/quant.sf"
    else
        echo "ERROR: No se generó quant.sf para $FASTQ"
        echo "Revisa el log: salmon_results/${FASTQ}_salmon.log"
    fi
done

#--------------------------------------------------------------------------------------------------------------------------------

### [ANÁLISIS] SECCIÓN 3: IMPORTACIÓN Y ANÁLISIS EN PYTHON
# A partir de aquí comienza el código Python puro.
# El objetivo es colapsar la información de transcritos (isoformas) a nivel de gen para análisis de expresión diferencial.

#Importación con Pymportx y creación de AnnData (Python)

# Propósito:

# - Leer los resultados de Salmon a nivel de transcrito.

# - Resumirlos a nivel de gen usando tx2gene.

# - Crear un objeto AnnData con las matrices de expresión.

# - Exportar conteos a h5ad y CSV para análisis posteriores.

# Lógica:

# - Define WORKING_DIR y RESULTS_DIR.

# - Usa salmon.read_salmon() de pymportx:

# * folders: rutas a carpetas de resultados de Salmon (4 muestras).

# * tx_out=False: salida a nivel de gen (no transcritos). 
#   IMPORTANTE: Para análisis funcional y DGE robusto, solemos trabajar a nivel de gen.

# * tx2gene: CSV con mapeo transcrito-gen.

# * countsFromAbundance='no': usar conteos sin reescalar por abundancia.
#   EXPLICACIÓN: 'no' significa que toma la columna 'NumReads' estimada por Salmon.
#   Otras opciones como 'lengthScaledTPM' ajustarían los conteos por la longitud del transcrito,
#   lo cual es útil si comparas genes entre sí, pero para comparar Muestras (DGE), los conteos brutos son preferibles.

# - Guarda AnnData en gene_counts.h5ad.

# - Extrae la matriz de conteos (adata.X) a un DataFrame.

# - Calcula estadísticas básicas (top genes, conteos por muestra).

# - Exporta counts_matrix.csv, tpm_matrix.csv (si existe 'abundance') y

# sample_stats.csv.

# Parámetros clave:

# - folders: lista de carpetas "SRR348037i_salmon_quant_sample".

# - tx2gene: gencode_tx2gene.csv (preparado previamente).

# - countsFromAbundance='no': configuración de tximport-like.

#-------------------------------------------------------------------------------------------------

# Realizamos el pymportx en base a los quant comprimidos (Python).

# ANTES DE HACER PYMPORTX, ES NECESARIO OBTENER EL ARCHIVO: genecode_tx2gene.csv ADAPTADO AL CHR21. A CONTINUACIÓN INCLUYO EL CÓDIGO COMENTADO.

from pathlib import Path

import pandas as pd

import gzip

import re

# %% 2. Organizar directorios
# Uso de pathlib para manejo de rutas agnóstico del sistema operativo (Windows/Linux/Mac).

WORKING_DIR = Path.home() / "Documentos" / "Transcriptomica" / "Trabajo_final_RNAseq" / "tx2gene"

# Ruta al archivo GTF de GENCODE (descargado)
# El archivo GTF (General Transfer Format) contiene la anotación genómica necesaria para saber qué transcrito pertenece a qué gen.

gtf_file = Path.home() / "Documentos" / "Transcriptomica" / "Trabajo_final_RNAseq" / "ref" / "GRCh38.gencode.v38.annotation.for.chr21.gtf"

# %% 3. Extración de campos par generar el archivo tx2gene

def extract_tx2gene_from_gencode(gtf_path, out_put_path=None):

    """Extrae mapeo transcript_id -> gene_id de archivo GTF de GENCODE"""
    
    print(f"Extrayendo tx2gene de GENCODE: {gtf_path}")
    
    tx2gene_data = set() # Usar set para evitar duplicados desde el principio (optimización de memoria)
    
    # Contadores
    
    line_count = 0
    
    transcript_lines = 0
    
    # Apertura inteligente del archivo (gzip o texto plano según extensión)
    with gzip.open(gtf_path, 'rt') if str(gtf_path).endswith('.gz') else open(gtf_path, 'r') as f:
    
        for line in f:
        
            line_count += 1
            
            if line_count % 1000000 == 0:
            
                print(f" Procesadas {line_count:,} líneas...")
            
            # Saltar comentarios del GTF
            if line.startswith('#'):
            
                continue
            
            parts = line.strip().split('\t')
            
            if len(parts) < 9:
            
                continue
            
            # Solo líneas de tipo 'transcript' (más eficiente)
            # El GTF tiene líneas para 'gene', 'exon', 'CDS', etc. Solo nos interesan las que definen un transcrito 
            # para vincularlo a su gen padre.
            if parts[2] == 'transcript':
            
                transcript_lines += 1
                
                # La columna 9 (índice 8) contiene los atributos: gene_id "..."; transcript_id "...";
                info = parts[8]
                
                # Parsear atributos con regex (más robusto para GENCODE)
                
                # GENCODE tiene atributos como: gene_id "ENSG00000223972.5"; gene_name "DDX11L1";
                
                attributes = {}
                
                # Buscar todos los pares clave-valor
                # Regex explica: (\w+) captura la clave, \s+ el espacio, "([^"]+)" captura el valor entre comillas.
                pattern = r'(\w+)\s+"([^"]+)"'
                
                matches = re.findall(pattern, info)
                
                for key, value in matches:
                
                    attributes[key] = value
                
                # Extraer transcript_id y gene_id
                
                transcript_id = attributes.get('transcript_id')
                
                gene_id = attributes.get('gene_id')
                
                if transcript_id and gene_id:
                    # NOTA: Los IDs de GENCODE suelen tener versión (ej. ENSG000001.5). 
                    # Si Salmon usó un índice sin versiones, esto podría causar problemas de mapeo.
                    # Aquí se extraen tal cual vienen en el GTF.
                    tx2gene_data.add((transcript_id, gene_id))
    
    # Crear DataFrame
    
    df = pd.DataFrame(list(tx2gene_data), columns=['TXNAME', 'GENEID'])
    
    print("\n📊 Estadísticas:")
    
    print(f" Líneas totales procesadas: {line_count:,}")
    
    print(f" Líneas de transcript: {transcript_lines:,}")
    
    print(f" Mapeos únicos encontrados: {len(df):,}")
    
    # Manejar la ruta de salida
    
    if out_put_path is None:
    
        out_put_path = WORKING_DIR
    
    else:
    
        out_put_path = Path(out_put_path)
    
    # Crear directorio si no existe
    
    out_put_path.mkdir(parents=True, exist_ok=True)
    
    # Extraer nombre base del archivo GTF
    
    gtf_name = Path(gtf_path).stem # Quita la extensión
    
    if gtf_name.endswith('.gtf'):
    
        gtf_name = Path(gtf_name).stem # Quita .gtf también si existe
    
    # Definir nombre del archivo de salida
    
    output_file = out_put_path / f"{gtf_name}_tx2gene.csv"
    
    # Guardar CSV
    
    df.to_csv(output_file, index=False)
    
    print(f"💾 Archivo CSV guardado como: {output_file}")
    
    return df


# Uso (manteniendo tu código original)

tx2_gene = extract_tx2gene_from_gencode(gtf_file)

# El DataFrame ya está guardado automáticamente en CSV

# Puedes usar tx2_gene normalmente en tu análisis

print(tx2_gene.head())

# %% 4. Después de crear tx2gene, verifica con tus datos

# Después de crear tx2gene, verifica con tus datos

def verify_tx2gene_with_data(tx2gene_path, sample_quant_path):

    """Verifica cuántos transcritos de tus datos están en tx2gene"""
    
    print("\n🔍 Verificando tx2gene con tus datos...")
    
    # Leer tx2gene
    
    tx2gene = pd.read_csv(tx2gene_path)
    
    # Leer un archivo quant.sf de muestra
    # Se usa compression='gzip' porque en el paso de Bash comprimimos los quant.sf
    df_quant = pd.read_csv(sample_quant_path, sep='\t', compression='gzip', nrows=1000)
    
    # IDs de transcritos en tus datos
    
    data_transcripts = set(df_quant['Name'])
    
    # IDs en tx2gene
    
    tx2gene_transcripts = set(tx2gene['TXNAME'])
    
    # Coincidencias
    # Intersección de conjuntos para ver si los IDs de Salmon coinciden con los del GTF.
    matches = data_transcripts.intersection(tx2gene_transcripts)
    
    print(f" Transcritos en datos (muestra): {len(data_transcripts)}")
    
    print(f" Transcritos en tx2gene: {len(tx2gene_transcripts)}")
    
    print(f" Coincidencias encontradas: {len(matches)}")
    
    print(f" Porcentaje de coincidencia: {len(matches)/len(data_transcripts)*100:.1f}%")
    
    if len(matches) == 0:
    
        print("\n⚠️ ¡CERO coincidencias! Mostrando ejemplos...")
        
        print("\nEjemplos en tus datos (primeros 5):")
        
        for tx in list(data_transcripts)[:5]:
        
            print(f" • {tx}")
        
        print("\nEjemplos en tx2gene (primeros 5):")
        
        for tx in list(tx2gene_transcripts)[:5]:
        
            print(f" • {tx}")
        
        # Verificar si es problema de versión
        # A menudo Salmon indexa "ENST00001" y el GTF tiene "ENST00001.2". 
        # Esta lógica detecta si quitar el ".2" solucionaría el problema.
        print("\n🔧 Probando sin números de versión...")
        
        data_no_version = {tx.split('.')[0] for tx in data_transcripts}
        
        tx2gene_no_version = {tx.split('.')[0] for tx in tx2gene_transcripts}
        
        matches_no_version = data_no_version.intersection(tx2gene_no_version)
        
        print(f" Coincidencias (sin versión): {len(matches_no_version)}")
        
        if len(matches_no_version) > 0:
        
            print(" ✅ ¡Mejor! Usa la versión sin números de versión")
    
    return len(matches) > 0

# Ejecutar verificación

sample_quant = WORKING_DIR.parent / "salmon_results" / "SRR3480371_salmon_quant_sample" / "quant.sf.gz"

tx2gene_file = WORKING_DIR / "GRCh38.gencode.v38.annotation.for.chr21_tx2gene.csv"

if sample_quant.exists() and tx2gene_file.exists():

    verification_ok = verify_tx2gene_with_data(tx2gene_file, sample_quant)
    
    if verification_ok:
    
        print("\n✅ tx2gene parece compatible con tus datos")
    
    else:
    
        print("\n❌ tx2gene NO es compatible con tus datos")

else:

    print("No se pueden encontrar los archivos para verificación")

# AHORA PODEMOS EJECUTAR ESTE SCRIPT DE PYMPORTX

# 1. Importar Paquetes

from pathlib import Path
import pandas as pd
import anndata as ad
from pymportx import salmon

# %% 2. Definir las rutas

WORKING_DIR = Path.home() / "Documentos" /"Transcriptomica" / "Trabajo_final_RNAseq"

RESULTS_DIR = WORKING_DIR / "Pymportx_results"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# %% 3. Importar los datos a Pymportx y generamos el objeto AnnData...

# Genero el archivo genecode_tx2gene.csv utilizando el script tx2gene.py utilizando como referencia el archivo: GRCh38.gencode.v38.annotation.for.chr21.gtf
# salmon.read_salmon es la función clave. Lee los archivos quant.sf, usa tx2gene para sumarizar a nivel de gen
# y devuelve un objeto AnnData (formato estándar para single-cell y transcriptómica en Python).
adata = salmon.read_salmon(
    folders=[str(WORKING_DIR / "salmon_results" / f"SRR348037{i}_salmon_quant_sample") for i in range(1, 5)],
    tx_out=False, # False = colapsar a genes. True = mantener a nivel de transcrito.
    tx2gene=str(WORKING_DIR / "tx2gene" / "genecode_tx2gene.csv"),
    countsFromAbundance='no' # 'no' usa la columna NumReads directamente. Esencial para DESeq2.
)

# 5. Guarda los resultados

# Guardar
# h5ad es un formato HDF5 optimizado para matrices grandes (genes x muestras).
adata.write(RESULTS_DIR / 'gene_counts.h5ad')

print("✅ Datos guardados en 'gene_counts.h5ad'")

# %% 4. Extraer datos del objeto AnnData

# Extraer datos

# - Si adata.X es una matriz dispersa (sparse), se convierte a densa con toarray().
#   AnnData usa matrices dispersas para ahorrar memoria si hay muchos ceros, pero para 4 muestras es mejor densa para visualización.

# - index=adata.obs_names: nombres de muestras (SRR...).

# - columns=adata.var_names: nombres de genes (ENSG...).

counts_df = pd.DataFrame(
    adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X,
    index=adata.obs_names,
    columns=adata.var_names
)

# Verificar datos

print("="*50)

print("VERIFICACIÓN DE DATOS")

print("="*50)

print(f"Matriz: {counts_df.shape}")

print(f"Suma total: {counts_df.sum().sum():,.0f}")

print(f"Conteos > 0: {(counts_df > 0).sum().sum():,}")

print("\nTop 10 genes más expresados:")
# Suma sobre el eje 0 (muestras) para ver qué genes tienen más lecturas en total.
top_genes = counts_df.sum(axis=0).sort_values(ascending=False).head(10)

for gene, count in top_genes.items():

    print(f" {gene}: {count:,.0f}")

print("\nConteos por muestra:")

for sample in counts_df.index:

    total = counts_df.loc[sample].sum()
    
    print(f" {sample}: {total:,.0f}")

# %% Extracion de resultados en formato CSV

# Guardar SOLO CSV

csv_files = {
    'counts': RESULTS_DIR / "counts_matrix.csv",
    'tpm': RESULTS_DIR / "tpm_matrix.csv" if 'abundance' in adata.layers else None,
    'stats': RESULTS_DIR / "sample_stats.csv"
}

# Guardar conteos (matriz de conteos a nivel de gen; filas = muestras, columnas = genes)

counts_df.to_csv(csv_files['counts'])

print(f"\n💾 Conteos: {csv_files['counts']}")

# Guardar TPM si existe la capa 'abundance'
# TPM (Transcripts Per Million) normaliza por longitud del gen y profundidad de secuenciación.
# Es útil para comparar la expresión de dos genes diferentes DENTRO de la misma muestra, pero NO para DGE.
if 'abundance' in adata.layers:

    tpm_df = pd.DataFrame(
        adata.layers['abundance'].toarray() if hasattr(adata.layers['abundance'], 'toarray') else adata.layers['abundance'],
        index=adata.obs_names,
        columns=adata.var_names
    )
    
    tpm_df.to_csv(csv_files['tpm'])
    
    print(f"💾 TPM: {csv_files['tpm']}")

# Guardar estadísticas por muestra

# - total_counts: suma total de conteos por muestra (Library Size).

# - expressed_genes: número de genes con conteo > 0 por muestra.

# - mean_expression: media de conteos por gen en cada muestra.

stats_df = pd.DataFrame({
    'total_counts': counts_df.sum(axis=1),
    'expressed_genes': (counts_df > 0).sum(axis=1),
    'mean_expression': counts_df.mean(axis=1)
})

stats_df.to_csv(csv_files['stats'])

print(f"💾 Estadísticas: {csv_files['stats']}")

print(f"\n✅ ¡Todos los archivos guardados en {RESULTS_DIR}!")


