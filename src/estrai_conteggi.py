import pandas as pd
import sys
import os

# Percorso del file compresso contenente i metadati cellulari GEO (GSE278688)
file_metadati = "data/raw/GSE278688/GSE278688_sc_metadata.csv.gz"

if not os.path.exists(file_metadati):
    print(f"ERRORE: Non trovo il file {file_metadati}.")
    sys.exit(1)

print("Caricamento metadati...")
obs = pd.read_csv(file_metadati)

# Esclusione delle popolazioni non target o a bassa numerosita' (Platelet ed Endocrine)
EXCLUDE_ALWAYS = ["Platelet", "Endocrine"]
obs = obs[~obs["all_celltype"].isin(EXCLUDE_ALWAYS)].copy()

# Identificazione dei 5 pazienti dotati di campioni tissutali appaiati (Tumor e Adjacent_normal)
tis = pd.crosstab(obs["patients"], obs["tissue"])
paired = sorted(
    tis.index[(tis.get("Tumor", 0) > 0) & (tis.get("Adjacent_normal", 0) > 0)]
)

# Selezione delle cellule appartenenti ai due compartimenti solidi per la coorte appaiata
t1all = obs[obs["tissue"].isin(["Tumor", "Adjacent_normal"])]
t1p = t1all[t1all["patients"].isin(paired)]

# Aggregazione e ordinamento dei conteggi cellulari per singolo paziente
conteggi = t1p["patients"].value_counts().sort_index()

print("\n--- RIEPILOGO CONTEGGI PER PAZIENTE (COORTE APPAIATA) ---")
print("Paziente\tCellule totali")
print("-" * 30)
for paziente, n_celle in conteggi.items():
    print(f"{paziente}\t\t{n_celle:,}".replace(",", "."))
print("-" * 30)
print(f"TOTALE\t\t{conteggi.sum():,}".replace(",", "."))
