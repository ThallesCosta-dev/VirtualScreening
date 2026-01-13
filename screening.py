import os
import subprocess
import re
import requests
import polars as pl
import numpy as np
import shutil
from rdkit import Chem
from rdkit.Chem import Descriptors

def check_tool(path):
    """Verifica se o executável existe e é funcional."""
    if not path:
        return False
    return shutil.which(path) is not None or os.path.exists(path)

# --- INTERFACE DE ENTRADA (PROMPTS) ---
print("=== CONFIGURAÇÃO DO PIPELINE DE BIOINFORMÁTICA ===")
protein_input = input("Digite o ID do PDB (ex: 6LU7) ou o caminho do arquivo .pdb: ").strip()

# Caminhos dos Binários
path_vina = input("Caminho do executável Vina-GPU (ou deixe vazio): ").strip()
path_gnina = input("Caminho do executável Gnina (ou deixe vazio): ").strip()
path_smina = input("Caminho do executável Smina (ou deixe vazio): ").strip()

# Configuração do GridBox (Baseado no seu AutoDock Tools)
print("\n--- Configurações do GridBox (Vindas do AutoDock Tools) ---")
c_x = input("Center X: ")
c_y = input("Center Y: ")
c_z = input("Center Z: ")
s_x = input("Size X (padrão 20): ") or "20"
s_y = input("Size Y (padrão 20): ") or "20"
s_z = input("Size Z (padrão 20): ") or "20"

# Validação das Ferramentas
tools = {
    "Vina-GPU": {"path": path_vina, "active": check_tool(path_vina)},
    "Gnina": {"path": path_gnina, "active": check_tool(path_gnina)},
    "Smina": {"path": path_smina, "active": check_tool(path_smina)}
}

for name, info in tools.items():
    if not info["active"]:
        print(f"[!] {name} OFF - Realize o download e configure o path corretamente para usar este módulo.")

# --- CONFIGURAÇÕES DE DIRETÓRIO ---
LIGANDS_DIR = "ligands_raw"
WORK_DIR = "screening_workdir"
os.makedirs(f"{WORK_DIR}/results", exist_ok=True)
os.makedirs(f"{WORK_DIR}/ligands", exist_ok=True)

# --- FUNÇÕES DE SUPORTE ---
def prepare_protein(target):
    if os.path.exists(target):
        return target
    print(f"[*] Baixando PDB: {target}...")
    r = requests.get(f"https://files.rcsb.org/download/{target}.pdb")
    path = f"{WORK_DIR}/{target}.pdb"
    with open(path, "wb") as f:
        f.write(r.content)
    # Converter para PDBQT
    path_qt = path + "qt"
    subprocess.run(["obabel", path, "-O", path_qt, "-xr", "-p", "7.4"], capture_output=True)
    return path_qt

def parse_generic_score(log_path):
    try:
        with open(log_path, "r") as f:
            content = f.read()
            match = re.search(r"\s+1\s+(-?\d+\.\d+)", content)
            return float(match.group(1)) if match else None
    except: return None

# --- LOOP PRINCIPAL ---
def run_pipeline():
    receptor_qt = prepare_protein(protein_input)
    results_list = []

    # Verificar se temos ao menos um arquivo SDF
    if not os.path.exists(LIGANDS_DIR):
        print(f"[ERRO] Pasta {LIGANDS_DIR} não encontrada!")
        return

    for sdf_file in os.listdir(LIGANDS_DIR):
        if not sdf_file.endswith(".sdf"): continue
        
        suppl = Chem.SDMolSupplier(os.path.join(LIGANDS_DIR, sdf_file))
        for i, mol in enumerate(suppl):
            if mol is None: continue
            
            mol_id = mol.GetProp("_Name") if mol.HasProp("_Name") else f"mol_{i}"
            smiles = Chem.MolToSmiles(mol)
            
            # Filtro Lipinski
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            if mw > 500 or logp > 5: continue # Filtro rígido para exemplo

            # Preparar Ligante
            lig_sdf = f"{WORK_DIR}/ligands/{mol_id}.sdf"
            lig_qt = f"{WORK_DIR}/ligands/{mol_id}.pdbqt"
            writer = Chem.SDWriter(lig_sdf)
            writer.write(mol)
            writer.close()
            subprocess.run(["obabel", lig_sdf, "-O", lig_qt, "-p", "7.4"], capture_output=True)

            res = {"ID": mol_id, "SMILES": smiles, "MW": mw, "LogP": logp}

            # Execução condicional
            if tools["Vina-GPU"]["active"]:
                log = f"{WORK_DIR}/results/{mol_id}_vina.log"
                subprocess.run([tools["Vina-GPU"]["path"], "--receptor", receptor_qt, "--ligand", lig_qt,
                                "--center_x", c_x, "--center_y", c_y, "--center_z", c_z,
                                "--size_x", s_x, "--size_y", s_y, "--size_z", s_z, "--log", log], capture_output=True)
                res["Vina_Score"] = parse_generic_score(log)

            if tools["Smina"]["active"]:
                log = f"{WORK_DIR}/results/{mol_id}_smina.log"
                subprocess.run([tools["Smina"]["path"], "-r", receptor_qt, "-l", lig_qt,
                                "--center_x", c_x, "--center_y", c_y, "--center_z", c_z,
                                "--size_x", s_x, "--size_y", s_y, "--size_z", s_z], stdout=open(log, "w"))
                res["Smina_Score"] = parse_generic_score(log)

            results_list.append(res)

    # Gerar Excel
    df = pl.DataFrame(results_list)
    df.write_excel("Resultado_Virtual_Screening.xlsx")
    print("\n[✓] Processo concluído! Relatório salvo em Resultado_Virtual_Screening.xlsx")

if __name__ == "__main__":
    if any(info["active"] for info in tools.values()):
        run_pipeline()
    else:
        print("[X] Nenhuma ferramenta de docking ativa. Abortando.")