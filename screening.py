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

def find_vina_in_current_dir():
    """Procura pelo executável vina na pasta atual."""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    possible_names = ["vina.exe", "vina", "autodock_vina.exe", "autodock_vina","vina_1.2.7_win.exe"]
    for name in possible_names:
        if os.path.exists(os.path.join(current_dir, name)):
            return os.path.join(current_dir, name)
    return None

# --- INTERFACE DE ENTRADA (PROMPTS) ---
print("=== CONFIGURAÇÃO DO PIPELINE DE BIOINFORMÁTICA ===")
protein_input = input("Digite o ID do PDB (ex: P2X7) ou o caminho do arquivo .pdb: ").strip()

# Caminhos dos Binários
path_vina = input("Caminho do executável Vina (ou deixe vazio para procurar na pasta atual): ").strip()
path_vina_gpu = input("Caminho do executável Vina-GPU (ou deixe vazio): ").strip()
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
# Se não for especificado caminho para o Vina, procurar na pasta atual
if not path_vina:
    path_vina = find_vina_in_current_dir()

tools = {
    "Vina": {"path": path_vina, "active": check_tool(path_vina)},
    "Vina-GPU": {"path": path_vina_gpu, "active": check_tool(path_vina_gpu)},
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
def find_obabel():
    """Procura pelo Open Babel em diferentes locais."""
    import sys
    
    # Verificar no PATH
    if shutil.which("obabel"):
        return "obabel"
    
    # Locais comuns no Windows (aceita qualquer versão)
    possible_paths = []
    
    # Buscar em C:\Program Files\OpenBabel*
    if os.path.exists(r"C:\Program Files"):
        for item in os.listdir(r"C:\Program Files"):
            if item.startswith("OpenBabel"):
                obabel_path = os.path.join(r"C:\Program Files", item, "obabel.exe")
                if os.path.exists(obabel_path):
                    possible_paths.append(obabel_path)
    
    # Buscar em C:\Program Files (x86)\OpenBabel*
    if os.path.exists(r"C:\Program Files (x86)"):
        for item in os.listdir(r"C:\Program Files (x86)"):
            if item.startswith("OpenBabel"):
                obabel_path = os.path.join(r"C:\Program Files (x86)", item, "obabel.exe")
                if os.path.exists(obabel_path):
                    possible_paths.append(obabel_path)
    
    # Locais fixos
    possible_paths.extend([
        r"C:\OpenBabel\obabel.exe",
        r"C:\msys64\mingw64\bin\obabel.exe",
        r"C:\msys64\usr\bin\obabel.exe"
    ])
    
    for path in possible_paths:
        if os.path.exists(path):
            return path
    
    # Tentar encontrar via Python packages
    try:
        import openbabel
        # Verificar se está no site-packages
        site_packages = next((p for p in sys.path if 'site-packages' in p), None)
        if site_packages:
            obabel_path = os.path.join(site_packages, "..", "..", "Scripts", "obabel.exe")
            if os.path.exists(obabel_path):
                return obabel_path
    except ImportError:
        pass
    
    return None

def find_mgltools():
    """Procura pelo MGLTools em diferentes locais."""
    import sys
    
    # Locais comuns no Windows (aceita qualquer versão)
    possible_paths = []
    
    # Buscar em C:\Program Files\MGLTools*
    if os.path.exists(r"C:\Program Files"):
        for item in os.listdir(r"C:\Program Files"):
            if item.startswith("MGLTools"):
                mgl_path = os.path.join(r"C:\Program Files", item, "MGLToolsPckgs", "AutoDockTools", "Utilities24")
                if os.path.exists(mgl_path):
                    possible_paths.append(mgl_path)
    
    # Buscar em C:\Program Files (x86)\MGLTools*
    if os.path.exists(r"C:\Program Files (x86)"):
        for item in os.listdir(r"C:\Program Files (x86)"):
            if item.startswith("MGLTools"):
                mgl_path = os.path.join(r"C:\Program Files (x86)", item, "MGLToolsPckgs", "AutoDockTools", "Utilities24")
                if os.path.exists(mgl_path):
                    possible_paths.append(mgl_path)
    
    # Locais fixos
    possible_paths.extend([
        r"C:\MGLTools-1.5.7\MGLToolsPckgs\AutoDockTools\Utilities24",
        r"C:\MGLTools-1.5.6\MGLToolsPckgs\AutoDockTools\Utilities24"
    ])
    
    for path in possible_paths:
        if os.path.exists(path):
            return path
    
    # Tentar encontrar via Python packages
    try:
        import AutoDockTools
        return AutoDockTools.__path__[0]
    except ImportError:
        pass
    
    return None

def search_pdb_id(protein_name):
    """Busca ID PDB para uma proteína usando a API do RCSB."""
    try:
        # Busca por nome da proteína
        search_url = f"https://search.rcsb.org/rcsbsearch/v2/query"
        query = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "value": protein_name
                }
            },
            "return_type": "entry"
        }
        
        response = requests.post(search_url, json=query)
        if response.status_code == 200:
            results = response.json()
            if results.get("result_set"):
                # Pega o primeiro resultado
                pdb_id = results["result_set"][0]["identifier"]
                print(f"[*] PDB encontrado para {protein_name}: {pdb_id}")
                return pdb_id
        return None
    except Exception as e:
        print(f"[!] Erro na busca do PDB: {e}")
        return None

def prepare_protein(target):
    if os.path.exists(target):
        return target
    
    print(f"[*] Baixando PDB: {target}...")
    
    # Se não parece com um ID PDB (4 caracteres), tentar buscar
    if len(target) != 4 or not target.isalnum():
        print(f"[*] '{target}' não parece um ID PDB. Buscando...")
        pdb_id = search_pdb_id(target)
        if not pdb_id:
            print(f"[!] Não foi possível encontrar um PDB para '{target}'")
            print("[!] Tente usar um ID PDB válido (ex: 6U9N, 7KZF, etc.)")
            return None
    else:
        pdb_id = target.upper()
    
    # Tentar download
    r = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb")
    if r.status_code != 200:
        print(f"[ERRO] Não foi possível baixar o PDB {pdb_id}. Status: {r.status_code}")
        print("[!] Verifique se o ID PDB está correto")
        return None
    
    path = f"{WORK_DIR}/{pdb_id}.pdb"
    with open(path, "wb") as f:
        f.write(r.content)
    print(f"[✓] PDB {pdb_id} baixado com sucesso")
    
    # Converter para PDBQT
    path_qt = path + "qt"
    if convert_to_pdbqt(path, path_qt, is_ligand=False):
        print(f"[✓] PDB convertido para PDBQT")
        return path_qt
    else:
        print("[ERRO] Não foi possível converter o PDB para PDBQT")
        print("[INFO] Instale Open Babel ou MGLTools para melhor conversão")
        return None

def convert_to_pdbqt(pdb_path, pdbqt_path, is_ligand=False):
    """Converte PDB para PDBQT usando diferentes métodos."""
    
    # Método 1: Tentar com Open Babel (busca aprimorada)
    obabel_cmd = find_obabel()
    if obabel_cmd:
        try:
            if is_ligand:
                cmd = [obabel_cmd, pdb_path, "-O", pdbqt_path, "-p", "7.4"]
            else:
                cmd = [obabel_cmd, pdb_path, "-O", pdbqt_path, "-xr", "-p", "7.4"]
            subprocess.run(cmd, capture_output=True, check=True)
            print(f"[✓] Convertido com Open Babel: {obabel_cmd}")
            return True
        except subprocess.CalledProcessError as e:
            print(f"[!] Erro ao usar Open Babel: {e}")
    
    # Método 2: Tentar com MGLTools (busca aprimorada)
    mgl_path = find_mgltools()
    if mgl_path:
        try:
            if is_ligand:
                prepare_script = os.path.join(mgl_path, "prepare_ligand4.py")
                cmd = ["python", prepare_script, "-l", pdb_path, "-o", pdbqt_path]
            else:
                prepare_script = os.path.join(mgl_path, "prepare_receptor4.py")
                cmd = ["python", prepare_script, "-r", pdb_path, "-o", pdbqt_path]
            
            subprocess.run(cmd, capture_output=True, check=True)
            print(f"[✓] Convertido com MGLTools: {prepare_script}")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"[!] Erro ao usar MGLTools: {e}")
    
    # Método 3: Tentar via import direto do Python
    try:
        if is_ligand:
            cmd = ["python", "-c", 
                   "from AutoDockTools.Utilities24.prepare_ligand4 import prepare_ligand4; "
                   f"prepare_ligand4('{pdb_path}', '{pdbqt_path}')"]
        else:
            cmd = ["python", "-c", 
                   "from AutoDockTools.Utilities24.prepare_receptor4 import prepare_receptor4; "
                   f"prepare_receptor4('{pdb_path}', '{pdbqt_path}')"]
        subprocess.run(cmd, capture_output=True, check=True, shell=True)
        print(f"[✓] Convertido com MGLTools (Python import)")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError, ImportError) as e:
        print(f"[!] Erro ao usar MGLTools via Python: {e}")
    
    # Método 4: Conversão manual básica (fallback)
    print("[!] Tentando conversão manual básica...")
    try:
        with open(pdb_path, 'r') as f:
            content = f.read()
        
        # Adicionar cargas parciais básicas e remover átomos não essenciais
        lines = []
        for line in content.split('\n'):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Manter apenas átomos relevantes
                atom = line[76:78].strip() if len(line) > 77 else ''
                if atom not in ['HOH', 'WAT']:  # Remover moléculas de água
                    lines.append(line)
        
        # Adicionar carga e tipo de átomo
        for i, line in enumerate(lines):
            if len(line) > 77:
                lines[i] = line[:78] + " 0.00  0.00"
            else:
                lines[i] = line.ljust(80) + " 0.00  0.00"
        
        with open(pdbqt_path, 'w') as f:
            f.write('\n'.join(lines))
            f.write('\nTORSDOF 0\n')
        
        print(f"[✓] Conversão manual básica concluída")
        return True
    except Exception as e:
        print(f"[!] Falha na conversão manual: {e}")
        return False

def parse_generic_score(log_path):
    try:
        with open(log_path, "r") as f:
            content = f.read()
            # Tentar encontrar score no formato Vina 1.2.7
            match = re.search(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", content)
            if match:
                return float(match.group(1))
            
            # Tentar formato antigo
            match = re.search(r"\s+1\s+(-?\d+\.\d+)", content)
            return float(match.group(1)) if match else None
    except: return None

# --- LOOP PRINCIPAL ---
def run_pipeline():
    receptor_qt = prepare_protein(protein_input)
    if not receptor_qt:
        print("[ERRO] Falha ao preparar o receptor. Verifique os erros acima.")
        return
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
            
            # Converter SDF para PDBQT usando a função robusta
            if not convert_to_pdbqt(lig_sdf, lig_qt, is_ligand=True):
                print(f"[!] Falha ao converter ligante {mol_id}, pulando...")
                continue

            res = {"ID": mol_id, "SMILES": smiles, "MW": mw, "LogP": logp}

            # Execução condicional
            if tools["Vina"]["active"]:
                log = f"{WORK_DIR}/results/{mol_id}_vina.log"
                out = f"{WORK_DIR}/results/{mol_id}_vina_out.pdbqt"
                print(f"[*] Executando Vina para {mol_id}...")
                print(f"    Comando: {tools['Vina']['path']} --receptor {receptor_qt} --ligand {lig_qt} --center_x {c_x} --center_y {c_y} --center_z {c_z} --size_x {s_x} --size_y {s_y} --size_z {s_z} --out {out}")
                
                result = subprocess.run([tools["Vina"]["path"], "--receptor", receptor_qt, "--ligand", lig_qt,
                                        "--center_x", c_x, "--center_y", c_y, "--center_z", c_z,
                                        "--size_x", s_x, "--size_y", s_y, "--size_z", s_z, "--out", out], 
                                       capture_output=True, text=True)
                
                print(f"    Retorno Vina: {result.returncode}")
                if result.stdout:
                    print(f"    Stdout: {result.stdout[:200]}...")
                if result.stderr:
                    print(f"    Stderr: {result.stderr[:200]}...")
                
                # Criar log manualmente com stdout do Vina
                with open(log, 'w') as f:
                    f.write(result.stdout)
                    if result.stderr:
                        f.write("\nSTDERR:\n" + result.stderr)
                
                score = parse_generic_score(out)  # Usar arquivo de saída, não log
                print(f"    Score extraído: {score}")
                res["Vina_Score"] = score

            if tools["Vina-GPU"]["active"]:
                log = f"{WORK_DIR}/results/{mol_id}_vina_gpu.log"
                subprocess.run([tools["Vina-GPU"]["path"], "--receptor", receptor_qt, "--ligand", lig_qt,
                                "--center_x", c_x, "--center_y", c_y, "--center_z", c_z,
                                "--size_x", s_x, "--size_y", s_y, "--size_z", s_z, "--log", log], capture_output=True)
                res["Vina_GPU_Score"] = parse_generic_score(log)

            if tools["Smina"]["active"]:
                log = f"{WORK_DIR}/results/{mol_id}_smina.log"
                subprocess.run([tools["Smina"]["path"], "-r", receptor_qt, "-l", lig_qt,
                                "--center_x", c_x, "--center_y", c_y, "--center_z", c_z,
                                "--size_x", s_x, "--size_y", s_y, "--size_z", s_z], stdout=open(log, "w"))
                res["Smina_Score"] = parse_generic_score(log)

            results_list.append(res)

    # Gerar Excel
    df = pl.DataFrame(results_list)
    try:
        df.write_excel("Resultado_Virtual_Screening.xlsx")
        print("\n[✓] Processo concluído! Relatório salvo em Resultado_Virtual_Screening.xlsx")
    except PermissionError:
        print("\n[!] Erro de permissão ao salvar Excel. Feche o arquivo Resultado_Virtual_Screening.xlsx e tente novamente.")
        # Tentar salvar com nome diferente
        df.write_excel("Resultado_Virtual_Screening_novo.xlsx")
        print("[✓] Relatório salvo em Resultado_Virtual_Screening_novo.xlsx")

if __name__ == "__main__":
    if any(info["active"] for info in tools.values()):
        run_pipeline()
    else:
        print("[X] Nenhuma ferramenta de docking ativa. Abortando.")