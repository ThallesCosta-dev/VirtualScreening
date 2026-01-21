#!/usr/bin/env python3
"""
Instalador automático para Open Babel e MGLTools
"""

import os
import sys
import subprocess
import urllib.request
import zipfile
import shutil

def download_and_extract(url, extract_to):
    """Baixa e extrai um arquivo ZIP"""
    print(f"[*] Baixando {url}...")
    try:
        urllib.request.urlretrieve(url, "temp.zip")
        print("[✓] Download concluído")
        
        with zipfile.ZipFile("temp.zip", 'r') as zip_ref:
            zip_ref.extractall(extract_to)
        print("[✓] Extração concluída")
        
        os.remove("temp.zip")
        return True
    except Exception as e:
        print(f"[!] Erro: {e}")
        return False

def install_openbabel():
    """Instala Open Babel para Windows"""
    print("\n=== Instalando Open Babel ===")
    
    # Download Open Babel Windows
    url = "https://sourceforge.net/projects/openbabel/files/openbabel/3.1.1/OpenBabel-3.1.1-win64.exe/download"
    
    if download_and_extract(url, "openbabel_temp"):
        # Mover para pasta local
        if os.path.exists("openbabel_temp"):
            for item in os.listdir("openbabel_temp"):
                s = os.path.join("openbabel_temp", item)
                d = os.path.join(os.getcwd(), item)
                if os.path.isdir(s):
                    if os.path.exists(d):
                        shutil.rmtree(d)
                    shutil.copytree(s, d)
                else:
                    shutil.copy2(s, d)
            shutil.rmtree("openbabel_temp")
            print("[✓] Open Babel instalado na pasta local")
            return True
    return False

def install_mgltools():
    """Instala MGLTools via pip"""
    print("\n=== Instalando MGLTools ===")
    
    try:
        # Tentar instalar via conda primeiro
        subprocess.run([sys.executable, "-m", "pip", "install", "conda"], capture_output=True)
        
        # Instalar MGLTools
        result = subprocess.run([
            sys.executable, "-m", "pip", "install", 
            "https://github.com/ccsb-scripps/AutoDockTools/archive/refs/tags/v1.5.7.tar.gz"
        ], capture_output=True, text=True)
        
        if result.returncode == 0:
            print("[✓] MGLTools instalado via pip")
            return True
        else:
            print(f"[!] Erro na instalação: {result.stderr}")
            return False
            
    except Exception as e:
        print(f"[!] Erro ao instalar MGLTools: {e}")
        return False

def test_installations():
    """Testa se as ferramentas estão funcionando"""
    print("\n=== Testando instalações ===")
    
    # Testar Open Babel
    obabel_paths = [
        "obabel.exe",
        "OpenBabel/obabel.exe",
        "openbabel-3.1.1/obabel.exe"
    ]
    
    obabel_found = False
    for path in obabel_paths:
        if os.path.exists(path):
            try:
                result = subprocess.run([path, "--version"], capture_output=True, text=True)
                if result.returncode == 0:
                    print(f"[✓] Open Babel encontrado: {path}")
                    print(f"    Versão: {result.stdout.strip()}")
                    obabel_found = True
                    break
            except:
                continue
    
    if not obabel_found:
        print("[!] Open Babel não encontrado ou não funcional")
    
    # Testar MGLTools
    try:
        result = subprocess.run([
            sys.executable, "-c", 
            "from AutoDockTools.Utilities24.prepare_receptor4 import prepare_receptor4; print('MGLTools OK')"
        ], capture_output=True, text=True)
        
        if result.returncode == 0 and "MGLTools OK" in result.stdout:
            print("[✓] MGLTools funcionando")
            return True
        else:
            print(f"[!] MGLTools com problemas: {result.stderr}")
    except Exception as e:
        print(f"[!] Erro ao testar MGLTools: {e}")
    
    return obabel_found

def main():
    print("=== Instalador Open Babel + MGLTools ===")
    print("Este script instalará as ferramentas necessárias para conversão PDB→PDBQT")
    
    # Instalar Open Babel
    obabel_ok = install_openbabel()
    
    # Instalar MGLTools
    mgl_ok = install_mgltools()
    
    # Testar instalações
    test_installations()
    
    if obabel_ok or mgl_ok:
        print("\n[✓] Instalação concluída! Execute o screening.py novamente.")
    else:
        print("\n[!] Alguns problemas ocorreram. Verifique os erros acima.")

if __name__ == "__main__":
    main()
