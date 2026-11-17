"""
Boltz Co-folding Job Builder
============================
Interface gráfica (roda no navegador) para gerar arquivos de input (.yaml)
para o Boltz (co-folding proteína-ligante), nos modos:

  1) Proteína única + Ligante único
  2) Screening (1 proteína x vários ligantes)
  3) Multidocking (várias proteínas x vários ligantes)

Como usar:
    pip install -r requirements.txt
    python app.py
Depois abra o link que aparecer no terminal (geralmente http://127.0.0.1:7860).

O app NÃO roda o Boltz em si — ele gera uma pasta .zip com os arquivos .yaml
prontos, um script run.sh/run_all.sh para rodar tudo de uma vez, e (nos modos
screening/multidocking) um script compile_results.py que junta os resultados
num único CSV ranqueado.
"""

import os
import re
import csv
import json
import shutil
import tempfile
import zipfile
from itertools import product

import gradio as gr
import pandas as pd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def sanitize(name: str) -> str:
    """Transforma um nome livre em algo seguro para nome de arquivo/pasta."""
    name = (name or "").strip()
    name = re.sub(r"\s+", "_", name)
    name = re.sub(r"[^A-Za-z0-9_\-]", "", name)
    return name or "sem_nome"


def clean_sequence(seq: str) -> str:
    """Remove espaços, quebras de linha e números (ex: colado de um FASTA numerado)."""
    seq = (seq or "").strip()
    # remove cabeçalho fasta se o usuário colou por engano
    lines = [l for l in seq.splitlines() if not l.startswith(">")]
    seq = "".join(lines)
    seq = re.sub(r"\s+", "", seq)
    return seq.upper()


def clean_smiles(smiles: str) -> str:
    return (smiles or "").strip()


HEADER_WORDS = {
    "nome", "name", "id", "ligante", "ligand", "proteina", "protein",
    "sequencia", "sequence", "seq", "smiles",
}


def parse_pairs(text: str, csv_file) -> list[tuple[str, str]]:
    """
    Extrai pares (nome, valor) de um textarea colado (formato 'nome,valor' por
    linha) ou de um CSV enviado (procura colunas por nome, senão usa as duas
    primeiras colunas). Ignora linhas de cabeçalho.
    """
    pairs: list[tuple[str, str]] = []

    if csv_file is not None:
        path = csv_file.name if hasattr(csv_file, "name") else csv_file
        df = pd.read_csv(path)
        cols_lower = {c.lower().strip(): c for c in df.columns}
        name_candidates = ["nome", "name", "id", "protein_name", "ligand_name"]
        val_candidates = ["smiles", "sequence", "sequencia", "seq", "ligand", "ligante"]
        name_col = next((cols_lower[c] for c in name_candidates if c in cols_lower), df.columns[0])
        val_col = next((cols_lower[c] for c in val_candidates if c in cols_lower and cols_lower[c] != name_col),
                        df.columns[1] if len(df.columns) > 1 else df.columns[0])
        for _, row in df.iterrows():
            n = str(row[name_col]).strip()
            v = str(row[val_col]).strip()
            if n and v and n.lower() != "nan" and v.lower() != "nan":
                pairs.append((n, v))
        return pairs

    for raw_line in (text or "").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if "," in line:
            n, v = line.split(",", 1)
        elif "\t" in line:
            n, v = line.split("\t", 1)
        elif ";" in line:
            n, v = line.split(";", 1)
        else:
            continue
        n, v = n.strip(), v.strip()
        if n.lower() in HEADER_WORDS and v.lower() in HEADER_WORDS:
            continue  # linha de cabeçalho colada por engano
        if n and v:
            pairs.append((n, v))
    return pairs


# ---------------------------------------------------------------------------
# Geração do YAML do Boltz
# ---------------------------------------------------------------------------

def build_yaml(protein_seq: str, ligand_value: str, ligand_is_ccd: bool,
                use_msa_server: bool, compute_affinity: bool) -> str:
    lines = ["version: 1", "sequences:", "  - protein:", "      id: A",
             f"      sequence: {protein_seq}"]
    if not use_msa_server:
        lines.append("      msa: empty")
    lines.append("  - ligand:")
    lines.append("      id: B")
    if ligand_is_ccd:
        lines.append(f"      ccd: {ligand_value}")
    else:
        lines.append(f"      smiles: '{ligand_value}'")
    if compute_affinity:
        lines.append("properties:")
        lines.append("  - affinity:")
        lines.append("      binder: B")
    return "\n".join(lines) + "\n"


README_TEMPLATE = """\
COMO RODAR O BOLTZ
===================

1) Instale o Boltz (uma vez só, precisa de GPU NVIDIA para ser rápido):
   pip install boltz -U

2) Descompacte este .zip e entre na pasta pelo terminal:
   cd {folder_hint}

3) Rode o comando abaixo (ele já está salvo no arquivo {run_script}):
   bash {run_script}

   - O comando usa "--use_msa_server", que gera o alinhamento (MSA) da
     proteína automaticamente por um servidor público (precisa de internet).
   - Os resultados aparecem na pasta "resultados/predictions/".
   - Cada predição tem um arquivo .cif (estrutura 3D) e um .json de confiança
     (e de afinidade, se você marcou essa opção).

{extra}
Dúvidas comuns:
- "boltz: command not found" -> rode "pip install boltz -U" de novo e confira
  se o ambiente/virtualenv certo está ativado.
- Sem GPU -> o comando funciona mas fica bem mais lento (pode usar
  --accelerator cpu, mas não é recomendado para proteínas grandes).
"""

COMPILE_RESULTS_SCRIPT = '''\
"""
Junta os resultados de várias predições do Boltz num único CSV ranqueado.
Rode isso DEPOIS que o "boltz predict" terminar, de dentro desta mesma pasta:

    python compile_results.py

Gera o arquivo resumo_resultados.csv com uma linha por job, já ordenado
(pelo valor de afinidade quando disponível, senão pelo confidence_score).
"""
import csv
import json
import os

BASE = os.path.join("resultados", "predictions")
MANIFEST = "manifest.csv"


def load_manifest():
    with open(MANIFEST, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main():
    if not os.path.isdir(BASE):
        print(f"Pasta '{BASE}' não encontrada. Rode 'boltz predict' primeiro.")
        return

    manifest = load_manifest()
    rows_out = []

    for row in manifest:
        job = row["job_id"]
        pred_dir = os.path.join(BASE, job)
        conf_path = os.path.join(pred_dir, f"confidence_{job}_model_0.json")
        aff_path = os.path.join(pred_dir, f"affinity_{job}.json")

        result = dict(row)
        result["confidence_score"] = None
        result["iptm"] = None
        result["ligand_iptm"] = None
        result["affinity_pred_value"] = None
        result["affinity_probability_binary"] = None

        if os.path.exists(conf_path):
            with open(conf_path) as cf:
                conf = json.load(cf)
            result["confidence_score"] = conf.get("confidence_score")
            result["iptm"] = conf.get("iptm")
            result["ligand_iptm"] = conf.get("ligand_iptm")
        else:
            print(f"[aviso] sem resultado de confiança para: {job}")

        if os.path.exists(aff_path):
            with open(aff_path) as af:
                aff = json.load(af)
            result["affinity_pred_value"] = aff.get("affinity_pred_value")
            result["affinity_probability_binary"] = aff.get("affinity_probability_binary")

        rows_out.append(result)

    def sort_key(r):
        if r.get("affinity_pred_value") is not None:
            return (0, r["affinity_pred_value"])
        return (1, -(r.get("confidence_score") or 0))

    rows_out.sort(key=sort_key)

    if not rows_out:
        print("Nenhum resultado encontrado.")
        return

    fieldnames = list(rows_out[0].keys())
    with open("resumo_resultados.csv", "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"Pronto! Resumo salvo em resumo_resultados.csv ({len(rows_out)} jobs).")
    print("Dica: abra esse CSV no Excel/LibreOffice e ordene pela coluna de afinidade ou confidence_score.")


if __name__ == "__main__":
    main()
'''


def make_zip(tmpdir: str, zip_name: str) -> str:
    zip_path = os.path.join(tempfile.gettempdir(), zip_name)
    if os.path.exists(zip_path):
        os.remove(zip_path)
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for root, _, files in os.walk(tmpdir):
            for f in files:
                full = os.path.join(root, f)
                rel = os.path.relpath(full, tmpdir)
                zf.write(full, rel)
    return zip_path


# ---------------------------------------------------------------------------
# Modo 1: Proteína única + Ligante único
# ---------------------------------------------------------------------------

def gerar_single(job_name, protein_seq, ligand_type, ligand_value,
                  use_msa_server, compute_affinity):
    protein_seq = clean_sequence(protein_seq)
    ligand_value = clean_smiles(ligand_value)
    job = sanitize(job_name) or "job"

    if not protein_seq:
        return None, "⚠️ Cole a sequência da proteína antes de gerar."
    if not ligand_value:
        return None, "⚠️ Informe o SMILES ou código CCD do ligante."

    is_ccd = ligand_type == "Código CCD"
    yaml_str = build_yaml(protein_seq, ligand_value, is_ccd, use_msa_server, compute_affinity)

    tmpdir = tempfile.mkdtemp()
    with open(os.path.join(tmpdir, f"{job}.yaml"), "w") as f:
        f.write(yaml_str)

    flag = " --use_msa_server" if use_msa_server else ""
    run_sh = (
        "#!/usr/bin/env bash\n"
        "set -e\n"
        f"boltz predict {job}.yaml --out_dir resultados{flag}\n"
        "echo 'Pronto! Veja a pasta resultados/predictions/'\n"
    )
    with open(os.path.join(tmpdir, "run.sh"), "w") as f:
        f.write(run_sh)

    readme = README_TEMPLATE.format(
        folder_hint="(pasta que você descompactou)",
        run_script="run.sh",
        extra="" if use_msa_server else (
            "ATENÇÃO: você desmarcou o uso do servidor de MSA automático.\n"
            "Isso roda em 'modo sequência única', o que reduz a precisão.\n"
            "Recomendado deixar essa opção marcada.\n\n"
        ),
    )
    with open(os.path.join(tmpdir, "LEIA-ME.txt"), "w") as f:
        f.write(readme)

    zip_path = make_zip(tmpdir, f"boltz_single_{job}.zip")
    shutil.rmtree(tmpdir, ignore_errors=True)
    return zip_path, f"✅ Gerado! 1 arquivo YAML ({job}.yaml) pronto no zip."


# ---------------------------------------------------------------------------
# Modo 2: Screening (1 proteína x vários ligantes)
# ---------------------------------------------------------------------------

def gerar_screening(protein_name, protein_seq, ligantes_text, ligantes_csv,
                     use_msa_server, compute_affinity):
    protein_seq = clean_sequence(protein_seq)
    protein_name = sanitize(protein_name) or "proteina"

    if not protein_seq:
        return None, "⚠️ Cole a sequência da proteína antes de gerar."

    pairs = parse_pairs(ligantes_text, ligantes_csv)
    if not pairs:
        return None, "⚠️ Nenhum ligante encontrado. Cole linhas no formato 'nome,SMILES' ou envie um CSV."

    tmpdir = tempfile.mkdtemp()
    manifest_rows = []
    used_names = set()

    for lig_name, lig_smiles in pairs:
        job = sanitize(f"{protein_name}__{lig_name}")
        base_job = job
        i = 2
        while job in used_names:
            job = f"{base_job}_{i}"
            i += 1
        used_names.add(job)

        yaml_str = build_yaml(protein_seq, clean_smiles(lig_smiles), False,
                               use_msa_server, compute_affinity)
        with open(os.path.join(tmpdir, f"{job}.yaml"), "w") as f:
            f.write(yaml_str)

        manifest_rows.append({
            "job_id": job,
            "protein_name": protein_name,
            "ligand_name": lig_name,
            "ligand_smiles": lig_smiles,
        })

    with open(os.path.join(tmpdir, "manifest.csv"), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["job_id", "protein_name", "ligand_name", "ligand_smiles"])
        writer.writeheader()
        writer.writerows(manifest_rows)

    flag = " --use_msa_server" if use_msa_server else ""
    run_all = (
        "#!/usr/bin/env bash\n"
        "set -e\n"
        f"boltz predict . --out_dir resultados{flag}\n"
        "echo 'Pronto! Rode agora: python compile_results.py'\n"
    )
    with open(os.path.join(tmpdir, "run_all.sh"), "w") as f:
        f.write(run_all)

    with open(os.path.join(tmpdir, "compile_results.py"), "w") as f:
        f.write(COMPILE_RESULTS_SCRIPT)

    readme = README_TEMPLATE.format(
        folder_hint="(pasta que você descompactou)",
        run_script="run_all.sh",
        extra=(
            f"Este pacote tem {len(manifest_rows)} ligantes contra a proteína "
            f"'{protein_name}'.\n\n"
            "4) Depois que o boltz terminar, rode:\n"
            "   python compile_results.py\n"
            "   Isso cria 'resumo_resultados.csv' com todos os ligantes "
            "ranqueados\n   (por afinidade, se marcada, ou por confiança).\n\n"
        ),
    )
    with open(os.path.join(tmpdir, "LEIA-ME.txt"), "w") as f:
        f.write(readme)

    zip_path = make_zip(tmpdir, f"boltz_screening_{protein_name}.zip")
    shutil.rmtree(tmpdir, ignore_errors=True)
    return zip_path, f"✅ Gerados {len(manifest_rows)} arquivos YAML (1 proteína x {len(manifest_rows)} ligantes)."


# ---------------------------------------------------------------------------
# Modo 3: Multidocking (várias proteínas x vários ligantes)
# ---------------------------------------------------------------------------

def gerar_multidocking(proteinas_text, proteinas_csv, ligantes_text, ligantes_csv,
                        use_msa_server, compute_affinity):
    proteinas = parse_pairs(proteinas_text, proteinas_csv)
    ligantes = parse_pairs(ligantes_text, ligantes_csv)

    if not proteinas:
        return None, "⚠️ Nenhuma proteína encontrada. Cole linhas 'nome,sequência' ou envie um CSV."
    if not ligantes:
        return None, "⚠️ Nenhum ligante encontrado. Cole linhas 'nome,SMILES' ou envie um CSV."

    tmpdir = tempfile.mkdtemp()
    manifest_rows = []
    used_names = set()

    for (prot_name, prot_seq), (lig_name, lig_smiles) in product(proteinas, ligantes):
        job = sanitize(f"{prot_name}__{lig_name}")
        base_job = job
        i = 2
        while job in used_names:
            job = f"{base_job}_{i}"
            i += 1
        used_names.add(job)

        yaml_str = build_yaml(clean_sequence(prot_seq), clean_smiles(lig_smiles), False,
                               use_msa_server, compute_affinity)
        with open(os.path.join(tmpdir, f"{job}.yaml"), "w") as f:
            f.write(yaml_str)

        manifest_rows.append({
            "job_id": job,
            "protein_name": prot_name,
            "ligand_name": lig_name,
            "ligand_smiles": lig_smiles,
        })

    with open(os.path.join(tmpdir, "manifest.csv"), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["job_id", "protein_name", "ligand_name", "ligand_smiles"])
        writer.writeheader()
        writer.writerows(manifest_rows)

    flag = " --use_msa_server" if use_msa_server else ""
    run_all = (
        "#!/usr/bin/env bash\n"
        "set -e\n"
        f"boltz predict . --out_dir resultados{flag}\n"
        "echo 'Pronto! Rode agora: python compile_results.py'\n"
    )
    with open(os.path.join(tmpdir, "run_all.sh"), "w") as f:
        f.write(run_all)

    with open(os.path.join(tmpdir, "compile_results.py"), "w") as f:
        f.write(COMPILE_RESULTS_SCRIPT)

    readme = README_TEMPLATE.format(
        folder_hint="(pasta que você descompactou)",
        run_script="run_all.sh",
        extra=(
            f"Este pacote tem {len(proteinas)} proteína(s) x {len(ligantes)} "
            f"ligante(s) = {len(manifest_rows)} combinações.\n\n"
            "ATENÇÃO: cada combinação recalcula o MSA da proteína (se "
            "--use_msa_server estiver ligado), então isso pode demorar bem "
            "mais que o modo screening.\n\n"
            "4) Depois que o boltz terminar, rode:\n"
            "   python compile_results.py\n"
            "   Isso cria 'resumo_resultados.csv' com todas as combinações "
            "ranqueadas.\n\n"
        ),
    )
    with open(os.path.join(tmpdir, "LEIA-ME.txt"), "w") as f:
        f.write(readme)

    zip_path = make_zip(tmpdir, "boltz_multidocking.zip")
    shutil.rmtree(tmpdir, ignore_errors=True)
    return zip_path, (f"✅ Gerados {len(manifest_rows)} arquivos YAML "
                       f"({len(proteinas)} proteínas x {len(ligantes)} ligantes).")


# ---------------------------------------------------------------------------
# Interface Gradio
# ---------------------------------------------------------------------------

EXEMPLO_LIGANTES = "aspirina,CC(=O)Oc1ccccc1C(=O)O\nibuprofeno,CC(C)Cc1ccc(cc1)C(C)C(=O)O"
EXEMPLO_PROTEINAS = "minha_proteina,MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWELVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL"

with gr.Blocks(title="Boltz Co-folding Job Builder") as demo:
    gr.Markdown(
        "# 🧬 Boltz Co-folding — Gerador de Arquivos\n"
        "Gera os arquivos `.yaml` prontos para rodar o **Boltz** "
        "(predição de estrutura proteína-ligante). Escolha uma aba abaixo, "
        "preencha os campos e clique em **Gerar**. Você vai baixar um `.zip` "
        "com tudo pronto (incluindo o comando para rodar)."
    )

    with gr.Tabs():
        # ---------------- Aba 1: Single ----------------
        with gr.Tab("1️⃣ Proteína única + Ligante único"):
            gr.Markdown("Um complexo só: uma proteína e um ligante.")
            with gr.Row():
                with gr.Column():
                    s_job = gr.Textbox(label="Nome do job", value="meu_complexo", placeholder="ex: EGFR_gefitinib")
                    s_protein_seq = gr.Textbox(label="Sequência da proteína (uma letra por aminoácido)",
                                                lines=6, placeholder="MKTAYIAKQ...")
                    s_ligand_type = gr.Radio(["SMILES", "Código CCD"], value="SMILES", label="Tipo de ligante")
                    s_ligand_value = gr.Textbox(label="SMILES ou código CCD do ligante",
                                                 placeholder="CC(=O)Oc1ccccc1C(=O)O")
                    with gr.Accordion("Opções avançadas", open=False):
                        s_msa = gr.Checkbox(value=True, label="Gerar MSA automaticamente (--use_msa_server) [recomendado]")
                        s_aff = gr.Checkbox(value=True, label="Calcular afinidade de ligação (binding affinity)")
                    s_btn = gr.Button("🚀 Gerar arquivos", variant="primary")
                with gr.Column():
                    s_file = gr.File(label="Download do pacote (.zip)")
                    s_status = gr.Markdown()

            s_btn.click(gerar_single,
                        inputs=[s_job, s_protein_seq, s_ligand_type, s_ligand_value, s_msa, s_aff],
                        outputs=[s_file, s_status])

        # ---------------- Aba 2: Screening ----------------
        with gr.Tab("2️⃣ Screening (1 proteína x vários ligantes)"):
            gr.Markdown(
                "Uma proteína contra uma lista de ligantes. Cole a lista no formato "
                "`nome,SMILES` (um por linha) ou envie um CSV."
            )
            with gr.Row():
                with gr.Column():
                    sc_protein_name = gr.Textbox(label="Nome da proteína", value="minha_proteina")
                    sc_protein_seq = gr.Textbox(label="Sequência da proteína", lines=6,
                                                 placeholder="MKTAYIAKQ...")
                    sc_ligantes_text = gr.Textbox(
                        label="Lista de ligantes (nome,SMILES — um por linha)",
                        lines=8, placeholder=EXEMPLO_LIGANTES)
                    sc_ligantes_csv = gr.File(label="OU envie um CSV (colunas: nome, smiles)", file_types=[".csv"])
                    with gr.Accordion("Opções avançadas", open=False):
                        sc_msa = gr.Checkbox(value=True, label="Gerar MSA automaticamente (--use_msa_server) [recomendado]")
                        sc_aff = gr.Checkbox(value=True, label="Calcular afinidade de ligação (recomendado p/ ranquear ligantes)")
                    sc_btn = gr.Button("🚀 Gerar arquivos", variant="primary")
                with gr.Column():
                    sc_file = gr.File(label="Download do pacote (.zip)")
                    sc_status = gr.Markdown()

            sc_btn.click(gerar_screening,
                         inputs=[sc_protein_name, sc_protein_seq, sc_ligantes_text, sc_ligantes_csv, sc_msa, sc_aff],
                         outputs=[sc_file, sc_status])

        # ---------------- Aba 3: Multidocking ----------------
        with gr.Tab("3️⃣ Multidocking (várias proteínas x vários ligantes)"):
            gr.Markdown(
                "Todas as combinações entre uma lista de proteínas e uma lista de "
                "ligantes. Cole as listas ou envie CSVs."
            )
            with gr.Row():
                with gr.Column():
                    md_proteinas_text = gr.Textbox(
                        label="Lista de proteínas (nome,sequência — um por linha)",
                        lines=6, placeholder=EXEMPLO_PROTEINAS)
                    md_proteinas_csv = gr.File(label="OU envie um CSV (colunas: nome, sequencia)", file_types=[".csv"])
                    md_ligantes_text = gr.Textbox(
                        label="Lista de ligantes (nome,SMILES — um por linha)",
                        lines=6, placeholder=EXEMPLO_LIGANTES)
                    md_ligantes_csv = gr.File(label="OU envie um CSV (colunas: nome, smiles)", file_types=[".csv"])
                    with gr.Accordion("Opções avançadas", open=False):
                        md_msa = gr.Checkbox(value=True, label="Gerar MSA automaticamente (--use_msa_server) [recomendado]")
                        md_aff = gr.Checkbox(value=True, label="Calcular afinidade de ligação")
                    md_btn = gr.Button("🚀 Gerar arquivos (todas as combinações)", variant="primary")
                with gr.Column():
                    md_file = gr.File(label="Download do pacote (.zip)")
                    md_status = gr.Markdown()

            md_btn.click(gerar_multidocking,
                         inputs=[md_proteinas_text, md_proteinas_csv, md_ligantes_text, md_ligantes_csv, md_msa, md_aff],
                         outputs=[md_file, md_status])

    gr.Markdown(
        "---\n"
        "💡 **Depois de gerar:** descompacte o `.zip`, leia o `LEIA-ME.txt` e "
        "rode o script indicado. Precisa ter o Boltz instalado "
        "(`pip install boltz -U`) e, de preferência, uma GPU NVIDIA."
    )


if __name__ == "__main__":
    demo.launch()
