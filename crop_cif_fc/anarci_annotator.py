from __future__ import annotations

import re
import subprocess
import tempfile
from pathlib import Path


ANARCI_DIR = Path("/home/mullagaliamova/ANARCI")
VENV_ACTIVATE = Path("/home/mullagaliamova/envs/anarci-env/bin/activate")
HMMER_PATH = Path("/usr/bin")


def run_anarci(sequence: str) -> tuple[list[tuple[str, int, int]], str]:
    """
    Запускает ANARCI для одной аминокислотной последовательности.

    Возвращает:
        domains: список найденных доменов (chain_type, start, end),
                 может быть несколько (например, VH+VL в одной цепи);
        status: статус выполнения.
    """
    sequence = "".join(sequence.split()).upper()

    if not sequence:
        return [], "EMPTY_SEQUENCE"

    temp_file: str | None = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False, encoding="utf-8",
        ) as file:
            file.write(">temporary_sequence\n")
            file.write(f"{sequence}\n")
            temp_file = file.name

        command = (
            f"cd {ANARCI_DIR} && "
            f". {VENV_ACTIVATE} && "
            f"export PATH=/usr/bin:$PATH && "
            f"ANARCI "
            f"-i {temp_file} "
            f"--hmmerpath {HMMER_PATH} "
            f"-s imgt"
        )

        result = subprocess.run(
            command, shell=True, capture_output=True, text=True, timeout=60,
        )

        output = result.stdout + "\n" + result.stderr

        if result.returncode != 0:
            return [], f"ERROR: {output.strip()}"

        pattern = (
            r"#\|[^|]+\|([^|]+)\|[^|]+\|[^|]+\|"
            r"(\d+)\|(\d+)\|"
        )

        matches = list(re.finditer(pattern, output))

        if not matches:
            return [], "NO_DOMAIN"

        domains = [
            (match.group(1), int(match.group(2)), int(match.group(3)))
            for match in matches
        ]

        return domains, "OK"

    except subprocess.TimeoutExpired:
        return [], "TIMEOUT"

    except Exception as error:
        return [], f"ERROR: {error}"

    finally:
        if temp_file is not None:
            try:
                Path(temp_file).unlink()
            except OSError:
                pass


def is_antibody_chain_type(chain_type: str | None) -> bool:
    """Проверяет, относится ли тип цепи ANARCI к обычной цепи антитела."""
    return chain_type in {"H", "K", "L"}


def chain_role(chain_type: str | None) -> str:
    """Преобразует тип цепи ANARCI в роль цепи."""
    if chain_type == "H":
        return "heavy"

    if chain_type in {"K", "L"}:
        return "light"

    return "other"

