"""Check the complete release file set and SHA-256 digests (standard library)."""
import hashlib
from pathlib import Path, PurePosixPath
import subprocess

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = "verification/SHA256SUMS"
expected = {}
for row in (ROOT / MANIFEST).read_text(encoding="utf-8").splitlines():
    digest, separator, name = row.partition("  ")
    path = PurePosixPath(name)
    if not separator or path.is_absolute() or ".." in path.parts or "\\" in name or ":" in name or name in expected:
        raise ValueError(f"Unsafe manifest entry: {name}")
    data = ROOT.joinpath(*path.parts).read_bytes()
    if hashlib.sha256(data).hexdigest() != digest:
        raise ValueError(f"SHA-256 mismatch: {name}")
    expected[name] = digest
try:
    top = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], cwd=ROOT, stderr=subprocess.DEVNULL, text=True).strip()
except (OSError, subprocess.CalledProcessError):
    top = ""
if top and Path(top).resolve() == ROOT:
    names = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode("utf-8").split("\0")
    actual = set(filter(None, names)) - {MANIFEST}
else:
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob("*") if p.is_file() and ".git" not in p.relative_to(ROOT).parts and "__pycache__" not in p.parts} - {MANIFEST}
if actual != set(expected):
    raise ValueError(f"File inventory differs: missing={sorted(set(expected)-actual)}, extra={sorted(actual-set(expected))}")
if (ROOT / "paper.pdf").read_bytes() != (ROOT / "paper/paper.pdf").read_bytes():
    raise ValueError("Compatibility PDF differs")
print(f"Verified {len(expected)} release files and the compatibility PDF.")
