$ErrorActionPreference = "Stop"

Write-Host "DockCADD Windows bootstrap"
Write-Host "This script validates the presence of Java, Open Babel, and AutoDock Vina in PATH."

$missing = @()
foreach ($cmd in @("java", "obabel", "vina")) {
    if (-not (Get-Command $cmd -ErrorAction SilentlyContinue)) {
        $missing += $cmd
    }
}

if ($missing.Count -gt 0) {
    Write-Host "Missing tools: $($missing -join ', ')"
    Write-Host "Install them or use WSL for docking runs."
    exit 1
}

python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip install -e .

Write-Host "DockCADD Python package installed."
Write-Host "Optional: install OpenMM and PDBFixer separately if you want full receptor fixing support."
