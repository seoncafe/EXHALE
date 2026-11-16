"""Preserve input options that the graphical interface does not edit."""
from pathlib import Path

INTERFACE_KEYS = {
    'Planet name', 'Log10 lower boundary number density', 'Planet radius',
    'Planet mass', 'Equilibrium temperature', 'Orbital distance', 'Escape radius',
    'He/H number ratio', '2D approximate method', 'Parent star mass',
    'Spectrum type', 'Spectrum file', 'Power-law index', 'Photon energy',
    'Use only EUV', 'Log10 of X-ray luminosity', 'Log10 of EUV luminosity',
    'Grid type', 'Numerical flux', 'Reconstruction scheme', 'Include He23S',
    'Load IC', 'Do only PP', 'Force start',
}


def retained_input_options(path):
    """Return unchanged labeled options outside the GUI's editable fields."""
    path = Path(path)
    if not path.exists():
        return []
    retained = []
    for line in path.read_text().splitlines():
        text = line.strip()
        if not text or text.startswith(('#', '[E_low')):
            continue
        label = text.split(':', 1)[0].split('?', 1)[0].split('[', 1)[0].strip()
        if (':' in text or '?' in text) and label not in INTERFACE_KEYS:
            retained.append(line)
    return retained
