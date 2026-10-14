#!/usr/bin/env python3
"""Check owned C++ sources with the pinned formatter; never start a solver."""
from pathlib import Path
import argparse
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--formatter', default='clang-format', help='clang-format 18.1.8 executable')
    args = parser.parse_args()
    formatter = shutil.which(args.formatter)
    if formatter is None:
        parser.error('Provide clang-format 18.1.8 through PATH or --formatter.')
    version = subprocess.check_output([formatter, '--version'], text=True).strip()
    if version != 'clang-format version 18.1.8':
        parser.error('Expected clang-format version 18.1.8; found ' + version)
    files = sorted(
        p for p in (ROOT / 'src').rglob('*')
        if p.is_file() and p.suffix in {'.C', '.H', '.cpp'} and p.name != 'PatekNative.C'
    )
    subprocess.run([formatter, '--dry-run', '--Werror', *map(str, files)], cwd=ROOT, check=True)
    print(f'Format check passed: {len(files)} C++ files; {version}.')


if __name__ == '__main__':
    main()
