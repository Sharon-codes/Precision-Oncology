#!/usr/bin/env python3
"""Restore losslessly chunked GitHub artifacts, checking every SHA-256.

The complete downloadable ZIP already has normal files and needs no restoration.
This module uses only the Python standard library and performs no network access.
Existing outputs are preserved by default: rebuilt gzip timestamps may differ.
--verify-existing additionally checks the exact delivered snapshot bytes.
"""
from pathlib import Path
import hashlib
import json
import os
import tempfile

ROOT = Path(__file__).resolve().parent


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def within(root, name):
    path = (root / name).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError('Artifact path escapes package root')
    return path


def restore(root=ROOT, verify_existing=False):
    root = Path(root).resolve()
    manifest = root / 'large_artifact_manifest.json'
    if not manifest.exists():
        return []
    restored = []
    for item in json.loads(manifest.read_text())['artifacts']:
        target = within(root, item['path'])
        if target.exists():
            if verify_existing and (target.stat().st_size != item['bytes'] or digest(target) != item['sha256']):
                raise ValueError(f'Existing artifact differs from manifest: {item["path"]}')
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        fd, tmp_name = tempfile.mkstemp(prefix=target.name + '.', dir=target.parent)
        temporary = Path(tmp_name)
        try:
            with os.fdopen(fd, 'wb') as out:
                for part in item['parts']:
                    source = within(root, part['path'])
                    if source.stat().st_size != part['bytes'] or digest(source) != part['sha256']:
                        raise ValueError(f'Invalid artifact part: {part["path"]}')
                    with source.open('rb') as inp:
                        for block in iter(lambda: inp.read(1024 * 1024), b''):
                            out.write(block)
            if temporary.stat().st_size != item['bytes'] or digest(temporary) != item['sha256']:
                raise ValueError(f'Assembled artifact failed verification: {item["path"]}')
            temporary.replace(target)
            restored.append(item['path'])
            print(f'Restored and verified {item["path"]}', flush=True)
        finally:
            if temporary.exists():
                temporary.unlink()
    return restored


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify-existing', action='store_true',
                        help='Also require existing files to match the delivered byte snapshot')
    args = parser.parse_args()
    restore(verify_existing=args.verify_existing)
