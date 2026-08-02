"""Collect every PDF the derived case list names into one directory.

The cases live in two places - the reference PDFs come from
generate_reference_outputs.py, the rest from generate_review_pdfs.py - which makes
them awkward to read through. This copies them into review_cases/all/, one
sub-directory per section, and renames each to its case code, so the set sorts into
review order and every filename says what the case demonstrates:

    all/A/A1_most_like_set.pdf        A. choosing the comparable transcript
    all/A/A2_most_like_unset.pdf
    all/B/B1_reduction.pdf            B. domain resolution
    ...

The list itself comes from review_cases/cases.json, written by
generate_review_index.py, so it cannot drift from the index.

Prerequisites, in order:

    python3 tests/generate_reference_outputs.py   # tests/reference_outputs/**
    python3 tests/generate_review_pdfs.py         # review_cases/{rmats,majiq,leafcutter}
    python3 tests/generate_review_index.py        # cases.json + INDEX.md
    python3 tests/collect_review_pdfs.py          # this
"""
import glob
import json
import os
import shutil
import sys

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(TESTS_DIR, '..'))
REF = os.path.join(REPO, 'tests', 'reference_outputs')
CASES = os.path.join(REPO, 'review_cases')
OUT = os.path.join(CASES, 'all')

CASES_JSON = os.path.join(CASES, 'cases.json')


def load_cases():
    """The picks, as derived by tests/generate_review_index.py.

    Read from cases.json rather than restated here: a second list would drift from
    the index the moment the outputs moved.
    """
    if not os.path.exists(CASES_JSON):
        raise SystemExit(
            f'{os.path.relpath(CASES_JSON, REPO)} not found - run '
            f'python3 tests/generate_review_index.py first')
    with open(CASES_JSON) as handle:
        return json.load(handle)


def main():
    if os.path.isdir(OUT):
        shutil.rmtree(OUT)
    os.makedirs(OUT)

    copied, missing = 0, []
    for case in load_cases():
        code, label = case['code'], case['label']
        if case.get('missing'):
            missing.append(f'{code} {label}  (no case in the current outputs)')
            continue
        directory = os.path.join(REPO, case['source'])
        matches = sorted(glob.glob(os.path.join(directory, case['pattern'])))
        if not matches:
            missing.append(f"{code} {label}  (no {case['pattern']} in {case['source']})")
            continue
        # One sub-directory per section, the leading letter of the case code - the
        # same grouping the index renders under.
        section = os.path.join(OUT, code[0])
        os.makedirs(section, exist_ok=True)
        # A pattern may legitimately match several files (the MAJIQ set); everything
        # else resolves to one, and taking the first keeps the naming predictable.
        multiple = len(matches) > 1 and case['pattern'] == '*.pdf'
        for n, src in enumerate(matches if multiple else matches[:1]):
            suffix = f'_{n + 1}' if multiple else ''
            shutil.copy2(src, os.path.join(section, f'{code}{suffix}_{label}.pdf'))
            copied += 1

    print(f'copied {copied} PDFs to {os.path.relpath(OUT, REPO)}/')
    if missing:
        print(f'\n{len(missing)} case(s) had no source file:')
        for item in missing:
            print(f'  {item}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
