"""
Shared helpers for building and comparing a lightweight, text-only golden
reference for generated PDFs.

Raw PDF bytes aren't diffable as a golden reference: matplotlib embeds a
CreationDate (and similar metadata) in every PDF it writes, so re-running the
exact same analysis produces a byte-different file even when the visible
content is identical. Extracting each page's text with PyPDF2 sidesteps that
- gene/transcript/protein ids, event names, and domain labels are all real
text elements, so a manifest of that text is stable across reruns and still
catches real content regressions (wrong transcript compared, dropped domain,
mislabeled event, ...). It won't catch a purely visual/layout bug that
doesn't change the extracted text - see the README for that tradeoff.
"""
import json
import os

PDF_TEXT_MANIFEST_FILENAME = 'pdf_text_reference.json'


def extract_pdf_pages_text(pdf_path):
    """Return a list with one string per page of extracted text."""
    import PyPDF2
    reader = PyPDF2.PdfReader(pdf_path)
    return [page.extract_text() for page in reader.pages]


def build_pdf_text_manifest(directory):
    """Map every '*_junction_comparison.pdf' filename in directory to its per-page text."""
    manifest = {}
    for filename in sorted(os.listdir(directory)):
        if filename.endswith('_junction_comparison.pdf'):
            manifest[filename] = extract_pdf_pages_text(os.path.join(directory, filename))
    return manifest


def write_pdf_text_manifest(directory):
    """Build the manifest for directory and write it to PDF_TEXT_MANIFEST_FILENAME there."""
    manifest = build_pdf_text_manifest(directory)
    manifest_path = os.path.join(directory, PDF_TEXT_MANIFEST_FILENAME)
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
    return manifest_path
