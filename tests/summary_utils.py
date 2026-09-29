"""Making a run summary's input-source block checkout-independent.

A real run records the full path of every file it read - that is what the block
is for. A summary committed here as a golden reference cannot: the same run from
another checkout would write different paths, and the file is both committed and
compared against on every run.

Dropping the paths is not the answer either. Which fixture a case read, and for a
multi-file format which of its files were present, is exactly what the block
says. So only the part in front of the file name is replaced, and
/somewhere/DOMAS/tests/rmats/SE.MATS.JC.txt is stored as
<tests>/rmats/SE.MATS.JC.txt wherever it was produced.

The run's timestamp is the same problem with no useful part to keep: every run
writes a different one, so it is replaced wholesale. The version line above it is
left alone - it SHOULD fail the comparison when it changes, that being the point
of recording it.

Shared by generate_reference_outputs.py, which applies it to what it writes, and
test_flags.py, which applies it to both sides before comparing.
"""
import os

TESTS_ROOT = '<tests>'
DATE_LINE = 'Date: <run time>'


def portable_source_line(line, tests_dir):
    """One indented input-file line, rewritten relative to tests_dir. Left as it
    is when the path lies outside it - there is nothing stable to rewrite it
    against, and inventing one would hide where the file actually was."""
    path = line.strip()
    if not path.startswith(tests_dir):
        return line
    relative = os.path.relpath(path, tests_dir).replace(os.sep, '/')
    return f'    {TESTS_ROOT}/{relative}'


def portable_summary_lines(lines, tests_dir):
    """A summary's lines with every listed input file, and the run's timestamp,
    made checkout-independent.

    The input block runs from the 'Input source' heading, past its underline,
    through the indented paths, and ends at the blank line before the next
    section. The Date line is replaced wherever it appears.
    """
    out = []
    lines = iter(lines)
    for line in lines:
        if line.startswith('Date: '):
            out.append(DATE_LINE)
            continue
        out.append(line)
        if not line.startswith('Input source'):
            continue
        out.append(next(lines))                 # the ----- underline
        for line in lines:
            out.append(line if not line.strip() else portable_source_line(line, tests_dir))
            if not line.strip():                # blank line ends the block
                break
    return out


def make_summary_portable(path, tests_dir):
    """Rewrite a summary file in place so it can be committed."""
    with open(path) as handle:
        lines = portable_summary_lines(handle.read().splitlines(), tests_dir)
    with open(path, 'w') as handle:
        handle.write('\n'.join(lines) + '\n')
