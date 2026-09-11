#!/usr/bin/env python3
"""Generate docs/Update_EXHALE_stage2.tex from docs/Update_EXHALE_stage2.md.

The Markdown log is the file the entries are appended to; the TeX file is
its typeset twin in the memo class the earlier stage logs use
(docs/my_memo.cls). Pandoc converts the body; this script supplies the
preamble and the macros Pandoc's LaTeX writer expects.

The appendices are not in the Markdown log. They live in
docs/Update_EXHALE_appendix.tex, which is written by hand, and the postamble
below emits \\appendix and \\input{Update_EXHALE_appendix} after the converted
body; the path is relative to docs/, where latexmk runs. Regenerate after
appending to the log or editing the appendix:

    python3 src/utils/update_log_to_tex.py
    cd docs && latexmk -pdf Update_EXHALE_stage2.tex
"""
import os, re, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MD = os.path.join(ROOT, 'docs', 'Update_EXHALE_stage2.md')
TEX = os.path.join(ROOT, 'docs', 'Update_EXHALE_stage2.tex')

PREAMBLE = r"""%% EXHALE: running update log, stage 2 (from 2026-09-05).
%% GENERATED from Update_EXHALE_stage2.md by src/utils/update_log_to_tex.py;
%% edit the Markdown, not this file.
%% Author: Kwang-Il Seon
\documentclass[english,a4paper]{my_memo}
\synctex=1
\usepackage{listings}
\usepackage{xcolor}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{array}
\usepackage{calc}
\usepackage{url}
\urlstyle{tt}
\usepackage{soul}
\usepackage{babel}
\providecommand{\tightlist}{\setlength{\itemsep}{0pt}\setlength{\parskip}{0pt}}
\providecommand{\pandocbounded}[1]{#1}
\providecommand{\passthrough}[1]{#1}
\newlength{\cslhangindent}
\newlength{\csllabelwidth}
\setcounter{secnumdepth}{0}
\sloppy
\begin{document}

\title{EXHALE: Update Log (stage 2, from 2026-09-05)}
\author{Kwang-Il Seon}
\date{Last updated: \docmoddate}
\maketitle

\tableofcontents
\bigskip

"""
# \appendix then the hand-maintained appendix fragment, then the end of
# the document. The \input path is relative to docs/, the directory
# latexmk runs in.
POSTAMBLE = "\n\\appendix\n\\input{Update_EXHALE_appendix}\n\\end{document}\n"

def main():
    md = open(MD).read()
    # The first line is the document title, carried by \title above.
    md = re.sub(r'^# .*\n', '', md, count=1)
    body = subprocess.run(
        ['pandoc', '-f', 'gfm', '-t', 'latex', '--wrap=preserve',
         '--top-level-division=section'],
        input=md, capture_output=True, text=True, check=True).stdout
    # Pandoc's longtable column specs come as \real{...}; keep them but
    # make every table \small so wide rows fit the memo's text width.
    # A pipe table whose cells never wrapped in the Markdown comes out with
    # a bare l/c/r column spec, which cannot break lines; give any such
    # table with four or more columns equal-width paragraph columns instead.
    def spec(m):
        letters = m.group(1); n = len(letters)
        if n < 4: return m.group(0)
        col = ('>{\\raggedright\\arraybackslash}p{\\dimexpr(\\linewidth-'
               + str(2*n) + '\\tabcolsep)/' + str(n) + '\\relax}')
        return '\\begin{longtable}[]{@{}' + col*n + '@{}}'
    body = re.sub(r'\\begin\{longtable\}\[\]\{@\{\}([lcr]+)@\{\}\}', spec, body)
    body = body.replace('\\begin{longtable}', '{\\small\n\\begin{longtable}')
    body = body.replace('\\end{longtable}', '\\end{longtable}\n}')
    with open(TEX, 'w') as f:
        f.write(PREAMBLE + body + POSTAMBLE)
    print('wrote', TEX, len(body.splitlines()), 'body lines')

if __name__ == '__main__':
    main()
