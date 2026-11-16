#!/usr/bin/env python3
"""Generate docs/Update_EXHALE_stage<N>.tex from md/Update_EXHALE_stage<N>.md.

The Markdown log is the file the entries are appended to; the TeX file is
its typeset twin in the memo class the earlier stage logs use
(docs/my_memo.cls). Pandoc converts the body; this script supplies the
preamble and the macros Pandoc's LaTeX writer expects. One script for every
stage log: the stage is the argument, and STAGES below holds what differs
between them (the start date in the title, and the appendix).

The appendices are not in the Markdown logs. The code-size appendix lives in
docs/Update_EXHALE_appendix.tex, which is written by hand; for a stage whose
entry in STAGES names it, the postamble emits \\appendix and
\\input{Update_EXHALE_appendix} after the converted body (the path is
relative to docs/, where latexmk runs). It follows the current log: it was
input by the stage 2 twin until 2026-09-25 and by the stage 3 twin since.
Regenerate after appending to a log or editing the appendix:

    python3 src/utils/update_log_to_tex.py 3
    cd docs && latexmk -pdf Update_EXHALE_stage3.tex
"""
import os, re, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# What differs between the stage logs: the first day of each, which the title
# carries, and the appendix fragment its postamble inputs (None for none).
STAGES = {
    2: dict(start='2026-09-05', appendix=None),
    3: dict(start='2026-09-23', appendix='Update_EXHALE_appendix'),
}

PREAMBLE = r"""%% EXHALE: running update log, stage STAGE (from START).
%% GENERATED from Update_EXHALE_stageSTAGE.md by src/utils/update_log_to_tex.py;
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

\title{EXHALE: Update Log (stage STAGE, from START)}
\author{Kwang-Il Seon}
\date{Last updated: \docmoddate}
\maketitle

\tableofcontents
\bigskip

"""

def postamble(appendix):
    """The end of the document: \\appendix and the hand-maintained appendix
    fragment for a stage that carries one, then \\end{document}."""
    if appendix:
        return "\n\\appendix\n\\input{" + appendix + "}\n\\end{document}\n"
    return "\n\\end{document}\n"

def main():
    if len(sys.argv) != 2 or not sys.argv[1].isdigit() \
            or int(sys.argv[1]) not in STAGES:
        sys.exit('usage: python3 src/utils/update_log_to_tex.py <stage>, '
                 'stage one of ' + ', '.join(str(k) for k in sorted(STAGES)))
    stage = int(sys.argv[1])
    MD = os.path.join(ROOT, 'md', 'Update_EXHALE_stage%d.md' % stage)
    TEX = os.path.join(ROOT, 'docs', 'Update_EXHALE_stage%d.tex' % stage)
    head = (PREAMBLE.replace('STAGE', str(stage))
                    .replace('START', STAGES[stage]['start']))
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
        f.write(head + body + postamble(STAGES[stage]['appendix']))
    print('wrote', TEX, len(body.splitlines()), 'body lines')

if __name__ == '__main__':
    main()
