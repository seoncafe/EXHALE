"""Plot accepted regional histories and the solver's actual stopping measures."""
import argparse
import csv
from collections import deque
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


class RegionalHistory:
    """Read only complete frames, including while the solver appends a frame."""

    def __init__(self, path):
        self.path = Path(path)
        self.offset = 0
        self.pending = None
        self.first = None
        self.frames = deque(maxlen=32)
        self.trace = []
        self.identity = None

    def update(self):
        if not self.path.exists():
            return
        st = self.path.stat()
        identity = (st.st_dev, st.st_ino)
        if self.identity != identity or st.st_size < self.offset:
            self.__init__(self.path)
            self.identity = identity
        with self.path.open() as stream:
            stream.seek(self.offset)
            while True:
                position = stream.tell()
                line = stream.readline()
                if not line or not line.endswith('\n'):
                    stream.seek(position)
                    break
                if line.startswith('# FRAME '):
                    words = line.split()
                    self.pending = [int(words[2]), float(words[3]), int(words[4]), []]
                elif line.startswith('# END_FRAME') and self.pending is not None:
                    step, time, size, rows = self.pending
                    self.pending = None
                    data = np.asarray(rows)
                    if data.shape != (size, 8) or not np.isfinite(data).all():
                        raise ValueError('Invalid complete regional profile frame')
                    if self.frames and (time <= self.frames[-1][0] or
                                        not np.array_equal(data[:, :2], self.first[:, :2])):
                        raise ValueError('Regional clock or fixed radial regions changed')
                    if self.first is None:
                        self.first = data.copy()
                    self.frames.append((time, data))
                    scale = np.abs(self.first[:, 2:6]).copy()
                    scale[:, 1] = np.sqrt((5/3)*self.first[:, 4]/self.first[:, 2])
                    change = np.abs(data[:, 2:6]-self.first[:, 2:6])/scale
                    self.trace.append((time, {reg: change[data[:, 1] == reg].max(axis=0)
                                             for reg in (1, 2, 3) if np.any(data[:, 1] == reg)}))
                elif line.strip() and not line.startswith('#') and self.pending is not None:
                    self.pending[3].append([float(x) for x in line.split()])
            self.offset = stream.tell()


def criterion_history(path):
    groups = {}
    if not Path(path).exists():
        return groups
    with Path(path).open() as stream:
        for row in csv.DictReader(stream):
            try:
                key = (float(row['elapsed_seconds']), int(row['candidate']),
                       int(row['windows']), int(row['region']))
                measures = np.array([float(row['mean_change']) / float(row['mean_limit']),
                                     float(row['range_width']) / float(row['range_limit'])])
                if not np.isfinite(measures).all():
                    continue
                groups[key] = np.maximum(groups.get(key, np.zeros(2)), measures)
            except (TypeError, ValueError, KeyError, ZeroDivisionError):
                continue  # A writer may still be completing its last CSV line.
    return groups


def draw(history, output):
    figures = []
    if history.first is not None:
        regions = [reg for reg in (1, 2, 3) if np.any(history.first[:, 1] == reg)]
        names = {1: 'Lower atmosphere', 2: 'Transition region', 3: 'Upper atmosphere'}
        labels = [r'$\rho\ [\mathrm{g\ cm^{-3}}]$', r'$v\ [\mathrm{km\ s^{-1}}]$',
                  r'$p\ [\mathrm{erg\ cm^{-3}}]$', r'$T\ [\mathrm{K}]$']
        fig, axes = plt.subplots(4, len(regions), figsize=(5*len(regions), 11), squeeze=False)
        latest_time, latest = history.frames[-1]
        for column, reg in enumerate(regions):
            mask = latest[:, 1] == reg
            r = latest[mask, 0]
            for k in range(4):
                ax = axes[k, column]
                factor = 1e-5 if k == 1 else 1.0
                records = np.array([data[mask, k+2]*factor for _, data in history.frames])
                ax.fill_between(r, records.min(axis=0), records.max(axis=0),
                                color='C0', alpha=.18, label='Recent saved range')
                ax.plot(r, history.first[mask, k+2]*factor, '--', color='0.45', label='Initial state')
                ax.plot(r, latest[mask, k+2]*factor, color='C0', label='Latest state')
                if k != 1:
                    ax.set_yscale('log')
                ax.set_ylabel(labels[k])
                ax.grid(alpha=.2)
                if k == 0:
                    ax.set_title(names[reg])
                if k == 3:
                    ax.set_xlabel(r'$r/R_0$')
        axes[0, 0].legend(fontsize=8)
        fig.suptitle(f'Accepted physical states: t = {latest_time:.6g} s\n'
                     'Shaded range uses recent saved profiles; stopping measures use every accepted step')
        fig.tight_layout(rect=(0, 0, 1, .94))
        figures.append(('profiles', fig))
        evolution, panels = plt.subplots(2, 2, figsize=(10, 7))
        duration = np.array([t-history.trace[0][0] for t, _ in history.trace])
        for k, ax in enumerate(panels.flat):
            for reg in regions:
                ax.plot(duration, [100*entry[reg][k] for _, entry in history.trace],
                        label=names[reg], color=f'C{reg-1}')
            ax.set_title(['Density', 'Velocity', 'Pressure', 'Temperature'][k])
            ax.set_ylabel('Largest change from initial state [%]')
            ax.set_xlabel('Physical duration in this run [s]')
            ax.grid(alpha=.2)
        panels[0, 0].legend(fontsize=8)
        evolution.suptitle('Saved profile changes; velocity scaled by initial sound speed\n'
                           'The stopping criterion uses its own first-window reference')
        evolution.tight_layout(rect=(0, 0, 1, .92))
        figures.append(('evolution', evolution))
    groups = criterion_history(output/'regional_convergence.csv')
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for reg, name in ((1, 'Lower'), (2, 'Transition'), (3, 'Upper')):
        for candidate in sorted({key[1] for key in groups}):
            points = sorted((key[0], val) for key, val in groups.items()
                            if key[3] == reg and key[1] == candidate)
            if points:
                for k, ax in enumerate(axes):
                    ax.plot([p[0] for p in points], [p[1][k] for p in points],
                            '.-', color=f'C{reg-1}', label=f'{name}, candidate {candidate}')
    for ax, label in zip(axes, ['Largest mean change\n/ allowed change',
                                'Largest full range\n/ allowed range']):
        ax.axhline(1, color='k', linestyle='--', label='Configured limit')
        ax.set_yscale('symlog', linthresh=.01)
        ax.set_ylabel(label)
        ax.grid(alpha=.2)
    axes[-1].set_xlabel('Physical duration in this run [s]')
    if groups:
        axes[0].legend(fontsize=8, ncol=2)
    status = output/'regional_convergence_status.txt'
    text = status.read_text() if status.exists() else 'Waiting for regional monitoring output'
    status_line = next((x for x in text.splitlines() if x.startswith('status=')), 'status=collecting')
    resolved = dict(line.split('=', 1) for line in text.splitlines() if '=' in line)
    progress = 'Completed windows: '+resolved.get('completed_windows', '0')+'/'+resolved.get('confirmation_windows', '?')
    progress += '; physical duration: '+resolved.get('elapsed_seconds', '?').strip()+' s'
    sampling_line = next((x for x in text.splitlines() if x.startswith('samples_in_candidate_windows=')), '')
    level = resolved.get('convergence_level', 'not recorded')
    fig.suptitle(status_line+'; convergence level: '+level+'\n'+progress+'\n'+sampling_line+
                 '\nRegional stabilization check (stage 1); stationary convergence is not evaluated')
    fig.tight_layout(rect=(0, 0, 1, .89))
    figures.append(('convergence', fig))
    return figures


def main(arguments=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--regional', action='store_true')
    parser.add_argument('--live', nargs='?', type=float, const=4.0)
    parser.add_argument('--save', type=Path, help='PDF filename prefix; omit for interactive plots')
    parser.add_argument('--output', type=Path, default=Path('output'))
    args = parser.parse_args(arguments)
    if args.live is not None and args.live <= 0:
        parser.error('--live requires a positive refresh interval')
    history = RegionalHistory(args.output/'regional_profile_history.txt')
    while True:
        history.update()
        figures = draw(history, args.output)
        if args.save:
            args.save.parent.mkdir(parents=True, exist_ok=True)
            for suffix, fig in figures:
                fig.savefig(str(args.save)+'_'+suffix+'.pdf')
        if args.live is None:
            if not args.save:
                plt.show()
            break
        if not args.save:
            plt.show(block=False)
        plt.pause(args.live)
        for _, fig in figures:
            plt.close(fig)
