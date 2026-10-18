"""Read-only terminal monitor. Never starts fits or modifies checkpoints."""
import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import time


def read_json(path):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return None


def snapshot(root, total=200):
    done = 0
    partial = 0
    active = 0
    failed = 0
    pilot_folds = 0
    finished_times = []
    for rep_id in range(total):
        directory = root / f'rep_{rep_id:03d}'
        completion = read_json(directory / 'COMPLETE.json')
        if completion:
            result = completion.get('result', {})
            if result.get('rep_id') == rep_id and result.get('identity') is False:
                done += 1
                finished_times.append((directory / 'COMPLETE.json').stat().st_mtime)
                if rep_id == 0:
                    pilot_folds = 10
                continue
        attempts = sorted(directory.glob('attempt_*'),
                          key=lambda p: p.stat().st_mtime, reverse=True)
        if not attempts:
            continue
        attempt = attempts[0]
        failure = read_json(attempt / 'FAILURE.json')
        failed += int(failure is not None)
        progress = read_json(attempt / 'progress.json') or {}
        folds = max(0, min(10, int(progress.get('completed_folds', 0))))
        partial += folds
        if rep_id == 0:
            pilot_folds = folds
        lock = read_json(directory / 'CONTROLLER.lock')
        if lock:
            try:
                os.kill(int(lock['pid']), 0)
                active += 1
            except (OSError, ValueError, KeyError):
                pass
    phase = 'WAITING'
    elapsed = None
    eta = None
    lock = read_json(root / 'CONTROLLER.lock')
    if lock:
        try:
            pid = int(lock['pid'])
            cmd = Path(f'/proc/{pid}/cmdline').read_bytes().replace(b'\0', b' ')
            phase = 'PILOT' if b' pilot ' in cmd else 'FULL RUN'
            # Linux starttime is field 22; split after the parenthesized comm.
            fields = Path(f'/proc/{pid}/stat').read_text().rsplit(')', 1)[1].split()
            uptime = float(Path('/proc/uptime').read_text().split()[0])
            elapsed = max(0., uptime - float(fields[19]) / os.sysconf('SC_CLK_TCK'))
            started = time.time() - elapsed
            new_done = sum(t >= started for t in finished_times)
            if phase == 'FULL RUN' and new_done >= 2:
                eta = elapsed * (total - done) / new_done
        except (OSError, ValueError, KeyError, IndexError):
            pass
    summary = read_json(root / 'SUMMARY.json')
    complete = bool(summary and summary.get('status') == 'complete'
                    and summary.get('B') == total and done == total)
    if complete:
        phase = 'COMPLETE'
        eta = 0.
    return dict(phase=phase, done=done, total=total, folds=10*done+partial,
                total_folds=10*total, pilot_folds=pilot_folds, active=active,
                failed=failed, elapsed=elapsed, eta=eta, complete=complete)


def duration(seconds):
    if seconds is None:
        return '--'
    minutes, sec = divmod(int(seconds), 60)
    hours, minutes = divmod(minutes, 60)
    return f'{hours:02d}:{minutes:02d}:{sec:02d}'


def render(state):
    pilot = state['phase'] == 'PILOT'
    units = state['pilot_folds'] if pilot else state['folds']
    total = 10 if pilot else state['total_folds']
    fraction = units / total
    filled = min(16, int(fraction * 16))
    bar = '#' * filled + '-' * (16-filled)
    label = f'pilot folds {units}/10' if pilot else f'folds {units}/{total}'
    return (f"{time.strftime('%H:%M:%S')} {state['phase']} [{bar}] "
            f"{100*fraction:5.1f}% {label} | reps {state['done']}/{state['total']} "
            f"| active {state['active']} | failed {state['failed']} "
            f"| elapsed {duration(state['elapsed'])} | ETA~ {duration(state['eta'])}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('/root/pdgeo/private_data/geo_permutations_cloud_v1'))
    parser.add_argument('--interval', type=int, default=10)
    parser.add_argument('--once', action='store_true')
    args = parser.parse_args()
    if args.interval < 1 or not args.root.is_dir():
        parser.error('An existing run directory and positive interval are required')
    print('Read-only: equal-fold progress is not percent of runtime. ETA uses new full-run completions, excludes resumed checkpoints, and remains approximate. Ctrl-C closes only this monitor.', flush=True)
    try:
        while True:
            state = snapshot(args.root)
            line = render(state)
            if sys.stdout.isatty():
                width = max(20, shutil.get_terminal_size((160, 24)).columns-1)
                print('\r' + line[:width].ljust(width), end='', flush=True)
            else:
                print(line, flush=True)
            if args.once or state['complete']:
                if sys.stdout.isatty():
                    print()
                break
            time.sleep(args.interval)
    except KeyboardInterrupt:
        print('\nMonitor closed; computation is unchanged.')


if __name__ == '__main__':
    main()
