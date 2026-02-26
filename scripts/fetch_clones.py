#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import datetime as dt
import json
import os
import sys
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError


API_VERSION = "2022-11-28"


def iso_utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def http_get_json(url: str, token: str) -> tuple[int, dict, str]:
    req = Request(url, method="GET")
    req.add_header("Accept", "application/vnd.github+json")
    req.add_header("Authorization", f"Bearer {token}")
    req.add_header("X-GitHub-Api-Version", API_VERSION)

    try:
        with urlopen(req, timeout=30) as resp:
            body = resp.read().decode("utf-8")
            headers = dict(resp.headers.items())
            return resp.status, headers, body
    except HTTPError as e:
        body = e.read().decode("utf-8", errors="replace")
        headers = dict(e.headers.items()) if e.headers else {}
        return e.code, headers, body
    except URLError as e:
        raise RuntimeError(f"Network error: {e}") from e


def load_json(path: Path) -> dict:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)


def main() -> int:
    # Prefer TRAFFIC_TOKEN if provided, fallback to GITHUB_TOKEN
    token = (os.getenv("TRAFFIC_TOKEN") or "").strip() or (os.getenv("GITHUB_TOKEN") or "").strip()
    if not token:
        print("ERROR: No token found. Set TRAFFIC_TOKEN (recommended) or GITHUB_TOKEN.", file=sys.stderr)
        return 2

    repo_full = (os.getenv("GITHUB_REPOSITORY") or "").strip()  # "owner/repo"
    if not repo_full or "/" not in repo_full:
        print("ERROR: GITHUB_REPOSITORY is missing or invalid.", file=sys.stderr)
        return 2
    owner, repo = repo_full.split("/", 1)

    api_base = (os.getenv("GITHUB_API_URL") or "https://api.github.com").rstrip("/")
    url = f"{api_base}/repos/{owner}/{repo}/traffic/clones?per=day"

    status, headers, body = http_get_json(url, token)
    if status != 200:
        # Helpful hints: REST API may return required permissions in X-Accepted-GitHub-Permissions
        # https://docs.github.com/en/rest/using-the-rest-api/troubleshooting-the-rest-api
        hint = headers.get("X-Accepted-GitHub-Permissions") or headers.get("x-accepted-github-permissions")
        print(f"ERROR: GitHub API returned {status}", file=sys.stderr)
        if hint:
            print(f"Hint (X-Accepted-GitHub-Permissions): {hint}", file=sys.stderr)
        print("Response body:", file=sys.stderr)
        print(body, file=sys.stderr)
        return 1

    payload = json.loads(body)

    # payload example:
    # {
    #   "count": 173,
    #   "uniques": 128,
    #   "clones": [{"timestamp":"2016-10-10T00:00:00Z","count":2,"uniques":1}, ...]
    # }
    clones = payload.get("clones", [])
    fetched_at = iso_utc_now()

    # Upsert by date (UTC day)
    out_dir = Path("traffic")
    daily_json_path = out_dir / "clones_daily.json"
    daily_csv_path = out_dir / "clones_daily.csv"
    last14_path = out_dir / "clones_last14.json"

    daily = load_json(daily_json_path)  # { "YYYY-MM-DD": {"count":..., "uniques":..., ...}, ... }
    if not isinstance(daily, dict):
        daily = {}

    for item in clones:
        ts = item.get("timestamp", "")
        date_utc = ts[:10]  # YYYY-MM-DD from ISO8601 Z timestamp
        if len(date_utc) != 10:
            continue
        daily[date_utc] = {
            "count": int(item.get("count", 0)),
            "uniques": int(item.get("uniques", 0)),
            "timestamp_utc": ts,
            "updated_at_utc": fetched_at,
        }

    # Write daily JSON (sorted by date)
    ordered = {k: daily[k] for k in sorted(daily.keys())}
    atomic_write_text(
        daily_json_path,
        json.dumps(ordered, ensure_ascii=False, indent=2) + "\n",
    )

    # Write CSV
    rows = []
    for date_utc in sorted(ordered.keys()):
        r = ordered[date_utc]
        rows.append([date_utc, r["count"], r["uniques"], r.get("updated_at_utc", "")])

    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_csv = daily_csv_path.with_suffix(".csv.tmp")
    with tmp_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["date_utc", "count", "uniques", "updated_at_utc"])
        w.writerows(rows)
    tmp_csv.replace(daily_csv_path)

    # Store latest 14-day snapshot (useful for debugging)
    snapshot = {
        "fetched_at_utc": fetched_at,
        "source": url,
        "payload": payload,
    }
    atomic_write_text(last14_path, json.dumps(snapshot, ensure_ascii=False, indent=2) + "\n")

    print(f"OK: updated {daily_json_path}, {daily_csv_path}, {last14_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
